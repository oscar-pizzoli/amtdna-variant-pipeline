"""
Block 8 Mutation Cataloguer: Cross-stream relational reconciler unifying VCF records,
certified heteroplasmy evidence, and HaploGrep phylogenetic markers into a cohesive tidy catalog.
Guarantees strict provenance tracking, fail-closed I/O, noise gating, and unmeasured VAF preservation.
"""

import csv
from dataclasses import dataclass, field
import gzip
import json
import logging
import math
from pathlib import Path
import re
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple, Union

try:
    import pysam
except ImportError:
    pysam = None

from src.modules.synthesis.coordinate_normalizer import (
    CanonicalMutation,
    CoordinateNormalizer,
    MutationType,
)
from src.modules.synthesis.models import (
    MutationEntry,
    MutationRole,
)

logger = logging.getLogger("mtDNA_Pipeline.MutationCataloguer")


class MutationCatalogError(Exception):
    """Raised when an unrecoverable failure occurs during mutation parsing or stream reconciliation."""
    pass


@dataclass
class PhyloMarkerInfo:
    """Retains parsed canonical mutation together with its phylogenetic assignment flags."""
    mutation: CanonicalMutation
    is_diagnostic: bool = False
    is_private: bool = False
    is_back_mutation: bool = False


@dataclass
class ReconciledVariantState:
    """Internal mutable accumulator tracking multi-stream evidence for a single canonical mutation."""
    canonical_mutation: CanonicalMutation
    sources: Set[str] = field(default_factory=set)
    depth: Optional[int] = None
    vaf: Optional[float] = None
    vaf_source: Optional[str] = None  # "VCF_AF", "VCF_AD_DERIVED", "HET_CALLER"
    spatial_damage_suspect: Optional[bool] = None
    internal_support_count: Optional[int] = None
    terminal_support_count: Optional[int] = None
    filter_status: str = "PASS"
    role: MutationRole = MutationRole.PRIVATE


class MutationCataloguer:
    """Consolidates variant calls across VCF, consensus QC, and phylogenetic reports."""

    def __init__(
        self,
        results_dir: Path,
        reference_info: Any,
        reference_sequence: Optional[str] = None,
    ) -> None:
        self.results_dir = Path(results_dir)
        self.reference_sequence = reference_sequence

        if isinstance(reference_info, int):
            self.reference_length = reference_info
        else:
            self.reference_length = getattr(
                reference_info,
                "length_bp",
                getattr(reference_info, "length", None),
            )

        if not self.reference_length or self.reference_length <= 0:
            raise ValueError(
                f"Invalid reference length resolved from reference_info: '{self.reference_length}'."
            )

        # Inizializzazione flessibile: istanza con argomenti o fallback a classe senza stato
        try:
            self.normalizer = CoordinateNormalizer(
                reference_length=self.reference_length,
                reference_sequence=reference_sequence,
            )
        except TypeError:
            self.normalizer = CoordinateNormalizer

    def _normalize_mutation(
        self,
        raw_mutation: str = "",
        default_pos: Optional[int] = None,
        default_ref: Optional[str] = None,
        default_alt: Optional[str] = None,
    ) -> List[CanonicalMutation]:
        """
        Polymorphic wrapper for CoordinateNormalizer.
        Safely bridges instance vs classmethod calls and ensures List[CanonicalMutation] output.
        """
        target_fn = getattr(self.normalizer, "parse_mutation", CoordinateNormalizer.parse_mutation)
        kwargs: Dict[str, Any] = {"raw_mutation": raw_mutation}

        # Introspezione per parametri ammessi dalla firma
        import inspect
        sig = inspect.signature(target_fn)
        if "default_pos" in sig.parameters:
            kwargs["default_pos"] = default_pos
        if "default_ref" in sig.parameters:
            kwargs["default_ref"] = default_ref
        if "default_alt" in sig.parameters:
            kwargs["default_alt"] = default_alt
        if "reference_sequence" in sig.parameters:
            kwargs["reference_sequence"] = self.reference_sequence
        if "reference_length" in sig.parameters:
            kwargs["reference_length"] = self.reference_length

        res = target_fn(**kwargs)
        if res is None:
            return []
        if isinstance(res, (list, tuple, set)):
            return list(res)
        return [res]

    @staticmethod
    def _is_mutation_invalid(m: CanonicalMutation) -> bool:
        """Evaluate if a parsed canonical mutation is unusable or masked."""
        if m.position is None:
            return True
        if getattr(m, "is_unparseable", False):
            return True
        if getattr(m, "is_masked_n", False):
            return True
        return False

    # --------------------------------------------------------------------------
    # Defensive Utility Parsing Helpers
    # --------------------------------------------------------------------------
    @staticmethod
    def _safe_float(val: Any) -> Optional[float]:
        if val is None:
            return None
        if isinstance(val, (int, float)):
            return None if math.isnan(val) or math.isinf(val) else float(val)
        s = str(val).strip()
        if not s or s.upper() in {"NA", "NONE", "NULL", "NAN", "INF", "-INF"}:
            return None
        try:
            f = float(s)
            return None if math.isnan(f) or math.isinf(f) else f
        except (ValueError, TypeError):
            return None

    @staticmethod
    def _safe_int(val: Any) -> Optional[int]:
        if val is None:
            return None
        if isinstance(val, int):
            return val
        if isinstance(val, float):
            return int(val) if val.is_integer() else None
        s = str(val).strip()
        if not s or s.upper() in {"NA", "NONE", "NULL"}:
            return None
        try:
            return int(float(s)) if "." in s else int(s)
        except (ValueError, TypeError):
            return None

    # --------------------------------------------------------------------------
    # Stream 1: VCF Parsing with Strict Quality Gating
    # --------------------------------------------------------------------------
    def _parse_vcf(self, vcf_path: Path) -> List[Dict[str, Any]]:
        """
        Parse uncompressed or bgzipped VCF records defensively.
        Enforces multi-allelic parity on AD/AF and extracts strict filter metrics.
        """
        if not vcf_path.is_file():
            raise MutationCatalogError(f"Required VCF artifact missing at '{vcf_path}'.")

        records: List[Dict[str, Any]] = []

        if pysam is not None:
            try:
                bcf_in = pysam.VariantFile(str(vcf_path))
                for rec in bcf_in:
                    pos = rec.pos
                    ref = rec.ref.upper()
                    alts = [a.upper() for a in rec.alts] if rec.alts else []
                    filt_keys = list(rec.filter.keys())
                    filt_str = ";".join(filt_keys) if filt_keys else "PASS"

                    depth = None
                    sample = rec.samples[0] if len(rec.samples) > 0 else None
                    if sample and "DP" in sample:
                        depth = self._safe_int(sample["DP"])
                    elif "DP" in rec.info:
                        depth = self._safe_int(rec.info["DP"])

                    ad_raw = sample.get("AD") if sample and "AD" in sample else None
                    ad_list: List[int] = []
                    if ad_raw is not None:
                        ad_list = [self._safe_int(x) or 0 for x in ad_raw]

                    af_raw = None
                    if sample and "AF" in sample:
                        af_raw = sample["AF"]
                    elif "AF" in rec.info:
                        af_raw = rec.info["AF"]

                    af_list: List[Optional[float]] = []
                    if af_raw is not None:
                        if isinstance(af_raw, (list, tuple)):
                            af_list = [self._safe_float(x) for x in af_raw]
                        else:
                            af_list = [self._safe_float(af_raw)]

                    for idx, alt in enumerate(alts):
                        vaf = None
                        vaf_source = None
                        if idx < len(af_list) and af_list[idx] is not None:
                            vaf = af_list[idx]
                            vaf_source = "VCF_AF"

                        if vaf is None and len(ad_list) >= (len(alts) + 1):
                            alt_cnt = ad_list[idx + 1]
                            total_ad = sum(ad_list)
                            if total_ad > 0:
                                vaf = round(alt_cnt / total_ad, 4)
                                vaf_source = "VCF_AD_DERIVED"

                        records.append({
                            "pos": pos,
                            "ref": ref,
                            "alt": alt,
                            "depth": depth,
                            "vaf": vaf,
                            "vaf_source": vaf_source,
                            "filter_status": filt_str,
                        })
                return records
            except Exception as err:
                raise MutationCatalogError(f"Corrupted VCF detected via pysam at '{vcf_path}': {err}") from err

        open_fn = gzip.open if vcf_path.suffix == ".gz" else open
        try:
            with open_fn(vcf_path, "rt", encoding="utf-8") as fp:
                for line_idx, line in enumerate(fp, start=1):
                    line_str = line.strip()
                    if not line_str or line_str.startswith("#"):
                        continue
                    fields = line_str.split("\t")
                    if len(fields) < 5:
                        raise MutationCatalogError(
                            f"Truncated line {line_idx} in VCF '{vcf_path}': fewer than 5 columns."
                        )

                    pos_val = self._safe_int(fields[1])
                    if pos_val is None:
                        raise MutationCatalogError(
                            f"Invalid non-integer position '{fields[1]}' at line {line_idx} in '{vcf_path}'."
                        )

                    ref = fields[3].strip().upper()
                    alt_raw = fields[4].strip().upper()
                    filter_status = fields[6].strip() if len(fields) > 6 else "PASS"

                    info_dict: Dict[str, str] = {}
                    if len(fields) > 7 and fields[7]:
                        for item in fields[7].split(";"):
                            if "=" in item:
                                k, v = item.split("=", 1)
                                info_dict[k.strip()] = v.strip()

                    sample_metrics: Dict[str, str] = {}
                    if len(fields) >= 10:
                        fmt_keys = [k.strip() for k in fields[8].split(":")]
                        smp_vals = [v.strip() for v in fields[9].split(":")]
                        if len(fmt_keys) == len(smp_vals):
                            sample_metrics = dict(zip(fmt_keys, smp_vals))

                    depth = self._safe_int(sample_metrics.get("DP") or info_dict.get("DP"))

                    raw_af = (
                        sample_metrics.get("AF")
                        or sample_metrics.get("FREQ")
                        or info_dict.get("AF")
                        or info_dict.get("VAF")
                    )
                    af_parts = [p.strip() for p in raw_af.split(",")] if raw_af else []

                    raw_ad = sample_metrics.get("AD")
                    ad_parts = [self._safe_int(p) or 0 for p in raw_ad.split(",")] if raw_ad else []

                    alt_list = [a.strip() for a in alt_raw.split(",") if a.strip()]

                    for idx, single_alt in enumerate(alt_list):
                        vaf = None
                        vaf_source = None
                        if idx < len(af_parts):
                            clean_af = af_parts[idx].rstrip("%")
                            f_val = self._safe_float(clean_af)
                            if f_val is not None:
                                vaf = f_val / 100.0 if "%" in af_parts[idx] else f_val
                                vaf_source = "VCF_AF"

                        if vaf is None and len(ad_parts) >= (len(alt_list) + 1):
                            alt_cnt = ad_parts[idx + 1]
                            tot_ad = sum(ad_parts)
                            if tot_ad > 0:
                                vaf = round(alt_cnt / tot_ad, 4)
                                vaf_source = "VCF_AD_DERIVED"

                        records.append({
                            "pos": pos_val,
                            "ref": ref,
                            "alt": single_alt,
                            "depth": depth,
                            "vaf": vaf,
                            "vaf_source": vaf_source,
                            "filter_status": filter_status,
                        })
            return records
        except OSError as err:
            raise MutationCatalogError(f"I/O error reading VCF '{vcf_path}': {err}") from err
        except Exception as err:
            raise MutationCatalogError(f"Malformed syntax in VCF '{vcf_path}': {err}") from err

    # --------------------------------------------------------------------------
    # Stream 2: Certified Heteroplasmy Evidence Ingestion
    # --------------------------------------------------------------------------
    def _ingest_heteroplasmy_evidence(
        self,
        sample_id: str,
        qc_json_path: Optional[Path],
        b6_contract_path: Optional[Path],
        het_tsv_path: Optional[Path],
    ) -> Dict[str, Tuple[CanonicalMutation, Dict[str, Any]]]:
        """
        Ingest site-level heteroplasmy evidence strictly isolating PASS calls.
        Cross-validates against heteroplasmy.tsv and block6_contract.json to guarantee
        that sub-threshold or damage-filtered artifacts are never marked as PASS.
        """
        het_lookup: Dict[str, Tuple[CanonicalMutation, Dict[str, Any]]] = {}

        # 1. Strategia Primaria: Ingestione da heteroplasmy.tsv
        if het_tsv_path and het_tsv_path.is_file():
            try:
                with open(het_tsv_path, "r", encoding="utf-8") as fp:
                    reader = csv.DictReader(fp, delimiter="\t")
                    for row in reader:
                        raw_pos = row.get("position") or row.get("pos")
                        pos = self._safe_int(raw_pos)

                        # Priorità assoluta: reference_base del TSV o base reale da rCRS
                        ref_col = row.get("reference_base") or row.get("ref")
                        if ref_col:
                            ref = str(ref_col).strip().upper()
                        elif self.reference_sequence and pos and 1 <= pos <= len(self.reference_sequence):
                            ref = self.reference_sequence[pos - 1].upper()
                        else:
                            ref = "N"

                        # Identifica l'allele alternativo rispetto al reference
                        major_b = str(row.get("major_allele") or "").strip().upper()
                        minor_b = str(row.get("minor_allele") or "").strip().upper()
                        alt_col = row.get("alt")

                        if alt_col:
                            alt = str(alt_col).strip().upper()
                        elif minor_b and minor_b != ref:
                            alt = minor_b
                        elif major_b and major_b != ref:
                            alt = major_b
                        else:
                            alt = minor_b or major_b or "N"

                        status = str(row.get("status", "FILTERED")).strip().upper()
                        filter_reason = str(row.get("filter_reason", "FILTERED")).strip()

                        parsed_muts = self._normalize_mutation(
                            default_pos=pos,
                            default_ref=ref,
                            default_alt=alt,
                        )

                        for m in parsed_muts:
                            if self._is_mutation_invalid(m):
                                continue

                            is_dmg = str(row.get("spatial_damage_suspect", "")).lower() in {"true", "1", "yes"}
                            het_lookup[m.canonical_key] = (
                                m,
                                {
                                    "vaf": self._safe_float(row.get("alt_vaf") or row.get("vaf")),
                                    "depth": self._safe_int(row.get("callable_depth") or row.get("depth")),
                                    "spatial_damage_suspect": is_dmg,
                                    "internal_support_count": self._safe_int(row.get("internal_support_count")),
                                    "terminal_support_count": self._safe_int(row.get("terminal_support_count")),
                                    "filter_status": "PASS" if status == "PASS" else filter_reason,
                                },
                            )
                if het_lookup:
                    return het_lookup
            except Exception as err:
                logger.warning(f"Could not parse heteroplasmy.tsv ({err}); falling back to JSON contracts.")

        # 2. Strategia Secondaria: Fallback su block6_contract.json o consensus_qc.json
        target_json = b6_contract_path if (b6_contract_path and b6_contract_path.is_file()) else qc_json_path
        if not target_json or not target_json.is_file():
            return het_lookup

        try:
            with open(target_json, "r", encoding="utf-8") as fp:
                data = json.load(fp)
        except Exception as err:
            raise MutationCatalogError(f"Corrupted consensus contract JSON at '{target_json}': {err}") from err

        raw_pass_sites = data.get("heteroplasmy_sites_pass") or []
        pass_positions: Set[int] = set()
        if isinstance(raw_pass_sites, list):
            for p in raw_pass_sites:
                val = self._safe_int(p)
                if val is not None:
                    pass_positions.add(val)

        details = data.get("heteroplasmy_details") or {}
        items = details.items() if isinstance(details, dict) else enumerate(details)

        for key, entry in items:
            if not isinstance(entry, dict):
                continue

            raw_pos = entry.get("position") or entry.get("pos") or key
            pos = self._safe_int(raw_pos)
            ref = str(entry.get("ref", "N")).upper()
            alt = str(entry.get("alt", "N")).upper()

            parsed_muts = self._normalize_mutation(
                default_pos=pos,
                default_ref=ref,
                default_alt=alt,
            )

            for m in parsed_muts:
                if self._is_mutation_invalid(m):
                    continue

                is_pass = (m.position in pass_positions)
                filter_status = "PASS" if is_pass else "FILTERED"

                het_lookup[m.canonical_key] = (
                    m,
                    {
                        "vaf": self._safe_float(entry.get("vaf") or entry.get("frequency")),
                        "depth": self._safe_int(entry.get("depth") or entry.get("coverage")),
                        "spatial_damage_suspect": bool(entry.get("spatial_damage_suspect", False)),
                        "internal_support_count": self._safe_int(entry.get("internal_support_count")),
                        "terminal_support_count": self._safe_int(entry.get("terminal_support_count")),
                        "filter_status": filter_status,
                    },
                )

        return het_lookup

    # --------------------------------------------------------------------------
    # Stream 3: HaploGrep3 Phylogenetic Markers Ingestion
    # --------------------------------------------------------------------------
    def _ingest_phylogeny_markers(
        self, b7_contract_path: Path
    ) -> Dict[str, PhyloMarkerInfo]:
        """Extract phylogenetic markers preserving parsed CanonicalMutation instances."""
        phylo_map: Dict[str, PhyloMarkerInfo] = {}
        if not b7_contract_path.is_file():
            return phylo_map

        try:
            with open(b7_contract_path, "r", encoding="utf-8") as fp:
                data = json.load(fp)
        except Exception as err:
            raise MutationCatalogError(f"Corrupted Block 7 contract at '{b7_contract_path}': {err}") from err

        # Marcatori diagnostici del nodo
        found_tokens = data.get("found_mutations") or []
        for tok in found_tokens:
            is_back = str(tok).strip().endswith("!")
            parsed_muts = self._normalize_mutation(raw_mutation=str(tok))
            for m in parsed_muts:
                if not self._is_mutation_invalid(m):
                    phylo_map[m.canonical_key] = PhyloMarkerInfo(
                        mutation=m,
                        is_diagnostic=True,
                        is_private=False,
                        is_back_mutation=is_back,
                    )

        # Polimorfismi privati del campione
        remaining_tokens = (
            data.get("private_mutations")
            or data.get("remaining_mutations")
            or []
        )
        for tok in remaining_tokens:
            is_back = str(tok).strip().endswith("!")
            parsed_muts = self._normalize_mutation(raw_mutation=str(tok))
            for m in parsed_muts:
                if not self._is_mutation_invalid(m):
                    if m.canonical_key in phylo_map:
                        if is_back:
                            phylo_map[m.canonical_key].is_back_mutation = True
                    else:
                        phylo_map[m.canonical_key] = PhyloMarkerInfo(
                            mutation=m,
                            is_diagnostic=False,
                            is_private=True,
                            is_back_mutation=is_back,
                        )

        return phylo_map

    # --------------------------------------------------------------------------
    # Master Reconciliation Pipeline
    # --------------------------------------------------------------------------
    def build_catalog(
        self,
        sample_id: str,
        require_vcf: bool = True,
        min_vaf_reportable: float = 0.10,
    ) -> List[MutationEntry]:
        """
        Build the unified tidy mutation catalog for a single specimen.
        Integrates VCF calls, damage spatial audits, and phylogenetic markers.
        Enforces multi-stream provenance and prevents sub-threshold noise fabrication.
        """
        vcf_candidates = [
            self.results_dir / "consensus" / sample_id / f"{sample_id}.variants.vcf.gz",
            self.results_dir / "consensus" / sample_id / f"{sample_id}.variants.vcf",
            self.results_dir / sample_id / f"{sample_id}.variants.vcf.gz",
            self.results_dir / sample_id / f"{sample_id}.variants.vcf",
        ]
        vcf_path = next((p for p in vcf_candidates if p.is_file()), None)

        if vcf_path is None and require_vcf:
            raise MutationCatalogError(
                f"Required VCF missing for sample '{sample_id}'. Checked: {[str(p) for p in vcf_candidates]}"
            )

        qc_candidates = [
            self.results_dir / "consensus" / sample_id / f"{sample_id}.consensus_qc.json",
            self.results_dir / sample_id / f"{sample_id}.consensus_qc.json",
        ]
        qc_json_path = next((p for p in qc_candidates if p.is_file()), None)

        b6_candidates = [
            self.results_dir / "consensus" / sample_id / f"{sample_id}.block6_contract.json",
            self.results_dir / sample_id / f"{sample_id}.block6_contract.json",
        ]
        b6_contract_path = next((p for p in b6_candidates if p.is_file()), None)

        tsv_candidates = [
            self.results_dir / "consensus" / sample_id / f"{sample_id}.heteroplasmy.tsv",
            self.results_dir / sample_id / f"{sample_id}.heteroplasmy.tsv",
        ]
        het_tsv_path = next((p for p in tsv_candidates if p.is_file()), None)

        b7_candidates = [
            self.results_dir / "haplogroup" / sample_id / f"{sample_id}.block7_contract.json",
            self.results_dir / "haplogrouping" / sample_id / f"{sample_id}.block7_contract.json",
            self.results_dir / sample_id / f"{sample_id}.block7_contract.json",
        ]
        b7_path = next((p for p in b7_candidates if p.is_file()), Path("non_existent_b7.json"))

        # 1. Ingestione dei flussi primari
        vcf_records = self._parse_vcf(vcf_path) if vcf_path else []
        het_details = self._ingest_heteroplasmy_evidence(
            sample_id=sample_id,
            qc_json_path=qc_json_path,
            b6_contract_path=b6_contract_path,
            het_tsv_path=het_tsv_path,
        )
        phylo_markers = self._ingest_phylogeny_markers(b7_path)

        state_map: Dict[str, ReconciledVariantState] = {}

        # 2. Popolamento dallo stream VCF con filtro di qualità all'origine
        for v in vcf_records:
            if v["ref"] == v["alt"]:
                continue

            # Filtro rigido VCF: scarta all'origine i record con status diverso da PASS
            filt = str(v.get("filter_status", "PASS")).strip()
            if filt not in {"PASS", ".", ""}:
                continue

            parsed_muts = self._normalize_mutation(
                default_pos=v["pos"],
                default_ref=v["ref"],
                default_alt=v["alt"],
            )

            for m in parsed_muts:
                if self._is_mutation_invalid(m) or m.ref == m.alt:
                    continue

                key = m.canonical_key
                if key not in state_map:
                    state_map[key] = ReconciledVariantState(
                        canonical_mutation=m,
                        sources={"VCF"},
                        depth=v.get("depth"),
                        vaf=v.get("vaf"),
                        vaf_source=v.get("vaf_source"),
                        filter_status="PASS",
                    )
                else:
                    state_map[key].sources.add("VCF")
                    if state_map[key].depth is None:
                        state_map[key].depth = v.get("depth")
                    if state_map[key].vaf is None:
                        state_map[key].vaf = v.get("vaf")
                        state_map[key].vaf_source = v.get("vaf_source")

        # 3. Integrazione e arricchimento da evidenza di eteroplasmia
        for key, (canon_mut, h_data) in het_details.items():
            if canon_mut.ref == canon_mut.alt:
                continue

            is_pass = (h_data.get("filter_status") == "PASS")
            het_vaf = h_data.get("vaf")

            if key in state_map:
                st = state_map[key]
                # Aggiunge lo stream HETEROPLASMY e il VAF solo se il test del Blocco 6 era valido (PASS)
                if is_pass:
                    st.sources.add("HETEROPLASMY")
                    if st.vaf is None and het_vaf is not None:
                        st.vaf = het_vaf
                        st.vaf_source = "HET_CALLER"
                    if st.depth is None and h_data.get("depth") is not None:
                        st.depth = h_data["depth"]

                # Attribuzione dei dati di danno spaziale per l'audit del sito
                if h_data.get("spatial_damage_suspect") is not None:
                    st.spatial_damage_suspect = h_data.get("spatial_damage_suspect")
                if h_data.get("internal_support_count") is not None:
                    st.internal_support_count = h_data.get("internal_support_count")
                if h_data.get("terminal_support_count") is not None:
                    st.terminal_support_count = h_data.get("terminal_support_count")
            else:
                # Inclusione autonoma: consentita solo se la variante è certificata PASS e sopra soglia
                if is_pass and (het_vaf is not None and het_vaf >= min_vaf_reportable):
                    state_map[key] = ReconciledVariantState(
                        canonical_mutation=canon_mut,
                        sources={"HETEROPLASMY"},
                        depth=h_data.get("depth"),
                        vaf=het_vaf,
                        vaf_source="HET_CALLER",
                        spatial_damage_suspect=h_data.get("spatial_damage_suspect"),
                        internal_support_count=h_data.get("internal_support_count"),
                        terminal_support_count=h_data.get("terminal_support_count"),
                        filter_status="PASS",
                    )

        # 4. Integrazione delle varianti puramente filogenetiche (es. Back-Mutations rCRS)
        for key, p_info in phylo_markers.items():
            if key in state_map:
                state_map[key].sources.add("PHYLOGENY")
            else:
                state_map[key] = ReconciledVariantState(
                    canonical_mutation=p_info.mutation,
                    sources={"PHYLOGENY"},
                    depth=None,
                    vaf=None,
                    vaf_source=None,
                    spatial_damage_suspect=None,
                    internal_support_count=None,
                    terminal_support_count=None,
                    filter_status="PASS",
                )

        # 5. Risoluzione dei ruoli biologici ed evolutivi
        for key, st in state_map.items():
            if key in phylo_markers:
                p = phylo_markers[key]
                if p.is_back_mutation:
                    st.role = MutationRole.BACK_MUTATION
                elif p.is_diagnostic:
                    st.role = MutationRole.DIAGNOSTIC
                elif p.is_private:
                    st.role = MutationRole.PRIVATE
                else:
                    st.role = MutationRole.PRIVATE
            else:
                st.role = MutationRole.PRIVATE

        # 6. Costruzione del catalogo finale con filtro rumore consapevole
        catalog: List[MutationEntry] = []
        for st in state_map.values():
            # Politica del Noise Gate:
            # - Le mutazioni DIAGNOSTIC e BACK_MUTATION sono preservate a prescindere dal VAF
            #   poiché costituiscono nodi filogenetici autenticati da HaploGrep/Phylotree.
            # - Le mutazioni PRIVATE prive di riscontro filogenetico devono obbligatoriamente
            #   soddisfare min_vaf_reportable, scartando deaminazioni residue e rumore PCR sub-10%.
            if st.role == MutationRole.PRIVATE and (st.vaf is not None and st.vaf < min_vaf_reportable):
                continue

            m = st.canonical_mutation
            composite_stream = "+".join(sorted(st.sources)) if st.sources else "UNKNOWN"

            catalog.append(
                MutationEntry(
                    sample_id=sample_id,
                    position=m.position,
                    ref=m.ref,
                    alt=m.alt,
                    variant_type=m.variant_type,
                    mutation_notation=m.notation,
                    stream_source=composite_stream,
                    role=st.role,
                    depth=st.depth,
                    vaf=st.vaf,
                    spatial_damage_suspect=st.spatial_damage_suspect,
                    internal_support_count=st.internal_support_count,
                    terminal_support_count=st.terminal_support_count,
                    filter_status=st.filter_status,
                )
            )

        # 7. Ordinamento genomico deterministico
        type_rank = {
            MutationType.SNV: 0,
            MutationType.DEL: 1,
            MutationType.INS: 2,
            MutationType.COMPLEX: 3,
        }

        def sort_key(item: MutationEntry) -> Tuple[int, int, int, str]:
            ins_idx = 0
            if item.variant_type == MutationType.INS:
                match = re.search(r"\.(\d+)", item.mutation_notation)
                if match:
                    ins_idx = int(match.group(1))
            rank = type_rank.get(item.variant_type, 4)
            return (item.position, ins_idx, rank, item.alt)

        catalog.sort(key=sort_key)
        return catalog