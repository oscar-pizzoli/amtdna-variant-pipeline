"""
Block 8 Contract Aggregator: Robust ingestion, verification, and cross-stage
consolidation of upstream contracts (Blocks 2 -> 7) into SampleMasterRecord.
Guarantees strict type safety, string-to-flag parsing, and non-destructive failure handling.
"""

import csv
from dataclasses import dataclass
import hashlib
import json
import logging
import math
from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Tuple, Union

from src.modules.synthesis.models import (
    AlignmentSummary,
    AuthenticitySummary,
    CategorizedWarning,
    ConsensusSummary,
    ContaminationSummary,
    HaplogroupSummary,
    InvariantViolation,
    PreprocessingSummary,
    SampleBiologicalType,
    SampleMasterRecord,
    SequencingQCSummary,
    WarningCategory,
    WarningSeverity,
)

logger = logging.getLogger("mtDNA_Pipeline.ContractAggregator")


class ContractIngestionError(Exception):
    """Raised when an unrecoverable structural corruption or truncated file occurs during ingestion."""
    pass


@dataclass
class ContractFileSet:
    """Resolved physical file paths for a sample's upstream contracts."""
    sample_id: str
    block2: Optional[Path] = None
    block3: Optional[Path] = None
    block4: Optional[Path] = None
    block5: Optional[Path] = None
    block6: Optional[Path] = None
    block7: Optional[Path] = None
    # Companion analytical artifacts
    damage_json: Optional[Path] = None
    consensus_qc_json: Optional[Path] = None
    contamination_json: Optional[Path] = None
    haplogroup_json: Optional[Path] = None
    haplogrep_raw_report: Optional[Path] = None


class ContractAggregator:
    """Ingests, parses, and unifies upstream cryptographic JSON contracts."""

    CODE_CATEGORY_MAP: Dict[str, WarningCategory] = {
        # QC & Coverage
        "LOW_MEAN_DEPTH": WarningCategory.QC_COVERAGE,
        "COVERAGE_DROPOUT": WarningCategory.QC_COVERAGE,
        "EXCESS_MASKED_N": WarningCategory.QC_COVERAGE,
        "EXCESS_MISSING_DATA": WarningCategory.QC_COVERAGE,
        "LOW_BREADTH_COVERAGE": WarningCategory.QC_COVERAGE,
        "HIGH_DUPLICATION_RATE": WarningCategory.QC_COVERAGE,
        "INSUFFICIENT_DATA": WarningCategory.QC_COVERAGE,
        # Authenticity & aDNA Degradation
        "AUTHENTICITY_UNSUPPORTED": WarningCategory.AUTHENTICITY_DAMAGE,
        "DAMAGE_SUSPECT_PRIVATE_MUTATION": WarningCategory.AUTHENTICITY_DAMAGE,
        "TERMINAL_DEL_SUSPECT": WarningCategory.AUTHENTICITY_DAMAGE,
        "INDEL_SOFT_CLIP_SUSPECT": WarningCategory.AUTHENTICITY_DAMAGE,
        "EXCESS_PRIVATE_MUTATIONS": WarningCategory.AUTHENTICITY_DAMAGE,
        "PMD_FILTERING_REQUIRED": WarningCategory.AUTHENTICITY_DAMAGE,
        # Contamination & Exogenous Gating
        "HIGH_CONTAMINATION_CAUTION": WarningCategory.CONTAMINATION,
        "CONTAMINANT_AFFINITY_WARNING": WarningCategory.CONTAMINATION,
        "SUSPECT_HETEROPLASMY_SUBCLADE_REJECTED": WarningCategory.CONTAMINATION,
        "RESTRICT_HETEROPLASMY": WarningCategory.CONTAMINATION,
        "REJECT_CONTAMINATED": WarningCategory.CONTAMINATION,
        "CONTAMINATED_REJECT": WarningCategory.CONTAMINATION,
        # Phylogeny & Node Resolution
        "MISSING_DIAGNOSTIC_MUTATIONS": WarningCategory.PHYLOGENY,
        "MULTIPLE_TOP_HITS_AMBIGUOUS": WarningCategory.PHYLOGENY,
        "HAPLOGROUP_INPUT_DISCORDANCE": WarningCategory.PHYLOGENY,
        "LOW_HAPLOGROUP_CONFIDENCE": WarningCategory.PHYLOGENY,
        "SUBCLADE_HETEROPLASMY_REFINEMENT": WarningCategory.PHYLOGENY,
        "SUBCLADE_FASTA_SPECIFIC": WarningCategory.PHYLOGENY,
        "HETEROPLASMY_SUBCLADE_SUPPORTED": WarningCategory.PHYLOGENY,
        # Pipeline Integrity & Provenance
        "UPSTREAM_CONTRACT_ANOMALY": WarningCategory.PIPELINE_INTEGRITY,
        "MISSING_UPSTREAM_CONTRACT": WarningCategory.PIPELINE_INTEGRITY,
        "MISSING_CONTAMINATION_METRIC": WarningCategory.PIPELINE_INTEGRITY,
        "PROVENANCE_HASH_MISMATCH": WarningCategory.PIPELINE_INTEGRITY,
        "SCHEMA_CORRUPTION": WarningCategory.PIPELINE_INTEGRITY,
    }

    def __init__(
        self,
        operational_files_dir: Optional[Path] = None,
        results_dir: Optional[Path] = None,
    ) -> None:
        target = operational_files_dir or results_dir or Path("operational_files")
        self.operational_files_dir = Path(target)
        self.results_dir = Path(results_dir) if results_dir else self.operational_files_dir

    # --------------------------------------------------------------------------
    # Defensive Type Casting Helpers
    # --------------------------------------------------------------------------
    @staticmethod
    def _safe_float(val: Any) -> Optional[float]:
        """Defensively parse floating point values without raising exceptions on arbitrary strings."""
        if val is None:
            return None
        if isinstance(val, (int, float)):
            if math.isnan(val) or math.isinf(val):
                return None
            return float(val)
        s = str(val).strip()
        if not s or s.upper() in {"NA", "NONE", "NULL", "UNKNOWN", "N/A", "NAN", "INF", "-INF"}:
            return None
        try:
            f = float(s)
            if math.isnan(f) or math.isinf(f):
                return None
            return f
        except (ValueError, TypeError):
            return None

    @staticmethod
    def _safe_int(val: Any) -> Optional[int]:
        """Defensively parse integer values without raising exceptions on arbitrary strings."""
        if val is None:
            return None
        if isinstance(val, int):
            return val
        if isinstance(val, float):
            return int(val) if val.is_integer() else None
        s = str(val).strip()
        if not s or s.upper() in {"NA", "NONE", "NULL", "UNKNOWN", "N/A"}:
            return None
        try:
            if "." in s:
                f = float(s)
                return int(f) if f.is_integer() else None
            return int(s)
        except (ValueError, TypeError):
            return None

    @staticmethod
    def _safe_str(val: Any, default: str = "UNKNOWN") -> str:
        """Defensively resolve strings avoiding 'None' stringifications."""
        if val is None:
            return default
        s = str(val).strip()
        return s if s and s.upper() not in {"NONE", "NULL", "NA"} else default

    # --------------------------------------------------------------------------
    # Contract Discovery & JSON I/O
    # --------------------------------------------------------------------------
    def locate_contracts(self, sample_id: str) -> ContractFileSet:
        """Locate all primary and companion JSON contracts across block directories."""
        cset = ContractFileSet(sample_id=sample_id)

        # Block 2
        for p in [
            self.operational_files_dir / "preprocessing" / sample_id / f"{sample_id}.block2_contract.json",
            self.operational_files_dir / "alignment" / sample_id / f"{sample_id}.block2_contract.json",
            self.results_dir / "preprocessing" / sample_id / f"{sample_id}.block2_contract.json",
            self.results_dir / "alignment" / sample_id / f"{sample_id}.block2_contract.json",
        ]:
            if p.is_file():
                cset.block2 = p
                break

        # Block 3
        for p in [
            self.operational_files_dir / "alignment" / sample_id / f"{sample_id}.block3_contract.json",
            self.operational_files_dir / "alignment" / sample_id / f"{sample_id}.alignment_qc.json",
            self.results_dir / "alignment" / sample_id / f"{sample_id}.block3_contract.json",
            self.results_dir / "alignment" / sample_id / f"{sample_id}.alignment_qc.json",
        ]:
            if p.is_file():
                cset.block3 = p
                break

        # Block 4
        for p in [
            self.operational_files_dir / "damage" / sample_id / f"{sample_id}.block4_contract.json",
            self.results_dir / "damage" / sample_id / f"{sample_id}.block4_contract.json",
        ]:
            if p.is_file():
                cset.block4 = p
                break

        for p in [
            self.operational_files_dir / "damage" / sample_id / f"{sample_id}.damage.json",
            self.operational_files_dir / "damage" / sample_id / f"{sample_id}.damage_model.json",
            self.results_dir / "damage" / sample_id / f"{sample_id}.damage.json",
            self.results_dir / "damage" / sample_id / f"{sample_id}.damage_model.json",
        ]:
            if p.is_file():
                cset.damage_json = p
                break

        # Block 5
        for p in [
            self.operational_files_dir / "contamination" / sample_id / f"{sample_id}.block5_contract.json",
            self.results_dir / "contamination" / sample_id / f"{sample_id}.block5_contract.json",
        ]:
            if p.is_file():
                cset.block5 = p
                break

        for p in [
            self.operational_files_dir / "contamination" / sample_id / f"{sample_id}.schmutzi_summary.json",
            self.results_dir / "contamination" / sample_id / f"{sample_id}.schmutzi_summary.json",
        ]:
            if p.is_file():
                cset.contamination_json = p
                break

        # Block 6
        for p in [
            self.operational_files_dir / "consensus" / sample_id / f"{sample_id}.block6_contract.json",
            self.results_dir / "consensus" / sample_id / f"{sample_id}.block6_contract.json",
        ]:
            if p.is_file():
                cset.block6 = p
                break

        for p in [
            self.operational_files_dir / "consensus" / sample_id / f"{sample_id}.consensus_qc.json",
            self.results_dir / "consensus" / sample_id / f"{sample_id}.consensus_qc.json",
        ]:
            if p.is_file():
                cset.consensus_qc_json = p
                break

        # Block 7 Primary Contract
        for p in [
            self.operational_files_dir / "haplogroup" / sample_id / f"{sample_id}.block7_contract.json",
            self.results_dir / "haplogroup" / sample_id / f"{sample_id}.block7_contract.json",
        ]:
            if p.is_file():
                cset.block7 = p
                break

        # Block 7 Companion JSON
        for p in [
            self.operational_files_dir / "haplogroup" / sample_id / f"{sample_id}.haplogroup.json",
            self.results_dir / "haplogroup" / sample_id / f"{sample_id}.haplogroup.json",
        ]:
            if p.is_file():
                cset.haplogroup_json = p
                break

        # Block 7 Companion Raw TSV Report (HaploGrep3 output)
        for p in [
            self.operational_files_dir / "haplogroup" / sample_id / f"{sample_id}_fasta.txt",
            self.operational_files_dir / "haplogroup" / sample_id / f"{sample_id}.haplogrep.fasta.raw.txt",
            self.operational_files_dir / "haplogroup" / sample_id / f"{sample_id}.haplogrep.txt",
            self.results_dir / "haplogroups" / sample_id / f"{sample_id}.haplogroup_report.txt",
            self.results_dir / "haplogroup" / sample_id / f"{sample_id}_fasta.txt",
        ]:
            if p.is_file():
                cset.haplogrep_raw_report = p
                break

        return cset

    @staticmethod
    def _read_json_safely(path: Optional[Path]) -> Optional[Dict[str, Any]]:
        if path is None or not path.is_file():
            return None
        if path.stat().st_size == 0:
            raise ContractIngestionError(f"Contract artifact '{path}' exists on disk but is completely empty (0 bytes).")
        try:
            with open(path, "r", encoding="utf-8") as fp:
                data = json.load(fp)
        except json.JSONDecodeError as err:
            raise ContractIngestionError(
                f"Corrupted or truncated JSON in contract '{path}': {err}"
            ) from err
        except OSError as err:
            raise ContractIngestionError(
                f"I/O error reading contract '{path}': {err}"
            ) from err

        if not isinstance(data, dict):
            raise ContractIngestionError(
                f"Contract root in '{path}' must be a JSON object (dict), found '{type(data).__name__}'."
            )
        return data

    @staticmethod
    def calculate_file_sha256(path: Path) -> str:
        hasher = hashlib.sha256()
        with open(path, "rb") as fp:
            while chunk := fp.read(65536):
                hasher.update(chunk)
        return hasher.hexdigest()

    def compute_upstream_composite_hash(self, cset: ContractFileSet) -> Tuple[str, Dict[str, str]]:
        composite_hasher = hashlib.sha256()
        individual_hashes: Dict[str, str] = {}

        for b_name in ("block2", "block3", "block4", "block5", "block6", "block7"):
            path: Optional[Path] = getattr(cset, b_name)
            if path and path.is_file():
                h = self.calculate_file_sha256(path)
                individual_hashes[b_name] = h
                composite_hasher.update(f"{b_name}:{h}".encode("utf-8"))
            else:
                individual_hashes[b_name] = "MISSING"
                composite_hasher.update(f"{b_name}:MISSING".encode("utf-8"))

        return composite_hasher.hexdigest(), individual_hashes

    # --------------------------------------------------------------------------
    # Polymorphic Warning Classification
    # --------------------------------------------------------------------------
    def _classify_flag(self, raw_flag: Any, source_block: str = "UNKNOWN") -> CategorizedWarning:
        ctx: Dict[str, Any] = {}
        raw_code: str = "UPSTREAM_CONTRACT_ANOMALY"
        raw_sev: str = "WARNING"
        msg: str = "Unspecified upstream warning"

        if isinstance(raw_flag, str):
            s = raw_flag.strip()
            if s.startswith("{") and s.endswith("}"):
                try:
                    parsed_dict = json.loads(s)
                    if isinstance(parsed_dict, dict):
                        return self._classify_flag(parsed_dict, source_block=source_block)
                except Exception:
                    pass

            sev_prefix_match = re.match(r"^(CRITICAL|WARNING|INFO)\s*:\s*(.*)$", s, re.IGNORECASE)
            if sev_prefix_match:
                raw_sev = sev_prefix_match.group(1).upper()
                remainder = sev_prefix_match.group(2).strip()
                code_split = re.match(r"^([A-Za-z0-9_]+)\s*(?::\s*(.*))?$", remainder)
                if code_split:
                    raw_code = code_split.group(1).upper()
                    msg = code_split.group(2).strip() if code_split.group(2) else remainder
                else:
                    raw_code = remainder.split()[0].upper() if remainder else "UPSTREAM_CONTRACT_ANOMALY"
                    msg = remainder
            else:
                code_split = re.match(r"^([A-Za-z0-9_]+)\s*(?::\s*(.*))?$", s)
                if code_split and (code_split.group(1).upper() in self.CODE_CATEGORY_MAP or code_split.group(1).isupper()):
                    raw_code = code_split.group(1).upper()
                    raw_sev = "WARNING"
                    msg = code_split.group(2).strip() if code_split.group(2) else s
                else:
                    raw_code = "UPSTREAM_CONTRACT_ANOMALY"
                    raw_sev = "WARNING"
                    msg = s

        elif isinstance(raw_flag, dict):
            raw_code = str(raw_flag.get("code", "UPSTREAM_CONTRACT_ANOMALY")).strip()
            raw_sev = str(raw_flag.get("severity", "WARNING")).strip().upper()
            msg = str(raw_flag.get("message", "Unspecified upstream warning"))
            raw_ctx = raw_flag.get("context", {})
            ctx = raw_ctx if isinstance(raw_ctx, dict) else {}
            if "source_block" in raw_flag:
                source_block = str(raw_flag["source_block"]).upper()

        else:
            raw_code = str(getattr(raw_flag, "code", "UPSTREAM_CONTRACT_ANOMALY")).strip()
            raw_sev = str(getattr(raw_flag, "severity", "WARNING")).strip().upper()
            msg = str(getattr(raw_flag, "message", "Unspecified upstream warning"))
            raw_ctx = getattr(raw_flag, "context", {})
            ctx = raw_ctx if isinstance(raw_ctx, dict) else {}
            if hasattr(raw_flag, "source_block"):
                source_block = str(getattr(raw_flag, "source_block")).upper()

        if hasattr(raw_code, "value"):
            raw_code = str(raw_code.value)
        raw_code = str(raw_code).strip()

        if hasattr(raw_sev, "value"):
            raw_sev = str(raw_sev.value)
        raw_sev = str(raw_sev).strip().upper()

        if "CRIT" in raw_sev:
            severity = WarningSeverity.CRITICAL
        elif "INFO" in raw_sev:
            severity = WarningSeverity.INFO
        else:
            severity = WarningSeverity.WARNING

        category = self.CODE_CATEGORY_MAP.get(raw_code)

        if category is None:
            ctx["unmapped_code"] = True
            if any(term in raw_code for term in ("CONTRACT", "SCHEMA", "HASH", "PROVENANCE", "VERSION")):
                category = WarningCategory.PIPELINE_INTEGRITY
            elif any(term in raw_code for term in ("DEPTH", "COVERAGE", "QUAL", "FILTER", "DROPOUT", "MASK", "MISSING_DATA")):
                category = WarningCategory.QC_COVERAGE
            elif any(term in raw_code for term in ("CONTAM", "DISSENT", "HET", "EXOGENOUS")):
                category = WarningCategory.CONTAMINATION
            elif any(term in raw_code for term in ("DAMAGE", "ANCIENT", "PMD", "DEAM", "AUTHENTIC")):
                category = WarningCategory.AUTHENTICITY_DAMAGE
            elif any(term in raw_code for term in ("HAPLOGROUP", "PHYLO", "TREE", "CLADE", "SUBCLADE", "DISCORD")):
                category = WarningCategory.PHYLOGENY
            else:
                category = WarningCategory.PIPELINE_INTEGRITY

        return CategorizedWarning(
            category=category,
            code=raw_code,
            severity=severity,
            source_block=source_block,
            message=msg,
            context=ctx,
        )

    # --------------------------------------------------------------------------
    # Master Aggregation Pipeline
    # --------------------------------------------------------------------------
    def aggregate_sample(
        self,
        sample_id: str,
        contracts: Optional[ContractFileSet] = None,
    ) -> SampleMasterRecord:
        cset = contracts or self.locate_contracts(sample_id)
        composite_hash, contract_hashes = self.compute_upstream_composite_hash(cset)

        warnings: List[CategorizedWarning] = []
        exec_times: Dict[str, float] = {}

        # ----------------------------------------------------------------------
        # Ingest Contract 2 (Preprocessing & Sequencing QC)
        # ----------------------------------------------------------------------
        d_b2 = self._read_json_safely(cset.block2)
        prep_sum = PreprocessingSummary()
        seq_qc_sum = SequencingQCSummary()

        if d_b2:
            raw_r = self._safe_int(d_b2.get("raw_reads_total"))
            ret_r = self._safe_int(d_b2.get("retained_reads_total") or d_b2.get("clean_reads_total"))
            ret_f = self._safe_float(d_b2.get("retention_fraction"))
            col_f = self._safe_float(d_b2.get("collapsed_fraction"))
            out_p = self._safe_str(d_b2.get("output_profile"), "UNKNOWN")

            prep_sum.raw_reads_total = raw_r
            prep_sum.retained_reads_total = ret_r
            prep_sum.retention_fraction = ret_f
            prep_sum.collapsed_fraction = col_f
            prep_sum.output_profile = out_p

            seq_qc_sum.raw_reads_total = raw_r
            seq_qc_sum.clean_reads_total = ret_r
            seq_qc_sum.retention_fraction = ret_f
            seq_qc_sum.collapsed_fraction = col_f

            # Extract adapter trimming with automatic fallback for collapsed libraries
            adapt_pct = self._safe_float(d_b2.get("adapter_trimmed_reads_pct"))
            if adapt_pct is None and col_f is not None:
                adapt_pct = round(col_f * 100.0, 2)
            seq_qc_sum.adapter_trimmed_reads_pct = adapt_pct

            seq_qc_sum.q20_rate_post = self._safe_float(d_b2.get("q20_rate_post"))
            seq_qc_sum.q30_rate_post = self._safe_float(d_b2.get("q30_rate_post"))
            seq_qc_sum.gc_content_pct = self._safe_float(d_b2.get("gc_content_pct"))
            seq_qc_sum.mean_read_length_pre = self._safe_float(d_b2.get("mean_read_length_pre"))
            seq_qc_sum.mean_read_length_post = self._safe_float(d_b2.get("mean_read_length_post"))

            if "execution_time_seconds" in d_b2:
                exec_times["block2"] = self._safe_float(d_b2["execution_time_seconds"]) or 0.0

        # ----------------------------------------------------------------------
        # Ingest Contract 3 (Alignment & Dedup)
        # ----------------------------------------------------------------------
        d_b3 = self._read_json_safely(cset.block3)
        aln_sum = AlignmentSummary()
        if d_b3:
            aln_sum.mapped_reads_mt = self._safe_int(d_b3.get("mapped_reads_mt") or d_b3.get("mapped_reads"))
            raw_endo = self._safe_float(d_b3.get("endogenous_dna_fraction"))
            if raw_endo is not None:
                aln_sum.endogenous_dna_fraction = raw_endo
                aln_sum.endogenous_dna_fraction_source = "BLOCK3"
            aln_sum.mean_depth = self._safe_float(d_b3.get("mean_depth"))
            aln_sum.coverage_breadth_1x = self._safe_float(d_b3.get("coverage_breadth_1x") or d_b3.get("coverage_1x"))
            aln_sum.coverage_breadth_3x = self._safe_float(d_b3.get("coverage_breadth_3x") or d_b3.get("coverage_3x"))
            aln_sum.coverage_breadth_5x = self._safe_float(d_b3.get("coverage_breadth_5x") or d_b3.get("coverage_5x"))
            aln_sum.coverage_breadth_10x = self._safe_float(d_b3.get("coverage_breadth_10x") or d_b3.get("coverage_10x"))
            aln_sum.duplication_rate = self._safe_float(d_b3.get("duplication_rate"))
            if "execution_time_seconds" in d_b3:
                exec_times["block3"] = self._safe_float(d_b3["execution_time_seconds"]) or 0.0

        if aln_sum.endogenous_dna_fraction is None:
            raw_tot = prep_sum.raw_reads_total
            mapped_mt = aln_sum.mapped_reads_mt
            if raw_tot is not None and mapped_mt is not None and raw_tot > 0:
                aln_sum.endogenous_dna_fraction = round(mapped_mt / raw_tot, 6)
                aln_sum.endogenous_dna_fraction_source = "DERIVED_BLOCK2+BLOCK3"

        # ----------------------------------------------------------------------
        # Ingest Contract 4 (Damage & Authenticity)
        # ----------------------------------------------------------------------
        d_b4 = self._read_json_safely(cset.block4)
        auth_sum = AuthenticitySummary()
        if d_b4:
            auth_str = self._safe_str(d_b4.get("authentication_status"), "UNKNOWN").upper()
            auth_sum.authentication_status = auth_str
            if "AUTHENTIC_ANCIENT" in auth_str:
                auth_sum.biological_type = SampleBiologicalType.ANCIENT
            elif any(k in auth_str for k in ("MODERN", "ZERO_DAMAGE", "AUTHENTICITY_UNSUPPORTED")):
                auth_sum.biological_type = SampleBiologicalType.MODERN
            else:
                auth_sum.biological_type = SampleBiologicalType.UNKNOWN

            auth_sum.rescaling_performed = bool(d_b4.get("rescaling_performed", False))
            auth_sum.damage_model_id = self._safe_str(d_b4.get("damage_model_id"), "UNKNOWN")
            auth_sum.library_treatment = self._safe_str(d_b4.get("library_treatment"), "NON_UDG")

            # Primary direct mapping from contract
            auth_sum.median_fragment_length = self._safe_float(d_b4.get("median_fragment_length"))
            auth_sum.mean_fragment_length = self._safe_float(d_b4.get("mean_fragment_length"))
            auth_sum.pmd_fraction_ge_3 = self._safe_float(d_b4.get("pmd_fraction_ge_3"))
            auth_sum.c_to_t_5p_rate_pos1 = self._safe_float(d_b4.get("c_to_t_5p_rate_pos1"))
            auth_sum.g_to_a_3p_rate_pos1 = self._safe_float(d_b4.get("g_to_a_3p_rate_pos1"))

            if "execution_time_seconds" in d_b4:
                exec_times["block4"] = self._safe_float(d_b4["execution_time_seconds"]) or 0.0

        # Fallback inspection on damage_model.json / damage.json using verified keys
        d_dam = self._read_json_safely(cset.damage_json)
        if d_dam:
            fl_obj = d_dam.get("fragment_length") or d_dam.get("fragment_length_metrics") or {}
            if auth_sum.median_fragment_length is None:
                auth_sum.median_fragment_length = self._safe_float(fl_obj.get("median_length_bp") or fl_obj.get("median_length"))
            if auth_sum.mean_fragment_length is None:
                auth_sum.mean_fragment_length = self._safe_float(fl_obj.get("mean_length_bp") or fl_obj.get("mean_length"))

            pmd_obj = d_dam.get("pmd_metrics") or {}
            if auth_sum.pmd_fraction_ge_3 is None:
                auth_sum.pmd_fraction_ge_3 = self._safe_float(pmd_obj.get("fraction_pmd_ge_3"))

            td_obj = d_dam.get("terminal_damage") or d_dam.get("profile", {}).get("damage_5p", {})
            if auth_sum.c_to_t_5p_rate_pos1 is None:
                auth_sum.c_to_t_5p_rate_pos1 = self._safe_float(td_obj.get("c_to_t_5p_pos1_fraction") or td_obj.get("1"))
            if auth_sum.g_to_a_3p_rate_pos1 is None:
                td_3p = d_dam.get("terminal_damage") or d_dam.get("profile", {}).get("damage_3p", {})
                auth_sum.g_to_a_3p_rate_pos1 = self._safe_float(td_3p.get("g_to_a_3p_pos1_fraction") or td_3p.get("1"))

        # ----------------------------------------------------------------------
        # Ingest Contract 5 (Contamination)
        # ----------------------------------------------------------------------
        d_b5 = self._read_json_safely(cset.block5)
        contam_sum = ContaminationSummary()
        if d_b5:
            contam_sum.recommended_action = self._safe_str(d_b5.get("recommended_action"), "INSUFFICIENT_DATA")
            contam_sum.estimation_confidence = self._safe_str(d_b5.get("estimation_confidence"), "UNRELIABLE")
            contam_sum.estimate_source = self._safe_str(d_b5.get("contamination_estimate_source"), "NONE")

            # Primary direct mapping of confidence interval bounds from contract
            contam_sum.ci_95_lower = self._safe_float(d_b5.get("ci_95_lower"))
            contam_sum.ci_95_upper = self._safe_float(d_b5.get("ci_95_upper"))

            raw_cf = d_b5.get("contamination_fraction")
            parsed_cf = self._safe_float(raw_cf)
            if raw_cf is not None and parsed_cf is None:
                contam_sum.contamination_fraction = None
                contam_sum.recommended_action = "INSUFFICIENT_DATA"
                contam_sum.estimation_confidence = "UNRELIABLE"
                warnings.append(
                    CategorizedWarning(
                        category=WarningCategory.PIPELINE_INTEGRITY,
                        code="MISSING_CONTAMINATION_METRIC",
                        severity=WarningSeverity.CRITICAL,
                        source_block="BLOCK5",
                        message=f"Non-numeric contamination_fraction reported in Block 5: '{raw_cf}'.",
                        context={"raw_contamination_fraction": raw_cf},
                    )
                )
            else:
                contam_sum.contamination_fraction = parsed_cf

            chgs = d_b5.get("contaminant_haplogroups") or []
            contam_sum.contaminant_haplogroups = [str(c) for c in chgs] if isinstance(chgs, list) else [str(chgs)]
            if "execution_time_seconds" in d_b5:
                exec_times["block5"] = self._safe_float(d_b5["execution_time_seconds"]) or 0.0
        else:
            contam_sum.contamination_fraction = None
            contam_sum.recommended_action = "INSUFFICIENT_DATA"
            contam_sum.estimation_confidence = "UNRELIABLE"
            contam_sum.estimate_source = "NONE"
            warnings.append(
                CategorizedWarning(
                    category=WarningCategory.PIPELINE_INTEGRITY,
                    code="MISSING_UPSTREAM_CONTRACT",
                    severity=WarningSeverity.CRITICAL,
                    source_block="BLOCK5",
                    message=f"Block 5 contract missing for '{sample_id}'. Contamination unquantified.",
                    context={"missing_contract": "block5"},
                )
            )

        # Fallback to companion JSON if confidence bounds were not present in contract
        if contam_sum.ci_95_lower is None or contam_sum.ci_95_upper is None:
            d_contam_comp = self._read_json_safely(cset.contamination_json)
            if d_contam_comp:
                ci = d_contam_comp.get("ci_95") or d_contam_comp.get("contamination_ci_95")
                if isinstance(ci, (list, tuple)) and len(ci) == 2:
                    if contam_sum.ci_95_lower is None:
                        contam_sum.ci_95_lower = self._safe_float(ci[0])
                    if contam_sum.ci_95_upper is None:
                        contam_sum.ci_95_upper = self._safe_float(ci[1])

        # ----------------------------------------------------------------------
        # Ingest Contract 6 (Consensus & Heteroplasmy)
        # ----------------------------------------------------------------------
        d_b6 = self._read_json_safely(cset.block6)
        cons_sum = ConsensusSummary()
        if d_b6:
            cons_sum.mean_consensus_depth = self._safe_float(d_b6.get("mean_depth"))
            cons_sum.covered_fraction = self._safe_float(d_b6.get("covered_positions_fraction"))
            cons_sum.n_fraction_canonical = self._safe_float(d_b6.get("n_fraction_canonical"))
            cons_sum.consensus_fasta_sha256 = self._safe_str(d_b6.get("consensus_fasta_sha256"), "")
            pass_sites = d_b6.get("heteroplasmy_sites_pass") or []
            cons_sum.heteroplasmy_sites_pass = len(pass_sites) if isinstance(pass_sites, list) else (self._safe_int(pass_sites) or 0)
            if "execution_time_seconds" in d_b6:
                exec_times["block6"] = self._safe_float(d_b6["execution_time_seconds"]) or 0.0
        else:
            warnings.append(
                CategorizedWarning(
                    category=WarningCategory.PIPELINE_INTEGRITY,
                    code="MISSING_UPSTREAM_CONTRACT",
                    severity=WarningSeverity.CRITICAL,
                    source_block="BLOCK6",
                    message=f"Block 6 contract missing for '{sample_id}'. Consensus artifacts unverified.",
                    context={"missing_contract": "block6"},
                )
            )

        d_cqc = self._read_json_safely(cset.consensus_qc_json)
        if d_cqc:
            cons_sum.iupac_fraction = self._safe_float(d_cqc.get("iupac_fraction"))
            cons_sum.heteroplasmy_sites_total = self._safe_int(d_cqc.get("heteroplasmy_sites_total")) or 0
            cons_sum.heteroplasmy_sites_filtered = self._safe_int(d_cqc.get("heteroplasmy_sites_filtered")) or 0
            cons_sum.gating_action_applied = self._safe_str(d_cqc.get("gating_action"), cons_sum.gating_action_applied)

        # ----------------------------------------------------------------------
        # Ingest Contract 7 (Haplogrouping & Cross-Block Audit)
        # ----------------------------------------------------------------------
        d_b7 = self._read_json_safely(cset.block7)
        hg_sum = HaplogroupSummary()
        overall_status_str = "FAIL"

        if d_b7:
            hg_sum.primary_haplogroup = self._safe_str(d_b7.get("primary_haplogroup"), "UNKNOWN")
            hg_sum.resolved_haplogroup = self._safe_str(d_b7.get("resolved_haplogroup"), "UNKNOWN")
            hg_sum.quality_score = self._safe_float(d_b7.get("quality_score"))
            hg_sum.quality_margin = self._safe_float(d_b7.get("quality_margin"))
            hg_sum.haplogroup_status = self._safe_str(d_b7.get("haplogroup_status"), "UNRESOLVED")
            hg_sum.subclade_shift_detected = bool(d_b7.get("subclade_shift_detected", False))
            hg_sum.subclade_shift_supported = d_b7.get("subclade_shift_supported")
            hg_sum.contaminant_affinity = self._safe_str(d_b7.get("contaminant_affinity"), "NOT_APPLICABLE")
            overall_status_str = self._safe_str(d_b7.get("overall_validation_status"), "FAIL").upper()
            if "execution_time_seconds" in d_b7:
                exec_times["block7"] = self._safe_float(d_b7["execution_time_seconds"]) or 0.0

            # Level 1: Estrazione da block7_contract.json
            raw_found = d_b7.get("found_mutations") or []
            hg_sum.found_diagnostic_mutations = [str(m) for m in raw_found] if isinstance(raw_found, list) else []

            raw_missing = (
                d_b7.get("missing_diagnostic_mutations")
                or d_b7.get("missing_mutations")
                or d_b7.get("not_found_mutations")
                or d_b7.get("not_found_polys")
                or []
            )
            hg_sum.missing_diagnostic_mutations = [str(m) for m in raw_missing] if isinstance(raw_missing, list) else []

            # Level 2: Fallback su haplogroup.json
            if not hg_sum.missing_diagnostic_mutations and cset.haplogroup_json:
                d_hg = self._read_json_safely(cset.haplogroup_json)
                if d_hg:
                    miss_cand = d_hg.get("missing_mutations") or d_hg.get("missing_diagnostic_mutations") or []
                    if isinstance(miss_cand, list) and miss_cand:
                        hg_sum.missing_diagnostic_mutations = [str(m) for m in miss_cand]
                    if not hg_sum.found_diagnostic_mutations:
                        fnd_cand = d_hg.get("found_mutations") or []
                        if isinstance(fnd_cand, list):
                            hg_sum.found_diagnostic_mutations = [str(m) for m in fnd_cand]

            # Level 3: Fallback diretto dal report esteso di HaploGrep3
            if not hg_sum.missing_diagnostic_mutations and cset.haplogrep_raw_report and cset.haplogrep_raw_report.is_file():
                try:
                    with open(cset.haplogrep_raw_report, "r", encoding="utf-8") as fp:
                        reader = csv.DictReader(fp, delimiter="\t")
                        for r in reader:
                            sid = r.get("SampleID") or r.get("Sample ID") or r.get("Sample") or ""
                            if sid == sample_id or not sid:
                                raw_nfp = r.get("Not_Found_Polys") or ""
                                cleaned_nfp = re.sub(r"\s*\([^)]*\)", "", raw_nfp)
                                hg_sum.missing_diagnostic_mutations = [t.strip() for t in cleaned_nfp.split() if t.strip()]
                                if not hg_sum.found_diagnostic_mutations:
                                    raw_fp = r.get("Found_Polys") or ""
                                    cleaned_fp = re.sub(r"\s*\([^)]*\)", "", raw_fp)
                                    hg_sum.found_diagnostic_mutations = [t.strip() for t in cleaned_fp.split() if t.strip()]
                                break
                except Exception as err:
                    logger.debug(f"Could not read raw HaploGrep report for '{sample_id}': {err}")

            raw_flags = d_b7.get("flags") or []
            if isinstance(raw_flags, list):
                for f in raw_flags:
                    w = self._classify_flag(f, source_block="BLOCK7")
                    # Arricchisce il messaggio di warning specificando le varianti esatte mancanti
                    if w.code == "MISSING_DIAGNOSTIC_MUTATIONS" and hg_sum.missing_diagnostic_mutations:
                        missing_str = ", ".join(hg_sum.missing_diagnostic_mutations)
                        if missing_str not in w.message:
                            w.message = f"{w.message}: [{missing_str}]"
                        w.context["missing_diagnostic_mutations"] = hg_sum.missing_diagnostic_mutations
                    warnings.append(w)
        else:
            warnings.append(
                CategorizedWarning(
                    category=WarningCategory.PIPELINE_INTEGRITY,
                    code="MISSING_UPSTREAM_CONTRACT",
                    severity=WarningSeverity.CRITICAL,
                    source_block="BLOCK7",
                    message=f"Block 7 contract missing for '{sample_id}'. Lineage unassigned.",
                    context={"missing_contract": "block7"},
                )
            )

        # Cross-Contract Artifact Hash Integrity
        if d_b6 and d_b7:
            b6_sha = self._safe_str(d_b6.get("consensus_fasta_sha256"), "")
            b7_sha = self._safe_str(d_b7.get("fasta_sha256"), "")
            if b6_sha and b7_sha and b6_sha != b7_sha:
                warnings.append(
                    CategorizedWarning(
                        category=WarningCategory.PIPELINE_INTEGRITY,
                        code="PROVENANCE_HASH_MISMATCH",
                        severity=WarningSeverity.CRITICAL,
                        source_block="BLOCK7",
                        message=f"Consensus FASTA SHA256 mismatch: Block 6 ({b6_sha}) vs Block 7 ({b7_sha}).",
                        context={"block6_sha256": b6_sha, "block7_sha256": b7_sha},
                    )
                )

        if d_b7 and overall_status_str == "PASS":
            b7_crit = [w for w in warnings if w.source_block == "BLOCK7" and w.severity == WarningSeverity.CRITICAL]
            if b7_crit:
                warnings.append(
                    CategorizedWarning(
                        category=WarningCategory.PIPELINE_INTEGRITY,
                        code="UPSTREAM_CONTRACT_ANOMALY",
                        severity=WarningSeverity.CRITICAL,
                        source_block="BLOCK7",
                        message="Block 7 reported overall_validation_status='PASS' despite carrying critical flag(s).",
                        context={"critical_flags": [w.code for w in b7_crit]},
                    )
                )

        record = SampleMasterRecord(
            sample_id=sample_id,
            overall_status=overall_status_str,
            preprocessing=prep_sum,
            sequencing_qc=seq_qc_sum,
            alignment=aln_sum,
            authenticity=auth_sum,
            contamination=contam_sum,
            consensus=cons_sum,
            haplogroup=hg_sum,
            warnings=warnings,
            mutations=[],
            provenance=[],
            execution_times_seconds=exec_times,
            upstream_contracts_sha256=contract_hashes,
            upstream_composite_hash=composite_hash,
        )

        violations: List[InvariantViolation] = record.validate_invariants()
        for v in violations:
            record.warnings.append(
                CategorizedWarning(
                    category=WarningCategory.PIPELINE_INTEGRITY,
                    code=v.code,
                    severity=v.severity,
                    source_block="SYNTHESIS",
                    message=f"Cross-block invariant anomaly: {v.message}",
                    context=v.context,
                )
            )

        if any(w.severity == WarningSeverity.CRITICAL for w in record.warnings):
            record.overall_status = "FAIL"
        elif any(w.severity == WarningSeverity.WARNING for w in record.warnings):
            if record.overall_status != "FAIL":
                record.overall_status = "WARN"

        return record