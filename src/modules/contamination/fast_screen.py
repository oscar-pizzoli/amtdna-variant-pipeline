"""Fast Contamination Screening & Triage Engine (Track B / Pre-Screening).

Pure Python/pysam implementation for allele-aware diagnostic screening,
strict separation of genome vs diagnostic depth, robust multi-factor caching,
and conservative escalation triage (FAST_CLEAR, ESCALATE, INSUFFICIENT_DATA).
"""

from collections import Counter
import hashlib
import json
import logging
from pathlib import Path
import re
import statistics
import time
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import pysam

from src.modules.contamination.models import (
    FastScreeningPolicy,
    FastScreenResult,
    ScreeningDecision,
)

logger = logging.getLogger("mtDNA_Pipeline.FastContaminationScreen")


class FastScreenExecutionError(Exception):
    """Raised when BAM parsing, pileup iteration, or screening metrics calculation fails."""
    pass


class FastContaminationScreen:
    """Performs rapid allelic dissent profiling on diagnostic positions and conservative triage."""

    @classmethod
    def sanitize_sample_id(cls, sample_id: str) -> str:
        """Sanitizes sample identifiers for safe filesystem and filename operations."""
        return re.sub(r"[^\w\-.]", "_", sample_id)

    @classmethod
    def compute_sha256(cls, file_path: Path) -> Optional[str]:
        """Calculates SHA-256 digest of a target file."""
        if not file_path.is_file():
            return None
        hasher = hashlib.sha256()
        with open(file_path, "rb") as f:
            while chunk := f.read(65536):
                hasher.update(chunk)
        return hasher.hexdigest()

    @classmethod
    def compute_canonical_panel_hash(cls, panel_data: Union[Dict[int, Dict[str, str]], Set[int]]) -> str:
        """Computes a deterministic SHA-256 hash from in-memory positions and panel alleles."""
        if isinstance(panel_data, dict):
            items = sorted((int(k), f"{v.get('ref', '')}:{v.get('alt', '')}") for k, v in panel_data.items())
        else:
            items = sorted((int(k), "") for k in panel_data)
        canonical_str = "\n".join(f"{pos}:{val}" for pos, val in items)
        return hashlib.sha256(canonical_str.encode("utf-8")).hexdigest()

    @classmethod
    def load_diagnostic_panel(cls, panel_path: Optional[Path]) -> Tuple[Dict[int, Dict[str, str]], Optional[str]]:
        """Loads 1-based diagnostic positions, reference bases, and known alternative alleles from panel TSV."""
        if not panel_path:
            return {}, None

        if not panel_path.is_file():
            raise FastScreenExecutionError(f"Diagnostic panel file not found: {panel_path}")

        panel_sha = cls.compute_sha256(panel_path)
        panel_dict: Dict[int, Dict[str, str]] = {}

        try:
            with open(panel_path, "r", encoding="utf-8") as f:
                for line in f:
                    line_str = line.strip()
                    if not line_str or line_str.startswith("#"):
                        continue
                    parts = line_str.replace(",", "\t").split()
                    if parts and parts[0].isdigit():
                        pos = int(parts[0])
                        ref_b = parts[1].upper() if len(parts) >= 2 and parts[1].upper() in ("A", "C", "G", "T") else ""
                        alt_b = parts[2].upper() if len(parts) >= 3 and parts[2].upper() in ("A", "C", "G", "T") else ""
                        panel_dict[pos] = {"ref": ref_b, "alt": alt_b}
            logger.info(f"Loaded {len(panel_dict)} diagnostic sites from {panel_path.name} (SHA: {panel_sha[:8]}...)")
            return panel_dict, panel_sha
        except Exception as exc:
            raise FastScreenExecutionError(f"Failed loading diagnostic panel {panel_path.name}: {exc}") from exc

    @classmethod
    def compute_genome_mean_depth(cls, bam_path: Path, contig_name: str, contig_length: int) -> float:
        """Calculates true overall sequencing depth across the entire mitochondrial contig."""
        if not bam_path.is_file() or contig_length <= 0:
            return 0.0
        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                if contig_name not in bam.references:
                    return 0.0
                total_bases = sum(
                    col.nsegments
                    for col in bam.pileup(
                        contig=contig_name,
                        start=0,
                        stop=contig_length,
                        truncate=True,
                        stepper="all",  # Obbligatorio per non scartare mate orfani (FLAG 129)
                    )
                )
                return float(total_bases / contig_length)
        except Exception as exc:
            logger.warning(f"Failed calculating full contig depth: {exc}")
            return 0.0

    @classmethod
    def _get_cache_path(cls, cache_dir: Path, sample_id: str) -> Path:
        """Returns sanitized path for cached screening JSON."""
        safe_id = cls.sanitize_sample_id(sample_id)
        return cache_dir / f"{safe_id}.fast_screen.json"

    @classmethod
    def load_cached_result(
        cls,
        cache_path: Path,
        current_bam: Path,
        current_bam_sha: Optional[str],
        current_ref_sha: Optional[str],
        current_panel_sha: Optional[str],
        current_policy_sha: str,
        algorithm_version: str,
        cache_mode: str = "strict",
    ) -> Optional[FastScreenResult]:
        """Loads screening result from cache after strictly verifying provenance fingerprints."""
        if not cache_path.is_file():
            return None

        try:
            with open(cache_path, "r", encoding="utf-8") as f:
                data = json.load(f)

            if data.get("algorithm_version") != algorithm_version:
                logger.info(f"Cache invalidated for {cache_path.name}: Algorithm version mismatch.")
                return None

            if data.get("policy_sha256") != current_policy_sha:
                logger.info(f"Cache invalidated for {cache_path.name}: Policy configuration changed.")
                return None

            if current_panel_sha and data.get("diagnostic_panel_sha256") != current_panel_sha:
                logger.info(f"Cache invalidated for {cache_path.name}: Diagnostic panel SHA-256 mismatch.")
                return None

            if current_ref_sha and data.get("reference_sha256") != current_ref_sha:
                logger.info(f"Cache invalidated for {cache_path.name}: Reference SHA-256 mismatch.")
                return None

            if cache_mode == "strict":
                if current_bam_sha and data.get("input_bam_sha256") != current_bam_sha:
                    logger.info(f"Cache invalidated for {cache_path.name}: BAM SHA-256 mismatch.")
                    return None
            else:
                bam_stat = current_bam.stat()
                cached_mtime = data.get("bam_mtime")
                cached_size = data.get("bam_size")
                if cached_mtime is not None and cached_size is not None:
                    if bam_stat.st_mtime != cached_mtime or bam_stat.st_size != cached_size:
                        logger.info(f"Cache invalidated for {cache_path.name}: BAM mtime/size changed.")
                        return None

            result = FastScreenResult(
                sample_id=data["sample_id"],
                decision=ScreeningDecision(data["decision"]),
                bypass_eligible=bool(data.get("bypass_eligible", False)),
                mean_genome_depth=float(data.get("mean_genome_depth", data.get("mean_filtered_depth", 0.0))),
                mean_diagnostic_depth=float(data.get("mean_diagnostic_depth", data.get("mean_filtered_depth", 0.0))),
                panel_total_sites=int(data.get("panel_total_sites", 0)),
                evaluated_diagnostic_positions=int(data["evaluated_diagnostic_positions"]),
                panel_evaluation_fraction=float(data.get("panel_evaluation_fraction", 0.0)),
                dissenting_positions_count=int(data["dissenting_positions_count"]),
                total_informative_reads=int(data["total_informative_reads"]),
                total_alt_reads=int(data["total_alt_reads"]),
                global_dissent_fraction=float(data.get("global_dissent_fraction", 0.0)),
                max_dissent_vaf=float(data["max_dissent_vaf"]),
                mean_dissent_vaf=float(data["mean_dissent_vaf"]),
                median_dissent_vaf=float(data.get("median_dissent_vaf", 0.0)),
                median_depth_at_dissent_sites=(
                    float(data["median_depth_at_dissent_sites"])
                    if data.get("median_depth_at_dissent_sites") is not None
                    else None
                ),
                applied_thresholds=data.get("applied_thresholds", {}),
                interpretation=data.get("interpretation", "no escalation evidence detected"),
                diagnostic_panel_path=data.get("diagnostic_panel_path"),
                diagnostic_panel_sha256=data.get("diagnostic_panel_sha256"),
                diagnostic_panel_version=data.get("diagnostic_panel_version"),
                policy_sha256=data.get("policy_sha256"),
                algorithm_version=data.get("algorithm_version", algorithm_version),
                input_bam_sha256=data.get("input_bam_sha256"),
                reference_sha256=data.get("reference_sha256"),
                runtime_seconds=float(data.get("runtime_seconds", 0.0)),
            )

            if not result.validate_invariants():
                logger.warning(f"Cached FastScreenResult failed invariant validation: {cache_path.name}")
                return None

            logger.info(f"Reusing verified cached FastScreenResult from {cache_path.name}")
            return result
        except Exception as exc:
            logger.warning(f"Could not read cache file {cache_path.name}: {exc}")
            return None

    @classmethod
    def save_cache_result(
        cls,
        cache_path: Path,
        result: FastScreenResult,
        input_bam: Path,
    ) -> None:
        """Serializes screening result to disk cache with filesystem metadata."""
        try:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            data = result.to_dict()
            bam_stat = input_bam.stat()
            data["bam_mtime"] = bam_stat.st_mtime
            data["bam_size"] = bam_stat.st_size
            with open(cache_path, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2)
        except Exception as exc:
            logger.warning(f"Failed caching screen result to {cache_path.name}: {exc}")

    @classmethod
    def run_screen(
        cls,
        sample_id: str,
        input_bam: Path,
        reference_fasta: Path,
        expected_contig: str,
        policy: FastScreeningPolicy,
        output_dir: Path,
        target_positions: Optional[Union[Dict[int, Dict[str, str]], Set[int]]] = None,
    ) -> Tuple[FastScreenResult, Optional[Path]]:
        """Executes rapid position-by-position allelic dissent screening against diagnostic sites."""
        t0 = time.perf_counter()
        safe_id = cls.sanitize_sample_id(sample_id)
        out_sample_dir = output_dir / safe_id
        out_sample_dir.mkdir(parents=True, exist_ok=True)
        tsv_path = out_sample_dir / f"{safe_id}.dissent_screening.tsv"
        cache_file = cls._get_cache_path(out_sample_dir, sample_id)

        # 1. Risoluzione Pannello Diagnostico
        panel_dict: Optional[Dict[int, Dict[str, str]]] = None
        panel_sha: Optional[str] = None

        if policy.diagnostic_panel_path:
            panel_dict, panel_sha = cls.load_diagnostic_panel(policy.diagnostic_panel_path)
        elif target_positions is not None:
            if isinstance(target_positions, dict):
                panel_dict = target_positions
            else:
                panel_dict = {int(k): {"ref": "", "alt": ""} for k in target_positions}
            panel_sha = cls.compute_canonical_panel_hash(panel_dict)

        policy_sha = policy.compute_policy_sha256()
        bam_sha = cls.compute_sha256(input_bam) if policy.cache_mode == "strict" else None
        ref_sha = cls.compute_sha256(reference_fasta)

        # 2. Controllo Cache Deterministica
        if policy.cache_screen_results:
            cached = cls.load_cached_result(
                cache_path=cache_file,
                current_bam=input_bam,
                current_bam_sha=bam_sha,
                current_ref_sha=ref_sha,
                current_panel_sha=panel_sha,
                current_policy_sha=policy_sha,
                algorithm_version=policy.algorithm_version,
                cache_mode=policy.cache_mode,
            )
            if cached is not None:
                return cached, tsv_path if tsv_path.is_file() else None

        if not input_bam.is_file():
            raise FastScreenExecutionError(f"Input BAM not found: {input_bam}")

        if panel_dict is not None and len(panel_dict) == 0:
            logger.warning(f"[{sample_id}] FastScreen provided an empty diagnostic position set.")
            res = FastScreenResult(
                sample_id=sample_id,
                decision=ScreeningDecision.INSUFFICIENT_DATA,
                bypass_eligible=False,
                interpretation="Diagnostic panel is empty; safety escalation to full Bayesian model required.",
                diagnostic_panel_path=str(policy.diagnostic_panel_path) if policy.diagnostic_panel_path else None,
                diagnostic_panel_sha256=panel_sha,
                diagnostic_panel_version=policy.diagnostic_panel_version,
                policy_sha256=policy_sha,
                algorithm_version=policy.algorithm_version,
                input_bam_sha256=bam_sha,
                reference_sha256=ref_sha,
                runtime_seconds=time.perf_counter() - t0,
            )
            return res, None

        # 3. Calcolo Contig Length & True Genome Depth
        with pysam.FastaFile(str(reference_fasta)) as fa:
            ref_seq = fa.fetch(expected_contig)
            target_contig_length = len(ref_seq)

        true_genome_depth = cls.compute_genome_mean_depth(input_bam, expected_contig, target_contig_length)

        # 4. Scansione Pileup sui Siti Diagnostici
        evaluated_sites = 0
        dissenting_sites = 0
        total_diagnostic_bases = 0
        total_informative_reads = 0
        total_alt_reads = 0
        observed_vafs: List[float] = []
        dissent_depths: List[int] = []

        tsv_rows: List[str] = [
            "pos\tref\tpanel_alt\tmajor_base\tsite_usable_depth\tmajor_count\talt_count\tminor_vaf\tis_dissenting\n"
        ]

        try:
            with pysam.AlignmentFile(str(input_bam), "rb") as bam:
                if expected_contig not in bam.references:
                    raise FastScreenExecutionError(f"Contig '{expected_contig}' missing from BAM references.")

                for pileup_col in bam.pileup(
                    contig=expected_contig,
                    start=0,
                    stop=target_contig_length,
                    min_base_quality=policy.min_base_quality,
                    min_mapping_quality=policy.min_mapping_quality,
                    truncate=True,
                    stepper="all",  # Obbligatorio per catturare tutte le read informative mitocondriali
                ):
                    pos_1based = pileup_col.reference_pos + 1

                    if panel_dict is not None and pos_1based not in panel_dict:
                        continue

                    # Esclusione delezioni/refskip
                    base_counts: Counter = Counter()
                    for read in pileup_col.pileups:
                        if read.is_del or read.is_refskip:
                            continue
                        q_base = read.alignment.query_sequence[read.query_position].upper()
                        if q_base in ("A", "C", "G", "T"):
                            base_counts[q_base] += 1

                    site_usable_depth = sum(base_counts.values())
                    if site_usable_depth < policy.min_depth_per_site:
                        continue

                    evaluated_sites += 1
                    total_diagnostic_bases += site_usable_depth
                    total_informative_reads += site_usable_depth

                    ref_base = ref_seq[pos_1based - 1].upper() if pos_1based <= target_contig_length else "N"
                    panel_alt = panel_dict[pos_1based].get("alt", ".") if panel_dict else "."

                    # Rilevamento della base dominante del campione a questa posizione
                    most_common = base_counts.most_common()
                    major_base, major_count = most_common[0]

                    # Dissenso = letture con base minoritaria (eteroplasmia / contaminazione)
                    alt_count = site_usable_depth - major_count
                    alt_vaf = alt_count / float(site_usable_depth)

                    is_dissenting = (alt_vaf >= policy.noise_vaf_floor) and (alt_count >= policy.min_alt_reads_per_site)
                    if is_dissenting:
                        dissenting_sites += 1
                        total_alt_reads += alt_count
                        observed_vafs.append(alt_vaf)
                        dissent_depths.append(site_usable_depth)

                    tsv_rows.append(
                        f"{pos_1based}\t{ref_base}\t{panel_alt}\t{major_base}\t{site_usable_depth}\t"
                        f"{major_count}\t{alt_count}\t{alt_vaf:.4f}\t{int(is_dissenting)}\n"
                    )

        except Exception as exc:
            if isinstance(exc, FastScreenExecutionError):
                raise
            raise FastScreenExecutionError(f"Fast screen execution failed on {input_bam.name}: {exc}") from exc

        # 5. Calcolo Metriche
        mean_diagnostic_depth = (
            float(total_diagnostic_bases / evaluated_sites) if evaluated_sites > 0 else 0.0
        )
        panel_total_sites = len(panel_dict) if panel_dict is not None else evaluated_sites
        panel_eval_fraction = (
            float(evaluated_sites / panel_total_sites) if panel_total_sites > 0 else 0.0
        )

        max_vaf = max(observed_vafs) if observed_vafs else 0.0
        mean_vaf = (sum(observed_vafs) / len(observed_vafs)) if observed_vafs else 0.0
        median_vaf = statistics.median(observed_vafs) if observed_vafs else 0.0
        median_dissent_depth = statistics.median(dissent_depths) if dissent_depths else None
        global_dissent_fraction = (
            float(total_alt_reads / total_informative_reads) if total_informative_reads > 0 else 0.0
        )

        applied_thresholds = {
            "min_base_quality": policy.min_base_quality,
            "min_mapping_quality": policy.min_mapping_quality,
            "min_depth_per_site": policy.min_depth_per_site,
            "noise_vaf_floor": policy.noise_vaf_floor,
            "min_alt_reads_per_site": policy.min_alt_reads_per_site,
            "escalate_dissenting_sites_min": policy.escalate_dissenting_sites_min,
            "escalate_vaf_threshold": policy.escalate_vaf_threshold,
            "min_depth_for_screen": policy.min_depth_for_screen,
            "min_evaluated_sites": policy.min_evaluated_sites,
        }

        # 6. Triage Decisionale Asimmetrico
        decision = ScreeningDecision.INSUFFICIENT_DATA
        bypass_eligible = False
        interpretation = ""

        if mean_diagnostic_depth < policy.min_depth_for_screen or evaluated_sites < policy.min_evaluated_sites:
            decision = ScreeningDecision.INSUFFICIENT_DATA
            bypass_eligible = False
            interpretation = (
                f"Insufficient diagnostic depth ({mean_diagnostic_depth:.1f}x < {policy.min_depth_for_screen:.1f}x) "
                f"or evaluated diagnostic sites ({evaluated_sites} < {policy.min_evaluated_sites}); "
                "escalation to full Bayesian model required."
            )
        elif (
            dissenting_sites >= policy.escalate_dissenting_sites_min
            or max_vaf >= policy.escalate_vaf_threshold
        ):
            decision = ScreeningDecision.ESCALATE
            bypass_eligible = False
            interpretation = (
                f"Detectable allelic dissent above background threshold ({dissenting_sites} sites >= "
                f"{policy.escalate_dissenting_sites_min}, max VAF {max_vaf:.2%}); "
                "full Bayesian modeling warranted."
            )
        else:
            decision = ScreeningDecision.FAST_CLEAR
            bypass_eligible = True
            interpretation = (
                "No significant allelic dissent detected above background noise floor; "
                "no escalation evidence detected."
            )

        elapsed = time.perf_counter() - t0

        result = FastScreenResult(
            sample_id=sample_id,
            decision=decision,
            bypass_eligible=bypass_eligible,
            mean_genome_depth=true_genome_depth,
            mean_diagnostic_depth=mean_diagnostic_depth,
            panel_total_sites=panel_total_sites,
            evaluated_diagnostic_positions=evaluated_sites,
            panel_evaluation_fraction=panel_eval_fraction,
            dissenting_positions_count=dissenting_sites,
            total_informative_reads=total_informative_reads,
            total_alt_reads=total_alt_reads,
            global_dissent_fraction=global_dissent_fraction,
            max_dissent_vaf=max_vaf,
            mean_dissent_vaf=mean_vaf,
            median_dissent_vaf=median_vaf,
            median_depth_at_dissent_sites=median_dissent_depth,
            applied_thresholds=applied_thresholds,
            interpretation=interpretation,
            diagnostic_panel_path=str(policy.diagnostic_panel_path) if policy.diagnostic_panel_path else None,
            diagnostic_panel_sha256=panel_sha,
            diagnostic_panel_version=policy.diagnostic_panel_version,
            policy_sha256=policy_sha,
            algorithm_version=policy.algorithm_version,
            input_bam_sha256=bam_sha,
            reference_sha256=ref_sha,
            runtime_seconds=elapsed,
        )

        try:
            with open(tsv_path, "w", encoding="utf-8") as f_tsv:
                f_tsv.writelines(tsv_rows)
        except Exception as exc:
            logger.warning(f"Could not write dissent TSV: {exc}")

        if policy.cache_screen_results:
            cls.save_cache_result(cache_file, result, input_bam)

        logger.info(
            f"[{sample_id}] FastScreen complete in {elapsed:.2f}s | Decision: {decision.value} "
            f"(Diagnostic Depth: {mean_diagnostic_depth:.1f}x, Genome Depth: {true_genome_depth:.1f}x, "
            f"Evaluated Sites: {evaluated_sites}/{panel_total_sites}, Dissent Sites: {dissenting_sites}, "
            f"Global Dissent: {global_dissent_fraction:.2%}, Max VAF: {max_vaf:.2%})"
        )

        return result, tsv_path