"""Track B Engine: Read-level Post-Mortem Damage (PMD) scoring, MD-tag validation, and threshold filtering."""

from dataclasses import dataclass, field
import logging
from pathlib import Path
import re
import shutil
import subprocess
import time
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pysam

from src.modules.damage.models import (
    DamageErrorCode,
    DamagePolicy,
    DamageWarningCode,
    LibraryTreatment,
    PMDMetrics,
)

logger = logging.getLogger("mtDNA_Pipeline.PMDEngine")


class PMDScoringError(Exception):
    """Raised when PMD scoring, MD-tag injection, or threshold filtering encounters a domain-specific failure."""
    pass


@dataclass
class PMDExecutionMetadata:
    """Audit container tracking tool parameters, MD injection, and execution times."""

    tool_name: str = "pmdtools"
    tool_version: str = "unknown"
    commands_executed: List[str] = field(default_factory=list)
    md_tags_injected: bool = False
    total_reads_evaluated: int = 0
    reads_passing_threshold: Optional[int] = None
    applied_threshold: Optional[float] = None
    scores_tsv: Optional[Path] = None
    filtered_bam: Optional[Path] = None
    runtime_seconds: float = 0.0
    skipped: bool = False

    def add_command(self, cmd: str) -> None:
        """Appends an executed shell command to the audit trail."""
        self.commands_executed.append(cmd)

    def to_dict(self) -> Dict[str, Any]:
        """Serializes PMD execution metadata to a dictionary."""
        return {
            "tool_name": self.tool_name,
            "tool_version": self.tool_version,
            "commands_executed": self.commands_executed,
            "md_tags_injected": self.md_tags_injected,
            "total_reads_evaluated": self.total_reads_evaluated,
            "reads_passing_threshold": self.reads_passing_threshold,
            "applied_threshold": self.applied_threshold,
            "scores_tsv": str(self.scores_tsv) if self.scores_tsv else None,
            "filtered_bam": str(self.filtered_bam) if self.filtered_bam else None,
            "runtime_seconds": round(self.runtime_seconds, 2),
            "skipped": self.skipped,
        }


class PMDScoringEngine:
    """Orchestrates per-read ancient damage likelihood scoring and optional filtering via PMDtools."""

    @classmethod
    def get_binary_command(cls) -> List[str]:
        """Resolves whether pmdtools is installed as a standalone binary or Python script in PATH."""
        if shutil.which("pmdtools"):
            return ["pmdtools"]
        pmd_script = shutil.which("pmdtools.py")
        if pmd_script:
            return ["python3", pmd_script]
        return []

    @classmethod
    def get_version(cls) -> str:
        """Extracts installed pmdtools version string."""
        cmd = cls.get_binary_command()
        if not cmd:
            return "not_found"
        try:
            res = subprocess.run(cmd + ["--version"], capture_output=True, text=True, check=False)
            raw = f"{res.stdout}\n{res.stderr}".strip()
            match = re.search(r"(?:pmdtools|version)\s*([\d\.]+)", raw, re.IGNORECASE)
            return match.group(1) if match else (raw.split("\n")[0] if raw else "unknown")
        except Exception:
            return "unknown"

    @classmethod
    def verify_and_ensure_md_tags(
        cls,
        bam_path: Path,
        reference_fasta: Path,
        work_dir: Path,
        sample_check_limit: int = 1000,
    ) -> Tuple[Path, bool]:
        """Verifies uniform presence of SAM MD tags across sampled mapped reads, injecting them via calmd if incomplete.

        Args:
            bam_path: Input BAM file.
            reference_fasta: Reference FASTA path (rCRS).
            work_dir: Directory for intermediate calmd BAM storage.
            sample_check_limit: Number of mapped records to inspect.

        Returns:
            Tuple of (effective_bam_path, md_injected_flag).
        """
        mapped_inspected = 0
        md_tag_count = 0

        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                for read in bam:
                    if read.is_unmapped:
                        continue
                    mapped_inspected += 1
                    if read.has_tag("MD"):
                        md_tag_count += 1
                    if mapped_inspected >= sample_check_limit:
                        break
        except Exception as exc:
            raise PMDScoringError(f"Failed to inspect BAM records for MD tag validation: {exc}") from exc

        # If all sampled mapped reads have MD tags (or if BAM has 0 mapped reads), BAM is valid as-is
        if mapped_inspected > 0 and md_tag_count == mapped_inspected:
            return bam_path, False

        if mapped_inspected == 0:
            return bam_path, False

        # Incomplete MD tag coverage detected -> execute samtools calmd (single-threaded for universal compatibility)
        if shutil.which("samtools") is None:
            raise PMDScoringError("BAM has incomplete MD tags, but 'samtools' is not found in PATH to execute calmd.")

        work_dir.mkdir(parents=True, exist_ok=True)
        calmd_bam = work_dir / f"{bam_path.stem}.calmd.bam"
        logger.info(
            f"MD tags missing or incomplete ({md_tag_count}/{mapped_inspected} sampled reads). "
            f"Injecting MD tags via samtools calmd -> {calmd_bam.name}..."
        )

        calmd_cmd = [
            "samtools", "calmd",
            "-b",
            str(bam_path),
            str(reference_fasta),
        ]

        try:
            with open(calmd_bam, "wb") as out_f:
                res = subprocess.run(calmd_cmd, stdout=out_f, stderr=subprocess.PIPE, check=False)
            if res.returncode != 0:
                err_msg = res.stderr.decode("utf-8", "replace").strip()
                if calmd_bam.is_file():
                    calmd_bam.unlink(missing_ok=True)
                raise PMDScoringError(f"samtools calmd execution failed with exit code {res.returncode}: {err_msg}")

            pysam.index(str(calmd_bam))
            return calmd_bam, True

        except (OSError, subprocess.SubprocessError) as exc:
            if calmd_bam.is_file():
                calmd_bam.unlink(missing_ok=True)
            raise PMDScoringError(f"Failed to execute samtools calmd process: {exc}") from exc

    @classmethod
    def parse_pmd_scores_table(cls, scores_tsv: Path) -> PMDMetrics:
        """Parses pmdtools per-read TSV output and computes empirical distribution metrics with parsing diagnostics.

        Args:
            scores_tsv: Path to output TSV generated by pmdtools --print.

        Returns:
            Populated PMDMetrics dataclass instance.
        """
        if not scores_tsv.is_file() or scores_tsv.stat().st_size == 0:
            logger.warning(f"PMD scores file is missing or empty: {scores_tsv.name}")
            return PMDMetrics()

        scores: List[float] = []
        total_lines = 0
        unparsed_lines = 0

        try:
            with open(scores_tsv, "r", encoding="utf-8") as f:
                for line in f:
                    line_str = line.strip()
                    if not line_str or line_str.startswith("#"):
                        continue

                    total_lines += 1
                    parts = line_str.split("\t")

                    # Scan parts backwards/forwards for the numeric PMD score
                    # pmdtools formats: [read_id, score, ...] or [chr, pos, read_id, score] or single float
                    parsed_score: Optional[float] = None
                    if len(parts) >= 2:
                        # Standard format: second column is score
                        try:
                            parsed_score = float(parts[1])
                        except ValueError:
                            # Fallback: check last column if multi-column table
                            try:
                                parsed_score = float(parts[-1])
                            except ValueError:
                                pass
                    else:
                        try:
                            parsed_score = float(parts[0])
                        except ValueError:
                            pass

                    if parsed_score is not None:
                        scores.append(parsed_score)
                    else:
                        unparsed_lines += 1

            # Diagnostic warning if unparsed lines exceed 2% of total rows (ignoring expected 1 header line)
            if total_lines > 1 and (unparsed_lines / total_lines) > 0.02:
                logger.warning(
                    f"[{scores_tsv.name}] High unparsed row fraction: {unparsed_lines}/{total_lines} lines "
                    f"({unparsed_lines / total_lines:.1%}) could not be converted to float scores."
                )

            if not scores:
                logger.info(f"[{scores_tsv.name}] Zero valid numeric scores parsed from PMD table.")
                return PMDMetrics()

            arr = np.array(scores, dtype=np.float64)
            total = len(arr)

            return PMDMetrics(
                total_reads_scored=total,
                min_pmd_score=float(np.min(arr)),
                max_pmd_score=float(np.max(arr)),
                mean_pmd_score=float(np.mean(arr)),
                median_pmd_score=float(np.median(arr)),
                p75_pmd_score=float(np.percentile(arr, 75)),
                p90_pmd_score=float(np.percentile(arr, 90)),
                fraction_pmd_ge_0=float(np.sum(arr >= 0.0) / total),
                fraction_pmd_ge_3=float(np.sum(arr >= 3.0) / total),
                fraction_pmd_ge_5=float(np.sum(arr >= 5.0) / total),
            )

        except Exception as exc:
            logger.warning(f"Failed to calculate distribution metrics from {scores_tsv.name}: {exc}")
            return PMDMetrics()

    @classmethod
    def execute_pmd_pipeline(
        cls,
        sample_id: str,
        input_bam: Path,
        reference_fasta: Path,
        output_dir: Path,
        log_dir: Path,
        policy: DamagePolicy,
        threads: int = 2,
        timeout: Optional[int] = 3600,
        cleanup_on_failure: bool = True,
    ) -> Tuple[Optional[Path], Optional[Path], Optional[Path], PMDMetrics, PMDExecutionMetadata]:
        """Executes full PMD scoring stream, parses score distribution, and optionally generates filtered BAM.

        Args:
            sample_id: Biological sample identifier.
            input_bam: Input BAM file (dedup or rescaled).
            reference_fasta: Target reference FASTA path (rCRS).
            output_dir: Destination folder for PMD deliverables.
            log_dir: Destination folder for stdout/stderr logs.
            policy: Validated DamagePolicy container.
            threads: Number of worker CPU threads for samtools.
            timeout: Maximum allowed execution time in seconds for subprocess pipes.
            cleanup_on_failure: Whether to delete partial/corrupted files if execution fails.

        Returns:
            Tuple of (scores_tsv, pmd_filtered_bam, pmd_filtered_bai, PMDMetrics, PMDExecutionMetadata).

        Raises:
            PMDScoringError: If pipeline processes fail, timeouts expire, or BAM streams break.
        """
        meta = PMDExecutionMetadata()

        # 1. Policy check: Exit gracefully if PMD scoring is disabled
        if not policy.enable_pmd_scoring:
            logger.info(f"[{sample_id}] PMD scoring disabled in policy (enable_pmd_scoring=False). Skipping.")
            meta.skipped = True
            return None, None, None, PMDMetrics(), meta

        pmd_cmd_base = cls.get_binary_command()
        if not pmd_cmd_base:
            raise PMDScoringError("pmdtools executable not found in system PATH.")

        if shutil.which("samtools") is None:
            raise PMDScoringError("samtools not found in PATH (required for high-throughput BAM streaming).")

        meta.tool_version = cls.get_version()

        out_sample_dir = output_dir / sample_id
        out_sample_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)

        stderr_log = log_dir / f"{sample_id}.pmdtools.err"

        # 2. Pre-condition: Uniform MD tag verification and injection
        effective_bam, md_injected = cls.verify_and_ensure_md_tags(
            bam_path=input_bam,
            reference_fasta=reference_fasta,
            work_dir=out_sample_dir / "tmp_md",
        )
        meta.md_tags_injected = md_injected
        if md_injected:
            meta.add_command(f"samtools calmd -b {input_bam.name} {reference_fasta.name}")

        # 3. Resolve pmdtools flags supporting orthogonal single-stranded + UDG treatment
        pmd_args: List[str] = []
        is_ss = (
            policy.library_treatment == LibraryTreatment.SINGLE_STRANDED
            or getattr(policy, "is_single_stranded", False)
        )
        if is_ss:
            pmd_args.append("--ss")

        if policy.library_treatment == LibraryTreatment.UDG_PARTIAL:
            pmd_args.append("--UDGhalf")
        elif policy.library_treatment == LibraryTreatment.UDG_FULL:
            pmd_args.append("--CpG")

        # 4. Stream 1: Per-Read PMD Scoring (samtools view -h -@ N | pmdtools --print [flags] > scores_tsv)
        scores_tsv = out_sample_dir / f"{sample_id}.pmd_scores.tsv"
        scoring_cmd = pmd_cmd_base + ["--printDS"] + pmd_args
        full_scoring_cmd_str = (
            f"samtools view -h -@ {threads} {effective_bam.name} | "
            f"{' '.join(scoring_cmd)} > {scores_tsv.name}"
        )
        meta.add_command(full_scoring_cmd_str)

        logger.info(
            f"[{sample_id}] Executing PMD per-read scoring (treatment={policy.library_treatment.value}, "
            f"flags={pmd_args or 'none'})..."
        )

        t0 = time.perf_counter()
        try:
            with open(stderr_log, "w", encoding="utf-8") as err_log_f:
                p_view = subprocess.Popen(
                    ["samtools", "view", "-h", "-@", str(threads), str(effective_bam)],
                    stdout=subprocess.PIPE,
                    stderr=err_log_f,
                )
                with open(scores_tsv, "w", encoding="utf-8") as out_tsv:
                    p_pmd = subprocess.Popen(
                        scoring_cmd,
                        stdin=p_view.stdout,
                        stdout=out_tsv,
                        stderr=err_log_f,
                    )

                # Allow SIGPIPE propagation by closing parent reference to pipe
                if p_view.stdout:
                    p_view.stdout.close()

                try:
                    p_pmd_code = p_pmd.wait(timeout=timeout)
                    p_view_code = p_view.wait(timeout=timeout)
                except subprocess.TimeoutExpired as exc:
                    p_pmd.kill()
                    p_view.kill()
                    raise PMDScoringError(f"Stream 1 (PMD Scoring) timed out after {timeout} seconds.") from exc

                if p_view_code != 0:
                    raise PMDScoringError(f"samtools view stream failed with exit code {p_view_code}. Check log: {stderr_log}")
                if p_pmd_code != 0:
                    raise PMDScoringError(f"pmdtools --print failed with exit code {p_pmd_code}. Check log: {stderr_log}")

        except (OSError, subprocess.SubprocessError) as exc:
            if cleanup_on_failure and scores_tsv.is_file():
                scores_tsv.unlink(missing_ok=True)
            raise PMDScoringError(f"Stream 1 (PMD Scoring) failed for sample '{sample_id}': {exc}") from exc
        except Exception:
            if cleanup_on_failure and scores_tsv.is_file():
                scores_tsv.unlink(missing_ok=True)
            raise
        finally:
            meta.runtime_seconds = time.perf_counter() - t0

        # Parse continuous score distribution
        pmd_metrics = cls.parse_pmd_scores_table(scores_tsv)
        meta.total_reads_evaluated = pmd_metrics.total_reads_scored
        meta.scores_tsv = scores_tsv

        # 5. Stream 2 (Optional): PMD Threshold Filtering (samtools view | pmdtools --threshold | samtools view -b)
        filtered_bam: Optional[Path] = None
        filtered_bai: Optional[Path] = None

        if policy.enable_pmd_filtering:
            threshold = policy.pmd_threshold
            meta.applied_threshold = threshold
            pmd_metrics.applied_threshold = threshold

            tag_str = f"{int(threshold)}" if float(threshold).is_integer() else f"{threshold}"
            filtered_bam = out_sample_dir / f"{sample_id}.pmd{tag_str}.bam"
            filtered_bai = Path(f"{filtered_bam}.bai")

            filter_cmd = pmd_cmd_base + ["--threshold", str(threshold)] + pmd_args
            full_filter_cmd_str = (
                f"samtools view -h -@ {threads} {effective_bam.name} | "
                f"{' '.join(filter_cmd)} | "
                f"samtools view -b -@ {threads} -o {filtered_bam.name}"
            )
            meta.add_command(full_filter_cmd_str)

            logger.info(f"[{sample_id}] Generating PMD-filtered BAM (threshold >= {threshold}) -> {filtered_bam.name}...")

            try:
                with open(stderr_log, "a", encoding="utf-8") as err_log_f:
                    p1 = subprocess.Popen(
                        ["samtools", "view", "-h", "-@", str(threads), str(effective_bam)],
                        stdout=subprocess.PIPE,
                        stderr=err_log_f,
                    )
                    p2 = subprocess.Popen(
                        filter_cmd,
                        stdin=p1.stdout,
                        stdout=subprocess.PIPE,
                        stderr=err_log_f,
                    )
                    if p1.stdout:
                        p1.stdout.close()

                    with open(filtered_bam, "wb") as b_out:
                        p3 = subprocess.Popen(
                            ["samtools", "view", "-b", "-@", str(threads), "-"],
                            stdin=p2.stdout,
                            stdout=b_out,
                            stderr=err_log_f,
                        )
                        if p2.stdout:
                            p2.stdout.close()

                        try:
                            p3_code = p3.wait(timeout=timeout)
                            p2_code = p2.wait(timeout=timeout)
                            p1_code = p1.wait(timeout=timeout)
                        except subprocess.TimeoutExpired as exc:
                            p3.kill()
                            p2.kill()
                            p1.kill()
                            raise PMDScoringError(f"Stream 2 (PMD Filtering) timed out after {timeout} seconds.") from exc

                if p1_code != 0:
                    raise PMDScoringError(f"samtools view in Stream 2 failed with exit code {p1_code}.")
                if p2_code != 0:
                    raise PMDScoringError(f"pmdtools --threshold failed with exit code {p2_code}.")
                if p3_code != 0:
                    raise PMDScoringError(f"samtools view BAM creation failed with exit code {p3_code}.")

                pysam.index(str(filtered_bam))

                with pysam.AlignmentFile(str(filtered_bam), "rb") as b_check:
                    retained_count = b_check.count(until_eof=True)

                pmd_metrics.reads_retained_post_pmd_filter = retained_count
                meta.reads_passing_threshold = retained_count
                meta.filtered_bam = filtered_bam

                logger.info(f"[{sample_id}] PMD Filtering Completed: {retained_count:,} reads retained post-filter.")

            except (OSError, subprocess.SubprocessError, pysam.SamtoolsError) as exc:
                if cleanup_on_failure:
                    if filtered_bam and filtered_bam.is_file():
                        filtered_bam.unlink(missing_ok=True)
                    if filtered_bai and filtered_bai.is_file():
                        filtered_bai.unlink(missing_ok=True)
                raise PMDScoringError(f"Stream 2 (PMD Filtering) failed for sample '{sample_id}': {exc}") from exc
            except Exception:
                if cleanup_on_failure:
                    if filtered_bam and filtered_bam.is_file():
                        filtered_bam.unlink(missing_ok=True)
                    if filtered_bai and filtered_bai.is_file():
                        filtered_bai.unlink(missing_ok=True)
                raise

        # Output logging
        if pmd_metrics.total_reads_scored > 0:
            logger.info(
                f"[{sample_id}] PMD Scoring Completed | Evaluated: {pmd_metrics.total_reads_scored:,} reads | "
                f"Median Score: {pmd_metrics.median_pmd_score:.2f} | "
                f"Fraction >= 0: {pmd_metrics.fraction_pmd_ge_0:.2%} | "
                f"Fraction >= 3: {pmd_metrics.fraction_pmd_ge_3:.2%}"
            )
        else:
            logger.warning(f"[{sample_id}] PMD Scoring Completed with 0 mapped reads evaluated.")

        return scores_tsv, filtered_bam, filtered_bai, pmd_metrics, meta
