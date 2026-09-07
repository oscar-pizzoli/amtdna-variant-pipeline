"""Alignment filtering engine: MAPQ thresholding, bitwise SAM flag filtering, multi-mapping resolution, and MAPQ profiling."""

from dataclasses import dataclass, field
from enum import Enum
import logging
from pathlib import Path
import shutil
import subprocess
from typing import Any, Dict, List, Optional, Tuple

import pysam

from src.modules.alignment.models import (
    AlignmentFilterPolicy,
    MAPQProfile,
    MultiMappingPolicy,
    MultiMappingStrategy,
)

logger = logging.getLogger("mtDNA_Pipeline.FilterEngine")


class FilterEngineError(Exception):
    """Raised when alignment filtering, verification, or MAPQ profiling fails."""
    pass


class ExclusionReason(str, Enum):
    """Mutually exclusive categories for read filtering."""

    RETAINED = "RETAINED"
    UNMAPPED = "UNMAPPED"
    SECONDARY = "SECONDARY"
    SUPPLEMENTARY = "SUPPLEMENTARY"
    LOW_MAPQ = "LOW_MAPQ"
    MULTIMAPPING_TAG = "MULTIMAPPING_TAG"


@dataclass
class FilterExclusionStats:
    """Exclusive first-failure read counts summing to total input reads."""

    total_input_reads: int = 0
    unmapped_reads_removed: int = 0
    secondary_reads_removed: int = 0
    supplementary_reads_removed: int = 0
    low_mapq_reads_removed: int = 0
    multimapping_tag_reads_removed: int = 0
    retained_reads: int = 0

    @property
    def total_removed(self) -> int:
        return (
            self.unmapped_reads_removed
            + self.secondary_reads_removed
            + self.supplementary_reads_removed
            + self.low_mapq_reads_removed
            + self.multimapping_tag_reads_removed
        )

    @property
    def retention_rate(self) -> float:
        return (self.retained_reads / self.total_input_reads) if self.total_input_reads > 0 else 0.0

    def to_dict(self) -> Dict[str, Any]:
        """Serializes exclusion statistics."""
        return {
            "total_input_reads": self.total_input_reads,
            "unmapped_removed": self.unmapped_reads_removed,
            "secondary_removed": self.secondary_reads_removed,
            "supplementary_removed": self.supplementary_reads_removed,
            "low_mapq_removed": self.low_mapq_reads_removed,
            "multimapping_tag_removed": self.multimapping_tag_reads_removed,
            "total_removed": self.total_removed,
            "retained_reads": self.retained_reads,
            "retention_rate": round(self.retention_rate, 4),
        }


@dataclass
class FilterExecutionResult:
    """Outcome, artifacts, and stratified metrics resulting from BAM filtering."""

    sample_id: str
    filtered_markdup_bam: Optional[Path] = None  # Filtered BAM retaining duplicate flags (0x400)
    filtered_markdup_bai: Optional[Path] = None
    filtered_dedup_bam: Optional[Path] = None    # Filtered BAM with duplicates physically removed
    filtered_dedup_bai: Optional[Path] = None
    primary_downstream_bam: Optional[Path] = None  # Recommended BAM for downstream consensus/calling
    filter_basis: str = "MARKDUP"                # Source dataset used for profiling ("MARKDUP" or "DEDUP")
    raw_mapq_profile: MAPQProfile = field(default_factory=MAPQProfile)
    retained_mapq_profile: MAPQProfile = field(default_factory=MAPQProfile)
    exclusion_stats: FilterExclusionStats = field(default_factory=FilterExclusionStats)
    command_lines: List[str] = field(default_factory=list)
    is_successful: bool = False
    error_message: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Serializes filtering outcome."""
        return {
            "sample_id": self.sample_id,
            "filtered_markdup_bam": str(self.filtered_markdup_bam) if self.filtered_markdup_bam else None,
            "filtered_markdup_bai": str(self.filtered_markdup_bai) if self.filtered_markdup_bai else None,
            "filtered_dedup_bam": str(self.filtered_dedup_bam) if self.filtered_dedup_bam else None,
            "filtered_dedup_bai": str(self.filtered_dedup_bai) if self.filtered_dedup_bai else None,
            "primary_downstream_bam": str(self.primary_downstream_bam) if self.primary_downstream_bam else None,
            "filter_basis": self.filter_basis,
            "raw_mapq_profile": self.raw_mapq_profile.to_dict(),
            "retained_mapq_profile": self.retained_mapq_profile.to_dict(),
            "exclusion_stats": self.exclusion_stats.to_dict(),
            "command_lines": self.command_lines,
            "is_successful": self.is_successful,
            "error_message": self.error_message,
        }


class AlignmentFilterEngine:
    """Applies standardized SAM bitwise flag filters, MAPQ thresholds, and generates stratified quality profiles."""

    @classmethod
    def index_bam(cls, bam_path: Path, threads: int = 2, use_csi: bool = False) -> Path:
        """Generates coordinate index (.bai or .csi) for a target BAM file."""
        if not shutil.which("samtools"):
            raise FilterEngineError("Cannot index BAM: 'samtools' executable not found in PATH.")

        bam = Path(bam_path).resolve()
        if not bam.is_file() or bam.stat().st_size == 0:
            raise FilterEngineError(f"Cannot index missing or empty BAM: {bam.name}")

        cmd = ["samtools", "index", "-@", str(max(1, threads))]
        if use_csi:
            cmd.append("-c")
        cmd.append(str(bam))

        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if res.returncode != 0:
            raise FilterEngineError(
                f"Failed to generate index for {bam.name} (exit code {res.returncode}): {res.stderr.strip()}"
            )

        expected_indices = [
            bam.with_name(bam.name + (".csi" if use_csi else ".bai")),
            bam.with_suffix(".csi" if use_csi else ".bai"),
        ]

        for idx_path in expected_indices:
            if idx_path.is_file() and idx_path.stat().st_size > 0:
                return idx_path

        raise FilterEngineError(f"Index generation completed but output index file was not found for {bam.name}")

    @classmethod
    def verify_bam_readability(cls, bam_path: Path) -> int:
        """Checks BAM readability, header validity, and quick structural integrity via samtools quickcheck.

        Returns:
            Total alignment records count in the BAM file.
        """
        if not bam_path.is_file() or bam_path.stat().st_size == 0:
            raise FilterEngineError(f"BAM file is missing or empty: {bam_path.name}")

        if shutil.which("samtools"):
            res_qc = subprocess.run(
                ["samtools", "quickcheck", "-v", str(bam_path)],
                capture_output=True,
                text=True,
                check=False,
            )
            if res_qc.returncode != 0:
                raise FilterEngineError(
                    f"Integrity check failed ('samtools quickcheck' code {res_qc.returncode}) for {bam_path.name}: {res_qc.stderr.strip()}"
                )

        try:
            with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam_f:
                if not bam_f.header:
                    raise FilterEngineError(f"Header missing in {bam_path.name}")
                return bam_f.count(until_eof=True)
        except FilterEngineError:
            raise
        except Exception as exc:
            raise FilterEngineError(f"Corrupted BAM file structure in {bam_path.name}: {exc}") from exc

    @classmethod
    def verify_filtered_bam_postconditions(
        cls,
        filtered_bam: Path,
        policy: AlignmentFilterPolicy,
    ) -> None:
        """Asserts that output filtered BAM contains zero records violating active filter policies."""
        min_mapq = policy.multimapping.min_mapq
        filter_unmapped = policy.multimapping.filter_unmapped
        filter_sec = policy.multimapping.filter_secondary
        filter_supp = policy.multimapping.filter_supplementary

        try:
            with pysam.AlignmentFile(str(filtered_bam), "rb", check_sq=False) as bam_f:
                for read in bam_f.fetch(until_eof=True):
                    if filter_unmapped and read.is_unmapped:
                        raise FilterEngineError(
                            f"Postcondition violation in {filtered_bam.name}: unmapped read '{read.query_name}' retained."
                        )
                    if filter_sec and read.is_secondary:
                        raise FilterEngineError(
                            f"Postcondition violation in {filtered_bam.name}: secondary read '{read.query_name}' retained."
                        )
                    if filter_supp and read.is_supplementary:
                        raise FilterEngineError(
                            f"Postcondition violation in {filtered_bam.name}: supplementary read '{read.query_name}' retained."
                        )
                    if not read.is_unmapped and read.mapping_quality < min_mapq:
                        raise FilterEngineError(
                            f"Postcondition violation in {filtered_bam.name}: read '{read.query_name}' has MAPQ={read.mapping_quality} < {min_mapq}."
                        )
        except FilterEngineError:
            raise
        except Exception as exc:
            raise FilterEngineError(f"Error during postcondition verification on {filtered_bam.name}: {exc}") from exc

    @classmethod
    def _has_multimapping_tag(cls, read: pysam.AlignedSegment) -> bool:
        """Inspects optional mapper-specific multi-mapping tags (e.g., BWA XT:A:R, XA:Z, or AS/XS delta)."""
        # BWA tag XT:A:R (Repetitive / non-unique)
        if read.has_tag("XT") and read.get_tag("XT") == "R":
            return True

        # BWA tag XA:Z (Alternative hits)
        if read.has_tag("XA"):
            return True

        # BWA tag X1:i (Number of suboptimal hits)
        if read.has_tag("X1") and read.get_tag("X1") > 0:
            return True

        # Modern aligners AS/XS score comparison
        if read.has_tag("AS") and read.has_tag("XS"):
            as_score = read.get_tag("AS")
            xs_score = read.get_tag("XS")
            if as_score == xs_score:
                return True

        return False

    @classmethod
    def compute_mapq_profiles_and_exclusions(
        cls,
        input_bam: Path,
        policy: AlignmentFilterPolicy,
    ) -> Tuple[MAPQProfile, MAPQProfile, FilterExclusionStats]:
        """Scans input BAM to compute raw profile, retained profile, and exclusive first-failure statistics.

        Returns:
            Tuple of (raw_mapq_profile, retained_mapq_profile, exclusion_stats).
        """
        min_mapq = policy.multimapping.min_mapq
        filter_unmapped = policy.multimapping.filter_unmapped
        filter_secondary = policy.multimapping.filter_secondary
        filter_supp = policy.multimapping.filter_supplementary
        inspect_tags = policy.multimapping.inspect_mapper_tags

        stats = FilterExclusionStats()

        # Raw mapped tracking
        raw_counts = [0] * 5  # [q0, q1_9, q10_19, q20_29, q30_plus]
        raw_mapq_sum = 0
        raw_mapped_total = 0

        # Retained mapped tracking
        ret_counts = [0] * 5
        ret_mapq_sum = 0
        ret_mapped_total = 0

        def update_bins(bin_list: List[int], q_val: int) -> None:
            if q_val == 0:
                bin_list[0] += 1
            elif 1 <= q_val <= 9:
                bin_list[1] += 1
            elif 10 <= q_val <= 19:
                bin_list[2] += 1
            elif 20 <= q_val <= 29:
                bin_list[3] += 1
            else:
                bin_list[4] += 1

        try:
            with pysam.AlignmentFile(str(input_bam), "rb", check_sq=False) as bam_f:
                for read in bam_f.fetch(until_eof=True):
                    stats.total_input_reads += 1

                    # Track raw mapped alignments before exclusion hierarchy
                    if not read.is_unmapped:
                        mapq = read.mapping_quality
                        raw_mapped_total += 1
                        raw_mapq_sum += mapq
                        update_bins(raw_counts, mapq)

                    # --- Exclusive First-Failure Decision Tree ---
                    # 1. Unmapped Check
                    if read.is_unmapped:
                        if filter_unmapped:
                            stats.unmapped_reads_removed += 1
                            continue
                        # If unmapped reads are explicitly not filtered, record in retained
                        stats.retained_reads += 1
                        continue

                    # 2. Secondary Alignment Check (0x100)
                    if filter_secondary and read.is_secondary:
                        stats.secondary_reads_removed += 1
                        continue

                    # 3. Supplementary Alignment Check (0x800)
                    if filter_supp and read.is_supplementary:
                        stats.supplementary_reads_removed += 1
                        continue

                    # 4. MAPQ Threshold Cutoff
                    if read.mapping_quality < min_mapq:
                        stats.low_mapq_reads_removed += 1
                        continue

                    # 5. Optional Mapper-Specific Multi-Mapping Tag Check
                    if inspect_tags and cls._has_multimapping_tag(read):
                        stats.multimapping_tag_reads_removed += 1
                        continue

                    # Read passed all filters -> RETAINED
                    stats.retained_reads += 1
                    ret_mapped_total += 1
                    ret_mapq_sum += read.mapping_quality
                    update_bins(ret_counts, read.mapping_quality)

        except Exception as exc:
            raise FilterEngineError(f"Failed to scan BAM records in {input_bam.name}: {exc}") from exc

        raw_mean = (raw_mapq_sum / raw_mapped_total) if raw_mapped_total > 0 else 0.0
        ret_mean = (ret_mapq_sum / ret_mapped_total) if ret_mapped_total > 0 else 0.0

        raw_profile = MAPQProfile(
            count_mapq_0=raw_counts[0],
            count_mapq_1_9=raw_counts[1],
            count_mapq_10_19=raw_counts[2],
            count_mapq_20_29=raw_counts[3],
            count_mapq_30_plus=raw_counts[4],
            mean_mapq=raw_mean,
        )

        retained_profile = MAPQProfile(
            count_mapq_0=ret_counts[0],
            count_mapq_1_9=ret_counts[1],
            count_mapq_10_19=ret_counts[2],
            count_mapq_20_29=ret_counts[3],
            count_mapq_30_plus=ret_counts[4],
            mean_mapq=ret_mean,
        )

        return raw_profile, retained_profile, stats

    @classmethod
    def filter_bam_file(
        cls,
        input_bam: Path,
        output_bam: Path,
        policy: AlignmentFilterPolicy,
        threads: int,
        log_dir: Path,
    ) -> str:
        """Applies SAM bitwise flag and MAPQ filters using samtools view."""
        output_bam.parent.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)

        log_err = log_dir / f"{output_bam.stem}.samtools_filter.err"

        # Calculate bitwise exclusion mask (-F)
        exclude_flag = 0
        if policy.multimapping.filter_unmapped:
            exclude_flag |= 0x4        # Segment unmapped
        if policy.multimapping.filter_secondary:
            exclude_flag |= 0x100      # Secondary alignment
        if policy.multimapping.filter_supplementary:
            exclude_flag |= 0x800      # Supplementary alignment

        cmd = [
            "samtools", "view",
            "-@", str(max(1, threads)),
            "-b",
            "-q", str(policy.multimapping.min_mapq),
        ]

        if exclude_flag > 0:
            cmd.extend(["-F", hex(exclude_flag)])

        cmd.extend([
            "-o", str(output_bam),
            str(input_bam),
        ])

        cmd_str = " ".join(cmd)
        logger.debug(f"Executing BAM filter: {cmd_str}")

        with open(log_err, "w", encoding="utf-8") as err_f:
            res = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=err_f, check=False)

        if res.returncode != 0 or not output_bam.is_file() or output_bam.stat().st_size == 0:
            raise FilterEngineError(
                f"'samtools view' filtering failed (exit code {res.returncode}) for {output_bam.name}. Check log: {log_err}"
            )

        return cmd_str

    @classmethod
    def execute_filtering(
        cls,
        sample_id: str,
        input_markdup_bam: Optional[Path],
        input_dedup_bam: Optional[Path],
        output_dir: Path,
        log_dir: Path,
        policy: Optional[AlignmentFilterPolicy] = None,
        threads: int = 4,
        fast_mode: bool = False,
        use_csi: bool = False,
        cleanup_on_failure: bool = True,
    ) -> FilterExecutionResult:
        """Applies configured quality filters to marked and/or deduplicated BAM files.

        Args:
            sample_id: Unique sample identifier.
            input_markdup_bam: Path to BAM with duplicates flagged (0x400).
            input_dedup_bam: Path to physically deduplicated BAM.
            output_dir: Target output directory for filtered BAM files.
            log_dir: Directory for process logs.
            policy: Filtering rules, thresholds, and multi-mapping policy.
            threads: Number of CPU threads allocated for filtering and compression.
            fast_mode: If True, skips Python-based granular MAPQ profiling to optimize speed.
            use_csi: If True, generates .csi indices instead of standard .bai.
            cleanup_on_failure: If True, unlinks partial files upon failure.

        Returns:
            FilterExecutionResult containing validated paths, indices, dual MAPQ profiles, and exclusion stats.
        """
        if not input_markdup_bam and not input_dedup_bam:
            raise FilterEngineError(f"[{sample_id}] No input BAM provided for filtering.")

        # Determine primary basis for profiling
        primary_input = input_markdup_bam or input_dedup_bam
        if not primary_input or not primary_input.is_file():
            raise FilterEngineError(f"[{sample_id}] Input BAM file does not exist: {primary_input}")

        filter_basis = "MARKDUP" if input_markdup_bam else "DEDUP"
        filter_policy = policy or AlignmentFilterPolicy()

        output_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)

        filtered_markdup_bam = output_dir / f"{sample_id}.filtered.markdup.bam" if input_markdup_bam else None
        filtered_markdup_bai: Optional[Path] = None
        filtered_dedup_bam = output_dir / f"{sample_id}.filtered.dedup.bam" if input_dedup_bam else None
        filtered_dedup_bai: Optional[Path] = None

        executed_cmds: List[str] = []

        logger.info(
            f"[{sample_id}] Starting alignment filtering on basis '{filter_basis}' (MAPQ >= {filter_policy.multimapping.min_mapq}, "
            f"Filter Unmapped: {filter_policy.multimapping.filter_unmapped}, Secondary: {filter_policy.multimapping.filter_secondary})..."
        )

        try:
            # 1. Compute granular MAPQ profiles and exclusive first-failure statistics
            if not fast_mode:
                raw_profile, retained_profile, exclusion_stats = cls.compute_mapq_profiles_and_exclusions(
                    input_bam=primary_input,
                    policy=filter_policy,
                )
                logger.info(
                    f"[{sample_id}] Profiling complete: Raw Mean MAPQ={raw_profile.mean_mapq:.2f} | "
                    f"Retained Mean MAPQ={retained_profile.mean_mapq:.2f} | "
                    f"Retained {exclusion_stats.retained_reads:,} / {exclusion_stats.total_input_reads:,} reads "
                    f"({exclusion_stats.retention_rate * 100:.2f}%)"
                )
            else:
                raw_profile = MAPQProfile()
                retained_profile = MAPQProfile()
                total_in = cls.verify_bam_readability(primary_input)
                exclusion_stats = FilterExclusionStats(total_input_reads=total_in)

            # 2. Filter marked duplicate BAM (.filtered.markdup.bam)
            if input_markdup_bam and filtered_markdup_bam:
                cmd_mark = cls.filter_bam_file(
                    input_bam=input_markdup_bam,
                    output_bam=filtered_markdup_bam,
                    policy=filter_policy,
                    threads=threads,
                    log_dir=log_dir,
                )
                executed_cmds.append(cmd_mark)

                cls.verify_bam_readability(filtered_markdup_bam)
                cls.verify_filtered_bam_postconditions(filtered_markdup_bam, filter_policy)
                filtered_markdup_bai = cls.index_bam(filtered_markdup_bam, threads=threads, use_csi=use_csi)
                executed_cmds.append(f"samtools index -@ {threads} {filtered_markdup_bam}")

            # 3. Filter physically deduplicated BAM (.filtered.dedup.bam)
            if input_dedup_bam and filtered_dedup_bam:
                cmd_dedup = cls.filter_bam_file(
                    input_bam=input_dedup_bam,
                    output_bam=filtered_dedup_bam,
                    policy=filter_policy,
                    threads=threads,
                    log_dir=log_dir,
                )
                executed_cmds.append(cmd_dedup)

                cls.verify_bam_readability(filtered_dedup_bam)
                cls.verify_filtered_bam_postconditions(filtered_dedup_bam, filter_policy)
                filtered_dedup_bai = cls.index_bam(filtered_dedup_bam, threads=threads, use_csi=use_csi)
                executed_cmds.append(f"samtools index -@ {threads} {filtered_dedup_bam}")

            # 4. Explicitly designate the primary downstream BAM for Block 4 and 5
            primary_downstream_bam = filtered_dedup_bam if filtered_dedup_bam else filtered_markdup_bam

            logger.info(
                f"[{sample_id}] Filtering completed successfully. Primary downstream BAM: "
                f"{primary_downstream_bam.name if primary_downstream_bam else 'None'}"
            )

            return FilterExecutionResult(
                sample_id=sample_id,
                filtered_markdup_bam=filtered_markdup_bam,
                filtered_markdup_bai=filtered_markdup_bai,
                filtered_dedup_bam=filtered_dedup_bam,
                filtered_dedup_bai=filtered_dedup_bai,
                primary_downstream_bam=primary_downstream_bam,
                filter_basis=filter_basis,
                raw_mapq_profile=raw_profile,
                retained_mapq_profile=retained_profile,
                exclusion_stats=exclusion_stats,
                command_lines=executed_cmds,
                is_successful=True,
            )

        except Exception as exc:
            logger.error(f"[{sample_id}] Filtering engine failed: {exc}")

            if cleanup_on_failure:
                for partial_file in (filtered_markdup_bam, filtered_markdup_bai, filtered_dedup_bam, filtered_dedup_bai):
                    if partial_file and partial_file.is_file():
                        try:
                            partial_file.unlink()
                        except OSError:
                            pass

            return FilterExecutionResult(
                sample_id=sample_id,
                filtered_markdup_bam=None,
                filtered_markdup_bai=None,
                filtered_dedup_bam=None,
                filtered_dedup_bai=None,
                primary_downstream_bam=None,
                filter_basis=filter_basis,
                raw_mapq_profile=MAPQProfile(),
                retained_mapq_profile=MAPQProfile(),
                exclusion_stats=FilterExclusionStats(),
                command_lines=executed_cmds,
                is_successful=False,
                error_message=str(exc),
            )
