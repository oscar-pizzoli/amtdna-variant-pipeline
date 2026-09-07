"""Quality control and statistical metrics collection engine: samtools flagstat parsing, depth profiling, and breadth calculation."""

from dataclasses import dataclass, field
import logging
import math
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any, Dict, List, Optional, Tuple, Union

import pysam

from src.modules.alignment.models import (
    AlignmentMetrics,
    CoverageMetrics,
    DedupMetrics,
    MAPQProfile,
    ProvenanceMetrics,
    ReadGroupSource,
    ReferenceInfo,
)
from src.modules.alignment.bwa_aligner import AlignmentBranchResults
from src.modules.alignment.dedup_engine import DedupExecutionResult
from src.modules.alignment.filter_engine import FilterExecutionResult

logger = logging.getLogger("mtDNA_Pipeline.MetricsCollector")


class MetricsCollectorError(Exception):
    """Raised when QC metrics extraction, depth calculation, or flagstat validation fails."""
    pass


@dataclass
class FlagstatMetrics:
    """Structured container holding parsed samtools flagstat counts and validation invariants."""

    total_reads: int = 0
    primary_reads: Optional[int] = None
    secondary_reads: int = 0
    supplementary_reads: int = 0
    duplicates: int = 0
    mapped_reads: int = 0
    mapped_pct: float = 0.0  # Mapping rate percentage [0.0 - 100.0%]
    paired_reads: int = 0
    read1: int = 0
    read2: int = 0
    properly_paired: int = 0
    with_itself_and_mate_mapped: int = 0
    singletons: int = 0

    def to_dict(self) -> Dict[str, Any]:
        """Serializes flagstat metrics."""
        return {
            "total_reads": self.total_reads,
            "primary_reads": self.primary_reads,
            "secondary_reads": self.secondary_reads,
            "supplementary_reads": self.supplementary_reads,
            "duplicates": self.duplicates,
            "mapped_reads": self.mapped_reads,
            "mapped_pct": round(self.mapped_pct, 2),
            "paired_reads": self.paired_reads,
            "read1": self.read1,
            "read2": self.read2,
            "properly_paired": self.properly_paired,
            "with_itself_and_mate_mapped": self.with_itself_and_mate_mapped,
            "singletons": self.singletons,
        }


@dataclass
class MetricsCollectionResult:
    """Top-level container grouping all intermediate flagstat reports, depth profiles, and consolidated metrics."""

    sample_id: str
    raw_flagstat: Optional[FlagstatMetrics] = None
    filtered_flagstat: Optional[FlagstatMetrics] = None
    dedup_flagstat: Optional[FlagstatMetrics] = None
    coverage: CoverageMetrics = field(default_factory=CoverageMetrics)
    provenance: ProvenanceMetrics = field(default_factory=ProvenanceMetrics)
    alignment_metrics: AlignmentMetrics = field(default_factory=AlignmentMetrics)
    depth_tsv: Optional[Path] = None
    command_lines: List[str] = field(default_factory=list)
    is_successful: bool = False
    error_message: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Serializes complete metrics collection result."""
        return {
            "sample_id": self.sample_id,
            "raw_flagstat": self.raw_flagstat.to_dict() if self.raw_flagstat else None,
            "filtered_flagstat": self.filtered_flagstat.to_dict() if self.filtered_flagstat else None,
            "dedup_flagstat": self.dedup_flagstat.to_dict() if self.dedup_flagstat else None,
            "coverage": self.coverage.to_dict(),
            "provenance": self.provenance.to_dict(),
            "alignment_metrics": self.alignment_metrics.to_dict(),
            "depth_tsv": str(self.depth_tsv) if self.depth_tsv else None,
            "command_lines": self.command_lines,
            "is_successful": self.is_successful,
            "error_message": self.error_message,
        }


class MetricsCollector:
    """Extracts alignment statistics, depth percentiles, and coverage breadth across pipeline stages."""

    @classmethod
    def _validate_flagstat_invariants(cls, metrics: FlagstatMetrics, source_name: str) -> None:
        """Validates mathematical consistency of parsed flagstat counts."""
        if metrics.mapped_reads > metrics.total_reads:
            raise MetricsCollectorError(
                f"Flagstat invariant violation in {source_name}: mapped reads ({metrics.mapped_reads:,}) > total reads ({metrics.total_reads:,})"
            )
        if metrics.primary_reads is not None and metrics.primary_reads > metrics.total_reads:
            raise MetricsCollectorError(
                f"Flagstat invariant violation in {source_name}: primary reads ({metrics.primary_reads:,}) > total reads ({metrics.total_reads:,})"
            )
        if metrics.duplicates > metrics.total_reads:
            raise MetricsCollectorError(
                f"Flagstat invariant violation in {source_name}: duplicates ({metrics.duplicates:,}) > total reads ({metrics.total_reads:,})"
            )
        if (metrics.secondary_reads + metrics.supplementary_reads) > metrics.total_reads:
            raise MetricsCollectorError(
                f"Flagstat invariant violation in {source_name}: secondary ({metrics.secondary_reads:,}) + supplementary ({metrics.supplementary_reads:,}) > total reads"
            )

    @classmethod
    def parse_flagstat_output(cls, flagstat_content_or_file: Union[Path, str]) -> FlagstatMetrics:
        """Parses samtools flagstat standard text output into a strongly typed FlagstatMetrics object."""
        text = ""
        source_name = "raw_string"
        if isinstance(flagstat_content_or_file, Path):
            source_name = flagstat_content_or_file.name
            if not flagstat_content_or_file.is_file():
                raise MetricsCollectorError(f"Flagstat file not found: {flagstat_content_or_file}")
            try:
                text = flagstat_content_or_file.read_text(encoding="utf-8")
            except OSError as exc:
                raise MetricsCollectorError(f"Could not read flagstat file {flagstat_content_or_file}: {exc}") from exc
        elif isinstance(flagstat_content_or_file, str):
            text = flagstat_content_or_file

        if not text.strip():
            raise MetricsCollectorError(f"Empty flagstat content provided from {source_name}")

        def extract_required(pattern: str, field_name: str) -> int:
            match = re.search(pattern, text, re.IGNORECASE)
            if not match:
                raise MetricsCollectorError(f"Required flagstat field '{field_name}' not found in {source_name}")
            return int(match.group(1))

        def extract_optional(pattern: str) -> Optional[int]:
            match = re.search(pattern, text, re.IGNORECASE)
            return int(match.group(1)) if match else None

        # Pattern matches: "<passed> + <failed> <description>"
        total_reads = extract_required(r"(\d+)\s*\+\s*\d+\s+in total", "in total")
        primary_reads = extract_optional(r"(\d+)\s*\+\s*\d+\s+primary")
        secondary_reads = extract_optional(r"(\d+)\s*\+\s*\d+\s+secondary") or 0
        supplementary_reads = extract_optional(r"(\d+)\s*\+\s*\d+\s+supplementary") or 0
        duplicates = extract_optional(r"(\d+)\s*\+\s*\d+\s+duplicates") or 0
        mapped_reads = extract_required(r"(\d+)\s*\+\s*\d+\s+mapped\s*\(", "mapped")
        paired_reads = extract_optional(r"(\d+)\s*\+\s*\d+\s+paired in sequencing") or 0
        read1 = extract_optional(r"(\d+)\s*\+\s*\d+\s+read1") or 0
        read2 = extract_optional(r"(\d+)\s*\+\s*\d+\s+read2") or 0
        properly_paired = extract_optional(r"(\d+)\s*\+\s*\d+\s+properly paired") or 0
        itself_and_mate = extract_optional(r"(\d+)\s*\+\s*\d+\s+with itself and mate mapped") or 0
        singletons = extract_optional(r"(\d+)\s*\+\s*\d+\s+singletons") or 0

        # If primary_reads was not explicitly reported, calculate it from total - secondary - supplementary
        if primary_reads is None and total_reads > 0:
            primary_reads = total_reads - secondary_reads - supplementary_reads

        mapped_pct = ((mapped_reads / total_reads) * 100.0) if total_reads > 0 else 0.0

        metrics = FlagstatMetrics(
            total_reads=total_reads,
            primary_reads=primary_reads,
            secondary_reads=secondary_reads,
            supplementary_reads=supplementary_reads,
            duplicates=duplicates,
            mapped_reads=mapped_reads,
            mapped_pct=mapped_pct,
            paired_reads=paired_reads,
            read1=read1,
            read2=read2,
            properly_paired=properly_paired,
            with_itself_and_mate_mapped=itself_and_mate,
            singletons=singletons,
        )

        cls._validate_flagstat_invariants(metrics, source_name)
        return metrics

    @classmethod
    def run_flagstat(cls, bam_path: Path, output_txt: Path, threads: int = 2) -> Tuple[FlagstatMetrics, str]:
        """Executes 'samtools flagstat' and saves the report file."""
        if not shutil.which("samtools"):
            raise MetricsCollectorError("Cannot run flagstat: 'samtools' executable not found in PATH.")

        bam = Path(bam_path).resolve()
        if not bam.is_file() or bam.stat().st_size == 0:
            raise MetricsCollectorError(f"Cannot run flagstat on missing or empty BAM: {bam.name}")

        output_txt.parent.mkdir(parents=True, exist_ok=True)
        cmd = ["samtools", "flagstat", "-@", str(max(1, threads)), str(bam)]
        cmd_str = " ".join(cmd)

        try:
            with open(output_txt, "w", encoding="utf-8") as out_f:
                res = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True, check=False)

            if res.returncode != 0:
                raise MetricsCollectorError(
                    f"'samtools flagstat' failed on {bam.name} (exit code {res.returncode}): {res.stderr.strip()}"
                )

            metrics = cls.parse_flagstat_output(output_txt)
            return metrics, cmd_str

        except OSError as exc:
            raise MetricsCollectorError(f"Failed to execute flagstat on {bam.name}: {exc}") from exc

    @classmethod
    def compute_depth_distribution(
        cls,
        bam_path: Path,
        reference_info: ReferenceInfo,
        output_depth_tsv: Optional[Path] = None,
        min_base_quality: int = 0,
        min_mapping_quality: int = 0,
    ) -> Tuple[List[int], float, float, int, float, float]:
        """Calculates per-base mitochondrial depth distribution across all positions using 'samtools depth -aa -r contig'.

        Args:
            bam_path: Path to the coordinate-sorted and indexed BAM file.
            reference_info: Validated reference genome metadata containing contig_name and length_bp.
            output_depth_tsv: Optional path to write position-by-position depth values.
            min_base_quality: Optional base quality cutoff (-q parameter).
            min_mapping_quality: Optional mapping quality cutoff (-Q parameter).

        Returns:
            Tuple of (depth_array, mean_depth, median_depth, max_depth, p95_depth, p99_depth).
        """
        if not shutil.which("samtools"):
            raise MetricsCollectorError("Cannot compute depth: 'samtools' executable not found in PATH.")

        bam = Path(bam_path).resolve()
        if not bam.is_file() or bam.stat().st_size == 0:
            raise MetricsCollectorError(f"Cannot compute depth on missing BAM: {bam.name}")

        ref_length = reference_info.length_bp
        contig_name = reference_info.contig_name or reference_info.sequence_id

        if not ref_length or ref_length <= 0:
            raise MetricsCollectorError(f"Invalid reference length ({ref_length}) provided in ReferenceInfo.")
        if not contig_name:
            raise MetricsCollectorError("Reference contig name is undefined in ReferenceInfo.")

        cmd = [
            "samtools", "depth",
            "-aa",
            "-r", str(contig_name),
        ]
        if min_base_quality > 0:
            cmd.extend(["-q", str(min_base_quality)])
        if min_mapping_quality > 0:
            cmd.extend(["-Q", str(min_mapping_quality)])
        cmd.append(str(bam))

        depth_array = [0] * ref_length
        positions_recorded = 0
        last_pos_seen = 0

        try:
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            tsv_file = open(output_depth_tsv, "w", encoding="utf-8") if output_depth_tsv else None

            if tsv_file:
                tsv_file.write("chrom\tpos\tdepth\n")

            if p.stdout:
                for line in p.stdout:
                    if tsv_file:
                        tsv_file.write(line)
                    parts = line.strip().split("\t")
                    if len(parts) >= 3:
                        chrom = parts[0]
                        pos = int(parts[1])
                        depth_val = int(parts[2])

                        if chrom != contig_name:
                            raise MetricsCollectorError(
                                f"Unexpected contig '{chrom}' found in depth stream (expected '{contig_name}')."
                            )
                        if pos < 1 or pos > ref_length:
                            raise MetricsCollectorError(
                                f"Depth position {pos} is out of bounds for reference length {ref_length:,} bp."
                            )
                        # Assert strict sequential coordinate progression
                        if pos != last_pos_seen + 1:
                            raise MetricsCollectorError(
                                f"Discontinuous or duplicate depth position encountered in {bam.name}: "
                                f"expected position {last_pos_seen + 1}, received position {pos}."
                            )

                        depth_array[pos - 1] = depth_val
                        positions_recorded += 1
                        last_pos_seen = pos

            if tsv_file:
                tsv_file.close()

            p.communicate()
            if p.returncode != 0:
                raise MetricsCollectorError(f"'samtools depth' failed on {bam.name} (exit code {p.returncode})")

        except Exception as exc:
            if output_depth_tsv and output_depth_tsv.is_file():
                output_depth_tsv.unlink()
            raise MetricsCollectorError(f"Error computing depth distribution on {bam.name}: {exc}") from exc

        # Assert full reference coverage streaming
        if positions_recorded != ref_length:
            raise MetricsCollectorError(
                f"Incomplete depth output for {bam.name}: recorded {positions_recorded:,} positions, "
                f"expected exactly {ref_length:,} bp for contig '{contig_name}'."
            )

        # Calculate exact percentiles and statistics
        sorted_depths = sorted(depth_array)
        total_depth_sum = sum(sorted_depths)
        mean_depth = total_depth_sum / ref_length
        max_depth = sorted_depths[-1]

        def get_percentile(p_rank: float) -> float:
            k = (ref_length - 1) * p_rank
            f = math.floor(k)
            c = math.ceil(k)
            if f == c:
                return float(sorted_depths[int(k)])
            return float(sorted_depths[int(f)] * (c - k) + sorted_depths[int(c)] * (k - f))

        median_depth = get_percentile(0.50)
        p95_depth = get_percentile(0.95)
        p99_depth = get_percentile(0.99)

        return depth_array, mean_depth, median_depth, max_depth, p95_depth, p99_depth

    @classmethod
    def calculate_coverage_metrics(
        cls,
        final_bam: Path,
        reference_info: ReferenceInfo,
        baseline_bam: Optional[Path] = None,
        output_depth_tsv: Optional[Path] = None,
        min_base_quality: int = 0,
        min_mapping_quality: int = 0,
    ) -> CoverageMetrics:
        """Computes comprehensive depth percentiles and breadth of coverage percentages for the target reference.

        Args:
            final_bam: Path to the primary downstream BAM file (e.g. filtered.dedup.bam).
            reference_info: Validated ReferenceInfo instance.
            baseline_bam: Optional path to pre-deduplication baseline BAM (e.g. prepared.sorted.bam).
            output_depth_tsv: Optional path to save tabular depth file.
            min_base_quality: Optional base quality threshold.
            min_mapping_quality: Optional mapping quality threshold.

        Returns:
            Fully populated CoverageMetrics container.
        """
        ref_length = reference_info.length_bp
        if not ref_length or ref_length <= 0:
            raise MetricsCollectorError(f"ReferenceInfo length is invalid: {ref_length}")

        # 1. Compute baseline mean depth before deduplication if baseline BAM is provided
        mean_pre: Optional[float] = None
        if baseline_bam and baseline_bam.is_file():
            _, mean_val, _, _, _, _ = cls.compute_depth_distribution(
                bam_path=baseline_bam,
                reference_info=reference_info,
                output_depth_tsv=None,
                min_base_quality=min_base_quality,
                min_mapping_quality=min_mapping_quality,
            )
            mean_pre = mean_val

        # 2. Compute final depth distribution and percentiles on the primary downstream BAM
        depths, mean_post, median_post, max_post, p95_post, p99_post = cls.compute_depth_distribution(
            bam_path=final_bam,
            reference_info=reference_info,
            output_depth_tsv=output_depth_tsv,
            min_base_quality=min_base_quality,
            min_mapping_quality=min_mapping_quality,
        )

        # Count covered bases at standard quality thresholds
        c_1x = sum(1 for d in depths if d >= 1)
        c_3x = sum(1 for d in depths if d >= 3)
        c_5x = sum(1 for d in depths if d >= 5)
        c_10x = sum(1 for d in depths if d >= 10)
        c_30x = sum(1 for d in depths if d >= 30)

        total_pos = float(ref_length)

        return CoverageMetrics(
            mean_depth_pre_dedup=mean_pre,
            mean_depth_post_dedup=mean_post,
            median_depth_post_dedup=median_post,
            max_depth_post_dedup=max_post,
            p95_depth_post_dedup=p95_post,
            p99_depth_post_dedup=p99_post,
            covered_bases_1x=c_1x,
            breadth_1x_pct=(c_1x / total_pos) * 100.0,
            breadth_3x_pct=(c_3x / total_pos) * 100.0,
            breadth_5x_pct=(c_5x / total_pos) * 100.0,
            breadth_10x_pct=(c_10x / total_pos) * 100.0,
            breadth_30x_pct=(c_30x / total_pos) * 100.0,
        )

    @classmethod
    def extract_branch_mapped_counts(
        cls,
        branch_results: AlignmentBranchResults,
    ) -> Dict[ReadGroupSource, int]:
        """Inspects branch BAMs using pysam to count mapped records per branch."""
        mapped_counts: Dict[ReadGroupSource, int] = {}
        for source, branch in branch_results.branches.items():
            if branch.is_successful and branch.output_bam.is_file():
                try:
                    with pysam.AlignmentFile(str(branch.output_bam), "rb", check_sq=False) as bam_f:
                        mapped = sum(1 for r in bam_f.fetch(until_eof=True) if not r.is_unmapped)
                        # For paired branch, convert mapped reads to pair count
                        if source == ReadGroupSource.PAIRED:
                            mapped_counts[source] = mapped // 2
                        else:
                            mapped_counts[source] = mapped
                except Exception as exc:
                    logger.warning(f"Could not count mapped reads in branch {source.value}: {exc}")
                    mapped_counts[source] = 0
            else:
                mapped_counts[source] = 0
        return mapped_counts

    @classmethod
    def validate_alignment_metrics_invariants(cls, metrics: AlignmentMetrics) -> None:
        """Asserts mathematical consistency across all final alignment metrics."""
        cov = metrics.coverage

        # 1. Breadth monotonic decrease invariant
        if not (cov.breadth_1x_pct >= cov.breadth_3x_pct >= cov.breadth_5x_pct >= cov.breadth_10x_pct >= cov.breadth_30x_pct):
            raise MetricsCollectorError(
                f"Breadth invariant violation: [1x={cov.breadth_1x_pct:.2f}%, 3x={cov.breadth_3x_pct:.2f}%, "
                f"5x={cov.breadth_5x_pct:.2f}%, 10x={cov.breadth_10x_pct:.2f}%, 30x={cov.breadth_30x_pct:.2f}%]"
            )

        # 2. Depth percentiles monotonic order invariant
        if not (cov.max_depth_post_dedup >= cov.p99_depth_post_dedup >= cov.p95_depth_post_dedup >= cov.median_depth_post_dedup >= 0):
            raise MetricsCollectorError(
                f"Depth percentiles invariant violation: [Max={cov.max_depth_post_dedup}, P99={cov.p99_depth_post_dedup}, "
                f"P95={cov.p95_depth_post_dedup}, Median={cov.median_depth_post_dedup}]"
            )

        # 3. Read retention invariants
        if metrics.retained_filtered_reads > metrics.total_input_reads:
            raise MetricsCollectorError(
                f"Retention invariant violation: retained filtered reads ({metrics.retained_filtered_reads:,}) > total input reads ({metrics.total_input_reads:,})"
            )
        if metrics.mapped_reads_raw is not None and metrics.mapped_reads_raw > metrics.total_input_reads:
            raise MetricsCollectorError(
                f"Raw mapped reads invariant violation: mapped raw ({metrics.mapped_reads_raw:,}) > total input ({metrics.total_input_reads:,})"
            )

    @classmethod
    def assemble_alignment_metrics(
        cls,
        branch_results: AlignmentBranchResults,
        filter_result: FilterExecutionResult,
        dedup_result: DedupExecutionResult,
        coverage_metrics: CoverageMetrics,
        raw_flagstat: Optional[FlagstatMetrics] = None,
        scan_branch_bams: bool = False,
    ) -> AlignmentMetrics:
        """Consolidates metrics from all alignment stages into the domain AlignmentMetrics schema."""
        # 1. Build provenance tracking
        prov = ProvenanceMetrics()
        mapped_by_branch = cls.extract_branch_mapped_counts(branch_results) if scan_branch_bams else {}

        for source, branch in branch_results.branches.items():
            if source == ReadGroupSource.MERGED:
                prov.input_merged_fragments = branch.input_units_count
                prov.mapped_merged_fragments = mapped_by_branch.get(source) if scan_branch_bams else None
            elif source == ReadGroupSource.MERGED_TRUNCATED:
                prov.input_merged_truncated_fragments = branch.input_units_count
                prov.mapped_merged_truncated_fragments = mapped_by_branch.get(source) if scan_branch_bams else None
            elif source == ReadGroupSource.SINGLETON:
                prov.input_singleton_reads = branch.input_units_count
                prov.mapped_singleton_reads = mapped_by_branch.get(source) if scan_branch_bams else None
            elif source == ReadGroupSource.PAIRED:
                prov.input_paired_pairs = branch.input_units_count
                prov.mapped_paired_pairs = mapped_by_branch.get(source) if scan_branch_bams else None

        # 2. Raw metrics (pre-filter)
        total_in = raw_flagstat.total_reads if raw_flagstat else filter_result.exclusion_stats.total_input_reads
        mapped_raw = raw_flagstat.mapped_reads if raw_flagstat else None
        unmapped_raw = (total_in - mapped_raw) if (raw_flagstat and mapped_raw is not None) else None
        sec_raw = raw_flagstat.secondary_reads if raw_flagstat else None
        supp_raw = raw_flagstat.supplementary_reads if raw_flagstat else None
        pre_dup = raw_flagstat.duplicates if raw_flagstat else None
        raw_map_pct = raw_flagstat.mapped_pct if raw_flagstat else None

        # 3. Filtering metrics
        ex = filter_result.exclusion_stats

        metrics = AlignmentMetrics(
            total_input_reads=total_in,
            mapped_reads_raw=mapped_raw,
            unmapped_reads_raw=unmapped_raw,
            secondary_reads_raw=sec_raw,
            supplementary_reads_raw=supp_raw,
            premarked_duplicates_raw=pre_dup,
            raw_mapping_rate_pct=raw_map_pct,
            filtered_unmapped_reads=ex.unmapped_reads_removed,
            filtered_secondary_reads=ex.secondary_reads_removed,
            filtered_supplementary_reads=ex.supplementary_reads_removed,
            filtered_low_mapq_reads=ex.low_mapq_reads_removed,
            retained_filtered_reads=ex.retained_reads,
            provenance=prov,
            mapq_profile=filter_result.retained_mapq_profile,
            dedup_metrics=dedup_result.metrics or DedupMetrics(),
            coverage=coverage_metrics,
        )

        cls.validate_alignment_metrics_invariants(metrics)
        return metrics

    @classmethod
    def collect_all_metrics(
        cls,
        sample_id: str,
        reference_info: ReferenceInfo,
        branch_results: AlignmentBranchResults,
        dedup_result: DedupExecutionResult,
        filter_result: FilterExecutionResult,
        output_dir: Path,
        log_dir: Path,
        prepared_bam: Optional[Path] = None,
        threads: int = 2,
        scan_branch_bams: bool = False,
        min_base_quality: int = 0,
        min_mapping_quality: int = 0,
    ) -> MetricsCollectionResult:
        """Orchestrates comprehensive flagstat extraction across stages, depth calculation, and metric assembly.

        Args:
            sample_id: Sample identifier.
            reference_info: Validated reference genome metadata.
            branch_results: Results from BWA multi-branch alignment.
            dedup_result: Results from samtools markdup deduplication.
            filter_result: Results from alignment filtering.
            output_dir: Directory where reports and depth TSVs are written.
            log_dir: Directory where execution logs are stored.
            prepared_bam: Optional path to pre-deduplication baseline BAM (prepared.sorted.bam).
            threads: Number of CPU threads for flagstat execution.
            scan_branch_bams: If True, performs detailed BAM traversal to count mapped reads per branch.
            min_base_quality: Base quality cutoff for depth distribution (-q).
            min_mapping_quality: Mapping quality cutoff for depth distribution (-Q).

        Returns:
            MetricsCollectionResult aggregating all structured stage metrics and reports.
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)

        executed_cmds: List[str] = []
        raw_flagstat: Optional[FlagstatMetrics] = None
        dedup_flagstat: Optional[FlagstatMetrics] = None
        filtered_flagstat: Optional[FlagstatMetrics] = None

        primary_final_bam = filter_result.primary_downstream_bam
        if not primary_final_bam or not primary_final_bam.is_file():
            raise MetricsCollectorError(f"[{sample_id}] Primary downstream BAM is missing for metrics collection.")

        depth_tsv = output_dir / f"{sample_id}.depth.tsv"

        try:
            # 1. Run flagstat on prepared raw/sorted BAM (Baseline)
            if prepared_bam and prepared_bam.is_file():
                raw_stat_txt = output_dir / f"{sample_id}.flagstat.prepared.txt"
                raw_flagstat, cmd1 = cls.run_flagstat(prepared_bam, raw_stat_txt, threads=threads)
                executed_cmds.append(cmd1)

            # 2. Run flagstat on post-dedup BAM
            if dedup_result.dedup_bam and dedup_result.dedup_bam.is_file():
                dedup_stat_txt = output_dir / f"{sample_id}.flagstat.dedup.txt"
                dedup_flagstat, cmd2 = cls.run_flagstat(dedup_result.dedup_bam, dedup_stat_txt, threads=threads)
                executed_cmds.append(cmd2)

            # 3. Run flagstat on final filtered BAM
            filtered_stat_txt = output_dir / f"{sample_id}.flagstat.filtered.txt"
            filtered_flagstat, cmd3 = cls.run_flagstat(primary_final_bam, filtered_stat_txt, threads=threads)
            executed_cmds.append(cmd3)

            # 4. Calculate coverage and depth percentiles on the primary downstream BAM
            coverage_metrics = cls.calculate_coverage_metrics(
                final_bam=primary_final_bam,
                reference_info=reference_info,
                baseline_bam=prepared_bam,
                output_depth_tsv=depth_tsv,
                min_base_quality=min_base_quality,
                min_mapping_quality=min_mapping_quality,
            )

            # 5. Assemble and mathematically validate consolidated AlignmentMetrics
            alignment_metrics = cls.assemble_alignment_metrics(
                branch_results=branch_results,
                filter_result=filter_result,
                dedup_result=dedup_result,
                coverage_metrics=coverage_metrics,
                raw_flagstat=raw_flagstat,
                scan_branch_bams=scan_branch_bams,
            )

            logger.info(
                f"[{sample_id}] Metrics collection completed: Mean Depth={coverage_metrics.mean_depth_post_dedup:.2f}x | "
                f"Breadth >= 3x: {coverage_metrics.breadth_3x_pct:.2f}% | "
                f"Retained Reads: {alignment_metrics.retained_filtered_reads:,}"
            )

            return MetricsCollectionResult(
                sample_id=sample_id,
                raw_flagstat=raw_flagstat,
                filtered_flagstat=filtered_flagstat,
                dedup_flagstat=dedup_flagstat,
                coverage=coverage_metrics,
                provenance=alignment_metrics.provenance,
                alignment_metrics=alignment_metrics,
                depth_tsv=depth_tsv if depth_tsv.is_file() else None,
                command_lines=executed_cmds,
                is_successful=True,
            )

        except Exception as exc:
            logger.error(f"[{sample_id}] Metrics collection failed: {exc}")
            return MetricsCollectionResult(
                sample_id=sample_id,
                raw_flagstat=None,
                filtered_flagstat=None,
                dedup_flagstat=None,
                coverage=CoverageMetrics(),
                provenance=ProvenanceMetrics(),
                alignment_metrics=AlignmentMetrics(),
                depth_tsv=None,
                command_lines=executed_cmds,
                is_successful=False,
                error_message=str(exc),
            )
