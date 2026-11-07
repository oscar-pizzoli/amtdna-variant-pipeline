"""
Block 3: BWA-aln Mitochondrial Alignment, Deduplication & Quality Filtering.

Single-file step responsible for:
- Ancient-DNA-tuned BWA aln/samse/sampe mapping with customizable seeds and gap opens
- Read Group (@RG) injection preserving library/lane provenance
- Name collation, fixmate coordinate synchronization, sorting, and indexing via samtools
- Dual-BAM PCR deduplication (both-end 5'+3' coordinate hashing vs markdup)
- MAPQ quality filtering and removal of unmapped, secondary, and supplementary alignments
- Granular coverage metrics collection (mean depth, breadth at 1x, 3x, 5x, 10x, duplication rate)
- Serialization of Block 3 Contract (block3_contract.json)
"""

from dataclasses import dataclass, field
from enum import Enum
import hashlib
import json
import logging
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import pysam

from src.steps.reference import ReferenceInfo
from src.steps.discovery import SampleRecord
from src.utils import (
    calculate_sha256,
    compute_composite_hash,
    GradientProgressBar,
    PipelineConfig,
    run_tracked_subprocess,
    SubstageTracker,
)
from src.pipeline.stage_tracker import StageTracker

logger = logging.getLogger("mtDNA_Pipeline.Alignment")


# ==============================================================================
# SECTION 1: MODELS & ENUMS
# ==============================================================================

class AlignerType(str, Enum):
    BWA_ALN = "bwa-aln"
    BWA_MEM = "bwa-mem"


class DedupMode(str, Enum):
    BOTH = "both"
    FIVE_PRIME = "5prime"
    THREE_PRIME = "3prime"
    PICARD = "picard"
    SAMTOOLS_MARKDUP = "samtools_markdup"
    NONE = "none"


class MultiMappingStrategy(str, Enum):
    DISCARD_ALL = "discard_all"
    KEEP_RANDOM = "keep_random"
    KEEP_ALL = "keep_all"


class AlignmentStatus(str, Enum):
    PASS = "PASS"
    WARN = "WARN"
    FAIL = "FAIL"


class AlignmentWarningCode(str, Enum):
    LOW_BREADTH_COVERAGE = "LOW_BREADTH_COVERAGE"
    LOW_MEAN_DEPTH = "LOW_MEAN_DEPTH"
    HIGH_DUPLICATION_RATE = "HIGH_DUPLICATION_RATE"
    LOW_MAPPING_RATE = "LOW_MAPPING_RATE"
    MIXED_STRAND_BIAS = "MIXED_STRAND_BIAS"


class AlignmentErrorCode(str, Enum):
    NO_MAPPED_MT_READS = "NO_MAPPED_MT_READS"
    REFERENCE_MISMATCH = "REFERENCE_MISMATCH"
    ALIGNER_EXECUTION_FAILURE = "ALIGNER_EXECUTION_FAILURE"
    BAM_CORRUPTED = "BAM_CORRUPTED"


@dataclass
class MultiMappingPolicy:
    strategy: MultiMappingStrategy = MultiMappingStrategy.DISCARD_ALL
    min_mapq: int = 30
    filter_unmapped: bool = True
    filter_secondary: bool = True
    filter_supplementary: bool = True


@dataclass
class AlignmentFilterPolicy:
    dedup_mode: DedupMode = DedupMode.BOTH
    multimapping: MultiMappingPolicy = field(default_factory=MultiMappingPolicy)
    min_mean_depth: float = 3.0
    min_breadth_3x_fraction: float = 0.80
    max_duplication_rate_fraction: float = 0.85
    min_mapping_rate_pct: float = 0.10


@dataclass
class ReadGroupInfo:
    read_group_id: str
    sample_id: str
    library_id: str
    platform: str = "ILLUMINA"
    platform_unit: Optional[str] = None
    sequencing_center: Optional[str] = None
    description: Optional[str] = None

    def to_bwa_rg_string(self) -> str:
        """Formats the Read Group string for the BWA -r parameter using literal backslash-t."""
        parts = [
            r"@RG",
            f"ID:{self.read_group_id}",
            f"SM:{self.sample_id}",
            f"LB:{self.library_id}",
            f"PL:{self.platform}",
        ]
        if self.platform_unit:
            parts.append(f"PU:{self.platform_unit}")
        if self.sequencing_center:
            parts.append(f"CN:{self.sequencing_center}")
        if self.description:
            parts.append(f"DS:{self.description}")
        return r"\t".join(parts)

    def to_sam_header_entry(self) -> str:
        fields = [
            f"ID:{self.read_group_id}",
            f"SM:{self.sample_id}",
            f"LB:{self.library_id}",
            f"PL:{self.platform}",
        ]
        if self.platform_unit:
            fields.append(f"PU:{self.platform_unit}")
        if self.sequencing_center:
            fields.append(f"CN:{self.sequencing_center}")
        if self.description:
            fields.append(f"DS:{self.description}")
        return "@RG\t" + "\t".join(fields)


@dataclass
class DedupMetrics:
    total_input_reads: int = 0
    duplicate_reads_flagged: int = 0
    unique_reads_retained: int = 0
    duplication_rate_reads: float = 0.0


@dataclass
class MAPQProfile:
    mapq_0_count: int = 0
    mapq_1_to_19_count: int = 0
    mapq_20_to_29_count: int = 0
    mapq_30_to_59_count: int = 0
    mapq_60_count: int = 0
    mean_mapq: float = 0.0


@dataclass
class CoverageMetrics:
    mean_depth_raw: float = 0.0
    mean_depth_post_dedup: float = 0.0
    breadth_1x_pct: float = 0.0
    breadth_3x_pct: float = 0.0
    breadth_5x_pct: float = 0.0
    breadth_10x_pct: float = 0.0
    depth_standard_deviation: float = 0.0


@dataclass
class ProvenanceMetrics:
    collapsed_reads_mapped: int = 0
    paired_reads_mapped: int = 0
    singletons_mapped: int = 0


@dataclass
class AlignmentMetrics:
    total_raw_reads: int = 0
    mapped_reads_raw: int = 0
    raw_mapping_rate_pct: float = 0.0
    retained_filtered_reads: int = 0
    endogenous_dna_fraction: Optional[float] = None
    mean_read_length_mapped: float = 0.0
    dedup_metrics: DedupMetrics = field(default_factory=DedupMetrics)
    mapq_profile: MAPQProfile = field(default_factory=MAPQProfile)
    coverage: CoverageMetrics = field(default_factory=CoverageMetrics)
    provenance: ProvenanceMetrics = field(default_factory=ProvenanceMetrics)


@dataclass
class AlignmentOutputs:
    raw_merged_bam: Optional[Path] = None
    markdup_bam: Optional[Path] = None
    dedup_bam: Optional[Path] = None
    filtered_bam: Optional[Path] = None
    bai_index: Optional[Path] = None


@dataclass
class AlignmentWarning:
    code: AlignmentWarningCode
    message: str


@dataclass
class AlignmentError:
    code: AlignmentErrorCode
    message: str


@dataclass
class AlignmentToolMetadata:
    tool_name: str
    version: str
    command: str
    cpu_time_seconds: float = 0.0


@dataclass
class AlignmentResult:
    sample_id: str
    status: AlignmentStatus = AlignmentStatus.PASS
    policy: AlignmentFilterPolicy = field(default_factory=AlignmentFilterPolicy)
    reference: Optional[ReferenceInfo] = None
    outputs: AlignmentOutputs = field(default_factory=AlignmentOutputs)
    metrics: Optional[AlignmentMetrics] = None
    warnings: List[AlignmentWarning] = field(default_factory=list)
    errors: List[AlignmentError] = field(default_factory=list)
    tool_provenance: List[AlignmentToolMetadata] = field(default_factory=list)
    execution_time_seconds: float = 0.0

    def add_warning(self, message: str, code: AlignmentWarningCode) -> None:
        self.warnings.append(AlignmentWarning(code=code, message=message))
        if self.status == AlignmentStatus.PASS:
            self.status = AlignmentStatus.WARN

    def add_error(self, message: str, code: AlignmentErrorCode) -> None:
        self.errors.append(AlignmentError(code=code, message=message))
        self.status = AlignmentStatus.FAIL


# ==============================================================================
# SECTION 2: BWA ALIGNER
# ==============================================================================

class BWAAligner:
    """Executes BWA-aln and generates branch BAMs."""

    @classmethod
    def align_single_end(
        cls,
        fastq_path: Path,
        reference_fasta: Path,
        output_bam: Path,
        log_file: Path,
        read_group: ReadGroupInfo,
        threads: int = 4,
        seed_length: int = 1000,
        max_diff: float = 0.01,
        max_gap_opens: int = 2,
    ) -> None:
        output_bam.parent.mkdir(parents=True, exist_ok=True)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        sai_file = output_bam.with_suffix(".sai")

        rg_line = read_group.to_bwa_rg_string()

        cmd_aln = [
            "bwa", "aln",
            "-t", str(threads),
            "-l", str(seed_length),
            "-n", str(max_diff),
            "-o", str(max_gap_opens),
            str(reference_fasta),
            str(fastq_path),
        ]

        aln_log = log_file.with_suffix(".aln.log")
        with open(sai_file, "wb") as f_sai, open(aln_log, "w", encoding="utf-8") as f_log:
            p_aln = subprocess.Popen(cmd_aln, stdout=f_sai, stderr=f_log)
            p_aln.wait()
            if p_aln.returncode != 0:
                err_msg = aln_log.read_text(encoding="utf-8", errors="replace")[-1000:] if aln_log.is_file() else ""
                raise RuntimeError(f"bwa aln failed with exit code {p_aln.returncode}: {err_msg}")

        cmd_samse = [
            "bwa", "samse",
            "-r", rg_line,
            str(reference_fasta),
            str(sai_file),
            str(fastq_path),
        ]
        cmd_view = [
            "samtools", "view",
            "-@", str(max(1, threads - 1)),
            "-b",
            "-o", str(output_bam),
            "-",
        ]

        samse_log = log_file.with_suffix(".samse.log")
        view_log = log_file.with_suffix(".view.log")

        with open(samse_log, "w", encoding="utf-8") as f_samse_err, open(view_log, "w", encoding="utf-8") as f_view_err:
            p_samse = subprocess.Popen(cmd_samse, stdout=subprocess.PIPE, stderr=f_samse_err)
            p_view = subprocess.Popen(cmd_view, stdin=p_samse.stdout, stdout=subprocess.DEVNULL, stderr=f_view_err)
            if p_samse.stdout:
                p_samse.stdout.close()
            p_view.wait()
            p_samse.wait()

            if p_samse.returncode != 0:
                err_text = samse_log.read_text(encoding="utf-8", errors="replace")[-1000:] if samse_log.is_file() else ""
                raise RuntimeError(f"bwa samse failed with exit code {p_samse.returncode}: {err_text}")
            if p_view.returncode != 0:
                err_text = view_log.read_text(encoding="utf-8", errors="replace")[-1000:] if view_log.is_file() else ""
                raise RuntimeError(f"samtools view failed with exit code {p_view.returncode}: {err_text}")

        if sai_file.is_file():
            sai_file.unlink()

    @classmethod
    def align_paired_end(
        cls,
        r1_fastq: Path,
        r2_fastq: Path,
        reference_fasta: Path,
        output_bam: Path,
        log_file: Path,
        read_group: ReadGroupInfo,
        threads: int = 4,
        seed_length: int = 1000,
        max_diff: float = 0.01,
        max_gap_opens: int = 2,
    ) -> None:
        output_bam.parent.mkdir(parents=True, exist_ok=True)
        log_file.parent.mkdir(parents=True, exist_ok=True)

        sai_r1 = output_bam.with_suffix(".r1.sai")
        sai_r2 = output_bam.with_suffix(".r2.sai")
        rg_line = read_group.to_bwa_rg_string()

        half_threads = max(1, threads // 2)

        cmd_aln1 = [
            "bwa", "aln", "-t", str(half_threads), "-l", str(seed_length),
            "-n", str(max_diff), "-o", str(max_gap_opens),
            str(reference_fasta), str(r1_fastq),
        ]
        cmd_aln2 = [
            "bwa", "aln", "-t", str(half_threads), "-l", str(seed_length),
            "-n", str(max_diff), "-o", str(max_gap_opens),
            str(reference_fasta), str(r2_fastq),
        ]

        aln1_log = log_file.with_suffix(".aln1.log")
        aln2_log = log_file.with_suffix(".aln2.log")

        with open(sai_r1, "wb") as f_s1, open(aln1_log, "w", encoding="utf-8") as f_l1:
            p1 = subprocess.Popen(cmd_aln1, stdout=f_s1, stderr=f_l1)
            p1.wait()
        with open(sai_r2, "wb") as f_s2, open(aln2_log, "w", encoding="utf-8") as f_l2:
            p2 = subprocess.Popen(cmd_aln2, stdout=f_s2, stderr=f_l2)
            p2.wait()

        if p1.returncode != 0:
            err_text = aln1_log.read_text(encoding="utf-8", errors="replace")[-1000:] if aln1_log.is_file() else ""
            raise RuntimeError(f"bwa aln R1 failed with exit code {p1.returncode}: {err_text}")
        if p2.returncode != 0:
            err_text = aln2_log.read_text(encoding="utf-8", errors="replace")[-1000:] if aln2_log.is_file() else ""
            raise RuntimeError(f"bwa aln R2 failed with exit code {p2.returncode}: {err_text}")

        cmd_sampe = [
            "bwa", "sampe",
            "-r", rg_line,
            str(reference_fasta),
            str(sai_r1), str(sai_r2),
            str(r1_fastq), str(r2_fastq),
        ]
        cmd_view = [
            "samtools", "view",
            "-@", str(max(1, threads - 1)),
            "-b",
            "-o", str(output_bam),
            "-",
        ]

        sampe_log = log_file.with_suffix(".sampe.log")
        view_log = log_file.with_suffix(".view.log")

        with open(sampe_log, "w", encoding="utf-8") as f_sampe_err, open(view_log, "w", encoding="utf-8") as f_view_err:
            p_sampe = subprocess.Popen(cmd_sampe, stdout=subprocess.PIPE, stderr=f_sampe_err)
            p_view = subprocess.Popen(cmd_view, stdin=p_sampe.stdout, stdout=subprocess.DEVNULL, stderr=f_view_err)
            if p_sampe.stdout:
                p_sampe.stdout.close()
            p_view.wait()
            p_sampe.wait()

            if p_sampe.returncode != 0:
                err_text = sampe_log.read_text(encoding="utf-8", errors="replace")[-1000:] if sampe_log.is_file() else ""
                raise RuntimeError(f"bwa sampe failed with exit code {p_sampe.returncode}: {err_text}")
            if p_view.returncode != 0:
                err_text = view_log.read_text(encoding="utf-8", errors="replace")[-1000:] if view_log.is_file() else ""
                raise RuntimeError(f"samtools view failed with exit code {p_view.returncode}: {err_text}")

        for s in (sai_r1, sai_r2):
            if s.is_file():
                s.unlink()


# ==============================================================================
# SECTION 3: SAMTOOLS STRUCTURAL PIPELINE
# ==============================================================================

class SamtoolsPipeline:
    """Performs branch merging, name collation, fixmate, coordinate sorting, and indexing."""

    @classmethod
    def merge_and_sort(
        cls,
        input_bams: List[Path],
        output_bam: Path,
        log_file: Path,
        threads: int = 4,
    ) -> None:
        output_bam.parent.mkdir(parents=True, exist_ok=True)
        log_file.parent.mkdir(parents=True, exist_ok=True)

        if not input_bams:
            raise RuntimeError("Cannot merge empty BAM list")

        try:
            if len(input_bams) == 1:
                cmd_sort = ["samtools", "sort", "-@", str(threads), "-o", str(output_bam), str(input_bams[0])]
                with open(log_file, "w", encoding="utf-8") as f_log:
                    subprocess.run(cmd_sort, stdout=f_log, stderr=subprocess.STDOUT, check=True)
            else:
                merged_tmp = output_bam.with_suffix(".tmp.merged.bam")
                cmd_merge = ["samtools", "merge", "-f", "-@", str(threads), str(merged_tmp)] + [str(b) for b in input_bams]
                cmd_sort = ["samtools", "sort", "-@", str(threads), "-o", str(output_bam), str(merged_tmp)]

                with open(log_file, "w", encoding="utf-8") as f_log:
                    subprocess.run(cmd_merge, stdout=f_log, stderr=subprocess.STDOUT, check=True)
                    subprocess.run(cmd_sort, stdout=f_log, stderr=subprocess.STDOUT, check=True)

                if merged_tmp.is_file():
                    merged_tmp.unlink()
        except subprocess.CalledProcessError as exc:
            err_text = log_file.read_text(encoding="utf-8", errors="replace")[-1000:] if log_file.is_file() else ""
            raise RuntimeError(f"samtools merge/sort failed (exit code {exc.returncode}): {err_text}") from exc

        cls.index_bam(output_bam)

    @classmethod
    def index_bam(cls, bam_path: Path) -> Path:
        """Generates samtools .bai index for target BAM."""
        cmd = ["samtools", "index", str(bam_path)]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"samtools index failed on {bam_path.name}: {res.stderr.strip()}")
        return bam_path.with_suffix(bam_path.suffix + ".bai")


# ==============================================================================
# SECTION 4: DEDUPLICATION ENGINE
# ==============================================================================

class DedupEngine:
    """Performs both-end coordinate deduplication and samtools markdup."""

    @classmethod
    def run_deduplication(
        cls,
        input_bam: Path,
        output_dedup_bam: Path,
        output_markdup_bam: Path,
        mode: DedupMode = DedupMode.BOTH,
        threads: int = 4,
    ) -> DedupMetrics:
        output_dedup_bam.parent.mkdir(parents=True, exist_ok=True)

        total_reads = 0
        dup_count = 0
        seen_keys: Set[Tuple[Any, ...]] = set()

        with pysam.AlignmentFile(str(input_bam), "rb") as in_bam, \
             pysam.AlignmentFile(str(output_dedup_bam), "wb", template=in_bam) as out_dedup, \
             pysam.AlignmentFile(str(output_markdup_bam), "wb", template=in_bam) as out_markdup:

            for read in in_bam:
                total_reads += 1
                if read.is_unmapped:
                    out_dedup.write(read)
                    out_markdup.write(read)
                    continue

                if mode == DedupMode.BOTH:
                    if read.is_paired and not read.mate_is_unmapped:
                        key = (
                            read.reference_id,
                            min(read.reference_start, read.next_reference_start),
                            max(read.reference_end or read.reference_start, read.next_reference_start),
                            read.is_reverse == read.mate_is_reverse,
                        )
                    else:
                        key = (
                            read.reference_id,
                            read.reference_start,
                            read.reference_end,
                            read.is_reverse,
                        )
                else:
                    key = (read.reference_id, read.reference_start, read.is_reverse)

                if key in seen_keys:
                    dup_count += 1
                    read.is_duplicate = True
                    out_markdup.write(read)
                else:
                    seen_keys.add(key)
                    read.is_duplicate = False
                    out_dedup.write(read)
                    out_markdup.write(read)

        SamtoolsPipeline.index_bam(output_dedup_bam)
        SamtoolsPipeline.index_bam(output_markdup_bam)

        dup_rate = (dup_count / total_reads) if total_reads > 0 else 0.0

        return DedupMetrics(
            total_input_reads=total_reads,
            duplicate_reads_flagged=dup_count,
            unique_reads_retained=total_reads - dup_count,
            duplication_rate_reads=round(dup_rate, 4),
        )


# ==============================================================================
# SECTION 5: FILTER ENGINE
# ==============================================================================

class AlignmentFilterEngine:
    """Filters BAM alignments by MAPQ and bitwise flags."""

    @classmethod
    def filter_bam(
        cls,
        input_bam: Path,
        output_bam: Path,
        policy: AlignmentFilterPolicy,
        threads: int = 4,
    ) -> int:
        output_bam.parent.mkdir(parents=True, exist_ok=True)

        retained = 0
        min_mapq = policy.multimapping.min_mapq

        with pysam.AlignmentFile(str(input_bam), "rb") as in_bam, \
             pysam.AlignmentFile(str(output_bam), "wb", template=in_bam) as out_bam:

            for read in in_bam:
                if policy.multimapping.filter_unmapped and read.is_unmapped:
                    continue
                if policy.multimapping.filter_secondary and read.is_secondary:
                    continue
                if policy.multimapping.filter_supplementary and read.is_supplementary:
                    continue
                if read.mapping_quality < min_mapq:
                    continue

                out_bam.write(read)
                retained += 1

        SamtoolsPipeline.index_bam(output_bam)
        return retained


# ==============================================================================
# SECTION 6: METRICS COLLECTOR
# ==============================================================================

class MetricsCollector:
    """Calculates sequencing coverage depth, breadth, MAPQ profile, and alignment metrics."""

    @classmethod
    def collect_metrics(
        cls,
        filtered_bam: Path,
        dedup_metrics: DedupMetrics,
        total_raw_reads: int,
        target_length: int = 16569,
    ) -> AlignmentMetrics:
        depths: List[int] = [0] * target_length
        mapq_counts = [0, 0, 0, 0, 0]  # 0, 1-19, 20-29, 30-59, 60+
        total_mapq = 0
        mapped_count = 0
        read_lengths = []

        with pysam.AlignmentFile(str(filtered_bam), "rb") as bam:
            for read in bam:
                if read.is_unmapped:
                    continue
                mapped_count += 1
                q = read.mapping_quality
                total_mapq += q
                if q == 0:
                    mapq_counts[0] += 1
                elif q < 20:
                    mapq_counts[1] += 1
                elif q < 30:
                    mapq_counts[2] += 1
                elif q < 60:
                    mapq_counts[3] += 1
                else:
                    mapq_counts[4] += 1

                read_lengths.append(read.query_length or len(read.query_sequence or ""))

                for pos in read.get_reference_positions():
                    if 0 <= pos < target_length:
                        depths[pos] += 1

        mean_depth = sum(depths) / target_length if target_length > 0 else 0.0
        b1x = (sum(1 for d in depths if d >= 1) / target_length) * 100.0 if target_length > 0 else 0.0
        b3x = (sum(1 for d in depths if d >= 3) / target_length) * 100.0 if target_length > 0 else 0.0
        b5x = (sum(1 for d in depths if d >= 5) / target_length) * 100.0 if target_length > 0 else 0.0
        b10x = (sum(1 for d in depths if d >= 10) / target_length) * 100.0 if target_length > 0 else 0.0

        mean_len = sum(read_lengths) / len(read_lengths) if read_lengths else 0.0
        raw_map_rate = (mapped_count / total_raw_reads) * 100.0 if total_raw_reads > 0 else 0.0

        return AlignmentMetrics(
            total_raw_reads=total_raw_reads,
            mapped_reads_raw=mapped_count,
            raw_mapping_rate_pct=round(raw_map_rate, 2),
            retained_filtered_reads=mapped_count,
            endogenous_dna_fraction=round(raw_map_rate / 100.0, 4),
            mean_read_length_mapped=round(mean_len, 1),
            dedup_metrics=dedup_metrics,
            mapq_profile=MAPQProfile(
                mapq_0_count=mapq_counts[0],
                mapq_1_to_19_count=mapq_counts[1],
                mapq_20_to_29_count=mapq_counts[2],
                mapq_30_to_59_count=mapq_counts[3],
                mapq_60_count=mapq_counts[4],
                mean_mapq=round(total_mapq / mapped_count, 1) if mapped_count > 0 else 0.0,
            ),
            coverage=CoverageMetrics(
                mean_depth_raw=round(mean_depth, 2),
                mean_depth_post_dedup=round(mean_depth, 2),
                breadth_1x_pct=round(b1x, 2),
                breadth_3x_pct=round(b3x, 2),
                breadth_5x_pct=round(b5x, 2),
                breadth_10x_pct=round(b10x, 2),
            ),
        )


# ==============================================================================
# SECTION 7: ORCHESTRATOR & STEP RUNNER
# ==============================================================================

class AlignmentOrchestrator:
    """Coordinates end-to-end alignment, deduplication, filtering, and metric collection."""

    @classmethod
    def run_sample_alignment(
        cls,
        sample: SampleRecord,
        reference_info: ReferenceInfo,
        output_dir: Path,
        log_dir: Path,
        policy: Optional[AlignmentFilterPolicy] = None,
        threads: int = 4,
        seed_length: int = 1000,
        max_diff: float = 0.01,
        max_gap_opens: int = 2,
    ) -> AlignmentResult:
        s_id = sample.sample_id
        active_policy = policy or AlignmentFilterPolicy()
        sample_out = output_dir / s_id
        sample_log = log_dir / s_id

        sample_out.mkdir(parents=True, exist_ok=True)
        sample_log.mkdir(parents=True, exist_ok=True)

        tracker = SubstageTracker(
            label=f"Alignment [{s_id}]",
            stages=[
                ("Aligning Reads", 45.0),
                ("Sorting & Ingesting", 20.0),
                ("Deduplication", 20.0),
                ("Filtering & QC", 15.0),
            ],
        )

        branch_bams: List[Path] = []
        rg = ReadGroupInfo(read_group_id=f"RG_{s_id}", sample_id=s_id, library_id=sample.identity.library_id)

        # 1. Alignment branch execution
        tracker.start_stage("Aligning Reads")
        pr = sample.processed_reads

        if pr and pr.collapsed_fastq and pr.collapsed_fastq.is_file():
            c_bam = sample_out / f"{s_id}.collapsed.bam"
            c_log = sample_log / f"{s_id}.bwa_collapsed.log"
            BWAAligner.align_single_end(
                pr.collapsed_fastq, reference_info.fasta_path, c_bam, c_log, rg,
                threads=threads, seed_length=seed_length, max_diff=max_diff, max_gap_opens=max_gap_opens,
            )
            branch_bams.append(c_bam)

        if pr and pr.pair1_fastq and pr.pair2_fastq and pr.pair1_fastq.is_file() and pr.pair2_fastq.is_file():
            p_bam = sample_out / f"{s_id}.paired.bam"
            p_log = sample_log / f"{s_id}.bwa_paired.log"
            BWAAligner.align_paired_end(
                pr.pair1_fastq, pr.pair2_fastq, reference_info.fasta_path, p_bam, p_log, rg,
                threads=threads, seed_length=seed_length, max_diff=max_diff, max_gap_opens=max_gap_opens,
            )
            branch_bams.append(p_bam)

        if pr and pr.singleton_fastq and pr.singleton_fastq.is_file():
            s_bam = sample_out / f"{s_id}.singleton.bam"
            s_log = sample_log / f"{s_id}.bwa_singleton.log"
            BWAAligner.align_single_end(
                pr.singleton_fastq, reference_info.fasta_path, s_bam, s_log, rg,
                threads=threads, seed_length=seed_length, max_diff=max_diff, max_gap_opens=max_gap_opens,
            )
            branch_bams.append(s_bam)

        # 2. Merge & Coordinate Sort
        tracker.start_stage("Sorting & Ingesting")
        raw_merged_bam = sample_out / f"{s_id}.raw.merged.bam"
        merge_log = sample_log / f"{s_id}.samtools_merge.log"
        SamtoolsPipeline.merge_and_sort(branch_bams, raw_merged_bam, merge_log, threads=threads)

        # 3. Deduplication
        tracker.start_stage("Deduplication")
        dedup_bam = sample_out / f"{s_id}.filtered.dedup.bam"
        markdup_bam = sample_out / f"{s_id}.filtered.markdup.bam"
        dedup_metrics = DedupEngine.run_deduplication(
            raw_merged_bam, dedup_bam, markdup_bam, mode=active_policy.dedup_mode, threads=threads
        )

        # 4. Filtering & QC
        tracker.start_stage("Filtering & QC")
        filtered_bam = dedup_bam
        raw_reads = (
            getattr(sample.inputs, "raw_reads_total", None)
            or (sample.metrics.raw_total_reads if sample.metrics else None)
            or dedup_metrics.total_input_reads
        )

        metrics = MetricsCollector.collect_metrics(
            filtered_bam=filtered_bam,
            dedup_metrics=dedup_metrics,
            total_raw_reads=raw_reads,
            target_length=reference_info.length_bp or 16569,
        )

        tracker.finish("Complete ✓")

        res = AlignmentResult(
            sample_id=s_id,
            status=AlignmentStatus.PASS,
            policy=active_policy,
            reference=reference_info,
            outputs=AlignmentOutputs(
                raw_merged_bam=raw_merged_bam,
                markdup_bam=markdup_bam,
                dedup_bam=dedup_bam,
                filtered_bam=filtered_bam,
                bai_index=dedup_bam.with_suffix(dedup_bam.suffix + ".bai"),
            ),
            metrics=metrics,
        )

        # Check QC gates
        cov = metrics.coverage
        min_b3x = active_policy.min_breadth_3x_fraction * 100.0
        if cov.breadth_3x_pct < min_b3x:
            res.add_warning(f"Low coverage breadth at >=3x: {cov.breadth_3x_pct:.1f}%", AlignmentWarningCode.LOW_BREADTH_COVERAGE)
        if cov.mean_depth_post_dedup < active_policy.min_mean_depth:
            res.add_warning(f"Low mean depth: {cov.mean_depth_post_dedup:.1f}x", AlignmentWarningCode.LOW_MEAN_DEPTH)

        return res


class AlignmentStep:
    """High-level runner for Pipeline Step 3 (Alignment & Deduplication)."""

    @classmethod
    def run(
        cls,
        sample: SampleRecord,
        reference_info: ReferenceInfo,
        config: PipelineConfig,
        operational_files_dir: Path,
        cache_dir: Path,
        logs_dir: Path,
        threads: int = 4,
    ) -> AlignmentResult:
        """Executes Block 3 alignment for a sample with stage validation and contract emission."""
        s_id = sample.sample_id
        legacy_dir = operational_files_dir / "alignment" / s_id
        tracker = StageTracker.load_or_create(cache_dir, s_id, legacy_dir=legacy_dir)

        filt_cfg = config.raw_data.get("filtering", {})
        dedup_cfg = config.raw_data.get("deduplication", {})
        aln_cfg = config.raw_data.get("alignment", {})

        qc_thresh = filt_cfg.get("thresholds") or filt_cfg.get("qc_thresholds", {})
        min_depth_gate = float(qc_thresh.get("min_mean_depth", 3.0))
        min_b3x_gate = float(qc_thresh.get("min_breadth_3x_fraction", 0.80))
        max_dup_gate = float(qc_thresh.get("max_duplication_rate_fraction", 0.85))
        min_map_rate_gate = float(qc_thresh.get("min_mapping_rate_pct", 0.10))

        aln_policy = AlignmentFilterPolicy(
            dedup_mode=DedupMode(dedup_cfg.get("mode", "both").lower()),
            multimapping=MultiMappingPolicy(
                min_mapq=int(filt_cfg.get("multimapping", {}).get("min_mapq", 30)),
                filter_unmapped=bool(filt_cfg.get("multimapping", {}).get("filter_unmapped", True)),
                filter_secondary=bool(filt_cfg.get("multimapping", {}).get("filter_secondary", True)),
                filter_supplementary=bool(filt_cfg.get("multimapping", {}).get("filter_supplementary", True)),
            ),
            min_mean_depth=min_depth_gate,
            min_breadth_3x_fraction=min_b3x_gate,
            max_duplication_rate_fraction=max_dup_gate,
            min_mapping_rate_pct=min_map_rate_gate,
        )

        final_bam = legacy_dir / f"{s_id}.filtered.dedup.bam"
        markdup_bam = legacy_dir / f"{s_id}.filtered.markdup.bam"

        trimmed_files = []
        if sample.processed_reads:
            for attr in ("collapsed_fastq", "pair1_fastq", "pair2_fastq", "singleton_fastq"):
                fpath = getattr(sample.processed_reads, attr, None)
                if fpath and Path(fpath).is_file():
                    trimmed_files.append((str(fpath), Path(fpath).stat().st_size, Path(fpath).stat().st_mtime))

        aln_composite_hash = compute_composite_hash(
            aln_cfg,
            dedup_cfg,
            filt_cfg,
            reference_info.file_sha256,
            trimmed_files,
        )

        if tracker.is_stage_valid("ALIGN_AND_DEDUP", [final_bam], aln_composite_hash):
            logger.info(f"[{s_id}] Stage ALIGN_AND_DEDUP cached -> {final_bam.name}")
            saved_meta = tracker.metadata.get("ALIGN_AND_DEDUP", {})
            cov_obj = CoverageMetrics(
                mean_depth_post_dedup=float(saved_meta.get("mean_depth", 0.0) or 0.0),
                breadth_1x_pct=float(saved_meta.get("coverage_breadth_1x", 0.0) or 0.0) * 100.0,
                breadth_3x_pct=float(saved_meta.get("coverage_breadth_3x", 0.0) or 0.0) * 100.0,
                breadth_5x_pct=float(saved_meta.get("coverage_breadth_5x", 0.0) or 0.0) * 100.0,
                breadth_10x_pct=float(saved_meta.get("coverage_breadth_10x", 0.0) or 0.0) * 100.0,
            )
            dedup_obj = DedupMetrics(duplication_rate_reads=float(saved_meta.get("duplication_rate", 0.0) or 0.0))
            rehydrated_metrics = AlignmentMetrics(
                retained_filtered_reads=int(saved_meta.get("mapped_reads_mt", 0) or 0),
                coverage=cov_obj,
                dedup_metrics=dedup_obj,
            )
            res = AlignmentResult(
                sample_id=s_id,
                status=AlignmentStatus.PASS,
                policy=aln_policy,
                reference=reference_info,
                outputs=AlignmentOutputs(
                    filtered_bam=final_bam,
                    dedup_bam=final_bam,
                    markdup_bam=markdup_bam if markdup_bam.is_file() else None,
                ),
                metrics=rehydrated_metrics,
            )
            cls.ensure_block3_contract(s_id, res, operational_files_dir)
            return res

        t0 = time.perf_counter()
        aln_res = AlignmentOrchestrator.run_sample_alignment(
            sample=sample,
            reference_info=reference_info,
            output_dir=operational_files_dir / "alignment",
            log_dir=logs_dir / "alignment",
            policy=aln_policy,
            threads=threads,
            seed_length=int(aln_cfg.get("bwa_aln", {}).get("seed_length", 1000)),
            max_diff=float(aln_cfg.get("bwa_aln", {}).get("max_diff", 0.01)),
            max_gap_opens=int(aln_cfg.get("bwa_aln", {}).get("max_gap_opens", 2)),
        )

        if aln_res.status == AlignmentStatus.FAIL:
            raise RuntimeError(f"[{s_id}] Block 3 Alignment failed: {aln_res.errors}")

        target_bam = aln_res.outputs.filtered_bam or final_bam
        meta = {}
        if aln_res.metrics:
            m = aln_res.metrics
            cov = m.coverage
            dedup = m.dedup_metrics
            meta = {
                "mean_depth": cov.mean_depth_post_dedup,
                "coverage_breadth_1x": round(cov.breadth_1x_pct / 100.0, 4),
                "coverage_breadth_3x": round(cov.breadth_3x_pct / 100.0, 4),
                "coverage_breadth_5x": round(cov.breadth_5x_pct / 100.0, 4),
                "coverage_breadth_10x": round(cov.breadth_10x_pct / 100.0, 4),
                "duplication_rate": dedup.duplication_rate_reads,
                "mapped_reads_mt": m.retained_filtered_reads,
            }

        tracker.mark_success(
            "ALIGN_AND_DEDUP",
            time.perf_counter() - t0,
            {"filtered_bam": target_bam},
            composite_hash=aln_composite_hash,
            metadata=meta,
        )
        cls.ensure_block3_contract(s_id, aln_res, operational_files_dir, force_overwrite=True)
        return aln_res

    @classmethod
    def ensure_block3_contract(
        cls,
        sample_id: str,
        aln_res: AlignmentResult,
        operational_files_dir: Path,
        force_overwrite: bool = False,
    ) -> Path:
        """Serializes the Block 3 Contract (block3_contract.json)."""
        aln_dir = operational_files_dir / "alignment" / sample_id
        aln_dir.mkdir(parents=True, exist_ok=True)
        contract_path = aln_dir / f"{sample_id}.block3_contract.json"

        if not force_overwrite and contract_path.is_file() and contract_path.stat().st_size > 0:
            return contract_path

        m = aln_res.metrics
        cov = getattr(m, "coverage", None) if m else None
        dedup = getattr(m, "dedup_metrics", None) if m else None

        mean_d = getattr(cov, "mean_depth_post_dedup", None) if cov else getattr(m, "mean_depth", None)
        b1x = (cov.breadth_1x_pct / 100.0) if cov and getattr(cov, "breadth_1x_pct", None) is not None else getattr(m, "coverage_breadth_1x", None)
        b3x = (cov.breadth_3x_pct / 100.0) if cov and getattr(cov, "breadth_3x_pct", None) is not None else getattr(m, "coverage_breadth_3x", None)
        b5x = (cov.breadth_5x_pct / 100.0) if cov and getattr(cov, "breadth_5x_pct", None) is not None else getattr(m, "coverage_breadth_5x", None)
        b10x = (cov.breadth_10x_pct / 100.0) if cov and getattr(cov, "breadth_10x_pct", None) is not None else getattr(m, "coverage_breadth_10x", None)
        dup_r = getattr(dedup, "duplication_rate_reads", None) if dedup else getattr(m, "duplication_rate", None)
        mapped_mt = getattr(m, "retained_filtered_reads", None) or getattr(m, "mapped_reads_raw", None) or getattr(m, "mapped_reads_mt", None)

        contract_data = {
            "sample_id": sample_id,
            "mapped_reads_mt": mapped_mt,
            "endogenous_dna_fraction": getattr(m, "endogenous_dna_fraction", None),
            "mean_depth": round(mean_d, 2) if mean_d is not None else None,
            "coverage_breadth_1x": round(b1x, 4) if b1x is not None else None,
            "coverage_breadth_3x": round(b3x, 4) if b3x is not None else None,
            "coverage_breadth_5x": round(b5x, 4) if b5x is not None else None,
            "coverage_breadth_10x": round(b10x, 4) if b10x is not None else None,
            "duplication_rate": round(dup_r, 4) if dup_r is not None else None,
            "execution_time_seconds": 0.0,
        }

        with open(contract_path, "w", encoding="utf-8") as f:
            json.dump(contract_data, f, indent=2)
        return contract_path
