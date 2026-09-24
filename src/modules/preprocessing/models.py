"""Data models and domain schemas for sample tracking, input fingerprinting, and preprocessing metrics."""

from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional


class DataType(str, Enum):
    """Supported input sequencing data types."""

    FASTQ = "FASTQ"
    BAM = "BAM"


class ValidationStatus(str, Enum):
    """Sample processing and quality validation status."""

    UNTESTED = "UNTESTED"
    PASS = "PASS"
    WARN = "WARN"
    FAIL = "FAIL"


class PreprocessingTool(str, Enum):
    """Supported adapter trimming and merging engines."""

    ADAPTER_REMOVAL = "AdapterRemoval"
    FASTP = "fastp"
    NONE = "none"


class OutputProfile(str, Enum):
    """Classifies the composition of reads surviving preprocessing."""

    EMPTY = "EMPTY"
    MERGED_ONLY = "MERGED_ONLY"
    PAIRED_ONLY = "PAIRED_ONLY"
    SINGLETON_ONLY = "SINGLETON_ONLY"
    MIXED = "MIXED"


@dataclass
class SampleIdentity:
    """Represents the biological and experimental provenance of a sample."""

    sample_id: str
    library_id: Optional[str] = None
    run_id: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Serializes identity fields to dictionary."""
        return {
            "sample_id": self.sample_id,
            "library_id": self.library_id or self.sample_id,
            "run_id": self.run_id or "run_01",
        }


@dataclass
class FileFingerprint:
    """Stores cryptographic integrity, file size, and read count for input/output files."""

    filename: str
    filepath: Path
    file_size_bytes: int
    sha256_checksum: Optional[str] = None
    read_count: Optional[int] = None
    is_input: bool = True

    def to_dict(self) -> Dict[str, Any]:
        """Serializes fingerprint metadata to dictionary."""
        return {
            "filename": self.filename,
            "filepath": str(self.filepath),
            "file_size_bytes": self.file_size_bytes,
            "sha256": self.sha256_checksum,
            "read_count": self.read_count,
            "is_input": self.is_input,
        }


@dataclass
class InputFiles:
    """Strongly validated container for input file paths with strict consistency checks."""

    data_type: DataType
    is_paired: bool = True
    r1_path: Optional[Path] = None
    r2_path: Optional[Path] = None
    bam_path: Optional[Path] = None

    def __post_init__(self) -> None:
        """Enforces domain consistency on input paths."""
        if self.data_type == DataType.BAM:
            if self.bam_path is None:
                raise ValueError("InputFiles configured as BAM but 'bam_path' is None.")
            if self.r1_path is not None or self.r2_path is not None:
                raise ValueError("InputFiles configured as BAM cannot have 'r1_path' or 'r2_path' assigned.")
            self.is_paired = False

        elif self.data_type == DataType.FASTQ:
            if self.bam_path is not None:
                raise ValueError("InputFiles configured as FASTQ cannot have 'bam_path' assigned.")
            if self.r1_path is None:
                raise ValueError("InputFiles configured as FASTQ requires at least 'r1_path'.")
            if self.is_paired:
                if self.r2_path is None:
                    raise ValueError("Paired-end FASTQ input requires both 'r1_path' and 'r2_path'.")
            else:
                if self.r2_path is not None:
                    raise ValueError("Single-end FASTQ input cannot have 'r2_path' assigned.")

    def to_dict(self) -> Dict[str, Any]:
        """Serializes input paths to dictionary."""
        return {
            "data_type": self.data_type.value,
            "is_paired": self.is_paired,
            "r1_path": str(self.r1_path) if self.r1_path else None,
            "r2_path": str(self.r2_path) if self.r2_path else None,
            "bam_path": str(self.bam_path) if self.bam_path else None,
        }


@dataclass
class ProcessedReads:
    """Granular container for post-trimming and merged sequence file paths."""

    merged_fastq: Optional[Path] = None
    paired_r1_fastq: Optional[Path] = None
    paired_r2_fastq: Optional[Path] = None
    singleton_fastq: Optional[Path] = None
    output_profile: OutputProfile = OutputProfile.EMPTY
    output_fingerprints: List[FileFingerprint] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Serializes processed read paths to dictionary."""
        return {
            "merged_fastq": str(self.merged_fastq) if self.merged_fastq else None,
            "paired_r1_fastq": str(self.paired_r1_fastq) if self.paired_r1_fastq else None,
            "paired_r2_fastq": str(self.paired_r2_fastq) if self.paired_r2_fastq else None,
            "singleton_fastq": str(self.singleton_fastq) if self.singleton_fastq else None,
            "output_profile": self.output_profile.value,
            "output_fingerprints": [fp.to_dict() for fp in self.output_fingerprints],
        }


@dataclass
class ToolMetadata:
    """Captures command-line invocations, raw version strings, and execution parameters."""

    tool_name: PreprocessingTool = PreprocessingTool.NONE
    tool_version: Optional[str] = None
    tool_raw_version: Optional[str] = None
    command_line: Optional[str] = None
    parameters: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Serializes tool metadata to dictionary."""
        return {
            "tool_name": self.tool_name.value,
            "tool_version": self.tool_version or "unknown",
            "tool_raw_version": self.tool_raw_version or "unknown",
            "command_line": self.command_line or "none",
            "parameters": self.parameters,
        }


@dataclass
class PreprocessingMetrics:
    """Detailed read counts, merging statistics, and sequencing quality metrics."""

    raw_read_pairs: Optional[int] = None
    raw_total_reads: Optional[int] = None
    surviving_original_reads: Optional[int] = None       # 2 * surviving_pairs + 2 * merged + singletons
    surviving_sequenced_fragments: Optional[int] = None  # surviving_pairs + merged + singletons
    surviving_read_pairs: Optional[int] = None           # Surviving unmerged PE pairs
    merged_fragments: Optional[int] = None               # Successfully collapsed pairs
    singleton_reads: Optional[int] = None                # Surviving single-end reads (discarded mate)
    discarded_reads: Optional[int] = None                # Discarded reads due to length/quality
    merged_fraction_raw: Optional[float] = None          # merged_fragments / raw_read_pairs
    merged_fraction_surviving: Optional[float] = None    # merged_fragments / (merged_fragments + surviving_read_pairs)
    adapter_trimmed_reads_pct: Optional[float] = None    # Percentage of reads containing adapters (0.0 - 100.0)
    q20_rate_post: Optional[float] = None                # Fraction of bases >= Q20 in surviving reads (0.0 - 1.0)
    q30_rate_post: Optional[float] = None                # Fraction of bases >= Q30 in surviving reads (0.0 - 1.0)
    gc_content_pct: Optional[float] = None               # GC base percentage in surviving reads (0.0 - 100.0)
    mean_read_length_pre: Optional[float] = None         # Mean read length before trimming
    mean_read_length_post: Optional[float] = None        # Mean read length of retained reads

    def to_dict(self) -> Dict[str, Any]:
        """Serializes metrics to dictionary."""
        return {
            "raw_read_pairs": self.raw_read_pairs,
            "raw_total_reads": self.raw_total_reads,
            "surviving_original_reads": self.surviving_original_reads,
            "surviving_sequenced_fragments": self.surviving_sequenced_fragments,
            "surviving_read_pairs": self.surviving_read_pairs,
            "merged_fragments": self.merged_fragments,
            "singleton_reads": self.singleton_reads,
            "discarded_reads": self.discarded_reads,
            "merged_fraction_raw": round(self.merged_fraction_raw, 4) if self.merged_fraction_raw is not None else None,
            "merged_fraction_surviving": round(self.merged_fraction_surviving, 4) if self.merged_fraction_surviving is not None else None,
            "adapter_trimmed_reads_pct": round(self.adapter_trimmed_reads_pct, 2) if self.adapter_trimmed_reads_pct is not None else None,
            "q20_rate_post": round(self.q20_rate_post, 4) if self.q20_rate_post is not None else None,
            "q30_rate_post": round(self.q30_rate_post, 4) if self.q30_rate_post is not None else None,
            "gc_content_pct": round(self.gc_content_pct, 2) if self.gc_content_pct is not None else None,
            "mean_read_length_pre": round(self.mean_read_length_pre, 2) if self.mean_read_length_pre is not None else None,
            "mean_read_length_post": round(self.mean_read_length_post, 2) if self.mean_read_length_post is not None else None,
        }


@dataclass
class ValidationResult:
    """Validation diagnostics and audit issues for the sample."""

    status: ValidationStatus = ValidationStatus.UNTESTED
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    notes: Optional[str] = None

    def add_warning(self, msg: str) -> None:
        """Adds a warning and updates status if not in FAIL."""
        self.warnings.append(msg)
        if self.status != ValidationStatus.FAIL:
            self.status = ValidationStatus.WARN

    def add_error(self, msg: str) -> None:
        """Adds an error and transitions status strictly to FAIL."""
        self.errors.append(msg)
        self.status = ValidationStatus.FAIL

    def to_dict(self) -> Dict[str, Any]:
        """Serializes validation result to dictionary."""
        return {
            "status": self.status.value,
            "warnings_count": len(self.warnings),
            "warnings": self.warnings,
            "errors_count": len(self.errors),
            "errors": self.errors,
            "notes": self.notes,
        }


@dataclass
class SampleRecord:
    """Top-level composable entity tracking a sample across all pipeline blocks."""

    identity: SampleIdentity
    inputs: InputFiles
    processed_reads: ProcessedReads = field(default_factory=ProcessedReads)
    fingerprints: List[FileFingerprint] = field(default_factory=list)
    metrics: PreprocessingMetrics = field(default_factory=PreprocessingMetrics)
    tool_meta: ToolMetadata = field(default_factory=ToolMetadata)
    validation: ValidationResult = field(default_factory=ValidationResult)
    processing_started: Optional[str] = None
    processing_finished: Optional[str] = None
    runtime_seconds: Optional[float] = None

    @property
    def sample_id(self) -> str:
        return self.identity.sample_id

    @property
    def data_type(self) -> DataType:
        return self.inputs.data_type

    @property
    def status(self) -> ValidationStatus:
        return self.validation.status

    def mark_started(self) -> None:
        """Records ISO start timestamp."""
        self.processing_started = datetime.now(timezone.utc).isoformat()

    def mark_finished(self) -> None:
        """Records ISO finish timestamp and calculates duration."""
        self.processing_finished = datetime.now(timezone.utc).isoformat()
        if self.processing_started:
            start_dt = datetime.fromisoformat(self.processing_started)
            end_dt = datetime.fromisoformat(self.processing_finished)
            self.runtime_seconds = round((end_dt - start_dt).total_seconds(), 2)

    def to_dict(self) -> Dict[str, Any]:
        """Recursively serializes entire sample entity into a structured JSON/TSV dictionary."""
        return {
            "identity": self.identity.to_dict(),
            "inputs": self.inputs.to_dict(),
            "processed_reads": self.processed_reads.to_dict(),
            "fingerprints": [fp.to_dict() for fp in self.fingerprints],
            "metrics": self.metrics.to_dict(),
            "tool_metadata": self.tool_meta.to_dict(),
            "validation": self.validation.to_dict(),
            "timestamps": {
                "started": self.processing_started,
                "finished": self.processing_finished,
                "runtime_seconds": self.runtime_seconds,
            },
        }