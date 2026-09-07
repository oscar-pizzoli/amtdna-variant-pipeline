"""Block 2 Preprocessing Orchestrator: bridges input ingestion and PreprocessingEngine execution."""

from dataclasses import dataclass
import logging
from pathlib import Path
from typing import Optional

from src.modules.preprocessing.models import (
    DataType,
    InputFiles,
    PreprocessingMetrics,
    ProcessedReads,
    SampleIdentity,
    SampleRecord,
    ValidationResult,
    ValidationStatus,
)
from src.modules.preprocessing.trimmer import PreprocessingEngine

logger = logging.getLogger("mtDNA_Pipeline.PreprocessingOrchestrator")


@dataclass
class RawSampleInput:
    """Convenience descriptor for paired-end or single-end FASTQ sample input."""

    sample_id: str
    r1_path: Path
    r2_path: Optional[Path] = None
    library_id: Optional[str] = None


class PreprocessingOrchestrator:
    """Orchestrates FASTQ ingestion, validation, and adapter trimming/merging via PreprocessingEngine."""

    @classmethod
    def run_sample(
        cls,
        sample_input: RawSampleInput,
        output_dir: Path,
        log_dir: Path,
        min_length: int = 30,
        min_overlap: int = 11,
        collapse_reads: bool = True,
        trim_ns: bool = True,
        threads: int = 4,
        integrity_mode: str = "fast",
    ) -> SampleRecord:
        """Executes full sample preprocessing, metrics extraction, and output verification."""
        sample_id = sample_input.sample_id
        is_paired = sample_input.r2_path is not None

        # 1. Build domain-compliant SampleRecord
        sample_record = SampleRecord(
            identity=SampleIdentity(
                sample_id=sample_id,
                library_id=sample_input.library_id or f"LIB_{sample_id}",
            ),
            inputs=InputFiles(
                data_type=DataType.FASTQ,
                is_paired=is_paired,
                r1_path=Path(sample_input.r1_path).resolve(),
                r2_path=Path(sample_input.r2_path).resolve() if is_paired else None,
            ),
            validation=ValidationResult(status=ValidationStatus.UNTESTED),
        )

        sample_out = output_dir / sample_id
        sample_log = log_dir / sample_id

        # 2. Delegate execution to the robust PreprocessingEngine
        PreprocessingEngine.run_adapterremoval(
            sample=sample_record,
            output_dir=sample_out,
            log_dir=sample_log,
            min_read_length=min_length,
            min_merge_overlap=min_overlap,
            collapse_reads=collapse_reads,
            trim_ns=trim_ns,
            threads=threads,
            integrity_mode=integrity_mode,
        )

        return sample_record