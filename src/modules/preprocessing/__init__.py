"""Preprocessing and Ingestion Engine package."""

from src.modules.preprocessing.bam_validator import (
    BAMValidationEngine,
    BAMValidationResult,
)
from src.modules.preprocessing.discovery import (
    DiscoveredFile,
    FilenameTokenizer,
    PairingKey,
    SampleDiscoveryEngine,
    SequencingUnit,
)
from src.modules.preprocessing.fastq_validator import (
    FastqValidationError,
    FastqValidator,
)
from src.modules.preprocessing.integrity import (
    IntegrityEngine,
    IntegrityError,
)
from src.modules.preprocessing.models import (
    DataType,
    FileFingerprint,
    InputFiles,
    OutputProfile,
    PreprocessingMetrics,
    PreprocessingTool,
    ProcessedReads,
    SampleIdentity,
    SampleRecord,
    ToolMetadata,
    ValidationResult,
    ValidationStatus,
)
from src.modules.preprocessing.trimmer import PreprocessingEngine
from src.modules.preprocessing.orchestrator import (
    PreprocessingOrchestrator,
    RawSampleInput,
)

__all__ = [
    # Data Models & Enums
    "DataType",
    "ValidationStatus",
    "PreprocessingTool",
    "OutputProfile",
    "SampleIdentity",
    "FileFingerprint",
    "InputFiles",
    "ProcessedReads",
    "ToolMetadata",
    "PreprocessingMetrics",
    "ValidationResult",
    "SampleRecord",
    # Intermediate Discovery Structures
    "PairingKey",
    "DiscoveredFile",
    "SequencingUnit",
    "FilenameTokenizer",
    # Engines & Validators
    "IntegrityEngine",
    "IntegrityError",
    "FastqValidator",
    "FastqValidationError",
    "BAMValidationEngine",
    "BAMValidationResult",
    "SampleDiscoveryEngine",
    "PreprocessingEngine",
    "PreprocessingOrchestrator",
    "RawSampleInput",
]