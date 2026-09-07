"""Contamination module public interface."""

from src.modules.contamination.models import (
    ConsensusGatingAction,
    ContaminationAuditMetadata,
    ContaminationCategory,
    ContaminationErrorCode,
    ContaminationEstimate,
    ContaminationInputContext,
    ContaminationMethod,
    ContaminationOutputFiles,
    ContaminationPolicy,
    ContaminationResult,
    ContaminationStatus,
    ContaminationWarningCode,
    EstimationConfidence,
    FastScreeningPolicy,
    FastScreenResult,
    GatingRecommendation,
    ModelQualityDiagnostics,
    SchmutziPolicy,
    ScreeningDecision,
    SubsamplingPolicy,
)
from src.modules.contamination.fast_screen import (
    FastContaminationScreen,
    FastScreenExecutionError,
)
from src.modules.contamination.schmutzi_engine import (
    SchmutziEngine,
    SchmutziExecutionError,
    SchmutziExecutionMetadata,
    SchmutziOutput,
)
from src.modules.contamination.dissent_engine import DissentEngine
from src.modules.contamination.qc_evaluator import ContaminationQCValidator
from src.modules.contamination.orchestrator import (
    ContaminationOrchestrator,
    ContaminationOrchestrationError,
)
from src.modules.contamination.batch_runner import (
    BatchSampleRecord,
    ContaminationBatchRunner,
    ResourceAllocation,
)

__all__ = [
    "ConsensusGatingAction",
    "ContaminationAuditMetadata",
    "ContaminationCategory",
    "ContaminationErrorCode",
    "ContaminationEstimate",
    "ContaminationInputContext",
    "ContaminationMethod",
    "ContaminationOutputFiles",
    "ContaminationPolicy",
    "ContaminationResult",
    "ContaminationStatus",
    "ContaminationWarningCode",
    "EstimationConfidence",
    "FastScreeningPolicy",
    "FastScreenResult",
    "GatingRecommendation",
    "ModelQualityDiagnostics",
    "SchmutziPolicy",
    "ScreeningDecision",
    "SubsamplingPolicy",
    "FastContaminationScreen",
    "FastScreenExecutionError",
    "SchmutziEngine",
    "SchmutziExecutionError",
    "SchmutziExecutionMetadata",
    "SchmutziOutput",
    "DissentEngine",
    "ContaminationQCValidator",
    "ContaminationOrchestrator",
    "ContaminationOrchestrationError",
    "ContaminationBatchRunner",
    "BatchSampleRecord",
    "ResourceAllocation",
]