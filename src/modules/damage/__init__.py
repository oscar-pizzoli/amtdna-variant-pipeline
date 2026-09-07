"""Block 4: aDNA Damage Authentication & Base Quality Rescaling Module."""

from src.modules.damage.damage_quantifier import (
    DamageQuantifierEngine,
    DamageQuantifierError,
)
from src.modules.damage.models import (
    AuthenticationVerdict,
    BayesianDamageParameters,
    DamageErrorCode,
    DamageInputBAMType,
    DamageModelProfile,
    DamageOutputFiles,
    DamagePolicy,
    DamageResult,
    DamageSignalLevel,
    DamageStatus,
    DamageWarningCode,
    FragmentLengthMetrics,
    LibraryTreatment,
    PMDMetrics,
    TerminalDamageFrequencies,
)
from src.modules.damage.orchestrator import (
    DamageOrchestrationError,
    DamageOrchestrator,
)
from src.modules.damage.pmd_engine import (
    PMDExecutionMetadata,
    PMDScoringEngine,
    PMDScoringError,
)
from src.modules.damage.qc_evaluator import DamageQCValidator
from src.modules.damage.rescale_engine import (
    DamageRescaleError,
    QualityRescaleEngine,
    RescaleExecutionMetadata,
)

__all__ = [
    # Orchestrator
    "DamageOrchestrator",
    "DamageOrchestrationError",
    # Engines
    "DamageQuantifierEngine",
    "DamageQuantifierError",
    "QualityRescaleEngine",
    "DamageRescaleError",
    "RescaleExecutionMetadata",
    "PMDScoringEngine",
    "PMDScoringError",
    "PMDExecutionMetadata",
    "DamageQCValidator",
    # Domain Models & Dataclasses
    "DamagePolicy",
    "DamageResult",
    "DamageModelProfile",
    "TerminalDamageFrequencies",
    "FragmentLengthMetrics",
    "BayesianDamageParameters",
    "PMDMetrics",
    "DamageOutputFiles",
    # Enums
    "LibraryTreatment",
    "DamageInputBAMType",
    "DamageSignalLevel",
    "AuthenticationVerdict",
    "DamageStatus",
    "DamageWarningCode",
    "DamageErrorCode",
]
