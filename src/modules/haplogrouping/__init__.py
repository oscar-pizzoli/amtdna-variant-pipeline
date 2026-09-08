"""
Block 7: Haplogrouping & Phylogenetic Validation Module.
Integrates sequence validation, dual-stream HaploGrep3 classification,
cross-block consistency audits, and formal contract emission.
"""

from src.modules.haplogrouping.cross_block_auditor import (
    CrossBlockAuditor,
    PhylotreeGraph,
)
from src.modules.haplogrouping.enums import (
    AuditCode,
    AuditSeverity,
    ContaminantAffinity,
    HaplogroupStatus,
    OverallValidationStatus,
    PrivateMutationDamageClass,
    SequenceType,
)
from src.modules.haplogrouping.haplogrep_engine import (
    HaplogrepEngine,
    HaplogrepExecutionError,
)
from src.modules.haplogrouping.models import (
    AuditFlag,
    Block7Contract,
    Block7JSONEncoder,
    CrossBlockAuditResult,
    HaplogrepHit,
    HaplogrepResult,
    HaplogrepRunConfig,
    PrivateMutationAssessment,
    RawHaplogrepResult,
    SequenceIntegrityMetrics,
    ValidationResult,
    to_json_serializable,
)
from src.modules.haplogrouping.orchestrator import HaplogroupOrchestrator
from src.modules.haplogrouping.sequence_validator import SequenceValidator

__all__ = [
    # Orchestrator
    "HaplogroupOrchestrator",
    # Engines & Validators
    "HaplogrepEngine",
    "HaplogrepExecutionError",
    "SequenceValidator",
    "CrossBlockAuditor",
    "PhylotreeGraph",
    # Enums
    "SequenceType",
    "HaplogroupStatus",
    "ContaminantAffinity",
    "AuditSeverity",
    "OverallValidationStatus",
    "PrivateMutationDamageClass",
    "AuditCode",
    # Contracts & Models
    "HaplogrepRunConfig",
    "HaplogrepHit",
    "RawHaplogrepResult",
    "PrivateMutationAssessment",
    "AuditFlag",
    "SequenceIntegrityMetrics",
    "CrossBlockAuditResult",
    "HaplogroupResult",
    "ValidationResult",
    "Block7Contract",
    "Block7JSONEncoder",
    "to_json_serializable",
]