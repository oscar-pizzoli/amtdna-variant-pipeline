"""
Block 6: Consensus sequence generation and intra-individual heteroplasmy calling.

Exposes observational pileup streaming, majority-rule consensus builders,
statistical strand-bias heteroplasmy callers, and pipeline orchestrators.
"""

from src.modules.consensus.consensus_builder import ConsensusBuilder
from src.modules.consensus.enums import (
    CallStatus,
    ConsensusCallReason,
    ConsensusCallStatus,
    ContaminationEstimateSource,
    FilterReason,
    GatingAction,
    InputRole,
    IUPACAmbiguity,
)
from src.modules.consensus.heteroplasmy_caller import HeteroplasmyCaller
from src.modules.consensus.models import (
    AlleleCounts,
    ConsensusPolicy,
    ConsensusQC,
    ConsensusSite,
    HeteroplasmyCall,
    InputResolution,
    PileupFilterMetrics,
    ReferenceInfo,
    RegionCoverage,
    SitePileup,
)
from src.modules.consensus.orchestrator import (
    ConsensusExecutionError,
    ConsensusOrchestrator,
    ContractValidationError,
    UnresolvedAuthenticityError,
)
from src.modules.consensus.pileup_engine import PileupEngine

__all__ = [
    # Engine & Processing Components
    "PileupEngine",
    "ConsensusBuilder",
    "HeteroplasmyCaller",
    "ConsensusOrchestrator",
    # Enums
    "InputRole",
    "GatingAction",
    "ContaminationEstimateSource",
    "CallStatus",
    "ConsensusCallStatus",
    "ConsensusCallReason",
    "FilterReason",
    "IUPACAmbiguity",
    # Data Models & Contracts
    "ReferenceInfo",
    "AlleleCounts",
    "PileupFilterMetrics",
    "SitePileup",
    "InputResolution",
    "ConsensusPolicy",
    "ConsensusSite",
    "HeteroplasmyCall",
    "RegionCoverage",
    "ConsensusQC",
    # Exceptions
    "ConsensusExecutionError",
    "ContractValidationError",
    "UnresolvedAuthenticityError",
]