"""
Modular Pipeline Steps for Ancient & Modern mtDNA Variant Analysis.

Separated into 1 clear Python file per analytical step:
- reference.py: Block 1 - Reference Genome Validation, Indexing & Verification
- discovery.py: Block 2 - Input Discovery, Preprocessing & Adapter Trimming
- alignment.py: Block 3 - BWA Alignment, Sorting, Deduplication & Quality Filtering
- damage.py: Block 4 - Ancient DNA Damage Quantification & Rescaling
- contamination.py: Block 5 - Two-Tier Contamination Estimation & Gating (Schmutzi & FastScreen)
- consensus.py: Block 6 - Mitochondrial Consensus Calling & Heteroplasmy Detection
- haplogroup.py: Block 7 - Phylogenetic Haplogroup Assignment & Cross-Block Audit
- synthesis.py: Block 8 - Forensic Synthesis, Provenance Audit & Cohort Delivery
"""

from src.steps.reference import (
    ReferenceInfo,
    ReferenceManager,
    ReferenceStep,
    ReferenceValidationResult,
)
from src.steps.discovery import (
    DataType,
    DiscoveryStep,
    DiscoveredFile,
    InputFiles,
    PreprocessingEngine,
    PreprocessingMetrics,
    ProcessedReads,
    SampleDiscoveryEngine,
    SampleRecord,
    ValidationStatus,
)
from src.steps.alignment import (
    AlignmentFilterPolicy,
    AlignmentMetrics,
    AlignmentOrchestrator,
    AlignmentResult,
    AlignmentStatus,
    AlignmentStep,
    CoverageMetrics,
    DedupEngine,
    DedupMetrics,
    DedupMode,
)
from src.steps.damage import (
    AuthenticationVerdict,
    DamageOrchestrator,
    DamagePolicy,
    DamageProfile,
    DamageQuantifierEngine,
    DamageResult,
    DamageStep,
    LibraryTreatment,
    PMDMetrics,
    PMDScoringEngine,
    TerminalDamage,
)
from src.steps.contamination import (
    ContaminationBatchRunner,
    ContaminationCategory,
    ContaminationMethod,
    ContaminationOrchestrator,
    ContaminationPolicy,
    ContaminationResult,
    ContaminationStep,
    ResourceAllocation,
    SchmutziEngine,
)
from src.steps.consensus import (
    CallStatus,
    ConsensusBuilder,
    ConsensusOrchestrator,
    ConsensusPolicy,
    ConsensusQC,
    ConsensusStep,
    GatingAction,
    HeteroplasmyCall,
    HeteroplasmyCaller,
    SitePileup,
)
from src.steps.haplogroup import (
    Block7Contract,
    ContaminantAffinity,
    HaplogrepEngine,
    HaplogroupOrchestrator,
    HaplogroupStatus,
    HaplogroupStep,
    OverallValidationStatus,
    SequenceValidator,
)
from src.steps.synthesis import (
    CategorizedWarning,
    FinalVerificationReport,
    ProvenanceItem,
    ProvenanceStatus,
    SampleMasterRecord,
    SynthesisOrchestrator,
    SynthesisStep,
    WarningCategory,
    WarningSeverity,
)

__all__ = [
    # Block 1
    "ReferenceInfo",
    "ReferenceManager",
    "ReferenceStep",
    "ReferenceValidationResult",
    # Block 2
    "DataType",
    "DiscoveryStep",
    "DiscoveredFile",
    "InputFiles",
    "PreprocessingEngine",
    "PreprocessingMetrics",
    "ProcessedReads",
    "SampleDiscoveryEngine",
    "SampleRecord",
    "ValidationStatus",
    # Block 3
    "AlignmentFilterPolicy",
    "AlignmentMetrics",
    "AlignmentOrchestrator",
    "AlignmentResult",
    "AlignmentStatus",
    "AlignmentStep",
    "CoverageMetrics",
    "DedupEngine",
    "DedupMetrics",
    "DedupMode",
    # Block 4
    "AuthenticationVerdict",
    "DamageOrchestrator",
    "DamagePolicy",
    "DamageProfile",
    "DamageQuantifierEngine",
    "DamageResult",
    "DamageStep",
    "LibraryTreatment",
    "PMDMetrics",
    "PMDScoringEngine",
    "TerminalDamage",
    # Block 5
    "ContaminationBatchRunner",
    "ContaminationCategory",
    "ContaminationMethod",
    "ContaminationOrchestrator",
    "ContaminationPolicy",
    "ContaminationResult",
    "ContaminationStep",
    "ResourceAllocation",
    "SchmutziEngine",
    # Block 6
    "CallStatus",
    "ConsensusBuilder",
    "ConsensusOrchestrator",
    "ConsensusPolicy",
    "ConsensusQC",
    "ConsensusStep",
    "GatingAction",
    "HeteroplasmyCall",
    "HeteroplasmyCaller",
    "SitePileup",
    # Block 7
    "Block7Contract",
    "ContaminantAffinity",
    "HaplogrepEngine",
    "HaplogroupOrchestrator",
    "HaplogroupStatus",
    "HaplogroupStep",
    "OverallValidationStatus",
    "SequenceValidator",
    # Block 8
    "CategorizedWarning",
    "FinalVerificationReport",
    "ProvenanceItem",
    "ProvenanceStatus",
    "SampleMasterRecord",
    "SynthesisOrchestrator",
    "SynthesisStep",
    "WarningCategory",
    "WarningSeverity",
]
