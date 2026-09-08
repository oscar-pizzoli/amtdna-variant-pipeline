"""
Block 7 Models: Type-safe contracts, run configurations, raw tool outputs,
multi-candidate hit containers, cross-block audit models, and validation contracts.
"""

import json
from dataclasses import asdict, dataclass, field, is_dataclass
from enum import Enum
from pathlib import Path
from types import MappingProxyType
from typing import Any, Dict, List, Mapping, Optional

from src.modules.haplogrouping.enums import (
    AuditCode,
    AuditSeverity,
    ContaminantAffinity,
    HaplogroupStatus,
    OverallValidationStatus,
    PrivateMutationDamageClass,
    SequenceType,
)


def to_json_serializable(val: Any) -> Any:
    """Recursively convert Enums, Paths, Sets, MappingProxy, and Dataclasses into JSON-safe objects."""
    if isinstance(val, Enum):
        return val.value
    if isinstance(val, Path):
        return str(val)
    if isinstance(val, (set, frozenset)):
        return sorted([to_json_serializable(x) for x in val], key=lambda x: str(x))
    if isinstance(val, tuple):
        return [to_json_serializable(x) for x in val]
    if isinstance(val, (dict, MappingProxyType)):
        return {str(k): to_json_serializable(v) for k, v in val.items()}
    if is_dataclass(val):
        return {k: to_json_serializable(v) for k, v in asdict(val).items()}
    if isinstance(val, list):
        return [to_json_serializable(x) for x in val]
    return val


class Block7JSONEncoder(json.JSONEncoder):
    """Custom JSON encoder handling Enum, Path, set, MappingProxyType, and dataclasses."""
    def default(self, o: Any) -> Any:
        return to_json_serializable(o)


@dataclass(frozen=True)
class HaplogrepRunConfig:
    """Execution parameters and provenance contract for HaploGrep3."""
    tool: str
    tool_version: str
    tree_id: str
    tree_version: str
    metric: str
    min_quality_score: float
    cross_check_vcf: bool
    timeout_seconds: int = 300
    retain_raw_output: bool = True


@dataclass(frozen=True)
class HaplogrepHit:
    """Represents an individual haplogroup classification hit candidate."""
    haplogroup: str
    quality_score: float


@dataclass(frozen=True)
class RawHaplogrepResult:
    """Parsed classification result from a single HaploGrep3 execution."""
    sample_id: str
    input_type: str
    haplogroup: str
    quality_score: float
    second_best_quality: Optional[float]
    quality_margin: Optional[float]
    co_optimal: bool
    top_hits: List[HaplogrepHit]
    expected_mutations: List[str]
    found_mutations: List[str]
    missing_mutations: List[str]
    private_mutations: List[str]
    mutation_parse_warnings: List[str]
    haplogrep_version: str
    tree_id: str
    tree_version: str
    metric: str
    raw_output_path: Optional[str] = None


@dataclass(frozen=True)
class PrivateMutationAssessment:
    """
    Evaluates private variant chemistry against authenticity triage.
    Spatial position within reads is audited in downstream BAM stages.
    """
    mutation: str
    position: int
    ref_base: str
    alt_base: str
    classification: PrivateMutationDamageClass
    evidence_source: str
    rationale: str


@dataclass(frozen=True)
class AuditFlag:
    """Immutable audit flag observation."""
    code: AuditCode
    severity: AuditSeverity
    message: str
    context: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.context, MappingProxyType):
            object.__setattr__(self, "context", MappingProxyType(dict(self.context)))

    def to_dict(self) -> Dict[str, Any]:
        """Convert flag to fully JSON-safe dictionary preserving context."""
        return {
            "code": self.code.value,
            "severity": self.severity.value,
            "message": self.message,
            "context": to_json_serializable(self.context),
        }


@dataclass(frozen=True)
class SequenceIntegrityMetrics:
    """Pre-flight sequence audit evaluating reference invariants, alphabet, and provenance."""
    sample_id: str
    sequence_type: SequenceType
    sequence_length: Optional[int]
    expected_length: int
    n_count: Optional[int]
    n_fraction: Optional[float]
    max_n_fraction: float
    canonical_fraction: Optional[float]
    iupac_fraction: Optional[float]
    normalization_applied: bool
    observed_sha256: str
    expected_sha256: Optional[str]
    is_length_valid: bool
    is_n_fraction_valid: bool
    is_provenance_valid: bool
    reference_compatible: bool
    has_invalid_bases: bool
    issues: List[str] = field(default_factory=list)

    @property
    def is_valid(self) -> bool:
        return (
            self.is_length_valid
            and self.is_n_fraction_valid
            and self.is_provenance_valid
            and self.reference_compatible
            and not self.has_invalid_bases
        )


@dataclass(frozen=True)
class CrossBlockAuditResult:
    """Multivariate consistency audit bridging Blocks 4, 5, 6, and 7."""
    overall_status: OverallValidationStatus
    fasta_vcf_concordant: bool
    fasta_haplogroup: str
    vcf_haplogroup: Optional[str]
    resolved_haplogroup: str
    quality_score: float
    quality_margin: Optional[float]
    contaminant_affinity: ContaminantAffinity
    contaminant_haplogroups: List[str]
    contamination_fraction: Optional[float]
    subclade_shift_detected: bool = False
    subclade_shift_supported: Optional[bool] = None
    private_mutations_assessment: List[PrivateMutationAssessment] = field(default_factory=list)
    flags: List[AuditFlag] = field(default_factory=list)

    @property
    def damage_suspect_mutations(self) -> List[str]:
        return [
            p.mutation
            for p in self.private_mutations_assessment
            if p.classification == PrivateMutationDamageClass.TRANSITION_REQUIRES_SPATIAL_AUDIT
        ]

    @property
    def transversion_mutations(self) -> List[str]:
        return [
            p.mutation
            for p in self.private_mutations_assessment
            if p.classification == PrivateMutationDamageClass.TRANSVERSION_NON_DAMAGE
        ]


@dataclass(frozen=True)
class HaplogroupResult:
    """Phylogenetic classification output retaining primary and candidate hits."""
    sample_id: str
    primary_haplogroup: str
    quality_score: float
    status: HaplogroupStatus
    top_hits: List[HaplogrepHit]
    vcf_haplogroup: Optional[str]
    vcf_quality_score: Optional[float]
    run_config: HaplogrepRunConfig
    fasta_details: RawHaplogrepResult
    vcf_details: Optional[RawHaplogrepResult] = None


@dataclass(frozen=True)
class ValidationResult:
    """Comprehensive biological and structural validation verdict."""
    sample_id: str
    overall_status: OverallValidationStatus
    sequence_integrity: SequenceIntegrityMetrics
    haplogroup_acceptable: bool
    cross_block_audit: CrossBlockAuditResult
    flags: List[AuditFlag] = field(default_factory=list)


@dataclass(frozen=True)
class Block7Contract:
    """Final immutable output contract schema for sample.block7_contract.json."""
    sample_id: str
    primary_haplogroup: str
    resolved_haplogroup: str
    quality_score: float
    quality_margin: Optional[float]
    haplogroup_status: HaplogroupStatus
    overall_validation_status: OverallValidationStatus
    tree_id: str
    tree_version: str
    metric: str
    contaminant_affinity: ContaminantAffinity
    contamination_fraction: Optional[float]
    vcf_haplogroup: Optional[str]
    vcf_concordant: Optional[bool]
    subclade_shift_detected: bool
    subclade_shift_supported: Optional[bool]
    fasta_sha256: str
    flags: List[Dict[str, Any]]
    execution_time_seconds: float

    def to_dict(self) -> Dict[str, Any]:
        """Generic serialization avoiding manual enum drift."""
        return to_json_serializable(asdict(self))