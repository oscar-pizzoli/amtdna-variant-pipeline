"""
Data structures, observational contracts, and QC containers for Block 6.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional

from src.modules.consensus.enums import (
    CallStatus,
    ConsensusCallReason,
    ConsensusCallStatus,
    ContaminationEstimateSource,
    FilterReason,
    GatingAction,
    InputRole,
)


class LibraryTreatment(str, Enum):
    """Molecular library preparation enzymatic treatment."""
    NON_UDG = "NON_UDG"
    UDG_HALF = "UDG_HALF"
    UDG_FULL = "UDG_FULL"


@dataclass(frozen=True)
class ReferenceInfo:
    """Reference genome contig provenance and verification dimensions."""
    contig: str
    length: int
    sequence_sha256: str


@dataclass(frozen=True)
class AlleleCounts:
    """Type-safe counts for canonical nucleotides and non-canonical events."""
    A: int = 0
    C: int = 0
    G: int = 0
    T: int = 0
    N: int = 0
    other: int = 0

    @property
    def canonical_total(self) -> int:
        """Total number of standard nucleotide reads (A, C, G, T)."""
        return self.A + self.C + self.G + self.T

    @property
    def total(self) -> int:
        """Sum of all observed reads including N and non-canonical symbols."""
        return self.canonical_total + self.N + self.other

    def get(self, base: str) -> int:
        """Retrieve count for a specific nucleotide symbol."""
        clean_base = base.upper()
        if clean_base in {"A", "C", "G", "T", "N"}:
            return getattr(self, clean_base)
        return self.other


@dataclass(frozen=True)
class PileupFilterMetrics:
    """Granular exclusion counters tracking read-level quality, alignment, flags, and overlap."""
    low_bq: int = 0
    low_mq: int = 0
    duplicates: int = 0
    unmapped: int = 0
    qc_fail: int = 0
    secondary: int = 0
    supplementary: int = 0
    deletions: int = 0
    refskips: int = 0
    no_query_position: int = 0
    noncanonical_bases: int = 0
    overlapping_mate: int = 0
    excluded_reads: int = 0


@dataclass(frozen=True)
class SitePileup:
    """Observational record for a single genomic locus capturing allelic and spatial metrics."""
    position: int
    reference_base: str
    depth_pileup_entries: int
    depth_filtered: int
    depth_callable: int
    counts: AlleleCounts
    forward_counts: AlleleCounts
    reverse_counts: AlleleCounts
    mean_bq: float
    mean_mq: float
    filter_metrics: PileupFilterMetrics
    allele_distances: Dict[str, List[int]] = field(default_factory=dict)


@dataclass(frozen=True)
class HeteroplasmySiteDetail:
    """Contract-level spatial damage authentication container exported to downstream blocks."""
    ref: str
    alt: str
    vaf: float
    internal_support_count: int
    terminal_support_count: int
    spatial_damage_suspect: bool


@dataclass(frozen=True)
class InputResolution:
    """Validated input specification consumed from upstream Block 4 and 5 contracts."""
    sample_id: str
    bam_path: str
    bam_sha256: str
    input_role: InputRole
    reference_fasta_path: str
    reference_info: ReferenceInfo
    damage_model_id: Optional[str]
    gating_action: GatingAction
    contamination_fraction: Optional[float]
    contamination_estimate_source: ContaminationEstimateSource
    library_treatment: LibraryTreatment = LibraryTreatment.NON_UDG

    def validate_for_execution(self) -> None:
        """Ensure the input resolution contract is valid for downstream processing."""
        if self.input_role == InputRole.UNKNOWN:
            raise ValueError(
                f"Sample '{self.sample_id}' has unresolved InputRole.UNKNOWN; execution blocked."
            )
        if self.gating_action == GatingAction.REJECT_CONTAMINATED:
            raise ValueError(
                f"Sample '{self.sample_id}' is flagged as REJECT_CONTAMINATED; execution blocked."
            )


@dataclass(frozen=True)
class ConsensusPolicy:
    """Complete, immutable parameter configuration governing Block 6 decisions."""
    policy_id: str = "consensus_strict_v1"
    policy_version: str = "1.0.0"
    min_base_quality: int = 30
    min_mapping_quality: int = 30
    min_consensus_depth: int = 10
    major_allele_fraction: float = 0.80
    ambiguity_min_vaf: float = 0.10
    heteroplasmy_min_depth: int = 15
    heteroplasmy_min_alt_reads: int = 3
    heteroplasmy_min_vaf: float = 0.10
    heteroplasmy_max_vaf: float = 0.80
    strand_bias_alpha: float = 0.05
    strand_bias_min_strand_reads: int = 1
    emit_iupac_fasta: bool = True
    mask_char: str = "N"
    restricted_heteroplasmy_min_depth: int = 20
    restricted_heteroplasmy_min_alt_reads: int = 5
    restricted_heteroplasmy_min_vaf: float = 0.25
    restricted_heteroplasmy_max_vaf: float = 0.75

    # Spatial Damage Authentication
    min_internal_fragments: int = 2
    spatial_end_distance_threshold: int = 5
    spatial_thresholds: Dict[str, int] = field(
        default_factory=lambda: {"NON_UDG": 5, "UDG_HALF": 2, "UDG_FULL": 0}
    )

    def __post_init__(self) -> None:
        if self.mask_char != "N":
            raise ValueError(f"ConsensusPolicy requires mask_char == 'N', observed '{self.mask_char}'")
        if not (0.50 <= self.major_allele_fraction <= 1.00):
            raise ValueError("major_allele_fraction must be in range [0.50, 1.00]")
        if not (0.0 <= self.ambiguity_min_vaf <= 0.50):
            raise ValueError("ambiguity_min_vaf must be in range [0.00, 0.50]")


@dataclass(frozen=True)
class ConsensusSite:
    """Consensus decision for a single site across canonical and IUPAC streams."""
    position: int
    reference_base: str
    canonical_base: str
    iupac_base: str
    callable_depth: int
    major_allele: str
    major_vaf: float
    call_status: ConsensusCallStatus
    call_reason: ConsensusCallReason

    @property
    def depth(self) -> int:
        """Backward-compatible alias for callable_depth."""
        return self.callable_depth


@dataclass(frozen=True)
class HeteroplasmyCall:
    """Evaluated heteroplasmy observation retaining full strand, spatial QC, and provenance."""
    position: int
    reference_base: str
    major_allele: str
    minor_allele: str
    callable_depth: int
    alt_reads: int
    alt_vaf: float
    forward_ref: int
    reverse_ref: int
    forward_alt: int
    reverse_alt: int
    mean_bq: Optional[float]
    mean_mq: Optional[float]
    strand_bias_pvalue: float
    status: CallStatus
    filter_reason: FilterReason
    gating_action: GatingAction
    internal_support_count: int = 0
    terminal_support_count: int = 0
    spatial_damage_suspect: bool = False

    @property
    def depth(self) -> int:
        """Backward-compatible alias for callable_depth."""
        return self.callable_depth


@dataclass(frozen=True)
class RegionCoverage:
    """Summary coverage metric container for specific intervals (e.g. HVS-I, HVS-II)."""
    name: str
    start: int
    end: int
    mean_depth: float
    covered_fraction: float


@dataclass
class ConsensusQC:
    """Global execution QC metrics, audit trail, and cross-block spatial contracts."""
    sample_id: str
    input_bam_path: str
    input_bam_sha256: str
    reference_sha256: str
    reference_length: int
    mean_depth: float
    covered_positions_fraction: float
    n_fraction_canonical: float
    iupac_fraction: float
    heteroplasmy_sites_total: int
    heteroplasmy_sites_pass: int
    heteroplasmy_sites_filtered: int
    regions: List[RegionCoverage]
    gating_action: GatingAction
    software_version: str
    policy_version: str
    execution_time_seconds: float
    heteroplasmy_details: Dict[str, HeteroplasmySiteDetail] = field(default_factory=dict)