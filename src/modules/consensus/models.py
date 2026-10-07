"""
Data structures, observational contracts, and QC containers for Block 6.
Extended with structural InDel tracking, homopolymer slippage definitions,
excluded blacklist intervals, strict input sanitization, and cryptographic provenance contracts.
"""

from dataclasses import dataclass, field
from enum import Enum
import re
from typing import Dict, List, Optional, Tuple

from src.modules.consensus.enums import (
    CallStatus,
    ConsensusCallReason,
    ConsensusCallStatus,
    ContaminationEstimateSource,
    FilterReason,
    GatingAction,
    InputRole,
)
from src.modules.synthesis.coordinate_normalizer import MutationType


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

    def __post_init__(self) -> None:
        if self.length <= 0:
            raise ValueError(f"ReferenceInfo length must be strictly positive, got {self.length}")
        if not re.fullmatch(r"[0-9a-fA-F]{64}", self.sequence_sha256):
            raise ValueError(
                f"ReferenceInfo sequence_sha256 must be a valid 64-char hex digest, got '{self.sequence_sha256}'"
            )


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
class IndelObservation:
    """
    Self-contained observational record of an insertion or deletion anchored at a locus.
    Carries spatial fragment boundary distributions and independent depth denominators.
    """
    event_type: MutationType
    anchor_position: int
    sequence: str
    length: int
    anchor_depth: int
    fwd_count: int = 0
    rev_count: int = 0
    read_end_distances: Tuple[int, ...] = field(default_factory=tuple)

    @property
    def total_count(self) -> int:
        """Total supporting molecules for this specific InDel event."""
        return self.fwd_count + self.rev_count

    @property
    def vaf(self) -> float:
        """Variant Allele Frequency computed over the local indel anchor depth."""
        if self.anchor_depth <= 0:
            return 0.0
        return self.total_count / self.anchor_depth


@dataclass(frozen=True)
class PileupFilterMetrics:
    """
    Granular exclusion counters strictly distinguishing filtered errors,
    structural genomic events, and deduplicated physical observations.
    """
    low_bq: int = 0
    low_mq: int = 0
    duplicates: int = 0
    unmapped: int = 0
    qc_fail: int = 0
    secondary: int = 0
    supplementary: int = 0
    refskips: int = 0
    no_query_position: int = 0
    noncanonical_bases: int = 0
    # Structural and molecular metrics (isolated from filtered errors)
    structural_deletion_spans: int = 0
    collapsed_mate_observations: int = 0
    raw_indel_observations: int = 0
    resolved_indel_observations: int = 0
    excluded_reads: int = 0


@dataclass(frozen=True)
class SitePileup:
    """
    Evidence-only container capturing single-base observations, structural InDels,
    and physical fragment terminus distances across deduplicated molecules.
    """
    position: int
    reference_base: str
    depth_pileup_entries: int
    depth_filtered: int
    depth_callable: int
    indel_anchor_depth: int
    counts: AlleleCounts
    forward_counts: AlleleCounts
    reverse_counts: AlleleCounts
    mean_bq: float
    mean_mq: float
    filter_metrics: PileupFilterMetrics
    read_end_distances: Dict[str, List[int]] = field(default_factory=dict)
    insertions: Dict[str, IndelObservation] = field(default_factory=dict)
    deletions: Dict[str, IndelObservation] = field(default_factory=dict)


@dataclass(frozen=True)
class HeteroplasmySiteDetail:
    """Self-contained spatial damage authentication container exported to downstream blocks."""
    position: int
    ref: str
    alt: str
    depth: int
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
        if self.contamination_fraction is not None:
            if not (0.0 <= self.contamination_fraction <= 1.0):
                raise ValueError(
                    f"Sample '{self.sample_id}' contamination_fraction out of bounds [0.0, 1.0]: {self.contamination_fraction}"
                )
        if not re.fullmatch(r"[0-9a-fA-F]{64}", self.bam_sha256):
            raise ValueError(
                f"Sample '{self.sample_id}' bam_sha256 must be a valid 64-char hex digest, got '{self.bam_sha256}'"
            )


@dataclass(frozen=True)
class HomopolymerRegion:
    """Region prone to polymerase slippage and alignment instability."""
    name: str
    start: int
    end: int
    homopolymer_base: str

    def __post_init__(self) -> None:
        if self.start > self.end or self.start <= 0:
            raise ValueError(f"Invalid HomopolymerRegion boundaries: [{self.start}, {self.end}]")


@dataclass(frozen=True)
class ConsensusPolicy:
    """Complete, immutable parameter configuration governing Block 6 decisions."""
    policy_id: str = "consensus_strict_v1"
    policy_version: str = "1.2.0"
    min_base_quality: int = 30
    min_mapping_quality: int = 30
    min_consensus_depth: int = 10
    major_allele_fraction: float = 0.80
    ambiguity_min_vaf: float = 0.10

    # SNV Heteroplasmy
    heteroplasmy_min_depth: int = 15
    heteroplasmy_min_alt_reads: int = 3
    heteroplasmy_min_vaf: float = 0.10
    heteroplasmy_max_vaf: float = 0.80
    strand_bias_alpha: float = 0.05
    strand_bias_min_strand_reads: int = 1

    # InDel Calling Parameters
    enable_indel_calling: bool = True
    indel_min_depth: int = 20
    indel_min_alt_reads: int = 4
    indel_min_vaf: float = 0.10
    indel_max_vaf: float = 0.80
    filter_homopolymer_slippage: bool = True

    emit_iupac_fasta: bool = True
    emit_subconsensus_vcf: bool = True
    mask_char: str = "N"

    restricted_heteroplasmy_min_depth: int = 20
    restricted_heteroplasmy_min_alt_reads: int = 5
    restricted_heteroplasmy_min_vaf: float = 0.25
    restricted_heteroplasmy_max_vaf: float = 0.75

    # Spatial Damage Authentication
    min_internal_fragments: int = 2
    spatial_end_distance_threshold: int = 5
    spatial_thresholds: Dict[LibraryTreatment, int] = field(
        default_factory=lambda: {
            LibraryTreatment.NON_UDG: 5,
            LibraryTreatment.UDG_HALF: 2,
            LibraryTreatment.UDG_FULL: 0,
        }
    )

    # Known slippage regions (HVS Poly-C tracts)
    slippage_regions: Tuple[HomopolymerRegion, ...] = field(
        default_factory=lambda: (
            HomopolymerRegion(name="HVS2_PolyC_303_315", start=303, end=315, homopolymer_base="C"),
            HomopolymerRegion(name="HVS3_PolyC_513_524", start=513, end=524, homopolymer_base="C"),
            HomopolymerRegion(name="HVS1_PolyC_16184_16193", start=16184, end=16193, homopolymer_base="C"),
        )
    )

    # Intervalli mascherati / blacklist per eteroplasmia (intervalli chiusi 1-based)
    excluded_intervals: Tuple[Tuple[int, int], ...] = field(
        default_factory=lambda: (
            (3105, 3108),      # Gap artificiale storico 3107N della rCRS
            (303, 315),        # Tratto omopolimerico poli-C HVS-II
            (513, 524),        # Tratto omopolimerico AC/C HVS-III
            (16184, 16193),    # Tratto omopolimerico poli-C HVS-I
        )
    )

    def is_position_excluded(self, pos: int) -> bool:
        """Restituisce True se la coordinata 1-based ricade in un intervallo escluso/blacklisted."""
        for start, end in self.excluded_intervals:
            if start <= pos <= end:
                return True
        return False

    def get_slippage_region(self, pos: int) -> Optional[HomopolymerRegion]:
        """Return the matching HomopolymerRegion if the locus falls into a tract, else None."""
        for r in self.slippage_regions:
            if r.start <= pos <= r.end:
                return r
        return None

    def __post_init__(self) -> None:
        if self.mask_char != "N":
            raise ValueError(f"ConsensusPolicy requires mask_char == 'N', observed '{self.mask_char}'")
        if not (0.50 <= self.major_allele_fraction <= 1.00):
            raise ValueError(f"major_allele_fraction must be in [0.50, 1.00], got {self.major_allele_fraction}")
        if not (0.0 <= self.ambiguity_min_vaf <= 0.50):
            raise ValueError(f"ambiguity_min_vaf must be in [0.00, 0.50], got {self.ambiguity_min_vaf}")

        for name, val in [
            ("min_consensus_depth", self.min_consensus_depth),
            ("heteroplasmy_min_depth", self.heteroplasmy_min_depth),
            ("restricted_heteroplasmy_min_depth", self.restricted_heteroplasmy_min_depth),
            ("indel_min_depth", self.indel_min_depth),
            ("heteroplasmy_min_alt_reads", self.heteroplasmy_min_alt_reads),
            ("restricted_heteroplasmy_min_alt_reads", self.restricted_heteroplasmy_min_alt_reads),
            ("indel_min_alt_reads", self.indel_min_alt_reads),
        ]:
            if val < 1:
                raise ValueError(f"{name} must be >= 1, got {val}")

        if not (0.0 <= self.heteroplasmy_min_vaf < self.heteroplasmy_max_vaf <= 1.0):
            raise ValueError(
                f"Invalid heteroplasmy VAF bounds: [{self.heteroplasmy_min_vaf}, {self.heteroplasmy_max_vaf}]"
            )
        if not (0.0 <= self.restricted_heteroplasmy_min_vaf < self.restricted_heteroplasmy_max_vaf <= 1.0):
            raise ValueError(
                f"Invalid restricted VAF bounds: [{self.restricted_heteroplasmy_min_vaf}, {self.restricted_heteroplasmy_max_vaf}]"
            )
        if not (0.0 <= self.indel_min_vaf < self.indel_max_vaf <= 1.0):
            raise ValueError(
                f"Invalid InDel VAF bounds: [{self.indel_min_vaf}, {self.indel_max_vaf}]"
            )

        if not (0.0 < self.strand_bias_alpha < 1.0):
            raise ValueError(f"strand_bias_alpha must be in (0.0, 1.0), got {self.strand_bias_alpha}")

        # Validazione limiti intervalli esclusi
        for start, end in self.excluded_intervals:
            if start > end or start <= 0:
                raise ValueError(f"Intervallo escluso non valido: [{start}, {end}]")


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
    variant_type: MutationType = MutationType.SNV
    internal_support_count: int = 0
    terminal_support_count: int = 0
    spatial_damage_suspect: bool = False
    is_slippage_suspect: bool = False
    slippage_region: Optional[str] = None

    @property
    def depth(self) -> int:
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
    reference_contig: str
    reference_sha256: str
    reference_length: int
    consensus_fasta_sha256: str
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