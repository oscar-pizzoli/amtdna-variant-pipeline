"""
Block 8 Models: Core dataclasses, formal schemas, multi-tier warning taxonomies,
and serialization primitives for holistic sample synthesis, sequencing QC, and cohort packaging.
"""

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from enum import Enum
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from src.modules.synthesis.coordinate_normalizer import MutationType


# ==============================================================================
# ENUMS
# ==============================================================================

class WarningCategory(str, Enum):
    """Five-tier categorical taxonomy separating biological warnings from technical anomalies."""
    QC_COVERAGE = "QC_COVERAGE"
    AUTHENTICITY_DAMAGE = "AUTHENTICITY_DAMAGE"
    CONTAMINATION = "CONTAMINATION"
    PHYLOGENY = "PHYLOGENY"
    PIPELINE_INTEGRITY = "PIPELINE_INTEGRITY"


class WarningSeverity(str, Enum):
    """Formal audit severity scale aligned across upstream and downstream blocks."""
    INFO = "INFO"
    WARNING = "WARNING"
    CRITICAL = "CRITICAL"


class MutationRole(str, Enum):
    """Phylogenetic assignment and biological role of an observed mutation."""
    DIAGNOSTIC = "DIAGNOSTIC"
    SUBCLADE_DIAGNOSTIC = "SUBCLADE_DIAGNOSTIC"
    PRIVATE = "PRIVATE"
    BACK_MUTATION = "BACK_MUTATION"
    UNKNOWN = "UNKNOWN"


class ProvenanceStatus(str, Enum):
    """Integrity and cryptographic verification state of a physical file."""
    VERIFIED = "VERIFIED"
    UNANCHORED = "UNANCHORED"
    MISMATCH = "MISMATCH"
    FILE_MISSING = "FILE_MISSING"
    IO_ERROR = "IO_ERROR"
    SEMANTIC_ERROR = "SEMANTIC_ERROR"
    GATED_OMITTED = "GATED_OMITTED"
    UNCHECKED = "UNCHECKED"


class SampleBiologicalType(str, Enum):
    """Authentication verdict derived from molecular degradation kinetic screening."""
    ANCIENT = "ANCIENT"
    MODERN = "MODERN"
    UNKNOWN = "UNKNOWN"


# ==============================================================================
# HELPER FORMATTERS
# ==============================================================================

def format_na(value: Any, float_format: Optional[str] = None) -> str:
    """Format any value, strictly guaranteeing emission of 'NA' when None, NaN, Inf, or empty."""
    if value is None:
        return "NA"
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return "NA"
        if float_format:
            return float_format.format(value)
        return f"{value:.4f}"
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    val_str = str(value).strip()
    return val_str if val_str and val_str.upper() not in {"NONE", "NULL", "NAN", "INF"} else "NA"


def format_pct_100(value: Optional[float], decimals: int = 2) -> str:
    """Format metrics already on a 0-100 scale (e.g. *_pct) with explicit percent sign."""
    if value is None or (isinstance(value, float) and (math.isnan(value) or math.isinf(value))):
        return "NA"
    return f"{value:.{decimals}f}%"


def format_pct_fraction(value: Optional[float], decimals: int = 2) -> str:
    """Format ratios and breadths on a 0.0-1.0 scale (e.g. *_fraction, *_rate) multiplying by 100."""
    if value is None or (isinstance(value, float) and (math.isnan(value) or math.isinf(value))):
        return "NA"
    return f"{value * 100.0:.{decimals}f}%"


# ==============================================================================
# CO-LOCATED MASTER TSV HEADERS
# ==============================================================================

SAMPLE_MASTER_TSV_HEADERS: List[str] = [
    # 1. General Identifiers
    "sample_id",
    "overall_status",
    "biological_type",
    # 2. Block 7 (Phylogeny & Lineage)
    "primary_haplogroup",
    "resolved_haplogroup",
    "quality_score",
    "quality_margin",
    "haplogroup_status",
    "subclade_shift_detected",
    "subclade_shift_supported",
    # 3. Sequencing QC & Trimming (Block 2)
    "raw_reads_total",
    "clean_reads_total",
    "retention_fraction",
    "collapsed_fraction",
    "adapter_trimmed_reads_pct",
    "q20_rate_post",
    "q30_rate_post",
    "gc_content_pct",
    "mean_read_length_pre",
    "mean_read_length_post",
    # 4. Block 3 (Alignment & Whole-mt Coverage)
    "mapped_reads_mt",
    "endogenous_dna_fraction",
    "endogenous_dna_source",
    "mean_depth",
    "coverage_breadth_1x",
    "coverage_breadth_3x",
    "coverage_breadth_5x",
    "duplication_rate",
    # 5. Mitochondrial Macro-Complexes Coverage (Block 6)
    "cov_dloop_mean",
    "cov_dloop_1x",
    "cov_rrna_mean",
    "cov_rrna_1x",
    "cov_cox_mean",
    "cov_cox_1x",
    "cov_nd_mean",
    "cov_nd_1x",
    "cov_cyb_mean",
    "cov_cyb_1x",
    # 6. Block 4 (Authenticity & Degradation Kinetics)
    "median_fragment_length",
    "pmd_fraction_ge_3",
    "c_to_t_5p_rate_pos1",
    "g_to_a_3p_rate_pos1",
    "rescaling_performed",
    # 7. Block 5 (Contamination Assessment)
    "contamination_fraction",
    "contamination_ci_95",
    "estimation_confidence",
    "recommended_action",
    # 8. Block 6 (Consensus & Heteroplasmy)
    "n_fraction_canonical",
    "heteroplasmy_sites_pass",
    # 9. Block 8 Synthesis & Audit
    "total_mutations_count",
    "total_warnings_count",
    "active_warning_codes",
]

MUTATION_TSV_HEADERS: List[str] = [
    "sample_id",
    "position",
    "ref",
    "alt",
    "variant_type",
    "mutation_notation",
    "stream_source",
    "role",
    "depth",
    "vaf",
    "spatial_damage_suspect",
    "internal_support_count",
    "terminal_support_count",
    "filter_status",
]


# ==============================================================================
# SUB-METRICS & COMPONENT RECORDS
# ==============================================================================

@dataclass
class InvariantViolation:
    """Diagnostic violation detected during cross-field invariant sanitization."""
    code: str
    severity: WarningSeverity
    message: str
    context: Dict[str, Any] = field(default_factory=dict)


@dataclass
class CategorizedWarning:
    """Atomic diagnostic flag classified within the 5-tier reporting taxonomy."""
    category: WarningCategory
    code: str
    severity: WarningSeverity
    message: str
    source_block: str = "UNKNOWN"
    context: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "category": self.category.value,
            "code": self.code,
            "severity": self.severity.value,
            "source_block": self.source_block,
            "message": self.message,
            "context": self.context,
        }


@dataclass
class ProvenanceItem:
    """Cryptographic chain-of-custody verification entry for an upstream artifact."""
    artifact_name: str
    file_path: str
    expected_sha256: str
    actual_sha256: Optional[str]
    file_size_bytes: Optional[int]
    status: ProvenanceStatus
    last_modified_utc: Optional[str] = None
    error_message: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "artifact_name": self.artifact_name,
            "file_path": self.file_path,
            "expected_sha256": self.expected_sha256,
            "actual_sha256": self.actual_sha256,
            "file_size_bytes": self.file_size_bytes,
            "status": self.status.value,
            "last_modified_utc": self.last_modified_utc,
            "error_message": self.error_message,
        }


@dataclass
class SemanticVerificationItem:
    """Detailed result of deep bioinformatic format and header inspection."""
    artifact_name: str
    file_path: str
    format_type: str  # BAM, FASTA, VCF, TSV
    is_valid: bool
    details: Dict[str, Any] = field(default_factory=dict)
    error_message: Optional[str] = None


@dataclass
class LineageVerificationItem:
    """Pairwise cross-block contract provenance identity check."""
    source_block: str
    target_block: str
    entity_name: str
    source_hash: str
    target_hash: str
    is_consistent: bool
    message: Optional[str] = None


@dataclass
class FinalVerificationReport:
    """Holistic forensic report assessing chain-of-custody and data integrity."""
    sample_id: str
    run_id: str
    overall_status: str  # PASS, WARN, FAIL
    manifest_anchored: bool
    contract_checks: List[ProvenanceItem] = field(default_factory=list)
    artifact_checks: List[ProvenanceItem] = field(default_factory=list)
    semantic_checks: List[SemanticVerificationItem] = field(default_factory=list)
    lineage_checks: List[LineageVerificationItem] = field(default_factory=list)
    critical_failures: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    verification_timestamp: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat()
    )


@dataclass
class MutationEntry:
    """Relational row for the tidy-format whole-cohort variant catalog."""
    sample_id: str
    position: int
    ref: str
    alt: str
    variant_type: MutationType
    mutation_notation: str
    stream_source: str  # "MAJOR_CONSENSUS", "HETEROPLASMY", "PHYLOGENY_ONLY"
    role: MutationRole
    depth: Optional[int] = None
    vaf: Optional[float] = None
    spatial_damage_suspect: Optional[bool] = None
    internal_support_count: Optional[int] = None
    terminal_support_count: Optional[int] = None
    filter_status: str = "PASS"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "position": self.position,
            "ref": self.ref,
            "alt": self.alt,
            "variant_type": self.variant_type.value,
            "mutation_notation": self.mutation_notation,
            "stream_source": self.stream_source,
            "role": self.role.value,
            "depth": self.depth,
            "vaf": round(self.vaf, 4) if self.vaf is not None else None,
            "spatial_damage_suspect": self.spatial_damage_suspect,
            "internal_support_count": self.internal_support_count,
            "terminal_support_count": self.terminal_support_count,
            "filter_status": self.filter_status,
        }

    def to_tsv_row(self) -> List[str]:
        row = [
            self.sample_id,
            str(self.position),
            self.ref,
            self.alt,
            self.variant_type.value,
            self.mutation_notation,
            self.stream_source,
            self.role.value,
            format_na(self.depth),
            format_na(self.vaf, "{:.4f}"),
            format_na(self.spatial_damage_suspect),
            format_na(self.internal_support_count),
            format_na(self.terminal_support_count),
            self.filter_status,
        ]
        if len(row) != len(MUTATION_TSV_HEADERS):
            raise ValueError(
                f"Mutation row length ({len(row)}) does not match MUTATION_TSV_HEADERS ({len(MUTATION_TSV_HEADERS)})."
            )
        return row


# ==============================================================================
# BLOCK-SPECIFIC SUMMARY STRUCTURES
# ==============================================================================

@dataclass
class PreprocessingSummary:
    """Block 2 Formal Summary."""
    raw_reads_total: Optional[int] = None
    retained_reads_total: Optional[int] = None
    retention_fraction: Optional[float] = None
    collapsed_fraction: Optional[float] = None
    output_profile: Optional[str] = None


@dataclass
class SequencingQCSummary:
    """Detailed sequencing quality and adapter trimming metrics (Block 2)."""
    raw_reads_total: Optional[int] = None
    clean_reads_total: Optional[int] = None
    retention_fraction: Optional[float] = None
    collapsed_fraction: Optional[float] = None
    adapter_trimmed_reads_pct: Optional[float] = None  # Scale 0.0 - 100.0
    q20_rate_post: Optional[float] = None              # Scale 0.0 - 1.0
    q30_rate_post: Optional[float] = None              # Scale 0.0 - 1.0
    gc_content_pct: Optional[float] = None             # Scale 0.0 - 100.0
    mean_read_length_pre: Optional[float] = None
    mean_read_length_post: Optional[float] = None
    per_base_quality_cycles: List[Dict[str, Any]] = field(default_factory=list)
    read_length_distribution: Dict[int, int] = field(default_factory=dict)


@dataclass
class AlignmentSummary:
    """Block 3 Formal Summary."""
    mapped_reads_mt: Optional[int] = None
    endogenous_dna_fraction: Optional[float] = None
    endogenous_dna_fraction_source: Optional[str] = None  # "BLOCK3" or "DERIVED_BLOCK2+BLOCK3"
    mean_depth: Optional[float] = None
    coverage_breadth_1x: Optional[float] = None
    coverage_breadth_3x: Optional[float] = None
    coverage_breadth_5x: Optional[float] = None
    coverage_breadth_10x: Optional[float] = None
    duplication_rate: Optional[float] = None


@dataclass
class MacroComplexCoverageSummary:
    """Coverage depth and breadth across mitochondrial functional complexes (Block 6)."""
    cov_dloop_mean: Optional[float] = None
    cov_dloop_1x: Optional[float] = None
    cov_rrna_mean: Optional[float] = None
    cov_rrna_1x: Optional[float] = None
    cov_cox_mean: Optional[float] = None
    cov_cox_1x: Optional[float] = None
    cov_nd_mean: Optional[float] = None
    cov_nd_1x: Optional[float] = None
    cov_cyb_mean: Optional[float] = None
    cov_cyb_1x: Optional[float] = None


@dataclass
class AuthenticitySummary:
    """Block 4 Formal Summary."""
    biological_type: SampleBiologicalType = SampleBiologicalType.UNKNOWN
    authentication_status: str = "UNKNOWN"
    median_fragment_length: Optional[float] = None
    mean_fragment_length: Optional[float] = None
    pmd_fraction_ge_3: Optional[float] = None
    c_to_t_5p_rate_pos1: Optional[float] = None
    g_to_a_3p_rate_pos1: Optional[float] = None
    rescaling_performed: bool = False
    damage_model_id: str = "UNKNOWN"
    library_treatment: str = "NON_UDG"


@dataclass
class ContaminationSummary:
    """Block 5 Formal Summary."""
    contamination_fraction: Optional[float] = None
    ci_95_lower: Optional[float] = None
    ci_95_upper: Optional[float] = None
    estimation_confidence: str = "UNRELIABLE"
    recommended_action: str = "EXCLUDE_OR_MANUAL_REVIEW"
    estimate_source: str = "NONE"
    contaminant_haplogroups: List[str] = field(default_factory=list)


@dataclass
class ConsensusSummary:
    """Block 6 Formal Summary."""
    mean_consensus_depth: Optional[float] = None
    covered_fraction: Optional[float] = None
    n_fraction_canonical: Optional[float] = None
    iupac_fraction: Optional[float] = None
    heteroplasmy_sites_total: int = 0
    heteroplasmy_sites_pass: int = 0
    heteroplasmy_sites_filtered: int = 0
    gating_action_applied: str = "RESTRICT_HETEROPLASMY"
    consensus_fasta_sha256: str = ""


@dataclass
class HaplogroupSummary:
    """Block 7 Formal Summary."""
    primary_haplogroup: str = "UNKNOWN"
    resolved_haplogroup: str = "UNKNOWN"
    quality_score: Optional[float] = None
    quality_margin: Optional[float] = None
    haplogroup_status: str = "UNRESOLVED"
    subclade_shift_detected: bool = False
    subclade_shift_supported: Optional[bool] = None
    contaminant_affinity: str = "NOT_APPLICABLE"
    found_diagnostic_mutations: List[str] = field(default_factory=list)
    missing_diagnostic_mutations: List[str] = field(default_factory=list)


# ==============================================================================
# MASTER SAMPLE RECORD & DELIVERABLE MODELS
# ==============================================================================

@dataclass
class SampleMasterRecord:
    """Unified single source of truth for an analyzed specimen across Blocks 1 through 7."""
    sample_id: str
    overall_status: str  # "PASS", "WARN", "FAIL"
    synthesis_timestamp: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat()
    )

    preprocessing: PreprocessingSummary = field(default_factory=PreprocessingSummary)
    sequencing_qc: SequencingQCSummary = field(default_factory=SequencingQCSummary)
    alignment: AlignmentSummary = field(default_factory=AlignmentSummary)
    macro_complexes: MacroComplexCoverageSummary = field(default_factory=MacroComplexCoverageSummary)
    authenticity: AuthenticitySummary = field(default_factory=AuthenticitySummary)
    contamination: ContaminationSummary = field(default_factory=ContaminationSummary)
    consensus: ConsensusSummary = field(default_factory=ConsensusSummary)
    haplogroup: HaplogroupSummary = field(default_factory=HaplogroupSummary)

    warnings: List[CategorizedWarning] = field(default_factory=list)
    mutations: List[MutationEntry] = field(default_factory=list)
    provenance: List[ProvenanceItem] = field(default_factory=list)

    execution_times_seconds: Dict[str, float] = field(default_factory=dict)
    upstream_contracts_sha256: Dict[str, str] = field(default_factory=dict)
    upstream_composite_hash: str = ""

    def validate_invariants(self) -> List[InvariantViolation]:
        """
        Validate cross-field logical invariants, sanitize contradictory states,
        and emit typed InvariantViolation items prior to packaging.
        """
        violations: List[InvariantViolation] = []

        # 1. Subclade shift coherence
        if not self.haplogroup.subclade_shift_detected:
            if self.haplogroup.subclade_shift_supported is not None:
                self.haplogroup.subclade_shift_supported = None
                violations.append(
                    InvariantViolation(
                        code="UPSTREAM_CONTRACT_ANOMALY",
                        severity=WarningSeverity.WARNING,
                        message="Normalized subclade_shift_supported to None because subclade_shift_detected is False.",
                        context={"sanitized_field": "subclade_shift_supported"},
                    )
                )

        # 2. Contamination confidence interval bounds
        c = self.contamination
        if c.ci_95_lower is not None and c.ci_95_upper is not None:
            if c.ci_95_lower > c.ci_95_upper:
                violations.append(
                    InvariantViolation(
                        code="UPSTREAM_CONTRACT_ANOMALY",
                        severity=WarningSeverity.CRITICAL,
                        message=f"Inverted contamination 95% CI: [{c.ci_95_lower}, {c.ci_95_upper}]. Lower bound exceeds upper bound.",
                        context={"ci_95_lower": c.ci_95_lower, "ci_95_upper": c.ci_95_upper},
                    )
                )
            if c.contamination_fraction is not None:
                if c.contamination_fraction < (c.ci_95_lower - 0.005) or c.contamination_fraction > (c.ci_95_upper + 0.005):
                    violations.append(
                        InvariantViolation(
                            code="UPSTREAM_CONTRACT_ANOMALY",
                            severity=WarningSeverity.WARNING,
                            message=f"Point estimate {c.contamination_fraction:.4f} lies outside 95% CI [{c.ci_95_lower:.4f}, {c.ci_95_upper:.4f}].",
                            context={"contamination_fraction": c.contamination_fraction, "ci_95": [c.ci_95_lower, c.ci_95_upper]},
                        )
                    )

        # 3. Derive endogenous DNA fraction if missing
        if self.alignment.endogenous_dna_fraction is None:
            raw = self.preprocessing.raw_reads_total
            mapped = self.alignment.mapped_reads_mt
            if raw and mapped is not None and raw > 0:
                self.alignment.endogenous_dna_fraction = round(mapped / raw, 6)
                self.alignment.endogenous_dna_fraction_source = "DERIVED_BLOCK2+BLOCK3"

        # 4. Probability bounds check [0.0, 1.0] for fractional fields
        bounded_fractions = [
            ("retention_fraction", self.preprocessing.retention_fraction),
            ("collapsed_fraction", self.preprocessing.collapsed_fraction),
            ("q20_rate_post", self.sequencing_qc.q20_rate_post),
            ("q30_rate_post", self.sequencing_qc.q30_rate_post),
            ("endogenous_dna_fraction", self.alignment.endogenous_dna_fraction),
            ("coverage_breadth_1x", self.alignment.coverage_breadth_1x),
            ("coverage_breadth_3x", self.alignment.coverage_breadth_3x),
            ("coverage_breadth_5x", self.alignment.coverage_breadth_5x),
            ("duplication_rate", self.alignment.duplication_rate),
            ("cov_dloop_1x", self.macro_complexes.cov_dloop_1x),
            ("cov_rrna_1x", self.macro_complexes.cov_rrna_1x),
            ("cov_cox_1x", self.macro_complexes.cov_cox_1x),
            ("cov_nd_1x", self.macro_complexes.cov_nd_1x),
            ("cov_cyb_1x", self.macro_complexes.cov_cyb_1x),
            ("pmd_fraction_ge_3", self.authenticity.pmd_fraction_ge_3),
            ("covered_fraction", self.consensus.covered_fraction),
            ("n_fraction_canonical", self.consensus.n_fraction_canonical),
            ("contamination_fraction", self.contamination.contamination_fraction),
        ]
        for name, val in bounded_fractions:
            if val is not None and not (-1e-6 <= val <= 1.0001):
                violations.append(
                    InvariantViolation(
                        code="UPSTREAM_CONTRACT_ANOMALY",
                        severity=WarningSeverity.CRITICAL,
                        message=f"Field '{name}' value {val} violates mathematical probability bounds [0.0, 1.0].",
                        context={"field": name, "value": val},
                    )
                )

        # 5. Percentage bounds check [0.0, 100.0] for fields with _pct suffix
        bounded_pcts = [
            ("adapter_trimmed_reads_pct", self.sequencing_qc.adapter_trimmed_reads_pct),
            ("gc_content_pct", self.sequencing_qc.gc_content_pct),
        ]
        for name, val in bounded_pcts:
            if val is not None and not (-1e-6 <= val <= 100.0001):
                violations.append(
                    InvariantViolation(
                        code="UPSTREAM_CONTRACT_ANOMALY",
                        severity=WarningSeverity.CRITICAL,
                        message=f"Field '{name}' value {val} violates percentage bounds [0.0, 100.0].",
                        context={"field": name, "value": val},
                    )
                )

        return violations

    def get_warning_codes(self, category: Optional[WarningCategory] = None) -> List[str]:
        """Retrieve deduplicated warning codes, optionally filtered by category."""
        if category is None:
            return sorted(list({w.code for w in self.warnings}))
        return sorted(list({w.code for w in self.warnings if w.category == category}))

    def get_warning_narrative(self) -> str:
        """Compose a natural-language executive narrative of all diagnostic findings."""
        if not self.warnings:
            return (
                f"Sample {self.sample_id} successfully passed all quality, authenticity, contamination, "
                f"and phylogenetic verification gates with zero anomalies detected."
            )

        sentences: List[str] = []
        crit_count = sum(1 for w in self.warnings if w.severity == WarningSeverity.CRITICAL)
        warn_count = sum(1 for w in self.warnings if w.severity == WarningSeverity.WARNING)

        if crit_count > 0:
            sentences.append(
                f"Sample {self.sample_id} failed final validation due to {crit_count} critical issue(s) "
                f"and {warn_count} warning(s)."
            )
        else:
            sentences.append(
                f"Sample {self.sample_id} passed validation with cautionary notices ({warn_count} warning(s))."
            )

        for cat, label in [
            (WarningCategory.PIPELINE_INTEGRITY, "Pipeline Integrity"),
            (WarningCategory.QC_COVERAGE, "Quality & Coverage"),
            (WarningCategory.AUTHENTICITY_DAMAGE, "Authenticity"),
            (WarningCategory.CONTAMINATION, "Contamination"),
            (WarningCategory.PHYLOGENY, "Phylogeny"),
        ]:
            cat_warns = [w for w in self.warnings if w.category == cat]
            if cat_warns:
                msgs = "; ".join(f"[{w.source_block}] {w.message}" for w in cat_warns)
                sentences.append(f"{label}: {msgs}")

        return " ".join(sentences)

    def to_dict(self) -> Dict[str, Any]:
        """Convert master record to structured dictionary with proper nested representations."""
        return {
            "sample_id": self.sample_id,
            "overall_status": self.overall_status,
            "synthesis_timestamp": self.synthesis_timestamp,
            "narrative": self.get_warning_narrative(),
            "preprocessing": asdict(self.preprocessing),
            "sequencing_qc": asdict(self.sequencing_qc),
            "alignment": asdict(self.alignment),
            "macro_complexes": asdict(self.macro_complexes),
            "authenticity": {
                **asdict(self.authenticity),
                "biological_type": self.authenticity.biological_type.value,
            },
            "contamination": asdict(self.contamination),
            "consensus": asdict(self.consensus),
            "haplogroup": asdict(self.haplogroup),
            "warning_counts": {
                cat.value: len(self.get_warning_codes(cat)) for cat in WarningCategory
            },
            "warnings": [w.to_dict() for w in self.warnings],
            "mutations_count": len(self.mutations),
            "mutations": [m.to_dict() for m in self.mutations],
            "provenance": [p.to_dict() for p in self.provenance],
            "execution_times_seconds": self.execution_times_seconds,
            "upstream_contracts_sha256": self.upstream_contracts_sha256,
            "upstream_composite_hash": self.upstream_composite_hash,
        }

    def to_summary_tsv_row(self) -> List[str]:
        """
        Emit flattened row strictly formatted and positional-parity verified
        against SAMPLE_MASTER_TSV_HEADERS (single source of truth for cohort exports).
        """
        active_warns = [w.code for w in self.warnings]
        warn_str = "; ".join(active_warns) if active_warns else "CLEAN"

        contam_str = format_na(self.contamination.contamination_fraction, "{:.4f}")
        ci_str = "NA"
        if self.contamination.ci_95_lower is not None and self.contamination.ci_95_upper is not None:
            ci_str = f"[{self.contamination.ci_95_lower:.4f}-{self.contamination.ci_95_upper:.4f}]"

        shift_supp_str = format_na(self.haplogroup.subclade_shift_supported)
        mc = self.macro_complexes
        sq = self.sequencing_qc

        row = [
            # 1. General Identifiers
            self.sample_id,
            self.overall_status,
            self.authenticity.biological_type.value,
            # 2. Block 7 (Phylogeny & Lineage)
            self.haplogroup.primary_haplogroup,
            self.haplogroup.resolved_haplogroup,
            format_na(self.haplogroup.quality_score, "{:.4f}"),
            format_na(self.haplogroup.quality_margin, "{:.4f}"),
            self.haplogroup.haplogroup_status,
            "TRUE" if self.haplogroup.subclade_shift_detected else "FALSE",
            shift_supp_str,
            # 3. Sequencing QC & Trimming (Block 2)
            format_na(sq.raw_reads_total if sq.raw_reads_total is not None else self.preprocessing.raw_reads_total),
            format_na(sq.clean_reads_total if sq.clean_reads_total is not None else self.preprocessing.retained_reads_total),
            format_na(self.preprocessing.retention_fraction, "{:.4f}"),
            format_na(self.preprocessing.collapsed_fraction, "{:.4f}"),
            format_pct_100(sq.adapter_trimmed_reads_pct),
            format_pct_fraction(sq.q20_rate_post),
            format_pct_fraction(sq.q30_rate_post),
            format_pct_100(sq.gc_content_pct),
            format_na(sq.mean_read_length_pre, "{:.1f}"),
            format_na(sq.mean_read_length_post, "{:.1f}"),
            # 4. Block 3 (Alignment & Whole-mt Coverage)
            format_na(self.alignment.mapped_reads_mt),
            format_na(self.alignment.endogenous_dna_fraction, "{:.4f}"),
            format_na(self.alignment.endogenous_dna_fraction_source),
            format_na(self.alignment.mean_depth, "{:.2f}"),
            format_na(self.alignment.coverage_breadth_1x, "{:.4f}"),
            format_na(self.alignment.coverage_breadth_3x, "{:.4f}"),
            format_na(self.alignment.coverage_breadth_5x, "{:.4f}"),
            format_na(self.alignment.duplication_rate, "{:.4f}"),
            # 5. Mitochondrial Macro-Complexes Coverage (Block 6)
            format_na(mc.cov_dloop_mean, "{:.1f}"),
            format_pct_fraction(mc.cov_dloop_1x),
            format_na(mc.cov_rrna_mean, "{:.1f}"),
            format_pct_fraction(mc.cov_rrna_1x),
            format_na(mc.cov_cox_mean, "{:.1f}"),
            format_pct_fraction(mc.cov_cox_1x),
            format_na(mc.cov_nd_mean, "{:.1f}"),
            format_pct_fraction(mc.cov_nd_1x),
            format_na(mc.cov_cyb_mean, "{:.1f}"),
            format_pct_fraction(mc.cov_cyb_1x),
            # 6. Block 4 (Authenticity & Degradation Kinetics)
            format_na(self.authenticity.median_fragment_length, "{:.1f}"),
            format_na(self.authenticity.pmd_fraction_ge_3, "{:.4f}"),
            format_na(self.authenticity.c_to_t_5p_rate_pos1, "{:.4f}"),
            format_na(self.authenticity.g_to_a_3p_rate_pos1, "{:.4f}"),
            "TRUE" if self.authenticity.rescaling_performed else "FALSE",
            # 7. Block 5 (Contamination Assessment)
            contam_str,
            ci_str,
            self.contamination.estimation_confidence,
            self.contamination.recommended_action,
            # 8. Block 6 (Consensus & Heteroplasmy)
            format_na(self.consensus.n_fraction_canonical, "{:.4f}"),
            str(self.consensus.heteroplasmy_sites_pass),
            # 9. Block 8 Synthesis & Audit
            str(len(self.mutations)),
            str(len(self.warnings)),
            warn_str,
        ]

        if len(row) != len(SAMPLE_MASTER_TSV_HEADERS):
            raise ValueError(
                f"Master TSV row length ({len(row)}) does not match SAMPLE_MASTER_TSV_HEADERS ({len(SAMPLE_MASTER_TSV_HEADERS)})."
            )
        return row


@dataclass
class PackagingConfig:
    """Settings controlling output formats, directory layout, and export options."""
    reports_dir: Path
    emit_tsv_summary: bool = True
    emit_mutations_catalog: bool = True
    emit_json_master: bool = True
    emit_markdown_sheets: bool = True
    enforce_strict_provenance: bool = True
    clean_na_tsv: bool = True