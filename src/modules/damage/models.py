"""Domain models, enums, dataclasses, and validation schemas for Block 4 (Damage & Authentication)."""

from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional


class LibraryTreatment(str, Enum):
    """Experimental enzymatic library treatment protocol."""

    NONE = "NONE"                     # Non-UDG treated standard ancient DNA library
    UDG_PARTIAL = "UDG_PARTIAL"       # Partial UDG / USER treatment (damage retained only at terminal nucleotides)
    UDG_FULL = "UDG_FULL"             # Full UDG / USER treatment (uracils enzymatically repaired)
    SINGLE_STRANDED = "SINGLE_STRANDED"  # Single-stranded DNA library (ssDNA)
    UNKNOWN = "UNKNOWN"               # Unspecified or unrecorded library preparation


class DamageInputBAMType(str, Enum):
    """Choice of upstream BAM file fed into damage estimation engines."""

    FILTERED_DEDUP = "FILTERED_DEDUP"      # Deduplicated BAM (prevents PCR duplicate bias)
    FILTERED_MARKDUP = "FILTERED_MARKDUP"  # Duplicate-flagged BAM (0x400 retained for full-library sampling)
    RAW_SORTED = "RAW_SORTED"              # Pre-filtered coordinate-sorted BAM


class DamageSignalLevel(str, Enum):
    """Categorical assessment of observable empirical deamination signal."""

    STRONG_SIGNAL = "STRONG_SIGNAL"                # Clear elevated terminal transitions supported by sufficient observations
    MODERATE_SIGNAL = "MODERATE_SIGNAL"            # Measurable and consistent damage above background divergence
    LOW_SIGNAL = "LOW_SIGNAL"                      # Baseline or flat transition rates near sequencing error floor
    NOISY_OR_INSUFFICIENT = "NOISY_OR_INSUFFICIENT"  # Insufficient informative sites or poor model convergence


class AuthenticationVerdict(str, Enum):
    """Protocol-aware authenticity diagnosis for ancient DNA authentication."""

    AUTHENTIC_ANCIENT_SUPPORTED = "AUTHENTIC_ANCIENT_SUPPORTED"  # Observed damage profile strictly matches protocol expectation
    AUTHENTICITY_UNSUPPORTED = "AUTHENTICITY_UNSUPPORTED"        # Damage profile does not provide positive ancient evidence
    UDG_DEPLETED_EXPECTED = "UDG_DEPLETED_EXPECTED"              # Nominal absence of deamination on UDG_FULL library
    INCONCLUSIVE = "INCONCLUSIVE"                                # Insufficient depth or ambiguous statistical metrics


class DamageStatus(str, Enum):
    """Overall execution status of Block 4."""

    UNTESTED = "UNTESTED"
    PASS = "PASS"
    WARN = "WARN"
    FAIL = "FAIL"


class DamageWarningCode(str, Enum):
    """Diagnostic warning codes for Block 4 quality control."""

    LOW_DAMAGE_ON_NON_UDG = "LOW_DAMAGE_ON_NON_UDG"
    UNEXPECTED_HIGH_DAMAGE_ON_UDG = "UNEXPECTED_HIGH_DAMAGE_ON_UDG"
    LOW_READ_COUNT_FOR_DAMAGE_MODEL = "LOW_READ_COUNT_FOR_DAMAGE_MODEL"
    INSUFFICIENT_TERMINAL_COVERAGE = "INSUFFICIENT_TERMINAL_COVERAGE"
    UNKNOWN_LIBRARY_TREATMENT = "UNKNOWN_LIBRARY_TREATMENT"
    BAYESIAN_CONVERGENCE_FAILURE = "BAYESIAN_CONVERGENCE_FAILURE"
    MISSING_MD_TAGS = "MISSING_MD_TAGS"
    LOW_PMD_FRACTION = "LOW_PMD_FRACTION"
    ASYMMETRIC_TERMINAL_DAMAGE = "ASYMMETRIC_TERMINAL_DAMAGE"


class DamageErrorCode(str, Enum):
    """Fatal error codes halting downstream processing."""

    MAPDAMAGE_EXECUTION_FAILED = "MAPDAMAGE_EXECUTION_FAILED"
    PMDTOOLS_EXECUTION_FAILED = "PMDTOOLS_EXECUTION_FAILED"
    RESCALING_FAILED = "RESCALING_FAILED"
    CORRUPTED_DAMAGE_OUTPUT = "CORRUPTED_DAMAGE_OUTPUT"
    INPUT_BAM_NOT_FOUND = "INPUT_BAM_NOT_FOUND"


@dataclass
class TerminalDamageFrequencies:
    """Empirical nucleotide substitution rates and supporting observation counts at terminal positions."""

    # 5' End C->T Mismatch Rates (Fractions in [0.0, 1.0])
    c_to_t_5p_pos1: float = 0.0
    c_to_t_5p_pos2: float = 0.0
    c_to_t_5p_profile: List[float] = field(default_factory=list)
    max_5p_c_to_t: float = 0.0

    # 3' End G->A Mismatch Rates (Fractions in [0.0, 1.0])
    g_to_a_3p_pos1: float = 0.0
    g_to_a_3p_pos2: float = 0.0
    g_to_a_3p_profile: List[float] = field(default_factory=list)
    max_3p_g_to_a: float = 0.0

    # Statistical Denominators & Supporting Counts
    n_c_to_t_5p_pos1_obs: int = 0
    n_c_ref_5p_pos1_total: int = 0
    n_g_to_a_3p_pos1_obs: int = 0
    n_g_ref_3p_pos1_total: int = 0

    c_to_t_5p_observations: List[int] = field(default_factory=list)
    c_to_t_5p_denominators: List[int] = field(default_factory=list)
    g_to_a_3p_observations: List[int] = field(default_factory=list)
    g_to_a_3p_denominators: List[int] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Serializes terminal mismatch rates and counts to a dictionary."""
        return {
            "c_to_t_5p_pos1_fraction": round(self.c_to_t_5p_pos1, 4),
            "c_to_t_5p_pos2_fraction": round(self.c_to_t_5p_pos2, 4),
            "g_to_a_3p_pos1_fraction": round(self.g_to_a_3p_pos1, 4),
            "g_to_a_3p_pos2_fraction": round(self.g_to_a_3p_pos2, 4),
            "max_5p_c_to_t_fraction": round(self.max_5p_c_to_t, 4),
            "max_3p_g_to_a_fraction": round(self.max_3p_g_to_a, 4),
            "c_to_t_5p_profile": [round(x, 4) for x in self.c_to_t_5p_profile],
            "g_to_a_3p_profile": [round(x, 4) for x in self.g_to_a_3p_profile],
            "n_c_to_t_5p_pos1_obs": self.n_c_to_t_5p_pos1_obs,
            "n_c_ref_5p_pos1_total": self.n_c_ref_5p_pos1_total,
            "n_g_to_a_3p_pos1_obs": self.n_g_to_a_3p_pos1_obs,
            "n_g_ref_3p_pos1_total": self.n_g_ref_3p_pos1_total,
            "c_to_t_5p_observations": self.c_to_t_5p_observations,
            "c_to_t_5p_denominators": self.c_to_t_5p_denominators,
            "g_to_a_3p_observations": self.g_to_a_3p_observations,
            "g_to_a_3p_denominators": self.g_to_a_3p_denominators,
        }


@dataclass
class FragmentLengthMetrics:
    """Summary statistics of aligned fragment length distribution."""

    total_fragments_sampled: int = 0
    mean_length: float = 0.0
    median_length: float = 0.0
    p25_length: float = 0.0
    p75_length: float = 0.0
    p90_length: float = 0.0
    min_length: int = 0
    max_length: int = 0

    def to_dict(self) -> Dict[str, Any]:
        """Serializes length distribution statistics to a dictionary."""
        return {
            "total_fragments_sampled": self.total_fragments_sampled,
            "mean_length_bp": round(self.mean_length, 2),
            "median_length_bp": round(self.median_length, 2),
            "p25_length_bp": round(self.p25_length, 2),
            "p75_length_bp": round(self.p75_length, 2),
            "p90_length_bp": round(self.p90_length, 2),
            "min_length_bp": self.min_length,
            "max_length_bp": self.max_length,
        }


@dataclass
class BayesianDamageParameters:
    """Statistical parameters estimated via mapDamage2 statistical model."""

    lambda_param: Optional[float] = None
    delta_s: Optional[float] = None
    delta_d: Optional[float] = None
    theta: Optional[float] = None
    model_log_likelihood: Optional[float] = None
    model_converged: bool = False
    iterations: int = 0

    def to_dict(self) -> Dict[str, Any]:
        """Serializes Bayesian model parameters to a dictionary."""
        return {
            "lambda": round(self.lambda_param, 4) if self.lambda_param is not None else None,
            "delta_s": round(self.delta_s, 4) if self.delta_s is not None else None,
            "delta_d": round(self.delta_d, 4) if self.delta_d is not None else None,
            "theta": round(self.theta, 6) if self.theta is not None else None,
            "model_log_likelihood": round(self.model_log_likelihood, 4) if self.model_log_likelihood is not None else None,
            "model_converged": self.model_converged,
            "iterations": self.iterations,
        }


@dataclass
class PMDMetrics:
    """Read-level Post-Mortem Damage (PMD) score distribution and filtering statistics."""

    total_reads_scored: int = 0
    min_pmd_score: float = 0.0
    max_pmd_score: float = 0.0
    mean_pmd_score: float = 0.0
    median_pmd_score: float = 0.0
    p75_pmd_score: float = 0.0
    p90_pmd_score: float = 0.0
    fraction_pmd_ge_0: float = 0.0
    fraction_pmd_ge_3: float = 0.0
    fraction_pmd_ge_5: float = 0.0
    applied_threshold: Optional[float] = None
    reads_retained_post_pmd_filter: Optional[int] = None

    def to_dict(self) -> Dict[str, Any]:
        """Serializes PMD metrics to a dictionary."""
        return {
            "total_reads_scored": self.total_reads_scored,
            "min_pmd_score": round(self.min_pmd_score, 4),
            "max_pmd_score": round(self.max_pmd_score, 4),
            "mean_pmd_score": round(self.mean_pmd_score, 4),
            "median_pmd_score": round(self.median_pmd_score, 4),
            "p75_pmd_score": round(self.p75_pmd_score, 4),
            "p90_pmd_score": round(self.p90_pmd_score, 4),
            "fraction_pmd_ge_0": round(self.fraction_pmd_ge_0, 4),
            "fraction_pmd_ge_3": round(self.fraction_pmd_ge_3, 4),
            "fraction_pmd_ge_5": round(self.fraction_pmd_ge_5, 4),
            "applied_threshold": self.applied_threshold,
            "reads_retained_post_pmd_filter": self.reads_retained_post_pmd_filter,
        }


@dataclass
class DamageOutputFiles:
    """Physical filesystem paths of Block 4 deliverables."""

    misincorporation_txt: Optional[Path] = None
    length_distribution_txt: Optional[Path] = None
    damage_plot_pdf: Optional[Path] = None
    rescaled_bam: Optional[Path] = None
    rescaled_bai: Optional[Path] = None
    pmd_scores_tsv: Optional[Path] = None
    pmd_filtered_bam: Optional[Path] = None
    pmd_filtered_bai: Optional[Path] = None
    damage_model_json: Optional[Path] = None

    def to_dict(self) -> Dict[str, Any]:
        """Serializes output paths to a dictionary."""
        return {
            "misincorporation_txt": str(self.misincorporation_txt) if self.misincorporation_txt else None,
            "length_distribution_txt": str(self.length_distribution_txt) if self.length_distribution_txt else None,
            "damage_plot_pdf": str(self.damage_plot_pdf) if self.damage_plot_pdf else None,
            "rescaled_bam": str(self.rescaled_bam) if self.rescaled_bam else None,
            "rescaled_bai": str(self.rescaled_bai) if self.rescaled_bai else None,
            "pmd_scores_tsv": str(self.pmd_scores_tsv) if self.pmd_scores_tsv else None,
            "pmd_filtered_bam": str(self.pmd_filtered_bam) if self.pmd_filtered_bam else None,
            "pmd_filtered_bai": str(self.pmd_filtered_bai) if self.pmd_filtered_bai else None,
            "damage_model_json": str(self.damage_model_json) if self.damage_model_json else None,
        }


@dataclass
class DamagePolicy:
    """Configurable execution parameters, heuristic threshold limits, and quality criteria for Block 4."""

    input_bam_type: DamageInputBAMType = DamageInputBAMType.FILTERED_DEDUP
    library_treatment: LibraryTreatment = LibraryTreatment.NONE
    enable_rescaling: bool = True
    enable_pmd_scoring: bool = True
    enable_pmd_filtering: bool = False
    pmd_threshold: float = 3.0

    # Policy Thresholds (Configurable empirical defaults)
    min_c_to_t_5p_strong_threshold: float = 0.10
    min_g_to_a_3p_strong_threshold: float = 0.10
    min_c_to_t_5p_moderate_threshold: float = 0.03
    min_g_to_a_3p_moderate_threshold: float = 0.03
    min_c_to_t_5p_observations: int = 50
    min_g_to_a_3p_observations: int = 50
    max_median_length_threshold_bp: int = 80
    min_reads_for_damage_model: int = 100
    terminal_positions_to_profile: int = 25

    # Bayesian Model & PMD Corroboration Thresholds
    bayesian_delta_s_strong_threshold: float = 0.20
    bayesian_delta_s_moderate_threshold: float = 0.05
    pmd_corroboration_fraction_ge_3: float = 0.10

    # Protocol-Specific Diagnostic Thresholds
    udg_partial_max_pos2_to_pos1_ratio: float = 0.70
    udg_partial_min_pos2_damage_floor: float = 0.03  # Noise threshold above which residual pos 2 damage triggers asymmetry warning
    terminal_asymmetry_max_difference: float = 0.15

    def to_dict(self) -> Dict[str, Any]:
        """Serializes policy configuration to a dictionary."""
        return {
            "input_bam_type": self.input_bam_type.value,
            "library_treatment": self.library_treatment.value,
            "enable_rescaling": self.enable_rescaling,
            "enable_pmd_scoring": self.enable_pmd_scoring,
            "enable_pmd_filtering": self.enable_pmd_filtering,
            "pmd_threshold": self.pmd_threshold,
            "min_c_to_t_5p_strong_threshold": self.min_c_to_t_5p_strong_threshold,
            "min_g_to_a_3p_strong_threshold": self.min_g_to_a_3p_strong_threshold,
            "min_c_to_t_5p_moderate_threshold": self.min_c_to_t_5p_moderate_threshold,
            "min_g_to_a_3p_moderate_threshold": self.min_g_to_a_3p_moderate_threshold,
            "min_c_to_t_5p_observations": self.min_c_to_t_5p_observations,
            "min_g_to_a_3p_observations": self.min_g_to_a_3p_observations,
            "max_median_length_threshold_bp": self.max_median_length_threshold_bp,
            "min_reads_for_damage_model": self.min_reads_for_damage_model,
            "terminal_positions_to_profile": self.terminal_positions_to_profile,
            "bayesian_delta_s_strong_threshold": self.bayesian_delta_s_strong_threshold,
            "bayesian_delta_s_moderate_threshold": self.bayesian_delta_s_moderate_threshold,
            "pmd_corroboration_fraction_ge_3": self.pmd_corroboration_fraction_ge_3,
            "udg_partial_max_pos2_to_pos1_ratio": self.udg_partial_max_pos2_to_pos1_ratio,
            "udg_partial_min_pos2_damage_floor": self.udg_partial_min_pos2_damage_floor,
            "terminal_asymmetry_max_difference": self.terminal_asymmetry_max_difference,
        }


@dataclass
class DamageModelProfile:
    """Consolidated statistical and biological damage profile for downstream consumption."""

    sample_id: str
    input_bam_type: DamageInputBAMType = DamageInputBAMType.FILTERED_DEDUP
    input_bam_path: Optional[Path] = None
    reference_contig: Optional[str] = None
    library_treatment: LibraryTreatment = LibraryTreatment.NONE
    terminal_damage: TerminalDamageFrequencies = field(default_factory=TerminalDamageFrequencies)
    fragment_length: FragmentLengthMetrics = field(default_factory=FragmentLengthMetrics)
    bayesian_parameters: BayesianDamageParameters = field(default_factory=BayesianDamageParameters)
    pmd_metrics: Optional[PMDMetrics] = None
    damage_signal_level: DamageSignalLevel = DamageSignalLevel.NOISY_OR_INSUFFICIENT
    authentication_verdict: AuthenticationVerdict = AuthenticationVerdict.INCONCLUSIVE
    rescaling_applied: bool = False

    def to_dict(self) -> Dict[str, Any]:
        """Serializes damage model profile to a structured dictionary."""
        return {
            "sample_id": self.sample_id,
            "input_bam_type": self.input_bam_type.value,
            "input_bam_path": str(self.input_bam_path) if self.input_bam_path else None,
            "reference_contig": self.reference_contig,
            "library_treatment": self.library_treatment.value,
            "terminal_damage": self.terminal_damage.to_dict(),
            "fragment_length": self.fragment_length.to_dict(),
            "bayesian_parameters": self.bayesian_parameters.to_dict(),
            "pmd_metrics": self.pmd_metrics.to_dict() if self.pmd_metrics else None,
            "damage_signal_level": self.damage_signal_level.value,
            "authentication_verdict": self.authentication_verdict.value,
            "rescaling_applied": self.rescaling_applied,
        }


@dataclass
class DamageResult:
    """Top-level container encapsulating execution results and QC diagnostics for Block 4."""

    sample_id: str
    status: DamageStatus = DamageStatus.UNTESTED
    policy: DamagePolicy = field(default_factory=DamagePolicy)
    profile: DamageModelProfile = field(default_factory=lambda: DamageModelProfile(sample_id=""))
    outputs: DamageOutputFiles = field(default_factory=DamageOutputFiles)
    warning_codes: List[DamageWarningCode] = field(default_factory=list)
    error_codes: List[DamageErrorCode] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    execution_started: Optional[str] = None
    execution_finished: Optional[str] = None
    runtime_seconds: Optional[float] = None

    def __post_init__(self) -> None:
        """Aligns sample_id across internal profile container."""
        if not self.profile.sample_id and self.sample_id:
            self.profile.sample_id = self.sample_id

    def mark_started(self) -> None:
        """Records ISO start timestamp."""
        self.execution_started = datetime.now(timezone.utc).isoformat()

    def mark_finished(self) -> None:
        """Records ISO finish timestamp and computes duration."""
        self.execution_finished = datetime.now(timezone.utc).isoformat()
        if self.execution_started:
            t0 = datetime.fromisoformat(self.execution_started)
            t1 = datetime.fromisoformat(self.execution_finished)
            self.runtime_seconds = round((t1 - t0).total_seconds(), 2)

    def add_warning(self, code: DamageWarningCode, message: str) -> None:
        """Appends a non-fatal warning code and message."""
        if code not in self.warning_codes:
            self.warning_codes.append(code)
        self.warnings.append(message)
        if self.status != DamageStatus.FAIL:
            self.status = DamageStatus.WARN

    def add_error(self, code: DamageErrorCode, message: str) -> None:
        """Appends a fatal error code and transitions status to FAIL."""
        if code not in self.error_codes:
            self.error_codes.append(code)
        self.errors.append(message)
        self.status = DamageStatus.FAIL

    def to_dict(self) -> Dict[str, Any]:
        """Serializes full Block 4 execution result into a structured dictionary."""
        return {
            "sample_id": self.sample_id,
            "status": self.status.value,
            "policy": self.policy.to_dict(),
            "profile": self.profile.to_dict(),
            "outputs": self.outputs.to_dict(),
            "warning_codes": [w.value for w in self.warning_codes],
            "error_codes": [e.value for e in self.error_codes],
            "warnings": self.warnings,
            "errors": self.errors,
            "timestamps": {
                "started": self.execution_started,
                "finished": self.execution_finished,
                "runtime_seconds": self.runtime_seconds,
            },
        }
