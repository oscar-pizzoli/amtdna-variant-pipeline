"""Domain models, enums, dataclasses, and validation schemas for Block 5 (mtDNA Contamination Estimation)."""

from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from src.modules.damage.models import LibraryTreatment


class ContaminationMethod(str, Enum):
    """Methodological approach selected for contamination estimation or triage."""

    SCHMUTZI = "SCHMUTZI"                        # Full iterative Bayesian mixture model (aDNA gold standard)
    FAST_SCREEN_TRIAGE = "FAST_SCREEN_TRIAGE"    # Conservative pre-screening triage (FAST_CLEAR)
    DISSENT_SCREENING = "DISSENT_SCREENING"      # Direct allelic dissent / empirical screening
    NONE = "NONE"                                # Contamination estimation disabled or bypassed


class ScreeningDecision(str, Enum):
    """Conservative triage decisions emitted by the FastContaminationScreen engine."""

    FAST_CLEAR = "FAST_CLEAR"                    # No significant allelic dissent detected; escalation unneeded
    ESCALATE = "ESCALATE"                        # Discordant allelic signal detected; mandatory Schmutzi execution
    INSUFFICIENT_DATA = "INSUFFICIENT_DATA"      # Insufficient depth or sites; safety escalation or manual review


class ContaminationStatus(str, Enum):
    """Overall execution status of Block 5."""

    UNTESTED = "UNTESTED"
    PASS = "PASS"
    WARN = "WARN"
    FAIL = "FAIL"


class ContaminationCategory(str, Enum):
    """Configured operational categorization of the estimated contamination fraction."""

    LOW = "LOW"                                  # Baseline noise / minimal contamination (typically <= 2%)
    MODERATE = "MODERATE"                        # Acceptable for standard haplogrouping (typically 2% - 5%)
    ELEVATED = "ELEVATED"                        # Elevated contamination requiring strict consensus gating (5% - 10%)
    CRITICAL = "CRITICAL"                        # Severely contaminated (typically > 10%)
    INSUFFICIENT_DATA = "INSUFFICIENT_DATA"      # Depth or informative sites too low for reliable categorization


class EstimationConfidence(str, Enum):
    """Statistical confidence level of the estimation based on depth, CI width, and convergence."""

    HIGH = "HIGH"                                # Adequate depth, narrow 95% CI, verified convergence
    MODERATE = "MODERATE"                        # Acceptable depth and convergence, moderate CI width
    LOW = "LOW"                                  # Low depth (< min coverage) or broad 95% CI
    UNRELIABLE = "UNRELIABLE"                    # MCMC failed to converge or contradictory signal metrics


class ConsensusGatingAction(str, Enum):
    """Recommended action and calling constraints passed downstream to Block 6 (Consensus & Variant Calling)."""

    STANDARD_CALLING = "STANDARD_CALLING"                          # Unrestricted standard consensus & variant calling
    RESTRICT_HETEROPLASMY = "RESTRICT_HETEROPLASMY"                # Suppress low-frequency heteroplasmy / raise VAF threshold
    USE_ENDOGENOUS_SCHMUTZI_CONSENSUS = "USE_ENDOGENOUS_SCHMUTZI_CONSENSUS"  # Use schmutzi-predicted endogenous FASTA
    PMD_FILTERING_REQUIRED = "PMD_FILTERING_REQUIRED"              # Require read partitioning by minimum PMD score
    EXCLUDE_OR_MANUAL_REVIEW = "EXCLUDE_OR_MANUAL_REVIEW"          # Flag sample as unreliable; manual review recommended


class ContaminationWarningCode(str, Enum):
    """Diagnostic warning codes emitted during contamination analysis."""

    LOW_COVERAGE_FOR_ESTIMATION = "LOW_COVERAGE_FOR_ESTIMATION"
    BROAD_CREDIBLE_INTERVAL = "BROAD_CREDIBLE_INTERVAL"
    MCMC_CONVERGENCE_FAILURE = "MCMC_CONVERGENCE_FAILURE"
    CONTAMINANT_DEAMINATION_DETECTED = "CONTAMINANT_DEAMINATION_DETECTED"
    MULTIPLE_CONTAMINANTS_SUSPECTED = "MULTIPLE_CONTAMINANTS_SUSPECTED"
    DATABASE_NOT_CONFIGURED = "DATABASE_NOT_CONFIGURED"
    DATABASE_NOT_FOUND = "DATABASE_NOT_FOUND"
    FALLBACK_TO_DISSENT = "FALLBACK_TO_DISSENT"
    LOW_INFORMATIVE_SITES = "LOW_INFORMATIVE_SITES"
    FAST_SCREEN_ESCALATION = "FAST_SCREEN_ESCALATION"
    FAST_SCREEN_INSUFFICIENT_DATA = "FAST_SCREEN_INSUFFICIENT_DATA"


class ContaminationErrorCode(str, Enum):
    """Fatal error codes halting downstream processing."""

    SCHMUTZI_EXECUTION_FAILED = "SCHMUTZI_EXECUTION_FAILED"
    DISSENT_EXECUTION_FAILED = "DISSENT_EXECUTION_FAILED"          # Standardizzato per BatchRunner / DissentEngine
    DISSENT_CALCULATION_FAILED = "DISSENT_EXECUTION_FAILED"        # Alias di retrocompatibilità
    FAST_SCREEN_FAILED = "FAST_SCREEN_FAILED"
    PREPARER_FAILED = "PREPARER_FAILED"
    CORRUPTED_OUTPUT = "CORRUPTED_OUTPUT"
    INPUT_BAM_NOT_FOUND = "INPUT_BAM_NOT_FOUND"
    REFERENCE_INCOMPATIBLE = "REFERENCE_INCOMPATIBLE"
    DATABASE_MISSING = "DATABASE_MISSING"


@dataclass
class SubsamplingPolicy:
    """Configuration policy for conditional depth subsampling prior to contamination estimation."""

    enabled: bool = True
    target_depth: float = 40.0
    seed: int = 42
    min_depth_threshold: float = 42.0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "enabled": self.enabled,
            "target_depth": self.target_depth,
            "seed": self.seed,
            "min_depth_threshold": self.min_depth_threshold,
        }


@dataclass
class FastScreeningPolicy:
    """Execution parameters and escalation criteria for the fast pre-screening triage engine."""

    enabled: bool = False
    algorithm_version: str = "fast-screen-v2.0"
    diagnostic_panel_path: Optional[Path] = None
    diagnostic_panel_version: str = "phylotree17_v1"
    escalate_dissenting_sites_min: int = 3
    escalate_vaf_threshold: float = 0.02
    min_alt_reads_per_site: int = 2
    min_depth_for_screen: float = 10.0
    min_genome_depth: float = 10.0
    min_evaluated_sites: int = 50
    min_panel_evaluation_fraction: float = 0.25
    min_base_quality: int = 30
    min_mapping_quality: int = 30
    min_depth_per_site: int = 10
    noise_vaf_floor: float = 0.015
    cache_screen_results: bool = True
    cache_mode: str = "strict"

    def compute_policy_sha256(self) -> str:
        data = self.to_dict()
        data.pop("cache_screen_results", None)
        data.pop("cache_mode", None)
        encoded = json.dumps(data, sort_keys=True).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()

    def to_dict(self) -> Dict[str, Any]:
        return {
            "enabled": self.enabled,
            "algorithm_version": self.algorithm_version,
            "diagnostic_panel_path": str(self.diagnostic_panel_path) if self.diagnostic_panel_path else None,
            "diagnostic_panel_version": self.diagnostic_panel_version,
            "escalate_dissenting_sites_min": self.escalate_dissenting_sites_min,
            "escalate_vaf_threshold": self.escalate_vaf_threshold,
            "min_alt_reads_per_site": self.min_alt_reads_per_site,
            "min_depth_for_screen": self.min_depth_for_screen,
            "min_genome_depth": self.min_genome_depth,
            "min_evaluated_sites": self.min_evaluated_sites,
            "min_panel_evaluation_fraction": self.min_panel_evaluation_fraction,
            "min_base_quality": self.min_base_quality,
            "min_mapping_quality": self.min_mapping_quality,
            "min_depth_per_site": self.min_depth_per_site,
            "noise_vaf_floor": self.noise_vaf_floor,
            "cache_screen_results": self.cache_screen_results,
            "cache_mode": self.cache_mode,
        }


@dataclass
class FastScreenResult:
    """Quantitative summary and conservative triage decision produced by FastContaminationScreen."""

    sample_id: str
    decision: ScreeningDecision = ScreeningDecision.INSUFFICIENT_DATA
    bypass_eligible: bool = False
    mean_genome_depth: float = 0.0
    mean_diagnostic_depth: float = 0.0
    panel_total_sites: int = 0
    evaluated_diagnostic_positions: int = 0
    panel_evaluation_fraction: float = 0.0
    dissenting_positions_count: int = 0
    total_informative_reads: int = 0
    total_alt_reads: int = 0
    global_dissent_fraction: float = 0.0
    max_dissent_vaf: float = 0.0
    mean_dissent_vaf: float = 0.0
    median_dissent_vaf: float = 0.0
    median_depth_at_dissent_sites: Optional[float] = None
    applied_thresholds: Dict[str, Any] = field(default_factory=dict)
    interpretation: str = "no escalation evidence detected"
    diagnostic_panel_path: Optional[str] = None
    diagnostic_panel_sha256: Optional[str] = None
    diagnostic_panel_version: Optional[str] = None
    policy_sha256: Optional[str] = None
    algorithm_version: str = "fast-screen-v2.0"
    input_bam_sha256: Optional[str] = None
    reference_sha256: Optional[str] = None
    runtime_seconds: float = 0.0

    def validate_invariants(self) -> bool:
        if not (0.0 <= self.max_dissent_vaf <= 1.0):
            return False
        if not (0.0 <= self.mean_dissent_vaf <= 1.0):
            return False
        if not (0.0 <= self.median_dissent_vaf <= 1.0):
            return False
        if self.mean_genome_depth < 0.0 or self.mean_diagnostic_depth < 0.0:
            return False
        if self.dissenting_positions_count < 0 or self.evaluated_diagnostic_positions < 0:
            return False
        if self.dissenting_positions_count > self.evaluated_diagnostic_positions:
            return False
        if self.total_alt_reads < 0 or self.total_informative_reads < 0:
            return False
        if self.total_alt_reads > self.total_informative_reads:
            return False
        if not (0.0 <= self.global_dissent_fraction <= 1.0):
            return False
        return True

    def to_dict(self) -> Dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "decision": self.decision.value,
            "bypass_eligible": self.bypass_eligible,
            "mean_genome_depth": round(self.mean_genome_depth, 2),
            "mean_diagnostic_depth": round(self.mean_diagnostic_depth, 2),
            "panel_total_sites": self.panel_total_sites,
            "evaluated_diagnostic_positions": self.evaluated_diagnostic_positions,
            "panel_evaluation_fraction": round(self.panel_evaluation_fraction, 4),
            "dissenting_positions_count": self.dissenting_positions_count,
            "total_informative_reads": self.total_informative_reads,
            "total_alt_reads": self.total_alt_reads,
            "global_dissent_fraction": round(self.global_dissent_fraction, 4),
            "max_dissent_vaf": round(self.max_dissent_vaf, 4),
            "mean_dissent_vaf": round(self.mean_dissent_vaf, 4),
            "median_dissent_vaf": round(self.median_dissent_vaf, 4),
            "median_depth_at_dissent_sites": (
                round(self.median_depth_at_dissent_sites, 2)
                if self.median_depth_at_dissent_sites is not None
                else None
            ),
            "applied_thresholds": self.applied_thresholds,
            "interpretation": self.interpretation,
            "diagnostic_panel_path": self.diagnostic_panel_path,
            "diagnostic_panel_sha256": self.diagnostic_panel_sha256,
            "diagnostic_panel_version": self.diagnostic_panel_version,
            "policy_sha256": self.policy_sha256,
            "algorithm_version": self.algorithm_version,
            "input_bam_sha256": self.input_bam_sha256,
            "reference_sha256": self.reference_sha256,
            "runtime_seconds": round(self.runtime_seconds, 2),
        }


@dataclass
class ContaminationInputContext:
    """Unified input context container passed to Block 5 orchestration."""

    sample_id: str
    input_bam_path: Path
    reference_fasta_path: Path
    reference_contig: str
    reference_length_bp: int
    contaminant_database_path: Optional[Path] = None
    is_unrescaled_bam_verified: bool = True
    is_ancient: bool = True

    def to_dict(self) -> Dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "input_bam_path": str(self.input_bam_path),
            "reference_fasta_path": str(self.reference_fasta_path),
            "reference_contig": self.reference_contig,
            "reference_length_bp": self.reference_length_bp,
            "contaminant_database_path": (
                str(self.contaminant_database_path) if self.contaminant_database_path else None
            ),
            "is_unrescaled_bam_verified": self.is_unrescaled_bam_verified,
            "is_ancient": self.is_ancient,
        }


@dataclass
class WrappingMetadata:
    """Metadata tracking circular mitochondrial origin padding and coordinate transformations."""

    applied: bool = False
    padding_bp: int = 0
    wrapped_contig_name: Optional[str] = None
    original_length_bp: Optional[int] = None
    extended_length_bp: Optional[int] = None
    remapping_offset: int = 0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "applied": self.applied,
            "padding_bp": self.padding_bp,
            "wrapped_contig_name": self.wrapped_contig_name,
            "original_length_bp": self.original_length_bp,
            "extended_length_bp": self.extended_length_bp,
            "remapping_offset": self.remapping_offset,
        }


@dataclass
class ContaminationEstimate:
    """Quantitative contamination metrics and statistical intervals (Level 1)."""

    point_estimate: float = 0.0
    ci_lower_95: Optional[float] = None
    ci_upper_95: Optional[float] = None
    method_used: ContaminationMethod = ContaminationMethod.NONE
    cont_deam_prior_estimate: Optional[float] = None
    log_likelihood: Optional[float] = None
    heterogeneity_fraction: Optional[float] = None

    def __post_init__(self) -> None:
        if not (0.0 <= self.point_estimate <= 1.0):
            raise ValueError(f"point_estimate must be in [0.0, 1.0], got {self.point_estimate}")

        if self.ci_lower_95 is not None and not (0.0 <= self.ci_lower_95 <= 1.0):
            raise ValueError(f"ci_lower_95 must be in [0.0, 1.0], got {self.ci_lower_95}")

        if self.ci_upper_95 is not None and not (0.0 <= self.ci_upper_95 <= 1.0):
            raise ValueError(f"ci_upper_95 must be in [0.0, 1.0], got {self.ci_upper_95}")

        if self.ci_lower_95 is not None and self.ci_upper_95 is not None:
            if self.ci_lower_95 > self.ci_upper_95 + 1e-5:
                raise ValueError(
                    f"ci_lower_95 ({self.ci_lower_95}) cannot be greater than ci_upper_95 ({self.ci_upper_95})"
                )
            if not (self.ci_lower_95 - 1e-5 <= self.point_estimate <= self.ci_upper_95 + 1e-5):
                raise ValueError(
                    f"point_estimate ({self.point_estimate}) must fall within [ci_lower_95 ({self.ci_lower_95}), "
                    f"ci_upper_95 ({self.ci_upper_95})]"
                )

        if self.heterogeneity_fraction is not None and not (0.0 <= self.heterogeneity_fraction <= 1.0):
            raise ValueError(
                f"heterogeneity_fraction must be in [0.0, 1.0], got {self.heterogeneity_fraction}"
            )

    @property
    def ci_width(self) -> Optional[float]:
        if self.ci_lower_95 is not None and self.ci_upper_95 is not None:
            return max(0.0, self.ci_upper_95 - self.ci_lower_95)
        return None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "point_estimate_fraction": round(self.point_estimate, 4),
            "ci_lower_95_fraction": round(self.ci_lower_95, 4) if self.ci_lower_95 is not None else None,
            "ci_upper_95_fraction": round(self.ci_upper_95, 4) if self.ci_upper_95 is not None else None,
            "ci_width": round(self.ci_width, 4) if self.ci_width is not None else None,
            "method_used": self.method_used.value,
            "cont_deam_prior_estimate_fraction": (
                round(self.cont_deam_prior_estimate, 4) if self.cont_deam_prior_estimate is not None else None
            ),
            "log_likelihood": round(self.log_likelihood, 4) if self.log_likelihood is not None else None,
            "heterogeneity_fraction": (
                round(self.heterogeneity_fraction, 4) if self.heterogeneity_fraction is not None else None
            ),
        }


@dataclass
class ModelQualityDiagnostics:
    """Model quality, convergence, depth metrics, CI width, and risk indicators (Level 2)."""

    mean_depth: float = 0.0
    informative_positions_evaluated: int = 0
    ci_width: Optional[float] = None
    model_converged: bool = False
    iterations_run: int = 0
    warn_contaminant_damage_pattern: bool = False
    possible_multiple_contaminants: bool = False
    estimation_confidence: EstimationConfidence = EstimationConfidence.UNRELIABLE
    contaminant_haplogroups: List[str] = field(default_factory=list)
    wrapping: WrappingMetadata = field(default_factory=WrappingMetadata)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "mean_depth": round(self.mean_depth, 2),
            "informative_positions_evaluated": self.informative_positions_evaluated,
            "ci_width": round(self.ci_width, 4) if self.ci_width is not None else None,
            "model_converged": self.model_converged,
            "iterations_run": self.iterations_run,
            "warn_contaminant_damage_pattern": self.warn_contaminant_damage_pattern,
            "possible_multiple_contaminants": self.possible_multiple_contaminants,
            "estimation_confidence": self.estimation_confidence.value,
            "contaminant_haplogroups": self.contaminant_haplogroups,
            "wrapping": self.wrapping.to_dict(),
        }


@dataclass
class GatingRecommendation:
    """Categorical interpretation and formal constraints for downstream Block 6 execution (Level 3)."""

    contamination_category: ContaminationCategory = ContaminationCategory.INSUFFICIENT_DATA
    recommended_action: ConsensusGatingAction = ConsensusGatingAction.EXCLUDE_OR_MANUAL_REVIEW
    preferred_consensus_source: Optional[str] = None
    required_pmd_threshold: Optional[float] = None
    summary_message: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "contamination_category": self.contamination_category.value,
            "recommended_action": self.recommended_action.value,
            "preferred_consensus_source": self.preferred_consensus_source,
            "required_pmd_threshold": self.required_pmd_threshold,
            "summary_message": self.summary_message,
        }


@dataclass
class ContaminationOutputFiles:
    """Physical filesystem paths of Block 5 deliverables."""

    schmutzi_ready_bam: Optional[Path] = None
    schmutzi_ready_bai: Optional[Path] = None
    cont_est_file: Optional[Path] = None
    cont_deam_file: Optional[Path] = None
    endogenous_fasta: Optional[Path] = None
    contaminant_fasta: Optional[Path] = None
    dissent_tsv: Optional[Path] = None
    contamination_json: Optional[Path] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "schmutzi_ready_bam": str(self.schmutzi_ready_bam) if self.schmutzi_ready_bam else None,
            "schmutzi_ready_bai": str(self.schmutzi_ready_bai) if self.schmutzi_ready_bai else None,
            "cont_est_file": str(self.cont_est_file) if self.cont_est_file else None,
            "cont_deam_file": str(self.cont_deam_file) if self.cont_deam_file else None,
            "endogenous_fasta": str(self.endogenous_fasta) if self.endogenous_fasta else None,
            "contaminant_fasta": str(self.contaminant_fasta) if self.contaminant_fasta else None,
            "dissent_tsv": str(self.dissent_tsv) if self.dissent_tsv else None,
            "contamination_json": str(self.contamination_json) if self.contamination_json else None,
        }


@dataclass
class ContaminationAuditMetadata:
    """Audit trail tracking input BAM provenance, reference hashes, tool versions, and commands."""

    input_bam_path: Optional[Path] = None
    input_bam_sha256: Optional[str] = None
    is_unrescaled_bam_verified: bool = False
    reference_fasta_path: Optional[Path] = None
    reference_fasta_sha256: Optional[str] = None
    wrapped_reference_path: Optional[Path] = None
    wrapped_reference_sha256: Optional[str] = None
    schmutzi_version: Optional[str] = None
    contaminant_database_path: Optional[Path] = None
    contaminant_database_sha256: Optional[str] = None
    executed_commands: List[str] = field(default_factory=list)
    step_runtimes_seconds: Dict[str, float] = field(default_factory=dict)

    # Subsampling Audit
    subsampling_applied: bool = False
    original_mean_depth: Optional[float] = None
    downsampled_mean_depth: Optional[float] = None
    subsampling_fraction: Optional[float] = None

    # Fast Pre-screening Audit
    fast_screen_applied: bool = False
    fast_screen_result: Optional[FastScreenResult] = None
    full_bayesian_estimation_performed: bool = False
    screened_without_full_quantification: bool = False

    def to_dict(self) -> Dict[str, Any]:
        return {
            "input_bam_path": str(self.input_bam_path) if self.input_bam_path else None,
            "input_bam_sha256": self.input_bam_sha256,
            "is_unrescaled_bam_verified": self.is_unrescaled_bam_verified,
            "reference_fasta_path": str(self.reference_fasta_path) if self.reference_fasta_path else None,
            "reference_fasta_sha256": self.reference_fasta_sha256,
            "wrapped_reference_path": str(self.wrapped_reference_path) if self.wrapped_reference_path else None,
            "wrapped_reference_sha256": self.wrapped_reference_sha256,
            "schmutzi_version": self.schmutzi_version,
            "contaminant_database_path": (
                str(self.contaminant_database_path) if self.contaminant_database_path else None
            ),
            "contaminant_database_sha256": self.contaminant_database_sha256,
            "executed_commands": self.executed_commands,
            "step_runtimes_seconds": {k: round(v, 2) for k, v in self.step_runtimes_seconds.items()},
            "subsampling_applied": self.subsampling_applied,
            "original_mean_depth": self.original_mean_depth,
            "downsampled_mean_depth": self.downsampled_mean_depth,
            "subsampling_fraction": self.subsampling_fraction,
            "fast_screen_applied": self.fast_screen_applied,
            "fast_screen_result": self.fast_screen_result.to_dict() if self.fast_screen_result else None,
            "full_bayesian_estimation_performed": self.full_bayesian_estimation_performed,
            "screened_without_full_quantification": self.screened_without_full_quantification,
        }


@dataclass
class SchmutziPolicy:
    """Execution parameters specific to schmutzi Bayesian mixture modeling."""

    contaminant_database: Optional[Path] = None
    use_wrapped_reference: bool = False
    padding_bp: int = 200
    max_iterations: int = 10
    threads: int = 4
    library_treatment: LibraryTreatment = LibraryTreatment.NONE

    def to_dict(self) -> Dict[str, Any]:
        return {
            "contaminant_database": str(self.contaminant_database) if self.contaminant_database else None,
            "use_wrapped_reference": self.use_wrapped_reference,
            "padding_bp": self.padding_bp,
            "max_iterations": self.max_iterations,
            "threads": self.threads,
            "library_treatment": self.library_treatment.value,
        }


@dataclass
class DissentPolicy:
    """Execution parameters specific to Track B allelic dissent and heterogeneity screening."""

    min_base_quality: int = 30
    min_mapping_quality: int = 30
    min_depth_per_site: int = 10
    noise_vaf_floor: float = 0.015
    base_min_dissenting_sites: int = 3
    dissent_high_confidence_min_sites: int = 50

    def to_dict(self) -> Dict[str, Any]:
        return {
            "min_base_quality": self.min_base_quality,
            "min_mapping_quality": self.min_mapping_quality,
            "min_depth_per_site": self.min_depth_per_site,
            "noise_vaf_floor": self.noise_vaf_floor,
            "base_min_dissenting_sites": self.base_min_dissenting_sites,
            "dissent_high_confidence_min_sites": self.dissent_high_confidence_min_sites,
        }


@dataclass
class ContaminationPolicy:
    """Configurable execution parameters, heuristic boundaries, and gating criteria for Block 5."""

    policy_id: str = "contamination_standard_v1"
    method: ContaminationMethod = ContaminationMethod.SCHMUTZI
    enabled: bool = True
    allow_fallback_to_dissent: bool = False
    schmutzi: SchmutziPolicy = field(default_factory=SchmutziPolicy)
    dissent: DissentPolicy = field(default_factory=DissentPolicy)
    subsampling: SubsamplingPolicy = field(default_factory=SubsamplingPolicy)
    fast_screening: FastScreeningPolicy = field(default_factory=FastScreeningPolicy)

    # Thresholds: Category Boundaries (Fractions in [0.0, 1.0])
    low_threshold: float = 0.02
    moderate_threshold: float = 0.05
    elevated_threshold: float = 0.10

    # Thresholds: Model Quality Criteria
    min_coverage_for_estimation: float = 5.0
    recommended_coverage: float = 10.0
    high_confidence_max_ci_width: float = 0.04
    moderate_confidence_max_ci_width: float = 0.08
    min_informative_sites: int = 10
    dissent_high_confidence_min_sites: int = 50
    pmd_filtering_threshold: float = 3.0

    def __post_init__(self) -> None:
        if not (0.0 <= self.low_threshold < self.moderate_threshold < self.elevated_threshold <= 1.0):
            raise ValueError(
                f"Thresholds must satisfy 0.0 <= low < moderate < elevated <= 1.0. "
                f"Got: low={self.low_threshold}, moderate={self.moderate_threshold}, elevated={self.elevated_threshold}"
            )
        if self.min_coverage_for_estimation <= 0.0:
            raise ValueError(
                f"min_coverage_for_estimation must be positive, got {self.min_coverage_for_estimation}"
            )
        if self.recommended_coverage < self.min_coverage_for_estimation:
            raise ValueError(
                f"recommended_coverage ({self.recommended_coverage}) cannot be lower than "
                f"min_coverage_for_estimation ({self.min_coverage_for_estimation})"
            )
        if not (0.0 < self.high_confidence_max_ci_width <= self.moderate_confidence_max_ci_width <= 1.0):
            raise ValueError(
                f"CI width thresholds must satisfy 0.0 < high <= moderate <= 1.0. "
                f"Got: high={self.high_confidence_max_ci_width}, moderate={self.moderate_confidence_max_ci_width}"
            )
        if self.min_informative_sites <= 0:
            raise ValueError(
                f"min_informative_sites must be a positive integer, got {self.min_informative_sites}"
            )
        if self.pmd_filtering_threshold < 0.0:
            raise ValueError(
                f"pmd_filtering_threshold must be non-negative, got {self.pmd_filtering_threshold}"
            )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "policy_id": self.policy_id,
            "method": self.method.value,
            "enabled": self.enabled,
            "allow_fallback_to_dissent": self.allow_fallback_to_dissent,
            "schmutzi": self.schmutzi.to_dict(),
            "dissent": self.dissent.to_dict(),
            "subsampling": self.subsampling.to_dict(),
            "fast_screening": self.fast_screening.to_dict(),
            "thresholds": {
                "category_limits": {
                    "low_threshold": self.low_threshold,
                    "moderate_threshold": self.moderate_threshold,
                    "elevated_threshold": self.elevated_threshold,
                },
                "model_quality": {
                    "min_coverage_for_estimation": self.min_coverage_for_estimation,
                    "recommended_coverage": self.recommended_coverage,
                    "high_confidence_max_ci_width": self.high_confidence_max_ci_width,
                    "moderate_confidence_max_ci_width": self.moderate_confidence_max_ci_width,
                    "min_informative_sites": self.min_informative_sites,
                    "dissent_high_confidence_min_sites": self.dissent_high_confidence_min_sites,
                    "pmd_filtering_threshold": self.pmd_filtering_threshold,
                },
            },
        }


@dataclass
class ContaminationResult:
    """Top-level deliverable container encapsulating Block 5 execution, quality, audit, and gating."""

    sample_id: str
    status: ContaminationStatus = ContaminationStatus.UNTESTED
    policy: ContaminationPolicy = field(default_factory=ContaminationPolicy)
    estimate: ContaminationEstimate = field(default_factory=ContaminationEstimate)
    diagnostics: ModelQualityDiagnostics = field(default_factory=ModelQualityDiagnostics)
    gating: GatingRecommendation = field(default_factory=GatingRecommendation)
    outputs: ContaminationOutputFiles = field(default_factory=ContaminationOutputFiles)
    audit: ContaminationAuditMetadata = field(default_factory=ContaminationAuditMetadata)
    fallback_reason: Optional[str] = None
    warning_codes: List[ContaminationWarningCode] = field(default_factory=list)
    error_codes: List[ContaminationErrorCode] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    execution_started: Optional[str] = None
    execution_finished: Optional[str] = None
    runtime_seconds: Optional[float] = None

    @property
    def started_at(self) -> Optional[datetime]:
        """Backward-compatibility property returning datetime for started_at."""
        if self.execution_started:
            try:
                return datetime.fromisoformat(self.execution_started)
            except Exception:
                return None
        return None

    @started_at.setter
    def started_at(self, dt: Optional[Union[datetime, str]]) -> None:
        """Backward-compatibility setter accepting datetime or ISO string."""
        if dt is None:
            self.execution_started = None
        elif isinstance(dt, datetime):
            self.execution_started = dt.isoformat()
        else:
            self.execution_started = str(dt)

    def mark_started(self) -> None:
        self.execution_started = datetime.now(timezone.utc).isoformat()

    def mark_finished(self) -> None:
        self.execution_finished = datetime.now(timezone.utc).isoformat()
        if self.execution_started:
            try:
                t0 = datetime.fromisoformat(self.execution_started)
                t1 = datetime.fromisoformat(self.execution_finished)
                self.runtime_seconds = round((t1 - t0).total_seconds(), 2)
            except Exception:
                pass

    def add_warning(self, code: ContaminationWarningCode, message: str) -> None:
        if code not in self.warning_codes:
            self.warning_codes.append(code)
        self.warnings.append(message)
        if self.status != ContaminationStatus.FAIL:
            self.status = ContaminationStatus.WARN

    def add_error(self, code: ContaminationErrorCode, message: str) -> None:
        if code not in self.error_codes:
            self.error_codes.append(code)
        self.errors.append(message)
        self.status = ContaminationStatus.FAIL

    def to_dict(self) -> Dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "status": self.status.value,
            "policy": self.policy.to_dict(),
            "estimate": self.estimate.to_dict(),
            "diagnostics": self.diagnostics.to_dict(),
            "gating": self.gating.to_dict(),
            "outputs": self.outputs.to_dict(),
            "audit": self.audit.to_dict(),
            "fallback_reason": self.fallback_reason,
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