"""Quality Control & Gating Engine for mtDNA Contamination and Heterogeneity.

Synthesizes quantitative estimates, credible intervals, diagnostic depth,
and FastScreen triage results to establish model confidence, assign categorical
risk tiers, and output formal gating directives for Block 6.
"""

from enum import Enum
import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from src.modules.contamination.models import (
    ConsensusGatingAction,
    ContaminationAuditMetadata,
    ContaminationCategory,
    ContaminationErrorCode,
    ContaminationEstimate,
    ContaminationMethod,
    ContaminationOutputFiles,
    ContaminationPolicy,
    ContaminationResult,
    ContaminationStatus,
    ContaminationWarningCode,
    EstimationConfidence,
    GatingRecommendation,
    ModelQualityDiagnostics,
    ScreeningDecision,
)

logger = logging.getLogger("mtDNA_Pipeline.ContaminationQCEvaluator")


def _contamination_json_serializer(obj: Any) -> Any:
    """JSON encoder helper ensuring clean Enum.value and string Path serialization."""
    if isinstance(obj, Enum):
        return obj.value
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


class ContaminationQCValidator:
    """Evaluates contamination model validity, assigns confidence tiers, and formulates Block 6 gating directives."""

    @classmethod
    def _deduplicate_items(cls, items: List[Any]) -> List[Any]:
        """Preserves emission order while removing duplicated diagnostic codes or identical string messages."""
        return list(dict.fromkeys(items))

    @classmethod
    def evaluate_model_confidence(
        cls,
        estimate: ContaminationEstimate,
        diagnostics: ModelQualityDiagnostics,
        policy: ContaminationPolicy,
    ) -> Tuple[EstimationConfidence, List[ContaminationWarningCode], List[str]]:
        """Evaluates statistical reliability based on method-specific quality metrics."""
        warnings: List[str] = []
        codes: List[ContaminationWarningCode] = []

        has_adequate_depth = diagnostics.mean_depth >= policy.recommended_coverage
        has_minimum_depth = diagnostics.mean_depth >= policy.min_coverage_for_estimation

        if not has_minimum_depth:
            codes.append(ContaminationWarningCode.LOW_COVERAGE_FOR_ESTIMATION)
            warnings.append(
                f"Mean depth ({diagnostics.mean_depth:.1f}x) is below the minimum threshold "
                f"({policy.min_coverage_for_estimation:.1f}x) required for reliable evaluation."
            )

        has_sufficient_sites = (
            diagnostics.informative_positions_evaluated >= policy.min_informative_sites
        )
        if not has_sufficient_sites and estimate.point_estimate > policy.low_threshold:
            codes.append(ContaminationWarningCode.LOW_INFORMATIVE_SITES)
            warnings.append(
                f"Informative diagnostic sites evaluated ({diagnostics.informative_positions_evaluated}) "
                f"is below the recommended minimum ({policy.min_informative_sites}) for non-zero contamination."
            )

        # Risk markers per Track A (Bayesian MCMC)
        if diagnostics.warn_contaminant_damage_pattern:
            codes.append(ContaminationWarningCode.CONTAMINANT_DEAMINATION_DETECTED)
            warnings.append(
                "Potential deamination detected on the contaminant fraction; "
                "schmutzi point estimate may underestimate actual modern contamination."
            )

        if diagnostics.possible_multiple_contaminants:
            codes.append(ContaminationWarningCode.MULTIPLE_CONTAMINANTS_SUSPECTED)
            warnings.append(
                "Divergent minor allele patterns suggest a mixture of multiple contaminant sources."
            )

        # --- Track A (SCHMUTZI): Bayesian Mixture Model ---
        if estimate.method_used == ContaminationMethod.SCHMUTZI:
            if not diagnostics.model_converged:
                codes.append(ContaminationWarningCode.MCMC_CONVERGENCE_FAILURE)
                warnings.append("schmutzi MCMC sampling did not achieve nominal convergence.")
                return EstimationConfidence.UNRELIABLE, codes, warnings

            if not has_minimum_depth:
                return EstimationConfidence.LOW, codes, warnings

            ci_width = diagnostics.ci_width
            if ci_width is None and estimate.ci_lower_95 is not None and estimate.ci_upper_95 is not None:
                ci_width = max(0.0, estimate.ci_upper_95 - estimate.ci_lower_95)

            if ci_width is None:
                codes.append(ContaminationWarningCode.BROAD_CREDIBLE_INTERVAL)
                warnings.append("Bayesian Credible Interval is missing or undefined.")
                return EstimationConfidence.LOW, codes, warnings

            if ci_width > policy.moderate_confidence_max_ci_width:
                codes.append(ContaminationWarningCode.BROAD_CREDIBLE_INTERVAL)
                warnings.append(
                    f"95% Credible Interval is excessively broad (width={ci_width:.2%}, "
                    f"threshold={policy.moderate_confidence_max_ci_width:.2%})."
                )
                return EstimationConfidence.LOW, codes, warnings

            if ci_width <= policy.high_confidence_max_ci_width and has_adequate_depth and has_sufficient_sites:
                return EstimationConfidence.HIGH, codes, warnings

            return EstimationConfidence.MODERATE, codes, warnings

        # --- Pre-Screening Triage (FAST_SCREEN_TRIAGE): Heuristic Triage Path ---
        if estimate.method_used == ContaminationMethod.FAST_SCREEN_TRIAGE:
            if not has_minimum_depth or not has_sufficient_sites:
                return EstimationConfidence.LOW, codes, warnings

            if has_adequate_depth and has_sufficient_sites:
                return EstimationConfidence.HIGH, codes, warnings

            return EstimationConfidence.MODERATE, codes, warnings

        # --- Track B (DISSENT_SCREENING): Non-MCMC Empirical Screen ---
        if estimate.method_used == ContaminationMethod.DISSENT_SCREENING:
            if not has_minimum_depth or not has_sufficient_sites:
                return EstimationConfidence.LOW, codes, warnings

            high_conf_sites = getattr(policy, "dissent_high_confidence_min_sites", 50)
            if has_adequate_depth and diagnostics.informative_positions_evaluated >= high_conf_sites:
                return EstimationConfidence.HIGH, codes, warnings

            return EstimationConfidence.MODERATE, codes, warnings

        return EstimationConfidence.UNRELIABLE, codes, warnings

    @classmethod
    def categorize_contamination(
        cls,
        estimate: ContaminationEstimate,
        confidence: EstimationConfidence,
        diagnostics: ModelQualityDiagnostics,
        policy: ContaminationPolicy,
    ) -> Tuple[ContaminationCategory, List[ContaminationWarningCode], List[str]]:
        """Classifies contamination rate incorporating point estimate and generalized Credible Interval bounds."""
        warnings: List[str] = []
        codes: List[ContaminationWarningCode] = []

        if confidence == EstimationConfidence.UNRELIABLE:
            return ContaminationCategory.INSUFFICIENT_DATA, codes, warnings

        if diagnostics.mean_depth < policy.min_coverage_for_estimation:
            return ContaminationCategory.INSUFFICIENT_DATA, codes, warnings

        # Per il Pre-Screening Triage, la categoria rispecchia l'assenza di segnale
        if estimate.method_used == ContaminationMethod.FAST_SCREEN_TRIAGE:
            if confidence == EstimationConfidence.LOW:
                return ContaminationCategory.LOW, codes, warnings
            return ContaminationCategory.LOW, codes, warnings

        point = estimate.point_estimate
        ci_high = estimate.ci_upper_95

        # Valutazione congiunta conservativa: Point Estimate + Upper Credible Bound
        if point <= policy.low_threshold:
            if ci_high is not None and ci_high > policy.elevated_threshold:
                codes.append(ContaminationWarningCode.BROAD_CREDIBLE_INTERVAL)
                warnings.append(
                    f"Risk category set to INSUFFICIENT_DATA because the upper 95% credible bound ({ci_high:.2%}) "
                    f"exceeds the elevated-risk threshold ({policy.elevated_threshold:.2%})."
                )
                return ContaminationCategory.INSUFFICIENT_DATA, codes, warnings

            if ci_high is not None and ci_high > policy.moderate_threshold:
                codes.append(ContaminationWarningCode.BROAD_CREDIBLE_INTERVAL)
                warnings.append(
                    f"Point estimate is low ({point:.2%}) but upper 95% credible bound ({ci_high:.2%}) "
                    f"spans into moderate risk. Conservative category assigned."
                )
                return ContaminationCategory.MODERATE, codes, warnings

            return ContaminationCategory.LOW, codes, warnings

        elif point <= policy.moderate_threshold:
            if ci_high is not None and ci_high > policy.elevated_threshold:
                codes.append(ContaminationWarningCode.BROAD_CREDIBLE_INTERVAL)
                warnings.append(
                    f"Point estimate is moderate ({point:.2%}) but upper 95% credible bound ({ci_high:.2%}) "
                    f"spans into elevated risk. Conservative category assigned."
                )
                return ContaminationCategory.ELEVATED, codes, warnings

            return ContaminationCategory.MODERATE, codes, warnings

        elif point <= policy.elevated_threshold:
            return ContaminationCategory.ELEVATED, codes, warnings

        else:
            return ContaminationCategory.CRITICAL, codes, warnings

    @classmethod
    def formulate_gating_recommendation(
        cls,
        category: ContaminationCategory,
        confidence: EstimationConfidence,
        estimate: ContaminationEstimate,
        outputs: ContaminationOutputFiles,
        policy: ContaminationPolicy,
    ) -> GatingRecommendation:
        """Formulates actionable gating constraints for Block 6 following a strict safety hierarchy."""
        # Tier 1: Dati Insufficienti o Modello Inaffidabile
        if category == ContaminationCategory.INSUFFICIENT_DATA or confidence == EstimationConfidence.UNRELIABLE:
            return GatingRecommendation(
                contamination_category=category,
                recommended_action=ConsensusGatingAction.EXCLUDE_OR_MANUAL_REVIEW,
                preferred_consensus_source=None,
                required_pmd_threshold=None,
                summary_message=(
                    "Contamination evaluation is unreliable or data is insufficient. "
                    "Manual inspection required before downstream variant calling."
                ),
            )

        # Tier 2: Modificatore di Bassa Confidenza (Si applica a TUTTI i metodi inclusi FastScreen)
        if confidence == EstimationConfidence.LOW:
            if category == ContaminationCategory.LOW:
                return GatingRecommendation(
                    contamination_category=category,
                    recommended_action=ConsensusGatingAction.RESTRICT_HETEROPLASMY,
                    preferred_consensus_source="STANDARD_CALLER",
                    required_pmd_threshold=None,
                    summary_message=(
                        f"Point estimate is low ({estimate.point_estimate:.2%}) but evaluation confidence is LOW. "
                        "Standard calling permitted with conservative suppression of low-frequency heteroplasmies."
                    ),
                )
            elif category == ContaminationCategory.MODERATE:
                return GatingRecommendation(
                    contamination_category=category,
                    recommended_action=ConsensusGatingAction.PMD_FILTERING_REQUIRED,
                    preferred_consensus_source="PMD_FILTERED_BAM",
                    required_pmd_threshold=policy.pmd_filtering_threshold,
                    summary_message=(
                        f"Moderate contamination ({estimate.point_estimate:.2%}) with LOW confidence. "
                        f"Strict read filtering (PMD >= {policy.pmd_filtering_threshold}) recommended."
                    ),
                )
            else:
                return GatingRecommendation(
                    contamination_category=category,
                    recommended_action=ConsensusGatingAction.EXCLUDE_OR_MANUAL_REVIEW,
                    preferred_consensus_source=None,
                    required_pmd_threshold=policy.pmd_filtering_threshold,
                    summary_message=(
                        f"High contamination ({estimate.point_estimate:.2%}) combined with LOW confidence. "
                        "Sample flagged for exclusion or expert manual review."
                    ),
                )

        # Tier 3: Pre-Screening Fast Triage con Alta/Moderata Confidenza (FAST_CLEAR)
        if estimate.method_used == ContaminationMethod.FAST_SCREEN_TRIAGE and category == ContaminationCategory.LOW:
            return GatingRecommendation(
                contamination_category=category,
                recommended_action=ConsensusGatingAction.STANDARD_CALLING,
                preferred_consensus_source="STANDARD_CALLER",
                required_pmd_threshold=None,
                summary_message=(
                    "Fast pre-screening detected no significant allelic dissent above background noise floor "
                    "(screened without full Bayesian quantification). Sample cleared as bypass candidate."
                ),
            )

        # Tier 4: Percorsi Standard ad Alta/Moderata Confidenza (Track A / Schmutzi)
        if category == ContaminationCategory.LOW:
            return GatingRecommendation(
                contamination_category=category,
                recommended_action=ConsensusGatingAction.STANDARD_CALLING,
                preferred_consensus_source="STANDARD_CALLER",
                required_pmd_threshold=None,
                summary_message=(
                    f"Low contamination fraction ({estimate.point_estimate:.2%}). "
                    "Sample cleared for unrestricted standard consensus and heteroplasmy analysis."
                ),
            )

        if category == ContaminationCategory.MODERATE:
            return GatingRecommendation(
                contamination_category=category,
                recommended_action=ConsensusGatingAction.RESTRICT_HETEROPLASMY,
                preferred_consensus_source="STANDARD_CALLER",
                required_pmd_threshold=None,
                summary_message=(
                    f"Moderate contamination fraction ({estimate.point_estimate:.2%}). "
                    "Consensus haplotyping is valid; low-frequency heteroplasmy calling must be suppressed."
                ),
            )

        if category == ContaminationCategory.ELEVATED:
            has_valid_endo_fa = outputs.endogenous_fasta is not None and outputs.endogenous_fasta.is_file()
            if estimate.method_used == ContaminationMethod.SCHMUTZI and has_valid_endo_fa:
                return GatingRecommendation(
                    contamination_category=category,
                    recommended_action=ConsensusGatingAction.USE_ENDOGENOUS_SCHMUTZI_CONSENSUS,
                    preferred_consensus_source="SCHMUTZI_ENDOGENOUS_FASTA",
                    required_pmd_threshold=None,
                    summary_message=(
                        f"Elevated contamination ({estimate.point_estimate:.2%}). "
                        "Block 6 should adopt the Schmutzi-predicted endogenous consensus FASTA (Bayesian reconstruction)."
                    ),
                )
            else:
                return GatingRecommendation(
                    contamination_category=category,
                    recommended_action=ConsensusGatingAction.PMD_FILTERING_REQUIRED,
                    preferred_consensus_source="PMD_FILTERED_BAM",
                    required_pmd_threshold=policy.pmd_filtering_threshold,
                    summary_message=(
                        f"Elevated contamination ({estimate.point_estimate:.2%}). "
                        f"Downstream caller must restrict input reads to damage-verified BAM (PMD >= {policy.pmd_filtering_threshold})."
                    ),
                )

        return GatingRecommendation(
            contamination_category=category,
            recommended_action=ConsensusGatingAction.EXCLUDE_OR_MANUAL_REVIEW,
            preferred_consensus_source=None,
            required_pmd_threshold=policy.pmd_filtering_threshold,
            summary_message=(
                f"Critical contamination detected ({estimate.point_estimate:.2%}). "
                "Exogenous modern sequences dominate the library; sample flagged for exclusion or per-read isolation."
            ),
        )

    @classmethod
    def determine_execution_status(
        cls,
        error_codes: List[ContaminationErrorCode],
        warning_codes: List[ContaminationWarningCode],
        confidence: EstimationConfidence,
        category: ContaminationCategory,
        gating_action: ConsensusGatingAction,
    ) -> ContaminationStatus:
        """Determines the overall execution status enforcing strict coherence with gating actions."""
        if error_codes:
            return ContaminationStatus.FAIL

        if gating_action in (
            ConsensusGatingAction.EXCLUDE_OR_MANUAL_REVIEW,
            ConsensusGatingAction.PMD_FILTERING_REQUIRED,
        ):
            return ContaminationStatus.WARN

        if (
            warning_codes
            or confidence in (EstimationConfidence.LOW, EstimationConfidence.UNRELIABLE)
            or category in (ContaminationCategory.CRITICAL, ContaminationCategory.INSUFFICIENT_DATA)
        ):
            return ContaminationStatus.WARN

        return ContaminationStatus.PASS

    @classmethod
    def serialize_deliverable_json(
        cls,
        result: ContaminationResult,
        output_json_dir: Path,
    ) -> Optional[Path]:
        """Serializes the full ContaminationResult into a structured JSON deliverable."""
        try:
            output_json_dir.mkdir(parents=True, exist_ok=True)
            json_file = output_json_dir / f"{result.sample_id}.contamination.json"
            with open(json_file, "w", encoding="utf-8") as jf:
                json.dump(result.to_dict(), jf, indent=2, default=_contamination_json_serializer)
            logger.info(f"[{result.sample_id}] Exported contamination JSON -> {json_file.name}")
            return json_file
        except Exception as exc:
            logger.warning(f"[{result.sample_id}] Failed serializing contamination JSON: {exc}")
            return None

    @classmethod
    def evaluate_sample(
        cls,
        sample_id: str,
        estimate: ContaminationEstimate,
        diagnostics: ModelQualityDiagnostics,
        outputs: ContaminationOutputFiles,
        audit: ContaminationAuditMetadata,
        policy: ContaminationPolicy,
        fallback_reason: Optional[str] = None,
        output_json_dir: Optional[Path] = None,
    ) -> ContaminationResult:
        """Orchestrates quality evaluation, gating formulation, and JSON deliverable export."""
        all_codes: List[ContaminationWarningCode] = []
        all_warnings: List[str] = []
        error_codes: List[ContaminationErrorCode] = []
        errors: List[str] = []

        # 1. Fatal Data Integrity Checks
        if estimate.method_used == ContaminationMethod.NONE:
            error_codes.append(ContaminationErrorCode.CORRUPTED_OUTPUT)
            errors.append(f"Contamination analysis method is NONE for sample '{sample_id}'.")

        # 2. Aggiornamento Flag di Audit Metodologico
        if estimate.method_used == ContaminationMethod.SCHMUTZI:
            audit.full_bayesian_estimation_performed = True
            audit.screened_without_full_quantification = False
        elif estimate.method_used == ContaminationMethod.FAST_SCREEN_TRIAGE:
            audit.full_bayesian_estimation_performed = False
            audit.screened_without_full_quantification = True

        # 3. Fast Screening Audit & Warning Extraction
        if audit.fast_screen_applied and audit.fast_screen_result:
            fs_res = audit.fast_screen_result
            if fs_res.decision == ScreeningDecision.ESCALATE:
                all_codes.append(ContaminationWarningCode.FAST_SCREEN_ESCALATION)
                all_warnings.append(
                    f"FastScreen detected allelic dissent ({fs_res.dissenting_positions_count} sites, "
                    f"max VAF {fs_res.max_dissent_vaf:.2%}); escalated to Bayesian model."
                )
            elif fs_res.decision == ScreeningDecision.INSUFFICIENT_DATA:
                all_codes.append(ContaminationWarningCode.FAST_SCREEN_INSUFFICIENT_DATA)
                all_warnings.append(
                    f"FastScreen had insufficient diagnostic depth ({fs_res.mean_diagnostic_depth:.1f}x); "
                    "escalated to full Bayesian model for safety."
                )

        # 4. Model Confidence Evaluation
        confidence, conf_codes, conf_warnings = cls.evaluate_model_confidence(
            estimate=estimate,
            diagnostics=diagnostics,
            policy=policy,
        )
        diagnostics.estimation_confidence = confidence
        all_codes.extend(conf_codes)
        all_warnings.extend(conf_warnings)

        # 5. Categorization (Point Estimate + Multi-Interval Credible Bounds)
        category, cat_codes, cat_warnings = cls.categorize_contamination(
            estimate=estimate,
            confidence=confidence,
            diagnostics=diagnostics,
            policy=policy,
        )
        all_codes.extend(cat_codes)
        all_warnings.extend(cat_warnings)

        # 6. Gating Recommendation Formulation
        gating = cls.formulate_gating_recommendation(
            category=category,
            confidence=confidence,
            estimate=estimate,
            outputs=outputs,
            policy=policy,
        )

        # 7. Deduplicazione Codici e Messaggi Diagnostici
        unique_warning_codes = cls._deduplicate_items(all_codes)
        unique_warnings = cls._deduplicate_items(all_warnings)
        unique_error_codes = cls._deduplicate_items(error_codes)
        unique_errors = cls._deduplicate_items(errors)

        # 8. Execution Status Determination
        status = cls.determine_execution_status(
            error_codes=unique_error_codes,
            warning_codes=unique_warning_codes,
            confidence=confidence,
            category=category,
            gating_action=gating.recommended_action,
        )

        result = ContaminationResult(
            sample_id=sample_id,
            status=status,
            policy=policy,
            estimate=estimate,
            diagnostics=diagnostics,
            gating=gating,
            outputs=outputs,
            audit=audit,
            fallback_reason=fallback_reason,
            warning_codes=unique_warning_codes,
            error_codes=unique_error_codes,
            warnings=unique_warnings,
            errors=unique_errors,
        )

        # 9. Esportazione JSON Deliverable con Gestione Fallimento Audit
        if output_json_dir:
            json_path = cls.serialize_deliverable_json(result, output_json_dir)
            if json_path:
                result.outputs.contamination_json = json_path
            else:
                result.add_warning(
                    ContaminationWarningCode.LOW_INFORMATIVE_SITES,
                    f"Failed writing contamination JSON deliverable to {output_json_dir}."
                )

        logger.info(
            f"[{sample_id}] Contamination QC Finished | Method: {estimate.method_used.value} | "
            f"Est: {estimate.point_estimate:.2%} | Category: {category.value} | "
            f"Confidence: {confidence.value} | Gating Action: {gating.recommended_action.value} | "
            f"Status: {status.value}"
        )

        return result