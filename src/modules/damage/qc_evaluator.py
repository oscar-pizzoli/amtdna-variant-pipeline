"""Damage Quality Control & Authenticity Evaluation Engine.

Synthesizes empirical deamination rates, Bayesian model parameters, and PMD distributions
into structured DamageModelProfile and DamageResult deliverables with strict protocol-aware validation.
"""

from enum import Enum
import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from src.modules.damage.models import (
    AuthenticationVerdict,
    BayesianDamageParameters,
    DamageErrorCode,
    DamageModelProfile,
    DamageOutputFiles,
    DamagePolicy,
    DamageResult,
    DamageSignalLevel,
    DamageStatus,
    DamageWarningCode,
    FragmentLengthMetrics,
    LibraryTreatment,
    PMDMetrics,
    TerminalDamageFrequencies,
)

logger = logging.getLogger("mtDNA_Pipeline.DamageQCEvaluator")


def _damage_json_serializer(obj: Any) -> Any:
    """Robust JSON encoder helper ensuring clean Enum.value, string Path, and numpy scalar exports."""
    if isinstance(obj, Enum):
        return obj.value
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


class DamageQCValidator:
    """Evaluates multi-criteria damage evidence and assigns protocol-aware authenticity verdicts."""

    @classmethod
    def evaluate_bayesian_convergence(
        cls,
        bayesian: BayesianDamageParameters,
    ) -> Tuple[bool, List[DamageWarningCode], List[str]]:
        """Audits mapDamage MCMC Bayesian parameter convergence prior to downstream usage.

        Args:
            bayesian: Parsed Bayesian parameters.

        Returns:
            Tuple of (is_usable, warning_codes, warning_messages).
        """
        warnings: List[str] = []
        codes: List[DamageWarningCode] = []

        # Case 1: Fitting failed to run or generate parameters
        if bayesian.lambda_param is None or bayesian.delta_s is None or bayesian.delta_d is None:
            codes.append(DamageWarningCode.BAYESIAN_CONVERGENCE_FAILURE)
            warnings.append(
                "mapDamage Bayesian statistical fitting failed or produced incomplete parameter estimates."
            )
            return False, codes, warnings

        # Case 2: Parameters exist but MCMC did not achieve nominal convergence
        if not bayesian.model_converged:
            codes.append(DamageWarningCode.BAYESIAN_CONVERGENCE_FAILURE)
            warnings.append(
                "mapDamage Bayesian MCMC statistical fit did not achieve nominal convergence; "
                "parameters excluded from primary signal weighting."
            )
            return False, codes, warnings

        return True, codes, warnings

    @classmethod
    def evaluate_damage_signal_level(
        cls,
        terminal: TerminalDamageFrequencies,
        bayesian: BayesianDamageParameters,
        policy: DamagePolicy,
        bayesian_usable: bool,
    ) -> Tuple[DamageSignalLevel, List[DamageWarningCode], List[str]]:
        """Determines empirical deamination signal strength supported by rigorous statistical denominators.

        Args:
            terminal: Terminal substitution rates and observation counts.
            bayesian: Bayesian model parameter container.
            policy: Configured damage thresholds and policy parameters.
            bayesian_usable: Flag indicating whether MCMC parameters are converged and trustworthy.

        Returns:
            Tuple of (DamageSignalLevel, warning_codes, warning_messages).
        """
        warnings: List[str] = []
        codes: List[DamageWarningCode] = []

        is_ss = policy.library_treatment == LibraryTreatment.SINGLE_STRANDED
        is_udg_full = policy.library_treatment == LibraryTreatment.UDG_FULL
        is_udg_partial = policy.library_treatment == LibraryTreatment.UDG_PARTIAL

        # 1. Statistical power check on terminal denominators
        has_sufficient_5p = (
            terminal.n_c_ref_5p_pos1_total >= policy.min_c_to_t_5p_observations
            and terminal.n_c_ref_5p_pos1_total > 0
        )

        has_sufficient_3p = True
        if not is_ss and not is_udg_full and not is_udg_partial:
            has_sufficient_3p = (
                terminal.n_g_ref_3p_pos1_total >= policy.min_g_to_a_3p_observations
                and terminal.n_g_ref_3p_pos1_total > 0
            )

        if not has_sufficient_5p or not has_sufficient_3p:
            codes.append(DamageWarningCode.INSUFFICIENT_TERMINAL_COVERAGE)
            warnings.append(
                f"Insufficient terminal observations (5' C={terminal.n_c_ref_5p_pos1_total}/"
                f"{policy.min_c_to_t_5p_observations}, 3' G={terminal.n_g_ref_3p_pos1_total}/"
                f"{policy.min_g_to_a_3p_observations}). Statistical power is insufficient for reliable damage estimation."
            )
            return DamageSignalLevel.NOISY_OR_INSUFFICIENT, codes, warnings

        c_to_t_5p = terminal.c_to_t_5p_pos1
        g_to_a_3p = terminal.g_to_a_3p_pos1
        delta_s = bayesian.delta_s if (bayesian_usable and bayesian.delta_s is not None) else None

        # 2. Terminal asymmetry check for standard double-stranded libraries
        if not is_ss and not is_udg_full and not is_udg_partial:
            if abs(c_to_t_5p - g_to_a_3p) > policy.terminal_asymmetry_max_difference:
                codes.append(DamageWarningCode.ASYMMETRIC_TERMINAL_DAMAGE)
                warnings.append(
                    f"Pronounced asymmetry detected between 5' C->T ({c_to_t_5p:.2%}) and 3' G->A ({g_to_a_3p:.2%})."
                )

        # 3. Categorize empirical signal level
        is_strong_both = (
            c_to_t_5p >= policy.min_c_to_t_5p_strong_threshold
            and g_to_a_3p >= policy.min_g_to_a_3p_strong_threshold
        )
        is_strong_5p_single = (is_ss or is_udg_partial) and (c_to_t_5p >= policy.min_c_to_t_5p_strong_threshold)
        is_strong_delta = delta_s is not None and delta_s >= policy.bayesian_delta_s_strong_threshold

        if is_strong_both or is_strong_5p_single or is_strong_delta:
            return DamageSignalLevel.STRONG_SIGNAL, codes, warnings

        # Moderate Signal: Measurable terminal deamination (>= 3%) at 5' or 3' end, or delta_S >= 0.05
        is_moderate_5p = c_to_t_5p >= policy.min_c_to_t_5p_moderate_threshold
        is_moderate_3p = g_to_a_3p >= policy.min_g_to_a_3p_moderate_threshold
        is_moderate_delta = delta_s is not None and delta_s >= policy.bayesian_delta_s_moderate_threshold
        is_elevated_5p = c_to_t_5p >= policy.min_c_to_t_5p_strong_threshold

        if is_moderate_5p or is_moderate_3p or is_moderate_delta or is_elevated_5p:
            return DamageSignalLevel.MODERATE_SIGNAL, codes, warnings

        return DamageSignalLevel.LOW_SIGNAL, codes, warnings

    @classmethod
    def determine_authenticity_verdict(
        cls,
        signal_level: DamageSignalLevel,
        terminal: TerminalDamageFrequencies,
        length_metrics: FragmentLengthMetrics,
        pmd_metrics: Optional[PMDMetrics],
        policy: DamagePolicy,
    ) -> Tuple[AuthenticationVerdict, List[DamageWarningCode], List[str]]:
        """Assigns protocol-aware authenticity verdict across all experimental library types."""
        warnings: List[str] = []
        codes: List[DamageWarningCode] = []

        # Universal Guardrail: If statistical power was insufficient, all protocols return INCONCLUSIVE
        if signal_level == DamageSignalLevel.NOISY_OR_INSUFFICIENT:
            return AuthenticationVerdict.INCONCLUSIVE, codes, warnings

        treatment = policy.library_treatment

        # Audit check: Unknown library preparation
        if treatment == LibraryTreatment.UNKNOWN:
            codes.append(DamageWarningCode.UNKNOWN_LIBRARY_TREATMENT)
            warnings.append("Library treatment protocol is UNKNOWN. Authentic ancient DNA evaluation is ambiguous.")
            return AuthenticationVerdict.INCONCLUSIVE, codes, warnings

        # Protocol 1: Full UDG/USER treatment (UDG_FULL)
        if treatment == LibraryTreatment.UDG_FULL:
            if signal_level in (DamageSignalLevel.STRONG_SIGNAL, DamageSignalLevel.MODERATE_SIGNAL):
                codes.append(DamageWarningCode.UNEXPECTED_HIGH_DAMAGE_ON_UDG)
                warnings.append(
                    f"Unexpectedly high terminal C->T damage ({terminal.c_to_t_5p_pos1:.2%}) on a declared UDG_FULL library."
                )
                return AuthenticationVerdict.INCONCLUSIVE, codes, warnings
            return AuthenticationVerdict.UDG_DEPLETED_EXPECTED, codes, warnings

        # Protocol 2: Partial UDG treatment (UDG_PARTIAL / UDG-half)
        if treatment == LibraryTreatment.UDG_PARTIAL:
            if signal_level in (DamageSignalLevel.STRONG_SIGNAL, DamageSignalLevel.MODERATE_SIGNAL):
                # Verify that damage is restricted to pos 1 and depleted at pos 2
                pos1 = terminal.c_to_t_5p_pos1
                pos2 = terminal.c_to_t_5p_pos2
                if (
                    pos1 > 0
                    and pos2 > (pos1 * policy.udg_partial_max_pos2_to_pos1_ratio)
                    and pos2 > policy.udg_partial_min_pos2_damage_floor
                ):
                    codes.append(DamageWarningCode.ASYMMETRIC_TERMINAL_DAMAGE)
                    warnings.append(
                        f"UDG_PARTIAL library exhibits elevated internal deamination at pos 2 ({pos2:.2%} vs pos 1 {pos1:.2%})."
                    )
                return AuthenticationVerdict.AUTHENTIC_ANCIENT_SUPPORTED, codes, warnings
            return AuthenticationVerdict.AUTHENTICITY_UNSUPPORTED, codes, warnings

        # Protocol 3: Single-stranded DNA library (ssDNA)
        if treatment == LibraryTreatment.SINGLE_STRANDED:
            if signal_level == DamageSignalLevel.STRONG_SIGNAL:
                return AuthenticationVerdict.AUTHENTIC_ANCIENT_SUPPORTED, codes, warnings
            if signal_level == DamageSignalLevel.MODERATE_SIGNAL:
                if (
                    length_metrics.median_length is not None
                    and 0 < length_metrics.median_length <= policy.max_median_length_threshold_bp
                ):
                    return AuthenticationVerdict.AUTHENTIC_ANCIENT_SUPPORTED, codes, warnings
                return AuthenticationVerdict.INCONCLUSIVE, codes, warnings
            return AuthenticationVerdict.AUTHENTICITY_UNSUPPORTED, codes, warnings

        # Protocol 4: Standard Non-UDG Double-Stranded DNA library (NONE)
        if signal_level == DamageSignalLevel.STRONG_SIGNAL:
            return AuthenticationVerdict.AUTHENTIC_ANCIENT_SUPPORTED, codes, warnings

        if signal_level == DamageSignalLevel.MODERATE_SIGNAL:
            has_short_fragments = (
                length_metrics.median_length is not None
                and 0 < length_metrics.median_length <= policy.max_median_length_threshold_bp
            )
            has_positive_pmd = (
                pmd_metrics is not None
                and pmd_metrics.fraction_pmd_ge_3 >= policy.pmd_corroboration_fraction_ge_3
            )

            if has_short_fragments or has_positive_pmd:
                return AuthenticationVerdict.AUTHENTIC_ANCIENT_SUPPORTED, codes, warnings
            return AuthenticationVerdict.INCONCLUSIVE, codes, warnings

        if signal_level == DamageSignalLevel.LOW_SIGNAL:
            codes.append(DamageWarningCode.LOW_DAMAGE_ON_NON_UDG)
            warnings.append(
                f"Low terminal damage ({terminal.c_to_t_5p_pos1:.2%}) observed on declared {treatment.value} library. "
                "Evidence is insufficient to confirm ancient authenticity."
            )
            return AuthenticationVerdict.AUTHENTICITY_UNSUPPORTED, codes, warnings

        return AuthenticationVerdict.INCONCLUSIVE, codes, warnings

    @classmethod
    def evaluate_sample(
        cls,
        sample_id: str,
        input_bam_path: Path,
        reference_contig: str,
        terminal_damage: TerminalDamageFrequencies,
        fragment_length: FragmentLengthMetrics,
        bayesian_params: BayesianDamageParameters,
        pmd_metrics: Optional[PMDMetrics],
        outputs: DamageOutputFiles,
        policy: DamagePolicy,
        rescaling_applied: bool = False,
        output_json_dir: Optional[Path] = None,
    ) -> DamageResult:
        """Synthesizes Block 4 outputs into DamageModelProfile, DamageResult, and exports damage_model.json."""
        all_codes: List[DamageWarningCode] = []
        all_warnings: List[str] = []
        error_codes: List[DamageErrorCode] = []
        errors: List[str] = []

        # 1. Fatal data integrity checks
        if fragment_length.total_fragments_sampled == 0 and (
            terminal_damage.n_c_ref_5p_pos1_total == 0 and terminal_damage.n_g_ref_3p_pos1_total == 0
        ):
            error_codes.append(DamageErrorCode.CORRUPTED_DAMAGE_OUTPUT)
            errors.append(f"Zero fragments or terminal bases sampled for sample '{sample_id}'. Damage outputs are empty.")

        # 2. Audit Bayesian fitting convergence before signal weighting
        bayesian_usable, bay_codes, bay_warnings = cls.evaluate_bayesian_convergence(bayesian_params)
        all_codes.extend(bay_codes)
        all_warnings.extend(bay_warnings)

        # 3. Evaluate empirical damage signal level
        signal_level, sig_codes, sig_warnings = cls.evaluate_damage_signal_level(
            terminal=terminal_damage,
            bayesian=bayesian_params,
            policy=policy,
            bayesian_usable=bayesian_usable,
        )
        all_codes.extend(sig_codes)
        all_warnings.extend(sig_warnings)

        # 4. Assign protocol-aware authenticity verdict
        verdict, ver_codes, ver_warnings = cls.determine_authenticity_verdict(
            signal_level=signal_level,
            terminal=terminal_damage,
            length_metrics=fragment_length,
            pmd_metrics=pmd_metrics,
            policy=policy,
        )
        all_codes.extend(ver_codes)
        all_warnings.extend(ver_warnings)

        # 5. Build consolidated DamageModelProfile
        profile = DamageModelProfile(
            sample_id=sample_id,
            input_bam_type=policy.input_bam_type,
            input_bam_path=input_bam_path,
            reference_contig=reference_contig,
            library_treatment=policy.library_treatment,
            terminal_damage=terminal_damage,
            fragment_length=fragment_length,
            bayesian_parameters=bayesian_params,
            pmd_metrics=pmd_metrics,
            damage_signal_level=signal_level,
            authentication_verdict=verdict,
            rescaling_applied=rescaling_applied,
        )

        # 6. Determine final execution status
        status = DamageStatus.PASS
        if error_codes or errors:
            status = DamageStatus.FAIL
        elif all_codes or all_warnings:
            status = DamageStatus.WARN

        # 7. Clean JSON artifact export with custom encoder
        if output_json_dir:
            output_json_dir.mkdir(parents=True, exist_ok=True)
            json_model_file = output_json_dir / f"{sample_id}.damage_model.json"
            json_legacy_file = output_json_dir / f"{sample_id}.damage.json"
            payload = profile.to_dict()
            try:
                for jf_path in (json_model_file, json_legacy_file):
                    with open(jf_path, "w", encoding="utf-8") as jf:
                        json.dump(payload, jf, indent=2, default=_damage_json_serializer)
                outputs.damage_model_json = json_model_file
                logger.info(f"[{sample_id}] Exported damage model profiles -> {json_model_file.name}, {json_legacy_file.name}")
            except Exception as exc:
                logger.warning(f"[{sample_id}] Failed to serialize damage model JSON: {exc}")

        result = DamageResult(
            sample_id=sample_id,
            status=status,
            policy=policy,
            profile=profile,
            outputs=outputs,
            warning_codes=all_codes,
            error_codes=error_codes,
            warnings=all_warnings,
            errors=errors,
        )

        logger.info(
            f"[{sample_id}] QC Evaluation Finished | Signal: {signal_level.value} | "
            f"Verdict: {verdict.value} | Status: {status.value}"
        )

        return result