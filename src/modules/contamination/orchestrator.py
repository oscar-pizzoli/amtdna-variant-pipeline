"""Orchestrator for Block 5: mtDNA Contamination Estimation & Two-Tier Gating.

Implements conservative two-tier routing:
  - Tier 1: Pre-Screening Triage (FastContaminationScreen, FAST_CLEAR bypass)
  - Tier 2: Subsampling, samtools calmd MD-tag injection & Full Bayesian / Dissent Modeling
  - Rigorous pre-flight checks, observed subsampling depth tracking, and audit trail
  - QC validation, gating action synthesis, and synchronized JSON delivery
"""

import hashlib
import json
import logging
from pathlib import Path
import shutil
import time
from typing import Any, Dict, List, Optional, Tuple

import pysam
import yaml

from src.modules.alignment.models import ReferenceInfo
from src.modules.contamination.dissent_engine import DissentEngine
from src.modules.contamination.fast_screen import (
    FastContaminationScreen,
    FastScreenExecutionError,
)
from src.modules.contamination.input_preparer import (
    SchmutziInputPreparer,
    SchmutziPreparationError,
)
from src.modules.contamination.models import (
    ConsensusGatingAction,
    ContaminationAuditMetadata,
    ContaminationCategory,
    ContaminationErrorCode,
    ContaminationEstimate,
    ContaminationInputContext,
    ContaminationMethod,
    ContaminationOutputFiles,
    ContaminationPolicy,
    ContaminationResult,
    ContaminationStatus,
    ContaminationWarningCode,
    DissentPolicy,
    EstimationConfidence,
    FastScreeningPolicy,
    FastScreenResult,
    GatingRecommendation,
    LibraryTreatment,
    ModelQualityDiagnostics,
    SchmutziPolicy,
    ScreeningDecision,
    SubsamplingPolicy,
)
from src.modules.contamination.qc_evaluator import ContaminationQCValidator
from src.modules.contamination.schmutzi_engine import (
    SchmutziEngine,
    SchmutziExecutionError,
)

logger = logging.getLogger("mtDNA_Pipeline.ContaminationOrchestrator")

PROJECT_ROOT = Path(__file__).resolve().parents[3]


class ContaminationOrchestrationError(Exception):
    """Raised when contamination workflow execution, pre-flight validation, or pipeline stages fail fatally."""
    pass


class ContaminationOrchestrator:
    """Orchestrates mtDNA contamination estimation, pre-screening triage, and downstream gating."""

    def __init__(
        self,
        config: Optional[Dict[str, Any]] = None,
        config_path: Optional[Path] = None,
        project_root: Optional[Path] = None,
    ) -> None:
        self.project_root = project_root or PROJECT_ROOT
        self.config = config or self._load_config(config_path)
        self.policy = self._build_policy()

    def _load_config(self, config_path: Optional[Path]) -> Dict[str, Any]:
        """Loads configuration from YAML file."""
        target_path = config_path or (self.project_root / "config" / "default_config.yaml")
        if not target_path.is_file():
            logger.warning(f"Config file not found at {target_path}, using built-in defaults.")
            return {}
        try:
            with open(target_path, "r", encoding="utf-8") as f:
                return yaml.safe_load(f) or {}
        except Exception as exc:
            logger.warning(f"Failed loading YAML config from {target_path}: {exc}")
            return {}

    def _resolve_path(self, path_str: Optional[str], default_rel: Path) -> Path:
        """Resolves path relative to project_root."""
        if not path_str:
            return default_rel if default_rel.is_absolute() else (self.project_root / default_rel)
        p = Path(path_str)
        return p if p.is_absolute() else (self.project_root / p)

    def _build_policy(self) -> ContaminationPolicy:
        """Constructs and validates ContaminationPolicy from configuration dictionary."""
        contam_cfg = self.config.get("contamination", {})
        schmutzi_cfg = contam_cfg.get("schmutzi", {})
        dissent_cfg = contam_cfg.get("dissent", {})
        subsample_cfg = contam_cfg.get("subsampling", {})
        fast_cfg = contam_cfg.get("fast_screening", {})
        thresh_cfg = contam_cfg.get("thresholds", {})

        db_path_str = schmutzi_cfg.get("database_path", "data/reference/schmutzi_db/197")
        db_path = self._resolve_path(db_path_str, self.project_root / "data" / "reference" / "schmutzi_db" / "197")

        panel_path_str = fast_cfg.get("diagnostic_panel_path", "data/reference/diagnostic_sites_phylotree17.tsv")
        panel_path = self._resolve_path(panel_path_str, self.project_root / "data" / "reference" / "diagnostic_sites_phylotree17.tsv")

        lib_treatment_val = schmutzi_cfg.get("library_treatment", "NONE")
        try:
            lib_treatment = LibraryTreatment(lib_treatment_val)
        except Exception:
            lib_treatment = LibraryTreatment.NONE

        return ContaminationPolicy(
            policy_id="contamination_standard_v2",
            method=ContaminationMethod(contam_cfg.get("method", "SCHMUTZI")),
            enabled=contam_cfg.get("enabled", True),
            allow_fallback_to_dissent=contam_cfg.get("allow_fallback_to_dissent", False),
            schmutzi=SchmutziPolicy(
                contaminant_database=db_path,
                threads=schmutzi_cfg.get("threads", 4),
                max_iterations=schmutzi_cfg.get("max_iterations", 10),
                use_wrapped_reference=schmutzi_cfg.get("use_wrapped_reference", False),
                library_treatment=lib_treatment,
            ),
            dissent=DissentPolicy(
                min_base_quality=dissent_cfg.get("min_base_quality", 30),
                min_mapping_quality=dissent_cfg.get("min_mapping_quality", 30),
                min_depth_per_site=dissent_cfg.get("min_depth_per_site", 10),
                noise_vaf_floor=dissent_cfg.get("noise_vaf_floor", 0.015),
                base_min_dissenting_sites=dissent_cfg.get("base_min_dissenting_sites", 3),
                dissent_high_confidence_min_sites=dissent_cfg.get("dissent_high_confidence_min_sites", 50),
            ),
            subsampling=SubsamplingPolicy(
                enabled=subsample_cfg.get("enabled", True),
                target_depth=subsample_cfg.get("target_depth", 40.0),
                seed=subsample_cfg.get("seed", 42),
                min_depth_threshold=subsample_cfg.get("min_depth_threshold", 42.0),
            ),
            fast_screening=FastScreeningPolicy(
                enabled=fast_cfg.get("enabled", False),
                algorithm_version=fast_cfg.get("algorithm_version", "fast-screen-v2.0"),
                diagnostic_panel_path=panel_path,
                diagnostic_panel_version=fast_cfg.get("diagnostic_panel_version", "phylotree17_v1"),
                escalate_dissenting_sites_min=fast_cfg.get("escalate_dissenting_sites_min", 3),
                escalate_vaf_threshold=fast_cfg.get("escalate_vaf_threshold", 0.02),
                min_alt_reads_per_site=fast_cfg.get("min_alt_reads_per_site", 2),
                min_depth_for_screen=fast_cfg.get("min_depth_for_screen", 10.0),
                min_genome_depth=fast_cfg.get("min_genome_depth", 10.0),
                min_evaluated_sites=fast_cfg.get("min_evaluated_sites", 50),
                min_panel_evaluation_fraction=fast_cfg.get("min_panel_evaluation_fraction", 0.25),
                min_base_quality=fast_cfg.get("min_base_quality", 30),
                min_mapping_quality=fast_cfg.get("min_mapping_quality", 30),
                min_depth_per_site=fast_cfg.get("min_depth_per_site", 10),
                noise_vaf_floor=fast_cfg.get("noise_vaf_floor", 0.015),
                cache_screen_results=fast_cfg.get("cache_screen_results", True),
                cache_mode=fast_cfg.get("cache_mode", "strict"),
            ),
            low_threshold=thresh_cfg.get("low_threshold", 0.02),
            moderate_threshold=thresh_cfg.get("moderate_threshold", 0.05),
            elevated_threshold=thresh_cfg.get("elevated_threshold", 0.10),
            min_coverage_for_estimation=thresh_cfg.get("min_coverage_for_estimation", 5.0),
            recommended_coverage=thresh_cfg.get("recommended_coverage", 10.0),
            high_confidence_max_ci_width=thresh_cfg.get("high_confidence_max_ci_width", 0.04),
            moderate_confidence_max_ci_width=thresh_cfg.get("moderate_confidence_max_ci_width", 0.08),
            min_informative_sites=thresh_cfg.get("min_informative_sites", 10),
            pmd_filtering_threshold=thresh_cfg.get("pmd_filtering_threshold", 3.0),
        )

    def _compute_sha256(self, file_path: Path) -> Optional[str]:
        """Calculates SHA-256 digest of a target file."""
        if not file_path.is_file():
            return None
        hasher = hashlib.sha256()
        with open(file_path, "rb") as f:
            while chunk := f.read(65536):
                hasher.update(chunk)
        return hasher.hexdigest()

    def _preflight_validation(self, context: ContaminationInputContext) -> Optional[Tuple[ContaminationErrorCode, str]]:
        """Performs structural integrity, contig, index, and header compatibility checks."""
        if not context.input_bam_path.is_file():
            return ContaminationErrorCode.INPUT_BAM_NOT_FOUND, f"Input BAM not found: {context.input_bam_path}"

        if not context.reference_fasta_path.is_file():
            return ContaminationErrorCode.REFERENCE_INCOMPATIBLE, f"Reference FASTA not found: {context.reference_fasta_path}"

        bai_path = context.input_bam_path.with_suffix(".bam.bai")
        csi_path = context.input_bam_path.with_suffix(".bam.csi")
        alt_bai = Path(f"{context.input_bam_path}.bai")
        if not (bai_path.is_file() or csi_path.is_file() or alt_bai.is_file()):
            logger.info(f"[{context.sample_id}] Index missing for BAM, indexing dynamically...")
            try:
                pysam.index(str(context.input_bam_path))
            except Exception as exc:
                return ContaminationErrorCode.CORRUPTED_OUTPUT, f"Failed creating BAM index: {exc}"

        try:
            with pysam.AlignmentFile(str(context.input_bam_path), "rb") as bam:
                if context.reference_contig not in bam.references:
                    return (
                        ContaminationErrorCode.REFERENCE_INCOMPATIBLE,
                        f"Target contig '{context.reference_contig}' missing from BAM headers: {bam.references}",
                    )
                bam_contig_len = bam.get_reference_length(context.reference_contig)
        except Exception as exc:
            return ContaminationErrorCode.CORRUPTED_OUTPUT, f"Input BAM header unreadable or corrupted: {exc}"

        try:
            with pysam.FastaFile(str(context.reference_fasta_path)) as fa:
                if context.reference_contig not in fa.references:
                    return (
                        ContaminationErrorCode.REFERENCE_INCOMPATIBLE,
                        f"Target contig '{context.reference_contig}' missing from FASTA: {fa.references}",
                    )
                fa_contig_len = fa.get_reference_length(context.reference_contig)
        except Exception as exc:
            return ContaminationErrorCode.REFERENCE_INCOMPATIBLE, f"Reference FASTA unreadable or missing .fai: {exc}"

        if bam_contig_len != fa_contig_len:
            return (
                ContaminationErrorCode.REFERENCE_INCOMPATIBLE,
                f"Contig length mismatch: BAM has {bam_contig_len} bp, FASTA has {fa_contig_len} bp for {context.reference_contig}",
            )

        if self.policy.method == ContaminationMethod.SCHMUTZI and context.is_ancient:
            db_dir = self.policy.schmutzi.contaminant_database
            if not db_dir or not db_dir.is_dir():
                return ContaminationErrorCode.DATABASE_MISSING, f"Schmutzi allele database not found at {db_dir}"

        return None

    def _save_contamination_json(self, result: ContaminationResult, output_path: Path) -> Path:
        """Serializes finalized ContaminationResult to disk ensuring complete timestamps."""
        output_path.parent.mkdir(parents=True, exist_ok=True)
        data = result.to_dict()
        with open(output_path, "w", encoding="utf-8") as f_json:
            json.dump(data, f_json, indent=2)
        return output_path

    def process_sample(
        self,
        context: ContaminationInputContext,
        output_dir: Path,
    ) -> ContaminationResult:
        """Executes complete contamination analysis adhering to two-tier screening architecture."""
        t_task_start = time.perf_counter()
        sample_id = context.sample_id
        sample_out_dir = output_dir / sample_id
        sample_out_dir.mkdir(parents=True, exist_ok=True)
        final_json_path = sample_out_dir / f"{sample_id}.contamination.json"

        result = ContaminationResult(sample_id=sample_id, policy=self.policy)
        result.mark_started()

        audit = ContaminationAuditMetadata(
            input_bam_path=context.input_bam_path,
            input_bam_sha256=self._compute_sha256(context.input_bam_path),
            is_unrescaled_bam_verified=context.is_unrescaled_bam_verified,
            reference_fasta_path=context.reference_fasta_path,
            reference_fasta_sha256=self._compute_sha256(context.reference_fasta_path),
        )
        outputs = ContaminationOutputFiles()

        preflight_error = self._preflight_validation(context)
        if preflight_error:
            err_code, err_msg = preflight_error
            logger.error(f"[{sample_id}] Pre-flight validation failed: {err_msg}")
            result.add_error(err_code, err_msg)
            result.audit = audit
            result.outputs = outputs
            result.mark_finished()
            result.runtime_seconds = round(time.perf_counter() - t_task_start, 2)
            self._save_contamination_json(result, final_json_path)
            return result

        # TIER 1: Pre-Screening Triage (Fast Screen)
        if self.policy.fast_screening.enabled:
            logger.info(f"[{sample_id}] Running Tier 1 Fast Contamination Screen...")
            t_fs0 = time.perf_counter()

            try:
                fs_res, dissent_tsv = FastContaminationScreen.run_screen(
                    sample_id=sample_id,
                    input_bam=context.input_bam_path,
                    reference_fasta=context.reference_fasta_path,
                    expected_contig=context.reference_contig,
                    policy=self.policy.fast_screening,
                    output_dir=sample_out_dir,
                )
            except Exception as exc:
                logger.error(f"[{sample_id}] FastScreen failed: {exc}")
                result.add_error(ContaminationErrorCode.FAST_SCREEN_FAILED, str(exc))
                result.audit = audit
                result.outputs = outputs
                result.mark_finished()
                result.runtime_seconds = round(time.perf_counter() - t_task_start, 2)
                self._save_contamination_json(result, final_json_path)
                return result

            t_fs = time.perf_counter() - t_fs0
            audit.fast_screen_applied = True
            audit.fast_screen_result = fs_res
            audit.step_runtimes_seconds["fast_screen"] = round(t_fs, 2)
            outputs.dissent_tsv = dissent_tsv

            if fs_res.decision == ScreeningDecision.FAST_CLEAR and fs_res.bypass_eligible:
                logger.info(f"[{sample_id}] FAST_CLEAR verified (0 discordant sites). Bypassing full Schmutzi MCMC.")
                audit.full_bayesian_estimation_performed = False
                audit.screened_without_full_quantification = True

                estimate = ContaminationEstimate(
                    point_estimate=0.0,
                    ci_lower_95=None,
                    ci_upper_95=None,
                    method_used=ContaminationMethod.FAST_SCREEN_TRIAGE,
                )

                diagnostics = ModelQualityDiagnostics(
                    mean_depth=fs_res.mean_diagnostic_depth,
                    informative_positions_evaluated=fs_res.evaluated_diagnostic_positions,
                    model_converged=False,
                    estimation_confidence=EstimationConfidence.LOW,
                )

                qc_result = ContaminationQCValidator.evaluate_sample(
                    sample_id=sample_id,
                    estimate=estimate,
                    diagnostics=diagnostics,
                    outputs=outputs,
                    audit=audit,
                    policy=self.policy,
                    output_json_dir=None,
                )

                qc_result.mark_finished()
                qc_result.runtime_seconds = round(time.perf_counter() - t_task_start, 2)
                qc_result.outputs.contamination_json = final_json_path
                self._save_contamination_json(qc_result, final_json_path)
                return qc_result

            elif fs_res.decision == ScreeningDecision.ESCALATE:
                logger.warning(
                    f"[{sample_id}] FAST_SCREEN escalated: {fs_res.dissenting_positions_count} dissenting sites detected "
                    f"(max VAF: {fs_res.max_dissent_vaf:.2%}). Routing to Phase 2 modeling."
                )
                result.add_warning(ContaminationWarningCode.FAST_SCREEN_ESCALATION, fs_res.interpretation)

            elif fs_res.decision == ScreeningDecision.INSUFFICIENT_DATA:
                logger.warning(f"[{sample_id}] FAST_SCREEN insufficient data: {fs_res.interpretation}. Routing to Phase 2 modeling.")
                result.add_warning(ContaminationWarningCode.FAST_SCREEN_INSUFFICIENT_DATA, fs_res.interpretation)

        # TIER 2: Modeling (Schmutzi for Ancient mtDNA, Direct Dissent for Modern mtDNA)
        is_ancient = getattr(context, "is_ancient", True)
        logger.info(
            f"[{sample_id}] Running Tier 2 Modeling "
            f"({'Schmutzi Bayesian MCMC' if is_ancient else 'Direct DissentEngine'})..."
        )
        audit.full_bayesian_estimation_performed = is_ancient
        audit.screened_without_full_quantification = False

        ref_info = ReferenceInfo(
            fasta_path=context.reference_fasta_path,
            contig_name=context.reference_contig,
            length_bp=context.reference_length_bp,
        )

        try:
            (
                schmutzi_ready_bam,
                schmutzi_ready_bai,
                was_subsampled,
                orig_depth,
                post_depth,
            ) = SchmutziInputPreparer.prepare_bam_for_schmutzi(
                sample_id=sample_id,
                input_bam=context.input_bam_path,
                reference_info=ref_info,
                work_dir=sample_out_dir,
                subsampling_policy=self.policy.subsampling,
                threads=self.policy.schmutzi.threads,
                allow_non_mt_contigs=True,
            )
            effective_bam = schmutzi_ready_bam
            outputs.schmutzi_ready_bam = schmutzi_ready_bam
            outputs.schmutzi_ready_bai = schmutzi_ready_bai

            if was_subsampled:
                audit.subsampling_applied = True
                audit.original_mean_depth = round(orig_depth, 2)
                audit.downsampled_mean_depth = round(post_depth if post_depth is not None else orig_depth, 2)
                if orig_depth > 0:
                    audit.subsampling_fraction = round(
                        self.policy.subsampling.target_depth / orig_depth, 4
                    )
        except Exception as prep_exc:
            result.add_error(
                ContaminationErrorCode.PREPARER_FAILED,
                f"BAM preparation with samtools calmd failed: {prep_exc}",
            )
            result.audit = audit
            result.outputs = outputs
            result.mark_finished()
            result.runtime_seconds = round(time.perf_counter() - t_task_start, 2)
            self._save_contamination_json(result, final_json_path)
            return result

        t_engine0 = time.perf_counter()

        if not is_ancient:
            logger.info(f"[{sample_id}] Direct routing to DissentEngine (Block 4: modern/zero-damage mtDNA certified).")
            try:
                estimate, diagnostics, tsv_file = DissentEngine.run(
                    sample_id=sample_id,
                    input_bam=effective_bam,
                    reference_fasta=context.reference_fasta_path,
                    output_dir=sample_out_dir,
                    policy=self.policy.dissent,
                    expected_contig=context.reference_contig,
                )
                audit.step_runtimes_seconds["dissent"] = round(time.perf_counter() - t_engine0, 2)
                outputs.dissent_tsv = tsv_file
            except Exception as d_exc:
                logger.error(f"[{sample_id}] DissentEngine execution failed: {d_exc}", exc_info=True)
                result.add_error(ContaminationErrorCode.DISSENT_EXECUTION_FAILED, f"DissentEngine failed: {d_exc}")
                result.audit = audit
                result.outputs = outputs
                result.mark_finished()
                result.runtime_seconds = round(time.perf_counter() - t_task_start, 2)
                self._save_contamination_json(result, final_json_path)
                return result
        else:
            try:
                schmutzi_out = SchmutziEngine.run(
                    sample_id=sample_id,
                    input_bam=effective_bam,
                    reference_fasta=context.reference_fasta_path,
                    database_dir=self.policy.schmutzi.contaminant_database,
                    output_dir=output_dir,
                    threads=self.policy.schmutzi.threads,
                    max_iterations=self.policy.schmutzi.max_iterations,
                    use_wrapped=self.policy.schmutzi.use_wrapped_reference,
                )
                audit.step_runtimes_seconds["schmutzi"] = round(time.perf_counter() - t_engine0, 2)
                outputs.cont_est_file = schmutzi_out.cont_est_file
                outputs.cont_deam_file = schmutzi_out.cont_deam_file
                outputs.endogenous_fasta = schmutzi_out.endogenous_fasta
                outputs.contaminant_fasta = schmutzi_out.contaminant_fasta

                estimate = schmutzi_out.estimate
                diagnostics = schmutzi_out.diagnostics

            except SchmutziExecutionError as exc:
                logger.error(f"[{sample_id}] Schmutzi Bayesian modeling execution error: {exc}")
                is_zero_damage_crash = any(k in str(exc).lower() for k in ("deam.config", "splitdeam", "-nan"))
                if self.policy.allow_fallback_to_dissent or is_zero_damage_crash:
                    logger.warning(f"[{sample_id}] Falling back to Track B Dissent Engine ({exc})...")
                    result.add_warning(ContaminationWarningCode.FALLBACK_TO_DISSENT, f"Schmutzi failed ({exc}); dissent used.")
                    try:
                        estimate, diagnostics, tsv_file = DissentEngine.run(
                            sample_id=sample_id,
                            input_bam=effective_bam,
                            reference_fasta=context.reference_fasta_path,
                            output_dir=sample_out_dir,
                            policy=self.policy.dissent,
                            expected_contig=context.reference_contig,
                        )
                        audit.step_runtimes_seconds["dissent"] = round(time.perf_counter() - t_engine0, 2)
                        audit.full_bayesian_estimation_performed = False
                        outputs.dissent_tsv = tsv_file
                        result.fallback_reason = str(exc)
                    except Exception as fb_exc:
                        result.add_error(ContaminationErrorCode.DISSENT_EXECUTION_FAILED, f"Fallback dissent failed: {fb_exc}")
                        result.audit = audit
                        result.outputs = outputs
                        result.mark_finished()
                        result.runtime_seconds = round(time.perf_counter() - t_task_start, 2)
                        self._save_contamination_json(result, final_json_path)
                        return result
                else:
                    result.add_error(ContaminationErrorCode.SCHMUTZI_EXECUTION_FAILED, str(exc))
                    result.audit = audit
                    result.outputs = outputs
                    result.mark_finished()
                    result.runtime_seconds = round(time.perf_counter() - t_task_start, 2)
                    self._save_contamination_json(result, final_json_path)
                    return result

        qc_result = ContaminationQCValidator.evaluate_sample(
            sample_id=sample_id,
            estimate=estimate,
            diagnostics=diagnostics,
            outputs=outputs,
            audit=audit,
            policy=self.policy,
            output_json_dir=None,
        )

        qc_result.mark_finished()
        qc_result.runtime_seconds = round(time.perf_counter() - t_task_start, 2)
        qc_result.outputs.contamination_json = final_json_path
        self._save_contamination_json(qc_result, final_json_path)

        logger.info(
            f"[{sample_id}] Block 5 contamination analysis COMPLETE in {qc_result.runtime_seconds}s "
            f"| Status: {qc_result.status.value} | Category: {qc_result.gating.contamination_category.value} "
            f"| Gating Action: {qc_result.gating.recommended_action.value}"
        )

        return qc_result