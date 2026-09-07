"""Block 4 Orchestrator: Ancient DNA Damage Quantification, Rescaling, and Authentication Pipeline."""

from dataclasses import dataclass, field
import logging
from pathlib import Path
import shutil
import subprocess
import time
from typing import Any, Dict, List, Optional, Tuple, Union
import pysam

from src.modules.alignment.models import ReferenceInfo
from src.modules.damage.damage_quantifier import (
    DamageQuantifierEngine,
    DamageQuantifierError,
)
from src.modules.damage.models import (
    AuthenticationVerdict,
    DamageErrorCode,
    DamageInputBAMType,
    DamageOutputFiles,
    DamagePolicy,
    DamageResult,
    DamageSignalLevel,
    DamageStatus,
    DamageWarningCode,
    LibraryTreatment,
)
from src.modules.damage.pmd_engine import (
    PMDExecutionMetadata,
    PMDScoringEngine,
    PMDScoringError,
)
from src.modules.damage.qc_evaluator import DamageQCValidator
from src.modules.damage.rescale_engine import (
    DamageRescaleError,
    QualityRescaleEngine,
    RescaleExecutionMetadata,
)

logger = logging.getLogger("mtDNA_Pipeline.DamageOrchestrator")


class DamageOrchestrationError(Exception):
    """Raised when an unrecoverable failure occurs during Block 4 orchestration."""
    pass


class DamageBAMValidationError(Exception):
    """Raised when pre-flight checks on input BAM fail."""
    pass


@dataclass
class DamageAuditTrail:
    """Comprehensive execution metadata consolidating tool versions, commands, and sub-engine metrics."""

    sample_id: str
    input_bam_resolved: Optional[Path] = None
    input_bam_type: Optional[str] = None
    input_mapped_reads: int = 0
    rescaling_requested: bool = False
    rescaling_applied: bool = False
    rescaling_failed: bool = False
    pmd_requested: bool = False
    pmd_completed: bool = False
    pmd_failed: bool = False
    rescale_metadata: Optional[RescaleExecutionMetadata] = None
    pmd_metadata: Optional[PMDExecutionMetadata] = None
    step_runtimes_seconds: Dict[str, float] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Serializes orchestration audit trail to a dictionary."""
        return {
            "sample_id": self.sample_id,
            "input_bam_resolved": str(self.input_bam_resolved) if self.input_bam_resolved else None,
            "input_bam_type": self.input_bam_type,
            "input_mapped_reads": self.input_mapped_reads,
            "rescaling_requested": self.rescaling_requested,
            "rescaling_applied": self.rescaling_applied,
            "rescaling_failed": self.rescaling_failed,
            "pmd_requested": self.pmd_requested,
            "pmd_completed": self.pmd_completed,
            "pmd_failed": self.pmd_failed,
            "rescale_metadata": self.rescale_metadata.to_dict() if self.rescale_metadata else None,
            "pmd_metadata": self.pmd_metadata.to_dict() if self.pmd_metadata else None,
            "step_runtimes_seconds": {k: round(v, 2) for k, v in self.step_runtimes_seconds.items()},
        }


class DamageOrchestrator:
    """End-to-end coordinator for Block 4 damage profiling, quality rescaling, and authenticity testing."""

    @classmethod
    def validate_input_bam_preflight(
        cls,
        bam_path: Path,
        expected_contig: str,
    ) -> int:
        """Executes full pre-flight integrity, header compatibility, and mapped read count checks on input BAM.

        Args:
            bam_path: Path to input BAM file.
            expected_contig: Target mitochondrial reference contig name.

        Returns:
            Count of mapped reads in input BAM.

        Raises:
            DamageBAMValidationError: If BAM is corrupt, truncated, lacks index, or has 0 mapped reads.
        """
        if not bam_path.is_file() or bam_path.stat().st_size == 0:
            raise DamageBAMValidationError(f"Input BAM does not exist or is empty: {bam_path}")

        # 1. samtools quickcheck validation
        if shutil.which("samtools") is not None:
            try:
                res = subprocess.run(["samtools", "quickcheck", str(bam_path)], capture_output=True, text=True, check=False)
                if res.returncode != 0:
                    raise DamageBAMValidationError(f"BAM failed 'samtools quickcheck': {bam_path.name}")
            except OSError as exc:
                logger.warning(f"Could not execute samtools quickcheck: {exc}")

        # 2. Check or generate BAI index
        bai_path = Path(f"{bam_path}.bai")
        if not bai_path.is_file() or bai_path.stat().st_size == 0:
            logger.info(f"Generating missing index for input BAM: {bam_path.name}...")
            try:
                pysam.index(str(bam_path))
            except Exception as exc:
                raise DamageBAMValidationError(f"Failed indexing input BAM '{bam_path.name}': {exc}") from exc

        # 3. Inspect header and mapped reads
        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                if expected_contig not in bam.references:
                    raise DamageBAMValidationError(
                        f"Expected contig '{expected_contig}' missing from BAM references: {bam.references}"
                    )

                mapped_reads = bam.count(contig=expected_contig, read_callback="all")
                if mapped_reads == 0:
                    raise DamageBAMValidationError(
                        f"Input BAM '{bam_path.name}' contains 0 mapped reads on contig '{expected_contig}'."
                    )
                return mapped_reads
        except Exception as exc:
            if isinstance(exc, DamageBAMValidationError):
                raise
            raise DamageBAMValidationError(f"Failed inspecting records in BAM '{bam_path.name}': {exc}") from exc

    @classmethod
    def resolve_input_bam(
        cls,
        sample_id: str,
        policy: DamagePolicy,
        dedup_bam: Optional[Path] = None,
        markdup_bam: Optional[Path] = None,
        custom_bam: Optional[Path] = None,
        allow_fallback: bool = False,
    ) -> Path:
        """Strictly resolves target input BAM path according to configured DamageInputBAMType policy.

        Args:
            sample_id: Biological sample identifier.
            policy: Validated DamagePolicy container.
            dedup_bam: Deduplicated BAM path from Block 3.
            markdup_bam: Duplicate-flagged BAM path from Block 3.
            custom_bam: Explicit BAM path override.
            allow_fallback: If True, allows falling back to alternate BAM if requested type is absent.

        Returns:
            Resolved Path to target BAM file.

        Raises:
            DamageOrchestrationError: If requested BAM is missing and fallback is disabled.
        """
        if custom_bam and custom_bam.is_file():
            logger.info(f"[{sample_id}] Using explicit custom input BAM: {custom_bam.name}")
            return custom_bam

        target_type = policy.input_bam_type

        if target_type == DamageInputBAMType.FILTERED_DEDUP:
            if dedup_bam and dedup_bam.is_file():
                return dedup_bam
            if allow_fallback and markdup_bam and markdup_bam.is_file():
                logger.warning(
                    f"[{sample_id}] Requested FILTERED_DEDUP BAM not found. Falling back to FILTERED_MARKDUP: {markdup_bam.name}"
                )
                return markdup_bam
            raise DamageOrchestrationError(
                f"Policy requested FILTERED_DEDUP, but dedup BAM is missing: {dedup_bam} (fallback=False)"
            )

        if target_type == DamageInputBAMType.FILTERED_MARKDUP:
            if markdup_bam and markdup_bam.is_file():
                return markdup_bam
            if allow_fallback and dedup_bam and dedup_bam.is_file():
                logger.warning(
                    f"[{sample_id}] Requested FILTERED_MARKDUP BAM not found. Falling back to FILTERED_DEDUP: {dedup_bam.name}"
                )
                return dedup_bam
            raise DamageOrchestrationError(
                f"Policy requested FILTERED_MARKDUP, but markdup BAM is missing: {markdup_bam} (fallback=False)"
            )

        if target_type == DamageInputBAMType.RAW_SORTED:
            if custom_bam and custom_bam.is_file():
                return custom_bam
            raise DamageOrchestrationError("Policy requested RAW_SORTED, but no valid sorted BAM was provided.")

        raise DamageOrchestrationError(f"Unsupported or unrecognized input BAM type: {target_type}")

    @classmethod
    def run_sample_damage_pipeline(
        cls,
        sample_id: str,
        reference_info: ReferenceInfo,
        output_dir: Path,
        log_dir: Path,
        policy: Optional[DamagePolicy] = None,
        dedup_bam: Optional[Path] = None,
        markdup_bam: Optional[Path] = None,
        custom_input_bam: Optional[Path] = None,
        threads: int = 2,
        allow_input_bam_fallback: bool = False,
        require_pmd_success: bool = False,
    ) -> DamageResult:
        """Executes full Block 4 workflow: mapDamage estimation -> Rescaling -> PMD scoring -> QC Evaluation.

        Args:
            sample_id: Biological sample identifier.
            reference_info: Target ReferenceInfo metadata container from Block 1/3.
            output_dir: Root output directory for damage deliverables.
            log_dir: Root directory for execution logs.
            policy: Configured DamagePolicy (defaults used if None).
            dedup_bam: Path to upstream deduplicated BAM from Block 3.
            markdup_bam: Path to upstream duplicate-flagged BAM from Block 3.
            custom_input_bam: Explicit BAM override if applicable.
            threads: CPU threads allocated for samtools operations.
            allow_input_bam_fallback: Whether to allow fallback to alternate BAM if requested type is absent.
            require_pmd_success: If True, failure in PMD scoring causes block execution failure.

        Returns:
            Fully populated and validated DamageResult container.
        """
        active_policy = policy or DamagePolicy()
        result = DamageResult(
            sample_id=sample_id,
            policy=active_policy,
        )
        result.mark_started()

        audit = DamageAuditTrail(
            sample_id=sample_id,
            input_bam_type=active_policy.input_bam_type.value,
            rescaling_requested=active_policy.enable_rescaling,
            pmd_requested=active_policy.enable_pmd_scoring,
        )

        logger.info("=" * 65)
        logger.info(f"[{sample_id}] Starting Block 4: aDNA Damage & Authentication Pipeline")
        logger.info(f"[{sample_id}] Target Contig: {reference_info.contig_name} ({reference_info.length_bp} bp)")
        logger.info("=" * 65)

        out_damage_dir = output_dir / sample_id
        out_damage_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)

        try:
            # 1. Input BAM Resolution
            input_bam = cls.resolve_input_bam(
                sample_id=sample_id,
                policy=active_policy,
                dedup_bam=dedup_bam,
                markdup_bam=markdup_bam,
                custom_bam=custom_input_bam,
                allow_fallback=allow_input_bam_fallback,
            )
            audit.input_bam_resolved = input_bam

            # 2. Pre-flight Validation
            mapped_count = cls.validate_input_bam_preflight(
                bam_path=input_bam,
                expected_contig=reference_info.contig_name,
            )
            audit.input_mapped_reads = mapped_count
            logger.info(
                f"[{sample_id}] Input BAM verified: {input_bam.name} ({mapped_count:,} mapped reads on {reference_info.contig_name})"
            )

            # 3. Track A: mapDamage2 Statistical Quantification (Terminal Frequencies & Bayesian Fit)
            t_start = time.perf_counter()
            logger.info(f"[{sample_id}] [1/4] Running mapDamage2 statistical quantification...")
            raw_outputs, terminal_freqs, length_metrics, bayesian_params = DamageQuantifierEngine.run_mapdamage_stats(
                sample_id=sample_id,
                input_bam=input_bam,
                reference_fasta=reference_info.fasta_path,
                output_dir=output_dir,
                log_dir=log_dir,
                policy=active_policy,
                threads=threads,
            )
            audit.step_runtimes_seconds["mapdamage_quantification"] = time.perf_counter() - t_start

            result.outputs.misincorporation_txt = raw_outputs.misincorporation_txt
            result.outputs.length_distribution_txt = raw_outputs.length_distribution_txt
            result.outputs.damage_plot_pdf = raw_outputs.damage_plot_pdf

            # 4. Track A: Base Quality Score Rescaling (--rescale-only)
            rescaling_applied = False
            if active_policy.enable_rescaling:
                t_start = time.perf_counter()
                logger.info(f"[{sample_id}] [2/4] Running Base Quality Rescaling (--rescale-only)...")
                try:
                    rescaled_bam, rescaled_bai, rescale_meta = QualityRescaleEngine.run_rescaling(
                        sample_id=sample_id,
                        input_bam=input_bam,
                        reference_fasta=reference_info.fasta_path,
                        damage_model_dir=out_damage_dir,
                        output_dir=output_dir,
                        log_dir=log_dir,
                        policy=active_policy,
                        expected_contig=reference_info.contig_name,
                    )
                    audit.rescale_metadata = rescale_meta
                    audit.rescaling_applied = (rescaled_bam is not None)
                    result.outputs.rescaled_bam = rescaled_bam
                    result.outputs.rescaled_bai = rescaled_bai
                    rescaling_applied = audit.rescaling_applied
                except DamageRescaleError as r_err:
                    audit.rescaling_failed = True
                    result.add_error(DamageErrorCode.RESCALING_FAILED, str(r_err))
                    logger.error(f"[{sample_id}] Quality Rescaling failed: {r_err}")
                finally:
                    audit.step_runtimes_seconds["quality_rescaling"] = time.perf_counter() - t_start
            else:
                logger.info(f"[{sample_id}] [2/4] Quality Rescaling skipped (enable_rescaling=False).")

            # 5. Track B: Read-Level PMD Scoring & Filtering (pmdtools on original input BAM)
            pmd_metrics = None
            if active_policy.enable_pmd_scoring:
                t_start = time.perf_counter()
                logger.info(f"[{sample_id}] [3/4] Running per-read PMD likelihood scoring (pmdtools)...")
                try:
                    scores_tsv, pmd_bam, pmd_bai, pmd_metrics, pmd_meta = PMDScoringEngine.execute_pmd_pipeline(
                        sample_id=sample_id,
                        input_bam=input_bam,
                        reference_fasta=reference_info.fasta_path,
                        output_dir=output_dir,
                        log_dir=log_dir,
                        policy=active_policy,
                        threads=threads,
                    )
                    audit.pmd_metadata = pmd_meta
                    audit.pmd_completed = True
                    result.outputs.pmd_scores_tsv = scores_tsv
                    result.outputs.pmd_filtered_bam = pmd_bam
                    result.outputs.pmd_filtered_bai = pmd_bai
                except PMDScoringError as p_err:
                    audit.pmd_failed = True
                    if require_pmd_success:
                        result.add_error(DamageErrorCode.PMDTOOLS_EXECUTION_FAILED, str(p_err))
                    else:
                        result.add_warning(
                            DamageWarningCode.LOW_READ_COUNT_FOR_DAMAGE_MODEL,
                            f"PMD scoring encountered an issue: {p_err}",
                        )
                    logger.warning(f"[{sample_id}] PMD Scoring warning: {p_err}")
                finally:
                    audit.step_runtimes_seconds["pmd_scoring"] = time.perf_counter() - t_start
            else:
                logger.info(f"[{sample_id}] [3/4] PMD scoring skipped (enable_pmd_scoring=False).")

            # 6. Damage QC Evaluation & JSON Deliverable Export
            t_start = time.perf_counter()
            logger.info(f"[{sample_id}] [4/4] Evaluating damage evidence and assigning authenticity verdict...")
            qc_result = DamageQCValidator.evaluate_sample(
                sample_id=sample_id,
                input_bam_path=input_bam,
                reference_contig=reference_info.contig_name,
                terminal_damage=terminal_freqs,
                fragment_length=length_metrics,
                bayesian_params=bayesian_params,
                pmd_metrics=pmd_metrics,
                outputs=result.outputs,
                policy=active_policy,
                rescaling_applied=rescaling_applied,
                output_json_dir=out_damage_dir,
            )
            audit.step_runtimes_seconds["qc_evaluation"] = time.perf_counter() - t_start

            # Harmonize top-level container status and metadata
            result.status = qc_result.status
            result.profile = qc_result.profile
            result.outputs = qc_result.outputs
            result.warning_codes.extend(qc_result.warning_codes)
            result.warnings.extend(qc_result.warnings)
            result.error_codes.extend(qc_result.error_codes)
            result.errors.extend(qc_result.errors)

        except DamageBAMValidationError as b_err:
            result.add_error(DamageErrorCode.INPUT_BAM_NOT_FOUND, str(b_err))
            logger.error(f"[{sample_id}] Input BAM pre-flight validation failed: {b_err}")
        except DamageOrchestrationError as o_err:
            result.add_error(DamageErrorCode.INPUT_BAM_NOT_FOUND, str(o_err))
            logger.error(f"[{sample_id}] BAM resolution error: {o_err}")
        except DamageQuantifierError as dq_err:
            result.add_error(DamageErrorCode.MAPDAMAGE_EXECUTION_FAILED, str(dq_err))
            logger.error(f"[{sample_id}] Damage quantification aborted: {dq_err}")
        except Exception as unexp_exc:
            result.add_error(DamageErrorCode.CORRUPTED_DAMAGE_OUTPUT, f"Unexpected pipeline failure: {unexp_exc}")
            logger.exception(f"[{sample_id}] Unexpected exception during Block 4 execution: {unexp_exc}")
        finally:
            result.mark_finished()

        # Safe summary logging across success and early-abortion states
        logger.info("=" * 65)
        logger.info(f"[{sample_id}] Block 4 Completed with Status: {result.status.value}")
        
        prof = result.profile
        if prof and prof.damage_signal_level:
            logger.info(f"[{sample_id}] Damage Signal Level : {prof.damage_signal_level.value}")
        if prof and prof.authentication_verdict:
            logger.info(f"[{sample_id}] Authenticity Verdict : {prof.authentication_verdict.value}")
        if prof and prof.terminal_damage and prof.terminal_damage.n_c_ref_5p_pos1_total > 0:
            td = prof.terminal_damage
            logger.info(
                f"[{sample_id}] 5' C->T (pos1)      : {td.c_to_t_5p_pos1:.2%} "
                f"({td.n_c_to_t_5p_pos1_obs}/{td.n_c_ref_5p_pos1_total} obs)"
            )
        if prof and prof.fragment_length and prof.fragment_length.total_fragments_sampled > 0:
            logger.info(f"[{sample_id}] Median Length       : {prof.fragment_length.median_length:.1f} bp")
        if result.outputs and result.outputs.rescaled_bam:
            logger.info(f"[{sample_id}] Rescaled BAM       : {result.outputs.rescaled_bam.name}")
        
        logger.info(f"[{sample_id}] Total Runtime       : {result.runtime_seconds:.2f}s")
        logger.info("=" * 65)

        return result
