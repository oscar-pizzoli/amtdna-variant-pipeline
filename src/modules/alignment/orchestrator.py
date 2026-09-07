"""Block 3 Alignment & Filtering Orchestrator: connects reference validation, BWA alignment, structural preparation, deduplication, filtering, and QC auditing."""

from dataclasses import dataclass, field
import logging
from pathlib import Path
import traceback
from typing import Any, Dict, List, Optional

from src.modules.alignment.models import (
    AlignerType,
    AlignmentErrorCode,
    AlignmentFilterPolicy,
    AlignmentMetrics,
    AlignmentOutputs,
    AlignmentResult,
    AlignmentStatus,
    AlignmentToolMetadata,
    AlignmentWarningCode,
    CoverageMetrics,
    DedupMode,
    MultiMappingPolicy,
    ReadGroupInfo,
    ReferenceInfo,
)
from src.modules.alignment.reference_manager import (
    ReferenceManager,
    ReferenceManagerError,
    ReferenceValidationResult,
)
from src.modules.alignment.bwa_aligner import (
    BWAAligner,
    BWAAlignerError,
    AlignmentBranchResults,
)
from src.modules.alignment.samtools_pipeline import (
    SamtoolsPipeline,
    SamtoolsPipelineError,
    PipelinePreparationResult,
)
from src.modules.alignment.dedup_engine import (
    DedupEngine,
    DedupEngineError,
    DedupExecutionResult,
)
from src.modules.alignment.filter_engine import (
    AlignmentFilterEngine,
    FilterEngineError,
    FilterExecutionResult,
)
from src.modules.alignment.metrics_collector import (
    MetricsCollector,
    MetricsCollectorError,
    MetricsCollectionResult,
)
from src.modules.preprocessing.models import SampleRecord

logger = logging.getLogger("mtDNA_Pipeline.AlignmentOrchestrator")


class AlignmentOrchestratorError(Exception):
    """Raised when an unrecoverable orchestration failure occurs in Block 3."""
    pass


class AlignmentOrchestrator:
    """Orchestrates end-to-end alignment, deduplication, filtering, and metric collection for ancient and modern mtDNA samples."""

    @classmethod
    def evaluate_qc_thresholds(
        cls,
        metrics: AlignmentMetrics,
        policy: AlignmentFilterPolicy,
        result: AlignmentResult,
    ) -> None:
        """Evaluates final alignment and coverage metrics against configured QC thresholds and assigns diagnostic warnings.

        Args:
            metrics: Consolidated AlignmentMetrics container.
            policy: AlignmentFilterPolicy defining operational cutoffs.
            result: Target AlignmentResult to record warnings and status updates.
        """
        cov: CoverageMetrics = metrics.coverage

        # 1. Breadth of coverage threshold check (>= 3x coverage)
        min_breadth_pct = policy.min_breadth_3x_fraction * 100.0
        if cov.breadth_3x_pct < min_breadth_pct:
            msg = (
                f"Low breadth of coverage at >=3x: {cov.breadth_3x_pct:.2f}% "
                f"(configured threshold: {min_breadth_pct:.2f}%)."
            )
            logger.warning(f"[{result.sample_id}] {msg}")
            result.add_warning(msg, AlignmentWarningCode.LOW_BREADTH_COVERAGE)

        # 2. Mean depth of coverage check
        if cov.mean_depth_post_dedup < policy.min_mean_depth:
            msg = (
                f"Low mean depth of coverage: {cov.mean_depth_post_dedup:.2f}x "
                f"(configured threshold: {policy.min_mean_depth:.2f}x)."
            )
            logger.warning(f"[{result.sample_id}] {msg}")
            result.add_warning(msg, AlignmentWarningCode.LOW_DEPTH_COVERAGE)

        # 3. PCR Duplication rate check
        dup_rate = metrics.dedup_metrics.duplication_rate_reads
        if dup_rate > policy.max_duplication_rate_fraction:
            msg = (
                f"High PCR duplication rate: {dup_rate * 100.0:.2f}% "
                f"(configured threshold: {policy.max_duplication_rate_fraction * 100.0:.2f}%)."
            )
            logger.warning(f"[{result.sample_id}] {msg}")
            result.add_warning(msg, AlignmentWarningCode.HIGH_DUPLICATION_RATE)
        if metrics.raw_mapping_rate_pct is not None and metrics.raw_mapping_rate_pct < 10.0:
            msg = f"Extremely low raw mapping rate: {metrics.raw_mapping_rate_pct:.2f}% (configured threshold: 10.00%)."
            logger.warning(f"[{result.sample_id}] {msg}")
            result.add_warning(msg, AlignmentWarningCode.LOW_MAPPING_RATE)

        # 5. Check if zero reads survived filtering
        if metrics.retained_filtered_reads == 0:
            msg = "Zero reads survived quality and multi-mapping filtering."
            logger.error(f"[{result.sample_id}] {msg}")
            result.add_error(msg, AlignmentErrorCode.NO_MAPPED_MT_READS)

    @classmethod
    def run_sample_alignment(
        cls,
        sample: SampleRecord,
        reference_info: ReferenceInfo,
        output_dir: Path,
        log_dir: Path,
        policy: Optional[AlignmentFilterPolicy] = None,
        threads: int = 4,
        parallel_aln: bool = True,
        keep_temporary_sai: bool = False,
        keep_raw_merged_bam: bool = True,
        use_csi: bool = False,
        fast_mode: bool = False,
        scan_branch_bams: bool = False,
        seed_length: int = 1000,
        max_diff: float = 0.01,
        max_gap_opens: int = 2,
        max_gap_extensions: int = -1,
        trim_quality_threshold: int = 0,
        platform: str = "ILLUMINA",
        platform_unit: Optional[str] = None,
        sequencing_center: Optional[str] = None,
        cleanup_on_failure: bool = True,
    ) -> AlignmentResult:
        """Executes the complete Block 3 pipeline for a single sample.

        Workflow:
            1. Validate reference genome and indices.
            2. Multi-branch BWA-aln alignment with branch-specific @RG injection.
            3. Merge branches, name collate, fixmate mate sync, coordinate sort, and index.
            4. PCR duplicate flagging and optional physical deduplication (Dual-BAM).
            5. Bitwise flag filtering and MAPQ cutoff thresholding.
            6. Flagstat extraction, depth profiling, breadth calculation, and QC validation.

        Args:
            sample: Validated SampleRecord containing preprocessed FASTQ files and metrics.
            reference_info: Fully validated ReferenceInfo instance.
            output_dir: Directory where BAMs and reports will be saved.
            log_dir: Directory where stdout/stderr process logs will be saved.
            policy: Filtering, deduplication, and quality thresholds policy.
            threads: Number of CPU threads allocated for alignment and samtools pipelines.
            parallel_aln: If True, aligns paired-end R1 and R2 concurrently.
            keep_temporary_sai: If True, preserves intermediate .sai index files.
            keep_raw_merged_bam: If True, preserves pre-sorted merged BAM on disk.
            use_csi: If True, generates .csi indices instead of standard .bai.
            fast_mode: If True, skips granular Python read inspection for maximum speed.
            scan_branch_bams: If True, inspects branch BAMs to count mapped reads per branch.
            seed_length: BWA aln seed length (-l).
            max_diff: BWA aln maximum edit distance / error rate (-n).
            max_gap_opens: BWA aln maximum gap opens (-o).
            max_gap_extensions: BWA aln maximum gap extensions (-e).
            trim_quality_threshold: BWA aln quality trimming threshold (-q).
            platform: Sequencing platform tag for @RG PL.
            platform_unit: Flowcell/run unit tag for @RG PU.
            sequencing_center: Sequencing center tag for @RG CN.
            cleanup_on_failure: If True, deletes partial files upon execution failure.

        Returns:
            Fully populated AlignmentResult object representing the outcome of Block 3.
        """
        sample_id = sample.sample_id
        active_policy = policy or AlignmentFilterPolicy()

        sample_out_dir = output_dir / sample_id
        sample_log_dir = log_dir / sample_id
        sample_out_dir.mkdir(parents=True, exist_ok=True)
        sample_log_dir.mkdir(parents=True, exist_ok=True)

        tool_meta = AlignmentToolMetadata(
            aligner_name=AlignerType.BWA_ALN,
            reference_file_sha256=reference_info.file_sha256,
            working_directory=str(sample_out_dir.resolve()),
            parameters={
                "threads": threads,
                "seed_length": seed_length,
                "max_diff": max_diff,
                "max_gap_opens": max_gap_opens,
                "max_gap_extensions": max_gap_extensions,
                "trim_quality_threshold": trim_quality_threshold,
                "dedup_mode": active_policy.dedup_mode.value,
                "min_mapq": active_policy.multimapping.min_mapq,
                "filter_unmapped": active_policy.multimapping.filter_unmapped,
                "filter_secondary": active_policy.multimapping.filter_secondary,
                "filter_supplementary": active_policy.multimapping.filter_supplementary,
            },
        )
        tool_meta.mark_started()

        result = AlignmentResult(
            sample_id=sample_id,
            status=AlignmentStatus.UNTESTED,
            policy=active_policy,
            reference=reference_info,
            tool_meta=tool_meta,
        )

        logger.info(f"[{sample_id}] ========================================================")
        logger.info(f"[{sample_id}] Starting Block 3: Alignment, Filtering, and QC Pipeline")
        logger.info(f"[{sample_id}] ========================================================")

        try:
            # ------------------------------------------------------------------
            # Step 1: Reference Genome Precondition Check
            # ------------------------------------------------------------------
            if not reference_info.has_bwa_index or not reference_info.has_fai_index:
                raise AlignmentOrchestratorError(
                    f"Reference genome indices are missing for {reference_info.fasta_path.name}."
                )

            # ------------------------------------------------------------------
            # Step 2: Multi-Branch BWA Alignment
            # ------------------------------------------------------------------
            logger.info(f"[{sample_id}] [Stage 1/5] Running multi-branch BWA alignment...")
            branch_results: AlignmentBranchResults = BWAAligner.align_sample_branches(
                sample=sample,
                reference_info=reference_info,
                output_dir=sample_out_dir / "raw_branches",
                log_dir=sample_log_dir / "bwa",
                threads=threads,
                seed_length=seed_length,
                max_diff=max_diff,
                max_gap_opens=max_gap_opens,
                max_gap_extensions=max_gap_extensions,
                trim_quality_threshold=trim_quality_threshold,
                platform=platform,
                platform_unit=platform_unit,
                sequencing_center=sequencing_center,
                keep_temporary_sai=keep_temporary_sai,
                parallel_aln=parallel_aln,
            )

            result.read_groups = branch_results.all_read_groups
            for cmd in branch_results.tool_meta.commands_executed:
                tool_meta.add_command(cmd)

            # ------------------------------------------------------------------
            # Step 3: Structural BAM Preparation (Merge -> Collate -> Fixmate -> Sort)
            # ------------------------------------------------------------------
            logger.info(f"[{sample_id}] [Stage 2/5] Merging branches and streaming collate -> fixmate -m -> sort...")
            prep_result: PipelinePreparationResult = SamtoolsPipeline.process_sample_branches(
                branch_results=branch_results,
                output_dir=sample_out_dir,
                log_dir=sample_log_dir / "samtools",
                threads=threads,
                keep_raw_merged_bam=keep_raw_merged_bam,
                cleanup_on_failure=cleanup_on_failure,
                expected_contig=reference_info.contig_name,
            )

            if not prep_result.is_successful or not prep_result.sorted_bam:
                raise AlignmentOrchestratorError(f"Structural BAM preparation failed: {prep_result.error_message}")

            for cmd in prep_result.command_lines:
                tool_meta.add_command(cmd)

            # ------------------------------------------------------------------
            # Step 4: PCR Deduplication (samtools markdup -s -> Dual-BAM)
            # ------------------------------------------------------------------
            logger.info(f"[{sample_id}] [Stage 3/5] Executing PCR duplicate marking and deduplication...")
            dedup_result: DedupExecutionResult = DedupEngine.execute_dedup(
                sample_id=sample_id,
                input_sorted_bam=prep_result.sorted_bam,
                output_dir=sample_out_dir,
                log_dir=sample_log_dir / "dedup",
                dedup_mode=active_policy.dedup_mode,
                threads=threads,
                use_csi=use_csi,
                cleanup_on_failure=cleanup_on_failure,
            )

            if not dedup_result.is_successful:
                raise AlignmentOrchestratorError(f"PCR deduplication failed: {dedup_result.error_message}")

            for cmd in dedup_result.command_lines:
                tool_meta.add_command(cmd)

            # ------------------------------------------------------------------
            # Step 5: Quality & Multi-Mapping Filtering
            # ------------------------------------------------------------------
            logger.info(f"[{sample_id}] [Stage 4/5] Applying quality, flag, and multi-mapping filters...")
            filter_result: FilterExecutionResult = AlignmentFilterEngine.execute_filtering(
                sample_id=sample_id,
                input_markdup_bam=dedup_result.markdup_bam,
                input_dedup_bam=dedup_result.dedup_bam,
                output_dir=sample_out_dir,
                log_dir=sample_log_dir / "filter",
                policy=active_policy,
                threads=threads,
                fast_mode=fast_mode,
                use_csi=use_csi,
                cleanup_on_failure=cleanup_on_failure,
            )

            if not filter_result.is_successful or not filter_result.primary_downstream_bam:
                raise AlignmentOrchestratorError(f"Alignment filtering failed: {filter_result.error_message}")

            for cmd in filter_result.command_lines:
                tool_meta.add_command(cmd)

            # ------------------------------------------------------------------
            # Step 6: Quality Metrics Collection & Whole-Genome Depth Profiling
            # ------------------------------------------------------------------
            logger.info(f"[{sample_id}] [Stage 5/5] Collecting flagstats, depth distribution, and coverage metrics...")
            metrics_result: MetricsCollectionResult = MetricsCollector.collect_all_metrics(
                sample_id=sample_id,
                reference_info=reference_info,
                branch_results=branch_results,
                dedup_result=dedup_result,
                filter_result=filter_result,
                output_dir=sample_out_dir / "reports",
                log_dir=sample_log_dir / "metrics",
                prepared_bam=prep_result.sorted_bam,
                threads=max(1, threads // 2),
                scan_branch_bams=scan_branch_bams,
                min_base_quality=0,
                min_mapping_quality=0,
            )

            if not metrics_result.is_successful:
                raise AlignmentOrchestratorError(f"Metrics collection failed: {metrics_result.error_message}")

            for cmd in metrics_result.command_lines:
                tool_meta.add_command(cmd)

            # ------------------------------------------------------------------
            # Step 7: Map Physical Outputs & Consolidated Metrics
            # ------------------------------------------------------------------
            result.metrics = metrics_result.alignment_metrics
            result.outputs = AlignmentOutputs(
                merged_raw_bam=prep_result.merged_raw_bam,
                sorted_bam=prep_result.sorted_bam,
                sorted_bai=prep_result.sorted_bai,
                markdup_bam=dedup_result.markdup_bam,
                markdup_bai=dedup_result.markdup_bai,
                dedup_bam=dedup_result.dedup_bam,
                dedup_bai=dedup_result.dedup_bai,
                filtered_bam=filter_result.primary_downstream_bam,
                filtered_bai=filter_result.filtered_dedup_bai or filter_result.filtered_markdup_bai,
                markdup_metrics_txt=dedup_result.markdup_metrics_txt,
                flagstat_pre_filter_txt=sample_out_dir / "reports" / f"{sample_id}.flagstat.prepared.txt",
                flagstat_post_dedup_txt=sample_out_dir / "reports" / f"{sample_id}.flagstat.dedup.txt" if dedup_result.dedup_bam else None,
                flagstat_post_filter_txt=sample_out_dir / "reports" / f"{sample_id}.flagstat.filtered.txt",
                samtools_stats_txt=prep_result.baseline_stats_txt,
                depth_tsv=metrics_result.depth_tsv,
            )

            # ------------------------------------------------------------------
            # Step 8: QC Threshold Evaluation & Final Status Determination
            # ------------------------------------------------------------------
            cls.evaluate_qc_thresholds(
                metrics=result.metrics,
                policy=active_policy,
                result=result,
            )

            if result.status == AlignmentStatus.UNTESTED:
                result.status = AlignmentStatus.PASS

            logger.info(
                f"[{sample_id}] Block 3 alignment completed successfully with status: {result.status.value}"
            )

        except Exception as exc:
            err_msg = f"Alignment orchestration failed for sample {sample_id}: {exc}"
            logger.error(err_msg)
            logger.debug(traceback.format_exc())
            result.add_error(err_msg, AlignmentErrorCode.ALIGNMENT_EXECUTION_FAILED)

        finally:
            tool_meta.mark_finished()
            result.tool_meta = tool_meta

        return result
