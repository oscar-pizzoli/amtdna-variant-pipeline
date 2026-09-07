"""Batch Runner for Block 5 (mtDNA Contamination Estimation & Two-Phase Fast-Gating).

Implements a two-phase cohort execution model:
  - Phase 1: High-concurrency Fast Screening pool (1 thread/sample, instant triage)
  - Phase 2: Bounded-concurrency Bayesian/Allelic modeling pool (Schmutzi for aDNA, DissentEngine for modern)
  - Strict host CPU budget enforcement (default: 50% CPU allocation ceiling)
  - Work conservation: 0 compute slots occupied by FAST_CLEAR samples
  - Multi-factor audit summaries tracking stage-specific runtimes, bypass rates, and gating
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from datetime import datetime, timezone
import inspect
import json
import logging
import math
import os
from pathlib import Path
import re
import sys
import time
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import pysam
import yaml

# Risoluzione dinamica della radice del progetto
PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

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
    EstimationConfidence,
    FastScreeningPolicy,
    FastScreenResult,
    GatingRecommendation,
    ModelQualityDiagnostics,
    SchmutziPolicy,
    ScreeningDecision,
    SubsamplingPolicy,
)
from src.modules.contamination.orchestrator import (
    ContaminationOrchestrator,
    ContaminationOrchestrationError,
)
from src.modules.contamination.qc_evaluator import ContaminationQCValidator
from src.modules.contamination.schmutzi_engine import (
    SchmutziEngine,
    SchmutziExecutionError,
    SchmutziExecutionMetadata,
    SchmutziOutput,
)

logger = logging.getLogger("mtDNA_Pipeline.BatchRunner")


@dataclass
class ResourceAllocation:
    """Quantitative allocation of system threads and multi-phase worker limits."""

    total_host_cores: int
    cpu_budget_fraction: float
    total_allocated_threads: int
    threads_per_schmutzi_instance: int
    max_fast_screen_workers: int
    max_schmutzi_workers: int
    oversubscription_allowed: bool = False

    def to_dict(self) -> Dict[str, Any]:
        return {
            "total_host_cores": self.total_host_cores,
            "cpu_budget_fraction": self.cpu_budget_fraction,
            "total_allocated_threads": self.total_allocated_threads,
            "threads_per_schmutzi_instance": self.threads_per_schmutzi_instance,
            "max_fast_screen_workers": self.max_fast_screen_workers,
            "max_schmutzi_workers": self.max_schmutzi_workers,
            "oversubscription_allowed": self.oversubscription_allowed,
        }


@dataclass
class BatchSampleRecord:
    """Summary record for a single processed sample in the cohort."""

    sample_id: str
    status: str
    screening_decision: str
    bypass_applied: bool
    schmutzi_executed: bool
    contamination_category: str
    recommended_action: str
    point_estimate: float
    ci_lower_95: Optional[float]
    ci_upper_95: Optional[float]
    method_used: str
    diagnostic_depth: float
    genome_depth: float
    dissenting_sites: int
    fast_screen_runtime_seconds: float
    schmutzi_runtime_seconds: float
    total_runtime_seconds: float
    subsampling_applied: bool
    error_stage: str = ""
    error_summary: str = ""

    def to_tsv_row(self) -> str:
        ci_lower = f"{self.ci_lower_95:.4f}" if self.ci_lower_95 is not None else "NA"
        ci_upper = f"{self.ci_upper_95:.4f}" if self.ci_upper_95 is not None else "NA"
        return (
            f"{self.sample_id}\t{self.status}\t{self.screening_decision}\t"
            f"{self.bypass_applied}\t{self.schmutzi_executed}\t{self.contamination_category}\t"
            f"{self.recommended_action}\t{self.point_estimate:.4f}\t{ci_lower}\t{ci_upper}\t"
            f"{self.method_used}\t{self.diagnostic_depth:.1f}\t{self.genome_depth:.1f}\t"
            f"{self.dissenting_sites}\t{self.fast_screen_runtime_seconds:.2f}\t"
            f"{self.schmutzi_runtime_seconds:.2f}\t{self.total_runtime_seconds:.2f}\t"
            f"{self.subsampling_applied}\t{self.error_stage}\t{self.error_summary}\n"
        )


class ContaminationBatchRunner:
    """Coordinates two-phase cohort contamination estimation and strict resource management."""

    @classmethod
    def load_yaml_config(cls, config_path: Optional[Path] = None) -> Dict[str, Any]:
        """Loads configuration from YAML file."""
        target_path = config_path or (PROJECT_ROOT / "config" / "default_config.yaml")
        if not target_path.is_file():
            logger.warning(f"Config file not found at {target_path}, using built-in defaults.")
            return {}
        try:
            with open(target_path, "r", encoding="utf-8") as f:
                return yaml.safe_load(f) or {}
        except Exception as exc:
            logger.warning(f"Failed loading YAML config from {target_path}: {exc}")
            return {}

    @classmethod
    def resolve_path(cls, path_str: Optional[str], default_rel: Path) -> Path:
        """Resolves path string relative to PROJECT_ROOT."""
        if not path_str:
            return default_rel if default_rel.is_absolute() else (PROJECT_ROOT / default_rel)
        p = Path(path_str)
        return p if p.is_absolute() else (PROJECT_ROOT / p)

    @classmethod
    def detect_available_cores(cls, host_cores_override: Optional[int] = None) -> int:
        """Determines host CPU cores respecting affinity, containers, or explicit overrides."""
        if host_cores_override is not None and host_cores_override > 0:
            return host_cores_override

        try:
            return len(os.sched_getaffinity(0))
        except (AttributeError, OSError):
            return os.cpu_count() or 4

    @classmethod
    def calculate_resource_allocation(
        cls,
        threads_per_sample: int = 4,
        cpu_budget_fraction: float = 0.50,
        explicit_max_workers: Optional[int] = None,
        host_cores_override: Optional[int] = None,
        allow_oversubscription: bool = False,
    ) -> ResourceAllocation:
        """Calculates safe concurrent worker limits strictly bounding host core allocation."""
        host_cores = cls.detect_available_cores(host_cores_override)
        budget_pct = max(0.05, min(1.0, cpu_budget_fraction))
        allocated_threads = max(1, math.floor(host_cores * budget_pct))

        budget_schmutzi_workers = max(1, math.floor(allocated_threads / max(1, threads_per_sample)))

        if explicit_max_workers is not None and explicit_max_workers > 0:
            if allow_oversubscription:
                schmutzi_workers = explicit_max_workers
                logger.warning(f"Oversubscription allowed: Schmutzi worker count forced to {explicit_max_workers}.")
            else:
                schmutzi_workers = min(budget_schmutzi_workers, explicit_max_workers)
        else:
            schmutzi_workers = budget_schmutzi_workers

        fast_screen_workers = max(1, min(allocated_threads, host_cores))

        alloc = ResourceAllocation(
            total_host_cores=host_cores,
            cpu_budget_fraction=budget_pct,
            total_allocated_threads=allocated_threads,
            threads_per_schmutzi_instance=threads_per_sample,
            max_fast_screen_workers=fast_screen_workers,
            max_schmutzi_workers=schmutzi_workers,
            oversubscription_allowed=allow_oversubscription,
        )

        logger.info(
            f"Resource allocation initialized: {alloc.total_host_cores} host cores -> "
            f"{alloc.total_allocated_threads} threads budget ({alloc.cpu_budget_fraction:.0%}) | "
            f"FastScreen concurrency: {alloc.max_fast_screen_workers} workers | "
            f"Schmutzi concurrency: {alloc.max_schmutzi_workers} workers ({alloc.threads_per_schmutzi_instance} threads/instance)"
        )
        return alloc

    @classmethod
    def discover_bam_inputs(
        cls,
        input_dir: Path,
        pattern: str = "*.bam",
        sample_manifest: Optional[Path] = None,
    ) -> List[Tuple[str, Path]]:
        """Discovers input BAM files via structured manifest or defensive directory scanning."""
        discovered: List[Tuple[str, Path]] = []
        seen_sample_ids: Set[str] = set()

        if sample_manifest and sample_manifest.is_file():
            logger.info(f"Loading samples from manifest: {sample_manifest}")
            with open(sample_manifest, "r", encoding="utf-8") as f_m:
                for line in f_m:
                    line_str = line.strip()
                    if not line_str or line_str.startswith("#"):
                        continue
                    parts = line_str.split("\t") if "\t" in line_str else line_str.split(",")
                    if len(parts) >= 2:
                        s_id = parts[0].strip()
                        b_path = Path(parts[1].strip())
                        b_resolved = b_path if b_path.is_absolute() else (PROJECT_ROOT / b_path)
                        if b_resolved.is_file():
                            discovered.append((s_id, b_resolved))
                            seen_sample_ids.add(s_id)
            return discovered

        if not input_dir.is_dir():
            raise FileNotFoundError(f"Input directory not found: {input_dir}")

        all_bams = sorted(list(input_dir.glob(pattern)))
        excluded_suffixes = (
            ".unsorted.bam",
            ".tmp.bam",
            ".subsampled_40x.bam",
            ".cont_deam.bam",
            ".schmutzi.bam",
            ".schmutzi_ready.bam",
        )

        for b_path in all_bams:
            if any(b_path.name.endswith(sfx) for sfx in excluded_suffixes):
                continue

            stem = b_path.name
            cleaned_id = re.sub(r"(\.filtered|\.dedup|\.markdup|\.sorted)?\.bam$", "", stem)
            if cleaned_id in seen_sample_ids:
                logger.warning(f"Duplicate sample ID '{cleaned_id}' detected from {b_path.name}, skipping.")
                continue

            discovered.append((cleaned_id, b_path))
            seen_sample_ids.add(cleaned_id)

        logger.info(f"Discovered {len(discovered)} unique BAM files in {input_dir}")
        return discovered

    @classmethod
    def _run_phase1_screen_task(
        cls,
        sample_id: str,
        bam_path: Path,
        reference_fasta: Path,
        reference_contig: str,
        policy: FastScreeningPolicy,
        output_dir: Path,
    ) -> Tuple[str, Optional[FastScreenResult], Optional[Path], float, Optional[str]]:
        """Worker task executing Phase 1 Fast Screening."""
        t0 = time.perf_counter()
        sample_out = output_dir / sample_id
        sample_out.mkdir(parents=True, exist_ok=True)
        try:
            fs_res, dissent_tsv = FastContaminationScreen.run_screen(
                sample_id=sample_id,
                input_bam=bam_path,
                reference_fasta=reference_fasta,
                expected_contig=reference_contig,
                policy=policy,
                output_dir=sample_out,
            )
            return sample_id, fs_res, dissent_tsv, time.perf_counter() - t0, None
        except Exception as exc:
            return sample_id, None, None, time.perf_counter() - t0, str(exc)

    @classmethod
    def _invoke_dissent_engine(
        cls,
        sample_id: str,
        input_bam: Path,
        reference_fasta: Path,
        reference_contig: str,
        output_dir: Path,
        orchestrator: ContaminationOrchestrator,
    ) -> Tuple[ContaminationEstimate, ModelQualityDiagnostics, Optional[Path]]:
        """Safely invokes DissentEngine handling dynamic signatures and missing policy attributes."""
        sig = inspect.signature(DissentEngine.run)
        d_kwargs: Dict[str, Any] = {
            "sample_id": sample_id,
            "input_bam": input_bam,
            "reference_fasta": reference_fasta,
            "output_dir": output_dir,
        }

        if "expected_contig" in sig.parameters:
            d_kwargs["expected_contig"] = reference_contig
        elif "reference_contig" in sig.parameters:
            d_kwargs["reference_contig"] = reference_contig

        if "policy" in sig.parameters:
            pol = getattr(orchestrator.policy, "dissent", None)
            if pol is None:
                d_module = sys.modules.get("src.modules.contamination.dissent_engine")
                if d_module and hasattr(d_module, "DissentPolicy"):
                    pol = getattr(d_module, "DissentPolicy")()
            if pol is not None:
                d_kwargs["policy"] = pol
            elif sig.parameters["policy"].default is inspect.Parameter.empty:
                d_kwargs["policy"] = None

        return DissentEngine.run(**d_kwargs)

    @classmethod
    def _run_phase2_schmutzi_task(
        cls,
        context: ContaminationInputContext,
        output_dir: Path,
        config: Dict[str, Any],
        threads_per_sample: int,
        precomputed_fs: Optional[FastScreenResult],
        precomputed_tsv: Optional[Path],
        is_ancient: bool = True,
    ) -> ContaminationResult:
        """Worker task executing Phase 2 Bayesian modeling (or direct DissentEngine for modern mtDNA)."""
        t_task_start = time.perf_counter()
        orchestrator = ContaminationOrchestrator(config=config, project_root=PROJECT_ROOT)
        orchestrator.policy.schmutzi.threads = threads_per_sample
        orchestrator.policy.fast_screening.enabled = False

        sample_id = context.sample_id
        sample_out_dir = output_dir / sample_id
        sample_out_dir.mkdir(parents=True, exist_ok=True)
        final_json_path = sample_out_dir / f"{sample_id}.contamination.json"

        result = ContaminationResult(sample_id=sample_id, policy=orchestrator.policy)
        result.mark_started()

        audit = ContaminationAuditMetadata(
            input_bam_path=context.input_bam_path,
            input_bam_sha256=orchestrator._compute_sha256(context.input_bam_path),
            is_unrescaled_bam_verified=context.is_unrescaled_bam_verified,
            reference_fasta_path=context.reference_fasta_path,
            reference_fasta_sha256=orchestrator._compute_sha256(context.reference_fasta_path),
            fast_screen_applied=True if precomputed_fs else False,
            fast_screen_result=precomputed_fs,
            full_bayesian_estimation_performed=is_ancient,
            screened_without_full_quantification=False,
        )
        outputs = ContaminationOutputFiles(dissent_tsv=precomputed_tsv)

        preflight_error = orchestrator._preflight_validation(context)
        if preflight_error:
            err_code, err_msg = preflight_error
            result.add_error(err_code, err_msg)
            result.audit = audit
            result.outputs = outputs
            result.mark_finished()
            result.runtime_seconds = round(time.perf_counter() - t_task_start, 2)
            orchestrator._save_contamination_json(result, final_json_path)
            return result

        # Iniezione obbligatoria di samtools calmd e subsampling atomico con SchmutziInputPreparer
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
                subsampling_policy=orchestrator.policy.subsampling,
                threads=threads_per_sample,
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
                        orchestrator.policy.subsampling.target_depth / orig_depth, 4
                    )
        except Exception as prep_exc:
            result.add_error(
                ContaminationErrorCode.SCHMUTZI_EXECUTION_FAILED,
                f"BAM preparation with samtools calmd failed: {prep_exc}",
            )
            result.audit = audit
            result.outputs = outputs
            result.mark_finished()
            result.runtime_seconds = round(time.perf_counter() - t_task_start, 2)
            orchestrator._save_contamination_json(result, final_json_path)
            return result

        t_engine0 = time.perf_counter()

        # BIFORCAZIONE ARCHITETTURALE: Routing Diretto a DissentEngine per Campioni Moderni/Zero-Damage
        if not is_ancient:
            logger.info(
                f"[{sample_id}] Direct routing to DissentEngine (Block 4: modern/zero-damage mtDNA certified)."
            )
            try:
                estimate, diagnostics, tsv_file = cls._invoke_dissent_engine(
                    sample_id=sample_id,
                    input_bam=effective_bam,
                    reference_fasta=context.reference_fasta_path,
                    reference_contig=context.reference_contig,
                    output_dir=sample_out_dir,
                    orchestrator=orchestrator,
                )
                audit.step_runtimes_seconds["dissent"] = round(time.perf_counter() - t_engine0, 2)
                audit.step_runtimes_seconds["schmutzi"] = audit.step_runtimes_seconds["dissent"]
                outputs.dissent_tsv = tsv_file
            except Exception as d_exc:
                logger.error(f"[{sample_id}] DissentEngine execution failed: {d_exc}", exc_info=True)
                result.add_error(
                    ContaminationErrorCode.SCHMUTZI_EXECUTION_FAILED,
                    f"DissentEngine failed: {d_exc}",
                )
                result.audit = audit
                result.outputs = outputs
                result.mark_finished()
                result.runtime_seconds = round(time.perf_counter() - t_task_start, 2)
                orchestrator._save_contamination_json(result, final_json_path)
                return result
        else:
            # Esecuzione Bayesiana Completa Schmutzi per campioni archeologici autentici
            try:
                schmutzi_out = SchmutziEngine.run(
                    sample_id=sample_id,
                    input_bam=effective_bam,
                    reference_fasta=context.reference_fasta_path,
                    database_dir=orchestrator.policy.schmutzi.contaminant_database,
                    output_dir=output_dir,
                    threads=threads_per_sample,
                    max_iterations=orchestrator.policy.schmutzi.max_iterations,
                    use_wrapped=orchestrator.policy.schmutzi.use_wrapped_reference,
                )
                audit.step_runtimes_seconds["schmutzi"] = round(time.perf_counter() - t_engine0, 2)
                outputs.cont_est_file = schmutzi_out.cont_est_file
                outputs.cont_deam_file = schmutzi_out.cont_deam_file
                outputs.endogenous_fasta = schmutzi_out.endogenous_fasta
                outputs.contaminant_fasta = schmutzi_out.contaminant_fasta
                estimate = schmutzi_out.estimate
                diagnostics = schmutzi_out.diagnostics
            except SchmutziExecutionError as exc:
                is_zero_damage_crash = any(k in str(exc).lower() for k in ("deam.config", "splitdeam", "-nan"))
                if orchestrator.policy.allow_fallback_to_dissent or is_zero_damage_crash:
                    logger.warning(
                        f"[{sample_id}] Schmutzi execution failed ({exc}); auto-recovering with DissentEngine."
                    )
                    result.add_warning(
                        ContaminationWarningCode.FALLBACK_TO_DISSENT,
                        f"Schmutzi failed ({exc}); dissent used.",
                    )
                    try:
                        estimate, diagnostics, tsv_file = cls._invoke_dissent_engine(
                            sample_id=sample_id,
                            input_bam=effective_bam,
                            reference_fasta=context.reference_fasta_path,
                            reference_contig=context.reference_contig,
                            output_dir=sample_out_dir,
                            orchestrator=orchestrator,
                        )
                        audit.step_runtimes_seconds["dissent"] = round(time.perf_counter() - t_engine0, 2)
                        audit.step_runtimes_seconds["schmutzi"] = audit.step_runtimes_seconds["dissent"]
                        audit.full_bayesian_estimation_performed = False
                        outputs.dissent_tsv = tsv_file
                        result.fallback_reason = str(exc)
                    except Exception as fb_exc:
                        logger.error(f"[{sample_id}] Fallback DissentEngine failed: {fb_exc}", exc_info=True)
                        result.add_error(
                            ContaminationErrorCode.SCHMUTZI_EXECUTION_FAILED,
                            f"Schmutzi and fallback Dissent failed: {fb_exc}",
                        )
                        result.audit = audit
                        result.outputs = outputs
                        result.mark_finished()
                        result.runtime_seconds = round(time.perf_counter() - t_task_start, 2)
                        orchestrator._save_contamination_json(result, final_json_path)
                        return result
                else:
                    result.add_error(ContaminationErrorCode.SCHMUTZI_EXECUTION_FAILED, str(exc))
                    result.audit = audit
                    result.outputs = outputs
                    result.mark_finished()
                    result.runtime_seconds = round(time.perf_counter() - t_task_start, 2)
                    orchestrator._save_contamination_json(result, final_json_path)
                    return result

        qc_result = ContaminationQCValidator.evaluate_sample(
            sample_id=sample_id,
            estimate=estimate,
            diagnostics=diagnostics,
            outputs=outputs,
            audit=audit,
            policy=orchestrator.policy,
            output_json_dir=None,
        )
        qc_result.mark_finished()
        qc_result.runtime_seconds = round(time.perf_counter() - t_task_start, 2)
        qc_result.outputs.contamination_json = final_json_path
        orchestrator._save_contamination_json(qc_result, final_json_path)
        return qc_result

    @classmethod
    def run_cohort_pipeline(
        cls,
        discovered_samples: List[Tuple[Any, ...]],
        reference_fasta: Path,
        reference_contig: str,
        reference_length_bp: int,
        output_dir: Path,
        config: Dict[str, Any],
        resource_alloc: ResourceAllocation,
    ) -> List[ContaminationResult]:
        """Orchestrates two-phase cohort execution maximizing concurrency and work conservation."""
        output_dir.mkdir(parents=True, exist_ok=True)
        final_results: List[ContaminationResult] = []

        sample_map: Dict[str, Path] = {}
        ancient_map: Dict[str, bool] = {}

        # Mappatura e rilevamento automatico dello stato antico/moderno da Blocco 4
        for item in discovered_samples:
            s_id = str(item[0])
            b_path = Path(item[1])
            sample_map[s_id] = b_path

            if len(item) >= 3:
                ancient_map[s_id] = bool(item[2])
            else:
                damage_json = PROJECT_ROOT / "results" / "damage" / s_id / f"{s_id}.damage.json"
                if damage_json.is_file():
                    try:
                        with open(damage_json, "r", encoding="utf-8") as f_d:
                            dam_data = json.load(f_d)
                            v = (
                                dam_data.get("profile", {}).get("authentication_verdict")
                                or dam_data.get("verdict")
                                or ""
                            )
                            ancient_map[s_id] = "AUTHENTIC_ANCIENT_SUPPORTED" in str(v)
                    except Exception:
                        ancient_map[s_id] = True
                else:
                    ancient_map[s_id] = True

        orchestrator_ref = ContaminationOrchestrator(config=config, project_root=PROJECT_ROOT)
        fast_policy = orchestrator_ref.policy.fast_screening

        logger.info("=" * 80)
        logger.info(f"▶ PHASE 1: HIGH-CONCURRENCY FAST SCREENING TRIAGE ({len(discovered_samples)} samples)")
        logger.info(f"[*] Dispatching on {resource_alloc.max_fast_screen_workers} lightweight workers...")
        logger.info("=" * 80)

        t_p1_start = time.perf_counter()
        phase1_results: Dict[str, Tuple[Optional[FastScreenResult], Optional[Path], float, Optional[str]]] = {}

        with ThreadPoolExecutor(max_workers=resource_alloc.max_fast_screen_workers) as p1_pool:
            future_to_id = {
                p1_pool.submit(
                    cls._run_phase1_screen_task,
                    s_id,
                    sample_map[s_id],
                    reference_fasta,
                    reference_contig,
                    fast_policy,
                    output_dir,
                ): s_id
                for s_id in sample_map
            }
            for fut in as_completed(future_to_id):
                s_id, fs_res, tsv_p, r_time, err = fut.result()
                phase1_results[s_id] = (fs_res, tsv_p, r_time, err)

        t_p1_elapsed = time.perf_counter() - t_p1_start
        logger.info(f"Phase 1 Triage completed in {t_p1_elapsed:.2f}s.")

        fast_clear_samples: List[Tuple[str, FastScreenResult, Optional[Path], float]] = []
        escalation_queue: List[Tuple[str, Optional[FastScreenResult], Optional[Path]]] = []

        for s_id, (fs_res, tsv_p, r_time, err) in phase1_results.items():
            if err:
                err_res = ContaminationResult(sample_id=s_id, policy=orchestrator_ref.policy)
                err_res.mark_started()
                err_res.add_error(ContaminationErrorCode.FAST_SCREEN_FAILED, f"Phase 1 failed: {err}")
                err_res.mark_finished()
                err_res.runtime_seconds = round(r_time, 2)
                final_results.append(err_res)
            elif fs_res and fs_res.decision == ScreeningDecision.FAST_CLEAR and fs_res.bypass_eligible:
                fast_clear_samples.append((s_id, fs_res, tsv_p, r_time))
            else:
                escalation_queue.append((s_id, fs_res, tsv_p))

        logger.info(
            f"Triage Decision Summary: {len(fast_clear_samples)} FAST_CLEAR (bypassed) | "
            f"{len(escalation_queue)} Escalated/Insufficient Data (queued for Phase 2)"
        )

        for s_id, fs_res, tsv_p, r_time in fast_clear_samples:
            sample_out_dir = output_dir / s_id
            final_json_path = sample_out_dir / f"{s_id}.contamination.json"
            audit = ContaminationAuditMetadata(
                input_bam_path=sample_map[s_id],
                input_bam_sha256=orchestrator_ref._compute_sha256(sample_map[s_id]),
                reference_fasta_path=reference_fasta,
                reference_fasta_sha256=orchestrator_ref._compute_sha256(reference_fasta),
                fast_screen_applied=True,
                fast_screen_result=fs_res,
                full_bayesian_estimation_performed=False,
                screened_without_full_quantification=True,
                step_runtimes_seconds={"fast_screen": round(r_time, 2)},
            )
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
            qc_res = ContaminationQCValidator.evaluate_sample(
                sample_id=s_id,
                estimate=estimate,
                diagnostics=diagnostics,
                outputs=ContaminationOutputFiles(dissent_tsv=tsv_p, contamination_json=final_json_path),
                audit=audit,
                policy=orchestrator_ref.policy,
                output_json_dir=None,
            )
            qc_res.mark_finished()
            qc_res.runtime_seconds = round(r_time, 2)
            orchestrator_ref._save_contamination_json(qc_res, final_json_path)
            final_results.append(qc_res)

        if escalation_queue:
            logger.info("=" * 80)
            logger.info(f"▶ PHASE 2: BOUNDED MODELING POOL ({len(escalation_queue)} samples)")
            logger.info(
                f"[*] Dispatching on {resource_alloc.max_schmutzi_workers} concurrent workers "
                f"({resource_alloc.threads_per_schmutzi_instance} threads/instance)..."
            )
            logger.info("=" * 80)

            t_p2_start = time.perf_counter()
            with ThreadPoolExecutor(max_workers=resource_alloc.max_schmutzi_workers) as p2_pool:
                future_to_sample = {}
                for s_id, fs_res, tsv_p in escalation_queue:
                    is_anc = ancient_map.get(s_id, True)
                    context = ContaminationInputContext(
                        sample_id=s_id,
                        input_bam_path=sample_map[s_id],
                        reference_fasta_path=reference_fasta,
                        reference_contig=reference_contig,
                        reference_length_bp=reference_length_bp,
                    )
                    fut = p2_pool.submit(
                        cls._run_phase2_schmutzi_task,
                        context=context,
                        output_dir=output_dir,
                        config=config,
                        threads_per_sample=resource_alloc.threads_per_schmutzi_instance,
                        precomputed_fs=fs_res,
                        precomputed_tsv=tsv_p,
                        is_ancient=is_anc,
                    )
                    future_to_sample[fut] = s_id

                completed_p2 = 0
                for fut in as_completed(future_to_sample):
                    s_id = future_to_sample[fut]
                    completed_p2 += 1
                    try:
                        res = fut.result()
                        final_results.append(res)
                        runtime_str = f"{res.runtime_seconds:.2f}" if res.runtime_seconds is not None else "0.00"
                        logger.info(
                            f"[{completed_p2}/{len(escalation_queue)}] Phase 2 completed: {s_id} in {runtime_str}s "
                            f"| Method: {res.estimate.method_used.value} | Status: {res.status.value} "
                            f"| Est: {res.estimate.point_estimate:.2%}"
                        )
                    except Exception as exc:
                        logger.critical(f"[{completed_p2}/{len(escalation_queue)}] Fatal failure on {s_id}: {exc}")
                        err_res = ContaminationResult(sample_id=s_id, policy=orchestrator_ref.policy)
                        err_res.mark_started()
                        err_res.add_error(ContaminationErrorCode.SCHMUTZI_EXECUTION_FAILED, str(exc))
                        err_res.mark_finished()
                        final_results.append(err_res)

            t_p2_elapsed = time.perf_counter() - t_p2_start
            logger.info(f"Phase 2 modeling completed in {t_p2_elapsed:.2f}s.")

        return final_results

    @classmethod
    def export_cohort_summary(
        cls,
        results: List[ContaminationResult],
        output_dir: Path,
        resource_alloc: ResourceAllocation,
    ) -> Tuple[Path, Path]:
        """Generates audit-ready cohort-level TSV and JSON summaries with performance analytics."""
        output_dir.mkdir(parents=True, exist_ok=True)
        tsv_path = output_dir / "cohort_contamination_summary.tsv"
        json_path = output_dir / "cohort_contamination_summary.json"

        tsv_lines = [
            "# cohort_contamination_summary_v2.0\n",
            f"# date={datetime.now(timezone.utc).isoformat()}\n",
            f"# total_samples={len(results)}\n",
            f"# host_cores_allocated={resource_alloc.total_allocated_threads}/{resource_alloc.total_host_cores}\n",
            f"# schmutzi_max_workers={resource_alloc.max_schmutzi_workers}\n",
            "sample_id\tstatus\tscreening_decision\tbypass_applied\tschmutzi_executed\t"
            "contamination_category\trecommended_action\tpoint_estimate\tci_lower_95\tci_upper_95\t"
            "method_used\tdiagnostic_depth\tgenome_depth\tdissenting_sites\tfast_screen_runtime_sec\t"
            "schmutzi_runtime_sec\ttotal_runtime_sec\tsubsampling_applied\terror_stage\terrors\n",
        ]

        summary_records: List[Dict[str, Any]] = []
        bypassed_count = 0
        schmutzi_count = 0

        for r in sorted(results, key=lambda x: x.sample_id):
            fs_res = r.audit.fast_screen_result
            screening_decision = fs_res.decision.value if fs_res else "NOT_SCREENED"
            bypass_applied = bool(r.audit.screened_without_full_quantification)
            schmutzi_executed = bool(r.audit.full_bayesian_estimation_performed)

            if bypass_applied:
                bypassed_count += 1
            if schmutzi_executed:
                schmutzi_count += 1

            diag_depth = fs_res.mean_diagnostic_depth if fs_res else 0.0
            genome_depth = fs_res.mean_genome_depth if fs_res else r.diagnostics.mean_depth
            dissent_sites = fs_res.dissenting_positions_count if fs_res else 0

            fs_time = r.audit.step_runtimes_seconds.get("fast_screen", 0.0)
            sch_time = r.audit.step_runtimes_seconds.get("schmutzi", r.audit.step_runtimes_seconds.get("dissent", 0.0))

            err_stage = ""
            if r.status == ContaminationStatus.FAIL:
                if any("FAST_SCREEN" in e.value for e in r.error_codes):
                    err_stage = "FAST_SCREEN"
                elif any("SCHMUTZI" in e.value for e in r.error_codes):
                    err_stage = "SCHMUTZI"
                elif any("DISSENT" in e.value for e in r.error_codes):
                    err_stage = "DISSENT"
                else:
                    err_stage = "VALIDATION"

            rec = BatchSampleRecord(
                sample_id=r.sample_id,
                status=r.status.value,
                screening_decision=screening_decision,
                bypass_applied=bypass_applied,
                schmutzi_executed=schmutzi_executed,
                contamination_category=r.gating.contamination_category.value,
                recommended_action=r.gating.recommended_action.value,
                point_estimate=r.estimate.point_estimate,
                ci_lower_95=r.estimate.ci_lower_95,
                ci_upper_95=r.estimate.ci_upper_95,
                method_used=r.estimate.method_used.value,
                diagnostic_depth=diag_depth,
                genome_depth=genome_depth,
                dissenting_sites=dissent_sites,
                fast_screen_runtime_seconds=fs_time,
                schmutzi_runtime_seconds=sch_time,
                total_runtime_seconds=r.runtime_seconds or (fs_time + sch_time),
                subsampling_applied=r.audit.subsampling_applied,
                error_stage=err_stage,
                error_summary="; ".join(r.errors) if r.errors else "",
            )
            tsv_lines.append(rec.to_tsv_row())
            summary_records.append(r.to_dict())

        with open(tsv_path, "w", encoding="utf-8") as f_tsv:
            f_tsv.writelines(tsv_lines)

        bypass_pct = (bypassed_count / len(results) * 100.0) if results else 0.0

        cohort_payload = {
            "metadata": {
                "generated_at": datetime.now(timezone.utc).isoformat(),
                "total_samples": len(results),
                "fast_clear_bypassed": bypassed_count,
                "schmutzi_executed": schmutzi_count,
                "bypass_efficiency_pct": round(bypass_pct, 2),
                "resource_allocation": resource_alloc.to_dict(),
            },
            "samples": summary_records,
        }

        with open(json_path, "w", encoding="utf-8") as f_json:
            json.dump(cohort_payload, f_json, indent=2)

        logger.info(
            f"Cohort summaries generated -> TSV: {tsv_path.name} | JSON: {json_path.name} "
            f"| Bypass Efficiency: {bypassed_count}/{len(results)} ({bypass_pct:.1f}%)"
        )
        return tsv_path, json_path


def main() -> None:
    """CLI interface for running two-phase multi-sample cohort contamination estimation."""
    config_data = ContaminationBatchRunner.load_yaml_config()

    ref_cfg = config_data.get("reference", {})
    contam_cfg = config_data.get("contamination", {})
    schmutzi_cfg = contam_cfg.get("schmutzi", {})
    out_cfg = config_data.get("output", {})

    default_ref = ContaminationBatchRunner.resolve_path(
        ref_cfg.get("fasta_path"),
        PROJECT_ROOT / "data" / "reference" / "rCRS.fasta",
    )
    default_contig = ref_cfg.get("expected_contig", "NC_012920.1")
    default_length = ref_cfg.get("expected_length", 16569)

    default_input_dir = ContaminationBatchRunner.resolve_path(
        "results/alignment",
        PROJECT_ROOT / "results" / "alignment",
    )
    default_out_dir = ContaminationBatchRunner.resolve_path(
        out_cfg.get("results_dir", "results") + "/contamination",
        PROJECT_ROOT / "results" / "contamination",
    )

    parser = argparse.ArgumentParser(
        description="Run two-phase batch mtDNA contamination estimation with dynamic CPU budgeting."
    )
    parser.add_argument("--config", type=Path, default=PROJECT_ROOT / "config" / "default_config.yaml", help="Path to YAML configuration.")
    parser.add_argument("--input-dir", type=Path, default=default_input_dir, help=f"Directory with input BAM files (default: {default_input_dir}).")
    parser.add_argument("--bam-pattern", type=str, default="*.bam", help="Glob pattern for discovering input BAMs (default: *.bam).")
    parser.add_argument("--sample-manifest", type=Path, default=None, help="Optional TSV/CSV manifest mapping sample_id to BAM path.")
    parser.add_argument("--output-dir", type=Path, default=default_out_dir, help=f"Directory to save outputs (default: {default_out_dir}).")
    parser.add_argument("--reference-fasta", type=Path, default=default_ref, help=f"Reference FASTA path (default: {default_ref}).")
    parser.add_argument("--contig", type=str, default=default_contig, help=f"Reference contig name (default: {default_contig}).")
    parser.add_argument("--contig-length", type=int, default=default_length, help=f"Reference contig length bp (default: {default_length}).")
    parser.add_argument("--cpu-budget-fraction", type=float, default=0.50, help="Fraction of host CPU cores allocated to batch (default: 0.50).")
    parser.add_argument("--threads-per-sample", type=int, default=schmutzi_cfg.get("threads", 4), help="Threads per Schmutzi instance (default: 4).")
    parser.add_argument("--max-workers", type=int, default=None, help="Explicit maximum concurrent Schmutzi workers (bounded by budget).")
    parser.add_argument("--host-cores", type=int, default=None, help="Explicit host core override for container/HPC jobs.")
    parser.add_argument("--allow-oversubscription", action="store_true", help="Allow max-workers to exceed calculated CPU budget.")

    args = parser.parse_args()

    resource_alloc = ContaminationBatchRunner.calculate_resource_allocation(
        threads_per_sample=args.threads_per_sample,
        cpu_budget_fraction=args.cpu_budget_fraction,
        explicit_max_workers=args.max_workers,
        host_cores_override=args.host_cores,
        allow_oversubscription=args.allow_oversubscription,
    )

    discovered = ContaminationBatchRunner.discover_bam_inputs(
        input_dir=args.input_dir,
        pattern=args.bam_pattern,
        sample_manifest=args.sample_manifest,
    )

    if not discovered:
        logger.warning("No valid BAM files found. Exiting.")
        sys.exit(0)

    results = ContaminationBatchRunner.run_cohort_pipeline(
        discovered_samples=discovered,
        reference_fasta=args.reference_fasta,
        reference_contig=args.contig,
        reference_length_bp=args.contig_length,
        output_dir=args.output_dir,
        config=config_data,
        resource_alloc=resource_alloc,
    )

    ContaminationBatchRunner.export_cohort_summary(
        results=results,
        output_dir=args.output_dir,
        resource_alloc=resource_alloc,
    )


if __name__ == "__main__":
    main()