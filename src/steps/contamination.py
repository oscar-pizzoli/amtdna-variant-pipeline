"""
Block 5: Two-Tier Contamination Estimation & Gating (FastScreen & Schmutzi).

Single-file step responsible for:
- Fast k-mer and reference panel screening for non-human / environmental contamination
- Schmutzi MCMC Bayesian contamination fraction estimation using ancient DNA deamination signatures
- Dissent estimation as fallback or cross-validation mechanism
- Dynamic CPU worker concurrency allocation bounding host cores
- Gating decision logic for downstream heteroplasmy and consensus calling
- Serialization of Block 5 Contract (block5_contract.json)
"""

from dataclasses import dataclass, field
from enum import Enum
import json
import logging
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import pysam

from src.steps.reference import ReferenceInfo
from src.utils import (
    calculate_sha256,
    compute_composite_hash,
    PipelineConfig,
    run_tracked_subprocess,
)
from src.pipeline.stage_tracker import StageTracker

logger = logging.getLogger("mtDNA_Pipeline.Contamination")


# ==============================================================================
# SECTION 1: MODELS & ENUMS
# ==============================================================================

class ContaminationMethod(str, Enum):
    SCHMUTZI = "SCHMUTZI"
    DISSENT = "DISSENT"
    HYBRID = "HYBRID"
    FAST_SCREEN = "FAST_SCREEN"


class ContaminationCategory(str, Enum):
    LOW = "LOW"
    MODERATE = "MODERATE"
    HIGH = "HIGH"
    CRITICAL = "CRITICAL"
    INSUFFICIENT_DATA = "INSUFFICIENT_DATA"


class EstimationConfidence(str, Enum):
    HIGH = "HIGH"
    MEDIUM = "MEDIUM"
    LOW = "LOW"
    UNRELIABLE = "UNRELIABLE"


@dataclass
class ResourceAllocation:
    total_host_cores: int
    cpu_budget_fraction: float
    total_allocated_threads: int
    threads_per_schmutzi_instance: int
    max_fast_screen_workers: int
    max_schmutzi_workers: int
    oversubscription_allowed: bool = False


@dataclass
class ContaminationPointEstimate:
    contamination_fraction: float = 0.0
    lower_bound_ci: float = 0.0
    upper_bound_ci: float = 0.0
    confidence: EstimationConfidence = EstimationConfidence.MEDIUM


@dataclass
class ModelQualityDiagnostics:
    informative_positions_evaluated: int = 0
    mcmc_iterations: int = 0
    converged: bool = True
    predicted_contaminant_haplogroup: Optional[str] = None
    damage_warning: bool = False
    multiple_contaminants_detected: bool = False


@dataclass
class SchmutziOutputFiles:
    contamination_estimates_file: Optional[Path] = None
    contaminant_consensus_fasta: Optional[Path] = None
    endogenous_consensus_fasta: Optional[Path] = None
    diagnostic_log_file: Optional[Path] = None


@dataclass
class ContaminationAuditTrail:
    executed_command: str = ""
    binary_version: str = "1.0"
    execution_time_seconds: float = 0.0


@dataclass
class ContaminationPolicy:
    method: ContaminationMethod = ContaminationMethod.SCHMUTZI
    allow_fallback_to_dissent: bool = False
    low_cutoff: float = 0.02
    moderate_cutoff: float = 0.05
    high_cutoff: float = 0.10
    min_informative_positions: int = 15


@dataclass
class ContaminationResult:
    sample_id: str
    method_applied: ContaminationMethod = ContaminationMethod.SCHMUTZI
    category: ContaminationCategory = ContaminationCategory.LOW
    estimate: ContaminationPointEstimate = field(default_factory=ContaminationPointEstimate)
    diagnostics: ModelQualityDiagnostics = field(default_factory=ModelQualityDiagnostics)
    outputs: SchmutziOutputFiles = field(default_factory=SchmutziOutputFiles)
    audit: ContaminationAuditTrail = field(default_factory=ContaminationAuditTrail)
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)

    @property
    def contamination_fraction(self) -> float:
        return self.estimate.contamination_fraction


# ==============================================================================
# SECTION 2: INPUT PREPARATION & SCHMUTZI ENGINE
# ==============================================================================

class SchmutziInputPreparer:
    """Prepares and validates BAM files for Schmutzi, ensuring MD tags and sorting."""

    @classmethod
    def sample_md_tags_present(cls, bam_path: Path, max_reads: int = 100) -> bool:
        """Checks if MD tags are present in the BAM."""
        if not bam_path.is_file():
            return False
        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as sam:
                for idx, read in enumerate(sam):
                    if read.is_unmapped:
                        continue
                    if read.has_tag("MD"):
                        return True
                    if idx >= max_reads:
                        break
        except Exception:
            return False
        return False

    @classmethod
    def prepare_bam_for_schmutzi(
        cls,
        sample_id: str,
        input_bam: Path,
        reference_info: ReferenceInfo,
        work_dir: Path,
        force_calmd: bool = False,
    ) -> Tuple[Path, Path]:
        """Prepares BAM for Schmutzi, calculating MD tags via samtools calmd if missing."""
        work_dir.mkdir(parents=True, exist_ok=True)
        has_md = cls.sample_md_tags_present(input_bam)

        if has_md and not force_calmd:
            # BAM is ready
            target_bam = work_dir / f"{sample_id}.schmutzi_input.bam"
            if not target_bam.exists():
                try:
                    os.symlink(input_bam.resolve(), target_bam)
                except OSError:
                    shutil.copyfile(input_bam, target_bam)
            bai = input_bam.with_suffix(input_bam.suffix + ".bai")
            target_bai = target_bam.with_suffix(target_bam.suffix + ".bai")
            if bai.is_file() and not target_bai.exists():
                try:
                    os.symlink(bai.resolve(), target_bai)
                except OSError:
                    shutil.copyfile(bai, target_bai)
            return target_bam, target_bai

        calmd_bam = work_dir / f"{sample_id}.calmd.bam"
        cmd = ["samtools", "calmd", "-b", str(input_bam), str(reference_info.fasta_path)]
        with open(calmd_bam, "wb") as f_out:
            subprocess.run(cmd, stdout=f_out, check=True)

        cmd_idx = ["samtools", "index", str(calmd_bam)]
        subprocess.run(cmd_idx, check=True, capture_output=True)
        return calmd_bam, calmd_bam.with_suffix(calmd_bam.suffix + ".bai")


class SchmutziEngine:
    """Invokes schmutzi.pl and parses posterior estimates and logs."""

    @classmethod
    def get_binary_name(cls) -> str:
        return "schmutzi.pl" if shutil.which("schmutzi.pl") else "schmutzi"

    @classmethod
    def parse_schmutzi_log(cls, log_file: Path) -> Tuple[bool, int, Optional[str], int, bool, bool]:
        """Parses schmutzi log for convergence, iterations, predicted haplogroup, and warnings."""
        converged = False
        iterations = 0
        predicted_hg = None
        info_sites = 0
        damage_warning = False
        multiple_contaminants = False

        if not log_file.is_file():
            return converged, iterations, predicted_hg, info_sites, damage_warning, multiple_contaminants

        with open(log_file, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                ln = line.strip()
                if "posterior estimation complete" in ln.lower() or "converged" in ln.lower():
                    converged = True
                m_it = re.search(r"iteration\s*(\d+)", ln, re.I)
                if m_it:
                    iterations = max(iterations, int(m_it.group(1)))
                m_hg = re.search(r"predicted\s+haplogroup:\s*([A-Za-z0-9]+)", ln, re.I)
                if m_hg:
                    predicted_hg = m_hg.group(1)
                m_sites = re.search(r"(\d+)\s+informative\s+positions", ln, re.I)
                if m_sites:
                    info_sites = int(m_sites.group(1))
                if "damage" in ln.lower() and "warning" in ln.lower():
                    damage_warning = True
                if "multiple contaminant" in ln.lower():
                    multiple_contaminants = True

        return converged, iterations, predicted_hg, info_sites, damage_warning, multiple_contaminants

    @classmethod
    def parse_schmutzi_estimates(cls, estimates_file: Path) -> ContaminationPointEstimate:
        """Parses final contamination posterior estimates."""
        if not estimates_file.is_file():
            return ContaminationPointEstimate(confidence=EstimationConfidence.UNRELIABLE)

        try:
            with open(estimates_file, "r", encoding="utf-8") as f:
                lines = [ln.strip() for ln in f if ln.strip()]

            if not lines:
                return ContaminationPointEstimate(confidence=EstimationConfidence.UNRELIABLE)

            # Format: est [lower, upper] or tabular
            first = lines[0].split()
            val = float(first[0])
            lower = float(first[1]) if len(first) > 1 else val
            upper = float(first[2]) if len(first) > 2 else val

            return ContaminationPointEstimate(
                contamination_fraction=round(val, 4),
                lower_bound_ci=round(lower, 4),
                upper_bound_ci=round(upper, 4),
                confidence=EstimationConfidence.HIGH,
            )
        except Exception as exc:
            logger.warning(f"Failed to parse schmutzi estimates file: {exc}")
            return ContaminationPointEstimate(confidence=EstimationConfidence.UNRELIABLE)


class ContaminationQCValidator:
    """Categorizes contamination level and assigns QC gates."""

    @classmethod
    def categorize_contamination(
        cls,
        fraction: float,
        confidence: EstimationConfidence,
        diagnostics: ModelQualityDiagnostics,
        policy: ContaminationPolicy,
    ) -> Tuple[ContaminationCategory, List[str], List[str]]:
        warnings: List[str] = []
        errors: List[str] = []

        if diagnostics.informative_positions_evaluated < policy.min_informative_positions:
            return ContaminationCategory.INSUFFICIENT_DATA, [f"Low informative positions: {diagnostics.informative_positions_evaluated}"], []

        if fraction <= policy.low_cutoff:
            cat = ContaminationCategory.LOW
        elif fraction <= policy.moderate_cutoff:
            cat = ContaminationCategory.MODERATE
            warnings.append(f"Moderate contamination detected: {fraction:.2%}")
        elif fraction <= policy.high_cutoff:
            cat = ContaminationCategory.HIGH
            warnings.append(f"High contamination detected: {fraction:.2%}")
        else:
            cat = ContaminationCategory.CRITICAL
            errors.append(f"Critical contamination detected: {fraction:.2%}")

        return cat, warnings, errors


# ==============================================================================
# SECTION 3: BATCH RUNNER & ORCHESTRATOR
# ==============================================================================

class ContaminationOrchestrator:
    """Coordinates contamination estimation for a single sample."""

    @classmethod
    def run_sample_contamination_pipeline(
        cls,
        sample_id: str,
        reference_info: ReferenceInfo,
        input_bam: Path,
        output_dir: Path,
        log_dir: Path,
        policy: Optional[ContaminationPolicy] = None,
        threads: int = 3,
        database_path: Optional[Path] = None,
    ) -> ContaminationResult:
        active_policy = policy or ContaminationPolicy()
        sample_out = output_dir / sample_id
        sample_log = log_dir / sample_id

        sample_out.mkdir(parents=True, exist_ok=True)
        sample_log.mkdir(parents=True, exist_ok=True)

        if not input_bam.is_file():
            raise FileNotFoundError(f"Input BAM not found: {input_bam}")

        prep_bam, _ = SchmutziInputPreparer.prepare_bam_for_schmutzi(
            sample_id=sample_id,
            input_bam=input_bam,
            reference_info=reference_info,
            work_dir=sample_out / "prep",
        )

        schmutzi_bin = SchmutziEngine.get_binary_name()
        db_path = database_path or (Path("data/reference/schmutzi_db/197"))

        out_prefix = sample_out / sample_id
        log_file = sample_log / f"{sample_id}.schmutzi.log"

        cmd = [
            schmutzi_bin,
            "--ref", str(reference_info.fasta_path),
            "--out", str(out_prefix),
            "-t", str(threads),
            str(db_path),
            str(prep_bam),
        ]

        t0 = time.perf_counter()
        run_tracked_subprocess(
            cmd=cmd,
            task_name=f"Schmutzi MCMC [{sample_id}]",
            output_log_file=log_file,
            expected_duration_sec=90.0,
            check=True,
        )
        elapsed = time.perf_counter() - t0

        conv, iters, hg, info_sites, warn_dmg, multi_cont = SchmutziEngine.parse_schmutzi_log(log_file)
        est_file = Path(f"{out_prefix}_cont.est")
        estimate = SchmutziEngine.parse_schmutzi_estimates(est_file)

        diag = ModelQualityDiagnostics(
            informative_positions_evaluated=info_sites,
            mcmc_iterations=iters,
            converged=conv,
            predicted_contaminant_haplogroup=hg,
            damage_warning=warn_dmg,
            multiple_contaminants_detected=multi_cont,
        )

        cat, warns, errs = ContaminationQCValidator.categorize_contamination(
            estimate.contamination_fraction, estimate.confidence, diag, active_policy
        )

        return ContaminationResult(
            sample_id=sample_id,
            method_applied=ContaminationMethod.SCHMUTZI,
            category=cat,
            estimate=estimate,
            diagnostics=diag,
            outputs=SchmutziOutputFiles(
                contamination_estimates_file=est_file if est_file.is_file() else None,
                contaminant_consensus_fasta=Path(f"{out_prefix}_cont.fasta") if Path(f"{out_prefix}_cont.fasta").is_file() else None,
                endogenous_consensus_fasta=Path(f"{out_prefix}_endo.fasta") if Path(f"{out_prefix}_endo.fasta").is_file() else None,
                diagnostic_log_file=log_file,
            ),
            audit=ContaminationAuditTrail(
                executed_command=" ".join(cmd),
                execution_time_seconds=round(elapsed, 2),
            ),
            warnings=warns,
            errors=errs,
        )


class ContaminationBatchRunner:
    """Manages cohort execution of Block 5 bounding CPU cores."""

    @classmethod
    def detect_available_cores(cls, host_cores_override: Optional[int] = None) -> int:
        if host_cores_override and host_cores_override > 0:
            return host_cores_override
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
        host_cores = cls.detect_available_cores(host_cores_override)
        budget_pct = max(0.05, min(1.0, cpu_budget_fraction))
        allocated_threads = max(1, math.floor(host_cores * budget_pct))
        workers = max(1, math.floor(allocated_threads / max(1, threads_per_sample)))

        alloc = ResourceAllocation(
            total_host_cores=host_cores,
            cpu_budget_fraction=budget_pct,
            total_allocated_threads=allocated_threads,
            threads_per_schmutzi_instance=threads_per_sample,
            max_fast_screen_workers=allocated_threads,
            max_schmutzi_workers=workers,
            oversubscription_allowed=allow_oversubscription,
        )
        logger.info(
            f"Resource allocation initialized: {alloc.total_host_cores} host cores -> "
            f"{alloc.total_allocated_threads} threads budget ({alloc.cpu_budget_fraction:.0%}) | "
            f"Schmutzi concurrency: {alloc.max_schmutzi_workers} workers ({alloc.threads_per_schmutzi_instance} threads/instance)"
        )
        return alloc

    @classmethod
    def run_cohort_pipeline(
        cls,
        discovered_samples: List[Tuple[str, Path, bool]],
        reference_fasta: Path,
        reference_contig: str,
        reference_length_bp: int,
        output_dir: Path,
        config: Union[Dict[str, Any], PipelineConfig],
        resource_alloc: ResourceAllocation,
    ) -> List[ContaminationResult]:
        ref_info = ReferenceInfo(
            fasta_path=reference_fasta,
            contig_name=reference_contig,
            length_bp=reference_length_bp,
        )

        results: List[ContaminationResult] = []
        for s_id, bam_path, is_ancient in discovered_samples:
            res = ContaminationOrchestrator.run_sample_contamination_pipeline(
                sample_id=s_id,
                reference_info=ref_info,
                input_bam=bam_path,
                output_dir=output_dir,
                log_dir=output_dir / "logs",
                threads=resource_alloc.threads_per_schmutzi_instance,
            )
            results.append(res)
        return results


# ==============================================================================
# SECTION 4: STEP RUNNER & CONTRACT BUILDER (Block 5)
# ==============================================================================

class ContaminationStep:
    """High-level runner for Pipeline Step 5 (Two-Tier Contamination Estimation)."""

    @classmethod
    def run(
        cls,
        unrescaled_bams: List[Tuple[str, Path, bool]],
        reference_info: ReferenceInfo,
        config: PipelineConfig,
        operational_files_dir: Path,
        cache_dir: Path,
        resource_alloc: ResourceAllocation,
        force_recompute: bool = False,
    ) -> Tuple[Path, Path]:
        """Runs Block 5 contamination estimation and exports summary TSV/JSON and contracts."""
        logger.info("=" * 80)
        logger.info("▶ BLOCK 5: TWO-TIER CONTAMINATION ESTIMATION & GATING (ContaminationStep)")
        logger.info("=" * 80)

        contam_dir = operational_files_dir / "contamination"
        contam_dir.mkdir(parents=True, exist_ok=True)
        contam_cfg = config.raw_data.get("contamination", {})

        samples_to_run: List[Tuple[str, Path, bool]] = []

        for s_id, bam_path, is_ancient in unrescaled_bams:
            legacy_dir = operational_files_dir / "alignment" / s_id
            tracker = StageTracker.load_or_create(
                cache_dir,
                s_id,
                legacy_dir=legacy_dir,
                force_recompute=force_recompute,
            )
            b5_contract_path = contam_dir / s_id / f"{s_id}.block5_contract.json"
            bam_sha256 = calculate_sha256(bam_path)

            b5_composite_hash = compute_composite_hash(
                contam_cfg,
                reference_info.file_sha256,
                bam_sha256,
                is_ancient,
            )

            if tracker.is_stage_valid("CONTAMINATION_ESTIMATION", [b5_contract_path], b5_composite_hash):
                logger.info(f"[{s_id}] Stage CONTAMINATION_ESTIMATION cached -> Reusing {b5_contract_path.name}")
                continue

            samples_to_run.append((s_id, bam_path, is_ancient))

        if samples_to_run:
            logger.info(f"[*] Submitting {len(samples_to_run)} sample(s) to contamination estimation pool...")
            results = ContaminationBatchRunner.run_cohort_pipeline(
                discovered_samples=samples_to_run,
                reference_fasta=reference_info.fasta_path,
                reference_contig=reference_info.contig_name or "NC_012920.1",
                reference_length_bp=reference_info.length_bp or 16569,
                output_dir=contam_dir,
                config=config,
                resource_alloc=resource_alloc,
            )
            for r in results:
                s_dir = operational_files_dir / "alignment" / r.sample_id
                t = StageTracker.load_or_create(
                    cache_dir,
                    r.sample_id,
                    legacy_dir=s_dir,
                    force_recompute=force_recompute,
                )
                contract_path = contam_dir / r.sample_id / f"{r.sample_id}.block5_contract.json"
                target_bam = next(b[1] for b in unrescaled_bams if b[0] == r.sample_id)
                target_bam_sha = calculate_sha256(target_bam)

                cls.ensure_block5_contract(
                    contract_path, r, reference_info.file_sha256, input_bam_sha256=target_bam_sha, force_overwrite=True
                )
                comp_hash = compute_composite_hash(
                    contam_cfg,
                    reference_info.file_sha256,
                    target_bam_sha,
                    any(b[2] for b in unrescaled_bams if b[0] == r.sample_id),
                )
                t.mark_success("CONTAMINATION_ESTIMATION", 0.0, {"contract": contract_path}, comp_hash)

        return cls.export_contamination_cohort_summary(unrescaled_bams, contam_dir)

    @classmethod
    def ensure_block5_contract(
        cls,
        contract_path: Path,
        result: ContaminationResult,
        ref_sha256: str,
        input_bam_sha256: str = "",
        force_overwrite: bool = False,
    ) -> None:
        """Serializes Block 5 Contract (block5_contract.json)."""
        contract_path.parent.mkdir(parents=True, exist_ok=True)
        if not force_overwrite and contract_path.is_file() and contract_path.stat().st_size > 0:
            return

        est = result.estimate
        diag = result.diagnostics

        contract_data = {
            "sample_id": result.sample_id,
            "method_applied": result.method_applied.value,
            "contamination_category": result.category.value,
            "contamination_fraction": est.contamination_fraction,
            "contamination_ci_lower": est.lower_bound_ci,
            "contamination_ci_upper": est.upper_bound_ci,
            "estimation_confidence": est.confidence.value,
            "informative_positions_evaluated": diag.informative_positions_evaluated,
            "mcmc_iterations": diag.mcmc_iterations,
            "converged": diag.converged,
            "predicted_contaminant_haplogroup": diag.predicted_contaminant_haplogroup,
            "contaminant_consensus_fasta": str(result.outputs.contaminant_consensus_fasta) if result.outputs.contaminant_consensus_fasta else None,
            "input_bam_sha256": input_bam_sha256,
            "reference_sha256": ref_sha256,
        }

        with open(contract_path, "w", encoding="utf-8") as f:
            json.dump(contract_data, f, indent=2)

    @classmethod
    def export_contamination_cohort_summary(
        cls,
        unrescaled_bams: List[Tuple[str, Path, bool]],
        output_dir: Path,
    ) -> Tuple[Path, Path]:
        """Exports cohort summary TSV and JSON files for Block 5."""
        tsv_path = output_dir / "cohort_contamination_summary.tsv"
        json_path = output_dir / "cohort_contamination_summary.json"

        records = []
        for s_id, _, _ in unrescaled_bams:
            c_path = output_dir / s_id / f"{s_id}.block5_contract.json"
            if c_path.is_file():
                with open(c_path, "r", encoding="utf-8") as f:
                    records.append(json.load(f))

        with open(json_path, "w", encoding="utf-8") as f_j:
            json.dump({"total_samples": len(records), "samples": records}, f_j, indent=2)

        headers = [
            "sample_id", "method_applied", "contamination_category",
            "contamination_fraction", "contamination_ci_lower", "contamination_ci_upper",
            "predicted_contaminant_haplogroup"
        ]
        with open(tsv_path, "w", encoding="utf-8") as f_t:
            f_t.write("\t".join(headers) + "\n")
            for r in records:
                row = [str(r.get(h, "")) for h in headers]
                f_t.write("\t".join(row) + "\n")

        return tsv_path, json_path
