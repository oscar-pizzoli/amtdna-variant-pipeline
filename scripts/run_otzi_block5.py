#!/usr/bin/env python3
"""Strict Integration Runner for Block 5: mtDNA Contamination Estimation (Ötzi).

Executes schmutzi Bayesian joint mixture modeling in STRICT MODE (allow_fallback_to_dissent=False).
Performs pre-flight BAM/reference integrity checks, enforces unrescaled input verification,
and validates downstream JSON deliverables and physical artifacts.
"""

import json
import logging
from pathlib import Path
import shutil
import subprocess
import sys
import time

import pysam

# Add repository root to Python module path
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.modules.alignment.models import ReferenceInfo
from src.modules.contamination.models import (
    ConsensusGatingAction,
    ContaminationCategory,
    ContaminationMethod,
    ContaminationPolicy,
    ContaminationStatus,
    DissentPolicy,
    EstimationConfidence,
    SchmutziPolicy,
)
from src.modules.contamination.orchestrator import ContaminationOrchestrator
from src.modules.damage.models import LibraryTreatment

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("mtDNA_Pipeline.Block5_StrictRunner")


def preflight_check(bam_path: Path, ref_path: Path, db_path: Path) -> ReferenceInfo:
    """Performs strict pre-flight validation of input BAM, reference FASTA, and database."""
    # 1. Reference FASTA Check
    if not ref_path.is_file() or ref_path.stat().st_size == 0:
        raise FileNotFoundError(f"Reference FASTA missing or empty: {ref_path}")

    with pysam.FastaFile(str(ref_path)) as fa:
        contigs = fa.references
        if not contigs:
            raise ValueError(f"Reference FASTA has no sequences: {ref_path}")
        target_contig = contigs[0]
        ref_length = fa.get_reference_length(target_contig)

    logger.info(f"[Pre-flight] Reference verified: '{target_contig}' ({ref_length:,} bp)")

    # 2. BAM Check & Integrity
    if not bam_path.is_file() or bam_path.stat().st_size == 0:
        raise FileNotFoundError(f"Input BAM missing or empty: {bam_path}")

    # Run quickcheck if samtools is available
    samtools_bin = shutil.which("samtools")
    if samtools_bin:
        qc_res = subprocess.run([samtools_bin, "quickcheck", str(bam_path)], check=False)
        if qc_res.returncode != 0:
            raise ValueError(f"Input BAM failed 'samtools quickcheck': {bam_path}")

    # Check header and mapped reads
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        if target_contig not in bam.references:
            raise ValueError(
                f"Contig '{target_contig}' missing from BAM header. Found: {bam.references}"
            )
        bam_len = bam.get_reference_length(target_contig)
        if bam_len != ref_length:
            raise ValueError(
                f"Length mismatch: BAM @SQ length ({bam_len} bp) != FASTA length ({ref_length} bp)"
            )
        mapped_mt = bam.count(contig=target_contig, read_callback="all")
        if mapped_mt == 0:
            raise ValueError(f"BAM contains 0 mapped reads on target contig '{target_contig}'.")

    logger.info(f"[Pre-flight] BAM verified: '{bam_path.name}' ({mapped_mt:,} mapped mtDNA reads)")

    # 3. Database Check
    if not db_path.exists():
        raise FileNotFoundError(
            f"Schmutzi contaminant database not found at '{db_path}'. "
            "In strict mode, the database must be present on disk."
        )
    logger.info(f"[Pre-flight] Contaminant DB verified: '{db_path.name}'")

    return ReferenceInfo(fasta_path=ref_path, contig_name=target_contig, length_bp=ref_length)


def main() -> int:
    sample_id = "OTZI"
    ref_fasta = REPO_ROOT / "data" / "reference" / "rCRS.fasta"
    input_bam = REPO_ROOT / "results" / "alignment" / "OTZI" / "OTZI.dedup.bam"
    output_dir = REPO_ROOT / "results" / "contamination"
    log_dir = REPO_ROOT / "logs" / "contamination"
    schmutzi_db = REPO_ROOT / "data" / "reference" / "schmutzi_db" / "197"

    logger.info("=" * 70)
    logger.info(f"STARTING STRICT BLOCK 5 RUNNER: Sample '{sample_id}'")
    logger.info("Mode: STRICT SCHMUTZI (allow_fallback_to_dissent = False)")
    logger.info("=" * 70)

    # 1. Pre-flight Validation
    try:
        reference_info = preflight_check(
            bam_path=input_bam,
            ref_path=ref_fasta,
            db_path=schmutzi_db,
        )
    except Exception as exc:
        logger.error(f"Pre-flight check failed: {exc}")
        return 1

    # 2. Strict Block 5 Policy Configuration
    schmutzi_pol = SchmutziPolicy(
        contaminant_database=schmutzi_db,
        use_wrapped_reference=True,
        padding_bp=200,
        max_iterations=10,
        threads=4,
        library_treatment=LibraryTreatment.NONE,  # Standard double-stranded library
    )

    dissent_pol = DissentPolicy(
        min_base_quality=30,
        min_mapping_quality=30,
        min_depth_per_site=10,
        noise_vaf_floor=0.015,
        base_min_dissenting_sites=3,
        dissent_high_confidence_min_sites=50,
    )

    policy = ContaminationPolicy(
        policy_id="otzi_contamination_strict_v1",
        method=ContaminationMethod.SCHMUTZI,
        enabled=True,
        allow_fallback_to_dissent=False,  # Strict mode: fail immediately if schmutzi encounters an error
        schmutzi=schmutzi_pol,
        dissent=dissent_pol,
        low_threshold=0.02,
        moderate_threshold=0.05,
        elevated_threshold=0.10,
        min_coverage_for_estimation=5.0,
        recommended_coverage=10.0,
        high_confidence_max_ci_width=0.04,
        moderate_confidence_max_ci_width=0.08,
        min_informative_sites=10,
        pmd_filtering_threshold=3.0,
    )

    # 3. Pipeline Execution
    t0 = time.perf_counter()
    result = ContaminationOrchestrator.run_sample_contamination_pipeline(
        sample_id=sample_id,
        reference_info=reference_info,
        input_bam=input_bam,
        output_dir=output_dir,
        log_dir=log_dir,
        policy=policy,
        is_unrescaled_confirmed=True,
        threads=4,
    )
    elapsed = time.perf_counter() - t0

    # 4. Strict Deliverables & JSON Contract Validation
    json_path = output_dir / sample_id / f"{sample_id}.contamination.json"
    if not json_path.is_file():
        logger.error(f"Deliverable JSON was not generated: {json_path}")
        return 1

    try:
        with open(json_path, "r", encoding="utf-8") as jf:
            json_data = json.load(jf)

        # Assert mandatory JSON schema keys
        for key in ["sample_id", "status", "estimate", "diagnostics", "gating", "outputs", "audit"]:
            if key not in json_data:
                raise ValueError(f"Missing mandatory key in contamination JSON: '{key}'")

        # Verify probabilistic range
        pt_est = json_data["estimate"]["point_estimate_fraction"]
        if not (0.0 <= pt_est <= 1.0):
            raise ValueError(f"Contamination point estimate outside [0, 1]: {pt_est}")

        # Enforce that schmutzi was actually executed (strict mode requirement)
        method_used = json_data["estimate"]["method_used"]
        if method_used != "SCHMUTZI":
            logger.error(f"Strict mode violation: expected 'SCHMUTZI', but executed method was '{method_used}'.")
            return 1

    except Exception as exc:
        logger.error(f"Deliverable JSON validation failed: {exc}")
        return 1

    # 5. Format & Report Results
    ci_display = "N/A"
    if result.estimate.ci_lower_95 is not None and result.estimate.ci_upper_95 is not None:
        width_str = f"{result.diagnostics.ci_width:.2%}" if result.diagnostics.ci_width is not None else "N/A"
        ci_display = f"[{result.estimate.ci_lower_95:.2%}, {result.estimate.ci_upper_95:.2%}] (width={width_str})"

    logger.info("=" * 70)
    logger.info(f"STRICT BLOCK 5 EXECUTION SUMMARY: Sample '{sample_id}'")
    logger.info("=" * 70)
    logger.info(f"Status               : {result.status.value}")
    logger.info(f"Method Executed      : {result.estimate.method_used.value}")
    logger.info(f"Point Estimate       : {result.estimate.point_estimate:.2%}")
    logger.info(f"95% Credible Interval: {ci_display}")
    if result.estimate.cont_deam_prior_estimate is not None:
        logger.info(f"contDeam Prior Est   : {result.estimate.cont_deam_prior_estimate:.2%}")
    logger.info(f"Contamination Cat    : {result.gating.contamination_category.value}")
    logger.info(f"Model Confidence     : {result.diagnostics.estimation_confidence.value}")
    logger.info(f"MCMC Convergence     : {result.diagnostics.model_converged}")
    logger.info(f"Iterations Run       : {result.diagnostics.iterations_run}")
    logger.info(f"Informative Sites    : {result.diagnostics.informative_positions_evaluated}")
    logger.info(f"Mean Depth           : {result.diagnostics.mean_depth:.1f}x")
    logger.info(f"Downstream Action    : {result.gating.recommended_action.value}")
    logger.info(f"Preferred Source     : {result.gating.preferred_consensus_source}")
    logger.info(f"Gating Directives    : {result.gating.summary_message}")
    if result.warning_codes:
        logger.info(f"Warning Codes        : {[w.value for w in result.warning_codes]}")
    if result.error_codes:
        logger.info(f"Error Codes          : {[e.value for e in result.error_codes]}")
    logger.info(f"Deliverable JSON     : {json_path}")
    logger.info(f"Execution Wall Time  : {elapsed:.2f}s")
    logger.info("=" * 70)

    return 0 if result.status != ContaminationStatus.FAIL else 1


if __name__ == "__main__":
    sys.exit(main())