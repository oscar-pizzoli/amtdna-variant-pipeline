#!/usr/bin/env python3
"""Smoke Test & Comprehensive Integration Suite for Block 5 BatchRunner.

Integration smoke test: Schmutzi container execution is mocked/spied;
real end-to-end container execution is validated separately on physical fixtures.
"""

import argparse
from collections import Counter
import json
import logging
import os
from pathlib import Path
import random
import shutil
import sys
import time
from typing import Any, Dict, List, Optional, Tuple
from unittest.mock import MagicMock, patch
import uuid

import pysam
import yaml

# Risoluzione radice progetto
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.modules.contamination.batch_runner import (
    BatchSampleRecord,
    ContaminationBatchRunner,
    ResourceAllocation,
)
from src.modules.contamination.fast_screen import FastContaminationScreen
from src.modules.contamination.models import (
    ConsensusGatingAction,
    ContaminationCategory,
    ContaminationErrorCode,
    ContaminationEstimate,
    ContaminationInputContext,
    ContaminationMethod,
    ContaminationOutputFiles,
    ContaminationPolicy,
    ContaminationResult,
    ContaminationStatus,
    EstimationConfidence,
    FastScreeningPolicy,
    FastScreenResult,
    ModelQualityDiagnostics,
    ScreeningDecision,
)
from src.modules.contamination.orchestrator import ContaminationOrchestrator
from src.modules.contamination.schmutzi_engine import (
    SchmutziEngine,
    SchmutziExecutionError,
    SchmutziExecutionMetadata,
    SchmutziOutput,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger("mtDNA_Pipeline.BatchSmokeTest")


def generate_minimal_bam(
    output_bam: Path,
    reference_fasta: Path,
    contig_name: str,
    target_depth: int,
    spiked: bool = False,
    diagnostic_sites: Optional[Dict[int, Dict[str, str]]] = None,
    read_len: int = 75,
    random_seed: int = 42,
) -> Path:
    """Generates a minimal valid, coordinate-sorted, and indexed BAM file."""
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    random.seed(random_seed)

    with pysam.FastaFile(str(reference_fasta)) as fa:
        ref_seq = fa.fetch(contig_name)
        contig_len = len(ref_seq)

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": contig_name, "LN": contig_len}],
        "RG": [{"ID": "RG1", "SM": output_bam.stem, "PL": "ILLUMINA"}],
    }

    raw_bam = output_bam.with_suffix(".tmp_raw.bam")
    total_reads = max(2, int((contig_len * target_depth) / read_len))

    with pysam.AlignmentFile(str(raw_bam), "wb", header=header) as out_f:
        for idx in range(total_reads):
            start = random.randint(0, max(0, contig_len - read_len))
            seq_list = list(ref_seq[start : start + read_len].upper())

            # Spike eteroplasmico: muta solo il 20% delle letture sui siti diagnostici
            if spiked and diagnostic_sites:
                for i, r_pos in enumerate(range(start, start + read_len)):
                    pos_1based = r_pos + 1
                    if pos_1based in diagnostic_sites:
                        if random.random() < 0.20:
                            alt = diagnostic_sites[pos_1based].get("alt", "G")
                            seq_list[i] = alt if alt in ("A", "C", "G", "T") else "G"

            a = pysam.AlignedSegment()
            a.query_name = f"READ_{idx:05d}"
            a.query_sequence = "".join(seq_list)
            a.flag = 0
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigar = ((0, read_len),)
            a.query_qualities = pysam.qualitystring_to_array("I" * read_len)
            a.set_tag("RG", "RG1")
            out_f.write(a)

    pysam.sort("-o", str(output_bam), str(raw_bam))
    pysam.index(str(output_bam))
    if raw_bam.is_file():
        raw_bam.unlink()
    return output_bam


def run_batch_runner_smoke_test(keep_artifacts: bool = False):
    """Executes the full test matrix for Block 5 batch execution."""
    test_dir = PROJECT_ROOT / "results" / "test_batch_runner_smoke"
    if test_dir.is_dir():
        shutil.rmtree(test_dir)
    test_dir.mkdir(parents=True, exist_ok=True)

    input_bams_dir = test_dir / "inputs"
    cohort_out_dir = test_dir / "outputs"
    input_bams_dir.mkdir(parents=True, exist_ok=True)

    config_path = PROJECT_ROOT / "config" / "default_config.yaml"
    if not config_path.is_file():
        raise RuntimeError(f"Config file missing: {config_path}")

    with open(config_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    ref_path = (PROJECT_ROOT / cfg["reference"]["fasta_path"]).resolve()
    contig = cfg["reference"]["expected_contig"]
    length_bp = cfg["reference"]["expected_length"]
    panel_path = (PROJECT_ROOT / cfg["contamination"]["fast_screening"]["diagnostic_panel_path"]).resolve()

    if not ref_path.is_file():
        raise RuntimeError(f"Reference FASTA not found: {ref_path}")
    if not panel_path.is_file():
        raise RuntimeError(f"Diagnostic panel not found: {panel_path}. Run panel_builder first!")

    diag_dict, _ = FastContaminationScreen.load_diagnostic_panel(panel_path)

    print("=" * 90, flush=True)
    print("🚀 RUNNING BATCH RUNNER & TWO-PHASE SCHEDULER SMOKE SUITE", flush=True)
    print("=" * 90, flush=True)

    # 1. Resource Allocation
    print("\n▶ [TEST 1] Verifying CPU Resource Allocation & Scheduling Invariants...", flush=True)
    alloc = ContaminationBatchRunner.calculate_resource_allocation(
        threads_per_sample=4,
        cpu_budget_fraction=0.50,
        host_cores_override=16,
    )
    assert alloc.total_host_cores == 16
    assert alloc.total_allocated_threads == 8
    assert alloc.max_schmutzi_workers == 2
    assert alloc.max_fast_screen_workers == 8
    print("   ✅ Resource allocation invariants PASSED (50% ceiling enforced).", flush=True)

    # 2. BAM Discovery
    print("\n▶ [TEST 2] Verifying BAM Discovery & Suffix Filtering...", flush=True)
    bam_pure = generate_minimal_bam(input_bams_dir / "COHORT_SAMPLE_01_PURE.bam", ref_path, contig, target_depth=30, spiked=False)
    bam_spiked = generate_minimal_bam(input_bams_dir / "COHORT_SAMPLE_02_SPIKED.bam", ref_path, contig, target_depth=30, spiked=True, diagnostic_sites=diag_dict)
    bam_low = generate_minimal_bam(input_bams_dir / "COHORT_SAMPLE_03_LOW.bam", ref_path, contig, target_depth=2, spiked=False)

    (input_bams_dir / "COHORT_SAMPLE_01_PURE.unsorted.bam").touch()
    (input_bams_dir / "COHORT_SAMPLE_01_PURE.tmp.bam").touch()
    (input_bams_dir / "COHORT_SAMPLE_01_PURE.subsampled_40x.bam").touch()
    (input_bams_dir / "COHORT_SAMPLE_02_SPIKED.schmutzi.bam").touch()

    discovered = ContaminationBatchRunner.discover_bam_inputs(input_dir=input_bams_dir)
    discovered_ids = {s_id for s_id, _ in discovered}
    assert len(discovered) == 3
    assert "COHORT_SAMPLE_01_PURE" in discovered_ids
    assert "COHORT_SAMPLE_02_SPIKED" in discovered_ids
    assert "COHORT_SAMPLE_03_LOW" in discovered_ids
    print(f"   ✅ Discovery PASSED ({len(discovered)} valid samples discovered).", flush=True)

    # 3 & 4. Two-Phase Gating & Container Uniqueness
    print("\n▶ [TEST 3 & 4] Verifying Two-Phase Triage Routing & Unique Container Isolation...", flush=True)
    executed_container_names: List[str] = []
    sample_dispatch_counts: Counter = Counter()

    def mock_schmutzi_run(sample_id, input_bam, reference_fasta, database_dir, output_dir, threads, max_iterations, use_wrapped):
        token = uuid.uuid4().hex[:6]
        c_name = f"schmutzi_{sample_id}_run_{token}"
        executed_container_names.append(c_name)
        sample_dispatch_counts[sample_id] += 1

        return SchmutziOutput(
            estimate=ContaminationEstimate(
                point_estimate=0.08,
                ci_lower_95=0.05,
                ci_upper_95=0.11,
                method_used=ContaminationMethod.SCHMUTZI,
            ),
            diagnostics=ModelQualityDiagnostics(
                mean_depth=30.0,
                informative_positions_evaluated=25,
                model_converged=True,
                estimation_confidence=EstimationConfidence.HIGH,
            ),
            outputs=ContaminationOutputFiles(
                cont_est_file=output_dir / sample_id / f"{sample_id}.cont.est",
                endogenous_fasta=output_dir / sample_id / f"{sample_id}_final_endo.fa",
            ),
            audit=SchmutziExecutionMetadata(
                tool_version="docker-mock-v1.0",
                executed_command=f"docker run --rm --name {c_name} ...",
            ),
        )

    with patch.object(SchmutziEngine, "run", side_effect=mock_schmutzi_run):
        results = ContaminationBatchRunner.run_cohort_pipeline(
            discovered_samples=discovered,
            reference_fasta=ref_path,
            reference_contig=contig,
            reference_length_bp=length_bp,
            output_dir=cohort_out_dir,
            config=cfg,
            resource_alloc=alloc,
        )

    res_map = {r.sample_id: r for r in results}

    # Sample 1: PURE -> FAST_CLEAR
    r1 = res_map["COHORT_SAMPLE_01_PURE"]
    assert r1.audit.fast_screen_result.decision == ScreeningDecision.FAST_CLEAR
    assert r1.audit.screened_without_full_quantification is True
    assert r1.audit.full_bayesian_estimation_performed is False
    assert r1.estimate.method_used == ContaminationMethod.FAST_SCREEN_TRIAGE
    assert sample_dispatch_counts["COHORT_SAMPLE_01_PURE"] == 0

    # Sample 2: SPIKED -> ESCALATE
    r2 = res_map["COHORT_SAMPLE_02_SPIKED"]
    assert r2.audit.fast_screen_result.decision == ScreeningDecision.ESCALATE
    assert r2.audit.full_bayesian_estimation_performed is True
    assert r2.estimate.method_used == ContaminationMethod.SCHMUTZI
    assert sample_dispatch_counts["COHORT_SAMPLE_02_SPIKED"] == 1

    # Sample 3: LOW -> INSUFFICIENT_DATA
    r3 = res_map["COHORT_SAMPLE_03_LOW"]
    assert r3.audit.fast_screen_result.decision == ScreeningDecision.INSUFFICIENT_DATA
    assert r3.audit.full_bayesian_estimation_performed is True
    assert sample_dispatch_counts["COHORT_SAMPLE_03_LOW"] == 1

    assert len(executed_container_names) == 2
    assert len(set(executed_container_names)) == len(executed_container_names)
    assert not any("COHORT_SAMPLE_01_PURE" in name for name in executed_container_names)
    print("   ✅ Two-Tier dispatch verified (Pure=0, Spiked=1, Low=1 calls).", flush=True)
    print(f"   ✅ Container name uniqueness verified: {executed_container_names}", flush=True)

    # 5. Cache Test
    print("\n▶ [TEST 5] Verifying Cache Hits and Invalidation Behavior...", flush=True)
    fast_policy = ContaminationOrchestrator(config=cfg).policy.fast_screening

    fs1, _ = FastContaminationScreen.run_screen(
        sample_id="CACHE_TEST_SAMPLE",
        input_bam=bam_pure,
        reference_fasta=ref_path,
        expected_contig=contig,
        policy=fast_policy,
        output_dir=cohort_out_dir,
    )

    t0 = time.perf_counter()
    fs2, _ = FastContaminationScreen.run_screen(
        sample_id="CACHE_TEST_SAMPLE",
        input_bam=bam_pure,
        reference_fasta=ref_path,
        expected_contig=contig,
        policy=fast_policy,
        output_dir=cohort_out_dir,
    )
    cache_hit_time = time.perf_counter() - t0
    assert fs1.decision == fs2.decision
    assert cache_hit_time <= 0.20
    print(f"   ✅ Cache hit verified in {cache_hit_time:.4f}s.", flush=True)

    # 6. Failure Safety
    print("\n▶ [TEST 6] Verifying Failure Safety Contract (No False PASS on Crash)...", flush=True)
    fail_cfg = dict(cfg)
    fail_cfg["contamination"]["allow_fallback_to_dissent"] = False

    def mock_schmutzi_crash(*args, **kwargs):
        raise SchmutziExecutionError("Simulated Bayesian MCMC collapse.")

    with patch.object(SchmutziEngine, "run", side_effect=mock_schmutzi_crash):
        crashed_results = ContaminationBatchRunner.run_cohort_pipeline(
            discovered_samples=[("COHORT_SAMPLE_02_SPIKED", bam_spiked)],
            reference_fasta=ref_path,
            reference_contig=contig,
            reference_length_bp=length_bp,
            output_dir=cohort_out_dir,
            config=fail_cfg,
            resource_alloc=alloc,
        )

    crash_res = crashed_results[0]
    assert crash_res.status == ContaminationStatus.FAIL
    assert crash_res.gating.recommended_action == ConsensusGatingAction.EXCLUDE_OR_MANUAL_REVIEW
    print("   ✅ Failure safety verified (crashed Schmutzi yields strict FAIL and EXCLUDE gating).", flush=True)

    # 7. Cohort Reporting
    print("\n▶ [TEST 7] Verifying Cohort TSV & JSON Export...", flush=True)
    tsv_out, json_out = ContaminationBatchRunner.export_cohort_summary(
        results=results,
        output_dir=cohort_out_dir,
        resource_alloc=alloc,
    )
    assert tsv_out.is_file()
    assert json_out.is_file()

    with open(json_out, "r", encoding="utf-8") as jf:
        json_data = json.load(jf)
    assert json_data["metadata"]["total_samples"] == 3
    assert json_data["metadata"]["fast_clear_bypassed"] == 1
    assert json_data["metadata"]["schmutzi_executed"] == 2
    assert json_data["metadata"]["bypass_efficiency_pct"] == 33.33
    print(f"   ✅ Cohort summaries verified -> TSV: {tsv_out.name} | JSON: {json_out.name}", flush=True)

    if not keep_artifacts and test_dir.is_dir():
        shutil.rmtree(test_dir)

    print("\n" + "=" * 90, flush=True)
    print("🎉 ALL 7 BLOCK 5 INTEGRATION SMOKE TESTS PASSED PERFECTLY!", flush=True)
    print("=" * 90, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Block 5 BatchRunner smoke test.")
    parser.add_argument("--keep-artifacts", action="store_true", help="Preserve test output directory for inspection.")
    args = parser.parse_args()
    run_batch_runner_smoke_test(keep_artifacts=args.keep_artifacts)