#!/usr/bin/env python3
"""Benchmark and Empirical Validation Suite for mtDNA Fast-Gating Triage.

Validates the screening and QC gating contracts across 6 test scenarios:
  1. Pure Sample (Scenario A) -> FAST_CLEAR (Bypass Eligible, under runtime budget)
  2. Borderline Sub-Threshold (Scenario B) -> FAST_CLEAR (1 site < 3)
  3. Spiked Dissent Sample (Scenario C) -> ESCALATE (7.68% spike on diagnostic sites)
  4. Low-Depth Sample (Scenario D) -> INSUFFICIENT_DATA (Safety escalation on <5x depth)
  5. Deterministic Cache Hit (Scenario E) -> Verified cache retrieval & strict runtime budget
  6. Algorithmic Determinism (Scenario F) -> Bitwise identical result across independent runs

Also validates the two-tier routing contract (0 downstream invocations on FAST_CLEAR,
1 downstream invocation on ESCALATE and INSUFFICIENT_DATA).
"""

import argparse
from dataclasses import dataclass
import json
import logging
from pathlib import Path
import random
import shutil
import sys
import time
from typing import Any, Dict, List, Optional, Tuple

import pysam
import yaml

# Risoluzione dinamica della radice del progetto
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.modules.contamination.fast_screen import FastContaminationScreen
from src.modules.contamination.models import (
    ConsensusGatingAction,
    ContaminationAuditMetadata,
    ContaminationCategory,
    ContaminationEstimate,
    ContaminationMethod,
    ContaminationOutputFiles,
    ContaminationPolicy,
    FastScreeningPolicy,
    FastScreenResult,
    ModelQualityDiagnostics,
    ScreeningDecision,
)
from src.modules.contamination.qc_evaluator import ContaminationQCValidator

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger("mtDNA_Pipeline.Benchmark")


@dataclass
class PerformanceBudget:
    """Configurable performance budgets for benchmark assertions."""

    max_fast_clear_runtime_seconds: float = 5.0
    max_cache_hit_runtime_seconds: float = 0.20


class TwoTierRoutingSpy:
    """Simulates and verifies the two-tier orchestration dispatch contract."""

    def __init__(self) -> None:
        self.escalation_call_count: int = 0

    def dispatch(self, decision: ScreeningDecision) -> bool:
        """Evaluates whether the decision triggers downstream escalation."""
        if decision in (ScreeningDecision.ESCALATE, ScreeningDecision.INSUFFICIENT_DATA):
            self.escalation_call_count += 1
            return True
        return False


class FastGateBenchmarkRunner:
    """Orchestrates empirical validation of FastContaminationScreen and QCValidator."""

    @classmethod
    def load_yaml_config(cls, config_path: Path) -> Dict[str, Any]:
        """Loads configuration from YAML with explicit error handling."""
        if not config_path.is_file():
            raise FileNotFoundError(f"Configuration file not found: {config_path}")
        try:
            with open(config_path, "r", encoding="utf-8") as f:
                return yaml.safe_load(f) or {}
        except Exception as exc:
            raise RuntimeError(f"Failed loading YAML configuration from {config_path}: {exc}") from exc

    @classmethod
    def resolve_path(cls, path_str: Optional[str], default_relative: Path) -> Path:
        """Resolves path string relative to PROJECT_ROOT."""
        if not path_str:
            return default_relative if default_relative.is_absolute() else (PROJECT_ROOT / default_relative)
        p = Path(path_str)
        return p if p.is_absolute() else (PROJECT_ROOT / p)

    @classmethod
    def measure_actual_mean_depth(cls, bam_path: Path, contig_name: str, contig_length: int) -> float:
        """Measures true average sequencing depth from generated BAM."""
        if not bam_path.is_file() or contig_length <= 0:
            return 0.0
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            if contig_name not in bam.references:
                return 0.0
            total_bases = sum(col.nsegments for col in bam.pileup(contig=contig_name, truncate=True))
            return float(total_bases / contig_length)

    @classmethod
    def generate_synthetic_bam(
        cls,
        output_bam: Path,
        reference_fasta: Path,
        contig_name: str,
        target_depth: int,
        spiked_error_rate: float = 0.0,
        specific_spiked_sites: Optional[Dict[int, str]] = None,
        diagnostic_sites: Optional[Dict[int, Dict[str, str]]] = None,
        read_length: int = 75,
        random_seed: int = 42,
    ) -> Tuple[Path, float]:
        """Generates a coordinate-sorted and indexed BAM, returning path and actual depth."""
        output_bam.parent.mkdir(parents=True, exist_ok=True)
        random.seed(random_seed)

        with pysam.FastaFile(str(reference_fasta)) as fa:
            if contig_name not in fa.references:
                raise RuntimeError(f"Contig '{contig_name}' missing from reference FASTA.")
            ref_seq = fa.fetch(contig_name)
            contig_len = len(ref_seq)

        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": contig_name, "LN": contig_len}],
            "RG": [{"ID": "BENCH_RG", "SM": output_bam.stem, "PL": "ILLUMINA"}],
        }

        raw_bam_path = output_bam.with_suffix(".unsorted.bam")
        total_reads = max(1, int((contig_len * target_depth) / read_length))

        with pysam.AlignmentFile(str(raw_bam_path), "wb", header=header) as out_bam:
            for read_idx in range(total_reads):
                start_pos = random.randint(0, max(0, contig_len - read_length))
                end_pos = start_pos + read_length
                seq_bases = list(ref_seq[start_pos:end_pos].upper())

                for i, r_pos_0based in enumerate(range(start_pos, end_pos)):
                    pos_1based = r_pos_0based + 1

                    if specific_spiked_sites and pos_1based in specific_spiked_sites:
                        seq_bases[i] = specific_spiked_sites[pos_1based]

                    elif spiked_error_rate > 0.0 and diagnostic_sites and pos_1based in diagnostic_sites:
                        if random.random() < spiked_error_rate:
                            expected_alt = diagnostic_sites[pos_1based].get("alt", "")
                            if expected_alt and expected_alt in ("A", "C", "G", "T"):
                                seq_bases[i] = expected_alt
                            else:
                                alt_pool = [b for b in ("A", "C", "G", "T") if b != seq_bases[i]]
                                seq_bases[i] = random.choice(alt_pool)

                read_seq = "".join(seq_bases)
                qual_str = "".join([chr(30 + 33)] * read_length)

                a = pysam.AlignedSegment()
                a.query_name = f"SYNTH_{read_idx:07d}"
                a.query_sequence = read_seq
                a.flag = 0
                a.reference_id = 0
                a.reference_start = start_pos
                a.mapping_quality = 60
                a.cigar = ((0, read_length),)
                a.query_qualities = pysam.qualitystring_to_array(qual_str)
                a.set_tag("RG", "BENCH_RG")
                out_bam.write(a)

        pysam.sort("-o", str(output_bam), str(raw_bam_path))
        pysam.index(str(output_bam))

        if raw_bam_path.is_file():
            raw_bam_path.unlink()

        actual_depth = cls.measure_actual_mean_depth(output_bam, contig_name, contig_len)
        return output_bam, actual_depth


def run_benchmark() -> None:
    """CLI and orchestration entrypoint for the complete benchmark suite."""
    parser = argparse.ArgumentParser(description="Empirical benchmark for mtDNA Fast-Gating triage.")
    parser.add_argument(
        "--config",
        type=Path,
        default=PROJECT_ROOT / "config" / "default_config.yaml",
        help="Path to YAML configuration file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "benchmark_fastgate",
        help="Directory to store benchmark outputs and temporary BAMs.",
    )
    parser.add_argument(
        "--max-runtime",
        type=float,
        default=5.0,
        help="Maximum allowed runtime in seconds for FAST_CLEAR execution (budget: 5.0s).",
    )
    parser.add_argument(
        "--max-cache-runtime",
        type=float,
        default=0.20,
        help="Maximum allowed runtime in seconds for cache hit retrieval (budget: 0.20s).",
    )
    parser.add_argument(
        "--keep-bams",
        action="store_true",
        help="Keep generated synthetic BAMs after benchmark completion.",
    )

    args = parser.parse_args()
    budget = PerformanceBudget(
        max_fast_clear_runtime_seconds=args.max_runtime,
        max_cache_hit_runtime_seconds=args.max_cache_runtime,
    )

    out_dir = args.output_dir
    # Reset della cartella del benchmark per eliminare cache obsolete
    if out_dir.is_dir():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 90)
    print("🔬 STARTING mtDNA FAST-GATING TRIAGE BENCHMARK SUITE")
    print("=" * 90)

    # 1. Caricamento configurazione e policy
    cfg = FastGateBenchmarkRunner.load_yaml_config(args.config)
    ref_cfg = cfg.get("reference", {})
    contam_cfg = cfg.get("contamination", {})
    fast_cfg = contam_cfg.get("fast_screening", {})

    ref_fasta = FastGateBenchmarkRunner.resolve_path(
        ref_cfg.get("fasta_path"),
        PROJECT_ROOT / "data" / "reference" / "rCRS.fasta",
    )
    contig_name = ref_cfg.get("expected_contig", "NC_012920.1")

    panel_path = FastGateBenchmarkRunner.resolve_path(
        fast_cfg.get("diagnostic_panel_path"),
        PROJECT_ROOT / "data" / "reference" / "diagnostic_sites_phylotree17.tsv",
    )

    if not ref_fasta.is_file():
        raise RuntimeError(f"Reference FASTA not found at {ref_fasta}")
    if not panel_path.is_file():
        raise RuntimeError(f"Diagnostic panel not found at {panel_path}. Run panel_builder.py first!")

    fast_policy = FastScreeningPolicy(
        enabled=True,
        algorithm_version=fast_cfg.get("algorithm_version", "fast-screen-v2.0"),
        diagnostic_panel_path=panel_path,
        diagnostic_panel_version=fast_cfg.get("diagnostic_panel_version", "phylotree17_v1"),
        escalate_dissenting_sites_min=fast_cfg.get("escalate_dissenting_sites_min", 3),
        escalate_vaf_threshold=fast_cfg.get("escalate_vaf_threshold", 0.02),
        min_alt_reads_per_site=fast_cfg.get("min_alt_reads_per_site", 2),
        min_depth_for_screen=fast_cfg.get("min_depth_for_screen", 5.0),
        min_evaluated_sites=fast_cfg.get("min_evaluated_sites", 10),
        min_base_quality=fast_cfg.get("min_base_quality", 30),
        min_mapping_quality=fast_cfg.get("min_mapping_quality", 30),
        min_depth_per_site=fast_cfg.get("min_depth_per_site", 10),
        noise_vaf_floor=fast_cfg.get("noise_vaf_floor", 0.015),
        cache_screen_results=True,
        cache_mode="strict",
    )

    diag_dict, _ = FastContaminationScreen.load_diagnostic_panel(panel_path)
    print(f"[*] Reference: {ref_fasta.name} ({contig_name})")
    print(f"[*] Panel: {panel_path.name} ({len(diag_dict)} sites)")
    print(f"[*] Performance Budget: FAST_CLEAR <= {budget.max_fast_clear_runtime_seconds:.1f}s | Cache <= {budget.max_cache_hit_runtime_seconds:.2f}s")
    print("-" * 90)

    # 2. Generazione dataset sintetici
    synthetic_bams_dir = out_dir / "synthetic_inputs"
    synthetic_bams_dir.mkdir(parents=True, exist_ok=True)

    print("[*] Generating Synthetic Validation Datasets...")

    # Dataset A: Puro (~40x)
    bam_pure, depth_pure = FastGateBenchmarkRunner.generate_synthetic_bam(
        output_bam=synthetic_bams_dir / "SYNTH_PURE_40X.bam",
        reference_fasta=ref_fasta,
        contig_name=contig_name,
        target_depth=40,
        spiked_error_rate=0.0,
        diagnostic_sites=diag_dict,
    )

    # Dataset B: Borderline sub-threshold (~40x, esattamente 1 sito divergente con alt >= 2 letture)
    borderline_sites = {}
    if diag_dict:
        first_diag_pos = sorted(diag_dict.keys())[0]
        alt_base = "G" if diag_dict[first_diag_pos].get("ref") != "G" else "A"
        borderline_sites[first_diag_pos] = alt_base

    bam_borderline, depth_borderline = FastGateBenchmarkRunner.generate_synthetic_bam(
        output_bam=synthetic_bams_dir / "SYNTH_BORDERLINE_1SITE_40X.bam",
        reference_fasta=ref_fasta,
        contig_name=contig_name,
        target_depth=40,
        specific_spiked_sites=borderline_sites,
        diagnostic_sites=diag_dict,
    )

    # Dataset C: Spiked con dissenso allelico diffuso (~40x, 7.68% spike probabilistico su siti diagnostici)
    bam_spiked, depth_spiked = FastGateBenchmarkRunner.generate_synthetic_bam(
        output_bam=synthetic_bams_dir / "SYNTH_SPIKED_7.68pct_40X.bam",
        reference_fasta=ref_fasta,
        contig_name=contig_name,
        target_depth=40,
        spiked_error_rate=0.0768,
        diagnostic_sites=diag_dict,
    )

    # Dataset D: Bassa profondità (~2x)
    bam_low, depth_low = FastGateBenchmarkRunner.generate_synthetic_bam(
        output_bam=synthetic_bams_dir / "SYNTH_LOW_DEPTH_2X.bam",
        reference_fasta=ref_fasta,
        contig_name=contig_name,
        target_depth=2,
        spiked_error_rate=0.0,
        diagnostic_sites=diag_dict,
    )

    print(f"    - Pure Sample: requested 40x -> actual {depth_pure:.1f}x")
    print(f"    - Borderline Sample: requested 40x -> actual {depth_borderline:.1f}x")
    print(f"    - Spiked Sample: requested 40x -> actual {depth_spiked:.1f}x")
    print(f"    - Low-Depth Sample: requested 2x -> actual {depth_low:.1f}x")
    print("✅ Synthetic BAM datasets prepared successfully.")
    print("-" * 90)

    # 3. Esecuzione Test Suite
    results_summary: List[Dict[str, Any]] = []
    routing_spy = TwoTierRoutingSpy()

    # TEST 1: Pure Sample
    print("\n▶ [TEST 1] Pure Sample (~40x) -> Expect: FAST_CLEAR & 0 Escalations")
    t0 = time.perf_counter()
    fs_pure, _ = FastContaminationScreen.run_screen(
        sample_id="SYNTH_PURE_40X",
        input_bam=bam_pure,
        reference_fasta=ref_fasta,
        expected_contig=contig_name,
        policy=fast_policy,
        output_dir=out_dir,
    )
    t1_pure = time.perf_counter() - t0
    esc_pure = routing_spy.dispatch(fs_pure.decision)

    pass_1 = (
        fs_pure.decision == ScreeningDecision.FAST_CLEAR
        and fs_pure.bypass_eligible is True
        and not esc_pure
        and t1_pure <= budget.max_fast_clear_runtime_seconds
    )
    results_summary.append({
        "Test": "1. Pure (40x)",
        "Runtime": f"{t1_pure:.3f}s",
        "Decision": fs_pure.decision.value,
        "Expected": ScreeningDecision.FAST_CLEAR.value,
        "Bypass": str(fs_pure.bypass_eligible),
        "Depth": f"{fs_pure.mean_diagnostic_depth:.1f}x",
        "Dissent": fs_pure.dissenting_positions_count,
        "Max VAF": f"{fs_pure.max_dissent_vaf:.2%}",
        "Escalations": "0",
        "Result": "PASS" if pass_1 else "FAIL",
    })

    # TEST 2: Borderline Sub-Threshold Sample
    print("\n▶ [TEST 2] Borderline Sub-Threshold (1 dissenting site < 3) -> Expect: FAST_CLEAR & 0 Escalations")
    t0 = time.perf_counter()
    fs_bord, _ = FastContaminationScreen.run_screen(
        sample_id="SYNTH_BORDERLINE_1SITE_40X",
        input_bam=bam_borderline,
        reference_fasta=ref_fasta,
        expected_contig=contig_name,
        policy=fast_policy,
        output_dir=out_dir,
    )
    t1_bord = time.perf_counter() - t0
    esc_bord = routing_spy.dispatch(fs_bord.decision)

    pass_2 = (
        fs_bord.decision == ScreeningDecision.FAST_CLEAR
        and fs_bord.bypass_eligible is True
        and fs_bord.dissenting_positions_count < fast_policy.escalate_dissenting_sites_min
        and not esc_bord
    )
    results_summary.append({
        "Test": "2. Borderline (1 Site)",
        "Runtime": f"{t1_bord:.3f}s",
        "Decision": fs_bord.decision.value,
        "Expected": ScreeningDecision.FAST_CLEAR.value,
        "Bypass": str(fs_bord.bypass_eligible),
        "Depth": f"{fs_bord.mean_diagnostic_depth:.1f}x",
        "Dissent": fs_bord.dissenting_positions_count,
        "Max VAF": f"{fs_bord.max_dissent_vaf:.2%}",
        "Escalations": "0",
        "Result": "PASS" if pass_2 else "FAIL",
    })

    # TEST 3: Spiked Dissent Sample
    print("\n▶ [TEST 3] Spiked Dissent Sample (7.68% spike) -> Expect: ESCALATE & 1 Escalation")
    t0 = time.perf_counter()
    fs_spiked, _ = FastContaminationScreen.run_screen(
        sample_id="SYNTH_SPIKED_7.68pct_40X",
        input_bam=bam_spiked,
        reference_fasta=ref_fasta,
        expected_contig=contig_name,
        policy=fast_policy,
        output_dir=out_dir,
    )
    t1_spiked = time.perf_counter() - t0
    esc_spiked = routing_spy.dispatch(fs_spiked.decision)

    pass_3 = (
        fs_spiked.decision == ScreeningDecision.ESCALATE
        and fs_spiked.bypass_eligible is False
        and esc_spiked is True
        and fs_spiked.dissenting_positions_count >= fast_policy.escalate_dissenting_sites_min
    )
    results_summary.append({
        "Test": "3. Spiked (7.68%)",
        "Runtime": f"{t1_spiked:.3f}s",
        "Decision": fs_spiked.decision.value,
        "Expected": ScreeningDecision.ESCALATE.value,
        "Bypass": str(fs_spiked.bypass_eligible),
        "Depth": f"{fs_spiked.mean_diagnostic_depth:.1f}x",
        "Dissent": fs_spiked.dissenting_positions_count,
        "Max VAF": f"{fs_spiked.max_dissent_vaf:.2%}",
        "Escalations": "1",
        "Result": "PASS" if pass_3 else "FAIL",
    })

    # TEST 4: Low-Depth Sample
    print("\n▶ [TEST 4] Low-Depth Sample (~2x < 5x) -> Expect: INSUFFICIENT_DATA & 1 Escalation")
    t0 = time.perf_counter()
    fs_low, _ = FastContaminationScreen.run_screen(
        sample_id="SYNTH_LOW_DEPTH_2X",
        input_bam=bam_low,
        reference_fasta=ref_fasta,
        expected_contig=contig_name,
        policy=fast_policy,
        output_dir=out_dir,
    )
    t1_low = time.perf_counter() - t0
    esc_low = routing_spy.dispatch(fs_low.decision)

    pass_4 = (
        fs_low.decision == ScreeningDecision.INSUFFICIENT_DATA
        and fs_low.bypass_eligible is False
        and esc_low is True
        and fs_low.mean_diagnostic_depth < fast_policy.min_depth_for_screen
    )
    results_summary.append({
        "Test": "4. Low Depth (2x)",
        "Runtime": f"{t1_low:.3f}s",
        "Decision": fs_low.decision.value,
        "Expected": ScreeningDecision.INSUFFICIENT_DATA.value,
        "Bypass": str(fs_low.bypass_eligible),
        "Depth": f"{fs_low.mean_diagnostic_depth:.1f}x",
        "Dissent": fs_low.dissenting_positions_count,
        "Max VAF": f"{fs_low.max_dissent_vaf:.2%}",
        "Escalations": "1",
        "Result": "PASS" if pass_4 else "FAIL",
    })

    # TEST 5: Cache Verification
    print(f"\n▶ [TEST 5] Deterministic Cache Retrieval -> Expect: FAST_CLEAR & Runtime <= {budget.max_cache_hit_runtime_seconds:.2f}s")
    t0 = time.perf_counter()
    fs_cache, _ = FastContaminationScreen.run_screen(
        sample_id="SYNTH_PURE_40X",
        input_bam=bam_pure,
        reference_fasta=ref_fasta,
        expected_contig=contig_name,
        policy=fast_policy,
        output_dir=out_dir,
    )
    t1_cache = time.perf_counter() - t0

    pass_5 = (
        fs_cache.decision == ScreeningDecision.FAST_CLEAR
        and fs_cache.bypass_eligible is True
        and t1_cache <= budget.max_cache_hit_runtime_seconds
    )
    results_summary.append({
        "Test": "5. Cache Hit",
        "Runtime": f"{t1_cache:.4f}s",
        "Decision": fs_cache.decision.value,
        "Expected": ScreeningDecision.FAST_CLEAR.value,
        "Bypass": str(fs_cache.bypass_eligible),
        "Depth": f"{fs_cache.mean_diagnostic_depth:.1f}x",
        "Dissent": fs_cache.dissenting_positions_count,
        "Max VAF": f"{fs_cache.max_dissent_vaf:.2%}",
        "Escalations": "0",
        "Result": "PASS" if pass_5 else "FAIL",
    })

    # TEST 6: Determinism Verification
    print("\n▶ [TEST 6] Algorithmic Determinism -> Expect: Identical result dictionary on repeat run")
    no_cache_policy = FastScreeningPolicy(
        enabled=True,
        algorithm_version=fast_policy.algorithm_version,
        diagnostic_panel_path=fast_policy.diagnostic_panel_path,
        diagnostic_panel_version=fast_policy.diagnostic_panel_version,
        cache_screen_results=False,
    )
    fs_det_1, _ = FastContaminationScreen.run_screen(
        sample_id="SYNTH_SPIKED_DET_1",
        input_bam=bam_spiked,
        reference_fasta=ref_fasta,
        expected_contig=contig_name,
        policy=no_cache_policy,
        output_dir=out_dir,
    )
    fs_det_2, _ = FastContaminationScreen.run_screen(
        sample_id="SYNTH_SPIKED_DET_2",
        input_bam=bam_spiked,
        reference_fasta=ref_fasta,
        expected_contig=contig_name,
        policy=no_cache_policy,
        output_dir=out_dir,
    )

    d1 = fs_det_1.to_dict()
    d2 = fs_det_2.to_dict()
    for k in ("sample_id", "runtime_seconds"):
        d1.pop(k, None)
        d2.pop(k, None)

    pass_6 = (d1 == d2)
    results_summary.append({
        "Test": "6. Determinism",
        "Runtime": f"{fs_det_2.runtime_seconds:.3f}s",
        "Decision": fs_det_2.decision.value,
        "Expected": fs_det_1.decision.value,
        "Bypass": str(fs_det_2.bypass_eligible),
        "Depth": f"{fs_det_2.mean_diagnostic_depth:.1f}x",
        "Dissent": fs_det_2.dissenting_positions_count,
        "Max VAF": f"{fs_det_2.max_dissent_vaf:.2%}",
        "Escalations": "N/A",
        "Result": "PASS" if pass_6 else "FAIL",
    })

    # 4. Tabella Riassuntiva Finale
    print("\n" + "=" * 90)
    print("📊 BENCHMARK EMPIRICAL SUMMARY REPORT")
    print("=" * 90)

    header = f"{'Test Case':<22} | {'Runtime':<9} | {'Decision':<18} | {'Depth':<7} | {'Dissent':<8} | {'Max VAF':<8} | {'Escalate':<8} | {'Status'}"
    print(header)
    print("-" * len(header))
    for r in results_summary:
        row = (
            f"{r['Test']:<22} | "
            f"{r['Runtime']:<9} | "
            f"{r['Decision']:<18} | "
            f"{r['Depth']:<7} | "
            f"{str(r['Dissent']):<8} | "
            f"{r['Max VAF']:<8} | "
            f"{r['Escalations']:<8} | "
            f"{r['Result']}"
        )
        print(row)
    print("=" * len(header))

    if not args.keep_bams and synthetic_bams_dir.is_dir():
        shutil.rmtree(synthetic_bams_dir)
        print("[*] Cleaned up temporary synthetic BAMs.")

    all_passed = all(r["Result"] == "PASS" for r in results_summary)
    if all_passed:
        print("\n🎉 ALL 6 BENCHMARK SCENARIOS & PERFORMANCE BUDGETS PASSED!")
    else:
        print("\n❌ BENCHMARK FAILED: One or more assertions did not pass.")
        sys.exit(1)


if __name__ == "__main__":
    run_benchmark()