#!/usr/bin/env python3
"""Test Suite for Block 4 (aDNA Damage Authentication & Quality Rescaling)."""

import array
from dataclasses import dataclass
import json
from pathlib import Path
import sys
import tempfile
import time
from typing import Callable, List, Tuple
import pysam

# Project Root Setup
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.modules.alignment.models import ReferenceInfo
from src.modules.damage.damage_quantifier import DamageQuantifierEngine
from src.modules.damage.models import (
    AuthenticationVerdict,
    BayesianDamageParameters,
    DamageErrorCode,
    DamageInputBAMType,
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
from src.modules.damage.orchestrator import DamageOrchestrator
from src.modules.damage.pmd_engine import PMDScoringEngine
from src.modules.damage.qc_evaluator import DamageQCValidator


def test_models_and_serialization() -> None:
    """Tests instantiation and JSON serialization of all data models and dataclasses."""
    policy = DamagePolicy(
        input_bam_type=DamageInputBAMType.FILTERED_DEDUP,
        library_treatment=LibraryTreatment.UDG_PARTIAL,
        min_c_to_t_5p_strong_threshold=0.10,
    )
    
    terminal = TerminalDamageFrequencies(
        c_to_t_5p_pos1=0.154,
        g_to_a_3p_pos1=0.142,
        n_c_to_t_5p_pos1_obs=77,
        n_c_ref_5p_pos1_total=500,
        c_to_t_5p_profile=[0.154, 0.082, 0.021],
    )
    
    profile = DamageModelProfile(
        sample_id="SAMPLE_TEST_01",
        library_treatment=policy.library_treatment,
        terminal_damage=terminal,
        damage_signal_level=DamageSignalLevel.STRONG_SIGNAL,
        authentication_verdict=AuthenticationVerdict.AUTHENTIC_ANCIENT_SUPPORTED,
    )
    
    res = DamageResult(
        sample_id="SAMPLE_TEST_01",
        status=DamageStatus.PASS,
        policy=policy,
        profile=profile,
    )
    
    data_dict = res.to_dict()
    json_str = json.dumps(data_dict, indent=2)
    assert "SAMPLE_TEST_01" in json_str
    assert data_dict["policy"]["library_treatment"] == "UDG_PARTIAL"
    assert data_dict["profile"]["terminal_damage"]["n_c_ref_5p_pos1_total"] == 500


def test_mapdamage_matrix_parser() -> None:
    """Tests misincorporation, length distribution, and MCMC parser on synthetic tables."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        
        misinc_file = tmp_path / "misincorporation.txt"
        with open(misinc_file, "w", encoding="utf-8") as f:
            f.write("# Chr\tEnd\tStd\tPos\tBase\tA\tC\tG\tT\tTotal\tG>A\tC>T\n")
            f.write("NC_012920.1\t5p\t+\t1\tC\t0\t85\t0\t15\t100\t0\t15\n")
            f.write("NC_012920.1\t5p\t+\t2\tC\t0\t96\t0\t4\t100\t0\t4\n")
            f.write("NC_012920.1\t3p\t+\t1\tG\t12\t0\t88\t0\t100\t12\t0\n")
            f.write("NC_012920.1\t3p\t+\t2\tG\t3\t0\t97\t0\t100\t3\t0\n")

        lg_file = tmp_path / "lgdistribution.txt"
        with open(lg_file, "w", encoding="utf-8") as f:
            f.write("Length\tCount\n")
            f.write("40\t10\n")
            f.write("50\t40\n")
            f.write("60\t30\n")
            f.write("70\t20\n")

        mcmc_file = tmp_path / "Stats_out_MCMC_iter_summ_stat.csv"
        with open(mcmc_file, "w", encoding="utf-8") as f:
            f.write("Parameter,Mean,Median,2.5%,97.5%,Std\n")
            f.write("Lambda,0.325,0.324,0.28,0.37,0.02\n")
            f.write("Delta_S,0.245,0.244,0.19,0.30,0.03\n")
            f.write("Delta_D,0.012,0.011,0.005,0.02,0.003\n")
            f.write("Theta,0.0015,0.0014,0.0008,0.002,0.0003\n")

        freqs, _ = DamageQuantifierEngine.parse_misincorporation_file(misinc_file, window_length=10)
        lg_metrics = DamageQuantifierEngine.parse_length_distribution_file(lg_file)
        mcmc_params = DamageQuantifierEngine.parse_mcmc_summary_file(mcmc_file)

        assert abs(freqs.c_to_t_5p_pos1 - 0.15) < 1e-4
        assert freqs.n_c_to_t_5p_pos1_obs == 15
        assert freqs.n_c_ref_5p_pos1_total == 100
        assert abs(freqs.g_to_a_3p_pos1 - 0.12) < 1e-4
        assert lg_metrics.total_fragments_sampled == 100
        assert 50 <= lg_metrics.median_length <= 60
        assert mcmc_params.model_converged is True
        assert abs(mcmc_params.delta_s - 0.245) < 1e-4


def test_pmd_scores_table_parser() -> None:
    """Tests PMD scoring table parser, handling header text, comments, and score percentiles."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        tsv_path = Path(tmp_dir) / "sample.pmd_scores.tsv"
        with open(tsv_path, "w", encoding="utf-8") as f:
            f.write("# PMDtools scores output\n")
            f.write("read_id\tscore\tlength\n")
            f.write("read_01\t5.82\t54\n")
            f.write("read_02\t5.12\t48\n")
            f.write("read_03\t3.91\t60\n")
            f.write("read_04\t3.10\t45\n")
            f.write("read_05\t1.50\t70\n")
            f.write("read_06\t0.80\t65\n")
            f.write("read_07\t0.10\t58\n")
            f.write("read_08\t-1.20\t80\n")
            f.write("read_09\t-2.50\t95\n")
            f.write("read_10\t-3.10\t110\n")

        metrics = PMDScoringEngine.parse_pmd_scores_table(tsv_path)

        assert metrics.total_reads_scored == 10
        assert abs(metrics.fraction_pmd_ge_0 - 0.70) < 1e-4
        assert abs(metrics.fraction_pmd_ge_3 - 0.40) < 1e-4
        assert abs(metrics.fraction_pmd_ge_5 - 0.20) < 1e-4
        assert metrics.min_pmd_score == -3.10
        assert metrics.max_pmd_score == 5.82


def test_qc_evaluator_decision_matrix() -> None:
    """Tests the multi-scenario biological decision matrix in DamageQCValidator."""
    policy_standard = DamagePolicy(
        library_treatment=LibraryTreatment.NONE,
        min_c_to_t_5p_strong_threshold=0.10,
        min_g_to_a_3p_strong_threshold=0.10,
        min_c_to_t_5p_observations=50,
        min_g_to_a_3p_observations=50,
    )
    
    # 1. Authentic Non-UDG Sample (Strong Signal)
    t_strong = TerminalDamageFrequencies(
        c_to_t_5p_pos1=0.14,
        g_to_a_3p_pos1=0.13,
        n_c_to_t_5p_pos1_obs=70,
        n_c_ref_5p_pos1_total=500,
        n_g_to_a_3p_pos1_obs=65,
        n_g_ref_3p_pos1_total=500,
    )
    bay_strong = BayesianDamageParameters(delta_s=0.22, model_converged=True, lambda_param=0.3, delta_d=0.01)
    len_short = FragmentLengthMetrics(total_fragments_sampled=500, median_length=55.0)

    res1 = DamageQCValidator.evaluate_sample(
        sample_id="SAMPLE_ANCIENT",
        input_bam_path=Path("sample.bam"),
        reference_contig="NC_012920.1",
        terminal_damage=t_strong,
        fragment_length=len_short,
        bayesian_params=bay_strong,
        pmd_metrics=PMDMetrics(total_reads_scored=500, fraction_pmd_ge_3=0.35),
        outputs=DamageOutputFiles(),
        policy=policy_standard,
    )
    assert res1.profile.damage_signal_level == DamageSignalLevel.STRONG_SIGNAL
    assert res1.profile.authentication_verdict == AuthenticationVerdict.AUTHENTIC_ANCIENT_SUPPORTED
    assert res1.status == DamageStatus.PASS

    # 2. Modern / Flat Sample
    t_flat = TerminalDamageFrequencies(
        c_to_t_5p_pos1=0.01,
        g_to_a_3p_pos1=0.01,
        n_c_to_t_5p_pos1_obs=5,
        n_c_ref_5p_pos1_total=500,
        n_g_to_a_3p_pos1_obs=5,
        n_g_ref_3p_pos1_total=500,
    )
    bay_flat = BayesianDamageParameters(delta_s=0.005, model_converged=True, lambda_param=0.3, delta_d=0.001)

    res2 = DamageQCValidator.evaluate_sample(
        sample_id="SAMPLE_MODERN",
        input_bam_path=Path("sample.bam"),
        reference_contig="NC_012920.1",
        terminal_damage=t_flat,
        fragment_length=FragmentLengthMetrics(total_fragments_sampled=500, median_length=120.0),
        bayesian_params=bay_flat,
        pmd_metrics=PMDMetrics(total_reads_scored=500, fraction_pmd_ge_3=0.01),
        outputs=DamageOutputFiles(),
        policy=policy_standard,
    )
    assert res2.profile.damage_signal_level == DamageSignalLevel.LOW_SIGNAL
    assert res2.profile.authentication_verdict == AuthenticationVerdict.AUTHENTICITY_UNSUPPORTED
    assert DamageWarningCode.LOW_DAMAGE_ON_NON_UDG in res2.warning_codes

    # 3. UDG_FULL Sample
    policy_udg_full = DamagePolicy(library_treatment=LibraryTreatment.UDG_FULL)
    res3 = DamageQCValidator.evaluate_sample(
        sample_id="SAMPLE_UDG_FULL",
        input_bam_path=Path("sample.bam"),
        reference_contig="NC_012920.1",
        terminal_damage=t_flat,
        fragment_length=len_short,
        bayesian_params=bay_flat,
        pmd_metrics=None,
        outputs=DamageOutputFiles(),
        policy=policy_udg_full,
    )
    assert res3.profile.authentication_verdict == AuthenticationVerdict.UDG_DEPLETED_EXPECTED

    # 4. Insufficient Power (n_c_ref < 50)
    t_insufficient = TerminalDamageFrequencies(
        c_to_t_5p_pos1=0.50,
        n_c_to_t_5p_pos1_obs=1,
        n_c_ref_5p_pos1_total=2,
    )
    res4 = DamageQCValidator.evaluate_sample(
        sample_id="SAMPLE_LOW_DEPTH",
        input_bam_path=Path("sample.bam"),
        reference_contig="NC_012920.1",
        terminal_damage=t_insufficient,
        fragment_length=FragmentLengthMetrics(),
        bayesian_params=BayesianDamageParameters(),
        pmd_metrics=None,
        outputs=DamageOutputFiles(),
        policy=policy_standard,
    )
    assert res4.profile.damage_signal_level == DamageSignalLevel.NOISY_OR_INSUFFICIENT
    assert res4.profile.authentication_verdict == AuthenticationVerdict.INCONCLUSIVE
    assert DamageWarningCode.INSUFFICIENT_TERMINAL_COVERAGE in res4.warning_codes


def test_synthetic_bam_and_preflight() -> None:
    """Tests input BAM resolution and pre-flight validator using a synthetic in-memory BAM."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        bam_path = tmp_path / "synthetic.dedup.bam"
        
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "NC_012920.1", "LN": 16569}],
            "RG": [{"ID": "SAMPLE1", "SM": "SAMPLE1", "PL": "ILLUMINA"}],
        }
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out_bam:
            for i in range(5):
                a = pysam.AlignedSegment()
                a.query_name = f"read_{i+1}"
                a.query_sequence = "ACGTACGTACGTACGT"
                a.flag = 0
                a.reference_id = 0
                a.reference_start = 100 + (i * 20)
                a.mapping_quality = 60
                a.cigar = ((0, 16),)
                a.query_qualities = array.array("B", [35] * 16)
                a.set_tag("MD", "16")
                a.set_tag("RG", "SAMPLE1")
                out_bam.write(a)

        mapped_count = DamageOrchestrator.validate_input_bam_preflight(
            bam_path=bam_path,
            expected_contig="NC_012920.1",
        )
        assert mapped_count == 5
        assert Path(f"{bam_path}.bai").is_file()

        policy = DamagePolicy(input_bam_type=DamageInputBAMType.FILTERED_DEDUP)
        resolved_bam = DamageOrchestrator.resolve_input_bam(
            sample_id="SYNTH",
            policy=policy,
            dedup_bam=bam_path,
            allow_fallback=False,
        )
        assert resolved_bam.resolve() == bam_path.resolve()