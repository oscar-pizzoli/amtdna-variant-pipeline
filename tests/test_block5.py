#!/usr/bin/env python3
"""Comprehensive Unit, Integration & State-Machine Test Suite for Block 5.

Covers domain invariant boundaries, parser edge-cases, synthetic BAM pileup
profiling, BAM preparation/MD tag checks, the full 2D QC gating matrix,
and orchestrator fallback state-machine execution paths.
"""

import copy
import json
import logging
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pysam

# Add repository root to Python path
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.modules.alignment.models import ReferenceInfo
from src.modules.contamination.dissent_engine import (
    DEFAULT_MT_MASKED_REGIONS,
    DissentEngine,
    DissentExecutionError,
    DissentExecutionMetadata,
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
    GatingRecommendation,
    ModelQualityDiagnostics,
    SchmutziPolicy,
    WrappingMetadata,
)
from src.modules.contamination.orchestrator import (
    ContaminationOrchestrationError,
    ContaminationOrchestrator,
)
from src.modules.contamination.qc_evaluator import (
    ContaminationQCValidator,
    _contamination_json_serializer,
)
from src.modules.contamination.schmutzi_engine import (
    SchmutziEngine,
    SchmutziExecutionError,
    SchmutziExecutionMetadata,
)
from src.modules.damage.models import LibraryTreatment

logging.basicConfig(level=logging.ERROR)


def create_synthetic_reference(fa_path: Path, contig_name: str = "NC_012920.1", length_bp: int = 1000) -> ReferenceInfo:
    """Creates a synthetic mitochondrial FASTA reference and builds its samtools index."""
    # Pattern: 1000 bp alternating standard bases
    seq = ("ACGT" * (length_bp // 4 + 1))[:length_bp]
    fa_path.write_text(f">{contig_name}\n{seq}\n", encoding="utf-8")
    pysam.faidx(str(fa_path))
    return ReferenceInfo(fasta_path=fa_path, contig_name=contig_name, length_bp=length_bp)


def create_synthetic_bam(
    bam_path: Path,
    reference_info: ReferenceInfo,
    read_count: int = 50,
    with_md: bool = True,
    heterogeneity_site: Optional[int] = None,
    heterogeneity_vaf: float = 0.0,
) -> Path:
    """Builds a valid, coordinate-sorted, and indexed BAM with controlled alignments."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": reference_info.contig_name, "LN": reference_info.length_bp}],
    }

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam_out:
        for i in range(read_count):
            rec = pysam.AlignedSegment()
            rec.query_name = f"read_{i:04d}"
            rec.reference_id = 0
            # Posizionamento monotonico crescente per rispettare l'ordinamento per coordinate
            rec.reference_start = 50 + (i // 5)
            rec.mapping_quality = 60
            rec.flag = 0  # Primary mapped single-end

            seq_len = 50
            seq = "ACGT" * 13
            seq = seq[:seq_len]

            # Iniezione allele minoritario sintetico
            if heterogeneity_site is not None and heterogeneity_vaf > 0.0:
                pos_in_read = heterogeneity_site - rec.reference_start - 1
                if 0 <= pos_in_read < seq_len:
                    if (i / read_count) < heterogeneity_vaf:
                        seq_list = list(seq)
                        seq_list[pos_in_read] = "T"
                        seq = "".join(seq_list)

            rec.query_sequence = seq
            rec.cigar = ((0, seq_len),)
            rec.query_qualities = [35] * seq_len

            if with_md:
                rec.set_tag("MD", "50", value_type="Z")

            bam_out.write(rec)

    pysam.index(str(bam_path))
    return bam_path


# ==============================================================================
# 1. UNIT TESTS: Domain Models, Boundaries & Invariants
# ==============================================================================

class TestBlock5ModelsAndInvariants(unittest.TestCase):
    """Verifies domain model invariants, boundary conditions, and JSON roundtrip serialization."""

    def test_estimate_probability_bounds_and_intervals(self):
        # Valid estimate
        est = ContaminationEstimate(point_estimate=0.02, ci_lower_95=0.01, ci_upper_95=0.03)
        self.assertAlmostEqual(est.ci_width, 0.02)

        # Invariant 1: Out of range [0.0, 1.0]
        with self.assertRaises(ValueError):
            ContaminationEstimate(point_estimate=1.05)
        with self.assertRaises(ValueError):
            ContaminationEstimate(point_estimate=-0.01)

        # Invariant 2: Lower CI > Upper CI
        with self.assertRaises(ValueError):
            ContaminationEstimate(point_estimate=0.04, ci_lower_95=0.06, ci_upper_95=0.02)

        # Invariant 3: Point estimate strictly outside CI bounds
        with self.assertRaises(ValueError):
            ContaminationEstimate(point_estimate=0.10, ci_lower_95=0.01, ci_upper_95=0.05)

    def test_non_bayesian_none_ci_handling(self):
        est = ContaminationEstimate(
            point_estimate=0.025,
            ci_lower_95=None,
            ci_upper_95=None,
            method_used=ContaminationMethod.DISSENT_SCREENING,
        )
        self.assertIsNone(est.ci_width)
        d = est.to_dict()
        self.assertIsNone(d["ci_lower_95_fraction"])
        self.assertIsNone(d["ci_upper_95_fraction"])
        self.assertIsNone(d["ci_width"])
        self.assertEqual(d["point_estimate_fraction"], 0.025)

    def test_policy_threshold_validations(self):
        # Inverted category boundaries
        with self.assertRaises(ValueError):
            ContaminationPolicy(low_threshold=0.05, moderate_threshold=0.02)

        # Recommended depth lower than minimum floor
        with self.assertRaises(ValueError):
            ContaminationPolicy(min_coverage_for_estimation=10.0, recommended_coverage=5.0)

        # Inverted CI width criteria
        with self.assertRaises(ValueError):
            ContaminationPolicy(high_confidence_max_ci_width=0.08, moderate_confidence_max_ci_width=0.04)

        # Non-positive informative sites
        with self.assertRaises(ValueError):
            ContaminationPolicy(min_informative_sites=0)

    def test_json_deliverable_serialization_roundtrip(self):
        """Ensures that ContaminationResult serializes to valid JSON without enum/Path leaks."""
        policy = ContaminationPolicy()
        est = ContaminationEstimate(point_estimate=0.015, ci_lower_95=0.010, ci_upper_95=0.020, method_used=ContaminationMethod.SCHMUTZI)
        diag = ModelQualityDiagnostics(mean_depth=12.5, informative_positions_evaluated=40, model_converged=True, ci_width=0.010)
        gating = GatingRecommendation(
            contamination_category=ContaminationCategory.LOW,
            recommended_action=ConsensusGatingAction.STANDARD_CALLING,
            preferred_consensus_source="STANDARD_CALLER",
        )
        audit = ContaminationAuditMetadata(input_bam_path=Path("/tmp/sample.bam"), schmutzi_version="v1.5.6")
        outputs = ContaminationOutputFiles(endogenous_fasta=Path("/tmp/sample_endo.fasta"))

        res = ContaminationResult(
            sample_id="SYNTH_01",
            status=ContaminationStatus.PASS,
            policy=policy,
            estimate=est,
            diagnostics=diag,
            gating=gating,
            outputs=outputs,
            audit=audit,
        )

        res.mark_started()
        res.mark_finished()

        raw_dict = res.to_dict()
        json_str = json.dumps(raw_dict, indent=2, default=_contamination_json_serializer)
        reloaded = json.loads(json_str)

        self.assertEqual(reloaded["sample_id"], "SYNTH_01")
        self.assertEqual(reloaded["status"], "PASS")
        self.assertEqual(reloaded["estimate"]["point_estimate_fraction"], 0.015)
        self.assertEqual(reloaded["gating"]["recommended_action"], "STANDARD_CALLING")
        self.assertIsInstance(reloaded["timestamps"]["runtime_seconds"], (int, float))


# ==============================================================================
# 2. UNIT TESTS: Schmutzi Parsers & Diagnostics
# ==============================================================================

class TestSchmutziParsers(unittest.TestCase):
    """Tests robustness of schmutzi parsers on nominal, percentage, scientific, and corrupt formats."""

    def setUp(self):
        self.test_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_parse_cont_est_scientific_and_percentages(self):
        cont_est = self.test_dir / "sample.cont.est"
        # Percentages: 3.2% [1.1%, 5.4%] with scientific log-likelihood
        cont_est.write_text("3.2 1.1 5.4 -1.25e+03\n", encoding="utf-8")

        point, low, high, log_lik = SchmutziEngine.parse_cont_est_file(cont_est)
        self.assertAlmostEqual(point, 0.032)
        self.assertAlmostEqual(low, 0.011)
        self.assertAlmostEqual(high, 0.054)
        self.assertAlmostEqual(log_lik, -1250.0)

    def test_parse_cont_est_gross_statistical_inconsistency_fails(self):
        cont_est = self.test_dir / "sample.cont.est"
        # Point estimate (0.75) is radically outside [0.01, 0.03] -> must raise error
        cont_est.write_text("0.75 0.01 0.03 -500.0\n", encoding="utf-8")

        with self.assertRaises(SchmutziExecutionError) as ctx:
            SchmutziEngine.parse_cont_est_file(cont_est)
        self.assertIn("Gross statistical inconsistency", str(ctx.exception))

    def test_parse_cont_est_corrupted_empty_file_fails(self):
        empty_file = self.test_dir / "empty.cont.est"
        empty_file.write_text("", encoding="utf-8")

        with self.assertRaises(SchmutziExecutionError):
            SchmutziEngine.parse_cont_est_file(empty_file)

    def test_parse_schmutzi_log_convergence_and_failure(self):
        log_pass = self.test_dir / "pass.log"
        log_pass.write_text(
            "Starting MCMC sampler\n"
            "iteration 1\n"
            "iteration 2\n"
            "posterior estimation complete\n"
            "predicted haplogroup: U5b1\n"
            "35 informative positions\n",
            encoding="utf-8",
        )

        conv, iters, hgs, info_sites, warn_dmg, multi_cont = SchmutziEngine.parse_schmutzi_log(log_pass)
        self.assertTrue(conv)
        self.assertEqual(iters, 2)
        self.assertIn("U5b1", hgs)
        self.assertEqual(info_sites, 35)
        self.assertFalse(warn_dmg)
        self.assertFalse(multi_cont)

        # Log without explicit convergence keywords
        log_fail = self.test_dir / "fail.log"
        log_fail.write_text("Starting iteration 1\nUnexpected process halt.\n", encoding="utf-8")
        conv_f, iters_f, _, _, _, _ = SchmutziEngine.parse_schmutzi_log(log_fail)
        self.assertFalse(conv_f)


# ==============================================================================
# 3. INTEGRATION TESTS: Input Preparer & MD Tags
# ==============================================================================

class TestSchmutziInputPreparer(unittest.TestCase):
    """Tests BAM sanitization, MD tag verification, and reference compatibility."""

    def setUp(self):
        self.test_dir = Path(tempfile.mkdtemp())
        self.ref_fa = self.test_dir / "ref.fasta"
        self.ref_info = create_synthetic_reference(self.ref_fa, contig_name="NC_012920.1", length_bp=1000)

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_check_md_tags_presence(self):
        bam_with_md = create_synthetic_bam(self.test_dir / "with_md.bam", self.ref_info, with_md=True)
        bam_without_md = create_synthetic_bam(self.test_dir / "no_md.bam", self.ref_info, with_md=False)

        self.assertTrue(SchmutziInputPreparer.sample_md_tags_present(bam_with_md))
        self.assertFalse(SchmutziInputPreparer.sample_md_tags_present(bam_without_md))

    def test_prepare_bam_symlink_when_md_present(self):
        bam_with_md = create_synthetic_bam(self.test_dir / "valid.bam", self.ref_info, with_md=True)
        work_dir = self.test_dir / "prep_work"

        prep_bam, prep_bai = SchmutziInputPreparer.prepare_bam_for_schmutzi(
            sample_id="SAMPLE_MD",
            input_bam=bam_with_md,
            reference_info=self.ref_info,
            work_dir=work_dir,
            force_calmd=False,
        )

        self.assertTrue(prep_bam.is_file())
        self.assertTrue(prep_bai.is_file())
        self.assertEqual(prep_bam.name, "SAMPLE_MD.schmutzi_ready.bam")

    def test_reference_contig_mismatch_fails_fast(self):
        # Reference expects NC_012920.1, BAM contains chrM
        other_ref = ReferenceInfo(fasta_path=self.ref_fa, contig_name="chrM", length_bp=1000)
        bam_nc = create_synthetic_bam(self.test_dir / "nc.bam", self.ref_info, with_md=True)

        with self.assertRaises(SchmutziPreparationError) as ctx:
            SchmutziInputPreparer.prepare_bam_for_schmutzi(
                sample_id="MISMATCH",
                input_bam=bam_nc,
                reference_info=other_ref,
                work_dir=self.test_dir / "prep_work",
            )
        self.assertIn("missing from BAM header", str(ctx.exception))


# ==============================================================================
# 4. INTEGRATION TESTS: Dissent Engine End-to-End on Synthetic BAM
# ==============================================================================

class TestDissentEngineEndToEnd(unittest.TestCase):
    """Tests Track B pileup scanning, minor allele calculations, and boundary masking."""

    def setUp(self):
        self.test_dir = Path(tempfile.mkdtemp())
        self.ref_fa = self.test_dir / "ref.fasta"
        self.ref_info = create_synthetic_reference(self.ref_fa, contig_name="NC_012920.1", length_bp=1000)

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_homopolymer_exact_boundaries(self):
        # Default region: (303, 315)
        # Position 302: outside -> False
        self.assertFalse(DissentEngine.is_masked_position(302, DEFAULT_MT_MASKED_REGIONS))
        # Position 303: boundary start -> True
        self.assertTrue(DissentEngine.is_masked_position(303, DEFAULT_MT_MASKED_REGIONS))
        # Position 310: inside -> True
        self.assertTrue(DissentEngine.is_masked_position(310, DEFAULT_MT_MASKED_REGIONS))
        # Position 315: boundary end -> True
        self.assertTrue(DissentEngine.is_masked_position(315, DEFAULT_MT_MASKED_REGIONS))
        # Position 316: outside -> False
        self.assertFalse(DissentEngine.is_masked_position(316, DEFAULT_MT_MASKED_REGIONS))

    def test_dissent_synthetic_pileup_minor_allele_detection(self):
        # Create BAM with 100 reads, and 10% VAF at position 60
        bam_path = self.test_dir / "hetero.bam"
        create_synthetic_bam(
            bam_path=bam_path,
            reference_info=self.ref_info,
            read_count=100,
            with_md=True,
            heterogeneity_site=60,
            heterogeneity_vaf=0.10,
        )

        out_dir = self.test_dir / "dissent_out"
        log_dir = self.test_dir / "dissent_log"

        est, diag, out_files, meta = DissentEngine.execute_dissent_screening(
            sample_id="SYNTH_DISSENT",
            input_bam=bam_path,
            reference_fasta=self.ref_fa,
            output_dir=out_dir,
            log_dir=log_dir,
            expected_contig=self.ref_info.contig_name,
            min_depth_per_site=5,
            noise_vaf_floor=0.015,
            base_min_dissenting_sites=1,
            target_positions={60},  # Targeted scan on site 60
        )

        self.assertEqual(est.method_used, ContaminationMethod.DISSENT_SCREENING)
        self.assertIsNone(est.ci_lower_95)
        self.assertIsNone(est.ci_upper_95)
        self.assertGreater(est.point_estimate, 0.05)
        self.assertTrue(out_files.dissent_tsv.is_file())
        self.assertEqual(meta.informative_positions_evaluated, 1)
        self.assertEqual(meta.dissenting_positions_count, 1)


# ==============================================================================
# 5. UNIT TESTS: Full 2D QC Gating Matrix (Category x Confidence)
# ==============================================================================

class TestQCValidatorFullMatrix(unittest.TestCase):
    """Exhaustively validates all cells of the 2D Gating Matrix (Category x Confidence)."""

    def setUp(self):
        self.policy = ContaminationPolicy(
            low_threshold=0.02,
            moderate_threshold=0.05,
            elevated_threshold=0.10,
            pmd_filtering_threshold=3.0,
        )

    def test_full_decision_grid(self):
        grid_cases = [
            # (Category, Confidence, ExpectedAction, ExpectedSource)
            (ContaminationCategory.LOW, EstimationConfidence.HIGH, ConsensusGatingAction.STANDARD_CALLING, "STANDARD_CALLER"),
            (ContaminationCategory.LOW, EstimationConfidence.MODERATE, ConsensusGatingAction.STANDARD_CALLING, "STANDARD_CALLER"),
            (ContaminationCategory.LOW, EstimationConfidence.LOW, ConsensusGatingAction.RESTRICT_HETEROPLASMY, "STANDARD_CALLER"),
            (ContaminationCategory.LOW, EstimationConfidence.UNRELIABLE, ConsensusGatingAction.EXCLUDE_OR_MANUAL_REVIEW, None),

            (ContaminationCategory.MODERATE, EstimationConfidence.HIGH, ConsensusGatingAction.RESTRICT_HETEROPLASMY, "STANDARD_CALLER"),
            (ContaminationCategory.MODERATE, EstimationConfidence.MODERATE, ConsensusGatingAction.RESTRICT_HETEROPLASMY, "STANDARD_CALLER"),
            (ContaminationCategory.MODERATE, EstimationConfidence.LOW, ConsensusGatingAction.PMD_FILTERING_REQUIRED, "PMD_FILTERED_BAM"),
            (ContaminationCategory.MODERATE, EstimationConfidence.UNRELIABLE, ConsensusGatingAction.EXCLUDE_OR_MANUAL_REVIEW, None),

            (ContaminationCategory.ELEVATED, EstimationConfidence.HIGH, ConsensusGatingAction.USE_ENDOGENOUS_SCHMUTZI_CONSENSUS, "SCHMUTZI_ENDOGENOUS_FASTA"),
            (ContaminationCategory.ELEVATED, EstimationConfidence.LOW, ConsensusGatingAction.EXCLUDE_OR_MANUAL_REVIEW, None),

            (ContaminationCategory.CRITICAL, EstimationConfidence.HIGH, ConsensusGatingAction.EXCLUDE_OR_MANUAL_REVIEW, None),
            (ContaminationCategory.CRITICAL, EstimationConfidence.UNRELIABLE, ConsensusGatingAction.EXCLUDE_OR_MANUAL_REVIEW, None),

            (ContaminationCategory.INSUFFICIENT_DATA, EstimationConfidence.LOW, ConsensusGatingAction.EXCLUDE_OR_MANUAL_REVIEW, None),
            (ContaminationCategory.INSUFFICIENT_DATA, EstimationConfidence.UNRELIABLE, ConsensusGatingAction.EXCLUDE_OR_MANUAL_REVIEW, None),
        ]

        dummy_endo = Path(tempfile.gettempdir()) / "dummy_endo.fasta"
        dummy_endo.write_text(">NC_012920.1\nACGT\n", encoding="utf-8")
        outputs = ContaminationOutputFiles(endogenous_fasta=dummy_endo)

        for cat, conf, exp_action, exp_source in grid_cases:
            with self.subTest(category=cat.value, confidence=conf.value):
                est = ContaminationEstimate(point_estimate=0.01, method_used=ContaminationMethod.SCHMUTZI)
                rec = ContaminationQCValidator.formulate_gating_recommendation(
                    category=cat,
                    confidence=conf,
                    estimate=est,
                    outputs=outputs,
                    policy=self.policy,
                )
                self.assertEqual(rec.recommended_action, exp_action)
                self.assertEqual(rec.preferred_consensus_source, exp_source)

    def test_exact_threshold_boundaries(self):
        diag = ModelQualityDiagnostics(informative_positions_evaluated=20)
        # Boundary 0.0200 -> LOW
        c1, _, _ = ContaminationQCValidator.categorize_contamination(0.0200, EstimationConfidence.HIGH, diag, self.policy)
        self.assertEqual(c1, ContaminationCategory.LOW)

        # Boundary 0.0201 -> MODERATE
        c2, _, _ = ContaminationQCValidator.categorize_contamination(0.0201, EstimationConfidence.HIGH, diag, self.policy)
        self.assertEqual(c2, ContaminationCategory.MODERATE)

        # Boundary 0.0500 -> MODERATE
        c3, _, _ = ContaminationQCValidator.categorize_contamination(0.0500, EstimationConfidence.HIGH, diag, self.policy)
        self.assertEqual(c3, ContaminationCategory.MODERATE)

        # Boundary 0.0501 -> ELEVATED
        c4, _, _ = ContaminationQCValidator.categorize_contamination(0.0501, EstimationConfidence.HIGH, diag, self.policy)
        self.assertEqual(c4, ContaminationCategory.ELEVATED)

        # Boundary 0.1000 -> ELEVATED
        c5, _, _ = ContaminationQCValidator.categorize_contamination(0.1000, EstimationConfidence.HIGH, diag, self.policy)
        self.assertEqual(c5, ContaminationCategory.ELEVATED)

        # Boundary 0.1001 -> CRITICAL
        c6, _, _ = ContaminationQCValidator.categorize_contamination(0.1001, EstimationConfidence.HIGH, diag, self.policy)
        self.assertEqual(c6, ContaminationCategory.CRITICAL)


# ==============================================================================
# 6. SYSTEM TESTS: Orchestrator State-Machine & Fallback Pathways
# ==============================================================================

class TestOrchestratorStateMachine(unittest.TestCase):
    """Tests the full orchestration state machine: nominal execution, database failures, and fallbacks."""

    def setUp(self):
        self.test_dir = Path(tempfile.mkdtemp())
        self.out_dir = self.test_dir / "results"
        self.log_dir = self.test_dir / "logs"

        self.ref_fa = self.test_dir / "ref.fasta"
        self.ref_info = create_synthetic_reference(self.ref_fa, contig_name="NC_012920.1", length_bp=1000)
        self.bam_path = create_synthetic_bam(self.test_dir / "sample.bam", self.ref_info, read_count=60, with_md=True)

        self.valid_db = self.test_dir / "valid_db.fa"
        self.valid_db.write_text(">hg1\nACGT\n>hg2\nTGCA\n", encoding="utf-8")

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_missing_input_bam_halts_fast(self):
        res = ContaminationOrchestrator.run_sample_contamination_pipeline(
            sample_id="ERR_BAM",
            reference_info=self.ref_info,
            input_bam=self.test_dir / "non_existent.bam",
            output_dir=self.out_dir,
            log_dir=self.log_dir,
            policy=ContaminationPolicy(),
        )
        self.assertEqual(res.status, ContaminationStatus.FAIL)
        self.assertIn(ContaminationErrorCode.INPUT_BAM_NOT_FOUND, res.error_codes)

    @patch("src.modules.contamination.schmutzi_engine.SchmutziEngine.get_binary_name", return_value="schmutzi")
    @patch("src.modules.contamination.schmutzi_engine.SchmutziEngine.run_schmutzi_pipeline")
    def test_schmutzi_nominal_success(self, mock_schmutzi, mock_bin):
        # Mock successful schmutzi execution
        est = ContaminationEstimate(point_estimate=0.012, ci_lower_95=0.008, ci_upper_95=0.016, method_used=ContaminationMethod.SCHMUTZI)
        diag = ModelQualityDiagnostics(mean_depth=15.0, informative_positions_evaluated=30, ci_width=0.008, model_converged=True)
        out_f = ContaminationOutputFiles()
        audit = SchmutziExecutionMetadata(tool_version="v1.5.6", executed_command="schmutzi ...")
        mock_schmutzi.return_value = (est, diag, out_f, audit)

        policy = ContaminationPolicy(
            method=ContaminationMethod.SCHMUTZI,
            schmutzi=SchmutziPolicy(contaminant_database=self.valid_db),
        )

        res = ContaminationOrchestrator.run_sample_contamination_pipeline(
            sample_id="SAMPLE_SUCCESS",
            reference_info=self.ref_info,
            input_bam=self.bam_path,
            output_dir=self.out_dir,
            log_dir=self.log_dir,
            policy=policy,
        )

        self.assertEqual(res.status, ContaminationStatus.PASS)
        self.assertEqual(res.estimate.method_used, ContaminationMethod.SCHMUTZI)
        self.assertEqual(res.gating.recommended_action, ConsensusGatingAction.STANDARD_CALLING)
        self.assertTrue(res.outputs.contamination_json.is_file())
        self.assertIsNotNone(res.runtime_seconds)

    @patch("src.modules.contamination.schmutzi_engine.SchmutziEngine.get_binary_name", return_value="schmutzi")
    @patch("src.modules.contamination.schmutzi_engine.SchmutziEngine.run_schmutzi_pipeline")
    def test_schmutzi_runtime_crash_triggers_fallback_when_allowed(self, mock_schmutzi, mock_bin):
        # Simulate schmutzi crashing with an exception
        mock_schmutzi.side_effect = SchmutziExecutionError("Process returned exit code 139 (SIGSEGV)")

        policy = ContaminationPolicy(
            method=ContaminationMethod.SCHMUTZI,
            allow_fallback_to_dissent=True,
            schmutzi=SchmutziPolicy(contaminant_database=self.valid_db),
        )

        res = ContaminationOrchestrator.run_sample_contamination_pipeline(
            sample_id="SAMPLE_CRASH_FALLBACK",
            reference_info=self.ref_info,
            input_bam=self.bam_path,
            output_dir=self.out_dir,
            log_dir=self.log_dir,
            policy=policy,
        )

        self.assertEqual(res.estimate.method_used, ContaminationMethod.DISSENT_SCREENING)
        self.assertIn(ContaminationWarningCode.FALLBACK_TO_DISSENT, res.warning_codes)
        self.assertIn("SIGSEGV", res.fallback_reason)

    @patch("src.modules.contamination.schmutzi_engine.SchmutziEngine.get_binary_name", return_value="schmutzi")
    @patch("src.modules.contamination.schmutzi_engine.SchmutziEngine.run_schmutzi_pipeline")
    def test_schmutzi_runtime_crash_fails_when_fallback_forbidden(self, mock_schmutzi, mock_bin):
        mock_schmutzi.side_effect = SchmutziExecutionError("MCMC internal divergence")

        policy = ContaminationPolicy(
            method=ContaminationMethod.SCHMUTZI,
            allow_fallback_to_dissent=False,
            schmutzi=SchmutziPolicy(contaminant_database=self.valid_db),
        )

        res = ContaminationOrchestrator.run_sample_contamination_pipeline(
            sample_id="SAMPLE_STRICT_FAIL",
            reference_info=self.ref_info,
            input_bam=self.bam_path,
            output_dir=self.out_dir,
            log_dir=self.log_dir,
            policy=policy,
        )

        self.assertEqual(res.status, ContaminationStatus.FAIL)
        self.assertIn(ContaminationErrorCode.SCHMUTZI_EXECUTION_FAILED, res.error_codes)


if __name__ == "__main__":
    unittest.main(verbosity=2)