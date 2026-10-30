"""
Exhaustive smoke and integration test suite for Block 6 (Consensus & Heteroplasmy Pipeline).
Organized into three strict operational levels:
- Level 1: Unit Tests (Enums, Models, Deterministic Tie-Breaking, Policy Invariants)
- Level 2: Component Integration Tests (PileupEngine, ConsensusBuilder, HeteroplasmyCaller)
- Level 3: End-to-End Pipeline & Orchestrator Integration Tests
"""

import gzip
import json
import os
from typing import Tuple
import pytest
import pysam

from src.modules.consensus.enums import (
    CallStatus,
    ConsensusCallReason,
    ConsensusCallStatus,
    ContaminationEstimateSource,
    FilterReason,
    GatingAction,
    InputRole,
    IUPACAmbiguity,
)
from src.modules.consensus.models import (
    AlleleCounts,
    ConsensusPolicy,
    InputResolution,
    ReferenceInfo,
    SitePileup,
)
from src.modules.consensus.pileup_engine import PileupEngine
from src.modules.consensus.consensus_builder import ConsensusBuilder
from src.modules.consensus.heteroplasmy_caller import HeteroplasmyCaller
from src.modules.consensus.orchestrator import (
    ConsensusOrchestrator,
    ContractValidationError,
    UnresolvedAuthenticityError,
)


# ============================================================================
# Synthetic Test Fixtures: Genome & Alignment Construction
# ============================================================================

@pytest.fixture(scope="module")
def synthetic_environment(tmp_path_factory) -> Tuple[str, str, str]:
    """
    Build a realistic 16,569 bp synthetic mitochondrial reference FASTA (with rCRS pos 3107 'N')
    and a targeted indexed BAM containing designated quality, flag, and overlap test cases.
    """
    tmp_dir = tmp_path_factory.mktemp("block6_smoke_env")
    fasta_path = str(tmp_dir / "rCRS_synthetic.fasta")
    bam_path = str(tmp_dir / "sample_synthetic.bam")
    contig_name = "NC_012920.1"
    ref_length = 16569

    # 1. Synthesize reference sequence: mostly 'A' with legacy 'N' at Cambridge position 3107
    ref_chars = ["A"] * ref_length
    ref_chars[3106] = "N"  # Position 3107 (1-based) is index 3106
    # Anchors for calling verification
    ref_chars[9] = "C"     # Pos 10: Reference C (majority consensus test)
    ref_chars[19] = "G"    # Pos 20: Reference G (heteroplasmy test)
    ref_chars[29] = "T"    # Pos 30: Reference T (strand bias test)
    ref_chars[39] = "A"    # Pos 40: Reference A (low-frequency minor allele test)

    ref_seq = "".join(ref_chars)
    with open(fasta_path, "w", encoding="utf-8") as f:
        f.write(f">{contig_name} synthetic rCRS\n")
        for i in range(0, len(ref_seq), 60):
            f.write(ref_seq[i : i + 60] + "\n")

    pysam.faidx(fasta_path)

    # 2. Build synthetic BAM with 8 bp reads so test loci (spaced 10 bp apart) do not bleed
    header = {
        "HD": {"VN": "1.0", "SO": "coordinate"},
        "SQ": [{"LN": ref_length, "SN": contig_name}],
    }

    sam_path = str(tmp_dir / "temp.sam")
    with pysam.AlignmentFile(sam_path, "wh", header=header) as sam:
        # Position 10 (index 9): High Depth Majority Consensus (20 reads 'C')
        for i in range(20):
            r = pysam.AlignedSegment()
            r.query_name = f"read_pos10_{i}"
            r.query_sequence = "C" * 8
            r.flag = 0
            r.reference_id = 0
            r.reference_start = 9
            r.mapping_quality = 40
            r.cigar = ((0, 8),)
            r.query_qualities = pysam.qualitystring_to_array("I" * 8)
            sam.write(r)

        # Position 20 (index 19): Balanced Heteroplasmy: 12 'G', 8 'A' (4 FWD, 4 REV)
        for i in range(12):
            r = pysam.AlignedSegment()
            r.query_name = f"read_pos20_maj_{i}"
            r.query_sequence = "G" * 8
            r.flag = 0 if i % 2 == 0 else 16
            r.reference_id = 0
            r.reference_start = 19
            r.mapping_quality = 40
            r.cigar = ((0, 8),)
            r.query_qualities = pysam.qualitystring_to_array("I" * 8)
            sam.write(r)

        for i in range(8):
            r = pysam.AlignedSegment()
            r.query_name = f"read_pos20_min_{i}"
            r.query_sequence = "A" * 8
            r.flag = 0 if i < 4 else 16
            r.reference_id = 0
            r.reference_start = 19
            r.mapping_quality = 40
            r.cigar = ((0, 8),)
            r.query_qualities = pysam.qualitystring_to_array("I" * 8)
            sam.write(r)

        # Position 30 (index 29): Strand Bias: 15 'T' (ref), 5 'C' (alt strictly FWD)
        for i in range(15):
            r = pysam.AlignedSegment()
            r.query_name = f"read_pos30_maj_{i}"
            r.query_sequence = "T" * 8
            r.flag = 0 if i % 2 == 0 else 16
            r.reference_id = 0
            r.reference_start = 29
            r.mapping_quality = 40
            r.cigar = ((0, 8),)
            r.query_qualities = pysam.qualitystring_to_array("I" * 8)
            sam.write(r)

        for i in range(5):
            r = pysam.AlignedSegment()
            r.query_name = f"read_pos30_min_{i}"
            r.query_sequence = "C" * 8
            r.flag = 0
            r.reference_id = 0
            r.reference_start = 29
            r.mapping_quality = 40
            r.cigar = ((0, 8),)
            r.query_qualities = pysam.qualitystring_to_array("I" * 8)
            sam.write(r)

        # Position 40 (index 39): Dominant Majority (18 'A', 2 'G', VAF=0.90)
        for i in range(18):
            r = pysam.AlignedSegment()
            r.query_name = f"read_pos40_maj_{i}"
            r.query_sequence = "A" * 8
            r.flag = 0 if i % 2 == 0 else 16
            r.reference_id = 0
            r.reference_start = 39
            r.mapping_quality = 40
            r.cigar = ((0, 8),)
            r.query_qualities = pysam.qualitystring_to_array("I" * 8)
            sam.write(r)

        for i in range(2):
            r = pysam.AlignedSegment()
            r.query_name = f"read_pos40_min_{i}"
            r.query_sequence = "G" * 8
            r.flag = 0 if i == 0 else 16
            r.reference_id = 0
            r.reference_start = 39
            r.mapping_quality = 40
            r.cigar = ((0, 8),)
            r.query_qualities = pysam.qualitystring_to_array("I" * 8)
            sam.write(r)

        # Position 50 (index 49): Overlapping Mate Pair Test (Mate 1 BQ=35, Mate 2 BQ=28 >= min_bq 25)
        # Position 50 (index 49): Overlapping Mate Pair Test (Mate 1 BQ=35, Mate 2 BQ=30 >= min_bq 25)
        r1 = pysam.AlignedSegment()
        r1.query_name = "mate_pair_overlap_1"
        r1.query_sequence = "T" * 8
        r1.flag = 99
        r1.reference_id = 0
        r1.reference_start = 49
        r1.mapping_quality = 40
        r1.cigar = ((0, 8),)
        r1.query_qualities = pysam.qualitystring_to_array("D" * 8)  # ord('D')=68 -> BQ=35
        sam.write(r1)

        r2 = pysam.AlignedSegment()
        r2.query_name = "mate_pair_overlap_1"
        r2.query_sequence = "T" * 8
        r2.flag = 147
        r2.reference_id = 0
        r2.reference_start = 49
        r2.mapping_quality = 40
        r2.cigar = ((0, 8),)
        r2.query_qualities = pysam.qualitystring_to_array("?" * 8)  # ord('?')=63 -> BQ=30 (>= 25)
        sam.write(r2)

        # Position 60 (index 59): Filter Exclusions (Valid, QC-Fail, Duplicate)
        r_valid = pysam.AlignedSegment()
        r_valid.query_name = "r_valid"
        r_valid.query_sequence = "A" * 8
        r_valid.flag = 0
        r_valid.reference_id = 0
        r_valid.reference_start = 59
        r_valid.mapping_quality = 40
        r_valid.cigar = ((0, 8),)
        r_valid.query_qualities = pysam.qualitystring_to_array("I" * 8)
        sam.write(r_valid)

        r_qcfail = pysam.AlignedSegment()
        r_qcfail.query_name = "r_qcfail"
        r_qcfail.query_sequence = "A" * 8
        r_qcfail.flag = 512
        r_qcfail.reference_id = 0
        r_qcfail.reference_start = 59
        r_qcfail.mapping_quality = 40
        r_qcfail.cigar = ((0, 8),)
        r_qcfail.query_qualities = pysam.qualitystring_to_array("I" * 8)
        sam.write(r_qcfail)

        r_dup = pysam.AlignedSegment()
        r_dup.query_name = "r_dup"
        r_dup.query_sequence = "A" * 8
        r_dup.flag = 1024
        r_dup.reference_id = 0
        r_dup.reference_start = 59
        r_dup.mapping_quality = 40
        r_dup.cigar = ((0, 8),)
        r_dup.query_qualities = pysam.qualitystring_to_array("I" * 8)
        sam.write(r_dup)

    pysam.sort("-o", bam_path, sam_path)
    pysam.index(bam_path)
    os.remove(sam_path)

    return fasta_path, bam_path, contig_name


# ============================================================================
# Level 1: Unit Tests (Models, Enums, Deterministic Tie-Breaking)
# ============================================================================

def test_unit_iupac_ambiguity_resolution():
    """Verify deterministic IUPAC resolution without collapsing identical bases to N."""
    assert IUPACAmbiguity.resolve("A", "A") == "A"
    assert IUPACAmbiguity.resolve("C", "C") == "C"
    assert IUPACAmbiguity.resolve("A", "G") == "R"
    assert IUPACAmbiguity.resolve("G", "A") == "R"
    assert IUPACAmbiguity.resolve("C", "T") == "Y"
    assert IUPACAmbiguity.resolve("G", "C") == "S"
    assert IUPACAmbiguity.resolve("A", "T") == "W"
    assert IUPACAmbiguity.resolve("G", "T") == "K"
    assert IUPACAmbiguity.resolve("A", "C") == "M"
    assert IUPACAmbiguity.resolve("N", "A") == "N"


def test_unit_deterministic_tie_break_ranking():
    """Verify ranking priorities: (1) count descending, (2) match ref, (3) alphabetical."""
    # Tie between A and G, reference is G -> G must win
    counts = AlleleCounts(A=10, G=10, C=0, T=0)
    ranked = ConsensusBuilder.rank_canonical_bases(counts, ref_base="G")
    assert ranked[0][0] == "G"
    assert ranked[1][0] == "A"

    # Tie between A and C, reference is T -> A must win by priority order (A > C > G > T)
    ranked_no_ref = ConsensusBuilder.rank_canonical_bases(counts, ref_base="T")
    assert ranked_no_ref[0][0] == "A"
    assert ranked_no_ref[1][0] == "G"


def test_unit_consensus_policy_invariants():
    """Verify ConsensusPolicy validates boundaries and parameter ranges."""
    with pytest.raises(ValueError, match="mask_char == 'N'"):
        ConsensusPolicy(mask_char="X")

    with pytest.raises(ValueError, match="major_allele_fraction"):
        ConsensusPolicy(major_allele_fraction=0.45)

    with pytest.raises(ValueError, match="ambiguity_min_vaf"):
        ConsensusPolicy(ambiguity_min_vaf=0.55)


# ============================================================================
# Level 2: Component Integration Tests
# ============================================================================

def test_component_pileup_engine(synthetic_environment):
    """Verify gapless streaming, single-molecule mate overlap deduplication, and QC filters."""
    fasta_path, bam_path, contig_name = synthetic_environment
    policy = ConsensusPolicy(min_base_quality=25, min_mapping_quality=30)
    engine = PileupEngine(
        bam_path=bam_path,
        reference_fasta_path=fasta_path,
        policy=policy,
        target_contig=contig_name,
        expected_reference_length=16569,
    )

    sites = list(engine.iter_sites())
    assert len(sites) == 16569
    assert [s.position for s in sites[:5]] == [1, 2, 3, 4, 5]

    # Verify mate overlap deduplication at position 50
    site_50 = sites[49]
    assert site_50.position == 50
    assert site_50.depth_pileup_entries == 2
    assert site_50.depth_filtered == 1
    assert site_50.depth_callable == 1
    assert site_50.counts.T == 1
    assert site_50.filter_metrics.overlapping_mate == 1

    # Verify QC filter tracking at position 60
    site_60 = sites[59]
    assert site_60.depth_pileup_entries == 3
    assert site_60.depth_callable == 1
    assert site_60.filter_metrics.qc_fail == 1
    assert site_60.filter_metrics.duplicates == 1


def test_component_consensus_builder(synthetic_environment):
    """Verify consensus builder across majority, contested, low-frequency alt, and rCRS N."""
    fasta_path, bam_path, contig_name = synthetic_environment
    ref_info = ReferenceInfo(contig=contig_name, length=16569, sequence_sha256="dummy")
    policy = ConsensusPolicy(min_consensus_depth=10, major_allele_fraction=0.80)
    builder = ConsensusBuilder(policy=policy, sample_id="COMP_TEST", reference_info=ref_info)

    engine = PileupEngine(bam_path, fasta_path, policy, contig_name, 16569)
    sites = list(engine.iter_sites())
    consensus_sites, canonical_seq, iupac_seq = builder.build_consensus(sites)

    assert len(canonical_seq) == 16569
    assert len(iupac_seq) == 16569

    # Pos 10: 100% C -> CALLED_MAJORITY
    s10 = consensus_sites[9]
    assert s10.canonical_base == "C"
    assert s10.iupac_base == "C"
    assert s10.call_reason == ConsensusCallReason.CALLED_MAJORITY

    # Pos 20: 60% G / 40% A -> Contested, canonical 'N', IUPAC 'R'
    s20 = consensus_sites[19]
    assert s20.canonical_base == "N"
    assert s20.iupac_base == "R"
    assert s20.call_reason == ConsensusCallReason.AMBIGUOUS_IUPAC

    # Pos 40: 18 A, 2 G (VAF A=0.90 >= 0.80) -> Dominant majority, canonical 'A', IUPAC 'A'
    s40 = consensus_sites[39]
    assert s40.canonical_base == "A"
    assert s40.iupac_base == "A"
    assert s40.call_reason == ConsensusCallReason.CALLED_MAJORITY

    # Pos 3107: Legacy rCRS placeholder 'N' -> Masked on both streams
    s3107 = consensus_sites[3106]
    assert s3107.canonical_base == "N"
    assert s3107.iupac_base == "N"
    assert s3107.call_reason == ConsensusCallReason.REFERENCE_N_PLACEHOLDER


def test_component_heteroplasmy_caller(synthetic_environment):
    """Verify Fisher exact strand bias test, floor thresholds, and INSUFFICIENT_DATA gating."""
    fasta_path, bam_path, contig_name = synthetic_environment
    ref_info = ReferenceInfo(contig=contig_name, length=16569, sequence_sha256="dummy")
    policy = ConsensusPolicy()

    engine = PileupEngine(bam_path, fasta_path, policy, contig_name, 16569)
    sites = list(engine.iter_sites())

    # 1. Standard Calling Action
    caller_std = HeteroplasmyCaller(policy, "COMP_TEST", ref_info, GatingAction.STANDARD_CALLING)
    calls_std = caller_std.call_heteroplasmies(sites)
    calls_by_pos = {c.position: c for c in calls_std}

    # Pos 20: Balanced heteroplasmy (4 FWD, 4 REV 'A') -> PASS
    c20 = calls_by_pos[20]
    assert c20.status == CallStatus.PASS
    assert c20.filter_reason == FilterReason.PASS
    assert c20.minor_allele == "A"
    assert c20.strand_bias_pvalue >= policy.strand_bias_alpha

    # Pos 30: Severe strand bias (5 FWD, 0 REV 'C') -> FILTER_STRAND_BIAS
    c30 = calls_by_pos[30]
    assert c30.status == CallStatus.FILTERED
    assert c30.filter_reason == FilterReason.FILTER_STRAND_BIAS
    assert c30.reverse_alt == 0

    # 2. INSUFFICIENT_DATA Gating Action -> Must suppress calling and tag reason explicitly
    caller_insufficient = HeteroplasmyCaller(policy, "COMP_TEST", ref_info, GatingAction.INSUFFICIENT_DATA)
    c20_insufficient = caller_insufficient.evaluate_site(sites[19])
    assert c20_insufficient is not None
    assert c20_insufficient.status == CallStatus.FILTERED
    assert c20_insufficient.filter_reason == FilterReason.FILTER_INSUFFICIENT_DATA


# ============================================================================
# Level 3: End-to-End Pipeline & Orchestrator Integration Tests
# ============================================================================

def test_e2e_orchestrator_full_pipeline(synthetic_environment, tmp_path):
    """Verify input contract resolution, pessimistic gating downgrade, and complete artifact emission."""
    fasta_path, bam_path, contig_name = synthetic_environment
    output_dir = str(tmp_path / "block6_results")
    sample_id = "E2E_SAMPLE_001"

    b4_json = str(tmp_path / "block4_contract.json")
    with open(b4_json, "w", encoding="utf-8") as f:
        json.dump(
            {
                "authentication_status": "AUTHENTIC_ANCIENT_SUPPORTED",
                "rescaling_performed": True,
                "rescaled_bam": bam_path,
                "damage_model_id": "model_ancient_v1",
            },
            f,
        )

    b5_json = str(tmp_path / "block5_contract.json")
    with open(b5_json, "w", encoding="utf-8") as f:
        json.dump(
            {
                "recommended_action": "STANDARD_CALLING",
                "estimation_confidence": "LOW",  # Must trigger downgrade to RESTRICT_HETEROPLASMY
                "contamination_fraction": 0.045,
                "contamination_estimate_source": "SCHMUTZI",
            },
            f,
        )

    orchestrator = ConsensusOrchestrator()

    resolution = orchestrator.resolve_inputs(
        sample_id=sample_id,
        block4_json_path=b4_json,
        block5_json_path=b5_json,
        reference_fasta_path=fasta_path,
        target_contig=contig_name,
    )
    assert resolution.input_role == InputRole.RESCALED_ANCIENT
    assert resolution.gating_action == GatingAction.RESTRICT_HETEROPLASMY
    assert resolution.reference_fasta_path == fasta_path

    qc = orchestrator.run(resolution, output_dir)
    assert qc.sample_id == sample_id
    assert qc.reference_length == 16569

    # Verify physical existence of all 5 primary artifacts
    canonical_fa = os.path.join(output_dir, f"{sample_id}.consensus.fasta")
    iupac_fa = os.path.join(output_dir, f"{sample_id}.consensus_iupac.fasta")
    het_tsv = os.path.join(output_dir, f"{sample_id}.heteroplasmy.tsv")
    vcf_gz = os.path.join(output_dir, f"{sample_id}.variants.vcf.gz")
    qc_json = os.path.join(output_dir, f"{sample_id}.consensus_qc.json")

    assert os.path.exists(canonical_fa)
    assert os.path.exists(iupac_fa)
    assert os.path.exists(het_tsv)
    assert os.path.exists(vcf_gz)
    assert os.path.exists(qc_json)

    # Verify FASTA invariant length
    with pysam.FastaFile(canonical_fa) as fa:
        assert fa.get_reference_length(sample_id) == 16569
    with pysam.FastaFile(iupac_fa) as fa:
        assert fa.get_reference_length(sample_id) == 16569

    # Structural VCF inspection using pysam.VariantFile
    with pysam.VariantFile(vcf_gz) as vcf:
        assert contig_name in vcf.header.contigs
        records = list(vcf.fetch())
        assert len(records) > 0

        pos_map = {rec.pos: rec for rec in records}

        # Verify Pos 20 (Heteroplasmy variant emitted)
        assert 20 in pos_map
        rec_20 = pos_map[20]
        assert "A" in rec_20.alts
        assert "PASS" in rec_20.filter.keys()
        assert rec_20.info["VARTYPE"] == "HETEROPLASMY"

    # Verify QC JSON structure and region metrics
    with open(qc_json, "r", encoding="utf-8") as f:
        qc_data = json.load(f)
        assert qc_data["sample_id"] == sample_id
        assert qc_data["gating_action"] == "RESTRICT_HETEROPLASMY"
        region_names = [r["name"] for r in qc_data["regions"]]
        assert "HVS-I" in region_names
        assert "HVS-II" in region_names


def test_e2e_orchestrator_authenticity_guardrails(tmp_path, synthetic_environment):
    """Verify orchestrator aborts execution on inconclusive or inconsistent authenticity contracts."""
    fasta_path, bam_path, contig_name = synthetic_environment
    orchestrator = ConsensusOrchestrator()

    # Guardrail 1: INCONCLUSIVE status
    b4_inconclusive = str(tmp_path / "b4_inconclusive.json")
    with open(b4_inconclusive, "w", encoding="utf-8") as f:
        json.dump({"authentication_status": "INCONCLUSIVE"}, f)

    b5_dummy = str(tmp_path / "b5_dummy.json")
    with open(b5_dummy, "w", encoding="utf-8") as f:
        json.dump({"recommended_action": "STANDARD_CALLING"}, f)

    with pytest.raises(UnresolvedAuthenticityError):
        orchestrator.resolve_inputs(
            sample_id="SAMPLE_ERR_1",
            block4_json_path=b4_inconclusive,
            block5_json_path=b5_dummy,
            reference_fasta_path=fasta_path,
            target_contig=contig_name,
        )

    # Guardrail 2: AUTHENTIC_ANCIENT_SUPPORTED but rescaling_performed is False
    b4_no_rescale = str(tmp_path / "b4_no_rescale.json")
    with open(b4_no_rescale, "w", encoding="utf-8") as f:
        json.dump(
            {
                "authentication_status": "AUTHENTIC_ANCIENT_SUPPORTED",
                "rescaling_performed": False,
                "rescaled_bam": bam_path,
            },
            f,
        )

    with pytest.raises(ContractValidationError, match="lacks rescaled BAM or rescaling flag"):
        orchestrator.resolve_inputs(
            sample_id="SAMPLE_ERR_2",
            block4_json_path=b4_no_rescale,
            block5_json_path=b5_dummy,
            reference_fasta_path=fasta_path,
            target_contig=contig_name,
        )