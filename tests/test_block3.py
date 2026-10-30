"""Comprehensive stress-testing and quality assurance suite for Block 3 (Alignment, Deduplication, Filtering, and QC)."""

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union
from unittest.mock import MagicMock, patch
import pytest
import pysam

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
    DedupMetrics,
    DedupMode,
    MAPQProfile,
    MultiMappingPolicy,
    MultiMappingStrategy,
    ProvenanceMetrics,
    ReadGroupInfo,
    ReadGroupSource,
    ReferenceInfo,
)
from src.modules.alignment.reference_manager import (
    ContigAliasResolver,
    ReferenceManager,
    ReferenceManagerError,
    ReferenceValidationResult,
)
from src.modules.alignment.bwa_aligner import (
    AlignmentBranchResult,
    AlignmentBranchResults,
    BWAAligner,
    BWAAlignerError,
)
from src.modules.alignment.samtools_pipeline import (
    PipelinePreparationResult,
    SamtoolsPipeline,
    SamtoolsPipelineError,
)
from src.modules.alignment.dedup_engine import (
    DedupEngine,
    DedupEngineError,
    DedupExecutionResult,
)
from src.modules.alignment.filter_engine import (
    AlignmentFilterEngine,
    ExclusionReason,
    FilterEngineError,
    FilterExclusionStats,
    FilterExecutionResult,
)
from src.modules.alignment.metrics_collector import (
    FlagstatMetrics,
    MetricsCollectionResult,
    MetricsCollector,
    MetricsCollectorError,
)
from src.modules.alignment.orchestrator import (
    AlignmentOrchestrator,
    AlignmentOrchestratorError,
)
from src.modules.preprocessing.models import (
    PreprocessingMetrics,
    ProcessedReads,
    SampleIdentity,
    SampleRecord,
)


# ==============================================================================
# FIXTURES & SYNTHETIC DATA GENERATORS
# ==============================================================================

@pytest.fixture
def synthetic_reference(tmp_path: Path) -> ReferenceInfo:
    """Generates a valid single-contig mitochondrial reference FASTA with exact coordinates."""
    fasta_path = tmp_path / "synthetic_rCRS.fasta"
    bases = "ACGT" * 125  # 500 bp cyclic reference
    fasta_content = f">NC_012920.1 Homo sapiens mitochondrion, complete genome\n{bases}\n"
    fasta_path.write_text(fasta_content, encoding="utf-8")

    fai_path = tmp_path / "synthetic_rCRS.fasta.fai"
    fai_path.write_text("NC_012920.1\t500\t56\t500\t501\n", encoding="utf-8")

    for ext in ReferenceManager.BWA_EXTENSIONS:
        (tmp_path / f"synthetic_rCRS.fasta{ext}").write_text("DUMMY_BWA_INDEX", encoding="utf-8")

    file_sha = ReferenceManager.compute_file_sha256(fasta_path)
    seq_md5, seq_id, desc, seq_len, _, _ = ReferenceManager.parse_and_hash_sequence(fasta_path)

    return ReferenceInfo(
        fasta_path=fasta_path,
        sequence_id=seq_id,
        contig_name=seq_id,
        description=desc,
        length_bp=seq_len,
        is_circular=True,
        file_sha256=file_sha,
        sequence_md5=seq_md5,
        has_bwa_index=True,
        has_fai_index=True,
    )


@pytest.fixture
def synthetic_sample(tmp_path: Path) -> Any:
    """Constructs a mock SampleRecord linking all active read branches without constructor coupling."""
    fastq_dir = tmp_path / "fastqs"
    fastq_dir.mkdir(parents=True, exist_ok=True)

    m_fq = fastq_dir / "sample1.collapsed.gz"
    mt_fq = fastq_dir / "sample1.collapsed.truncated.gz"
    s_fq = fastq_dir / "sample1.singleton.truncated.gz"
    r1_fq = fastq_dir / "sample1.pair1.truncated.gz"
    r2_fq = fastq_dir / "sample1.pair2.truncated.gz"

    dummy_content = "@READ1\nACGTACGT\n+\nIIIIIIII\n"
    for p in (m_fq, mt_fq, s_fq, r1_fq, r2_fq):
        p.write_text(dummy_content, encoding="utf-8")

    processed = ProcessedReads(
        merged_fastq=m_fq,
        singleton_fastq=s_fq,
        paired_r1_fastq=r1_fq,
        paired_r2_fastq=r2_fq,
    )
    processed.merged_truncated_fastq = mt_fq

    sample = MagicMock(spec=SampleRecord)
    sample.sample_id = "NA12878"
    sample.identity = SampleIdentity(sample_id="NA12878", library_id="LIB_ANCIENT_01")
    sample.processed_reads = processed
    sample.metrics = PreprocessingMetrics(
        merged_fragments=1000,
        singleton_reads=200,
        surviving_read_pairs=500,
    )
    return sample


def create_mock_bam(
    bam_path: Path,
    reference_name: str,
    reference_length: int,
    records: List[Dict[str, Any]],
    read_groups: Optional[List[Dict[str, str]]] = None,
) -> Path:
    """Helper creating physical BAM files with custom flags, MAPQ scores, and header tags via pysam."""
    bam_path.parent.mkdir(parents=True, exist_ok=True)

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": reference_name, "LN": reference_length}],
    }
    if read_groups:
        header["RG"] = read_groups

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam_out:
        for idx, rec in enumerate(records):
            a = pysam.AlignedSegment(bam_out.header)
            a.query_name = rec.get("query_name", f"READ_{idx}")
            a.query_sequence = rec.get("seq", "ACGTACGT")
            a.flag = rec.get("flag", 0)
            a.reference_id = 0 if not a.is_unmapped else -1
            a.reference_start = rec.get("pos", 0)
            a.mapping_quality = rec.get("mapq", 60)
            a.cigar = rec.get("cigar", ((0, len(a.query_sequence)),))
            a.query_qualities = pysam.qualitystring_to_array("I" * len(a.query_sequence))

            tags = rec.get("tags", [])
            if "rg_id" in rec:
                tags.append(("RG", rec["rg_id"]))
            if tags:
                a.set_tags(tags)

            bam_out.write(a)

    return bam_path


# ==============================================================================
# 1. REFERENCE MANAGER STRESS TESTS
# ==============================================================================

@pytest.mark.unit
class TestReferenceManagerStress:
    """Stress-tests reference FASTA parsing, multi-contig rejection, and cryptographic validation."""

    def test_strict_multi_contig_rejection(self, tmp_path: Path):
        """Asserts that multi-contig FASTA files are immediately rejected in strict mode."""
        multi_fasta = tmp_path / "multi_contig.fasta"
        multi_fasta.write_text(">NC_012920.1 mt\nACGT\n>chr1 nuclear\nTGCA\n", encoding="utf-8")

        _, _, _, _, _, errors = ReferenceManager.parse_and_hash_sequence(multi_fasta)
        assert len(errors) > 0
        assert any("Multiple contigs detected" in e for e in errors)

    def test_invalid_iupac_characters_rejected(self, tmp_path: Path):
        """Asserts that ambiguous or non-canonical nucleotide characters are rejected under strict ACGT rules."""
        corrupt_fasta = tmp_path / "corrupt_bases.fasta"
        corrupt_fasta.write_text(">NC_012920.1\nACGTXRYK\n", encoding="utf-8")

        _, _, _, _, _, errors = ReferenceManager.parse_and_hash_sequence(
            corrupt_fasta, allow_n=False
        )
        assert len(errors) > 0
        assert any("Invalid nucleotide characters" in e for e in errors)

    def test_whitespace_sanitization_and_md5(self, tmp_path: Path):
        """Asserts that internal whitespace and formatting do not alter the biological sequence MD5 hash."""
        f1 = tmp_path / "ref_clean.fasta"
        f2 = tmp_path / "ref_spaces.fasta"

        f1.write_text(">NC_012920.1\nACGTACGT\n", encoding="utf-8")
        f2.write_text(">NC_012920.1\n  ACGT  \n  ACGT \n", encoding="utf-8")

        md5_1, _, _, len_1, _, _ = ReferenceManager.parse_and_hash_sequence(f1)
        md5_2, _, _, len_2, _, _ = ReferenceManager.parse_and_hash_sequence(f2)

        assert md5_1 == md5_2
        assert len_1 == len_2 == 8

    def test_alias_resolution(self):
        """Tests that common nomenclature aliases resolve to canonical RefSeq identifiers."""
        assert ContigAliasResolver.are_aliases("chrM", "NC_012920.1") is True
        assert ContigAliasResolver.are_aliases("MT", "rCRS") is True
        assert ContigAliasResolver.are_aliases("chr1", "NC_012920.1") is False
        assert ContigAliasResolver.resolve_canonical("CHRM") == "NC_012920.1"


# ==============================================================================
# 2. DATA MODELS & READ GROUP SAFETY TESTS
# ==============================================================================

@pytest.mark.unit
class TestModelsAndReadGroupSafety:
    """Stress-tests defensive dataclass validation, header serialization, and collision resistance."""

    def test_read_group_illegal_character_injection(self):
        """Asserts that tab and newline characters in Read Group fields raise ValueError."""
        with pytest.raises(ValueError, match="illegal tab or newline characters"):
            ReadGroupInfo(
                read_group_id="RG1\tID",
                sample_id="SAMPLE1",
                library_id="LIB1",
                platform="ILLUMINA",
            )

    def test_read_group_invalid_platform_rejection(self):
        """Asserts that unsupported sequencing platforms raise ValueError."""
        with pytest.raises(ValueError, match="invalid sequencing platform 'SANGER_UNKNOWN'"):
            ReadGroupInfo(
                read_group_id="RG1",
                sample_id="SAMPLE1",
                library_id="LIB1",
                platform="SANGER_UNKNOWN",
            )

    def test_multi_library_collision_resistant_rg_id(self):
        """Asserts that distinct libraries for the same sample receive unique Read Group IDs."""
        rg1 = BWAAligner._build_read_group(
            sample_id="BSP01",
            library_id="LIB_A",
            source=ReadGroupSource.MERGED,
        )
        rg2 = BWAAligner._build_read_group(
            sample_id="BSP01",
            library_id="LIB_B",
            source=ReadGroupSource.MERGED,
        )
        assert rg1.read_group_id != rg2.read_group_id
        assert "LIB_A" in rg1.read_group_id
        assert "LIB_B" in rg2.read_group_id


# ==============================================================================
# 3. DEDUPLICATION ENGINE TESTS
# ==============================================================================

@pytest.mark.unit
class TestDedupEngineStress:
    """Stress-tests samtools markdup stats parsing, dual-BAM generation, and postcondition assertions."""

    def test_markdup_stats_parsing_mutation_resistant(self):
        """Tests key-value extraction with distinct field counts to prevent false-positive parsing."""
        raw_stats = """
        # Program: samtools markdup
        COMMAND: samtools markdup -s in.bam out.bam
        READ: 10500
        WRITTEN: 10200
        EXAMINED: 9800
        PAIRED: 7600
        SINGLE: 2200
        DUPLICATE PAIR: 400
        DUPLICATE SINGLE: 150
        DUPLICATE TOTAL: 950
        """
        metrics = DedupEngine.parse_markdup_stats(raw_stats)
        assert metrics is not None
        assert metrics.examined_reads == 9800
        assert metrics.examined_pairs == 3800  # 7600 / 2
        assert metrics.examined_single_fragments == 2200
        assert metrics.duplicate_pairs == 400
        assert metrics.duplicate_single_fragments == 150
        assert metrics.duplicate_reads == 950
        assert metrics.retained_reads == 8850  # 9800 - 950
        assert math.isclose(metrics.duplication_rate_reads, 950 / 9800, rel_tol=1e-4)

    def test_markdup_unparseable_stats_returns_none(self):
        """Asserts that corrupted stats content returns None rather than zero-filled objects."""
        corrupt_stats = "RANDOM ERROR: Segmentation fault\n"
        assert DedupEngine.parse_markdup_stats(corrupt_stats) is None

    def test_verify_dedup_bam_duplicate_leak_detection(self, tmp_path: Path):
        """Asserts that verify_dedup_bam raises DedupEngineError if 0x400 duplicate flags are present."""
        leaky_bam = tmp_path / "leaky_dedup.bam"
        records = [
            {"query_name": "READ_CLEAN", "flag": 0, "pos": 10},
            {"query_name": "READ_DUP_LEAK", "flag": 0x400, "pos": 20},
        ]
        create_mock_bam(leaky_bam, "NC_012920.1", 500, records)

        with pytest.raises(DedupEngineError, match="retains duplicate flag 0x400"):
            DedupEngine.verify_dedup_bam(leaky_bam, markdup_record_count=2)


# ==============================================================================
# 4. FILTER ENGINE & EXCLUSION HIERARCHY TESTS
# ==============================================================================

@pytest.mark.unit
class TestFilterEngineStress:
    """Stress-tests exclusive first-failure categorization, MAPQ profiling, and postcondition verification."""

    @pytest.mark.parametrize(
        "records,expected_counts",
        [
            (
                [
                    {"query_name": "R_UNMAPPED", "flag": 0x4, "mapq": 0},
                    {"query_name": "R_SECONDARY", "flag": 0x100, "mapq": 60},
                    {"query_name": "R_SUPP", "flag": 0x800, "mapq": 60},
                    {"query_name": "R_LOW_MAPQ", "flag": 0, "mapq": 15},
                    {"query_name": "R_MULTIMAP_TAG", "flag": 0, "mapq": 35, "tags": [("XT", "R")]},
                    {"query_name": "R_PASS_1", "flag": 0, "mapq": 60},
                    {"query_name": "R_PASS_2", "flag": 0, "mapq": 40},
                ],
                {
                    "total": 7,
                    "unmapped": 1,
                    "secondary": 1,
                    "supplementary": 1,
                    "low_mapq": 1,
                    "multimap_tag": 1,
                    "retained": 2,
                },
            ),
            (
                [
                    {"query_name": "R_SEC_AND_LOW_MAPQ", "flag": 0x100, "mapq": 5},
                    {"query_name": "R_PASS", "flag": 0, "mapq": 60},
                ],
                {
                    "total": 2,
                    "unmapped": 0,
                    "secondary": 1,
                    "supplementary": 0,
                    "low_mapq": 0,
                    "multimap_tag": 0,
                    "retained": 1,
                },
            ),
        ],
    )
    def test_exclusive_first_failure_accounting_sums_to_total(
        self, tmp_path: Path, records: List[Dict[str, Any]], expected_counts: Dict[str, int]
    ):
        """Asserts that every record falls into exactly one exclusive category summing to total input reads."""
        test_bam = tmp_path / f"test_{len(records)}_alignments.bam"
        create_mock_bam(test_bam, "NC_012920.1", 500, records)

        policy = AlignmentFilterPolicy(
            multimapping=MultiMappingPolicy(
                min_mapq=30,
                filter_unmapped=True,
                filter_secondary=True,
                filter_supplementary=True,
                inspect_mapper_tags=True,
            )
        )

        _, ret_prof, stats = AlignmentFilterEngine.compute_mapq_profiles_and_exclusions(
            input_bam=test_bam,
            policy=policy,
        )

        assert stats.total_input_reads == expected_counts["total"]
        assert stats.unmapped_reads_removed == expected_counts["unmapped"]
        assert stats.secondary_reads_removed == expected_counts["secondary"]
        assert stats.supplementary_reads_removed == expected_counts["supplementary"]
        assert stats.low_mapq_reads_removed == expected_counts["low_mapq"]
        assert stats.multimapping_tag_reads_removed == expected_counts["multimap_tag"]
        assert stats.retained_reads == expected_counts["retained"]
        assert stats.total_removed + stats.retained_reads == stats.total_input_reads

    def test_postcondition_violation_detection(self, tmp_path: Path):
        """Asserts that verify_filtered_bam_postconditions raises error if low MAPQ reads leak into filtered BAM."""
        leaky_filtered_bam = tmp_path / "leaky_filtered.bam"
        records = [
            {"query_name": "R_GOOD", "flag": 0, "mapq": 60},
            {"query_name": "R_BAD", "flag": 0, "mapq": 10},
        ]
        create_mock_bam(leaky_filtered_bam, "NC_012920.1", 500, records)

        policy = AlignmentFilterPolicy(multimapping=MultiMappingPolicy(min_mapq=30))

        with pytest.raises(FilterEngineError, match="has MAPQ=10 < 30"):
            AlignmentFilterEngine.verify_filtered_bam_postconditions(leaky_filtered_bam, policy)


# ==============================================================================
# 5. METRICS COLLECTOR & MATHEMATICAL INVARIANTS TESTS
# ==============================================================================

@pytest.mark.unit
class TestMetricsCollectorStress:
    """Stress-tests flagstat parsing, depth coverage extraction, and mathematical invariants."""

    def test_flagstat_invariant_rejection(self):
        """Asserts that mathematically impossible flagstat counts raise MetricsCollectorError."""
        corrupted_flagstat = """
        100 + 0 in total (QC-passed reads + QC-failed reads)
        150 + 0 mapped (150.00% : N/A)
        """
        with pytest.raises(MetricsCollectorError, match="mapped reads .* > total reads"):
            MetricsCollector.parse_flagstat_output(corrupted_flagstat)

    def test_alignment_metrics_invariants_breadth_monotonicity(self):
        """Asserts that non-monotonic breadth percentages trigger invariant failure."""
        invalid_cov = CoverageMetrics(
            breadth_1x_pct=80.0,
            breadth_3x_pct=85.0,  # Impossible: 3x coverage cannot exceed 1x coverage
            breadth_5x_pct=70.0,
            breadth_10x_pct=60.0,
            breadth_30x_pct=50.0,
            max_depth_post_dedup=100,
            p99_depth_post_dedup=90.0,
            p95_depth_post_dedup=80.0,
            median_depth_post_dedup=50.0,
        )
        invalid_metrics = AlignmentMetrics(
            total_input_reads=1000,
            retained_filtered_reads=500,
            coverage=invalid_cov,
        )

        with pytest.raises(MetricsCollectorError, match="Breadth invariant violation"):
            MetricsCollector.validate_alignment_metrics_invariants(invalid_metrics)

    def test_alignment_metrics_invariants_depth_percentiles_order(self):
        """Asserts that out-of-order depth percentiles trigger invariant failure."""
        invalid_cov = CoverageMetrics(
            breadth_1x_pct=90.0,
            breadth_3x_pct=80.0,
            breadth_5x_pct=70.0,
            breadth_10x_pct=60.0,
            breadth_30x_pct=50.0,
            max_depth_post_dedup=100,
            p99_depth_post_dedup=50.0,
            p95_depth_post_dedup=80.0,  # Impossible: P95 cannot exceed P99
            median_depth_post_dedup=20.0,
        )
        invalid_metrics = AlignmentMetrics(
            total_input_reads=1000,
            retained_filtered_reads=500,
            coverage=invalid_cov,
        )

        with pytest.raises(MetricsCollectorError, match="Depth percentiles invariant violation"):
            MetricsCollector.validate_alignment_metrics_invariants(invalid_metrics)


# ==============================================================================
# 6. ORCHESTRATOR INTEGRATION & QC THRESHOLD EVALUATION TESTS
# ==============================================================================

class TestAlignmentOrchestratorIntegration:
    """Stress-tests top-level orchestration, warning code triggers, and fail status determinations."""

    @pytest.mark.unit
    def test_evaluate_qc_thresholds_triggers_warnings(self):
        """Asserts that sub-threshold breadth, depth, and duplication rates assign structured warning codes."""
        result = AlignmentResult(sample_id="SAMPLE_LOW_QC")
        policy = AlignmentFilterPolicy(
            min_breadth_3x_fraction=0.80,         # Requires >= 80% breadth at 3x
            max_duplication_rate_fraction=0.50,   # Requires <= 50% duplication
            min_mean_depth=10.0,                  # Requires >= 10x mean depth
        )

        metrics = AlignmentMetrics(
            total_input_reads=10000,
            retained_filtered_reads=1000,
            raw_mapping_rate_pct=5.0,  # < 10.0% triggers LOW_MAPPING_RATE
            dedup_metrics=DedupMetrics(
                examined_reads=10000,
                duplicate_reads=6000,
                duplication_rate_reads=0.60,  # > 50% duplication triggers HIGH_DUPLICATION_RATE
            ),
            coverage=CoverageMetrics(
                mean_depth_post_dedup=5.0,  # < 10x triggers LOW_DEPTH_COVERAGE
                breadth_3x_pct=65.0,        # < 80% triggers LOW_BREADTH_COVERAGE
            ),
        )

        AlignmentOrchestrator.evaluate_qc_thresholds(metrics, policy, result)

        assert result.status == AlignmentStatus.WARN
        assert AlignmentWarningCode.LOW_BREADTH_COVERAGE in result.warning_codes
        assert AlignmentWarningCode.LOW_DEPTH_COVERAGE in result.warning_codes
        assert AlignmentWarningCode.HIGH_DUPLICATION_RATE in result.warning_codes
        assert AlignmentWarningCode.LOW_MAPPING_RATE in result.warning_codes

    @pytest.mark.unit
    def test_evaluate_qc_thresholds_zero_retained_reads_fails(self):
        """Asserts that zero retained reads sets status strictly to FAIL with NO_MAPPED_MT_READS."""
        result = AlignmentResult(sample_id="SAMPLE_EMPTY")
        policy = AlignmentFilterPolicy()
        metrics = AlignmentMetrics(
            total_input_reads=5000,
            retained_filtered_reads=0,
            coverage=CoverageMetrics(),
        )

        AlignmentOrchestrator.evaluate_qc_thresholds(metrics, policy, result)

        assert result.status == AlignmentStatus.FAIL
        assert AlignmentErrorCode.NO_MAPPED_MT_READS in result.error_codes

    @pytest.mark.integration
    def test_run_sample_alignment_end_to_end_mocked(
        self, tmp_path: Path, synthetic_sample: Any, synthetic_reference: ReferenceInfo
    ):
        """Tests the full orchestration pipeline with mocked binary execution, exercising all stages."""
        out_dir = tmp_path / "results"
        log_dir = tmp_path / "logs"

        raw_bam = tmp_path / "mock_raw.bam"
        sorted_bam = tmp_path / "mock_sorted.bam"
        markdup_bam = tmp_path / "mock_markdup.bam"
        dedup_bam = tmp_path / "mock_dedup.bam"
        filtered_bam = tmp_path / "mock_filtered.bam"

        rg_list = [{"ID": "NA12878_MERGED", "SM": "NA12878", "LB": "LIB_ANCIENT_01", "PL": "ILLUMINA"}]
        records = [{"query_name": "R1", "flag": 0, "pos": 100, "mapq": 60}]

        for b in (raw_bam, sorted_bam, markdup_bam, dedup_bam, filtered_bam):
            create_mock_bam(b, "NC_012920.1", 500, records, read_groups=rg_list)

        rg_info = ReadGroupInfo(
            read_group_id="NA12878_MERGED",
            sample_id="NA12878",
            library_id="LIB_ANCIENT_01",
            platform="ILLUMINA",
        )

        mock_branch_results = AlignmentBranchResults(
            sample_id="NA12878",
            branches={
                ReadGroupSource.MERGED: AlignmentBranchResult(
                    source=ReadGroupSource.MERGED,
                    read_group=rg_info,
                    output_bam=raw_bam,
                    input_units_count=1000,
                    is_successful=True,
                )
            },
        )

        mock_prep_result = PipelinePreparationResult(
            sample_id="NA12878",
            sorted_bam=sorted_bam,
            sorted_bai=sorted_bam.with_suffix(".bam.bai"),
            is_successful=True,
        )

        mock_dedup_result = DedupExecutionResult(
            sample_id="NA12878",
            markdup_bam=markdup_bam,
            dedup_bam=dedup_bam,
            metrics=DedupMetrics(examined_reads=1000, duplicate_reads=50, duplication_rate_reads=0.05),
            is_successful=True,
        )

        mock_filter_result = FilterExecutionResult(
            sample_id="NA12878",
            filtered_dedup_bam=filtered_bam,
            primary_downstream_bam=filtered_bam,
            retained_mapq_profile=MAPQProfile(count_mapq_30_plus=1, mean_mapq=60.0),
            exclusion_stats=FilterExclusionStats(total_input_reads=1, retained_reads=1),
            is_successful=True,
        )

        mock_metrics_result = MetricsCollectionResult(
            sample_id="NA12878",
            alignment_metrics=AlignmentMetrics(
                total_input_reads=1000,
                mapped_reads_raw=1000,
                raw_mapping_rate_pct=100.0,
                retained_filtered_reads=950,
                coverage=CoverageMetrics(
                    mean_depth_post_dedup=25.0,
                    breadth_1x_pct=95.0,
                    breadth_3x_pct=90.0,
                    breadth_5x_pct=85.0,
                    breadth_10x_pct=80.0,
                    breadth_30x_pct=50.0,
                ),
            ),
            is_successful=True,
        )

        with patch.object(BWAAligner, "align_sample_branches", return_value=mock_branch_results), \
             patch.object(SamtoolsPipeline, "process_sample_branches", return_value=mock_prep_result), \
             patch.object(DedupEngine, "execute_dedup", return_value=mock_dedup_result), \
             patch.object(AlignmentFilterEngine, "execute_filtering", return_value=mock_filter_result), \
             patch.object(MetricsCollector, "collect_all_metrics", return_value=mock_metrics_result):

            result = AlignmentOrchestrator.run_sample_alignment(
                sample=synthetic_sample,
                reference_info=synthetic_reference,
                output_dir=out_dir,
                log_dir=log_dir,
                threads=2,
            )

            assert result.status == AlignmentStatus.PASS
            assert result.sample_id == "NA12878"
            assert result.outputs.filtered_bam == filtered_bam
            assert result.metrics.coverage.breadth_3x_pct == 90.0
            assert len(result.errors) == 0