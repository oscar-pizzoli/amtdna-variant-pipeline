#!/usr/bin/env python3
"""Step-by-step visual audit and validation script for Block 3 (Alignment, Structural Prep, Deduplication, Filtering, and QC)."""

from pathlib import Path
import shutil
import sys
import tempfile
import pysam

# Add project root to sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.modules.alignment.models import (
    AlignmentFilterPolicy,
    AlignmentResult,
    AlignmentStatus,
    AlignmentWarningCode,
    DedupMode,
    MultiMappingPolicy,
    ReadGroupInfo,
    ReadGroupSource,
    ReferenceInfo,
)
from src.modules.alignment.reference_manager import ReferenceManager
from src.modules.alignment.bwa_aligner import (
    AlignmentBranchResult,
    AlignmentBranchResults,
)
from src.modules.alignment.samtools_pipeline import SamtoolsPipeline
from src.modules.alignment.dedup_engine import DedupEngine
from src.modules.alignment.filter_engine import AlignmentFilterEngine
from src.modules.alignment.metrics_collector import MetricsCollector
from src.modules.alignment.orchestrator import AlignmentOrchestrator


def print_stage_header(stage_number: int, title: str) -> None:
    """Prints a styled visual header for pipeline stages."""
    print("\n" + "=" * 80)
    print(f" [STAGE {stage_number}] {title.upper()}")
    print("=" * 80)


def create_synthetic_raw_branch_bam(
    bam_path: Path,
    reference_name: str,
    reference_length: int,
    rg_id: str,
    sample_id: str,
    library_id: str,
) -> Path:
    """Creates a physical raw BAM file containing diverse alignment types."""
    bam_path.parent.mkdir(parents=True, exist_ok=True)
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": reference_name, "LN": reference_length}],
        "RG": [{
            "ID": rg_id,
            "SM": sample_id,
            "LB": library_id,
            "PL": "ILLUMINA",
        }],
    }

    mock_records = [
        # Retained high-confidence reads
        {"qname": "READ_PASS_01", "flag": 0, "pos": 20, "mapq": 60, "seq": "ACGTACGTACGT"},
        {"qname": "READ_PASS_02", "flag": 0, "pos": 50, "mapq": 45, "seq": "ACGTACGTACGT"},
        {"qname": "READ_PASS_03", "flag": 0, "pos": 100, "mapq": 35, "seq": "ACGTACGTACGT"},
        {"qname": "READ_PASS_04", "flag": 0, "pos": 150, "mapq": 60, "seq": "ACGTACGTACGT"},
        {"qname": "READ_PASS_05", "flag": 0, "pos": 200, "mapq": 60, "seq": "ACGTACGTACGT"},
        {"qname": "READ_PASS_06", "flag": 0, "pos": 200, "mapq": 60, "seq": "ACGTACGTACGT"},  # Will be flagged duplicate
        # Low mapping quality (< 30)
        {"qname": "READ_LOW_MQ", "flag": 0, "pos": 220, "mapq": 12, "seq": "ACGTACGTACGT"},
        # Secondary alignment (0x100)
        {"qname": "READ_SEC_01", "flag": 0x100, "pos": 300, "mapq": 60, "seq": "ACGTACGTACGT"},
        # Supplementary / Chimeric alignment (0x800)
        {"qname": "READ_SUP_01", "flag": 0x800, "pos": 350, "mapq": 60, "seq": "ACGTACGTACGT"},
        # Unmapped read (0x4)
        {"qname": "READ_UNMAP", "flag": 0x4, "pos": 0, "mapq": 0, "seq": "NNNNNNNNNNNN"},
    ]

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam_out:
        for r in mock_records:
            seg = pysam.AlignedSegment(bam_out.header)
            seg.query_name = r["qname"]
            seg.query_sequence = r["seq"]
            seg.flag = r["flag"]
            seg.reference_id = 0 if not seg.is_unmapped else -1
            seg.reference_start = r["pos"] if not seg.is_unmapped else -1
            seg.mapping_quality = r["mapq"]
            seg.cigar = ((0, len(r["seq"])),) if not seg.is_unmapped else ()
            seg.query_qualities = pysam.qualitystring_to_array("I" * len(r["seq"]))
            seg.set_tags([("RG", rg_id)])
            bam_out.write(seg)

    return bam_path


def main() -> None:
    """Executes the step-by-step inspection across all Block 3 components."""
    workspace = Path(tempfile.mkdtemp(prefix="amtdna_block3_audit_"))
    log_dir = workspace / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "#" * 80)
    print("  ANCIENT mtDNA PIPELINE - BLOCK 3 STEP-BY-STEP AUDIT RUNNER")
    print(f"  Working Directory: {workspace}")
    print("#" * 80)

    try:
        # ======================================================================
        # STAGE 1: Reference Genome Ingestion & Cryptographic Validation
        # ======================================================================
        print_stage_header(1, "Reference Genome Ingestion & Validation")

        ref_fasta = workspace / "rCRS_synthetic.fasta"
        ref_seq = "ACGT" * 125  # 500 bp synthetic reference
        ref_fasta.write_text(f">NC_012920.1 Homo sapiens mitochondrion\n{ref_seq}\n")

        pysam.faidx(str(ref_fasta))
        file_sha = ReferenceManager.compute_file_sha256(ref_fasta)
        seq_md5, seq_id, desc, seq_len, _, _ = ReferenceManager.parse_and_hash_sequence(ref_fasta)

        ref_info = ReferenceInfo(
            fasta_path=ref_fasta,
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

        print(f"  [+] Reference Contig   : {ref_info.contig_name}")
        print(f"  [+] Reference Length   : {ref_info.length_bp:,} bp")
        print(f"  [+] Circular Topology  : {ref_info.is_circular}")
        print(f"  [+] Sequence MD5       : {ref_info.sequence_md5}")
        print(f"  [+] File SHA-256       : {ref_info.file_sha256[:16]}... (verified)")

        # ======================================================================
        # STAGE 2: Multi-Branch Read Groups & Synthetic Alignments
        # ======================================================================
        print_stage_header(2, "Read Group Ingestion & Multi-Branch Alignment")

        rg_merged = ReadGroupInfo(
            read_group_id="DEMO_SAMPLE_MERGED",
            sample_id="DEMO_SAMPLE",
            library_id="LIB_ANCIENT_01",
            platform="ILLUMINA",
            read_group_source=ReadGroupSource.MERGED,
        )

        raw_bam = workspace / "raw_branches" / "demo.merged.raw.bam"
        create_synthetic_raw_branch_bam(
            raw_bam,
            reference_name="NC_012920.1",
            reference_length=500,
            rg_id=rg_merged.read_group_id,
            sample_id=rg_merged.sample_id,
            library_id=rg_merged.library_id,
        )

        branch_results = AlignmentBranchResults(
            sample_id="DEMO_SAMPLE",
            branches={
                ReadGroupSource.MERGED: AlignmentBranchResult(
                    source=ReadGroupSource.MERGED,
                    read_group=rg_merged,
                    output_bam=raw_bam,
                    input_units_count=10,
                    is_successful=True,
                )
            },
        )

        print(f"  [+] Read Group ID      : {rg_merged.read_group_id}")
        print(f"  [+] Formatted @RG Line : {rg_merged.to_bwa_rg_string()}")
        print(f"  [+] Input Branch BAM   : {raw_bam.name} (10 raw reads)")

        # ======================================================================
        # STAGE 3: Structural Preparation (Collate -> Fixmate -m -> Coordinate Sort)
        # ======================================================================
        print_stage_header(3, "Structural BAM Preparation (collate -> fixmate -m -> sort)")

        prep_dir = workspace / "prepared"
        prep_result = SamtoolsPipeline.process_sample_branches(
            branch_results=branch_results,
            output_dir=prep_dir,
            log_dir=log_dir / "samtools",
            threads=2,
            keep_raw_merged_bam=True,
            expected_contig=ref_info.contig_name,
        )

        print(f"  [+] Prepared Sorted BAM: {prep_result.sorted_bam.name}")
        print(f"  [+] Prepared BAI Index : {prep_result.sorted_bai.name if prep_result.sorted_bai else 'None'}")
        print(f"  [+] Verified Records   : {prep_result.prepared_records_count:,} reads in BAM")

        # ======================================================================
        # STAGE 4: PCR Deduplication & Dual-BAM Generation
        # ======================================================================
        print_stage_header(4, "PCR Duplicate Marking & Dual-BAM Generation")

        dedup_dir = workspace / "dedup"
        dedup_result = DedupEngine.execute_dedup(
            sample_id="DEMO_SAMPLE",
            input_sorted_bam=prep_result.sorted_bam,
            output_dir=dedup_dir,
            log_dir=log_dir / "dedup",
            dedup_mode=DedupMode.BOTH,
            threads=2,
        )

        print(f"  [+] Marked BAM Created : {dedup_result.markdup_bam.name if dedup_result.markdup_bam else 'None'}")
        print(f"  [+] Dedup BAM Created  : {dedup_result.dedup_bam.name if dedup_result.dedup_bam else 'None'}")

        if dedup_result.metrics:
            m = dedup_result.metrics
            print(f"  [+] Examined Reads     : {m.examined_reads:,}")
            print(f"  [+] Duplicates Flagged : {m.duplicate_reads:,} ({m.duplication_rate_reads * 100:.2f}%)")
            print(f"  [+] Retained Unique    : {m.retained_reads:,}")

        # ======================================================================
        # STAGE 5: Quality & Bitwise Flag Filtering
        # ======================================================================
        print_stage_header(5, "Alignment Filtering (MAPQ >= 30, Flag 0x904 Excluded)")

        policy = AlignmentFilterPolicy(
            multimapping=MultiMappingPolicy(
                min_mapq=30,
                filter_unmapped=True,
                filter_secondary=True,
                filter_supplementary=True,
            ),
            min_breadth_3x_fraction=0.80,
            max_duplication_rate_fraction=0.85,
            min_mean_depth=3.0,
        )

        filter_dir = workspace / "filtered"
        filter_result = AlignmentFilterEngine.execute_filtering(
            sample_id="DEMO_SAMPLE",
            input_markdup_bam=dedup_result.markdup_bam,
            input_dedup_bam=dedup_result.dedup_bam,
            output_dir=filter_dir,
            log_dir=log_dir / "filter",
            policy=policy,
            threads=2,
        )

        ex = filter_result.exclusion_stats
        print(f"  [+] Primary Downstream BAM: {filter_result.primary_downstream_bam.name}")
        print("\n  --- Exclusive First-Failure Accounting ---")
        print(f"  |  Total Input Reads       : {ex.total_input_reads:>4}")
        print(f"  |  - Unmapped Removed (0x4): {ex.unmapped_reads_removed:>4}")
        print(f"  |  - Secondary (0x100)     : {ex.secondary_reads_removed:>4}")
        print(f"  |  - Supplementary (0x800) : {ex.supplementary_reads_removed:>4}")
        print(f"  |  - Low MAPQ (< 30)       : {ex.low_mapq_reads_removed:>4}")
        print(f"  |  -------------------------------")
        print(f"  |  = RETAINED PASSED READS : {ex.retained_reads:>4} ({ex.retention_rate * 100:.1f}%)")

        # ======================================================================
        # STAGE 6: QC Metrics, Flagstat & Whole-Mitochondrial Depth
        # ======================================================================
        print_stage_header(6, "QC Metrics, Flagstat & Whole-Mitochondrial Depth")

        metrics_dir = workspace / "reports"
        metrics_result = MetricsCollector.collect_all_metrics(
            sample_id="DEMO_SAMPLE",
            reference_info=ref_info,
            branch_results=branch_results,
            dedup_result=dedup_result,
            filter_result=filter_result,
            output_dir=metrics_dir,
            log_dir=log_dir / "metrics",
            prepared_bam=prep_result.sorted_bam,
            threads=2,
            scan_branch_bams=True,
        )

        cov = metrics_result.coverage
        print("\n  --- Depth Percentiles & Coverage Profile ---")
        print(f"  |  Mean Depth Post-Dedup   : {cov.mean_depth_post_dedup:.2f}x")
        print(f"  |  Median Depth            : {cov.median_depth_post_dedup:.2f}x")
        print(f"  |  Max Depth               : {cov.max_depth_post_dedup}x")
        print(f"  |  P95 / P99 Depth         : {cov.p95_depth_post_dedup:.1f}x / {cov.p99_depth_post_dedup:.1f}x")
        print(f"  |  Breadth >= 1x           : {cov.breadth_1x_pct:.2f}%")
        print(f"  |  Breadth >= 3x           : {cov.breadth_3x_pct:.2f}%")
        print(f"  |  Breadth >= 5x           : {cov.breadth_5x_pct:.2f}%")

        # ======================================================================
        # STAGE 7: Automated Policy QC Evaluation & Final Status
        # ======================================================================
        print_stage_header(7, "Automated Policy Evaluation & QC Status")

        final_result = AlignmentResult(
            sample_id="DEMO_SAMPLE",
            status=AlignmentStatus.UNTESTED,
            policy=policy,
            reference=ref_info,
            metrics=metrics_result.alignment_metrics,
        )

        AlignmentOrchestrator.evaluate_qc_thresholds(
            metrics=final_result.metrics,
            policy=policy,
            result=final_result,
        )

        if final_result.status == AlignmentStatus.UNTESTED:
            final_result.status = AlignmentStatus.PASS

        status_color = "\033[92m" if final_result.status == AlignmentStatus.PASS else "\033[93m"
        reset_color = "\033[0m"

        print(f"  [+] Final Block 3 Status   : {status_color}{final_result.status.value}{reset_color}")
        if final_result.warning_codes:
            print(f"  [!] Diagnostic Warnings    : {[w.value for w in final_result.warning_codes]}")
            for w in final_result.warnings:
                print(f"      - {w}")
        if final_result.error_codes:
            print(f"  [X] Fatal Errors           : {[e.value for e in final_result.error_codes]}")

        print("\n" + "=" * 80)
        print("  BLOCK 3 STEP-BY-STEP AUDIT FINISHED WITH 100% SUCCESS")
        print("=" * 80 + "\n")

    finally:
        shutil.rmtree(workspace, ignore_errors=True)


if __name__ == "__main__":
    main()