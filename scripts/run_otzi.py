#!/usr/bin/env python3
"""Run complete pipeline for sample: Preprocessing -> Alignment -> BAM Filtering."""

import logging
from pathlib import Path
import subprocess
import sys
import pysam

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.modules.alignment.models import (
    AlignmentFilterPolicy,
    DedupMode,
    ReferenceInfo,
)
from src.modules.alignment.reference_manager import ReferenceManager
from src.modules.alignment.orchestrator import AlignmentOrchestrator
from src.modules.preprocessing.orchestrator import (
    PreprocessingOrchestrator,
    RawSampleInput,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("mtDNA_Pipeline.Runner")


def prepare_reference(fasta_path: Path) -> ReferenceInfo:
    """Validates FASTA reference, creates .fai and BWA indexes if missing, and returns ReferenceInfo."""
    ref_path = Path(fasta_path).resolve()
    if not ref_path.is_file():
        raise FileNotFoundError(f"Reference FASTA not found: {ref_path}")

    # 1. Indicizzazione FAI (samtools / pysam)
    fai_path = Path(f"{ref_path}.fai")
    if not fai_path.is_file():
        logger.info(f"Generating FASTA index (.fai) for {ref_path.name}...")
        pysam.faidx(str(ref_path))

    # 2. Indicizzazione BWA (.bwt, .pac, .ann, .amb, .sa)
    bwa_exts = [".bwt", ".pac", ".ann", ".amb", ".sa"]
    has_bwa = all(Path(f"{ref_path}{ext}").is_file() for ext in bwa_exts)
    if not has_bwa:
        logger.info(f"Generating BWA index for {ref_path.name}...")
        res = subprocess.run(["bwa", "index", str(ref_path)], capture_output=True, text=True, check=False)
        if res.returncode != 0:
            raise RuntimeError(f"bwa index failed: {res.stderr}")

    # 3. Parsing e hashing crittografico (allow_n=True per tollerare il placeholder N a pos 3107)
    file_sha = ReferenceManager.compute_file_sha256(ref_path)
    seq_md5, seq_id, desc, seq_len, _, errors = ReferenceManager.parse_and_hash_sequence(
        ref_path,
        allow_n=True,
    )

    if errors:
        raise ValueError(f"Reference validation errors: {'; '.join(errors)}")

    return ReferenceInfo(
        fasta_path=ref_path,
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


def main() -> None:
    sample_id = "OTZI"
    r1_path = Path("data/raw_fastq/OTZI_R1.fastq.gz")
    r2_path = Path("data/raw_fastq/OTZI_R2.fastq.gz")
    ref_fasta = Path("data/reference/rCRS.fasta")

    out_dir = Path("results")
    log_dir = Path("logs")

    # Step 1: Validazione e indicizzazione del genoma di riferimento rCRS
    logger.info("--- [1/3] Validating and indexing rCRS reference ---")
    ref_info = prepare_reference(ref_fasta)
    logger.info(f"Reference loaded: {ref_info.contig_name} ({ref_info.length_bp} bp)")

    # Step 2: Block 2 - Preprocessing e collasso letture sovrapposte (AdapterRemoval)
    logger.info("--- [2/3] Running Block 2: AdapterRemoval & Read Merging ---")
    raw_input = RawSampleInput(
        sample_id=sample_id,
        r1_path=r1_path,
        r2_path=r2_path,
    )
    sample_record = PreprocessingOrchestrator.run_sample(
        sample_input=raw_input,
        output_dir=out_dir / "preprocessed",
        log_dir=log_dir / "preprocessing",
    )

    # Step 3: Block 3 - BWA aDNA alignment, fixmate, markdup e quality filtering
    logger.info("--- [3/3] Running Block 3: Alignment, Deduplication & Filtering ---")
    policy = AlignmentFilterPolicy(
        dedup_mode=DedupMode.BOTH,
        min_breadth_3x_fraction=0.80,
    )

    result = AlignmentOrchestrator.run_sample_alignment(
        sample=sample_record,
        reference_info=ref_info,
        output_dir=out_dir / "alignment",
        log_dir=log_dir / "alignment",
        policy=policy,
        threads=4,
    )

    logger.info("=" * 65)
    logger.info(f"Execution completed with Status: {result.status.value}")
    logger.info(f"Filtered Dedup BAM     : {result.outputs.dedup_bam}")
    logger.info(f"Filtered Markdup BAM   : {result.outputs.markdup_bam}")
    logger.info(f"Mean Coverage Depth    : {result.metrics.coverage.mean_depth_post_dedup:.2f}x")
    logger.info(f"Breadth (>=3x coverage): {result.metrics.coverage.breadth_3x_pct:.2f}%")
    logger.info(f"Retained mtDNA Reads   : {result.metrics.retained_filtered_reads:,}")
    logger.info("=" * 65)


if __name__ == "__main__":
    main()