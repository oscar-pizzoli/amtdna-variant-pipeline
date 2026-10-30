#!/usr/bin/env python3
"""Run Block 4 (aDNA Damage Profiling, Quality Rescaling & Authenticity QC) on sample OTZI."""

import logging
from pathlib import Path
import subprocess
import sys
import pysam

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.modules.alignment.models import ReferenceInfo
from src.modules.alignment.reference_manager import ReferenceManager
from src.modules.damage.models import (
    DamageInputBAMType,
    DamagePolicy,
    LibraryTreatment,
)
from src.modules.damage.orchestrator import DamageOrchestrator

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("mtDNA_Pipeline.Block4_Runner")


def load_reference_info(fasta_path: Path) -> ReferenceInfo:
    """Loads and validates reference FASTA metadata and indices."""
    ref_path = Path(fasta_path).resolve()
    if not ref_path.is_file():
        raise FileNotFoundError(f"Reference FASTA not found: {ref_path}")

    # Ensure .fai index exists
    fai_path = Path(f"{ref_path}.fai")
    if not fai_path.is_file():
        logger.info(f"Generating FASTA index (.fai) for {ref_path.name}...")
        pysam.faidx(str(ref_path))

    # Parse and compute sequence hashes
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
    ref_fasta = Path("data/reference/rCRS.fasta")
    out_dir = Path("results")
    log_dir = Path("logs")

    # 1. Verifica la presenza dei file BAM generati a monte (Blocco 3)
    dedup_bam = out_dir / "alignment" / sample_id / f"{sample_id}.dedup.bam"
    markdup_bam = out_dir / "alignment" / sample_id / f"{sample_id}.markdup.bam"

    if not dedup_bam.is_file():
        logger.error(
            f"BAM deduplicato non trovato in: {dedup_bam}\n"
            "Esegui prima 'python scripts/run_otzi.py' per completare i Blocchi 1-3."
        )
        sys.exit(1)

    # 2. Carica i metadati del riferimento genomico (rCRS)
    logger.info("--- [1/2] Loading reference genome metadata ---")
    ref_info = load_reference_info(ref_fasta)
    logger.info(f"Reference: {ref_info.contig_name} ({ref_info.length_bp} bp)")

    # 3. Configura la policy per il DNA antico di Ötzi (libreria standard Non-UDG)
    damage_policy = DamagePolicy(
        input_bam_type=DamageInputBAMType.FILTERED_DEDUP,
        library_treatment=LibraryTreatment.NONE,
        enable_rescaling=True,
        enable_pmd_scoring=True,
        enable_pmd_filtering=False,
        terminal_positions_to_profile=25,
        min_c_to_t_5p_strong_threshold=0.10,
        min_g_to_a_3p_strong_threshold=0.10,
        min_c_to_t_5p_moderate_threshold=0.03,
        min_g_to_a_3p_moderate_threshold=0.03,
        min_c_to_t_5p_observations=50,
        min_g_to_a_3p_observations=50,
        max_median_length_threshold_bp=80,
        min_reads_for_damage_model=100,
    )

    # 4. Esegui il Blocco 4
    logger.info("--- [2/2] Running Block 4: Damage Profiling, Rescaling & Authentication ---")
    damage_result = DamageOrchestrator.run_sample_damage_pipeline(
        sample_id=sample_id,
        reference_info=ref_info,
        output_dir=out_dir / "damage",
        log_dir=log_dir / "damage",
        policy=damage_policy,
        dedup_bam=dedup_bam,
        markdup_bam=markdup_bam,
        threads=4,
    )

    # 5. Riepilogo finale
    logger.info("=" * 65)
    logger.info(f"EXECUTION SUMMARY | Status: {damage_result.status.value}")
    logger.info(f"Damage Signal Level  : {damage_result.profile.damage_signal_level.value}")
    logger.info(f"Authenticity Verdict : {damage_result.profile.authentication_verdict.value}")
    logger.info(
        f"5' C->T (pos 1)      : {damage_result.profile.terminal_damage.c_to_t_5p_pos1:.2%} "
        f"({damage_result.profile.terminal_damage.n_c_to_t_5p_pos1_obs}/{damage_result.profile.terminal_damage.n_c_ref_5p_pos1_total} obs)"
    )
    logger.info(
        f"3' G->A (pos 1)      : {damage_result.profile.terminal_damage.g_to_a_3p_pos1:.2%} "
        f"({damage_result.profile.terminal_damage.n_g_to_a_3p_pos1_obs}/{damage_result.profile.terminal_damage.n_g_ref_3p_pos1_total} obs)"
    )
    logger.info(f"Median Fragment Len  : {damage_result.profile.fragment_length.median_length:.1f} bp")
    if damage_result.outputs.rescaled_bam:
        logger.info(f"Rescaled BAM         : {damage_result.outputs.rescaled_bam}")
    if damage_result.outputs.damage_model_json:
        logger.info(f"Damage Model JSON    : {damage_result.outputs.damage_model_json}")
    logger.info("=" * 65)


if __name__ == "__main__":
    main()