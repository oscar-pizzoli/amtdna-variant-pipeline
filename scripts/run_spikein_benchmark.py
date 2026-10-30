#!/usr/bin/env python3
"""Synthetic Spike-in Benchmark for Block 5 Contamination Quantification.

Creates an isolated contaminated BAM (~8.5% synthetic modern contaminant)
and executes the complete Bayesian schmutzi pipeline to verify sensitivity.
"""

from pathlib import Path
import random
import shutil
import sys
import time

import pysam

# Progetto root setup
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.modules.alignment.models import ReferenceInfo
from src.modules.contamination.models import (
    ContaminationMethod,
    ContaminationPolicy,
    LibraryTreatment,
    SchmutziPolicy,
)
from src.modules.contamination.orchestrator import ContaminationOrchestrator


# Diagnostic SNVs of divergent haplogroup (1-based rCRS positions)
DIVERGENT_SNVS = {
    73: "G", 152: "C", 263: "G", 750: "G", 1438: "G", 2706: "G",
    4769: "G", 7028: "T", 8701: "G", 8860: "G", 9540: "C", 10398: "G",
    10873: "C", 11719: "A", 12705: "T", 14766: "T", 15326: "G", 16223: "T"
}


def create_spiked_bam(
    input_bam_path: Path,
    output_spiked_bam: Path,
    spike_fraction: float = 0.085,
) -> int:
    """Clones a fraction of reads and mutates them to simulate modern contamination."""
    output_spiked_bam.parent.mkdir(parents=True, exist_ok=True)
    random.seed(42)

    total_orig_reads = 0
    spiked_reads_count = 0

    with pysam.AlignmentFile(str(input_bam_path), "rb") as in_bam:
        with pysam.AlignmentFile(str(output_spiked_bam), "wb", template=in_bam) as out_bam:
            # 1. Scrivi tutte le read originali endogene
            for read in in_bam.fetch():
                total_orig_reads += 1
                out_bam.write(read)

                # 2. Clona e muta una percentuale per creare il contaminante sintetico
                if random.random() < spike_fraction:
                    spiked_reads_count += 1
                    contam_read = pysam.AlignedSegment(out_bam.header)
                    contam_read.query_name = f"SYNTH_CONTAM_{read.query_name}_{spiked_reads_count}"
                    contam_read.reference_name = read.reference_name
                    contam_read.reference_start = read.reference_start
                    contam_read.mapping_quality = read.mapping_quality
                    contam_read.cigar = read.cigar
                    contam_read.is_reverse = read.is_reverse

                    # Sostituzione allelica sui siti diagnostici
                    seq_list = list(read.query_sequence)
                    qual_list = list(read.query_qualities)
                    aligned_pairs = read.get_aligned_pairs(matches_only=True)

                    for qpos, rpos in aligned_pairs:
                        rcrs_pos = rpos + 1  # 1-based
                        if rcrs_pos in DIVERGENT_SNVS:
                            seq_list[qpos] = DIVERGENT_SNVS[rcrs_pos]
                            qual_list[qpos] = 37  # High quality modern base

                    contam_read.query_sequence = "".join(seq_list)
                    contam_read.query_qualities = qual_list
                    out_bam.write(contam_read)

    # Crea indice BAI per il BAM contaminato
    pysam.index(str(output_spiked_bam))
    return spiked_reads_count, total_orig_reads

def find_otzi_bam(root_dir: Path) -> Path:
    """Locates OTZI.dedup.bam across project directories."""
    candidates = [
        root_dir / "data" / "processed" / "OTZI.dedup.bam",
        root_dir / "data" / "processed" / "OTZI" / "OTZI.dedup.bam",
        root_dir / "data" / "bam" / "OTZI.dedup.bam",
        root_dir / "results" / "dedup" / "OTZI" / "OTZI.dedup.bam",
        root_dir / "results" / "alignment" / "OTZI" / "OTZI.dedup.bam",
    ]
    for cand in candidates:
        if cand.is_file():
            return cand
    
    # Fallback su ricerca ricorsiva
    matches = list(root_dir.glob("**/OTZI.dedup.bam"))
    if matches:
        return matches[0]
        
    raise FileNotFoundError(f"Impossibile trovare OTZI.dedup.bam sotto {root_dir}")


def main() -> int:
    sample_id = "OTZI_SPIKED"
    orig_bam = find_otzi_bam(PROJECT_ROOT)
    ref_fasta = PROJECT_ROOT / "data" / "reference" / "rCRS.fasta"
    contam_db = PROJECT_ROOT / "data" / "reference" / "schmutzi_db" / "197"
    spiked_bam = PROJECT_ROOT / "data" / "processed" / sample_id / "OTZI_spiked_10pct.bam"

    out_dir = PROJECT_ROOT / "results" / "contamination"
    log_dir = PROJECT_ROOT / "logs" / "contamination"

    print("=" * 70)
    print("STARTING POSITIVE CONTROL BENCHMARK: SYNTHETIC SPIKE-IN TEST")
    print(f"Original BAM: {orig_bam}")
    print(f"Target Spiked BAM: {spiked_bam}")
    print("=" * 70)

    # 1. Genera il BAM sintetico
    print("[1/3] Generating synthetic contaminated BAM...")
    t0_gen = time.perf_counter()
    n_spiked, n_orig = create_spiked_bam(orig_bam, spiked_bam, spike_fraction=0.085)
    injected_pct = (n_spiked / (n_orig + n_spiked)) * 100.0
    print(f"Spike-in complete in {time.perf_counter() - t0_gen:.2f}s:")
    print(f"  - Endogenous reads : {n_orig:,}")
    print(f"  - Contaminant reads: {n_spiked:,}")
    print(f"  - Injected Rate    : {injected_pct:.2f}%")

    # 2. Configura la Reference e Policy
    ref_info = ReferenceInfo(
        fasta_path=ref_fasta,
        contig_name="NC_012920.1",
        length_bp=16569,
    )

    schmutzi_pol = SchmutziPolicy(
        contaminant_database=contam_db,
        threads=4,
        max_iterations=10,
        library_treatment=LibraryTreatment.NONE,
        use_wrapped_reference=False,
    )

    policy = ContaminationPolicy(
        policy_id="spikein_test_policy_v1",
        method=ContaminationMethod.SCHMUTZI,
        enabled=True,
        allow_fallback_to_dissent=False,
        schmutzi=schmutzi_pol,
        low_threshold=0.02,
        moderate_threshold=0.05,
        elevated_threshold=0.10,
    )

    # 3. Esegui la pipeline di contaminazione
    print("\n[2/3] Running Schmutzi Bayesian Estimation on Spiked BAM...")
    result = ContaminationOrchestrator.run_sample_contamination_pipeline(
        sample_id=sample_id,
        reference_info=ref_info,
        input_bam=spiked_bam,
        output_dir=out_dir,
        log_dir=log_dir,
        policy=policy,
        threads=4,
    )

    # 4. Confronto e Verdetto Finale
    print("\n[3/3] BENCHMARK RESULTS COMPARISON")
    print("=" * 70)
    print(f"Injected Ground-Truth Contamination: {injected_pct:.2f}%")
    print(f"Schmutzi Point Estimate            : {result.estimate.point_estimate:.2%}")
    if result.estimate.ci_lower_95 is not None:
        print(f"Schmutzi 95% Credible Interval     : [{result.estimate.ci_lower_95:.2%}, {result.estimate.ci_upper_95:.2%}]")
    print(f"Informative Diagnostic Sites       : {result.diagnostics.informative_positions_evaluated}")
    print(f"Contamination Category             : {result.gating.contamination_category.value}")
    print(f"Downstream Gating Recommendation   : {result.gating.recommended_action.value}")
    print(f"Pipeline Execution Status          : {result.status.value}")
    print("=" * 70)

    return 0

if __name__ == "__main__":
    sys.exit(main())