"""Track B Engine: Direct allelic dissent and modern mtDNA heterogeneity screening."""

import csv
from dataclasses import dataclass
import logging
from pathlib import Path
import time
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pysam

from src.modules.contamination.models import (
    ContaminationEstimate,
    ContaminationMethod,
    ContaminationOutputFiles,
    ModelQualityDiagnostics,
    WrappingMetadata,
)

logger = logging.getLogger("mtDNA_Pipeline.DissentEngine")

# Standard mitochondrial homopolymer, tandem repeat, and legacy assembly artifact regions
DEFAULT_MT_MASKED_REGIONS: Set[Tuple[int, int]] = {
    (303, 315),      # D-loop Poly-C tract (HVS-II)
    (513, 526),      # D-loop CA/AC repeat region
    (3105, 3108),    # rCRS legacy 3107 'N' phantom deletion site
    (16182, 16194),  # HVS-I Poly-C tract
}

# SAM flag filter mask: UNMAP (0x4) | SECONDARY (0x100) | QCFAIL (0x200) | DUP (0x400) | SUPPLEMENTARY (0x800)
DEFAULT_SAM_FLAG_FILTER: int = 0x4 | 0x100 | 0x200 | 0x400 | 0x800  # 3844


class DissentExecutionError(Exception):
    """Raised when BAM reading, pileup generation, or dissent parsing encounters an unrecoverable failure."""
    pass


@dataclass
class DissentExecutionMetadata:
    """Audit metadata tracking dissent calculation parameters, filters, and execution runtimes."""

    min_base_quality: int = 30
    min_mapping_quality: int = 30
    min_depth_threshold: int = 10
    heterogeneity_vaf_noise_floor: float = 0.015
    is_targeted_scan: bool = False
    effective_min_dissenting_required: int = 3
    sam_flag_filter: int = DEFAULT_SAM_FLAG_FILTER
    total_positions_scanned: int = 0
    masked_positions_skipped: int = 0
    informative_positions_evaluated: int = 0
    dissenting_positions_count: int = 0
    multi_minor_positions_count: int = 0
    runtime_seconds: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        """Serializes dissent execution metadata to a dictionary."""
        return {
            "min_base_quality": self.min_base_quality,
            "min_mapping_quality": self.min_mapping_quality,
            "min_depth_threshold": self.min_depth_threshold,
            "heterogeneity_vaf_noise_floor": self.heterogeneity_vaf_noise_floor,
            "is_targeted_scan": self.is_targeted_scan,
            "effective_min_dissenting_required": self.effective_min_dissenting_required,
            "sam_flag_filter": self.sam_flag_filter,
            "total_positions_scanned": self.total_positions_scanned,
            "masked_positions_skipped": self.masked_positions_skipped,
            "informative_positions_evaluated": self.informative_positions_evaluated,
            "dissenting_positions_count": self.dissenting_positions_count,
            "multi_minor_positions_count": self.multi_minor_positions_count,
            "runtime_seconds": round(self.runtime_seconds, 2),
        }


class DissentEngine:
    """Calculates non-Bayesian allelic dissent and heterogeneity screening metrics across mtDNA sites."""

    @classmethod
    def is_masked_position(
        cls,
        pos_1based: int,
        masked_regions: Set[Tuple[int, int]],
    ) -> bool:
        """Checks if a 1-based genomic position falls within masked homopolymer/artifact regions."""
        for start, end in masked_regions:
            if start <= pos_1based <= end:
                return True
        return False

    @classmethod
    def calculate_effective_dissent_threshold(
        cls,
        num_evaluated: int,
        is_targeted: bool,
        base_min_dissenting: int = 3,
        genome_wide_dissent_rate_floor: float = 0.001,
    ) -> int:
        """Scales the required dissenting site count based on whether the scan is targeted or genome-wide."""
        if is_targeted:
            return max(base_min_dissenting, int(np.ceil(0.05 * num_evaluated))) if num_evaluated > 0 else base_min_dissenting
        return max(base_min_dissenting, int(np.ceil(genome_wide_dissent_rate_floor * num_evaluated)))

    @classmethod
    def evaluate_model_sufficiency(
        cls,
        num_evaluated: int,
        is_targeted: bool,
        target_set_size: Optional[int] = None,
        min_genome_wide_positions: int = 1000,
        min_targeted_positions: int = 10,
    ) -> bool:
        """Evaluates statistical validity and sufficiency of the evaluated position set."""
        if num_evaluated == 0:
            return False

        if is_targeted:
            expected = target_set_size if target_set_size is not None else min_targeted_positions
            required = min(min_targeted_positions, expected)
            return num_evaluated >= required

        return num_evaluated >= min_genome_wide_positions

    @classmethod
    def calculate_robust_point_estimate(
        cls,
        vafs_all: np.ndarray,
        vafs_dissenting: np.ndarray,
        effective_min_dissenting: int,
    ) -> float:
        """Calculates a continuous point estimate without step-function discontinuities."""
        if len(vafs_all) == 0:
            return 0.0

        n_dissent = len(vafs_dissenting)

        if n_dissent == 0:
            return float(np.mean(vafs_all))

        if n_dissent >= effective_min_dissenting:
            return float(np.median(vafs_dissenting))

        weight_dissent = n_dissent / float(effective_min_dissenting)
        base_estimate = float(np.mean(vafs_all))
        sparse_median = float(np.median(vafs_dissenting))

        interpolated_est = ((1.0 - weight_dissent) * base_estimate) + (weight_dissent * sparse_median)
        return float(interpolated_est)

    @classmethod
    def execute_dissent_screening(
        cls,
        sample_id: str,
        input_bam: Path,
        reference_fasta: Path,
        output_dir: Path,
        log_dir: Path,
        expected_contig: Optional[str] = None,
        min_base_qual: int = 30,
        min_mapq: int = 30,
        min_depth_per_site: int = 10,
        noise_vaf_floor: float = 0.015,
        base_min_dissenting_sites: int = 3,
        target_positions: Optional[Set[int]] = None,
        custom_masked_regions: Optional[Set[Tuple[int, int]]] = None,
        sam_flag_filter: int = DEFAULT_SAM_FLAG_FILTER,
    ) -> Tuple[ContaminationEstimate, ModelQualityDiagnostics, ContaminationOutputFiles, DissentExecutionMetadata]:
        """Scans BAM pileup with strict flag filters and homopolymer masking to quantify allelic heterogeneity."""
        if not input_bam.is_file() or input_bam.stat().st_size == 0:
            raise DissentExecutionError(f"Input BAM file not found or empty: {input_bam}")
        if not reference_fasta.is_file():
            raise DissentExecutionError(f"Reference FASTA missing: {reference_fasta}")

        out_sample_dir = output_dir / sample_id
        out_sample_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)

        dissent_tsv_path = out_sample_dir / f"{sample_id}.dissent.tsv"
        masked_regions = custom_masked_regions if custom_masked_regions is not None else DEFAULT_MT_MASKED_REGIONS
        is_targeted = target_positions is not None and len(target_positions) > 0

        logger.info(
            f"[{sample_id}] Starting Track B: Allelic Dissent & Heterogeneity Screening "
            f"(mode={'TARGETED' if is_targeted else 'GENOME_WIDE'}, min_BQ={min_base_qual}, "
            f"min_MAPQ={min_mapq}, min_depth={min_depth_per_site}x, flag_filter=0x{sam_flag_filter:X})..."
        )
        t0 = time.perf_counter()

        minor_vafs: List[float] = []
        informative_depths: List[int] = []
        multi_minor_positions = 0
        total_scanned_pos = 0
        skipped_masked_pos = 0

        try:
            with pysam.AlignmentFile(str(input_bam), "rb") as bam, pysam.FastaFile(str(reference_fasta)) as fa:
                target_contig = (
                    expected_contig if (expected_contig and expected_contig in bam.references) else bam.references[0]
                )
                target_contig_length = bam.get_reference_length(target_contig)

                with open(dissent_tsv_path, "w", newline="", encoding="utf-8") as tsv_f:
                    writer = csv.writer(tsv_f, delimiter="\t")
                    writer.writerow([
                        "pos",
                        "ref_base",
                        "depth_filtered",
                        "major_allele",
                        "major_count",
                        "minor_allele",
                        "minor_count",
                        "minor_vaf",
                        "is_informative",
                        "is_masked",
                    ])

                    # Inclusione di stepper="all" e limiti espliciti del contig per non escludere mate orfani
                    for col in bam.pileup(
                        contig=target_contig,
                        start=0,
                        stop=target_contig_length,
                        truncate=True,
                        stepper="all",
                        min_base_quality=min_base_qual,
                        min_mapping_quality=min_mapq,
                        flag_filter=sam_flag_filter,
                    ):
                        total_scanned_pos += 1
                        pos_1based = col.reference_pos + 1

                        if is_targeted and target_positions is not None and pos_1based not in target_positions:
                            continue

                        is_masked = cls.is_masked_position(pos_1based, masked_regions)
                        if is_masked:
                            skipped_masked_pos += 1
                            continue

                        ref_base = fa.fetch(target_contig, col.reference_pos, col.reference_pos + 1).upper()

                        base_counts = {"A": 0, "C": 0, "G": 0, "T": 0}
                        for read in col.pileups:
                            if not read.is_del and not read.is_refskip:
                                q_seq = read.alignment.query_sequence
                                if q_seq and read.query_position is not None:
                                    q_base = q_seq[read.query_position].upper()
                                    if q_base in base_counts:
                                        base_counts[q_base] += 1

                        depth = sum(base_counts.values())
                        if depth < min_depth_per_site:
                            writer.writerow([pos_1based, ref_base, depth, "N", 0, "N", 0, 0.0, False, is_masked])
                            continue

                        sorted_alleles = sorted(base_counts.items(), key=lambda x: x[1], reverse=True)
                        major_base, major_cnt = sorted_alleles[0]
                        minor_base, minor_cnt = sorted_alleles[1]

                        minor_vaf = float(minor_cnt / depth) if depth > 0 else 0.0
                        is_informative = minor_vaf >= noise_vaf_floor

                        if len(sorted_alleles) > 2 and (sorted_alleles[2][1] / depth) >= noise_vaf_floor:
                            multi_minor_positions += 1

                        minor_vafs.append(minor_vaf)
                        informative_depths.append(depth)

                        writer.writerow([
                            pos_1based,
                            ref_base,
                            depth,
                            major_base,
                            major_cnt,
                            minor_base,
                            minor_cnt,
                            round(minor_vaf, 5),
                            is_informative,
                            is_masked,
                        ])

        except Exception as exc:
            raise DissentExecutionError(f"Pileup profiling failed for {sample_id}: {exc}") from exc

        elapsed = time.perf_counter() - t0
        num_evaluated = len(informative_depths)

        if num_evaluated == 0:
            logger.warning(
                f"[{sample_id}] Zero positions satisfied minimum depth >= {min_depth_per_site}x "
                f"(masked positions skipped: {skipped_masked_pos})."
            )
            point_est = 0.0
            heterogeneity_frac = 0.0
            mean_cov = 0.0
            n_dissenting = 0
            effective_min_dissent = base_min_dissenting_sites
            has_multiple_contaminants = False
            model_converged = False
        else:
            vafs_arr = np.array(minor_vafs)
            mean_cov = float(np.mean(informative_depths))
            heterogeneity_frac = float(np.percentile(vafs_arr, 95))

            dissenting_arr = vafs_arr[vafs_arr >= noise_vaf_floor]
            n_dissenting = len(dissenting_arr)

            effective_min_dissent = cls.calculate_effective_dissent_threshold(
                num_evaluated=num_evaluated,
                is_targeted=is_targeted,
                base_min_dissenting=base_min_dissenting_sites,
            )

            point_est = cls.calculate_robust_point_estimate(
                vafs_all=vafs_arr,
                vafs_dissenting=dissenting_arr,
                effective_min_dissenting=effective_min_dissent,
            )

            multi_minor_fraction = multi_minor_positions / float(num_evaluated)
            has_multiple_contaminants = (multi_minor_fraction >= 0.01) and (multi_minor_positions >= 3)

            model_converged = cls.evaluate_model_sufficiency(
                num_evaluated=num_evaluated,
                is_targeted=is_targeted,
                target_set_size=len(target_positions) if target_positions else None,
            )

        estimate = ContaminationEstimate(
            point_estimate=point_est,
            ci_lower_95=None,
            ci_upper_95=None,
            method_used=ContaminationMethod.DISSENT_SCREENING,
            cont_deam_prior_estimate=None,
            log_likelihood=None,
            heterogeneity_fraction=heterogeneity_frac,
        )

        diagnostics = ModelQualityDiagnostics(
            mean_depth=mean_cov,
            informative_positions_evaluated=num_evaluated,
            ci_width=None,
            model_converged=model_converged,
            iterations_run=0,
            warn_contaminant_damage_pattern=False,
            possible_multiple_contaminants=has_multiple_contaminants,
            contaminant_haplogroups=[],
            wrapping=WrappingMetadata(),
        )

        outputs = ContaminationOutputFiles(
            schmutzi_ready_bam=None,
            schmutzi_ready_bai=None,
            dissent_tsv=dissent_tsv_path if dissent_tsv_path.is_file() else None,
        )

        metadata = DissentExecutionMetadata(
            min_base_quality=min_base_qual,
            min_mapping_quality=min_mapq,
            min_depth_threshold=min_depth_per_site,
            heterogeneity_vaf_noise_floor=noise_vaf_floor,
            is_targeted_scan=is_targeted,
            effective_min_dissenting_required=effective_min_dissent,
            sam_flag_filter=sam_flag_filter,
            total_positions_scanned=total_scanned_pos,
            masked_positions_skipped=skipped_masked_pos,
            informative_positions_evaluated=num_evaluated,
            dissenting_positions_count=n_dissenting,
            multi_minor_positions_count=multi_minor_positions,
            runtime_seconds=elapsed,
        )

        logger.info(
            f"[{sample_id}] Dissent Screening finished in {elapsed:.2f}s | "
            f"Evaluated: {num_evaluated:,} sites (mean depth={mean_cov:.1f}x, masked={skipped_masked_pos}, "
            f"mode={'TARGETED' if is_targeted else 'GENOME_WIDE'}) | "
            f"Dissenting sites: {n_dissenting}/{effective_min_dissent} req | Point Est: {point_est:.2%} | "
            f"Heterogeneity (p95): {heterogeneity_frac:.2%} | Sufficiency: {model_converged}"
        )

        return estimate, diagnostics, outputs, metadata

    @classmethod
    def run(
        cls,
        sample_id: str,
        input_bam: Path,
        reference_fasta: Path,
        output_dir: Path,
        policy: Optional[Any] = None,
        expected_contig: Optional[str] = None,
        log_dir: Optional[Path] = None,
        **kwargs: Any,
    ) -> Tuple[ContaminationEstimate, ModelQualityDiagnostics, Optional[Path]]:
        """Adapter matching the orchestrator interface: executes screening and returns (estimate, diagnostics, tsv_file)."""
        actual_log_dir = log_dir or (output_dir / sample_id / "logs")

        min_bq = getattr(policy, "min_base_quality", 30) if policy else 30
        min_mapq = getattr(policy, "min_mapping_quality", 30) if policy else 30
        min_depth = getattr(policy, "min_depth_per_site", 10) if policy else 10
        noise_floor = getattr(policy, "noise_vaf_floor", 0.015) if policy else 0.015
        base_min_dissent = getattr(policy, "base_min_dissenting_sites", 3) if policy else 3
        target_positions = getattr(policy, "target_positions", None) if policy else None

        estimate, diagnostics, outputs, _ = cls.execute_dissent_screening(
            sample_id=sample_id,
            input_bam=input_bam,
            reference_fasta=reference_fasta,
            output_dir=output_dir,
            log_dir=actual_log_dir,
            expected_contig=expected_contig,
            min_base_qual=min_bq,
            min_mapq=min_mapq,
            min_depth_per_site=min_depth,
            noise_vaf_floor=noise_floor,
            base_min_dissenting_sites=base_min_dissent,
            target_positions=target_positions,
        )

        return estimate, diagnostics, outputs.dissent_tsv