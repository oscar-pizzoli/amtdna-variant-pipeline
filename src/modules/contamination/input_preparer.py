"""Block 5 Input & Reference Preparation Engines for schmutzi.

Provides two specialized adapters:
1. SchmutziReferenceAdapter: Handles mitochondrial circular reference wrapping,
   flank padding, indexing, and coordinate metadata tracking.
2. SchmutziInputPreparer: Strictly validates mtDNA-only BAM structure, executes
   conditional depth subsampling (target 40x), runs mandatory atomic 'samtools calmd'
   MD tag injection, manages BAI indexing, and verifies integrity.
"""

import hashlib
import logging
from pathlib import Path
import shutil
import subprocess
import time
from typing import Optional, Tuple
import pysam

from src.modules.alignment.models import ReferenceInfo
from src.modules.contamination.models import SubsamplingPolicy, WrappingMetadata

logger = logging.getLogger("mtDNA_Pipeline.ContaminationInputPreparer")


class SchmutziPreparationError(Exception):
    """Raised when BAM validation, subsampling, calmd MD tag injection, or reference wrapping fails."""
    pass


class SchmutziReferenceAdapter:
    """Adapts mitochondrial reference genomes for circularity awareness in schmutzi."""

    @classmethod
    def compute_sha256(cls, file_path: Path) -> str:
        """Computes SHA-256 digest of a target file."""
        hasher = hashlib.sha256()
        with open(file_path, "rb") as f:
            while chunk := f.read(65536):
                hasher.update(chunk)
        return hasher.hexdigest()

    @classmethod
    def prepare_wrapped_reference(
        cls,
        reference_info: ReferenceInfo,
        work_dir: Path,
        padding_bp: int = 200,
    ) -> Tuple[Path, WrappingMetadata]:
        """Generates a wrapped (circularly extended) FASTA reference with flanking padding."""
        if not reference_info.fasta_path.is_file():
            raise SchmutziPreparationError(f"Canonical reference FASTA not found: {reference_info.fasta_path}")

        work_dir.mkdir(parents=True, exist_ok=True)
        wrapped_fa_path = work_dir / f"{reference_info.contig_name}_wrapped_{padding_bp}bp.fasta"
        fai_path = Path(f"{wrapped_fa_path}.fai")

        try:
            with pysam.FastaFile(str(reference_info.fasta_path)) as fa:
                if reference_info.contig_name not in fa.references:
                    raise SchmutziPreparationError(
                        f"Contig '{reference_info.contig_name}' missing from FASTA references: {fa.references}"
                    )
                seq = fa.fetch(reference_info.contig_name)

            orig_len = len(seq)
            if orig_len != reference_info.length_bp:
                raise SchmutziPreparationError(
                    f"Declared reference length ({reference_info.length_bp} bp) differs from actual FASTA length ({orig_len} bp)."
                )

            padding_seq = seq[:padding_bp]
            wrapped_seq = seq + padding_seq
            ext_len = len(wrapped_seq)

            tmp_fa = wrapped_fa_path.with_suffix(".tmp.fasta")
            with open(tmp_fa, "w", encoding="utf-8") as out_f:
                out_f.write(f">{reference_info.contig_name}\n")
                for i in range(0, ext_len, 60):
                    out_f.write(f"{wrapped_seq[i:i+60]}\n")
            tmp_fa.replace(wrapped_fa_path)

            if not fai_path.is_file() or fai_path.stat().st_mtime < wrapped_fa_path.stat().st_mtime:
                pysam.faidx(str(wrapped_fa_path))

            wrapping_meta = WrappingMetadata(
                applied=True,
                padding_bp=padding_bp,
                wrapped_contig_name=reference_info.contig_name,
                original_length_bp=orig_len,
                extended_length_bp=ext_len,
                remapping_offset=padding_bp,
            )

            logger.info(
                f"Generated circular wrapped reference -> {wrapped_fa_path.name} "
                f"({orig_len} bp + {padding_bp} bp padding = {ext_len} bp)."
            )
            return wrapped_fa_path, wrapping_meta

        except Exception as exc:
            if isinstance(exc, SchmutziPreparationError):
                raise
            raise SchmutziPreparationError(f"Failed to generate wrapped reference: {exc}") from exc


class SchmutziInputPreparer:
    """Prepares and strictly validates BAM files to satisfy schmutzi's input requirements."""

    @classmethod
    def run_samtools_quickcheck(cls, bam_path: Path) -> bool:
        """Runs 'samtools quickcheck' to verify file completeness and intact EOF blocks."""
        samtools_bin = shutil.which("samtools")
        if not samtools_bin:
            return bam_path.is_file() and bam_path.stat().st_size > 0

        res = subprocess.run([samtools_bin, "quickcheck", str(bam_path)], check=False)
        return res.returncode == 0

    @classmethod
    def ensure_bai_index(cls, bam_path: Path, reference_contig: Optional[str] = None, threads: int = 4) -> Path:
        """Ensures a valid, queryable BAI index exists."""
        bai_path = Path(f"{bam_path}.bai")
        needs_index = not bai_path.is_file() or bai_path.stat().st_mtime < bam_path.stat().st_mtime

        if needs_index:
            try:
                samtools_bin = shutil.which("samtools")
                if samtools_bin:
                    subprocess.run(
                        [samtools_bin, "index", "-@", str(threads), str(bam_path)],
                        capture_output=True,
                        check=True,
                    )
                else:
                    pysam.index(str(bam_path))
            except Exception as exc:
                raise SchmutziPreparationError(f"Failed indexing BAM '{bam_path.name}': {exc}") from exc

        if reference_contig:
            try:
                with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                    itr = bam.fetch(contig=reference_contig)
                    next(itr, None)
            except Exception as exc:
                raise SchmutziPreparationError(f"BAM index query verification failed for {bam_path.name}: {exc}") from exc

        return bai_path

    @classmethod
    def compute_contig_mean_depth(cls, bam_path: Path, contig_name: str, contig_length: int) -> float:
        """Calculates exact mean sequencing depth across target contig without dropping orphan mates."""
        if not bam_path.is_file():
            raise SchmutziPreparationError(f"BAM file not found: {bam_path}")

        cls.ensure_bai_index(bam_path)

        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                if contig_name not in bam.references:
                    raise SchmutziPreparationError(f"Contig '{contig_name}' missing from BAM references: {bam.references}")

                total_bases = sum(
                    col.nsegments
                    for col in bam.pileup(
                        contig=contig_name,
                        start=0,
                        stop=contig_length,
                        truncate=True,
                        stepper="all",  # Obbligatorio per non escludere read con mate nucleare/orfano (FLAG 129)
                    )
                )
                return float(total_bases / contig_length) if contig_length > 0 else 0.0
        except Exception as exc:
            raise SchmutziPreparationError(f"Failed to compute BAM depth on {contig_name}: {exc}") from exc

    @classmethod
    def subsample_bam(
        cls,
        input_bam: Path,
        output_bam: Path,
        target_fraction: float,
        seed: int = 42,
        threads: int = 4,
    ) -> None:
        """Executes reproducible pseudo-random subsampling using samtools view -s."""
        samtools_bin = shutil.which("samtools")
        if not samtools_bin:
            raise SchmutziPreparationError("samtools executable not found in PATH.")

        subsample_param = float(seed) + float(target_fraction)
        cmd = [
            samtools_bin,
            "view",
            "-@", str(threads),
            "-s", f"{subsample_param:.5f}",
            "-b",
            str(input_bam),
            "-o", str(output_bam),
        ]

        logger.info(f"Subsampling BAM: fraction={target_fraction:.4f} (seed={seed}) -> {output_bam.name}")
        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if res.returncode != 0:
            raise SchmutziPreparationError(f"samtools subsampling failed: {res.stderr.strip()}")

        cls.ensure_bai_index(output_bam, threads=threads)

    @classmethod
    def validate_bam_structure(
        cls,
        bam_path: Path,
        reference_info: ReferenceInfo,
        allow_non_mt_contigs: bool = False,
    ) -> None:
        """Validates BAM header, @SQ sequence length, and mtDNA-only contig constraints."""
        if not cls.run_samtools_quickcheck(bam_path):
            raise SchmutziPreparationError(f"BAM quickcheck failed or file is truncated: {bam_path}")

        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                if reference_info.contig_name not in bam.references:
                    raise SchmutziPreparationError(
                        f"Target contig '{reference_info.contig_name}' missing from BAM header. "
                        f"Available contigs: {bam.references}"
                    )

                bam_contig_len = bam.get_reference_length(reference_info.contig_name)
                if bam_contig_len != reference_info.length_bp:
                    raise SchmutziPreparationError(
                        f"Length mismatch for contig '{reference_info.contig_name}': "
                        f"BAM @SQ={bam_contig_len} bp, Reference FASTA={reference_info.length_bp} bp."
                    )

                if not allow_non_mt_contigs and len(bam.references) > 1:
                    raise SchmutziPreparationError(
                        f"schmutzi requires an mtDNA-only BAM. Found {len(bam.references)} contigs: {bam.references}"
                    )
        except Exception as exc:
            if isinstance(exc, SchmutziPreparationError):
                raise
            raise SchmutziPreparationError(f"BAM structure validation failed on {bam_path.name}: {exc}") from exc

    @classmethod
    def run_calmd_atomic(
        cls,
        input_bam: Path,
        reference_fasta: Path,
        output_bam: Path,
        threads: int = 4,
    ) -> Path:
        """Executes 'samtools calmd -b' writing atomically to guarantee valid MD and NM tags."""
        samtools_bin = shutil.which("samtools")
        if not samtools_bin:
            raise SchmutziPreparationError("samtools binary not found in system PATH.")

        tmp_bam = output_bam.with_suffix(".tmp.bam")
        cmd = [
            samtools_bin,
            "calmd",
            "-b",
            "-@", str(threads),
            str(input_bam),
            str(reference_fasta),
        ]

        logger.info(f"Injecting MD tags via samtools calmd -> {output_bam.name} (threads={threads})...")
        t0 = time.perf_counter()

        try:
            with open(tmp_bam, "wb") as out_f:
                res = subprocess.run(
                    cmd,
                    stdout=out_f,
                    stderr=subprocess.PIPE,
                    check=False,
                )

            if res.returncode != 0:
                err_msg = res.stderr.decode("utf-8", errors="replace").strip()
                if tmp_bam.exists():
                    tmp_bam.unlink()
                raise SchmutziPreparationError(f"samtools calmd failed (exit code {res.returncode}): {err_msg}")

            if not tmp_bam.exists() or tmp_bam.stat().st_size == 0:
                raise SchmutziPreparationError(f"samtools calmd produced an empty file: {tmp_bam}")

            if not cls.run_samtools_quickcheck(tmp_bam):
                if tmp_bam.exists():
                    tmp_bam.unlink()
                raise SchmutziPreparationError(f"Generated BAM failed quickcheck validation: {tmp_bam}")

            tmp_bam.replace(output_bam)

        except OSError as exc:
            if tmp_bam.exists():
                tmp_bam.unlink()
            raise SchmutziPreparationError(f"Failed launching samtools calmd process: {exc}") from exc

        elapsed = time.perf_counter() - t0
        logger.info(
            f"MD tags generated successfully in {elapsed:.2f}s "
            f"({output_bam.stat().st_size / (1024*1024):.2f} MB)."
        )
        return output_bam

    @classmethod
    def prepare_bam_for_schmutzi(
        cls,
        sample_id: str,
        input_bam: Path,
        reference_info: ReferenceInfo,
        work_dir: Path,
        subsampling_policy: Optional[SubsamplingPolicy] = None,
        threads: int = 4,
        allow_non_mt_contigs: bool = False,
    ) -> Tuple[Path, Path, bool, float, Optional[float]]:
        """Prepares and validates BAM strictly for schmutzi."""
        if not input_bam.is_file() or input_bam.stat().st_size == 0:
            raise SchmutziPreparationError(f"Input BAM not found or empty: {input_bam}")

        cls.validate_bam_structure(
            bam_path=input_bam,
            reference_info=reference_info,
            allow_non_mt_contigs=allow_non_mt_contigs,
        )

        work_dir.mkdir(parents=True, exist_ok=True)

        orig_depth = cls.compute_contig_mean_depth(
            bam_path=input_bam,
            contig_name=reference_info.contig_name,
            contig_length=reference_info.length_bp,
        )
        logger.info(f"[{sample_id}] Initial mtDNA coverage depth: {orig_depth:.2f}x")

        source_bam_for_calmd = input_bam
        was_subsampled = False
        post_depth = orig_depth

        sub_pol = subsampling_policy or SubsamplingPolicy(enabled=False)
        if sub_pol.enabled and orig_depth > sub_pol.min_depth_threshold:
            sampling_fraction = min(1.0, sub_pol.target_depth / orig_depth)
            subsampled_bam = work_dir / f"{sample_id}.subsampled_{sub_pol.target_depth:.0f}x.bam"

            cls.subsample_bam(
                input_bam=input_bam,
                output_bam=subsampled_bam,
                target_fraction=sampling_fraction,
                seed=sub_pol.seed,
                threads=threads,
            )

            post_depth = cls.compute_contig_mean_depth(
                bam_path=subsampled_bam,
                contig_name=reference_info.contig_name,
                contig_length=reference_info.length_bp,
            )
            logger.info(
                f"[{sample_id}] Subsampling applied: {orig_depth:.1f}x -> {post_depth:.1f}x "
                f"(fraction={sampling_fraction:.3f}, seed={sub_pol.seed})"
            )
            source_bam_for_calmd = subsampled_bam
            was_subsampled = True
        else:
            logger.info(
                f"[{sample_id}] Subsampling skipped (depth {orig_depth:.1f}x <= {sub_pol.min_depth_threshold:.1f}x or disabled)."
            )

        target_bam = work_dir / f"{sample_id}.schmutzi_ready.bam"
        cls.run_calmd_atomic(
            input_bam=source_bam_for_calmd,
            reference_fasta=reference_info.fasta_path,
            output_bam=target_bam,
            threads=threads,
        )

        target_bai = cls.ensure_bai_index(
            bam_path=target_bam,
            reference_contig=reference_info.contig_name,
            threads=threads,
        )

        try:
            with pysam.AlignmentFile(str(target_bam), "rb") as bam:
                mapped_count = bam.count(contig=reference_info.contig_name, read_callback="all")
                if mapped_count == 0:
                    raise SchmutziPreparationError(
                        f"Prepared BAM '{target_bam.name}' has 0 mapped reads on contig '{reference_info.contig_name}'."
                    )
        except Exception as exc:
            if isinstance(exc, SchmutziPreparationError):
                raise
            raise SchmutziPreparationError(f"Failed mapped read counting on {target_bam.name}: {exc}") from exc

        logger.info(
            f"[{sample_id}] schmutzi BAM preparation complete -> {target_bam.name} "
            f"({mapped_count:,} mapped reads verified, quickcheck PASS)."
        )
        return target_bam, target_bai, was_subsampled, orig_depth, post_depth