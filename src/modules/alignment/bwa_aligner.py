"""BWA-aln multi-branch alignment engine for ancient DNA: parallel SA indexing, library-safe @RG injection, and streamed BAM compilation."""

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
import logging
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any, Dict, List, Optional, Tuple

import pysam

from src.modules.alignment.models import (
    AlignerType,
    AlignmentToolMetadata,
    ReadGroupInfo,
    ReadGroupSource,
    ReferenceInfo,
)
from src.modules.preprocessing.models import ProcessedReads, SampleRecord, ValidationStatus
from src.modules.preprocessing.integrity import IntegrityEngine

logger = logging.getLogger("mtDNA_Pipeline.BWAAligner")


class BWAAlignerError(Exception):
    """Raised when BWA alignment, streaming execution, or BAM validation fails."""
    pass


@dataclass
class AlignmentBranchResult:
    """Outcome and traceability artifacts for an individual alignment branch."""

    source: ReadGroupSource
    read_group: ReadGroupInfo
    output_bam: Path
    sai_files: List[Path] = field(default_factory=list)
    temporary_sai_deleted: bool = False
    input_units_count: int = 0
    unit_type: str = "reads"  # "fragments", "reads", or "pairs"
    command_lines: List[str] = field(default_factory=list)
    is_successful: bool = False
    error_message: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Serializes branch alignment outcome."""
        return {
            "source": self.source.value,
            "read_group": self.read_group.to_dict(),
            "output_bam": str(self.output_bam),
            "sai_files": [str(p) for p in self.sai_files],
            "temporary_sai_deleted": self.temporary_sai_deleted,
            "input_units_count": self.input_units_count,
            "unit_type": self.unit_type,
            "command_lines": self.command_lines,
            "is_successful": self.is_successful,
            "error_message": self.error_message,
        }


@dataclass
class AlignmentBranchResults:
    """Container aggregating all executed alignment branches for a sample."""

    sample_id: str
    branches: Dict[ReadGroupSource, AlignmentBranchResult] = field(default_factory=dict)
    tool_meta: AlignmentToolMetadata = field(default_factory=AlignmentToolMetadata)

    @property
    def all_bams(self) -> List[Path]:
        """Returns physical paths to all successfully generated branch BAM files."""
        return [b.output_bam for b in self.branches.values() if b.is_successful and b.output_bam.is_file()]

    @property
    def all_read_groups(self) -> List[ReadGroupInfo]:
        """Returns all injected ReadGroupInfo objects."""
        return [b.read_group for b in self.branches.values()]

    def to_dict(self) -> Dict[str, Any]:
        """Serializes aggregated branch results."""
        return {
            "sample_id": self.sample_id,
            "branches": {k.value: v.to_dict() for k, v in self.branches.items()},
            "tool_metadata": self.tool_meta.to_dict(),
        }


class BWAAligner:
    """Executes conservative BWA-aln passes and streams SAM output directly into compressed BAM files."""

    @classmethod
    def get_bwa_version(cls) -> Tuple[str, str]:
        """Inspects and returns the installed BWA version and raw banner."""
        if not shutil.which("bwa"):
            return "not_found", "BWA executable not found in PATH"

        try:
            res = subprocess.run(["bwa"], capture_output=True, text=True, check=False)
            raw = f"{res.stdout}\n{res.stderr}".strip()
            match = re.search(r"Version:\s*([\d\.\-a-zA-Z]+)", raw, re.IGNORECASE)
            version_str = match.group(1) if match else "unknown"
            return version_str, raw
        except Exception as exc:
            logger.warning(f"Could not determine BWA version: {exc}")
            return "unknown", str(exc)

    @classmethod
    def get_samtools_version(cls) -> Tuple[str, str]:
        """Inspects and returns the installed samtools version and raw banner."""
        if not shutil.which("samtools"):
            return "not_found", "samtools executable not found in PATH"

        try:
            res = subprocess.run(["samtools", "--version"], capture_output=True, text=True, check=False)
            raw = f"{res.stdout}\n{res.stderr}".strip()
            first_line = raw.split("\n")[0]
            match = re.search(r"samtools\s+([\d\.]+)", first_line, re.IGNORECASE)
            version_str = match.group(1) if match else "unknown"
            return version_str, raw
        except Exception as exc:
            logger.warning(f"Could not determine samtools version: {exc}")
            return "unknown", str(exc)

    @classmethod
    def _resolve_read_count(cls, fastq_path: Path, cached_val: Optional[int]) -> int:
        """Returns cached metric if valid and positive; otherwise falls back to fast physical FASTQ counting."""
        if cached_val is not None and cached_val > 0:
            return cached_val

        logger.debug(f"Cached read count not available or zero for {fastq_path.name}. Executing fast physical count.")
        return IntegrityEngine.fast_count_fastq_reads(fastq_path)

    @classmethod
    def _build_read_group(
        cls,
        sample_id: str,
        library_id: str,
        source: ReadGroupSource,
        platform: str = "ILLUMINA",
        platform_unit: Optional[str] = None,
        sequencing_center: Optional[str] = None,
    ) -> ReadGroupInfo:
        """Constructs a collision-resistant, library-safe ReadGroupInfo."""
        source_suffix = source.value.upper()
        # Costruzione ID library-safe (es. SAMPLE1_LIB1_MERGED)
        if library_id and library_id != sample_id:
            rg_id = f"{sample_id}_{library_id}_{source_suffix}"
            pu_val = platform_unit or f"{sample_id}_{library_id}_{source_suffix}_RUN"
        else:
            rg_id = f"{sample_id}_{source_suffix}"
            pu_val = platform_unit or f"{sample_id}_{source_suffix}_RUN"

        return ReadGroupInfo(
            read_group_id=rg_id,
            sample_id=sample_id,
            library_id=library_id,
            platform=platform,
            platform_unit=pu_val,
            sequencing_center=sequencing_center,
            description=f"Branch: {source.value}",
            read_group_source=source,
        )

    @classmethod
    def _run_aln_pass(
        cls,
        fastq_path: Path,
        reference_path: Path,
        output_sai: Path,
        threads: int,
        seed_length: int,
        max_diff: float,
        max_gap_opens: int,
        max_gap_extensions: int,
        trim_quality_threshold: int,
        log_dir: Path,
    ) -> str:
        """Executes 'bwa aln' to generate suffix coordinates (.sai)."""
        output_sai.parent.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)

        log_err = log_dir / f"{output_sai.stem}.bwa_aln.err"

        cmd = [
            "bwa", "aln",
            "-t", str(max(1, threads)),
            "-l", str(seed_length),
            "-n", str(max_diff),
            "-o", str(max_gap_opens),
        ]
        if max_gap_extensions >= 0:
            cmd.extend(["-e", str(max_gap_extensions)])
        if trim_quality_threshold > 0:
            cmd.extend(["-q", str(trim_quality_threshold)])

        cmd.extend([
            str(reference_path),
            str(fastq_path),
        ])

        cmd_str = " ".join(cmd)
        logger.debug(f"Executing BWA aln: {cmd_str}")

        try:
            with open(output_sai, "wb") as out_f, open(log_err, "w", encoding="utf-8") as err_f:
                res = subprocess.run(cmd, stdout=out_f, stderr=err_f, check=False)
        except OSError as exc:
            raise BWAAlignerError(f"Failed to launch 'bwa aln' for {fastq_path.name}: {exc}") from exc

        if res.returncode != 0 or not output_sai.is_file() or output_sai.stat().st_size == 0:
            err_msg = f"'bwa aln' failed (exit code {res.returncode}) for {fastq_path.name}. Check log: {log_err}"
            raise BWAAlignerError(err_msg)

        return cmd_str

    @classmethod
    def _run_samse_pipe(
        cls,
        fastq_path: Path,
        sai_path: Path,
        reference_path: Path,
        rg_info: ReadGroupInfo,
        output_bam: Path,
        threads: int,
        log_dir: Path,
    ) -> Tuple[str, str]:
        """Streams 'bwa samse' directly into 'samtools view -b' creating a branch BAM."""
        output_bam.parent.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)

        samse_err = log_dir / f"{output_bam.stem}.bwa_samse.err"
        samtools_err = log_dir / f"{output_bam.stem}.samtools_view.err"

        cmd_samse = [
            "bwa", "samse",
            "-r", rg_info.to_bwa_rg_string(),
            str(reference_path),
            str(sai_path),
            str(fastq_path),
        ]

        cmd_samtools = [
            "samtools", "view",
            "-@", str(max(1, threads - 1)),
            "-b",
            "-o", str(output_bam),
            "-",
        ]

        cmd_samse_str = " ".join(cmd_samse)
        cmd_samtools_str = " ".join(cmd_samtools)
        logger.debug(f"Streaming: {cmd_samse_str} | {cmd_samtools_str}")

        with open(samse_err, "w", encoding="utf-8") as err_samse_f, open(samtools_err, "w", encoding="utf-8") as err_samtools_f:
            p_samse = subprocess.Popen(cmd_samse, stdout=subprocess.PIPE, stderr=err_samse_f)
            p_samtools = subprocess.Popen(
                cmd_samtools,
                stdin=p_samse.stdout,
                stdout=subprocess.DEVNULL,
                stderr=err_samtools_f,
            )

            if p_samse.stdout:
                p_samse.stdout.close()

            p_samtools.communicate()
            ret_samse = p_samse.wait()
            ret_samtools = p_samtools.returncode

        if ret_samse != 0:
            raise BWAAlignerError(f"'bwa samse' failed with code {ret_samse}. Log: {samse_err}")
        if ret_samtools != 0:
            raise BWAAlignerError(f"'samtools view' failed with code {ret_samtools}. Log: {samtools_err}")

        return cmd_samse_str, cmd_samtools_str

    @classmethod
    def _run_sampe_pipe(
        cls,
        r1_fastq: Path,
        r2_fastq: Path,
        r1_sai: Path,
        r2_sai: Path,
        reference_path: Path,
        rg_info: ReadGroupInfo,
        output_bam: Path,
        threads: int,
        log_dir: Path,
    ) -> Tuple[str, str]:
        """Streams 'bwa sampe' directly into 'samtools view -b' creating a paired branch BAM."""
        output_bam.parent.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)

        sampe_err = log_dir / f"{output_bam.stem}.bwa_sampe.err"
        samtools_err = log_dir / f"{output_bam.stem}.samtools_view.err"

        cmd_sampe = [
            "bwa", "sampe",
            "-r", rg_info.to_bwa_rg_string(),
            str(reference_path),
            str(r1_sai),
            str(r2_sai),
            str(r1_fastq),
            str(r2_fastq),
        ]

        cmd_samtools = [
            "samtools", "view",
            "-@", str(max(1, threads - 1)),
            "-b",
            "-o", str(output_bam),
            "-",
        ]

        cmd_sampe_str = " ".join(cmd_sampe)
        cmd_samtools_str = " ".join(cmd_samtools)
        logger.debug(f"Streaming: {cmd_sampe_str} | {cmd_samtools_str}")

        with open(sampe_err, "w", encoding="utf-8") as err_sampe_f, open(samtools_err, "w", encoding="utf-8") as err_samtools_f:
            p_sampe = subprocess.Popen(cmd_sampe, stdout=subprocess.PIPE, stderr=err_sampe_f)
            p_samtools = subprocess.Popen(
                cmd_samtools,
                stdin=p_sampe.stdout,
                stdout=subprocess.DEVNULL,
                stderr=err_samtools_f,
            )

            if p_sampe.stdout:
                p_sampe.stdout.close()

            p_samtools.communicate()
            ret_sampe = p_sampe.wait()
            ret_samtools = p_samtools.returncode

        if ret_sampe != 0:
            raise BWAAlignerError(f"'bwa sampe' failed with code {ret_sampe}. Log: {sampe_err}")
        if ret_samtools != 0:
            raise BWAAlignerError(f"'samtools view' failed with code {ret_samtools}. Log: {samtools_err}")

        return cmd_sampe_str, cmd_samtools_str

    @classmethod
    def _verify_branch_bam(
        cls,
        bam_path: Path,
        rg_info: ReadGroupInfo,
        expected_contig: Optional[str] = None,
    ) -> None:
        """Validates BAM physical integrity via samtools quickcheck and rigorously inspects ID/SM/LB/PL and contigs."""
        if not bam_path.is_file() or bam_path.stat().st_size == 0:
            raise BWAAlignerError(f"Branch BAM output is missing or empty: {bam_path.name}")

        # 1. Controllo profondo di integrità BGZF con samtools quickcheck
        if shutil.which("samtools"):
            res_qc = subprocess.run(["samtools", "quickcheck", "-v", str(bam_path)], capture_output=True, text=True, check=False)
            if res_qc.returncode != 0:
                raise BWAAlignerError(
                    f"BAM integrity verification failed ('samtools quickcheck' exit code {res_qc.returncode}) for {bam_path.name}: {res_qc.stderr.strip()}"
                )

        # 2. Ispezione logica e semantica degli header (@RG ed @SQ)
        try:
            with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam_f:
                header = bam_f.header.to_dict()

                rgs = header.get("RG", [])
                target_rg = None
                for rg in rgs:
                    if rg.get("ID") == rg_info.read_group_id:
                        target_rg = rg
                        break

                if target_rg is None:
                    found_ids = [rg.get("ID") for rg in rgs if "ID" in rg]
                    raise BWAAlignerError(
                        f"Expected Read Group ID '{rg_info.read_group_id}' not found in BAM header of {bam_path.name}. Found: {found_ids}"
                    )

                # Verifica semantica approfondita dei campi obbligatori
                if target_rg.get("SM") != rg_info.sample_id:
                    raise BWAAlignerError(
                        f"Read Group SM mismatch in {bam_path.name}: expected '{rg_info.sample_id}', found '{target_rg.get('SM')}'"
                    )
                if target_rg.get("LB") != rg_info.library_id:
                    raise BWAAlignerError(
                        f"Read Group LB mismatch in {bam_path.name}: expected '{rg_info.library_id}', found '{target_rg.get('LB')}'"
                    )
                if target_rg.get("PL") != rg_info.platform:
                    raise BWAAlignerError(
                        f"Read Group PL mismatch in {bam_path.name}: expected '{rg_info.platform}', found '{target_rg.get('PL')}'"
                    )

                if expected_contig:
                    sqs = header.get("SQ", [])
                    sq_names = [sq.get("SN") for sq in sqs if "SN" in sq]
                    if expected_contig not in sq_names:
                        raise BWAAlignerError(
                            f"Expected contig '{expected_contig}' missing from @SQ headers in {bam_path.name}."
                        )

        except BWAAlignerError:
            raise
        except Exception as exc:
            raise BWAAlignerError(f"Corrupted branch BAM header structure in {bam_path.name}: {exc}") from exc

    @classmethod
    def align_single_end_branch(
        cls,
        fastq_path: Path,
        source: ReadGroupSource,
        sample_id: str,
        library_id: str,
        reference_info: ReferenceInfo,
        output_dir: Path,
        log_dir: Path,
        cached_count: Optional[int] = None,
        threads: int = 4,
        seed_length: int = 1000,
        max_diff: float = 0.01,
        max_gap_opens: int = 2,
        max_gap_extensions: int = -1,
        trim_quality_threshold: int = 0,
        platform: str = "ILLUMINA",
        platform_unit: Optional[str] = None,
        sequencing_center: Optional[str] = None,
        keep_temporary_sai: bool = False,
    ) -> AlignmentBranchResult:
        """Executes full alignment lifecycle for a Single-End branch (Merged, Merged-Truncated, or Singleton)."""
        input_count = cls._resolve_read_count(fastq_path, cached_count)
        unit_type = "fragments" if source in (ReadGroupSource.MERGED, ReadGroupSource.MERGED_TRUNCATED) else "reads"

        rg_info = cls._build_read_group(
            sample_id=sample_id,
            library_id=library_id,
            source=source,
            platform=platform,
            platform_unit=platform_unit,
            sequencing_center=sequencing_center,
        )

        branch_tag = source.value
        sai_path = output_dir / "tmp" / f"{sample_id}_{library_id}_{branch_tag}.sai"
        bam_path = output_dir / f"{sample_id}_{library_id}_{branch_tag}.raw.bam"
        executed_cmds: List[str] = []

        logger.info(f"[{sample_id}] Aligning {branch_tag.upper()} branch ({input_count:,} {unit_type}) -> {bam_path.name}")

        try:
            # 1. BWA aln pass
            cmd_aln = cls._run_aln_pass(
                fastq_path=fastq_path,
                reference_path=reference_info.fasta_path,
                output_sai=sai_path,
                threads=threads,
                seed_length=seed_length,
                max_diff=max_diff,
                max_gap_opens=max_gap_opens,
                max_gap_extensions=max_gap_extensions,
                trim_quality_threshold=trim_quality_threshold,
                log_dir=log_dir,
            )
            executed_cmds.append(cmd_aln)

            # 2. BWA samse streaming pass
            cmd_samse, cmd_view = cls._run_samse_pipe(
                fastq_path=fastq_path,
                sai_path=sai_path,
                reference_path=reference_info.fasta_path,
                rg_info=rg_info,
                output_bam=bam_path,
                threads=threads,
                log_dir=log_dir,
            )
            executed_cmds.extend([cmd_samse, cmd_view])

            # 3. Validazione BAM approfondita
            cls._verify_branch_bam(
                bam_path=bam_path,
                rg_info=rg_info,
                expected_contig=reference_info.contig_name,
            )

            # 4. Gestione file temporaneo .sai
            sai_deleted = False
            if not keep_temporary_sai and sai_path.is_file():
                sai_path.unlink()
                sai_deleted = True

            return AlignmentBranchResult(
                source=source,
                read_group=rg_info,
                output_bam=bam_path,
                sai_files=[sai_path],
                temporary_sai_deleted=sai_deleted,
                input_units_count=input_count,
                unit_type=unit_type,
                command_lines=executed_cmds,
                is_successful=True,
            )

        except Exception as exc:
            logger.error(f"[{sample_id}] Branch {branch_tag} failed: {exc}")
            return AlignmentBranchResult(
                source=source,
                read_group=rg_info,
                output_bam=bam_path,
                sai_files=[sai_path],
                temporary_sai_deleted=False,
                input_units_count=input_count,
                unit_type=unit_type,
                command_lines=executed_cmds,
                is_successful=False,
                error_message=str(exc),
            )

    @classmethod
    def align_paired_end_branch(
        cls,
        r1_fastq: Path,
        r2_fastq: Path,
        sample_id: str,
        library_id: str,
        reference_info: ReferenceInfo,
        output_dir: Path,
        log_dir: Path,
        cached_pairs_count: Optional[int] = None,
        threads: int = 4,
        seed_length: int = 1000,
        max_diff: float = 0.01,
        max_gap_opens: int = 2,
        max_gap_extensions: int = -1,
        trim_quality_threshold: int = 0,
        platform: str = "ILLUMINA",
        platform_unit: Optional[str] = None,
        sequencing_center: Optional[str] = None,
        keep_temporary_sai: bool = False,
        parallel_aln: bool = True,
    ) -> AlignmentBranchResult:
        """Executes full alignment lifecycle for intact unmerged Paired-End mates using parallel SA indexing."""
        if cached_pairs_count is not None and cached_pairs_count > 0:
            pairs_count = cached_pairs_count
        else:
            r1_count = cls._resolve_read_count(r1_fastq, None)
            r2_count = cls._resolve_read_count(r2_fastq, None)
            if r1_count != r2_count:
                raise BWAAlignerError(
                    f"Paired-end input count mismatch for {sample_id}: R1 has {r1_count:,} reads, R2 has {r2_count:,} reads."
                )
            pairs_count = r1_count

        rg_info = cls._build_read_group(
            sample_id=sample_id,
            library_id=library_id,
            source=ReadGroupSource.PAIRED,
            platform=platform,
            platform_unit=platform_unit,
            sequencing_center=sequencing_center,
        )

        r1_sai = output_dir / "tmp" / f"{sample_id}_{library_id}_paired_R1.sai"
        r2_sai = output_dir / "tmp" / f"{sample_id}_{library_id}_paired_R2.sai"
        bam_path = output_dir / f"{sample_id}_{library_id}_paired.raw.bam"
        executed_cmds: List[str] = []

        logger.info(f"[{sample_id}] Aligning PAIRED branch ({pairs_count:,} pairs) -> {bam_path.name}")

        try:
            # 1. BWA aln su R1 ed R2 in parallelo
            if parallel_aln and threads >= 2:
                threads_per_mate = max(1, threads // 2)
                with ThreadPoolExecutor(max_workers=2) as executor:
                    future_r1 = executor.submit(
                        cls._run_aln_pass,
                        fastq_path=r1_fastq,
                        reference_path=reference_info.fasta_path,
                        output_sai=r1_sai,
                        threads=threads_per_mate,
                        seed_length=seed_length,
                        max_diff=max_diff,
                        max_gap_opens=max_gap_opens,
                        max_gap_extensions=max_gap_extensions,
                        trim_quality_threshold=trim_quality_threshold,
                        log_dir=log_dir,
                    )
                    future_r2 = executor.submit(
                        cls._run_aln_pass,
                        fastq_path=r2_fastq,
                        reference_path=reference_info.fasta_path,
                        output_sai=r2_sai,
                        threads=threads_per_mate,
                        seed_length=seed_length,
                        max_diff=max_diff,
                        max_gap_opens=max_gap_opens,
                        max_gap_extensions=max_gap_extensions,
                        trim_quality_threshold=trim_quality_threshold,
                        log_dir=log_dir,
                    )
                    cmd_aln1 = future_r1.result()
                    cmd_aln2 = future_r2.result()
                    executed_cmds.extend([cmd_aln1, cmd_aln2])
            else:
                cmd_aln1 = cls._run_aln_pass(
                    fastq_path=r1_fastq,
                    reference_path=reference_info.fasta_path,
                    output_sai=r1_sai,
                    threads=threads,
                    seed_length=seed_length,
                    max_diff=max_diff,
                    max_gap_opens=max_gap_opens,
                    max_gap_extensions=max_gap_extensions,
                    trim_quality_threshold=trim_quality_threshold,
                    log_dir=log_dir,
                )
                cmd_aln2 = cls._run_aln_pass(
                    fastq_path=r2_fastq,
                    reference_path=reference_info.fasta_path,
                    output_sai=r2_sai,
                    threads=threads,
                    seed_length=seed_length,
                    max_diff=max_diff,
                    max_gap_opens=max_gap_opens,
                    max_gap_extensions=max_gap_extensions,
                    trim_quality_threshold=trim_quality_threshold,
                    log_dir=log_dir,
                )
                executed_cmds.extend([cmd_aln1, cmd_aln2])

            # 2. BWA sampe streaming pass
            cmd_sampe, cmd_view = cls._run_sampe_pipe(
                r1_fastq=r1_fastq,
                r2_fastq=r2_fastq,
                r1_sai=r1_sai,
                r2_sai=r2_sai,
                reference_path=reference_info.fasta_path,
                rg_info=rg_info,
                output_bam=bam_path,
                threads=threads,
                log_dir=log_dir,
            )
            executed_cmds.extend([cmd_sampe, cmd_view])

            # 3. Validazione BAM approfondita
            cls._verify_branch_bam(
                bam_path=bam_path,
                rg_info=rg_info,
                expected_contig=reference_info.contig_name,
            )

            # 4. Gestione file temporanei .sai
            sai_deleted = False
            if not keep_temporary_sai:
                if r1_sai.is_file():
                    r1_sai.unlink()
                if r2_sai.is_file():
                    r2_sai.unlink()
                sai_deleted = True

            return AlignmentBranchResult(
                source=ReadGroupSource.PAIRED,
                read_group=rg_info,
                output_bam=bam_path,
                sai_files=[r1_sai, r2_sai],
                temporary_sai_deleted=sai_deleted,
                input_units_count=pairs_count,
                unit_type="pairs",
                command_lines=executed_cmds,
                is_successful=True,
            )

        except Exception as exc:
            logger.error(f"[{sample_id}] Paired branch alignment failed: {exc}")
            return AlignmentBranchResult(
                source=ReadGroupSource.PAIRED,
                read_group=rg_info,
                output_bam=bam_path,
                sai_files=[r1_sai, r2_sai],
                temporary_sai_deleted=False,
                input_units_count=pairs_count,
                unit_type="pairs",
                command_lines=executed_cmds,
                is_successful=False,
                error_message=str(exc),
            )

    @classmethod
    def align_sample_branches(
        cls,
        sample: SampleRecord,
        reference_info: ReferenceInfo,
        output_dir: Path,
        log_dir: Path,
        threads: int = 4,
        seed_length: int = 1000,
        max_diff: float = 0.01,
        max_gap_opens: int = 2,
        max_gap_extensions: int = -1,
        trim_quality_threshold: int = 0,
        platform: str = "ILLUMINA",
        platform_unit: Optional[str] = None,
        sequencing_center: Optional[str] = None,
        keep_temporary_sai: bool = False,
        parallel_aln: bool = True,
    ) -> AlignmentBranchResults:
        """Dispatches all active read branches (merged, merged_truncated, singleton, paired) for a SampleRecord.

        Args:
            sample: Target SampleRecord containing populated ProcessedReads and PreprocessingMetrics.
            reference_info: Validated ReferenceInfo instance.
            output_dir: Output directory for raw branch BAMs.
            log_dir: Output directory for stdout/stderr logs.
            threads: CPU threads allocated for alignment.
            seed_length: Seed length (-l parameter).
            max_diff: Maximum edit distance / difference rate (-n parameter).
            max_gap_opens: Maximum open gaps (-o parameter).
            max_gap_extensions: Maximum gap extensions (-e parameter).
            trim_quality_threshold: Quality trimming threshold (-q parameter).
            platform: Platform tag for @RG PL.
            platform_unit: Platform unit tag for @RG PU.
            sequencing_center: Sequencing center tag for @RG CN.
            keep_temporary_sai: If True, retains intermediate .sai index files on disk.
            parallel_aln: If True, aligns PE mates R1/R2 concurrently.

        Returns:
            AlignmentBranchResults aggregating all aligned branch BAMs and audit metadata.
        """
        sample_id = sample.sample_id
        library_id = sample.identity.library_id or sample_id
        processed: ProcessedReads = sample.processed_reads

        bwa_ver, _ = cls.get_bwa_version()
        sam_ver, _ = cls.get_samtools_version()

        tool_meta = AlignmentToolMetadata(
            aligner_name=AlignerType.BWA_ALN,
            bwa_version=bwa_ver,
            samtools_version=sam_ver,
            reference_file_sha256=reference_info.file_sha256,
            working_directory=str(output_dir.resolve()),
            parameters={
                "seed_length": seed_length,
                "max_diff": max_diff,
                "max_gap_opens": max_gap_opens,
                "max_gap_extensions": max_gap_extensions,
                "trim_quality_threshold": trim_quality_threshold,
                "threads": threads,
            },
        )
        tool_meta.mark_started()

        branch_results: Dict[ReadGroupSource, AlignmentBranchResult] = {}

        # 1. Ramo Frammenti Merged (Full-length collapsed)
        if processed.merged_fastq and processed.merged_fastq.is_file() and processed.merged_fastq.stat().st_size > 0:
            cached_m = sample.metrics.merged_fragments if sample.metrics else None
            res_merged = cls.align_single_end_branch(
                fastq_path=processed.merged_fastq,
                source=ReadGroupSource.MERGED,
                sample_id=sample_id,
                library_id=library_id,
                reference_info=reference_info,
                output_dir=output_dir,
                log_dir=log_dir,
                cached_count=cached_m,
                threads=threads,
                seed_length=seed_length,
                max_diff=max_diff,
                max_gap_opens=max_gap_opens,
                max_gap_extensions=max_gap_extensions,
                trim_quality_threshold=trim_quality_threshold,
                platform=platform,
                platform_unit=platform_unit,
                sequencing_center=sequencing_center,
                keep_temporary_sai=keep_temporary_sai,
            )
            branch_results[ReadGroupSource.MERGED] = res_merged
            for cmd in res_merged.command_lines:
                tool_meta.add_command(cmd)

        # 2. Ramo Frammenti Merged Truncated (Collapsed con trimming terminale)
        merged_trunc_path = getattr(processed, "merged_truncated_fastq", None)
        if merged_trunc_path and merged_trunc_path.is_file() and merged_trunc_path.stat().st_size > 0:
            res_trunc = cls.align_single_end_branch(
                fastq_path=merged_trunc_path,
                source=ReadGroupSource.MERGED_TRUNCATED,
                sample_id=sample_id,
                library_id=library_id,
                reference_info=reference_info,
                output_dir=output_dir,
                log_dir=log_dir,
                threads=threads,
                seed_length=seed_length,
                max_diff=max_diff,
                max_gap_opens=max_gap_opens,
                max_gap_extensions=max_gap_extensions,
                trim_quality_threshold=trim_quality_threshold,
                platform=platform,
                platform_unit=platform_unit,
                sequencing_center=sequencing_center,
                keep_temporary_sai=keep_temporary_sai,
            )
            branch_results[ReadGroupSource.MERGED_TRUNCATED] = res_trunc
            for cmd in res_trunc.command_lines:
                tool_meta.add_command(cmd)

        # 3. Ramo Singletons (Letture orfane)
        if processed.singleton_fastq and processed.singleton_fastq.is_file() and processed.singleton_fastq.stat().st_size > 0:
            cached_s = sample.metrics.singleton_reads if sample.metrics else None
            res_singleton = cls.align_single_end_branch(
                fastq_path=processed.singleton_fastq,
                source=ReadGroupSource.SINGLETON,
                sample_id=sample_id,
                library_id=library_id,
                reference_info=reference_info,
                output_dir=output_dir,
                log_dir=log_dir,
                cached_count=cached_s,
                threads=threads,
                seed_length=seed_length,
                max_diff=max_diff,
                max_gap_opens=max_gap_opens,
                max_gap_extensions=max_gap_extensions,
                trim_quality_threshold=trim_quality_threshold,
                platform=platform,
                platform_unit=platform_unit,
                sequencing_center=sequencing_center,
                keep_temporary_sai=keep_temporary_sai,
            )
            branch_results[ReadGroupSource.SINGLETON] = res_singleton
            for cmd in res_singleton.command_lines:
                tool_meta.add_command(cmd)

        # 4. Ramo Coppie Paired-End non collassate
        has_r1 = processed.paired_r1_fastq and processed.paired_r1_fastq.is_file() and processed.paired_r1_fastq.stat().st_size > 0
        has_r2 = processed.paired_r2_fastq and processed.paired_r2_fastq.is_file() and processed.paired_r2_fastq.stat().st_size > 0

        if has_r1 and has_r2:
            cached_p = sample.metrics.surviving_read_pairs if sample.metrics else None
            res_paired = cls.align_paired_end_branch(
                r1_fastq=processed.paired_r1_fastq,  # type: ignore[arg-type]
                r2_fastq=processed.paired_r2_fastq,  # type: ignore[arg-type]
                sample_id=sample_id,
                library_id=library_id,
                reference_info=reference_info,
                output_dir=output_dir,
                log_dir=log_dir,
                cached_pairs_count=cached_p,
                threads=threads,
                seed_length=seed_length,
                max_diff=max_diff,
                max_gap_opens=max_gap_opens,
                max_gap_extensions=max_gap_extensions,
                trim_quality_threshold=trim_quality_threshold,
                platform=platform,
                platform_unit=platform_unit,
                sequencing_center=sequencing_center,
                keep_temporary_sai=keep_temporary_sai,
                parallel_aln=parallel_aln,
            )
            branch_results[ReadGroupSource.PAIRED] = res_paired
            for cmd in res_paired.command_lines:
                tool_meta.add_command(cmd)

        tool_meta.mark_finished()

        successful_branches = [b for b in branch_results.values() if b.is_successful]
        if not successful_branches:
            failed_msgs = [f"[{b.source.value}]: {b.error_message}" for b in branch_results.values() if not b.is_successful]
            raise BWAAlignerError(
                f"No alignment branches succeeded for {sample_id}. Failures: {'; '.join(failed_msgs)}"
            )

        return AlignmentBranchResults(
            sample_id=sample_id,
            branches=branch_results,
            tool_meta=tool_meta,
        )
