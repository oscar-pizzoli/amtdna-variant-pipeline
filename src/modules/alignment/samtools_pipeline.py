"""Structural BAM preparation engine: branch merging, name-collation, mate synchronization via fixmate -m, and coordinate sorting."""

from dataclasses import dataclass, field
import logging
from pathlib import Path
import shutil
import subprocess
from typing import Any, Dict, List, Optional, Tuple

import pysam

from src.modules.alignment.models import ReadGroupInfo
from src.modules.alignment.bwa_aligner import AlignmentBranchResults

logger = logging.getLogger("mtDNA_Pipeline.SamtoolsPipeline")


class SamtoolsPipelineError(Exception):
    """Raised when BAM merging, streaming execution, sorting, indexing, or verification fails."""
    pass


@dataclass
class PipelinePreparationResult:
    """Output artifacts and audit records resulting from structural BAM preparation."""

    sample_id: str
    merged_raw_bam: Optional[Path] = None
    sorted_bam: Optional[Path] = None           # sample.prepared.sorted.bam (fixmate + coordinate sort)
    sorted_bai: Optional[Path] = None
    baseline_stats_txt: Optional[Path] = None   # samtools stats pre-filtro
    input_branch_records: int = 0               # Somma record dichiarati dai rami di input
    prepared_records_count: int = 0             # Record fisici nel BAM ordinato
    command_lines: List[str] = field(default_factory=list)
    is_successful: bool = False
    error_message: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Serializes preparation outcome."""
        return {
            "sample_id": self.sample_id,
            "merged_raw_bam": str(self.merged_raw_bam) if self.merged_raw_bam else None,
            "sorted_bam": str(self.sorted_bam) if self.sorted_bam else None,
            "sorted_bai": str(self.sorted_bai) if self.sorted_bai else None,
            "baseline_stats_txt": str(self.baseline_stats_txt) if self.baseline_stats_txt else None,
            "input_branch_records": self.input_branch_records,
            "prepared_records_count": self.prepared_records_count,
            "command_lines": self.command_lines,
            "is_successful": self.is_successful,
            "error_message": self.error_message,
        }


class SamtoolsPipeline:
    """Orchestrates BAM branch merging, UNIX pipe collate -> fixmate -m -> sort, indexing, and baseline stats."""

    @classmethod
    def index_bam(cls, bam_path: Path, threads: int = 2) -> Path:
        """Generates coordinate index (.bai) for a target BAM file."""
        if not shutil.which("samtools"):
            raise SamtoolsPipelineError("Cannot index BAM: 'samtools' executable not found in system PATH.")

        bam = Path(bam_path).resolve()
        if not bam.is_file() or bam.stat().st_size == 0:
            raise SamtoolsPipelineError(f"Cannot index missing or empty BAM: {bam.name}")

        bai_path = bam.with_name(bam.name + ".bai")
        cmd = ["samtools", "index", "-@", str(max(1, threads)), str(bam)]

        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if res.returncode != 0 or not bai_path.is_file():
            raise SamtoolsPipelineError(
                f"Failed to generate index for {bam.name} (exit code {res.returncode}): {res.stderr.strip()}"
            )

        return bai_path

    @classmethod
    def generate_samtools_stats(cls, bam_path: Path, output_stats_txt: Path, threads: int = 2) -> str:
        """Runs 'samtools stats' on the prepared BAM to generate comprehensive baseline metrics."""
        output_stats_txt.parent.mkdir(parents=True, exist_ok=True)
        cmd = ["samtools", "stats", "-@", str(max(1, threads)), str(bam_path)]
        cmd_str = " ".join(cmd)

        try:
            with open(output_stats_txt, "w", encoding="utf-8") as out_f:
                res = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True, check=False)

            if res.returncode != 0:
                raise SamtoolsPipelineError(
                    f"'samtools stats' failed on {bam_path.name} (exit code {res.returncode}): {res.stderr.strip()}"
                )
            return cmd_str
        except OSError as exc:
            raise SamtoolsPipelineError(f"Failed to write samtools stats to {output_stats_txt.name}: {exc}") from exc

    @classmethod
    def verify_prepared_bam(
        cls,
        bam_path: Path,
        expected_read_groups: List[ReadGroupInfo],
        expected_contig: Optional[str] = None,
    ) -> int:
        """Performs fast integrity/truncation check, validates header consistency (SO:coordinate, @RG, @SQ), and returns total records.

        Returns:
            Total alignment records count in the BAM file.
        """
        if not bam_path.is_file() or bam_path.stat().st_size == 0:
            raise SamtoolsPipelineError(f"Prepared BAM is missing or empty: {bam_path.name}")

        # 1. Controllo rapido di integrità e assenza di troncamento BGZF
        if shutil.which("samtools"):
            res_qc = subprocess.run(
                ["samtools", "quickcheck", "-v", str(bam_path)],
                capture_output=True,
                text=True,
                check=False,
            )
            if res_qc.returncode != 0:
                raise SamtoolsPipelineError(
                    f"samtools quickcheck failed (exit code {res_qc.returncode}) on {bam_path.name}: {res_qc.stderr.strip()}"
                )

        # 2. Ispezione degli header e conteggio record con pysam
        total_records = 0
        try:
            with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam_f:
                # Sanity check: dichiarazione dell'ordinamento nell'header
                so_flag = bam_f.header.get("HD", {}).get("SO")
                if so_flag != "coordinate":
                    logger.warning(
                        f"BAM header SO flag is '{so_flag}' instead of 'coordinate' in {bam_path.name}"
                    )

                header = bam_f.header.to_dict()
                existing_rgs = {rg["ID"]: rg for rg in header.get("RG", []) if "ID" in rg}

                for expected_rg in expected_read_groups:
                    rg_id = expected_rg.read_group_id
                    if rg_id not in existing_rgs:
                        raise SamtoolsPipelineError(
                            f"Missing Read Group '{rg_id}' in prepared BAM header {bam_path.name}. "
                            f"Present RGs: {list(existing_rgs.keys())}"
                        )
                    actual_rg = existing_rgs[rg_id]
                    if actual_rg.get("SM") != expected_rg.sample_id:
                        raise SamtoolsPipelineError(
                            f"Corrupted SM in Read Group '{rg_id}': expected '{expected_rg.sample_id}', got '{actual_rg.get('SM')}'"
                        )
                    if actual_rg.get("LB") != expected_rg.library_id:
                        raise SamtoolsPipelineError(
                            f"Corrupted LB in Read Group '{rg_id}': expected '{expected_rg.library_id}', got '{actual_rg.get('LB')}'"
                        )
                    if actual_rg.get("PL") != expected_rg.platform:
                        raise SamtoolsPipelineError(
                            f"Corrupted PL in Read Group '{rg_id}': expected '{expected_rg.platform}', got '{actual_rg.get('PL')}'"
                        )
                    if expected_rg.platform_unit and actual_rg.get("PU") != expected_rg.platform_unit:
                        raise SamtoolsPipelineError(
                            f"Corrupted PU in Read Group '{rg_id}': expected '{expected_rg.platform_unit}', got '{actual_rg.get('PU')}'"
                        )

                if expected_contig:
                    sq_names = [sq["SN"] for sq in header.get("SQ", []) if "SN" in sq]
                    if expected_contig not in sq_names:
                        raise SamtoolsPipelineError(
                            f"Expected contig '{expected_contig}' missing from @SQ headers in {bam_path.name}."
                        )

                total_records = bam_f.count(until_eof=True)

        except SamtoolsPipelineError:
            raise
        except Exception as exc:
            raise SamtoolsPipelineError(f"Corrupted BAM structure in {bam_path.name}: {exc}") from exc

        return total_records

    @classmethod
    def merge_branch_bams(
        cls,
        branch_bams: List[Path],
        output_merged_bam: Path,
        threads: int,
        log_dir: Path,
    ) -> str:
        """Merges multiple raw branch BAM files into a single unified BAM preserving all distinct Read Groups."""
        output_merged_bam.parent.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)

        log_err = log_dir / f"{output_merged_bam.stem}.samtools_merge.err"

        if len(branch_bams) == 1:
            # Singolo ramo attivo: copia atomica
            shutil.copyfile(branch_bams[0], output_merged_bam)
            return f"cp {branch_bams[0]} {output_merged_bam}"

        cmd = [
            "samtools", "merge",
            "-@", str(max(1, threads)),
            "-f",
            str(output_merged_bam),
        ] + [str(b) for b in branch_bams]

        cmd_str = " ".join(cmd)
        logger.debug(f"Merging {len(branch_bams)} branch BAMs: {cmd_str}")

        with open(log_err, "w", encoding="utf-8") as err_f:
            res = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=err_f, check=False)

        if res.returncode != 0 or not output_merged_bam.is_file() or output_merged_bam.stat().st_size == 0:
            raise SamtoolsPipelineError(
                f"'samtools merge' failed (exit code {res.returncode}) for {output_merged_bam.name}. Check log: {log_err}"
            )

        return cmd_str

    @classmethod
    def run_collate_fixmate_sort_stream(
        cls,
        input_bam: Path,
        output_sorted_bam: Path,
        threads: int,
        log_dir: Path,
    ) -> List[str]:
        """Executes streamed UNIX pipe: 'collate -u -O' -> 'fixmate -m -u' -> 'sort -@ sort_threads' without intermediate disk files."""
        output_sorted_bam.parent.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)

        sample_stem = output_sorted_bam.stem
        collate_err = log_dir / f"{sample_stem}.samtools_collate.err"
        fixmate_err = log_dir / f"{sample_stem}.samtools_fixmate.err"
        sort_err = log_dir / f"{sample_stem}.samtools_sort.err"

        # Ripartizione esplicita dei thread
        sort_threads = max(1, threads)
        collate_threads = max(1, min(4, threads // 2))

        cmd_collate = [
            "samtools", "collate",
            "-@", str(collate_threads),
            "-u",
            "-O",
            str(input_bam),
            "-",
        ]

        cmd_fixmate = [
            "samtools", "fixmate",
            "-m",
            "-u",
            "-",
            "-",
        ]

        cmd_sort = [
            "samtools", "sort",
            "-@", str(sort_threads),
            "-o", str(output_sorted_bam),
            "-",
        ]

        cmds_str = [" ".join(cmd_collate), " ".join(cmd_fixmate), " ".join(cmd_sort)]
        logger.debug(f"Streaming pipe: {' | '.join(cmds_str)}")

        with open(collate_err, "w", encoding="utf-8") as f_collate, \
             open(fixmate_err, "w", encoding="utf-8") as f_fixmate, \
             open(sort_err, "w", encoding="utf-8") as f_sort:

            p_collate = subprocess.Popen(cmd_collate, stdout=subprocess.PIPE, stderr=f_collate)
            p_fixmate = subprocess.Popen(cmd_fixmate, stdin=p_collate.stdout, stdout=subprocess.PIPE, stderr=f_fixmate)
            p_sort = subprocess.Popen(cmd_sort, stdin=p_fixmate.stdout, stdout=subprocess.DEVNULL, stderr=f_sort)

            # Inoltro corretto dei segnali SIGPIPE
            if p_collate.stdout:
                p_collate.stdout.close()
            if p_fixmate.stdout:
                p_fixmate.stdout.close()

            p_sort.communicate()

            ret_collate = p_collate.wait()
            ret_fixmate = p_fixmate.wait()
            ret_sort = p_sort.returncode

        if ret_collate != 0:
            raise SamtoolsPipelineError(f"'samtools collate' failed (code {ret_collate}). Log: {collate_err}")
        if ret_fixmate != 0:
            raise SamtoolsPipelineError(f"'samtools fixmate' failed (code {ret_fixmate}). Log: {fixmate_err}")
        if ret_sort != 0 or not output_sorted_bam.is_file() or output_sorted_bam.stat().st_size == 0:
            raise SamtoolsPipelineError(f"'samtools sort' failed (code {ret_sort}). Log: {sort_err}")

        return cmds_str

    @classmethod
    def process_sample_branches(
        cls,
        branch_results: AlignmentBranchResults,
        output_dir: Path,
        log_dir: Path,
        threads: int = 4,
        keep_raw_merged_bam: bool = True,
        cleanup_on_failure: bool = True,
        expected_contig: Optional[str] = "NC_012920.1",
    ) -> PipelinePreparationResult:
        """Dispatches full BAM merge, streaming fixmate sync, coordinate sort, indexing, stats, and validation.

        Args:
            branch_results: Aggregated output container from BWAAligner.
            output_dir: Output folder for prepared BAMs.
            log_dir: Folder for stdout/stderr logs.
            threads: Number of CPU threads allocated.
            keep_raw_merged_bam: If True, preserves sample.merged.raw.bam on disk for auditing.
            cleanup_on_failure: If True, automatically unlinks partial output files if the pipeline fails.
            expected_contig: Contig identifier expected in headers.

        Returns:
            PipelinePreparationResult container with verified paths and audit logs.
        """
        sample_id = branch_results.sample_id
        valid_bams = branch_results.all_bams
        read_groups = branch_results.all_read_groups

        if not valid_bams:
            raise SamtoolsPipelineError(f"[{sample_id}] No valid branch BAMs available to process.")

        output_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)

        merged_raw_bam = output_dir / f"{sample_id}.merged.raw.bam"
        sorted_bam = output_dir / f"{sample_id}.prepared.sorted.bam"
        sorted_bai: Optional[Path] = None
        stats_txt = output_dir / f"{sample_id}.prepared.stats.txt"
        executed_cmds: List[str] = []

        # Calcolo totale record attesi dai rami BWA
        expected_records = sum(
            (b.input_units_count * 2 if b.unit_type == "pairs" else b.input_units_count)
            for b in branch_results.branches.values()
            if b.is_successful
        )

        logger.info(f"[{sample_id}] Starting structural BAM preparation from {len(valid_bams)} branch(es)...")

        try:
            # 1. Fusione fisica dei rami grezzi
            cmd_merge = cls.merge_branch_bams(
                branch_bams=valid_bams,
                output_merged_bam=merged_raw_bam,
                threads=threads,
                log_dir=log_dir,
            )
            executed_cmds.append(cmd_merge)

            # 2. Pipeline in streaming UNIX pipe: collate -> fixmate -m -> coordinate sort
            cmds_stream = cls.run_collate_fixmate_sort_stream(
                input_bam=merged_raw_bam,
                output_sorted_bam=sorted_bam,
                threads=threads,
                log_dir=log_dir,
            )
            executed_cmds.extend(cmds_stream)

            # 3. Indicizzazione del BAM ordinato per coordinate
            cmd_idx = f"samtools index -@ {threads} {sorted_bam}"
            sorted_bai = cls.index_bam(sorted_bam, threads=threads)
            executed_cmds.append(cmd_idx)

            # 4. Validazione e conteggio record
            actual_records = cls.verify_prepared_bam(
                bam_path=sorted_bam,
                expected_read_groups=read_groups,
                expected_contig=expected_contig,
            )

            # 5. Generazione statistiche baseline (pre-filtro)
            cmd_stats = cls.generate_samtools_stats(
                bam_path=sorted_bam,
                output_stats_txt=stats_txt,
                threads=max(1, threads // 2),
            )
            executed_cmds.append(cmd_stats)

            # 6. Pulizia del BAM grezzo se richiesto dalla policy (production mode)
            if not keep_raw_merged_bam and merged_raw_bam.is_file():
                merged_raw_bam.unlink()
                merged_raw_bam = None

            logger.info(
                f"[{sample_id}] Structural preparation completed -> {sorted_bam.name} "
                f"({actual_records:,} records verified, Stats generated)"
            )

            return PipelinePreparationResult(
                sample_id=sample_id,
                merged_raw_bam=merged_raw_bam,
                sorted_bam=sorted_bam,
                sorted_bai=sorted_bai,
                baseline_stats_txt=stats_txt,
                input_branch_records=expected_records,
                prepared_records_count=actual_records,
                command_lines=executed_cmds,
                is_successful=True,
            )

        except Exception as exc:
            logger.error(f"[{sample_id}] Structural BAM preparation failed: {exc}")

            # Pulizia difensiva dei file parziali su errore
            if cleanup_on_failure:
                for partial_file in (sorted_bam, sorted_bai, stats_txt):
                    if partial_file and partial_file.is_file():
                        try:
                            partial_file.unlink()
                        except OSError:
                            pass
                if not keep_raw_merged_bam and merged_raw_bam and merged_raw_bam.is_file():
                    try:
                        merged_raw_bam.unlink()
                    except OSError:
                        pass

            return PipelinePreparationResult(
                sample_id=sample_id,
                merged_raw_bam=merged_raw_bam if (merged_raw_bam and merged_raw_bam.is_file()) else None,
                sorted_bam=None,
                sorted_bai=None,
                baseline_stats_txt=None,
                input_branch_records=expected_records,
                prepared_records_count=0,
                command_lines=executed_cmds,
                is_successful=False,
                error_message=str(exc),
            )
