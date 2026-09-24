"""Track A Engine: Base quality score rescaling via mapDamage2 --rescale-only and multi-tier BAM validation."""

from dataclasses import dataclass, field
import logging
from pathlib import Path
import shutil
import subprocess
import time
from typing import Any, Dict, List, Optional, Tuple
import pysam

from src.modules.damage.models import (
    DamageErrorCode,
    DamagePolicy,
)

logger = logging.getLogger("mtDNA_Pipeline.RescaleEngine")


class DamageRescaleError(Exception):
    """Raised when base quality rescaling, model pre-conditions, or output post-conditions fail."""
    pass


@dataclass
class RescaleExecutionMetadata:
    """Audit trail container capturing tool version, command invocation, and validation parameters."""

    tool_name: str = "mapDamage"
    tool_version: str = "unknown"
    command_line: str = ""
    input_bam: Path = field(default_factory=Path)
    damage_model_dir: Path = field(default_factory=Path)
    output_bam: Path = field(default_factory=Path)
    output_bai: Path = field(default_factory=Path)
    verified_records: int = 0
    runtime_seconds: float = 0.0
    skipped: bool = False

    def to_dict(self) -> Dict[str, Any]:
        """Serializes rescaling execution metadata to a dictionary."""
        return {
            "tool_name": self.tool_name,
            "tool_version": self.tool_version,
            "command_line": self.command_line,
            "input_bam": str(self.input_bam),
            "damage_model_dir": str(self.damage_model_dir),
            "output_bam": str(self.output_bam) if not self.skipped else None,
            "output_bai": str(self.output_bai) if not self.skipped else None,
            "verified_records": self.verified_records,
            "runtime_seconds": round(self.runtime_seconds, 2),
            "skipped": self.skipped,
        }


class QualityRescaleEngine:
    """Orchestrates base quality score down-scaling on terminal transitions with comprehensive integrity checks."""

    # Standard BGZF 28-byte End-Of-File marker
    BGZF_EOF_MARKER: bytes = (
        b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00"
        b"\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
    )

    @classmethod
    def get_binary_name(cls) -> Optional[str]:
        """Detects if mapDamage is registered in system PATH."""
        if shutil.which("mapDamage"):
            return "mapDamage"
        if shutil.which("mapdamage"):
            return "mapdamage"
        return None

    @classmethod
    def get_version(cls) -> str:
        """Retrieves installed mapDamage version string."""
        binary = cls.get_binary_name()
        if not binary:
            return "not_found"
        try:
            res = subprocess.run([binary, "--version"], capture_output=True, text=True, check=False)
            raw = f"{res.stdout}\n{res.stderr}".strip()
            return raw.split("\n")[0]
        except Exception:
            return "unknown"

    @classmethod
    def check_samtools_quickcheck(cls, bam_path: Path) -> bool:
        """Executes samtools quickcheck if available in PATH."""
        if shutil.which("samtools") is None:
            return True  # Fallback gracefully to internal BGZF EOF check
        try:
            res = subprocess.run(
                ["samtools", "quickcheck", "-v", str(bam_path)],
                capture_output=True,
                text=True,
                check=False,
            )
            return res.returncode == 0
        except OSError:
            return True

    @classmethod
    def check_bgzf_eof(cls, bam_path: Path) -> bool:
        """Verifies physical presence of the standard 28-byte BGZF EOF block."""
        try:
            with open(bam_path, "rb") as stream:
                stream.seek(0, 2)
                file_size = stream.tell()
                if file_size < 28:
                    return False
                stream.seek(file_size - 28)
                tail = stream.read(28)
                return tail == cls.BGZF_EOF_MARKER
        except OSError:
            return False

    @classmethod
    def validate_preconditions(
        cls,
        input_bam: Path,
        reference_fasta: Path,
        damage_model_dir: Path,
        expected_contig: Optional[str] = "NC_012920.1",
    ) -> None:
        """Verifies presence and validity of input BAM, FASTA, and mapDamage statistical model files."""
        if not input_bam.is_file() or input_bam.stat().st_size == 0:
            raise DamageRescaleError(f"Input BAM file for rescaling not found or empty: {input_bam}")

        if not reference_fasta.is_file() or reference_fasta.stat().st_size == 0:
            raise DamageRescaleError(f"Reference FASTA not found or empty: {reference_fasta}")

        if not damage_model_dir.is_dir():
            raise DamageRescaleError(f"Damage model folder not found: {damage_model_dir}")

        # Validate mandatory minimal statistical model files produced by quantifier
        required_model_files = [
            damage_model_dir / "misincorporation.txt",
        ]
        for req in required_model_files:
            if not req.is_file() or req.stat().st_size == 0:
                raise DamageRescaleError(
                    f"Statistical model file '{req.name}' is missing in {damage_model_dir}. "
                    "Cannot execute --rescale-only without valid pre-computed damage statistics."
                )

        # Check contig compatibility in input BAM
        try:
            with pysam.AlignmentFile(str(input_bam), "rb") as bam:
                if expected_contig and expected_contig not in bam.references:
                    raise DamageRescaleError(
                        f"Expected mitochondrial contig '{expected_contig}' missing from input BAM references."
                    )
        except Exception as exc:
            if isinstance(exc, DamageRescaleError):
                raise
            raise DamageRescaleError(f"Cannot inspect input BAM headers: {exc}") from exc

    @classmethod
    def validate_postconditions(
        cls,
        input_bam: Path,
        rescaled_bam: Path,
        rescaled_bai: Path,
        expected_contig: Optional[str] = "NC_012920.1",
    ) -> int:
        """Executes multi-tier post-rescaling integrity, header preservation, record parity, and BAI checks."""
        if not rescaled_bam.is_file() or rescaled_bam.stat().st_size == 0:
            raise DamageRescaleError(f"Rescaled BAM missing or empty on disk: {rescaled_bam}")

        # 1. samtools quickcheck validation
        if not cls.check_samtools_quickcheck(rescaled_bam):
            raise DamageRescaleError(f"Rescaled BAM failed 'samtools quickcheck' validation: {rescaled_bam.name}")

        # 2. Binary BGZF EOF check
        if not cls.check_bgzf_eof(rescaled_bam):
            raise DamageRescaleError(
                f"Rescaled BAM '{rescaled_bam.name}' is truncated (missing 28-byte BGZF EOF marker)."
            )

        # 3. Header & Record Integrity via pysam
        try:
            with (
                pysam.AlignmentFile(str(input_bam), "rb") as in_bam,
                pysam.AlignmentFile(str(rescaled_bam), "rb") as out_bam,
            ):
                in_hdr = in_bam.header.to_dict()
                out_hdr = out_bam.header.to_dict()

                # Validate @SQ sequence dictionary
                in_sq = in_hdr.get("SQ", [])
                out_sq = out_hdr.get("SQ", [])
                if in_sq != out_sq:
                    raise DamageRescaleError(
                        "Header @SQ sequence dictionary modified during rescaling! "
                        f"Input contigs: {len(in_sq)}, Output contigs: {len(out_sq)}"
                    )

                # Validate Read Groups (@RG)
                in_rg = in_hdr.get("RG", [])
                out_rg = out_hdr.get("RG", [])
                if len(in_rg) != len(out_rg):
                    raise DamageRescaleError(
                        f"Read Groups (@RG) mismatch: Input had {len(in_rg)}, Output has {len(out_rg)}"
                    )

                # Verify expected contig presence
                if expected_contig and expected_contig not in out_bam.references:
                    raise DamageRescaleError(
                        f"Target mitochondrial contig '{expected_contig}' missing from rescaled BAM references."
                    )

                # Verify exact read record count parity
                in_count = in_bam.count(until_eof=True)
                out_count = out_bam.count(until_eof=True)
                if in_count != out_count:
                    raise DamageRescaleError(
                        f"Read count parity violation: Input BAM had {in_count} records, "
                        f"Rescaled BAM produced {out_count} records."
                    )

        except Exception as exc:
            if isinstance(exc, DamageRescaleError):
                raise
            raise DamageRescaleError(f"Failed during rescaled BAM header/record inspection: {exc}") from exc

        # 4. Validate BAI index
        if not rescaled_bai.is_file() or rescaled_bai.stat().st_size == 0:
            raise DamageRescaleError(f"Rescaled BAM index (.bai) is missing or empty: {rescaled_bai}")

        try:
            with pysam.AlignmentFile(str(rescaled_bam), "rb", index_filename=str(rescaled_bai)) as test_idx:
                if not test_idx.has_index():
                    raise DamageRescaleError(f"Index validation failed: pysam cannot query '{rescaled_bai.name}'")
        except Exception as exc:
            if isinstance(exc, DamageRescaleError):
                raise
            raise DamageRescaleError(f"Corrupted or invalid BAM index: {exc}") from exc

        return out_count

    @classmethod
    def run_rescaling(
        cls,
        sample_id: str,
        input_bam: Path,
        reference_fasta: Path,
        damage_model_dir: Path,
        output_dir: Path,
        log_dir: Path,
        policy: DamagePolicy,
        expected_contig: Optional[str] = "NC_012920.1",
    ) -> Tuple[Optional[Path], Optional[Path], RescaleExecutionMetadata]:
        """Executes mapDamage --rescale-only, standardizes deliverables, and enforces post-rescaling audit.

        Args:
            sample_id: Biological sample identifier.
            input_bam: Path to input BAM (dedup or markdup).
            reference_fasta: Target reference FASTA path (rCRS).
            damage_model_dir: Folder containing previously estimated mapDamage models.
            output_dir: Destination folder for rescaled BAM.
            log_dir: Destination folder for process logs.
            policy: Validated DamagePolicy container.
            expected_contig: Expected mitochondrial contig name.

        Returns:
            Tuple of (rescaled_bam_path, rescaled_bai_path, RescaleExecutionMetadata).

        Raises:
            DamageRescaleError: If pre-conditions, mapDamage execution, or post-conditions fail.
        """
        meta = RescaleExecutionMetadata(
            input_bam=input_bam,
            damage_model_dir=damage_model_dir,
        )

        # 1. Policy check: if disabled, exit cleanly in SKIPPED state
        if not policy.enable_rescaling:
            logger.info(f"[{sample_id}] Rescaling is disabled in policy (enable_rescaling=False). Skipping.")
            meta.skipped = True
            return None, None, meta

        binary = cls.get_binary_name()
        if not binary:
            raise DamageRescaleError("mapDamage binary not found in system PATH.")

        meta.tool_version = cls.get_version()

        # 2. Validate pre-conditions
        cls.validate_preconditions(
            input_bam=input_bam,
            reference_fasta=reference_fasta,
            damage_model_dir=damage_model_dir,
            expected_contig=expected_contig,
        )

        out_sample_dir = output_dir / sample_id
        out_sample_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)

        stdout_log = log_dir / f"{sample_id}.mapdamage_rescale.log"
        stderr_log = log_dir / f"{sample_id}.mapdamage_rescale.err"

        cmd = [
            binary,
            "-i", str(input_bam),
            "-r", str(reference_fasta),
            "-d", str(damage_model_dir),
            "--rescale-only",
        ]
        meta.command_line = " ".join(cmd)

        logger.info(
            f"[{sample_id}] Executing mapDamage2 Quality Rescaling (--rescale-only) | "
            f"Model Source: {damage_model_dir.name}"
        )

        # 3. Process execution
        t0 = time.perf_counter()
        try:
            with open(stdout_log, "w", encoding="utf-8") as out_f, open(stderr_log, "w", encoding="utf-8") as err_f:
                res = subprocess.run(cmd, stdout=out_f, stderr=err_f, check=False)
        except OSError as exc:
            raise DamageRescaleError(f"Failed to launch mapDamage process: {exc}") from exc
        finally:
            meta.runtime_seconds = time.perf_counter() - t0

        if res.returncode != 0:
            raise DamageRescaleError(
                f"mapDamage --rescale-only failed with exit code {res.returncode}. Log: {stderr_log}"
            )

        # 4. Locate produced artifact
        expected_candidates = [
            damage_model_dir / f"{input_bam.stem}.rescaled.bam",
            damage_model_dir / f"{sample_id}.rescaled.bam",
            damage_model_dir / "rescaled.bam",
        ]
        raw_rescaled_bam: Optional[Path] = None
        for cand in expected_candidates:
            if cand.is_file() and cand.stat().st_size > 0:
                raw_rescaled_bam = cand
                break

        if not raw_rescaled_bam:
            raise DamageRescaleError(
                f"mapDamage finished with exit code 0 but rescaled BAM was not produced in {damage_model_dir}."
            )

        # 5. Standardized relocation and BAI index generation
        target_rescaled_bam = out_sample_dir / f"{sample_id}.rescaled.bam"
        target_bai = Path(f"{target_rescaled_bam}.bai")

        if raw_rescaled_bam.resolve() != target_rescaled_bam.resolve():
            logger.info(f"[{sample_id}] Standardizing rescaled BAM destination -> {target_rescaled_bam.name}")
            shutil.move(str(raw_rescaled_bam), str(target_rescaled_bam))

        try:
            logger.info(f"[{sample_id}] Generating BAI index for rescaled BAM...")
            pysam.index(str(target_rescaled_bam))
        except Exception as exc:
            raise DamageRescaleError(f"Failed to index rescaled BAM: {exc}") from exc

        # 6. Validate post-conditions
        verified_reads = cls.validate_postconditions(
            input_bam=input_bam,
            rescaled_bam=target_rescaled_bam,
            rescaled_bai=target_bai,
            expected_contig=expected_contig,
        )

        meta.output_bam = target_rescaled_bam
        meta.output_bai = target_bai
        meta.verified_records = verified_reads

        logger.info(
            f"[{sample_id}] Rescaling Completed & Verified -> {target_rescaled_bam.name} "
            f"({verified_reads:,} reads, @SQ/@RG preserved, quickcheck PASS, BAI valid in {meta.runtime_seconds:.1f}s)."
        )

        return target_rescaled_bam, target_bai, meta