"""
Block 1: Reference Genome Validation, Indexing & Cryptographic Auditing.

Single-file step responsible for:
- Reference FASTA sequence inspection (contig ID, length, strict base checks)
- Deterministic cryptographic hashing (File SHA-256 and Sequence MD5)
- Automated index verification and generation (samtools faidx, bwa index, picard dict)
- Contig alias resolution (NC_012920.1 <-> chrM <-> MT)
"""

from dataclasses import dataclass, field
import hashlib
import logging
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from src.utils import calculate_sha256, EnvironmentValidator, PipelineConfig

logger = logging.getLogger("mtDNA_Pipeline.Reference")


class ReferenceManagerError(Exception):
    """Raised strictly on unrecoverable I/O or external binary execution failures."""
    pass


@dataclass
class ReferenceInfo:
    """Mitochondrial reference genome metadata dynamically resolved against input FASTA."""

    fasta_path: Path
    sequence_id: Optional[str] = None         # Header ID extracted from FASTA (e.g. NC_012920.1)
    contig_name: Optional[str] = None         # Contig name mapped in BAM (@SQ SN:)
    description: Optional[str] = None         # Full FASTA header comment
    length_bp: Optional[int] = None           # True base pair count parsed from FASTA
    is_circular: bool = True                  # Biological circularity flag
    file_sha256: Optional[str] = None         # SHA-256 of physical FASTA file on disk
    sequence_md5: Optional[str] = None        # MD5 computed strictly on uppercase biological sequence bases
    has_bwa_index: bool = False
    has_fai_index: bool = False

    def to_dict(self) -> Dict[str, Any]:
        """Serializes reference metadata."""
        return {
            "fasta_path": str(self.fasta_path),
            "sequence_id": self.sequence_id,
            "contig_name": self.contig_name,
            "description": self.description,
            "length_bp": self.length_bp,
            "is_circular": self.is_circular,
            "file_sha256": self.file_sha256,
            "sequence_md5": self.sequence_md5,
            "has_bwa_index": self.has_bwa_index,
            "has_fai_index": self.has_fai_index,
        }


class ContigAliasResolver:
    """Standardizes mitochondrial contig nomenclature across RefSeq, UCSC, and Ensembl conventions."""

    CANONICAL_RCRS = "NC_012920.1"

    KNOWN_MT_ALIASES: Dict[str, Set[str]] = {
        CANONICAL_RCRS: {
            "NC_012920.1",
            "NC_012920",
            "CHRM",
            "MT",
            "CHRMT",
            "RCRS",
            "GI|251831106|REF|NC_012920.1|",
        }
    }

    @classmethod
    def are_aliases(cls, contig_a: str, contig_b: str) -> bool:
        """Checks if two contig names represent equivalent naming conventions."""
        a_norm = contig_a.strip().upper()
        b_norm = contig_b.strip().upper()

        if a_norm == b_norm:
            return True

        for canonical, aliases in cls.KNOWN_MT_ALIASES.items():
            norm_aliases = {alias.upper() for alias in aliases}
            if a_norm in norm_aliases and b_norm in norm_aliases:
                return True

        return False

    @classmethod
    def resolve_canonical(cls, contig_name: str) -> str:
        """Resolves a contig name to its canonical RefSeq identifier if registered."""
        norm = contig_name.strip().upper()
        for canonical, aliases in cls.KNOWN_MT_ALIASES.items():
            if norm in {alias.upper() for alias in aliases}:
                return canonical
        return contig_name.strip()


@dataclass
class ReferenceValidationResult:
    """Summary of comprehensive integrity, sequence, and indexing checks for a reference genome."""

    is_valid: bool
    reference_info: ReferenceInfo
    checks_passed: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, object]:
        """Serializes validation summary."""
        return {
            "is_valid": self.is_valid,
            "reference": self.reference_info.to_dict(),
            "checks_passed": self.checks_passed,
            "warnings": self.warnings,
            "errors": self.errors,
        }


class ReferenceManager:
    """Manages validation, cryptographic auditing, cross-consistency, and indexing of the reference FASTA."""

    CHUNK_SIZE = 65536
    BWA_EXTENSIONS = (".bwt", ".pac", ".ann", ".amb", ".sa")
    STRICT_BASES = frozenset(b"ACGT")
    BASES_WITH_N = frozenset(b"ACGTN")

    @classmethod
    def compute_file_sha256(cls, file_path: Path) -> str:
        """Computes the SHA-256 cryptographic hash across the physical reference file on disk."""
        path = Path(file_path).resolve()
        if not path.is_file():
            raise ReferenceManagerError(f"Reference file does not exist: {path}")
        return calculate_sha256(path)

    @classmethod
    def parse_and_hash_sequence(
        cls,
        file_path: Path,
        allow_n: bool = False,
    ) -> Tuple[str, str, str, int, List[str], List[str]]:
        """Parses the FASTA stream, extracting sequence ID, description, uppercase bases, and MD5."""
        path = Path(file_path).resolve()
        if not path.is_file():
            raise ReferenceManagerError(f"Reference file not found: {path}")

        warnings: List[str] = []
        errors: List[str] = []

        seq_id = ""
        description = ""
        total_length = 0
        contig_count = 0
        md5_hasher = hashlib.md5()
        valid_charset = cls.BASES_WITH_N if allow_n else cls.STRICT_BASES

        with open(path, "rb") as f:
            for line in f:
                line_stripped = line.strip()
                if not line_stripped:
                    continue

                if line_stripped.startswith(b">"):
                    contig_count += 1
                    if contig_count == 1:
                        header_str = line_stripped[1:].decode("utf-8", errors="replace")
                        parts = header_str.split(None, 1)
                        seq_id = parts[0]
                        description = parts[1] if len(parts) > 1 else ""
                    else:
                        errors.append(
                            f"Multiple contigs detected in mitochondrial reference: {contig_count}. "
                            "Target mitochondrial reference must contain exactly one contig."
                        )
                else:
                    if contig_count == 0:
                        errors.append("Invalid FASTA format: Sequence content preceding contig header.")
                        break

                    clean_line = line_stripped.upper()
                    line_chars = frozenset(clean_line)
                    invalid_chars = line_chars - valid_charset

                    if invalid_chars:
                        invalid_str = "".join(chr(c) for c in sorted(invalid_chars))
                        errors.append(
                            f"Non-standard nucleotide bases found in sequence: '{invalid_str}'. "
                            f"(allow_n={allow_n}, allowed={sorted([chr(b) for b in valid_charset])})"
                        )

                    total_length += len(clean_line)
                    md5_hasher.update(clean_line)

        if contig_count == 0:
            errors.append(f"No FASTA contigs found in {path.name}.")

        seq_md5 = md5_hasher.hexdigest() if total_length > 0 else ""
        return seq_md5, seq_id, description, total_length, warnings, errors

    @classmethod
    def check_bwa_index(cls, fasta_path: Path) -> Tuple[bool, str]:
        """Checks if all 5 BWA-aln index files are present and strictly newer than the FASTA."""
        path = Path(fasta_path).resolve()
        if not path.is_file():
            return False, f"FASTA file does not exist: {path}"

        fasta_mtime = path.stat().st_mtime
        for ext in cls.BWA_EXTENSIONS:
            idx_file = path.parent / f"{path.name}{ext}"
            if not idx_file.is_file():
                return False, f"Missing BWA index file: {idx_file.name}"
            if idx_file.stat().st_size == 0:
                return False, f"Empty BWA index file: {idx_file.name}"
            if idx_file.stat().st_mtime < fasta_mtime:
                return False, f"BWA index file {idx_file.name} is older than reference FASTA."

        return True, "All BWA index files verified."

    @classmethod
    def parse_and_validate_fai(
        cls,
        fasta_path: Path,
        expected_contig: str,
        expected_length: int,
    ) -> Tuple[bool, str]:
        """Validates the samtools FAI index file for single-contig mitochondrial reference."""
        path = Path(fasta_path).resolve()
        fai_file = path.parent / f"{path.name}.fai"

        if not fai_file.is_file():
            return False, f"FAI index missing: {fai_file.name}"

        if fai_file.stat().st_size == 0:
            return False, f"FAI index is empty: {fai_file.name}"

        if fai_file.stat().st_mtime < path.stat().st_mtime:
            return False, f"FAI index {fai_file.name} is older than reference FASTA."

        with open(fai_file, "r", encoding="utf-8") as f:
            lines = [ln.strip() for ln in f if ln.strip()]

        if len(lines) != 1:
            return False, f"FAI index contains {len(lines)} contigs, expected exactly 1."

        columns = lines[0].split("\t")
        if len(columns) < 2:
            return False, "Invalid FAI format: fewer than 2 columns."

        fai_contig = columns[0]
        try:
            fai_len = int(columns[1])
        except ValueError:
            return False, f"Invalid length in FAI index: '{columns[1]}'."

        if fai_contig != expected_contig:
            if not ContigAliasResolver.are_aliases(fai_contig, expected_contig):
                return False, f"Contig in FAI '{fai_contig}' does not match expected '{expected_contig}'."

        if fai_len != expected_length:
            return False, f"Length in FAI ({fai_len}) does not match expected FASTA length ({expected_length})."

        return True, "FAI index cross-validated successfully."

    @classmethod
    def generate_bwa_index(cls, fasta_path: Path, log_dir: Optional[Path] = None) -> None:
        """Executes 'bwa index' to generate reference index files."""
        path = Path(fasta_path).resolve()
        cmd = ["bwa", "index", str(path)]

        log_file = None
        if log_dir:
            log_dir.mkdir(parents=True, exist_ok=True)
            log_file = log_dir / f"{path.name}.bwa_index.log"

        logger.info(f"Generating BWA index for {path.name}...")
        try:
            with open(log_file, "w", encoding="utf-8") if log_file else subprocess.DEVNULL as log_fp:
                subprocess.run(
                    cmd,
                    stdout=log_fp,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
        except subprocess.CalledProcessError as exc:
            raise ReferenceManagerError(f"Failed to generate BWA index: {exc}") from exc
        except FileNotFoundError as exc:
            raise ReferenceManagerError("bwa executable not found in PATH.") from exc

    @classmethod
    def generate_fai_index(cls, fasta_path: Path, log_dir: Optional[Path] = None) -> None:
        """Executes 'samtools faidx' to generate the .fai index."""
        path = Path(fasta_path).resolve()
        cmd = ["samtools", "faidx", str(path)]

        log_file = None
        if log_dir:
            log_dir.mkdir(parents=True, exist_ok=True)
            log_file = log_dir / f"{path.name}.faidx.log"

        logger.info(f"Generating samtools FAI index for {path.name}...")
        try:
            with open(log_file, "w", encoding="utf-8") if log_file else subprocess.DEVNULL as log_fp:
                subprocess.run(
                    cmd,
                    stdout=log_fp,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
        except subprocess.CalledProcessError as exc:
            raise ReferenceManagerError(f"Failed to generate samtools FAI index: {exc}") from exc
        except FileNotFoundError as exc:
            raise ReferenceManagerError("samtools executable not found in PATH.") from exc

    @classmethod
    def prepare_and_validate_reference(
        cls,
        fasta_path: Path,
        expected_contig: str = "NC_012920.1",
        strict_contig: bool = True,
        expected_length: Optional[int] = 16569,
        expected_sequence_md5: Optional[str] = None,
        is_circular: bool = True,
        allow_n: bool = False,
        auto_index: bool = True,
        log_dir: Optional[Path] = None,
    ) -> ReferenceValidationResult:
        """Orchestrates reference verification, cryptographic auditing, alias resolution, and indexing."""
        path = Path(fasta_path).resolve()
        passed_checks: List[str] = []
        warnings: List[str] = []
        errors: List[str] = []

        if not path.is_file():
            raise ReferenceManagerError(f"Reference FASTA does not exist: {path}")

        file_sha = cls.compute_file_sha256(path)
        passed_checks.append("FILE_SHA256_COMPUTED")

        seq_md5, seq_id, desc, seq_len, parse_warn, parse_err = cls.parse_and_hash_sequence(
            file_path=path,
            allow_n=allow_n,
        )
        warnings.extend(parse_warn)
        errors.extend(parse_err)

        if not parse_err:
            passed_checks.append("FASTA_STRUCTURE_AND_BASES_VALID")

        contig_name = seq_id
        if seq_id == expected_contig:
            passed_checks.append("CONTIG_EXACT_MATCH")
        elif not strict_contig and ContigAliasResolver.are_aliases(seq_id, expected_contig):
            canonical = ContigAliasResolver.resolve_canonical(seq_id)
            warnings.append(
                f"Contig '{seq_id}' in FASTA matches expected '{expected_contig}' via naming alias (Canonical: {canonical})."
            )
            passed_checks.append("CONTIG_ALIAS_MATCH")
        else:
            errors.append(
                f"Contig identifier mismatch in {path.name}: found '{seq_id}', "
                f"expected '{expected_contig}' (strict_contig={strict_contig})."
            )

        if expected_length is not None and seq_len > 0:
            if seq_len == expected_length:
                passed_checks.append("LENGTH_EXACT_MATCH")
            else:
                errors.append(
                    f"Sequence length mismatch: expected {expected_length:,} bp, but parsed {seq_len:,} bp."
                )

        if expected_sequence_md5 is not None and seq_md5:
            if seq_md5.lower() == expected_sequence_md5.lower():
                passed_checks.append("SEQUENCE_MD5_IDENTITY_VERIFIED")
            else:
                errors.append(
                    f"Biological sequence MD5 mismatch: expected '{expected_sequence_md5}', got '{seq_md5}'."
                )

        has_bwa, bwa_reason = cls.check_bwa_index(path)
        if not has_bwa:
            if auto_index and not errors:
                cls.generate_bwa_index(path, log_dir=log_dir)
                has_bwa, bwa_reason = cls.check_bwa_index(path)
                if not has_bwa:
                    errors.append(f"BWA auto-indexing failed for {path.name}: {bwa_reason}")
                else:
                    passed_checks.append("BWA_INDEX_AUTO_GENERATED")
            else:
                errors.append(f"BWA index invalid or outdated: {bwa_reason}")
        else:
            passed_checks.append("BWA_INDEX_VERIFIED")

        has_fai, fai_reason = cls.parse_and_validate_fai(
            fasta_path=path,
            expected_contig=seq_id,
            expected_length=seq_len,
        )
        if not has_fai:
            if auto_index and not errors:
                cls.generate_fai_index(path, log_dir=log_dir)
                has_fai, fai_reason = cls.parse_and_validate_fai(
                    fasta_path=path,
                    expected_contig=seq_id,
                    expected_length=seq_len,
                )
                if not has_fai:
                    errors.append(f"FAI auto-indexing / cross-validation failed: {fai_reason}")
                else:
                    passed_checks.append("FAI_INDEX_AUTO_GENERATED")
            else:
                errors.append(f"FAI index invalid: {fai_reason}")
        else:
            passed_checks.append("FAI_INDEX_CROSS_VALIDATED")

        ref_info = ReferenceInfo(
            fasta_path=path,
            sequence_id=seq_id,
            contig_name=contig_name,
            description=desc,
            length_bp=seq_len,
            is_circular=is_circular,
            file_sha256=file_sha,
            sequence_md5=seq_md5,
            has_bwa_index=has_bwa,
            has_fai_index=has_fai,
        )

        is_valid = len(errors) == 0 and has_bwa and has_fai

        logger.info(
            f"Reference Validation {'PASSED' if is_valid else 'FAILED'} for {path.name} | "
            f"Contig: {contig_name} | Length: {seq_len:,} bp | "
            f"Seq-MD5: {seq_md5} | File-SHA256: {file_sha[:12]}..."
        )

        return ReferenceValidationResult(
            is_valid=is_valid,
            reference_info=ref_info,
            checks_passed=passed_checks,
            warnings=warnings,
            errors=errors,
        )


class ReferenceStep:
    """High-level runner for Pipeline Step 1 (Reference Genome Validation)."""

    @classmethod
    def run(
        cls,
        fasta_path: Union[Path, str, Any],
        expected_contig: str = "NC_012920.1",
        expected_length: int = 16569,
        expected_sequence_md5: Optional[str] = None,
        auto_index: bool = True,
        strict_contig: bool = True,
        allow_n: bool = True,
        log_dir: Optional[Path] = None,
    ) -> ReferenceInfo:
        """Executes reference validation and index verification, returning ReferenceInfo."""
        # Handle PipelineConfig instance if passed directly
        if hasattr(fasta_path, "reference_fasta") or hasattr(fasta_path, "reference"):
            if hasattr(fasta_path, "reference_fasta"):
                raw_fasta = getattr(fasta_path, "reference_fasta")
            else:
                ref_sec = getattr(fasta_path, "reference", {})
                if isinstance(ref_sec, dict):
                    expected_contig = ref_sec.get("expected_contig", expected_contig)
                    expected_length = ref_sec.get("expected_length", expected_length)
                    auto_index = ref_sec.get("auto_index", auto_index)
                    raw_fasta = ref_sec.get("fasta_path")
                else:
                    raw_fasta = getattr(ref_sec, "fasta_path", None)
            from src.utils import resolve_path
            fasta_path = resolve_path(str(raw_fasta)) if raw_fasta else Path("data/reference/rCRS.fasta")

        logger.info("=" * 80)
        logger.info("▶ BLOCK 1: REFERENCE GENOME VALIDATION (ReferenceStep)")
        logger.info("=" * 80)

        res = ReferenceManager.prepare_and_validate_reference(
            fasta_path=fasta_path,
            expected_contig=expected_contig,
            strict_contig=strict_contig,
            expected_length=expected_length,
            expected_sequence_md5=expected_sequence_md5,
            allow_n=allow_n,
            auto_index=auto_index,
            log_dir=log_dir,
        )

        if not res.is_valid:
            raise RuntimeError(f"Reference validation failed: {res.errors}")

        logger.info(
            f"✅ Reference verified: {res.reference_info.fasta_path.name} "
            f"({res.reference_info.contig_name}, {res.reference_info.length_bp} bp) "
            f"[SHA256: {res.reference_info.file_sha256[:12]}...]"
        )
        return res.reference_info
