"""Mitochondrial reference genome management, sequence validation, cryptographic auditing, and indexing."""

from dataclasses import dataclass, field
import hashlib
import logging
from pathlib import Path
import re
import shutil
import subprocess
from typing import Dict, List, Optional, Set, Tuple

from src.modules.alignment.models import ReferenceInfo

logger = logging.getLogger("mtDNA_Pipeline.ReferenceManager")


class ReferenceManagerError(Exception):
    """Raised strictly on unrecoverable I/O or external binary execution failures."""
    pass


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

        sha256 = hashlib.sha256()
        try:
            with open(path, "rb") as stream:
                for block in iter(lambda: stream.read(cls.CHUNK_SIZE), b""):
                    sha256.update(block)
            return sha256.hexdigest()
        except OSError as exc:
            raise ReferenceManagerError(f"Failed to read reference file for hashing {path.name}: {exc}") from exc

    @classmethod
    def parse_and_hash_sequence(
        cls,
        file_path: Path,
        allow_n: bool = False,
    ) -> Tuple[str, str, Optional[str], int, List[str], List[str]]:
        """Parses FASTA reference, validates single-contig structure, and computes biological sequence MD5.

        Note:
            Sequence MD5 is computed on concatenated uppercase nucleotide bases with all whitespace and
            newlines stripped, strictly adhering to the standard SAM/BAM @SQ M5 specification.

        Returns:
            Tuple of (sequence_md5, sequence_id, description, sequence_length_bp, validation_warnings, validation_errors).
        """
        path = Path(file_path).resolve()
        if not path.is_file():
            raise ReferenceManagerError(f"FASTA reference file not found: {path}")

        valid_bases = cls.BASES_WITH_N if allow_n else cls.STRICT_BASES
        sequence_id: Optional[str] = None
        description: Optional[str] = None
        seq_md5 = hashlib.md5()
        total_len = 0
        contig_count = 0
        warnings: List[str] = []
        errors: List[str] = []

        try:
            with open(path, "rb") as stream:
                for line_idx, line in enumerate(stream, start=1):
                    # Rimozione di tutti gli spazi e ritorni a capo
                    line_clean = re.sub(rb"\s+", b"", line)
                    if not line_clean:
                        continue

                    if line_clean.startswith(b">"):
                        contig_count += 1
                        if contig_count > 1:
                            errors.append(
                                f"Multiple contigs detected in {path.name} at line {line_idx}. "
                                "Strict mitochondrial mode requires an isolated single-contig reference."
                            )
                            break

                        header_line = line.strip()[1:].decode("utf-8", errors="replace").strip()
                        parts = header_line.split(maxsplit=1)
                        sequence_id = parts[0]
                        description = parts[1] if len(parts) > 1 else None
                    else:
                        upper_line = line_clean.upper()
                        # Controllo caratteri nucleotidici non ammessi
                        invalid_chars = set(upper_line) - valid_bases
                        if invalid_chars:
                            invalid_repr = "".join(chr(c) for c in sorted(invalid_chars))
                            errors.append(
                                f"Invalid nucleotide characters '{invalid_repr}' found at line {line_idx}. "
                                f"Allowed bases: {''.join(chr(b) for b in sorted(valid_bases))}"
                            )
                            break

                        seq_md5.update(upper_line)
                        total_len += len(upper_line)

        except OSError as exc:
            raise ReferenceManagerError(f"I/O error reading FASTA sequence from {path.name}: {exc}") from exc

        if contig_count == 0 or sequence_id is None:
            errors.append(f"FASTA file {path.name} contains no sequence headers (missing '>').")

        if total_len == 0 and not errors:
            errors.append(f"FASTA file {path.name} contains 0 base pairs.")

        return seq_md5.hexdigest(), sequence_id or "UNKNOWN", description, total_len, warnings, errors

    @classmethod
    def check_bwa_index(cls, fasta_path: Path) -> Tuple[bool, Optional[str]]:
        """Checks if all required BWA index files exist, are non-empty, and are not older than the FASTA."""
        fasta = Path(fasta_path).resolve()
        if not fasta.is_file():
            return False, "Reference FASTA not found."

        fasta_mtime = fasta.stat().st_mtime

        for ext in cls.BWA_EXTENSIONS:
            index_file = fasta.with_name(fasta.name + ext)
            if not index_file.is_file():
                return False, f"Missing index file: {index_file.name}"
            if index_file.stat().st_size == 0:
                return False, f"Empty index file: {index_file.name}"
            if index_file.stat().st_mtime < fasta_mtime:
                return False, f"Outdated index file: {index_file.name} is older than {fasta.name}"

        return True, None

    @classmethod
    def parse_and_validate_fai(
        cls,
        fasta_path: Path,
        expected_contig: str,
        expected_length: int,
    ) -> Tuple[bool, Optional[str]]:
        """Validates that samtools faidx index (.fai) exists, has 1 record, and matches contig and length exactly."""
        fasta = Path(fasta_path).resolve()
        fai_file = fasta.with_name(fasta.name + ".fai")

        if not fai_file.is_file():
            return False, f"FAI index file not found: {fai_file.name}"
        if fai_file.stat().st_size == 0:
            return False, f"FAI index file is empty: {fai_file.name}"
        if fai_file.stat().st_mtime < fasta.stat().st_mtime:
            return False, f"FAI index is older than FASTA: {fai_file.name}"

        try:
            with open(fai_file, "r", encoding="utf-8") as stream:
                lines = [line.strip() for line in stream if line.strip()]

            if not lines:
                return False, "FAI index contains no records."

            if len(lines) > 1:
                return False, f"FAI index contains {len(lines)} contigs; expected exactly 1."

            fields = lines[0].split("\t")
            if len(fields) < 2:
                return False, f"Malformed FAI entry: {lines[0]}"

            fai_contig = fields[0]
            fai_length = int(fields[1])

            if fai_contig != expected_contig:
                return False, f"Contig mismatch in FAI: expected '{expected_contig}', found '{fai_contig}'"

            if fai_length != expected_length:
                return False, f"Length mismatch in FAI: parser calculated {expected_length} bp, but .fai has {fai_length} bp"

            return True, None

        except (OSError, ValueError) as exc:
            return False, f"Error validating FAI file {fai_file.name}: {exc}"

    @classmethod
    def generate_bwa_index(cls, fasta_path: Path, log_dir: Optional[Path] = None) -> None:
        """Executes 'bwa index' to generate Burrows-Wheeler Transform index files."""
        if not shutil.which("bwa"):
            raise ReferenceManagerError("Cannot generate BWA index: 'bwa' executable not found in system PATH.")

        fasta = Path(fasta_path).resolve()
        logger.info(f"Generating BWA index for: {fasta.name}")

        cmd = ["bwa", "index", str(fasta)]

        if log_dir:
            log_dir.mkdir(parents=True, exist_ok=True)
            log_file = log_dir / "bwa_index.log"
            err_file = log_dir / "bwa_index.err"
            with open(log_file, "w", encoding="utf-8") as out_f, open(err_file, "w", encoding="utf-8") as err_f:
                res = subprocess.run(cmd, stdout=out_f, stderr=err_f, text=True, check=False)
        else:
            res = subprocess.run(cmd, capture_output=True, text=True, check=False)

        if res.returncode != 0:
            err_detail = res.stderr.strip() if not log_dir else f"Check {log_dir / 'bwa_index.err'}"
            raise ReferenceManagerError(
                f"Failed to generate BWA index for {fasta.name} (exit code {res.returncode}): {err_detail}"
            )

    @classmethod
    def generate_fai_index(cls, fasta_path: Path, log_dir: Optional[Path] = None) -> None:
        """Executes 'samtools faidx' to create sequence dictionary index."""
        if not shutil.which("samtools"):
            raise ReferenceManagerError("Cannot generate FAI index: 'samtools' executable not found in system PATH.")

        fasta = Path(fasta_path).resolve()
        logger.info(f"Generating samtools .fai index for: {fasta.name}")

        cmd = ["samtools", "faidx", str(fasta)]

        if log_dir:
            log_dir.mkdir(parents=True, exist_ok=True)
            log_file = log_dir / "samtools_faidx.log"
            err_file = log_dir / "samtools_faidx.err"
            with open(log_file, "w", encoding="utf-8") as out_f, open(err_file, "w", encoding="utf-8") as err_f:
                res = subprocess.run(cmd, stdout=out_f, stderr=err_f, text=True, check=False)
        else:
            res = subprocess.run(cmd, capture_output=True, text=True, check=False)

        if res.returncode != 0:
            err_detail = res.stderr.strip() if not log_dir else f"Check {log_dir / 'samtools_faidx.err'}"
            raise ReferenceManagerError(
                f"Failed to generate samtools faidx for {fasta.name} (exit code {res.returncode}): {err_detail}"
            )

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
        """Orchestrates reference verification, cryptographic auditing, alias resolution, indexing, and cross-consistency checks.

        Args:
            fasta_path: Path to the FASTA file.
            expected_contig: Target contig identifier (defaults to rCRS 'NC_012920.1').
            strict_contig: If True, header must match expected_contig exactly. If False, allows recognized naming aliases.
            expected_length: Expected base pair count (defaults to 16569 for rCRS).
            expected_sequence_md5: Expected biological MD5 sequence hash for absolute identity verification.
            is_circular: Biological circularity metadata.
            allow_n: Whether to allow 'N' bases (default: False, requires strict ACGT).
            auto_index: Automatically generate missing BWA and FAI indices.
            log_dir: Directory to save indexing logs.

        Returns:
            ReferenceValidationResult containing populated ReferenceInfo, checklist, warnings, and errors.
        """
        path = Path(fasta_path).resolve()
        passed_checks: List[str] = []
        warnings: List[str] = []
        errors: List[str] = []

        if not path.is_file():
            raise ReferenceManagerError(f"Reference FASTA does not exist: {path}")

        # 1. Hashing dell'intero file
        file_sha = cls.compute_file_sha256(path)
        passed_checks.append("FILE_SHA256_COMPUTED")

        # 2. Parsing sequenza, calcolo MD5 biologico, verifica caratteri e singolo contig
        seq_md5, seq_id, desc, seq_len, parse_warn, parse_err = cls.parse_and_hash_sequence(
            file_path=path,
            allow_n=allow_n,
        )
        warnings.extend(parse_warn)
        errors.extend(parse_err)

        if not parse_err:
            passed_checks.append("FASTA_STRUCTURE_AND_BASES_VALID")

        # 3. Verifica corrispondenza contig (Strict vs Alias Mode)
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

        # 4. Verifica biologica di lunghezza e MD5 atteso
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

        # 5. Verifica e generazione indici BWA
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

        # 6. Verifica e generazione indice FAI
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

        # 7. Costruzione del modello ReferenceInfo sincronizzato al 100%
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
