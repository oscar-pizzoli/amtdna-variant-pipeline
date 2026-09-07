"""Cryptographic hashing, file metadata extraction, and high-throughput read counting."""

import gzip
import hashlib
from pathlib import Path
from typing import Optional, Tuple
from src.modules.preprocessing.models import FileFingerprint


class IntegrityError(Exception):
    """Raised when file integrity, existence, or hash computation fails."""
    pass


class IntegrityEngine:
    """Computes SHA-256 hashes, detects file formats, and extracts file metadata."""

    CHUNK_SIZE = 65536

    @staticmethod
    def compute_sha256(file_path: Path) -> str:
        path = Path(file_path).resolve()
        if not path.is_file():
            raise IntegrityError(f"Input file does not exist: {path}")

        sha256 = hashlib.sha256()
        try:
            with open(path, "rb") as stream:
                for block in iter(lambda: stream.read(IntegrityEngine.CHUNK_SIZE), b""):
                    sha256.update(block)
            return sha256.hexdigest()
        except OSError as exc:
            raise IntegrityError(f"Failed to read file for hashing {path.name}: {exc}") from exc

    @staticmethod
    def detect_format_and_compression(file_path: Path) -> Tuple[str, str]:
        path = Path(file_path).resolve()
        name_lower = path.name.lower()

        compression = "none"
        if name_lower.endswith(".gz"):
            compression = "gzip"
        elif name_lower.endswith(".bgz"):
            compression = "bgzip"

        file_format = "UNKNOWN"
        if ".fastq" in name_lower or ".fq" in name_lower:
            file_format = "FASTQ"
        elif name_lower.endswith(".bam"):
            file_format = "BAM"
            compression = "bgzip"
        elif name_lower.endswith((".fasta", ".fa", ".fna")):
            file_format = "FASTA"

        return file_format, compression

    @classmethod
    def fast_count_fastq_reads(cls, file_path: Path) -> int:
        path = Path(file_path).resolve()
        _, compression = cls.detect_format_and_compression(path)
        open_fn = gzip.open if compression == "gzip" else open

        newline_count = 0
        try:
            with open_fn(path, "rb") as stream:
                for chunk in iter(lambda: stream.read(cls.CHUNK_SIZE * 4), b""):
                    newline_count += chunk.count(b"\n")
        except gzip.BadGzipFile as exc:
            raise IntegrityError(f"Corrupted gzip archive in {path.name}: {exc}") from exc
        except OSError as exc:
            raise IntegrityError(f"I/O error reading {path.name}: {exc}") from exc

        if newline_count % 4 != 0:
            raise IntegrityError(
                f"FASTQ file {path.name} is corrupted: contains {newline_count} lines (not a multiple of 4)."
            )

        return newline_count // 4

    @classmethod
    def generate_fingerprint(
        cls,
        file_path: Path,
        mode: str = "fast",
        is_input: bool = True,
    ) -> FileFingerprint:
        """Generates a complete metadata and provenance fingerprint for a file."""
        path = Path(file_path).resolve()
        if not path.is_file():
            raise IntegrityError(f"Cannot generate fingerprint: file does not exist: {path}")

        size_bytes = path.stat().st_size
        file_format, _ = cls.detect_format_and_compression(path)

        sha256_hash: Optional[str] = None
        read_count: Optional[int] = None

        if mode == "strict" or is_input:
            sha256_hash = cls.compute_sha256(path)
            if file_format == "FASTQ":
                read_count = cls.fast_count_fastq_reads(path)

        return FileFingerprint(
            filename=path.name,
            filepath=path,
            file_size_bytes=size_bytes,
            sha256_checksum=sha256_hash,
            read_count=read_count,
            is_input=is_input,
        )
