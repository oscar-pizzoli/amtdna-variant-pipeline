"""Structural syntax validation, quality encoding detection, and paired-end synchronization validator for FASTQ files."""

import gzip
import logging
from pathlib import Path
from typing import Any, BinaryIO, Dict, List, Optional, Tuple

logger = logging.getLogger("mtDNA_Pipeline.FastqValidator")

# ==============================================================================
# Domain Constants for Quality Detection & FASTQ Validation
# ==============================================================================
# Phred+33 distinctive characters: '!' (33) to ':' (58), including '#', '$', '%', '+' and digits
ASCII_PHRED33_THRESHOLD: int = 59

# Phred+64 (Solexa) baseline: ';' (59) to 'h' (104)
ASCII_SOLEXA_MIN: int = 59

# Phred+64 (Illumina 1.3-1.7): '@' (64) to 'h' (104)
ASCII_ILLUMINA64_MIN: int = 64
ASCII_ILLUMINA64_MAX: int = 104

# Operational upper bounds for AdapterRemoval --qualitymax
DEFAULT_PHRED33_QUALITY_MAX: int = 93  # Accommodates extended Sanger and Q60 ']' bins
DEFAULT_PHRED64_QUALITY_MAX: int = 50  # Safe envelope for Q0-Q40 on offset 64

IO_CHUNK_SIZE: int = 262144  # 256 KB buffer for high-throughput counting


class FastqValidationError(Exception):
    """Raised when FASTQ structural checks, line parity, or mate synchronization fail."""
    pass


class FastqValidator:
    """Validates 4-line record structure, base-quality length parity, quality encoding, and PE mate synchronization."""

    @staticmethod
    def _open_stream(file_path: Path) -> BinaryIO:
        """Opens standard or gzipped FASTQ streams in binary read mode with case-insensitive extension handling."""
        path = Path(file_path).resolve()
        name_lower = path.name.lower()
        if name_lower.endswith((".gz", ".gzip", ".bgz")):
            return gzip.open(path, "rb")
        return open(path, "rb")

    @classmethod
    def _fast_count_lines(cls, file_path: Path) -> int:
        """Fast byte-chunk counting of newline characters to verify total record integrity and count."""
        path = Path(file_path).resolve()
        if not path.is_file():
            raise FastqValidationError(f"File not found for line counting: {path}")

        total_lines = 0
        try:
            with cls._open_stream(path) as stream:
                for chunk in iter(lambda: stream.read(IO_CHUNK_SIZE), b""):
                    total_lines += chunk.count(b"\n")
        except OSError as exc:
            raise FastqValidationError(f"I/O error counting lines in {path.name}: {exc}") from exc

        return total_lines

    @classmethod
    def detect_quality_encoding(
        cls,
        fastq_path: Path,
        sample_records: int = 10000,
    ) -> Tuple[int, int]:
        """Auto-detects Phred quality score offset (Phred+33 vs Phred+64) and recommended qualitymax.

        Decision Invariants:
            - Any ASCII character < 59 -> Sanger / Illumina 1.8+ (Phred+33). Default max: 93.
            - All characters in [64, 104] -> Illumina 1.3-1.7 (Phred+64). Default max: 50.
            - All characters in [59, 104] -> Solexa / early Illumina (Phred+64). Default max: 50.

        Returns:
            Tuple of (quality_base: int, quality_max: int) -> e.g. (33, 93) or (64, 50).
        """
        path = Path(fastq_path).resolve()
        if not path.is_file():
            raise FastqValidationError(f"FASTQ file not found for quality detection: {path}")

        min_ascii = 127
        max_ascii = 0
        record_idx = 0

        try:
            with cls._open_stream(path) as stream:
                while True:
                    header = stream.readline()
                    if not header:
                        break
                    seq = stream.readline()
                    plus = stream.readline()
                    qual = stream.readline()

                    if not (seq and plus and qual):
                        break

                    record_idx += 1
                    qual_strip = qual.rstrip(b"\r\n")

                    for b in qual_strip:
                        if b < min_ascii:
                            min_ascii = b
                        if b > max_ascii:
                            max_ascii = b

                    if sample_records and record_idx >= sample_records:
                        break

        except OSError as exc:
            raise FastqValidationError(f"I/O error reading quality scores from {path.name}: {exc}") from exc

        if record_idx == 0:
            logger.warning(f"[{path.name}] FASTQ is empty during quality detection. Defaulting to Phred+33 (base=33, max=93).")
            return 33, DEFAULT_PHRED33_QUALITY_MAX

        # Exhaustive and non-overlapping decision matrix
        if min_ascii < ASCII_PHRED33_THRESHOLD:
            detected_base = 33
            detected_max = DEFAULT_PHRED33_QUALITY_MAX
            encoding_name = "Sanger / Illumina 1.8+ (Phred+33)"
        elif min_ascii >= ASCII_ILLUMINA64_MIN and max_ascii <= ASCII_ILLUMINA64_MAX:
            detected_base = 64
            detected_max = DEFAULT_PHRED64_QUALITY_MAX
            encoding_name = "Illumina 1.3-1.7 (Phred+64)"
        elif min_ascii >= ASCII_SOLEXA_MIN and max_ascii <= ASCII_ILLUMINA64_MAX:
            detected_base = 64
            detected_max = DEFAULT_PHRED64_QUALITY_MAX
            encoding_name = "Solexa / Early Illumina (Phred+64)"
        else:
            detected_base = 33
            detected_max = DEFAULT_PHRED33_QUALITY_MAX
            encoding_name = "Phred+33 (Out-of-Bounds Safe Fallback)"
            logger.warning(
                f"[{path.name}] Non-standard ASCII quality distribution observed: "
                f"min={min_ascii} ('{chr(min_ascii)}'), max={max_ascii} ('{chr(max_ascii)}'). "
                f"Applying default Sanger Phred+33 fallback."
            )

        logger.info(
            f"[{path.name}] Quality Encoding Detected: {encoding_name} "
            f"(ASCII Range: {min_ascii} ['{chr(min_ascii)}'] -> {max_ascii} ['{chr(max_ascii)}']) "
            f"=> qualitybase={detected_base}, qualitymax={detected_max}"
        )
        return detected_base, detected_max

    @classmethod
    def validate_fastq_structure(
        cls,
        file_path: Path,
        max_records_to_inspect: Optional[int] = 10000,
    ) -> Dict[str, Any]:
        """Inspects structural syntax, 4-line invariant, and seq-qual length parity on a FASTQ file."""
        path = Path(file_path).resolve()
        record_idx = 0
        min_len = float("inf")
        max_len = 0

        try:
            with cls._open_stream(path) as stream:
                while True:
                    header = stream.readline()
                    if not header:
                        break

                    seq = stream.readline()
                    plus = stream.readline()
                    qual = stream.readline()

                    record_idx += 1

                    if not (seq and plus and qual):
                        raise FastqValidationError(
                            f"Truncated FASTQ record at record #{record_idx} in {path.name}."
                        )

                    if not header.startswith(b"@"):
                        raise FastqValidationError(
                            f"Malformed header at record #{record_idx} in {path.name}: expected '@', got {header[:1]!r}."
                        )
                    if not plus.startswith(b"+"):
                        raise FastqValidationError(
                            f"Malformed separator at record #{record_idx} in {path.name}: expected '+', got {plus[:1]!r}."
                        )

                    seq_strip = seq.rstrip(b"\r\n")
                    qual_strip = qual.rstrip(b"\r\n")

                    if len(seq_strip) != len(qual_strip):
                        raise FastqValidationError(
                            f"Sequence length ({len(seq_strip)}) does not match quality length ({len(qual_strip)}) "
                            f"at record #{record_idx} in {path.name}."
                        )

                    read_len = len(seq_strip)
                    if read_len < min_len:
                        min_len = read_len
                    if read_len > max_len:
                        max_len = read_len

                    if max_records_to_inspect and record_idx >= max_records_to_inspect:
                        break

        except OSError as exc:
            raise FastqValidationError(f"I/O error validating {path.name}: {exc}") from exc

        return {
            "inspected_records": record_idx,
            "min_read_length": min_len if min_len != float("inf") else 0,
            "max_read_length": max_len,
        }

    @classmethod
    def validate_paired_end_sync(
        cls,
        r1_path: Path,
        r2_path: Path,
        sample_check_records: int = 5000,
        verify_total_counts: bool = True,
    ) -> bool:
        """Validates that R1 and R2 have equal total read counts, valid structure, and matched mate identifiers."""
        r1 = Path(r1_path).resolve()
        r2 = Path(r2_path).resolve()

        if not r1.is_file() or not r2.is_file():
            raise FastqValidationError(f"Missing paired-end file: R1={r1.is_file()}, R2={r2.is_file()}")

        # 1. Whole-file line parity check (guarantees exact equal read counts across full files)
        if verify_total_counts:
            r1_lines = cls._fast_count_lines(r1)
            r2_lines = cls._fast_count_lines(r2)

            if r1_lines % 4 != 0:
                raise FastqValidationError(
                    f"R1 file {r1.name} is corrupt/truncated: total lines ({r1_lines:,}) not divisible by 4."
                )
            if r2_lines % 4 != 0:
                raise FastqValidationError(
                    f"R2 file {r2.name} is corrupt/truncated: total lines ({r2_lines:,}) not divisible by 4."
                )
            if r1_lines != r2_lines:
                raise FastqValidationError(
                    f"Paired-end read count mismatch between mate files:\n"
                    f"  {r1.name}: {r1_lines // 4:,} reads ({r1_lines:,} lines)\n"
                    f"  {r2.name}: {r2_lines // 4:,} reads ({r2_lines:,} lines)"
                )

        # 2. Single-pass structural parity and header synchronization check
        try:
            with cls._open_stream(r1) as s1, cls._open_stream(r2) as s2:
                for idx in range(1, sample_check_records + 1):
                    h1 = s1.readline()
                    h2 = s2.readline()

                    if not h1 or not h2:
                        if h1 != h2:
                            raise FastqValidationError(
                                f"Paired-end files truncated prematurely at record #{idx}: {r1.name} vs {r2.name}."
                            )
                        break

                    seq1 = s1.readline()
                    seq2 = s2.readline()
                    p1 = s1.readline()
                    p2 = s2.readline()
                    q1 = s1.readline()
                    q2 = s2.readline()

                    if not (seq1 and p1 and q1 and seq2 and p2 and q2):
                        raise FastqValidationError(f"Truncated FASTQ block encountered at record #{idx}.")

                    # Structure checks
                    if not h1.startswith(b"@") or not h2.startswith(b"@"):
                        raise FastqValidationError(f"Malformed '@' header at record #{idx}.")
                    if not p1.startswith(b"+") or not p2.startswith(b"+"):
                        raise FastqValidationError(f"Malformed '+' separator at record #{idx}.")
                    if len(seq1.rstrip(b"\r\n")) != len(q1.rstrip(b"\r\n")):
                        raise FastqValidationError(f"Length mismatch in R1 at record #{idx}.")
                    if len(seq2.rstrip(b"\r\n")) != len(q2.rstrip(b"\r\n")):
                        raise FastqValidationError(f"Length mismatch in R2 at record #{idx}.")

                    # Header matching
                    id1 = h1.split()[0]
                    id2 = h2.split()[0]

                    if id1.endswith(b"/1"):
                        id1 = id1[:-2]
                    if id2.endswith(b"/2"):
                        id2 = id2[:-2]

                    if id1 != id2:
                        raise FastqValidationError(
                            f"Mismatched paired read headers at record #{idx}:\n"
                            f"  R1: {h1.decode('ascii', 'replace').strip()}\n"
                            f"  R2: {h2.decode('ascii', 'replace').strip()}"
                        )

            logger.info(f"Paired-end structural sync validation PASSED for: {r1.name} <-> {r2.name}")
            return True

        except OSError as exc:
            raise FastqValidationError(f"I/O error during synchronization check: {exc}") from exc