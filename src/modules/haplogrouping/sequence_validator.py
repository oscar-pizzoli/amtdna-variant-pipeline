"""
Block 7 Sequence Validator: Pre-flight genomic integrity checks, reference compatibility,
single-record strict enforcement, N-masking thresholds, and SHA-256 provenance verification.
"""

import hashlib
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from src.modules.alignment.models import ReferenceInfo
from src.modules.haplogrouping.enums import AuditCode, AuditSeverity, SequenceType
from src.modules.haplogrouping.models import AuditFlag, SequenceIntegrityMetrics

logger = logging.getLogger("mtDNA_Pipeline.SequenceValidator")


class SequenceValidator:
    """
    Verifies consensus sequence compliance prior to phylogenetic haplogrouping.
    Protects against artifact staleness, corruption, and reference drift.
    """

    CANONICAL_BASES: Set[str] = {"A", "C", "G", "T", "N"}
    IUPAC_BASES: Set[str] = {"A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "N"}

    def __init__(
        self,
        max_n_fraction: float = 0.15,
        sequence_type: SequenceType = SequenceType.CANONICAL,
        allowed_bases_override: Optional[Set[str]] = None,
    ):
        if not (0.0 <= max_n_fraction <= 1.0):
            raise ValueError(f"max_n_fraction must be in range [0.0, 1.0], observed {max_n_fraction}")

        self.max_n_fraction = max_n_fraction
        self.sequence_type = sequence_type
        self.allowed_bases = allowed_bases_override or (
            self.CANONICAL_BASES if sequence_type == SequenceType.CANONICAL else self.IUPAC_BASES
        )

    @staticmethod
    def calculate_sha256(file_path: Path, chunk_size: int = 65536) -> str:
        """Compute hex SHA-256 digest of target file."""
        hasher = hashlib.sha256()
        with open(file_path, "rb") as f:
            while chunk := f.read(chunk_size):
                hasher.update(chunk)
        return hasher.hexdigest()

    @staticmethod
    def _is_valid_hex_sha256(hash_str: str) -> bool:
        return len(hash_str) == 64 and all(c in "0123456789abcdefABCDEF" for c in hash_str)

    @classmethod
    def stream_fasta_sequence(
        cls,
        path: Path,
        max_length: int = 30000,
    ) -> Tuple[Optional[str], List[str], bool]:
        """
        Public reusable streaming FASTA parser with circuit-breaker.
        Stops immediately upon detecting multiple contigs or size > max_length.

        :return: Tuple of (assembled_sequence_or_None, headers_list, is_multi_record_triggered)
        """
        if not path.is_file():
            return None, [], False

        headers: List[str] = []
        raw_fragments: List[str] = []
        header_count = 0
        total_len = 0
        multi_record_triggered = False

        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                stripped = line.strip()
                if not stripped:
                    continue
                if stripped.startswith(">"):
                    header_count += 1
                    headers.append(stripped[1:].strip())
                    if header_count > 1:
                        multi_record_triggered = True
                        raw_fragments.clear()
                        break
                else:
                    if not multi_record_triggered:
                        raw_fragments.append(stripped)
                        total_len += len(stripped)
                        if total_len > max_length:
                            logger.error(
                                f"Sequence in '{path}' exceeded {max_length} bp. Circuit breaker triggered."
                            )
                            raw_fragments.clear()
                            return None, headers, True

        if multi_record_triggered or not headers:
            return None, headers, multi_record_triggered

        return "".join(raw_fragments), headers, False

    def load_and_verify_fasta(
        self,
        fasta_path: Path,
        expected_length: int,
        expected_sha256: Optional[str] = None,
        sample_id: Optional[str] = None,
        reference_info: Optional[ReferenceInfo] = None,
    ) -> Tuple[SequenceIntegrityMetrics, List[AuditFlag], str]:
        """
        Evaluate FASTA record integrity with streaming circuit-breaker.
        """
        if not fasta_path.is_file():
            raise FileNotFoundError(f"Consensus FASTA not found at: {fasta_path}")

        resolved_sample_id = sample_id or fasta_path.stem.split(".")[0]
        observed_sha256 = self.calculate_sha256(fasta_path)
        issues: List[str] = []
        flags: List[AuditFlag] = []

        # 1. Provenance Hash Check
        is_provenance_valid = True
        if expected_sha256 is not None and expected_sha256.strip():
            clean_expected_sha = expected_sha256.strip().lower()
            if not self._is_valid_hex_sha256(clean_expected_sha):
                raise ValueError(f"Invalid expected_sha256 format for '{resolved_sample_id}': {expected_sha256}")

            if observed_sha256.lower() != clean_expected_sha:
                is_provenance_valid = False
                msg = (
                    f"Consensus FASTA SHA-256 mismatch for '{resolved_sample_id}': "
                    f"observed {observed_sha256}, expected {clean_expected_sha}."
                )
                issues.append(msg)
                flags.append(
                    AuditFlag(
                        code=AuditCode.PROVENANCE_HASH_MISMATCH,
                        severity=AuditSeverity.CRITICAL,
                        message=msg,
                        context={"observed_sha256": observed_sha256, "expected_sha256": clean_expected_sha},
                    )
                )

        # 2. Reference Compatibility Check
        reference_compatible = True
        if reference_info is not None:
            if reference_info.length_bp != expected_length:
                reference_compatible = False
                msg = (
                    f"Reference dimension mismatch: ReferenceInfo length ({reference_info.length_bp} bp) "
                    f"differs from expected ({expected_length} bp)."
                )
                issues.append(msg)
                flags.append(
                    AuditFlag(
                        code=AuditCode.REFERENCE_INCOMPATIBLE,
                        severity=AuditSeverity.CRITICAL,
                        message=msg,
                        context={
                            "reference_contig": reference_info.contig_name,
                            "reference_length": reference_info.length_bp,
                            "expected_length": expected_length,
                        },
                    )
                )

        # 3. Stream FASTA using public circuit-breaker
        raw_sequence, headers, multi_record_triggered = self.stream_fasta_sequence(fasta_path)

        if multi_record_triggered:
            msg = f"Multi-record FASTA detected ({len(headers)}+ records). Parse aborted to prevent memory exhaustion."
            issues.append(msg)
            flags.append(
                AuditFlag(
                    code=AuditCode.MULTI_RECORD_FASTA,
                    severity=AuditSeverity.CRITICAL,
                    message=msg,
                    context={"headers": headers, "record_count": len(headers)},
                )
            )
            return (
                SequenceIntegrityMetrics(
                    sample_id=resolved_sample_id,
                    sequence_type=self.sequence_type,
                    sequence_length=None,
                    expected_length=expected_length,
                    n_count=None,
                    n_fraction=None,
                    max_n_fraction=self.max_n_fraction,
                    canonical_fraction=None,
                    iupac_fraction=None,
                    normalization_applied=False,
                    observed_sha256=observed_sha256,
                    expected_sha256=expected_sha256,
                    is_length_valid=False,
                    is_n_fraction_valid=False,
                    is_provenance_valid=is_provenance_valid,
                    reference_compatible=reference_compatible,
                    has_invalid_bases=True,
                    issues=issues,
                ),
                flags,
                "",
            )

        if not headers or raw_sequence is None:
            msg = f"Consensus FASTA lacks a valid sequence or header: {fasta_path}"
            issues.append(msg)
            flags.append(
                AuditFlag(
                    code=AuditCode.MISSING_FASTA_HEADER,
                    severity=AuditSeverity.CRITICAL,
                    message=msg,
                    context={"fasta_path": str(fasta_path)},
                )
            )
            return (
                SequenceIntegrityMetrics(
                    sample_id=resolved_sample_id,
                    sequence_type=self.sequence_type,
                    sequence_length=None,
                    expected_length=expected_length,
                    n_count=None,
                    n_fraction=None,
                    max_n_fraction=self.max_n_fraction,
                    canonical_fraction=None,
                    iupac_fraction=None,
                    normalization_applied=False,
                    observed_sha256=observed_sha256,
                    expected_sha256=expected_sha256,
                    is_length_valid=False,
                    is_n_fraction_valid=False,
                    is_provenance_valid=is_provenance_valid,
                    reference_compatible=reference_compatible,
                    has_invalid_bases=True,
                    issues=issues,
                ),
                flags,
                "",
            )

        # 4. Header Matching
        primary_header_id = headers[0].split()[0]
        if sample_id is not None and primary_header_id != sample_id:
            msg = f"FASTA header ID mismatch: found '{primary_header_id}', expected '{sample_id}'."
            issues.append(msg)
            flags.append(
                AuditFlag(
                    code=AuditCode.SAMPLE_ID_MISMATCH,
                    severity=AuditSeverity.CRITICAL,
                    message=msg,
                    context={"header_id": primary_header_id, "sample_id": sample_id},
                )
            )

        # 5. Sequence Assembly & Auditing
        normalization_applied = any(c.islower() for c in raw_sequence)
        sequence = raw_sequence.upper()
        seq_length = len(sequence)

        is_length_valid = seq_length == expected_length
        if not is_length_valid:
            msg = f"Sequence length mismatch for '{resolved_sample_id}': observed {seq_length} bp, expected {expected_length} bp."
            issues.append(msg)
            flags.append(
                AuditFlag(
                    code=AuditCode.SEQUENCE_LENGTH_MISMATCH,
                    severity=AuditSeverity.CRITICAL,
                    message=msg,
                    context={"observed_length": seq_length, "expected_length": expected_length},
                )
            )

        n_count = sequence.count("N")
        n_fraction = round(n_count / seq_length, 4) if seq_length > 0 else 1.0
        canonical_count = sum(1 for c in sequence if c in {"A", "C", "G", "T"})
        canonical_fraction = round(canonical_count / seq_length, 4) if seq_length > 0 else 0.0
        iupac_count = sum(1 for c in sequence if c in {"R", "Y", "S", "W", "K", "M"})
        iupac_fraction = round(iupac_count / seq_length, 4) if seq_length > 0 else 0.0

        is_n_fraction_valid = n_fraction <= self.max_n_fraction
        if not is_n_fraction_valid:
            msg = f"Excessive missing data for '{resolved_sample_id}': N-fraction {n_fraction:.2%} exceeds limit {self.max_n_fraction:.2%}."
            issues.append(msg)
            flags.append(
                AuditFlag(
                    code=AuditCode.EXCESS_MISSING_DATA,
                    severity=AuditSeverity.WARNING,
                    message=msg,
                    context={"n_fraction": n_fraction, "threshold": self.max_n_fraction},
                )
            )

        unique_bases = set(sequence)
        invalid_bases = unique_bases - self.allowed_bases
        has_invalid_bases = len(invalid_bases) > 0

        if has_invalid_bases:
            msg = f"Non-canonical characters for {self.sequence_type.value} detected: {sorted(invalid_bases)}"
            issues.append(msg)
            flags.append(
                AuditFlag(
                    code=AuditCode.NON_CANONICAL_CHARACTERS,
                    severity=AuditSeverity.CRITICAL,
                    message=msg,
                    context={"invalid_bases": sorted(list(invalid_bases))},
                )
            )

        metrics = SequenceIntegrityMetrics(
            sample_id=resolved_sample_id,
            sequence_type=self.sequence_type,
            sequence_length=seq_length,
            expected_length=expected_length,
            n_count=n_count,
            n_fraction=n_fraction,
            max_n_fraction=self.max_n_fraction,
            canonical_fraction=canonical_fraction,
            iupac_fraction=iupac_fraction,
            normalization_applied=normalization_applied,
            observed_sha256=observed_sha256,
            expected_sha256=expected_sha256,
            is_length_valid=is_length_valid,
            is_n_fraction_valid=is_n_fraction_valid,
            is_provenance_valid=is_provenance_valid,
            reference_compatible=reference_compatible,
            has_invalid_bases=has_invalid_bases,
            issues=issues,
        )

        return metrics, flags, sequence