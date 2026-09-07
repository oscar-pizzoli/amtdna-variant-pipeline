"""Robust BAM validation engine with structural integrity, contig resolution, and reference compatibility checks."""

from dataclasses import dataclass, field
import logging
from pathlib import Path
import shutil
import subprocess
from typing import Any, Dict, List, Optional, Set, Tuple
import pysam

from src.modules.preprocessing.models import SampleRecord, ValidationStatus

logger = logging.getLogger("mtDNA_Pipeline.BAMValidator")


@dataclass
class BAMValidationResult:
    """Comprehensive diagnostic evaluation of a BAM alignment file."""

    is_valid: bool = False
    status: ValidationStatus = ValidationStatus.UNTESTED
    reasons: List[str] = field(default_factory=list)
    eof_valid: bool = False
    index_present: bool = False
    index_type: Optional[str] = None  # "BAI", "CSI", or None
    declared_sort_order: str = "unknown"
    actual_sorting_verified: bool = False
    detected_mt_contig: Optional[str] = None
    target_contig_length: Optional[int] = None
    total_mapped_reads: int = 0
    mt_mapped_reads: int = 0
    reference_compatible: bool = False

    def add_failure(self, reason: str) -> None:
        """Registers a critical failure, strictly transitioning status to FAIL."""
        self.reasons.append(f"FAIL: {reason}")
        self.status = ValidationStatus.FAIL
        self.is_valid = False

    def add_warning(self, reason: str) -> None:
        """Registers a non-critical anomaly, setting status to WARN if not already in FAIL."""
        self.reasons.append(f"WARN: {reason}")
        if self.status != ValidationStatus.FAIL:
            self.status = ValidationStatus.WARN
            self.is_valid = True

    def mark_passed(self) -> None:
        """Sets status to PASS if no failures or warnings were recorded."""
        if not self.reasons:
            self.status = ValidationStatus.PASS
            self.is_valid = True

    def to_dict(self) -> Dict[str, Any]:
        """Serializes validation diagnostics to a structured dictionary."""
        return {
            "is_valid": self.is_valid,
            "status": self.status.value,
            "reasons": self.reasons,
            "eof_valid": self.eof_valid,
            "index_present": self.index_present,
            "index_type": self.index_type,
            "declared_sort_order": self.declared_sort_order,
            "actual_sorting_verified": self.actual_sorting_verified,
            "detected_mt_contig": self.detected_mt_contig,
            "target_contig_length": self.target_contig_length,
            "total_mapped_reads": self.total_mapped_reads,
            "mt_mapped_reads": self.mt_mapped_reads,
            "reference_compatible": self.reference_compatible,
        }


class BAMValidationEngine:
    """Multi-tier validator for BAM file integrity, indexation, contig resolution, and alignment metrics."""

    # Standard BGZF 28-byte End-Of-File marker
    BGZF_EOF_MARKER: bytes = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"

    # Common canonical mitochondrial contig names
    COMMON_MT_ALIASES: List[str] = [
        "NC_012920.1",
        "chrM",
        "MT",
        "chrMT",
        "rCRS",
        "gi|251831106|ref|NC_012920.1|",
        "NC_001807.4",
    ]

    @classmethod
    def check_bgzf_eof(cls, bam_path: Path) -> bool:
        """Verifies presence of the standard 28-byte BGZF EOF block to detect truncation.

        Args:
            bam_path: Target BAM file.

        Returns:
            True if EOF marker is present or verified via samtools quickcheck, False otherwise.
        """
        # Primo tentativo: samtools quickcheck se disponibile nel PATH
        if shutil.which("samtools") is not None:
            try:
                res = subprocess.run(
                    ["samtools", "quickcheck", "-v", str(bam_path)],
                    capture_output=True,
                    text=True,
                )
                return res.returncode == 0
            except OSError:
                pass

        # Fallback nativo: lettura binaria degli ultimi 28 byte
        try:
            with open(bam_path, "rb") as stream:
                stream.seek(0, 2)  # Vai a fine file
                file_size = stream.tell()
                if file_size < 28:
                    return False
                stream.seek(file_size - 28)
                tail = stream.read(28)
                return tail == cls.BGZF_EOF_MARKER
        except OSError:
            return False

    @staticmethod
    def locate_index(bam_path: Path) -> Tuple[bool, Optional[str], Optional[Path]]:
        """Checks for standard .bai or .csi indices in both naming conventions.

        Args:
            bam_path: Path to BAM file.

        Returns:
            Tuple of (index_found, index_type, index_path).
        """
        candidates = [
            (bam_path.with_suffix(bam_path.suffix + ".bai"), "BAI"),
            (bam_path.with_suffix(".bai"), "BAI"),
            (bam_path.with_suffix(bam_path.suffix + ".csi"), "CSI"),
            (bam_path.with_suffix(".csi"), "CSI"),
        ]
        for candidate_path, idx_type in candidates:
            if candidate_path.is_file():
                return True, idx_type, candidate_path
        return False, None, None

    @classmethod
    def resolve_mt_contig(
        cls,
        bam_references: List[str],
        expected_contig: Optional[str] = None,
        auto_detect: bool = True,
    ) -> Optional[str]:
        """Resolves target mitochondrial contig in BAM using exact match or authorized alias list.

        Args:
            bam_references: List of reference contigs present in BAM header.
            expected_contig: User-configured contig name (strict if auto_detect=False).
            auto_detect: Whether to search through common aliases if expected_contig is missing.

        Returns:
            Resolved contig name or None if not found.
        """
        if expected_contig and expected_contig in bam_references:
            return expected_contig

        if auto_detect:
            for alias in cls.COMMON_MT_ALIASES:
                if alias in bam_references:
                    return alias

        return None

    @classmethod
    def verify_sorting_order(
        cls,
        bam: pysam.AlignmentFile,
        target_contig: str,
        sample_reads: int = 1000,
    ) -> bool:
        """Verifies actual coordinate monotonicity for the target contig by sampling records.

        Args:
            bam: Open pysam AlignmentFile.
            target_contig: Contig to inspect.
            sample_reads: Number of contiguous reads to check.

        Returns:
            True if reads appear strictly sorted by leftmost coordinate.
        """
        try:
            last_pos = -1
            checked = 0
            for read in bam.fetch(target_contig):
                if read.is_unmapped:
                    continue
                if read.reference_start < last_pos:
                    return False
                last_pos = read.reference_start
                checked += 1
                if checked >= sample_reads:
                    break
            return True
        except Exception:
            return False

    @classmethod
    def validate_and_inspect_bam(
        cls,
        sample: SampleRecord,
        expected_contig: str = "NC_012920.1",
        reference_fasta: Optional[Path] = None,
        auto_detect_contig: bool = True,
        auto_index: bool = True,
    ) -> BAMValidationResult:
        """Executes full diagnostic validation on input BAM.

        Args:
            sample: Target SampleRecord instance.
            expected_contig: Primary target mitochondrial contig name.
            reference_fasta: Optional reference FASTA to check contig and length compatibility.
            auto_detect_contig: Allow fallback to known mitochondrial aliases.
            auto_index: Attempt index generation if missing.

        Returns:
            Populated BAMValidationResult instance.
        """
        result = BAMValidationResult()
        bam_path = sample.inputs.bam_path

        # 1. File existence
        if not bam_path or not bam_path.is_file():
            result.add_failure(f"BAM file does not exist: {bam_path}")
            sample.validation.status = result.status
            sample.validation.add_error("BAM file missing on disk.")
            return result

        # 2. EOF Truncation Check
        result.eof_valid = cls.check_bgzf_eof(bam_path)
        if not result.eof_valid:
            result.add_failure("BAM file truncated or missing BGZF EOF block (corrupted file).")

        # 3. Index Verification / Generation
        idx_found, idx_type, _ = cls.locate_index(bam_path)
        if not idx_found:
            if auto_index:
                logger.info(f"[{sample.sample_id}] Index missing. Attempting automatic indexing...")
                try:
                    pysam.index(str(bam_path))
                    idx_found, idx_type, _ = cls.locate_index(bam_path)
                    result.index_present = idx_found
                    result.index_type = idx_type
                except Exception as exc:
                    result.add_failure(f"Failed to generate index for BAM: {exc}")
            else:
                result.add_failure("BAM index (.bai/.csi) is missing and auto_index is disabled.")
        else:
            result.index_present = True
            result.index_type = idx_type

        # 4. Header, Coordinate Sorting & Contig Inspection
        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                # Header sort order declaration
                header_dict = bam.header.to_dict()
                result.declared_sort_order = header_dict.get("HD", {}).get("SO", "unknown")
                if result.declared_sort_order != "coordinate":
                    result.add_warning(
                        f"Header declared sort order is '{result.declared_sort_order}' (expected 'coordinate')."
                    )

                # Target contig resolution
                references = list(bam.references)
                lengths = list(bam.lengths)
                ref_length_map = dict(zip(references, lengths))

                resolved_contig = cls.resolve_mt_contig(
                    bam_references=references,
                    expected_contig=expected_contig,
                    auto_detect=auto_detect_contig,
                )

                if resolved_contig is None:
                    result.add_failure(
                        f"Mitochondrial contig '{expected_contig}' not found in BAM references. "
                        f"Available contigs: {references[:8]}"
                    )
                else:
                    result.detected_mt_contig = resolved_contig
                    result.target_contig_length = ref_length_map.get(resolved_contig)

                    if resolved_contig != expected_contig:
                        result.add_warning(
                            f"Target contig auto-detected as '{resolved_contig}' (configured: '{expected_contig}')."
                        )

                    # Actual coordinate sorting check
                    if result.index_present:
                        result.actual_sorting_verified = cls.verify_sorting_order(bam, resolved_contig)
                        if not result.actual_sorting_verified:
                            result.add_failure(
                                f"Actual reads on contig '{resolved_contig}' are NOT sorted by coordinate."
                            )

                    # Read Count Metrics (Global & Contig-Specific)
                    result.total_mapped_reads = bam.mapped
                    if result.index_present:
                        stats = bam.get_index_statistics()
                        contig_stats = [s for s in stats if s.contig == resolved_contig]
                        result.mt_mapped_reads = contig_stats[0].mapped if contig_stats else 0
                    else:
                        result.mt_mapped_reads = 0

                    if result.mt_mapped_reads == 0:
                        result.add_warning(
                            f"Contig '{resolved_contig}' contains 0 mapped reads (possible blank/negative control)."
                        )

                # 5. Reference FASTA Compatibility Check
                if reference_fasta and Path(reference_fasta).is_file() and result.detected_mt_contig:
                    ref_path = Path(reference_fasta).resolve()
                    try:
                        with pysam.FastaFile(str(ref_path)) as ref_fa:
                            if result.detected_mt_contig not in ref_fa.references:
                                result.add_failure(
                                    f"Reference FASTA '{ref_path.name}' does not contain BAM contig '{result.detected_mt_contig}'."
                                )
                            else:
                                ref_contig_len = ref_fa.get_reference_length(result.detected_mt_contig)
                                if result.target_contig_length and ref_contig_len != result.target_contig_length:
                                    result.add_failure(
                                        f"Contig length mismatch for '{result.detected_mt_contig}': "
                                        f"BAM has {result.target_contig_length} bp, Reference has {ref_contig_len} bp."
                                    )
                                else:
                                    result.reference_compatible = True
                    except Exception as exc:
                        result.add_warning(f"Unable to verify reference FASTA compatibility: {exc}")

        except Exception as exc:
            result.add_failure(f"Corrupted or unreadable BAM file: {exc}")

        # Final Status Finalization
        if result.status == ValidationStatus.UNTESTED:
            result.mark_passed()

        # Update SampleRecord audit fields
        sample.validation.status = result.status
        for reason in result.reasons:
            if reason.startswith("FAIL:"):
                sample.validation.add_error(reason[6:])
            elif reason.startswith("WARN:"):
                sample.validation.add_warning(reason[6:])

        # Log appropriate conclusion
        if result.status == ValidationStatus.PASS:
            logger.info(
                f"[{sample.sample_id}] BAM Validation SUCCESS: {result.mt_mapped_reads} reads on contig '{result.detected_mt_contig}'."
            )
        elif result.status == ValidationStatus.WARN:
            logger.warning(
                f"[{sample.sample_id}] BAM Validation WARNING: {len(sample.validation.warnings)} non-critical issue(s) detected."
            )
        else:
            logger.error(
                f"[{sample.sample_id}] BAM Validation FAILED: {'; '.join(sample.validation.errors)}"
            )

        return result
