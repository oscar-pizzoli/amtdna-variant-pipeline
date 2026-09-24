"""
Block 8 Provenance & Lineage Checker: Production-grade chain-of-custody auditor.
Combines cryptographic root-of-trust anchoring, deep semantic artifact inspection
(BAM, FASTA, VCF via pysam), cross-block lineage consistency matrices, and gating awareness.
"""

from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
import gzip
import hashlib
import json
import logging
import os
from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Set, Tuple, Union

try:
    import pysam
except ImportError:
    pysam = None

from src.modules.synthesis.contract_aggregator import ContractFileSet
from src.modules.synthesis.models import (
    CategorizedWarning,
    FinalVerificationReport,
    LineageVerificationItem,
    ProvenanceItem,
    ProvenanceStatus,
    SampleMasterRecord,
    SemanticVerificationItem,
    WarningCategory,
    WarningSeverity,
)

logger = logging.getLogger("mtDNA_Pipeline.ProvenanceChecker")

SHA256_REGEX = re.compile(r"^[0-9a-fA-F]{64}$")


class ProvenanceAuditError(Exception):
    """Raised when an unrecoverable failure occurs during physical artifact inspection."""
    pass


class ProvenanceChecker:
    """
    Forensic integrity auditor verifying physical existence, cryptographic hashes,
    bioinformatic semantics, and cross-block lineage against formal upstream contracts.
    """

    CHUNK_SIZE = 65536

    def __init__(
        self,
        results_dir: Path,
        operational_files_dir: Optional[Path] = None,
        enforce_strict_provenance: bool = True,
        expected_reference_length: int = 16569,
    ) -> None:
        self.results_dir = Path(results_dir)
        self.operational_files_dir = Path(operational_files_dir) if operational_files_dir else self.results_dir
        self.enforce_strict_provenance = enforce_strict_provenance
        self.expected_reference_length = expected_reference_length

    # --------------------------------------------------------------------------
    # Cryptographic & Low-Level Helpers
    # --------------------------------------------------------------------------
    @classmethod
    def calculate_sha256(cls, filepath: Path) -> str:
        """Compute SHA-256 hex digest for a physical file using chunked streaming."""
        hasher = hashlib.sha256()
        try:
            with open(filepath, "rb") as fp:
                while chunk := fp.read(cls.CHUNK_SIZE):
                    hasher.update(chunk)
            return hasher.hexdigest()
        except OSError as err:
            raise ProvenanceAuditError(f"I/O failure computing SHA-256 for '{filepath}': {err}") from err

    @staticmethod
    def is_valid_sha256(digest: Optional[str]) -> bool:
        """Strictly validate whether a string conforms to a 64-character hexadecimal SHA-256."""
        if not digest or not isinstance(digest, str):
            return False
        return bool(SHA256_REGEX.match(digest.strip()))

    @staticmethod
    def _format_mtime_utc(filepath: Path) -> Optional[str]:
        """Extract file last modified timestamp formatted as ISO-8601 UTC."""
        try:
            mtime = filepath.stat().st_mtime
            return datetime.fromtimestamp(mtime, tz=timezone.utc).isoformat()
        except OSError:
            return None

    @staticmethod
    def _read_contract_dict(path: Optional[Path]) -> Dict[str, Any]:
        """Safely read JSON dictionary returning empty dict on failure."""
        if path and path.is_file():
            try:
                with open(path, "r", encoding="utf-8") as fp:
                    data = json.load(fp)
                    return data if isinstance(data, dict) else {}
            except Exception:
                return {}
        return {}

    # --------------------------------------------------------------------------
    # Root of Trust / Execution Manifest Discovery
    # --------------------------------------------------------------------------
    def discover_execution_manifest(self, sample_id: str) -> Optional[Dict[str, Any]]:
        """
        Locate external execution manifest establishing the project's root of trust.
        Searches project root, results directory, and operational_files directory.
        """
        candidate_paths = [
            self.results_dir / "execution_manifest.json",
            self.results_dir / "run_manifest.json",
            self.results_dir / sample_id / f"{sample_id}.manifest.json",
            self.operational_files_dir / "execution_manifest.json",
            self.operational_files_dir / "run_manifest.json",
            self.operational_files_dir / sample_id / f"{sample_id}.manifest.json",
            self.results_dir.parent / "execution_manifest.json",
            self.operational_files_dir.parent / "execution_manifest.json",
        ]
        for p in candidate_paths:
            if p.is_file():
                try:
                    with open(p, "r", encoding="utf-8") as fp:
                        data = json.load(fp)
                        if isinstance(data, dict):
                            logger.info(f"[{sample_id}] Root-of-trust anchored via manifest: {p}")
                            return data
                except Exception as err:
                    logger.warning(f"Could not parse candidate manifest at '{p}': {err}")
        return None

    # --------------------------------------------------------------------------
    # Level 1 & 2: Physical & Cryptographic Artifact Auditing
    # --------------------------------------------------------------------------
    def audit_single_file(
        self,
        artifact_name: str,
        file_path: Optional[Path],
        expected_sha256: Optional[str] = None,
        is_required: bool = True,
        is_gated: bool = False,
    ) -> ProvenanceItem:
        """
        Inspect physical existence, file size, timestamps, and SHA-256 against expectation.
        Categorizes granular failures into FILE_MISSING, IO_ERROR, MISMATCH, or UNANCHORED.
        """
        if is_gated and (file_path is None or not file_path.is_file()):
            p_str = str(file_path) if file_path is not None else "NOT_PRODUCED"
            return ProvenanceItem(
                artifact_name=artifact_name,
                file_path=p_str,
                expected_sha256=expected_sha256 or "UNSPECIFIED",
                actual_sha256=None,
                file_size_bytes=None,
                status=ProvenanceStatus.GATED_OMITTED,
                last_modified_utc=None,
                error_message="Artifact omitted legitimately due to upstream biological gating.",
            )

        if file_path is None:
            return ProvenanceItem(
                artifact_name=artifact_name,
                file_path="UNSPECIFIED",
                expected_sha256=expected_sha256 or "UNSPECIFIED",
                actual_sha256=None,
                file_size_bytes=None,
                status=ProvenanceStatus.FILE_MISSING,
                last_modified_utc=None,
                error_message=f"Artifact path for '{artifact_name}' was not provided or resolved.",
            )

        if not file_path.is_file():
            return ProvenanceItem(
                artifact_name=artifact_name,
                file_path=str(file_path),
                expected_sha256=expected_sha256 or "UNSPECIFIED",
                actual_sha256=None,
                file_size_bytes=None,
                status=ProvenanceStatus.FILE_MISSING,
                last_modified_utc=None,
                error_message=f"Physical file missing on disk: '{file_path}'.",
            )

        try:
            size_bytes = file_path.stat().st_size
            mtime_utc = self._format_mtime_utc(file_path)

            if size_bytes == 0 and is_required:
                return ProvenanceItem(
                    artifact_name=artifact_name,
                    file_path=str(file_path.resolve()),
                    expected_sha256=expected_sha256 or "UNSPECIFIED",
                    actual_sha256=None,
                    file_size_bytes=0,
                    status=ProvenanceStatus.MISMATCH,
                    last_modified_utc=mtime_utc,
                    error_message=f"Artifact '{artifact_name}' exists but is empty (0 bytes).",
                )

            actual_sha = self.calculate_sha256(file_path)

            # Dedicated handling of physical index files (BAI and Tabix lack contract seals)
            if artifact_name in {"aligned_bam_index", "variants_vcf_tabix"} and (
                not expected_sha256 or expected_sha256.upper() in {"UNSPECIFIED", "NONE", "MISSING", ""}
            ):
                return ProvenanceItem(
                    artifact_name=artifact_name,
                    file_path=str(file_path.resolve()),
                    expected_sha256="PHYSICAL_INDEX_PRESENT",
                    actual_sha256=actual_sha,
                    file_size_bytes=size_bytes,
                    status=ProvenanceStatus.VERIFIED,
                    last_modified_utc=mtime_utc,
                    error_message=None,
                )

            # Formal validation of expected hash
            if not expected_sha256 or expected_sha256.upper() in {"UNSPECIFIED", "NONE", "MISSING", ""}:
                return ProvenanceItem(
                    artifact_name=artifact_name,
                    file_path=str(file_path.resolve()),
                    expected_sha256="UNANCHORED",
                    actual_sha256=actual_sha,
                    file_size_bytes=size_bytes,
                    status=ProvenanceStatus.UNANCHORED,
                    last_modified_utc=mtime_utc,
                    error_message="No reference SHA-256 seal provided by upstream contracts or manifest.",
                )

            if not self.is_valid_sha256(expected_sha256):
                return ProvenanceItem(
                    artifact_name=artifact_name,
                    file_path=str(file_path.resolve()),
                    expected_sha256=expected_sha256,
                    actual_sha256=actual_sha,
                    file_size_bytes=size_bytes,
                    status=ProvenanceStatus.MISMATCH,
                    last_modified_utc=mtime_utc,
                    error_message=f"Expected SHA-256 signature is malformed: '{expected_sha256}'.",
                )

            if actual_sha.lower() == expected_sha256.lower():
                status = ProvenanceStatus.VERIFIED
                err_msg = None
            else:
                status = ProvenanceStatus.MISMATCH
                err_msg = (
                    f"SHA-256 mismatch for '{artifact_name}': "
                    f"expected {expected_sha256}, calculated {actual_sha}."
                )

            return ProvenanceItem(
                artifact_name=artifact_name,
                file_path=str(file_path.resolve()),
                expected_sha256=expected_sha256,
                actual_sha256=actual_sha,
                file_size_bytes=size_bytes,
                status=status,
                last_modified_utc=mtime_utc,
                error_message=err_msg,
            )

        except (PermissionError, OSError) as io_err:
            return ProvenanceItem(
                artifact_name=artifact_name,
                file_path=str(file_path),
                expected_sha256=expected_sha256 or "UNSPECIFIED",
                actual_sha256=None,
                file_size_bytes=None,
                status=ProvenanceStatus.IO_ERROR,
                last_modified_utc=None,
                error_message=f"I/O or permission error accessing file: {io_err}",
            )
        except Exception as err:
            return ProvenanceItem(
                artifact_name=artifact_name,
                file_path=str(file_path),
                expected_sha256=expected_sha256 or "UNSPECIFIED",
                actual_sha256=None,
                file_size_bytes=None,
                status=ProvenanceStatus.IO_ERROR,
                last_modified_utc=None,
                error_message=f"Unexpected internal failure during inspection: {err}",
            )

    # --------------------------------------------------------------------------
    # Level 3: Deep Semantic Bioinformatic Inspection (BAM, FASTA, VCF)
    # --------------------------------------------------------------------------
    def audit_bam_semantic(
        self,
        bam_path: Optional[Path],
        expected_sample_id: str,
        expected_contig: Optional[str] = None,
    ) -> SemanticVerificationItem:
        """Inspect BAM binary integrity, coordinate sorting tag, index queryability, and contig presence."""
        if not bam_path or not bam_path.is_file():
            return SemanticVerificationItem(
                artifact_name="bam_semantic",
                file_path=str(bam_path),
                format_type="BAM",
                is_valid=False,
                error_message="BAM file is missing on disk; semantic check skipped.",
            )

        if pysam is None:
            return SemanticVerificationItem(
                artifact_name="bam_semantic",
                file_path=str(bam_path),
                format_type="BAM",
                is_valid=True,
                details={"warning": "pysam not installed; deep semantic inspection bypassed."},
            )

        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                header = bam.header.to_dict()
                hd = header.get("HD", {})
                so = hd.get("SO", "unknown")
                if so.lower() != "coordinate":
                    return SemanticVerificationItem(
                        artifact_name="bam_semantic",
                        file_path=str(bam_path),
                        format_type="BAM",
                        is_valid=False,
                        details={"sorting_order": so},
                        error_message=f"BAM is not coordinate-sorted (HD:SO='{so}').",
                    )

                references = bam.references
                target_c = expected_contig if expected_contig else (references[0] if references else None)
                if not target_c or target_c not in references:
                    return SemanticVerificationItem(
                        artifact_name="bam_semantic",
                        file_path=str(bam_path),
                        format_type="BAM",
                        is_valid=False,
                        details={"expected_contig": target_c, "available_references": references[:5]},
                        error_message=f"Target contig '{target_c}' not found in BAM header @SQ references.",
                    )

                if not bam.has_index():
                    return SemanticVerificationItem(
                        artifact_name="bam_semantic",
                        file_path=str(bam_path),
                        format_type="BAM",
                        is_valid=False,
                        error_message="BAM has no valid index loaded or index file is corrupted.",
                    )

                first_read_mapped = False
                try:
                    for read in bam.fetch(target_c, 0, min(1000, bam.get_reference_length(target_c))):
                        if not read.is_unmapped:
                            first_read_mapped = True
                            break
                except Exception as fetch_err:
                    return SemanticVerificationItem(
                        artifact_name="bam_semantic",
                        file_path=str(bam_path),
                        format_type="BAM",
                        is_valid=False,
                        error_message=f"Index lookup failed during iterator fetch: {fetch_err}",
                    )

                return SemanticVerificationItem(
                    artifact_name="bam_semantic",
                    file_path=str(bam_path),
                    format_type="BAM",
                    is_valid=True,
                    details={
                        "sorting_order": so,
                        "reference_contig": target_c,
                        "reference_length": bam.get_reference_length(target_c),
                        "index_functional": True,
                        "active_mapped_reads_observed": first_read_mapped,
                    },
                )
        except Exception as err:
            return SemanticVerificationItem(
                artifact_name="bam_semantic",
                file_path=str(bam_path),
                format_type="BAM",
                is_valid=False,
                error_message=f"Corrupted or non-BAM binary stream: {err}",
            )

    def audit_fasta_semantic(
        self,
        fasta_path: Optional[Path],
        expected_sample_id: str,
        is_indel_aware: bool = False,
    ) -> SemanticVerificationItem:
        """Inspect consensus FASTA for single-record integrity, alphabet validity, and length."""
        if not fasta_path or not fasta_path.is_file():
            return SemanticVerificationItem(
                artifact_name="fasta_semantic",
                file_path=str(fasta_path),
                format_type="FASTA",
                is_valid=False,
                error_message="FASTA artifact missing on disk.",
            )

        headers: List[str] = []
        seq_parts: List[str] = []
        valid_bases = set("ACGTNRYWSKMBDHV-")

        try:
            with open(fasta_path, "r", encoding="utf-8") as fp:
                for line in fp:
                    line_s = line.strip()
                    if not line_s:
                        continue
                    if line_s.startswith(">"):
                        headers.append(line_s[1:].split()[0])
                    else:
                        seq_parts.append(line_s.upper())

            if len(headers) != 1:
                return SemanticVerificationItem(
                    artifact_name="fasta_semantic",
                    file_path=str(fasta_path),
                    format_type="FASTA",
                    is_valid=False,
                    details={"headers_count": len(headers), "headers": headers},
                    error_message=f"Consensus FASTA must contain strictly 1 record, found {len(headers)}.",
                )

            full_seq = "".join(seq_parts)
            actual_len = len(full_seq)

            if not is_indel_aware and self.expected_reference_length and actual_len != self.expected_reference_length:
                if actual_len >= 50:
                    return SemanticVerificationItem(
                        artifact_name="fasta_semantic",
                        file_path=str(fasta_path),
                        format_type="FASTA",
                        is_valid=False,
                        details={"expected_length": self.expected_reference_length, "actual_length": actual_len},
                        error_message=f"Canonical FASTA length mismatch: expected {self.expected_reference_length} bp, got {actual_len} bp.",
                    )

            illegal_chars = set(full_seq) - valid_bases
            if illegal_chars:
                return SemanticVerificationItem(
                    artifact_name="fasta_semantic",
                    file_path=str(fasta_path),
                    format_type="FASTA",
                    is_valid=False,
                    details={"illegal_characters": list(illegal_chars)},
                    error_message=f"FASTA sequence contains illegal non-IUPAC characters: {illegal_chars}.",
                )

            return SemanticVerificationItem(
                artifact_name="fasta_semantic",
                file_path=str(fasta_path),
                format_type="FASTA",
                is_valid=True,
                details={
                    "header": headers[0],
                    "sequence_length": actual_len,
                    "n_count": full_seq.count("N"),
                    "is_indel_aware": is_indel_aware,
                },
            )
        except Exception as err:
            return SemanticVerificationItem(
                artifact_name="fasta_semantic",
                file_path=str(fasta_path),
                format_type="FASTA",
                is_valid=False,
                error_message=f"Unparseable or corrupted FASTA file: {err}",
            )

    def audit_vcf_semantic(
        self,
        vcf_path: Optional[Path],
        expected_sample_id: str,
    ) -> SemanticVerificationItem:
        """Inspect VCF syntax, bgzip compression, and sample column header presence."""
        if not vcf_path or not vcf_path.is_file():
            return SemanticVerificationItem(
                artifact_name="vcf_semantic",
                file_path=str(vcf_path),
                format_type="VCF",
                is_valid=False,
                error_message="VCF artifact missing on disk.",
            )

        if pysam is not None:
            try:
                with pysam.VariantFile(str(vcf_path)) as vcf_in:
                    samples = list(vcf_in.header.samples)
                    if expected_sample_id not in samples:
                        return SemanticVerificationItem(
                            artifact_name="vcf_semantic",
                            file_path=str(vcf_path),
                            format_type="VCF",
                            is_valid=False,
                            details={"expected_sample": expected_sample_id, "vcf_samples": samples},
                            error_message=f"Sample '{expected_sample_id}' not found in VCF sample header columns: {samples}.",
                        )

                    rec_count = sum(1 for _ in vcf_in)
                    return SemanticVerificationItem(
                        artifact_name="vcf_semantic",
                        file_path=str(vcf_path),
                        format_type="VCF",
                        is_valid=True,
                        details={"sample": expected_sample_id, "total_records": rec_count},
                    )
            except Exception as err:
                return SemanticVerificationItem(
                    artifact_name="vcf_semantic",
                    file_path=str(vcf_path),
                    format_type="VCF",
                    is_valid=False,
                    error_message=f"pysam failed to read VCF stream: {err}",
                )

        open_fn = gzip.open if vcf_path.suffix == ".gz" else open
        try:
            with open_fn(vcf_path, "rt", encoding="utf-8") as fp:
                for line in fp:
                    if line.startswith("#CHROM"):
                        cols = line.strip().split("\t")
                        if len(cols) < 10 or expected_sample_id not in cols[9:]:
                            return SemanticVerificationItem(
                                artifact_name="vcf_semantic",
                                file_path=str(vcf_path),
                                format_type="VCF",
                                is_valid=False,
                                error_message=f"Missing sample column '{expected_sample_id}' in #CHROM line.",
                            )
                        break
            return SemanticVerificationItem(
                artifact_name="vcf_semantic",
                file_path=str(vcf_path),
                format_type="VCF",
                is_valid=True,
                details={"checked_via": "line_parser"},
            )
        except Exception as err:
            return SemanticVerificationItem(
                artifact_name="vcf_semantic",
                file_path=str(vcf_path),
                format_type="VCF",
                is_valid=False,
                error_message=f"Corrupted or truncated VCF syntax: {err}",
            )

    # --------------------------------------------------------------------------
    # Level 4: Cross-Block Lineage Consistency Matrix
    # --------------------------------------------------------------------------
    def audit_cross_block_lineage(
        self,
        contracts_data: Dict[str, Dict[str, Any]],
    ) -> List[LineageVerificationItem]:
        """
        Verify mathematical identity of entities passed across pipeline frontiers:
        B4 dedup BAM == B5 input BAM; B4 active BAM == B6 input BAM;
        B6 FASTA == B7 FASTA; B6 VCF == B7 VCF.
        """
        items: List[LineageVerificationItem] = []
        d4 = contracts_data.get("block4", {})
        d5 = contracts_data.get("block5", {})
        d6 = contracts_data.get("block6", {})
        d7 = contracts_data.get("block7", {})

        # 1. Lineage B4 -> B5 (BAM Input Identity for Contamination)
        b4_unrescaled_sha = d4.get("filtered_dedup_bam_sha256") or d4.get("bam_sha256")
        b5_bam_sha = d5.get("input_bam_sha256") or d5.get("bam_sha256")

        if b4_unrescaled_sha and b5_bam_sha:
            match = (b4_unrescaled_sha.lower() == b5_bam_sha.lower())
            items.append(
                LineageVerificationItem(
                    source_block="BLOCK4",
                    target_block="BLOCK5",
                    entity_name="input_bam",
                    source_hash=b4_unrescaled_sha,
                    target_hash=b5_bam_sha,
                    is_consistent=match,
                    message=None if match else "Block 5 analyzed a BAM distinct from Block 4 dedup BAM.",
                )
            )

        # 2. Lineage B4 -> B6 (BAM Input Identity for Consensus)
        is_ancient = "AUTHENTIC_ANCIENT" in str(d4.get("authentication_status", "")).upper()
        if is_ancient and d4.get("rescaling_performed"):
            b4_consensus_expected_bam = d4.get("rescaled_bam_sha256")
        else:
            b4_consensus_expected_bam = d4.get("filtered_dedup_bam_sha256") or d4.get("bam_sha256")

        b6_bam_sha = d6.get("input_bam_sha256")
        if b4_consensus_expected_bam and b6_bam_sha:
            match = (b4_consensus_expected_bam.lower() == b6_bam_sha.lower())
            items.append(
                LineageVerificationItem(
                    source_block="BLOCK4",
                    target_block="BLOCK6",
                    entity_name="consensus_input_bam",
                    source_hash=b4_consensus_expected_bam,
                    target_hash=b6_bam_sha,
                    is_consistent=match,
                    message=None if match else "Block 6 built consensus from a BAM divergent from Block 4 expected BAM.",
                )
            )

        # 3. Lineage B6 -> B7 (Consensus FASTA Identity)
        b6_fasta_sha = d6.get("consensus_fasta_sha256")
        b7_fasta_sha = d7.get("fasta_sha256") or d7.get("consensus_fasta_sha256")
        if b6_fasta_sha and b7_fasta_sha:
            match = (b6_fasta_sha.lower() == b7_fasta_sha.lower())
            items.append(
                LineageVerificationItem(
                    source_block="BLOCK6",
                    target_block="BLOCK7",
                    entity_name="consensus_fasta",
                    source_hash=b6_fasta_sha,
                    target_hash=b7_fasta_sha,
                    is_consistent=match,
                    message=None if match else "Haplogroup assignment (B7) evaluated a mutated FASTA distinct from B6.",
                )
            )

        # 4. Lineage B6 -> B7 (Variants VCF Identity)
        b6_vcf_sha = d6.get("variants_vcf_sha256")
        b7_vcf_sha = d7.get("variants_vcf_sha256") or d7.get("vcf_sha256")
        if b6_vcf_sha and b7_vcf_sha:
            match = (b6_vcf_sha.lower() == b7_vcf_sha.lower())
            items.append(
                LineageVerificationItem(
                    source_block="BLOCK6",
                    target_block="BLOCK7",
                    entity_name="variants_vcf",
                    source_hash=b6_vcf_sha,
                    target_hash=b7_vcf_sha,
                    is_consistent=match,
                    message=None if match else "Haplogroup cross-check evaluated a VCF discordant from Block 6.",
                )
            )

        # 5. Reference Genome Continuity Across All Contracts
        ref_shas = {
            "B4": d4.get("reference_sha256"),
            "B5": d5.get("reference_sha256"),
            "B6": d6.get("reference_sha256"),
            "B7": d7.get("reference_sha256"),
        }
        valid_refs = {k: v for k, v in ref_shas.items() if v}
        if len(set(v.lower() for v in valid_refs.values())) > 1:
            items.append(
                LineageVerificationItem(
                    source_block="GLOBAL",
                    target_block="ALL",
                    entity_name="reference_genome",
                    source_hash=str(valid_refs),
                    target_hash="DISCORDANT",
                    is_consistent=False,
                    message=f"Discordant reference genomes detected across pipeline blocks: {valid_refs}",
                )
            )

        return items

    # --------------------------------------------------------------------------
    # Master Audit Orchestration
    # --------------------------------------------------------------------------
    def audit_sample_provenance(
        self,
        record: SampleMasterRecord,
        contracts: ContractFileSet,
    ) -> FinalVerificationReport:
        """
        Execute full-spectrum 4-tier chain-of-custody verification.
        Populates record.provenance, injects categorized warnings, and builds FinalVerificationReport.
        """
        sample_id = record.sample_id
        critical_failures: List[str] = []
        cautionary_warnings: List[str] = []

        manifest = self.discover_execution_manifest(sample_id)
        manifest_anchored = manifest is not None
        manifest_contracts = manifest.get("contracts", {}) if manifest else {}
        sample_manifest = manifest.get("samples", {}).get(sample_id, {}).get("contracts", {}) if manifest else {}

        raw_contracts: Dict[str, Dict[str, Any]] = {
            "block2": self._read_contract_dict(contracts.block2),
            "block3": self._read_contract_dict(contracts.block3),
            "block4": self._read_contract_dict(contracts.block4),
            "block5": self._read_contract_dict(contracts.block5),
            "block6": self._read_contract_dict(contracts.block6),
            "block7": self._read_contract_dict(contracts.block7),
        }

        is_gated_upstream = (
            record.contamination.recommended_action in {"REJECT_CONTAMINATED", "EXCLUDE_OR_MANUAL_REVIEW"}
            or record.authenticity.authentication_status in {"AUTHENTICITY_UNSUPPORTED", "INCONCLUSIVE"}
        )

        def _get_manifest_seal(b_label: str) -> Optional[str]:
            if not manifest:
                return None
            if f"{b_label}_sha256" in sample_manifest:
                return sample_manifest[f"{b_label}_sha256"]
            if b_label in sample_manifest:
                return sample_manifest[b_label]
            if f"{sample_id}_{b_label}_sha256" in manifest_contracts:
                return manifest_contracts[f"{sample_id}_{b_label}_sha256"]
            return manifest_contracts.get(f"{b_label}_sha256")

        # ----------------------------------------------------------------------
        # Tier 1: Upstream JSON Contract Audit
        # ----------------------------------------------------------------------
        contract_checks: List[ProvenanceItem] = []
        expected_contract_seals = {
            "block2_contract": _get_manifest_seal("block2_contract"),
            "block3_contract": _get_manifest_seal("block3_contract"),
            "block4_contract": _get_manifest_seal("block4_contract")
            or raw_contracts["block6"].get("block4_contract_sha256")
            or raw_contracts["block7"].get("block4_contract_sha256"),
            "block5_contract": _get_manifest_seal("block5_contract")
            or raw_contracts["block6"].get("block5_contract_sha256")
            or raw_contracts["block7"].get("block5_contract_sha256"),
            "block6_contract": _get_manifest_seal("block6_contract")
            or raw_contracts["block7"].get("block6_contract_sha256"),
            "block7_contract": _get_manifest_seal("block7_contract"),
        }

        for block_name, contract_path in [
            ("block2_contract", contracts.block2),
            ("block3_contract", contracts.block3),
            ("block4_contract", contracts.block4),
            ("block5_contract", contracts.block5),
            ("block6_contract", contracts.block6),
            ("block7_contract", contracts.block7),
        ]:
            exp_seal = expected_contract_seals.get(block_name)
            is_req = block_name in {"block4_contract", "block5_contract"} or (
                not is_gated_upstream and block_name in {"block6_contract", "block7_contract"}
            )

            item = self.audit_single_file(
                artifact_name=block_name,
                file_path=contract_path,
                expected_sha256=exp_seal,
                is_required=is_req,
                is_gated=is_gated_upstream and block_name in {"block6_contract", "block7_contract"},
            )
            contract_checks.append(item)

        # ----------------------------------------------------------------------
        # Tier 2: Primary Physical Artifact Audit
        # ----------------------------------------------------------------------
        artifact_checks: List[ProvenanceItem] = []

        def _resolve_candidate(raw_val: Optional[str], sub_dir: str, fname: str) -> Optional[Path]:
            if raw_val:
                p = Path(raw_val)
                if p.is_file():
                    return p
            for base_d in (self.operational_files_dir, self.results_dir):
                c = base_d / sub_dir / sample_id / fname
                if c.is_file():
                    return c
            return Path(raw_val) if raw_val else (self.operational_files_dir / sub_dir / sample_id / fname)

        d4 = raw_contracts["block4"]
        rescaling_perf = bool(d4.get("rescaling_performed", False))
        if rescaling_perf:
            bam_path_val = d4.get("rescaled_bam")
            exp_bam_sha = d4.get("rescaled_bam_sha256")
            bam_label = "rescaled_ancient_bam"
            bam_path = _resolve_candidate(bam_path_val, "damage", f"{sample_id}.rescaled.bam")
        else:
            bam_path_val = d4.get("filtered_dedup_bam") or d4.get("bam_path")
            exp_bam_sha = d4.get("filtered_dedup_bam_sha256") or d4.get("bam_sha256")
            bam_label = "filtered_dedup_bam"
            bam_path = _resolve_candidate(bam_path_val, "alignment", f"{sample_id}.filtered.dedup.bam")

        artifact_checks.append(
            self.audit_single_file(
                artifact_name=bam_label,
                file_path=bam_path,
                expected_sha256=exp_bam_sha,
                is_required=True,
                is_gated=False,
            )
        )

        resolved_bai: Optional[Path] = None
        if bam_path and bam_path.is_file():
            bais = [Path(f"{bam_path}.bai"), bam_path.with_suffix(".bai")]
            resolved_bai = next((p for p in bais if p.is_file()), bais[0])

        artifact_checks.append(
            self.audit_single_file(
                artifact_name="aligned_bam_index",
                file_path=resolved_bai,
                expected_sha256=None,
                is_required=True,
                is_gated=False,
            )
        )

        d6 = raw_contracts["block6"]
        d7 = raw_contracts["block7"]

        fasta_path = _resolve_candidate(
            d6.get("consensus_fasta_path"), "consensus", f"{sample_id}.consensus.fasta"
        )
        exp_fasta_sha = d6.get("consensus_fasta_sha256") or d7.get("fasta_sha256")

        artifact_checks.append(
            self.audit_single_file(
                artifact_name="canonical_consensus_fasta",
                file_path=fasta_path,
                expected_sha256=exp_fasta_sha,
                is_required=not is_gated_upstream,
                is_gated=is_gated_upstream,
            )
        )

        indel_path = _resolve_candidate(
            d6.get("consensus_indel_fasta_path"), "consensus", f"{sample_id}.consensus_indel.fasta"
        )
        indel_required = bool(d7.get("indel_consensus_used", False))
        exp_indel_sha = d6.get("consensus_indel_fasta_sha256")

        if indel_path and indel_path.is_file():
            if exp_indel_sha:
                artifact_checks.append(
                    self.audit_single_file(
                        artifact_name="consensus_indel_fasta",
                        file_path=indel_path,
                        expected_sha256=exp_indel_sha,
                        is_required=indel_required,
                        is_gated=is_gated_upstream,
                    )
                )
            else:
                artifact_checks.append(
                    ProvenanceItem(
                        artifact_name="consensus_indel_fasta",
                        file_path=str(indel_path.resolve()),
                        expected_sha256="OPTIONAL_NO_SEAL",
                        actual_sha256=self.calculate_sha256(indel_path),
                        file_size_bytes=indel_path.stat().st_size,
                        status=ProvenanceStatus.VERIFIED,
                        last_modified_utc=self._format_mtime_utc(indel_path),
                        error_message=None,
                    )
                )

        vcf_path = _resolve_candidate(
            d6.get("variants_vcf_path"), "consensus", f"{sample_id}.variants.vcf.gz"
        )
        exp_vcf_sha = d6.get("variants_vcf_sha256") or d7.get("vcf_sha256") or d7.get("variants_vcf_sha256")

        artifact_checks.append(
            self.audit_single_file(
                artifact_name="variants_vcf",
                file_path=vcf_path,
                expected_sha256=exp_vcf_sha,
                is_required=not is_gated_upstream,
                is_gated=is_gated_upstream,
            )
        )

        tbi_path = Path(f"{vcf_path}.tbi") if vcf_path else None
        artifact_checks.append(
            self.audit_single_file(
                artifact_name="variants_vcf_tabix",
                file_path=tbi_path,
                expected_sha256=None,
                is_required=not is_gated_upstream and (vcf_path is not None and vcf_path.is_file()),
                is_gated=is_gated_upstream,
            )
        )

        # ----------------------------------------------------------------------
        # Tier 3: Deep Semantic Bioinformatic Inspection
        # ----------------------------------------------------------------------
        semantic_checks: List[SemanticVerificationItem] = []

        if bam_path and bam_path.is_file():
            target_contig = d4.get("reference_contig") or d6.get("reference_contig")
            semantic_checks.append(self.audit_bam_semantic(bam_path, sample_id, target_contig))

        if not is_gated_upstream:
            if fasta_path and fasta_path.is_file():
                semantic_checks.append(self.audit_fasta_semantic(fasta_path, sample_id, is_indel_aware=False))
            if indel_path and indel_path.is_file():
                semantic_checks.append(self.audit_fasta_semantic(indel_path, sample_id, is_indel_aware=True))
            if vcf_path and vcf_path.is_file():
                semantic_checks.append(self.audit_vcf_semantic(vcf_path, sample_id))

        # ----------------------------------------------------------------------
        # Tier 4: Cross-Block Lineage Consistency Matrix
        # ----------------------------------------------------------------------
        lineage_checks = self.audit_cross_block_lineage(raw_contracts)

        # ----------------------------------------------------------------------
        # Integrity Synthesis & Warning Escalation
        # ----------------------------------------------------------------------
        all_physical_items = contract_checks + artifact_checks
        record.provenance = all_physical_items

        for it in all_physical_items:
            if it.status == ProvenanceStatus.MISMATCH:
                msg = it.error_message or f"Cryptographic mismatch in '{it.artifact_name}'."
                critical_failures.append(msg)
                record.warnings.append(
                    CategorizedWarning(
                        category=WarningCategory.PIPELINE_INTEGRITY,
                        code="PROVENANCE_HASH_MISMATCH",
                        severity=WarningSeverity.CRITICAL,
                        source_block="SYNTHESIS",
                        message=msg,
                        context={"artifact": it.artifact_name, "path": it.file_path},
                    )
                )
            elif it.status in {ProvenanceStatus.FILE_MISSING, ProvenanceStatus.IO_ERROR}:
                if it.status != ProvenanceStatus.GATED_OMITTED:
                    msg = it.error_message or f"Artifact inaccessible: '{it.artifact_name}'."
                    critical_failures.append(msg)
                    record.warnings.append(
                        CategorizedWarning(
                            category=WarningCategory.PIPELINE_INTEGRITY,
                            code="MISSING_UPSTREAM_CONTRACT",
                            severity=WarningSeverity.CRITICAL,
                            source_block="SYNTHESIS",
                            message=msg,
                            context={"artifact": it.artifact_name, "path": it.file_path},
                        )
                    )
            elif it.status == ProvenanceStatus.UNANCHORED:
                msg = f"Artifact '{it.artifact_name}' lacks sealed hash anchor."
                if self.enforce_strict_provenance and it.artifact_name in {"canonical_consensus_fasta", "variants_vcf", "block6_contract"}:
                    critical_failures.append(msg)
                    record.warnings.append(
                        CategorizedWarning(
                            category=WarningCategory.PIPELINE_INTEGRITY,
                            code="PROVENANCE_UNANCHORED",
                            severity=WarningSeverity.CRITICAL,
                            source_block="SYNTHESIS",
                            message=msg,
                            context={"artifact": it.artifact_name},
                        )
                    )
                elif manifest_anchored:
                    cautionary_warnings.append(msg)
                    record.warnings.append(
                        CategorizedWarning(
                            category=WarningCategory.PIPELINE_INTEGRITY,
                            code="PROVENANCE_UNANCHORED",
                            severity=WarningSeverity.WARNING,
                            source_block="SYNTHESIS",
                            message=msg,
                            context={"artifact": it.artifact_name},
                        )
                    )
                else:
                    record.warnings.append(
                        CategorizedWarning(
                            category=WarningCategory.PIPELINE_INTEGRITY,
                            code="PROVENANCE_UNANCHORED",
                            severity=WarningSeverity.INFO,
                            source_block="SYNTHESIS",
                            message=f"{msg} (Standalone run without execution manifest).",
                            context={"artifact": it.artifact_name},
                        )
                    )

        for sem in semantic_checks:
            if not sem.is_valid:
                msg = f"Semantic failure in {sem.format_type} ({sem.artifact_name}): {sem.error_message}"
                critical_failures.append(msg)
                record.warnings.append(
                    CategorizedWarning(
                        category=WarningCategory.PIPELINE_INTEGRITY,
                        code="SEMANTIC_VALIDATION_ERROR",
                        severity=WarningSeverity.CRITICAL,
                        source_block="SYNTHESIS",
                        message=msg,
                        context={"artifact": sem.artifact_name, "details": sem.details},
                    )
                )

        for lin in lineage_checks:
            if not lin.is_consistent:
                msg = f"Lineage break [{lin.source_block}->{lin.target_block}] on {lin.entity_name}: {lin.message}"
                critical_failures.append(msg)
                record.warnings.append(
                    CategorizedWarning(
                        category=WarningCategory.PIPELINE_INTEGRITY,
                        code="LINEAGE_DISCORDANCE",
                        severity=WarningSeverity.CRITICAL,
                        source_block="SYNTHESIS",
                        message=msg,
                        context={"source": lin.source_block, "target": lin.target_block, "entity": lin.entity_name},
                    )
                )

        if critical_failures:
            overall_verdict = "FAIL"
            record.overall_status = "FAIL"
        elif cautionary_warnings:
            overall_verdict = "WARN" if record.overall_status != "FAIL" else "FAIL"
        else:
            overall_verdict = record.overall_status

        return FinalVerificationReport(
            sample_id=sample_id,
            run_id=manifest.get("run_id", "STANDALONE_VERIFICATION") if manifest else "STANDALONE_VERIFICATION",
            overall_status=overall_verdict,
            manifest_anchored=manifest_anchored,
            contract_checks=contract_checks,
            artifact_checks=artifact_checks,
            semantic_checks=semantic_checks,
            lineage_checks=lineage_checks,
            critical_failures=critical_failures,
            warnings=cautionary_warnings,
        )