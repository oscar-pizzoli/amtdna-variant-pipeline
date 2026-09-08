"""
Block 7 Haplogroup Orchestrator: End-to-end pipeline coordinator executing sequence
validation, dual-stream HaploGrep3 classification, multi-block audit checks, and contract emission.
"""

import json
import logging
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

from src.modules.alignment.models import ReferenceInfo
from src.modules.haplogrouping.cross_block_auditor import CrossBlockAuditor
from src.modules.haplogrouping.enums import (
    AuditCode,
    AuditSeverity,
    ContaminantAffinity,
    HaplogroupStatus,
    OverallValidationStatus,
    SequenceType,
)
from src.modules.haplogrouping.haplogrep_engine import (
    HaplogrepEngine,
    HaplogrepExecutionError,
)
from src.modules.haplogrouping.models import (
    AuditFlag,
    Block7Contract,
    Block7JSONEncoder,
    CrossBlockAuditResult,
    HaplogrepHit,
    HaplogrepResult,
    HaplogrepRunConfig,
    RawHaplogrepResult,
    SequenceIntegrityMetrics,
    ValidationResult,
    to_json_serializable,
)
from src.modules.haplogrouping.sequence_validator import SequenceValidator

logger = logging.getLogger("mtDNA_Pipeline.HaplogroupOrchestrator")


class HaplogroupOrchestrator:
    """Coordinates Block 7 execution for individual samples."""

    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.hg_cfg = config.get("haplogroup", {})
        self.val_cfg = config.get("validation", {})
        self.audit_cfg = config.get("cross_block_audit", {})
        self.ref_cfg = config.get("reference", {})

        self.run_config = HaplogrepRunConfig(
            tool=self.hg_cfg.get("tool", "haplogrep3"),
            tool_version=self.hg_cfg.get("tool_version", "unknown"),
            tree_id=self.hg_cfg.get("tree_id", "phylotree-fu-rcrs"),
            tree_version=str(self.hg_cfg.get("tree_version", "1.2")),
            metric=self.hg_cfg.get("metric", "kulczynski"),
            min_quality_score=float(self.hg_cfg.get("min_quality_score", 0.80)),
            cross_check_vcf=bool(self.hg_cfg.get("cross_check_vcf", True)),
            timeout_seconds=int(self.hg_cfg.get("timeout_seconds", 300)),
            retain_raw_output=bool(self.hg_cfg.get("retain_raw_output", True)),
        )
        self.engine = HaplogrepEngine(run_config=self.run_config)

        seq_type_str = str(self.val_cfg.get("sequence_type", "CANONICAL")).upper()
        seq_type = SequenceType.CANONICAL if seq_type_str == "CANONICAL" else SequenceType.IUPAC
        max_n_pct = float(self.val_cfg.get("max_n_percentage", 15.0))

        self.validator = SequenceValidator(
            max_n_fraction=max_n_pct / 100.0,
            sequence_type=seq_type,
        )

        self.auditor = CrossBlockAuditor(
            min_quality_score=self.run_config.min_quality_score,
            max_missing_diagnostic_mutations=int(self.audit_cfg.get("max_missing_diagnostic_mutations", 2)),
            max_private_mutations=int(self.audit_cfg.get("max_private_mutations", 5)),
            contam_affinity_threshold=float(self.audit_cfg.get("contam_affinity_threshold", 0.015)),
            subclade_contam_max_fraction=float(self.audit_cfg.get("subclade_contam_max_fraction", 0.03)),
            evaluate_contaminant_closeness=bool(self.audit_cfg.get("evaluate_contaminant_closeness", True)),
            flag_damage_in_private_mutations=bool(self.audit_cfg.get("flag_damage_in_private_mutations", True)),
        )

    def run_sample(
        self,
        sample_id: str,
        consensus_fasta: Path,
        sample_output_dir: Path,
        variants_vcf: Optional[Path] = None,
        reference_info: Optional[ReferenceInfo] = None,
        reference_fasta: Optional[Path] = None,
        expected_fasta_sha256: Optional[str] = None,
        block4_contract_path: Optional[Path] = None,
        block5_contract_path: Optional[Path] = None,
        block6_contract_path: Optional[Path] = None,
    ) -> Block7Contract:
        start_time = time.time()
        logger.info(f"[{sample_id}] Starting Block 7 Haplogrouping and Validation.")
        sample_output_dir.mkdir(parents=True, exist_ok=True)

        expected_length = reference_info.length_bp if reference_info else 16569

        # Stream reference sequence using SequenceValidator without redundant IO methods
        ref_seq: Optional[str] = None
        target_ref_path = reference_fasta or (
            Path(self.ref_cfg["fasta_path"]) if self.ref_cfg.get("fasta_path") else None
        )
        if target_ref_path and target_ref_path.is_file():
            seq, _, _ = SequenceValidator.stream_fasta_sequence(target_ref_path)
            ref_seq = seq.upper() if seq else None

        # 1. Hard-Stop on Missing Upstream Contracts
        b4_contract, err_b4 = self._load_mandatory_contract(block4_contract_path, "Block 4")
        b5_contract, err_b5 = self._load_mandatory_contract(block5_contract_path, "Block 5")
        b6_contract, err_b6 = self._load_mandatory_contract(block6_contract_path, "Block 6")

        contract_errors = [e for e in (err_b4, err_b5, err_b6) if e is not None]
        if contract_errors:
            logger.error(f"[{sample_id}] Missing mandatory upstream contract(s): {contract_errors}")
            rejection_flags = [
                AuditFlag(
                    code=AuditCode.MISSING_UPSTREAM_CONTRACT,
                    severity=AuditSeverity.CRITICAL,
                    message=err,
                    context={"sample_id": sample_id},
                )
                for err in contract_errors
            ]
            empty_metrics = SequenceIntegrityMetrics(
                sample_id=sample_id,
                sequence_type=self.validator.sequence_type,
                sequence_length=None,
                expected_length=expected_length,
                n_count=None,
                n_fraction=None,
                max_n_fraction=self.validator.max_n_fraction,
                canonical_fraction=None,
                iupac_fraction=None,
                normalization_applied=False,
                observed_sha256="",
                expected_sha256=expected_fasta_sha256,
                is_length_valid=False,
                is_n_fraction_valid=False,
                is_provenance_valid=False,
                reference_compatible=False,
                has_invalid_bases=True,
                issues=contract_errors,
            )
            return self._emit_early_rejection(
                sample_id=sample_id,
                sample_output_dir=sample_output_dir,
                metrics=empty_metrics,
                flags=rejection_flags,
                start_time=start_time,
            )

        if expected_fasta_sha256 is None and b6_contract:
            expected_fasta_sha256 = b6_contract.get("consensus_fasta_sha256")

        # 2. Pre-flight Sequence Validation
        integrity_metrics, validation_flags, valid_seq = self.validator.load_and_verify_fasta(
            fasta_path=consensus_fasta,
            expected_length=expected_length,
            expected_sha256=expected_fasta_sha256,
            sample_id=sample_id,
            reference_info=reference_info,
        )

        if not integrity_metrics.is_valid:
            logger.error(f"[{sample_id}] Sequence validation failed: {integrity_metrics.issues}")
            return self._emit_early_rejection(
                sample_id=sample_id,
                sample_output_dir=sample_output_dir,
                metrics=integrity_metrics,
                flags=validation_flags,
                start_time=start_time,
            )

        # Phantom Normalization Resolution: pass clean uppercase FASTA to tool
        input_fasta_for_engine = consensus_fasta
        if integrity_metrics.normalization_applied:
            normalized_fasta = sample_output_dir / f"{sample_id}.consensus.normalized.fasta"
            with open(normalized_fasta, "w", encoding="utf-8") as nf:
                nf.write(f">{sample_id}\n{valid_seq}\n")
            input_fasta_for_engine = normalized_fasta

        # 3. Dual-Stream HaploGrep3 Execution
        try:
            fasta_raw = self.engine.classify_fasta(
                fasta_path=input_fasta_for_engine,
                sample_id=sample_id,
                output_dir=sample_output_dir,
            )
        except HaplogrepExecutionError as exc:
            logger.error(f"[{sample_id}] HaploGrep3 FASTA classification failed: {exc}")
            flags = validation_flags + [
                AuditFlag(
                    code=AuditCode.LOW_HAPLOGROUP_CONFIDENCE,
                    severity=AuditSeverity.CRITICAL,
                    message=f"HaploGrep3 execution error: {exc}",
                )
            ]
            return self._emit_early_rejection(
                sample_id=sample_id,
                sample_output_dir=sample_output_dir,
                metrics=integrity_metrics,
                flags=flags,
                start_time=start_time,
            )

        vcf_raw: Optional[RawHaplogrepResult] = None
        if self.run_config.cross_check_vcf and variants_vcf and variants_vcf.is_file():
            try:
                vcf_raw = self.engine.classify_vcf(
                    vcf_path=variants_vcf,
                    sample_id=sample_id,
                    output_dir=sample_output_dir,
                )
            except HaplogrepExecutionError as exc:
                logger.warning(f"[{sample_id}] HaploGrep3 VCF cross-check failed (non-fatal): {exc}")
                validation_flags.append(
                    AuditFlag(
                        code=AuditCode.UPSTREAM_CONTRACT_ANOMALY,
                        severity=AuditSeverity.WARNING,
                        message=f"VCF classification failed: {exc}",
                        context={"vcf_path": str(variants_vcf)},
                    )
                )

        # 4. Integrated Cross-Block Audit
        audit_result = self.auditor.audit(
            sample_id=sample_id,
            fasta_result=fasta_raw,
            vcf_result=vcf_raw,
            block4_contract=b4_contract,
            block5_contract=b5_contract,
            block6_contract=b6_contract,
            reference_sequence=ref_seq,
        )

        all_flags = validation_flags + audit_result.flags

        final_overall_status = self._resolve_combined_status(
            seq_metrics=integrity_metrics,
            audit_status=audit_result.overall_status,
            flags=all_flags,
        )

        haplogroup_status = self._derive_haplogroup_status(
            raw_result=fasta_raw,
            audit_result=audit_result,
            integrity_metrics=integrity_metrics,
        )

        hg_result = HaplogroupResult(
            sample_id=sample_id,
            primary_haplogroup=fasta_raw.haplogroup,
            quality_score=fasta_raw.quality_score,
            status=haplogroup_status,
            top_hits=fasta_raw.top_hits,
            vcf_haplogroup=vcf_raw.haplogroup if vcf_raw else None,
            vcf_quality_score=vcf_raw.quality_score if vcf_raw else None,
            run_config=self.run_config,
            fasta_details=fasta_raw,
            vcf_details=vcf_raw,
        )

        val_result = ValidationResult(
            sample_id=sample_id,
            overall_status=final_overall_status,
            sequence_integrity=integrity_metrics,
            haplogroup_acceptable=(haplogroup_status in {
                HaplogroupStatus.ASSIGNED_CONFIDENT,
                HaplogroupStatus.ASSIGNED_WITH_WARNINGS,
            }),
            cross_block_audit=audit_result,
            flags=all_flags,
        )

        execution_seconds = round(time.time() - start_time, 2)

        # Context-Preserving Flags Serialization
        serialized_flags = [f.to_dict() for f in all_flags]

        contract = Block7Contract(
            sample_id=sample_id,
            primary_haplogroup=fasta_raw.haplogroup,
            resolved_haplogroup=audit_result.resolved_haplogroup,
            quality_score=fasta_raw.quality_score,
            quality_margin=fasta_raw.quality_margin,
            haplogroup_status=haplogroup_status,
            overall_validation_status=final_overall_status,
            tree_id=self.run_config.tree_id,
            tree_version=self.run_config.tree_version,
            metric=self.run_config.metric,
            contaminant_affinity=audit_result.contaminant_affinity,
            contamination_fraction=audit_result.contamination_fraction,
            vcf_haplogroup=vcf_raw.haplogroup if vcf_raw else None,
            vcf_concordant=audit_result.fasta_vcf_concordant if vcf_raw else None,
            subclade_shift_detected=audit_result.subclade_shift_detected,
            subclade_shift_supported=audit_result.subclade_shift_supported,
            fasta_sha256=integrity_metrics.observed_sha256,
            flags=serialized_flags,
            execution_time_seconds=execution_seconds,
        )

        self._save_artifacts(sample_output_dir, sample_id, hg_result, val_result, contract)

        logger.info(
            f"[{sample_id}] Block 7 Completed: Haplogroup={contract.primary_haplogroup} "
            f"({contract.haplogroup_status.value} / {contract.overall_validation_status.value}) "
            f"in {execution_seconds}s"
        )
        return contract

    @staticmethod
    def _load_mandatory_contract(path: Optional[Path], name: str) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
        """Safely load JSON contract, enforcing existence for mandatory upstream dependencies."""
        if path is None or not path.is_file():
            return None, f"Mandatory contract '{name}' does not exist on disk: {path}"
        try:
            with open(path, "r", encoding="utf-8") as f:
                data = json.load(f)
                if not isinstance(data, dict):
                    return None, f"Contract '{name}' is not a valid JSON dictionary: {path}"
                return data, None
        except Exception as exc:
            return None, f"Failed parsing contract '{name}' at {path}: {exc}"

    def _derive_haplogroup_status(
        self,
        raw_result: RawHaplogrepResult,
        audit_result: CrossBlockAuditResult,
        integrity_metrics: SequenceIntegrityMetrics,
    ) -> HaplogroupStatus:
        if not integrity_metrics.is_valid:
            return HaplogroupStatus.REJECTED
        if raw_result.haplogroup in {"UNKNOWN", "N/A"} or raw_result.quality_score == 0.0:
            return HaplogroupStatus.UNRESOLVED
        if raw_result.co_optimal:
            return HaplogroupStatus.AMBIGUOUS
        if not audit_result.fasta_vcf_concordant:
            return HaplogroupStatus.AMBIGUOUS
        if raw_result.quality_score < self.run_config.min_quality_score:
            return HaplogroupStatus.UNRESOLVED

        warning_codes = {
            AuditCode.EXCESS_MISSING_DATA,
            AuditCode.MISSING_DIAGNOSTIC_MUTATIONS,
            AuditCode.EXCESS_PRIVATE_MUTATIONS,
            AuditCode.CONTAMINANT_AFFINITY_WARNING,
            AuditCode.DAMAGE_SUSPECT_PRIVATE_MUTATION,
        }
        if any(f.code in warning_codes for f in audit_result.flags):
            return HaplogroupStatus.ASSIGNED_WITH_WARNINGS

        return HaplogroupStatus.ASSIGNED_CONFIDENT

    @staticmethod
    def _resolve_combined_status(
        seq_metrics: SequenceIntegrityMetrics,
        audit_status: OverallValidationStatus,
        flags: List[AuditFlag],
    ) -> OverallValidationStatus:
        if (
            not seq_metrics.is_valid
            or audit_status == OverallValidationStatus.FAIL
            or any(f.severity == AuditSeverity.CRITICAL for f in flags)
        ):
            return OverallValidationStatus.FAIL

        if (
            audit_status == OverallValidationStatus.WARN
            or any(f.severity == AuditSeverity.WARNING for f in flags)
        ):
            return OverallValidationStatus.WARN

        return OverallValidationStatus.PASS

    def _emit_early_rejection(
        self,
        sample_id: str,
        sample_output_dir: Path,
        metrics: SequenceIntegrityMetrics,
        flags: List[AuditFlag],
        start_time: float,
    ) -> Block7Contract:
        execution_seconds = round(time.time() - start_time, 2)
        serialized_flags = [f.to_dict() for f in flags]

        contract = Block7Contract(
            sample_id=sample_id,
            primary_haplogroup="REJECTED",
            resolved_haplogroup="REJECTED",
            quality_score=0.0,
            quality_margin=None,
            haplogroup_status=HaplogroupStatus.REJECTED,
            overall_validation_status=OverallValidationStatus.FAIL,
            tree_id=self.run_config.tree_id,
            tree_version=self.run_config.tree_version,
            metric=self.run_config.metric,
            contaminant_affinity=ContaminantAffinity.NOT_APPLICABLE,
            contamination_fraction=None,
            vcf_haplogroup=None,
            vcf_concordant=None,
            subclade_shift_detected=False,
            subclade_shift_supported=None,
            fasta_sha256=metrics.observed_sha256,
            flags=serialized_flags,
            execution_time_seconds=execution_seconds,
        )

        empty_raw = RawHaplogrepResult(
            sample_id=sample_id,
            input_type="fasta",
            haplogroup="REJECTED",
            quality_score=0.0,
            second_best_quality=None,
            quality_margin=None,
            co_optimal=False,
            top_hits=[],
            expected_mutations=[],
            found_mutations=[],
            missing_mutations=[],
            private_mutations=[],
            mutation_parse_warnings=[],
            haplogrep_version="none",
            tree_id=self.run_config.tree_id,
            tree_version=self.run_config.tree_version,
            metric=self.run_config.metric,
        )

        hg_result = HaplogroupResult(
            sample_id=sample_id,
            primary_haplogroup="REJECTED",
            quality_score=0.0,
            status=HaplogroupStatus.REJECTED,
            top_hits=[],
            vcf_haplogroup=None,
            vcf_quality_score=None,
            run_config=self.run_config,
            fasta_details=empty_raw,
            vcf_details=None,
        )

        val_result = ValidationResult(
            sample_id=sample_id,
            overall_status=OverallValidationStatus.FAIL,
            sequence_integrity=metrics,
            haplogroup_acceptable=False,
            cross_block_audit=CrossBlockAuditResult(
                overall_status=OverallValidationStatus.FAIL,
                fasta_vcf_concordant=False,
                fasta_haplogroup="REJECTED",
                vcf_haplogroup=None,
                resolved_haplogroup="REJECTED",
                quality_score=0.0,
                quality_margin=None,
                contaminant_affinity=ContaminantAffinity.NOT_APPLICABLE,
                contaminant_haplogroups=[],
                contamination_fraction=None,
                subclade_shift_detected=False,
                subclade_shift_supported=None,
                flags=flags,
            ),
            flags=flags,
        )

        self._save_artifacts(sample_output_dir, sample_id, hg_result, val_result, contract)
        return contract

    @classmethod
    def _save_artifacts(
        cls,
        output_dir: Path,
        sample_id: str,
        hg_result: HaplogroupResult,
        val_result: ValidationResult,
        contract: Block7Contract,
    ) -> None:
        with open(output_dir / f"{sample_id}.block7_contract.json", "w", encoding="utf-8") as f:
            json.dump(contract.to_dict(), f, indent=2, cls=Block7JSONEncoder)

        hg_payload = {
            "sample_id": hg_result.sample_id,
            "primary_haplogroup": hg_result.primary_haplogroup,
            "quality_score": hg_result.quality_score,
            "status": hg_result.status.value,
            "top_hits": [{"haplogroup": h.haplogroup, "quality": h.quality_score} for h in hg_result.top_hits],
            "vcf_haplogroup": hg_result.vcf_haplogroup,
            "vcf_quality_score": hg_result.vcf_quality_score,
            "expected_mutations": hg_result.fasta_details.expected_mutations,
            "found_mutations": hg_result.fasta_details.found_mutations,
            "missing_mutations": hg_result.fasta_details.missing_mutations,
            "private_mutations": hg_result.fasta_details.private_mutations,
            "mutation_parse_warnings": hg_result.fasta_details.mutation_parse_warnings,
        }
        with open(output_dir / f"{sample_id}.haplogroup.json", "w", encoding="utf-8") as f:
            json.dump(to_json_serializable(hg_payload), f, indent=2, cls=Block7JSONEncoder)

        val_payload = cls._serialize_validation_result(val_result)
        with open(output_dir / f"{sample_id}.validation.json", "w", encoding="utf-8") as f:
            json.dump(to_json_serializable(val_payload), f, indent=2, cls=Block7JSONEncoder)

    @staticmethod
    def _serialize_validation_result(val: ValidationResult) -> Dict[str, Any]:
        seq = val.sequence_integrity
        audit = val.cross_block_audit
        return {
            "sample_id": val.sample_id,
            "overall_status": val.overall_status.value,
            "haplogroup_acceptable": val.haplogroup_acceptable,
            "sequence_integrity": {
                "sample_id": seq.sample_id,
                "sequence_type": seq.sequence_type.value,
                "sequence_length": seq.sequence_length,
                "expected_length": seq.expected_length,
                "n_count": seq.n_count,
                "n_fraction": seq.n_fraction,
                "max_n_fraction": seq.max_n_fraction,
                "canonical_fraction": seq.canonical_fraction,
                "iupac_fraction": seq.iupac_fraction,
                "normalization_applied": seq.normalization_applied,
                "observed_sha256": seq.observed_sha256,
                "expected_sha256": seq.expected_sha256,
                "is_length_valid": seq.is_length_valid,
                "is_n_fraction_valid": seq.is_n_fraction_valid,
                "is_provenance_valid": seq.is_provenance_valid,
                "reference_compatible": seq.reference_compatible,
                "has_invalid_bases": seq.has_invalid_bases,
                "is_valid": seq.is_valid,
                "issues": seq.issues,
            },
            "cross_block_audit": {
                "overall_status": audit.overall_status.value,
                "fasta_vcf_concordant": audit.fasta_vcf_concordant,
                "fasta_haplogroup": audit.fasta_haplogroup,
                "vcf_haplogroup": audit.vcf_haplogroup,
                "resolved_haplogroup": audit.resolved_haplogroup,
                "quality_score": audit.quality_score,
                "quality_margin": audit.quality_margin,
                "contaminant_affinity": audit.contaminant_affinity.value,
                "contaminant_haplogroups": audit.contaminant_haplogroups,
                "contamination_fraction": audit.contamination_fraction,
                "subclade_shift_detected": audit.subclade_shift_detected,
                "subclade_shift_supported": audit.subclade_shift_supported,
                "private_mutations": [
                    {
                        "mutation": p.mutation,
                        "position": p.position,
                        "ref_base": p.ref_base,
                        "alt_base": p.alt_base,
                        "classification": p.classification.value,
                        "evidence_source": p.evidence_source,
                        "rationale": p.rationale,
                    }
                    for p in audit.private_mutations_assessment
                ],
            },
            "flags": [f.to_dict() for f in val.flags],
        }