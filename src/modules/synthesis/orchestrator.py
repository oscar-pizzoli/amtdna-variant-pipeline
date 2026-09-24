"""
Block 8 Synthesis Orchestrator: Master coordinator governing end-to-end specimen
synthesis, forensic provenance certification, procedure-based deliverable emission,
and cohort-wide aggregation (multi-FASTA, master QC TSV, unified variant matrices).
Enforces the "Verified Inputs Only" architectural paradigm.
"""

import csv
from datetime import datetime, timezone
from enum import Enum
import json
import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

try:
    import pysam
except ImportError:
    pysam = None

from src.modules.synthesis.contract_aggregator import (
    ContractAggregator,
    ContractFileSet,
)
from src.modules.synthesis.models import (
    MUTATION_TSV_HEADERS,
    SAMPLE_MASTER_TSV_HEADERS,
    CategorizedWarning,
    FinalVerificationReport,
    MacroComplexCoverageSummary,
    PackagingConfig,
    ProvenanceItem,
    ProvenanceStatus,
    SampleMasterRecord,
    SequencingQCSummary,
    WarningCategory,
    WarningSeverity,
    format_na,
)
from src.modules.synthesis.mutation_cataloguer import MutationCataloguer
from src.modules.synthesis.provenance_checker import ProvenanceChecker
from src.modules.synthesis.synthesis_packager import SynthesisPackager

logger = logging.getLogger("mtDNA_Pipeline.SynthesisOrchestrator")


class SynthesisOrchestrationError(Exception):
    """Raised when an unrecoverable failure occurs during sample or cohort synthesis."""
    pass


class MutationCatalogStatus(str, Enum):
    """Execution status of the mutation catalog reconciliation stream."""
    COMPLETE = "COMPLETE"
    EMPTY_VALID = "EMPTY_VALID"
    NOT_APPLICABLE = "NOT_APPLICABLE"
    FAILED = "FAILED"


class SequencingQCStatus(str, Enum):
    """Completeness status of upstream preprocessing QC metrics."""
    QC_PRESENT_VALID = "QC_PRESENT_VALID"
    QC_PARTIAL = "QC_PARTIAL"
    QC_MISSING = "QC_MISSING"
    QC_INVALID = "QC_INVALID"


class SynthesisOrchestrator:
    """
    Coordinates Block 8 execution under the 'Verified Inputs Only' principle:
    1. Aggregates upstream contracts from operational_files/ and verifies invariants.
    2. Executes 4-tier provenance audit on physical files BEFORE usage.
    3. Resolves only certified BAM and VCF artifacts for downstream computation.
    4. Computes macro-complex coverages and reconciles variant catalogs.
    5. Evaluates releaseability and packages deliverables into results/<category>/<sample_id>/.
    6. Aggregates cohort matrices (multi-FASTA, master TSVs, horizontal mutation summary).
    """

    DEFAULT_INTERVALS = {
        "dloop": [(16024, 16569), (1, 576)],
        "rrna": [(648, 1601), (1671, 3229)],
        "cox": [(5904, 7445), (7586, 8269), (9207, 9990)],
        "nd": [
            (3307, 4262),
            (4470, 5511),
            (10059, 10404),
            (10470, 10766),
            (10760, 12137),
            (12337, 14148),
            (14149, 14673),
        ],
        "cyb": [(14747, 15887)],
    }

    def __init__(
        self,
        results_dir: Path,
        reference_info: Any,
        operational_files_dir: Optional[Path] = None,
        reference_sequence: Optional[str] = None,
        config: Optional[Dict[str, Any]] = None,
    ) -> None:
        self.results_dir = Path(results_dir)
        self.operational_files_dir = Path(operational_files_dir) if operational_files_dir else self.results_dir
        self.config = config or {}
        self.reference_sequence = reference_sequence

        # Base procedural directories shortcuts for cohort aggregation and deliverables
        self.dir_master = self.results_dir / "master_summary"
        self.dir_fasta = self.results_dir / "fasta"
        self.dir_vcf = self.results_dir / "vcf"
        self.dir_haplo = self.results_dir / "haplogroups"
        self.dir_contam = self.results_dir / "contamination"
        self.dir_muts = self.results_dir / "mutations"
        self.dir_qc = self.results_dir / "qc_sequencing"
        self.dir_reports = self.results_dir / "reports_markdown"

        for d in [
            self.dir_master,
            self.dir_fasta,
            self.dir_vcf,
            self.dir_haplo,
            self.dir_contam,
            self.dir_muts,
            self.dir_qc,
            self.dir_reports,
        ]:
            d.mkdir(parents=True, exist_ok=True)

        # Preserve complete ReferenceInfo object without degradation
        self.reference_info = reference_info
        if hasattr(reference_info, "length_bp"):
            self.reference_length = reference_info.length_bp
        elif hasattr(reference_info, "length"):
            self.reference_length = reference_info.length
        elif isinstance(reference_info, int):
            self.reference_length = reference_info
        else:
            self.reference_length = 16569

        if hasattr(reference_info, "contig_name"):
            self.contig_name = str(reference_info.contig_name)
        elif hasattr(reference_info, "contig"):
            self.contig_name = str(reference_info.contig)
        else:
            self.contig_name = self.config.get("reference", {}).get("expected_contig", "NC_012920.1")

        self.pipeline_version = self.config.get("project", {}).get("version", "1.2.0")

        synth_cfg = self.config.get("synthesis", {})
        self.min_vaf_reportable = synth_cfg.get("min_vaf_reportable", 0.10)
        self.enforce_strict_provenance = synth_cfg.get("enforce_strict_provenance", True)

        consensus_cfg = self.config.get("consensus", {})
        self.min_base_quality = consensus_cfg.get("min_base_quality", 30)

        pkg_cfg = PackagingConfig(
            reports_dir=self.dir_master,
            emit_tsv_summary=synth_cfg.get("deliverables", {}).get("emit_tsv_summary", True),
            emit_mutations_catalog=synth_cfg.get("deliverables", {}).get("emit_mutations_catalog", True),
            emit_json_master=synth_cfg.get("deliverables", {}).get("emit_json_master", True),
            emit_markdown_sheets=synth_cfg.get("deliverables", {}).get("emit_markdown_sheets", True),
            enforce_strict_provenance=self.enforce_strict_provenance,
        )

        # Satellite modules reading from operational_files_dir and emitting to results_dir
        self.aggregator = ContractAggregator(self.operational_files_dir)
        self.cataloguer = MutationCataloguer(
            results_dir=self.operational_files_dir,
            reference_info=self.reference_info,
            reference_sequence=self.reference_sequence,
        )
        self.checker = ProvenanceChecker(
            results_dir=self.operational_files_dir,
            enforce_strict_provenance=self.enforce_strict_provenance,
            expected_reference_length=self.reference_length,
        )
        self.packager = SynthesisPackager(
            results_dir=self.results_dir,
            operational_files_dir=self.operational_files_dir,
            config=pkg_cfg,
        )

    # --------------------------------------------------------------------------
    # Centralized Gating Predicate (Single Source of Truth)
    # --------------------------------------------------------------------------
    @staticmethod
    def is_sample_gated(record: SampleMasterRecord) -> bool:
        """
        Evaluate if sample was legitimately excluded from downstream sequence generation.
        Strictly distinguishes authentic modern DNA (AUTHENTICITY_UNSUPPORTED) from failure.
        """
        contam_gated = record.contamination.recommended_action in {
            "REJECT_CONTAMINATED",
            "EXCLUDE_OR_MANUAL_REVIEW",
        }
        auth_gated = record.authenticity.authentication_status in {
            "INCONCLUSIVE",
            "NOISY_OR_INSUFFICIENT",
        }
        return bool(contam_gated or auth_gated)

    # --------------------------------------------------------------------------
    # Verified Macro-Complex Coverage Calculation (Fail-Closed Contig Resolution)
    # --------------------------------------------------------------------------
    def _compute_macro_complex_coverages(
        self, bam_path: Optional[Path]
    ) -> MacroComplexCoverageSummary:
        """
        Compute mean callable depth and breadth across mitochondrial complexes (Q >= min_base_quality).
        Fail-closed: rejects non-mitochondrial contigs without falling back to refs[0].
        """
        summary = MacroComplexCoverageSummary()
        if not bam_path or not bam_path.is_file() or pysam is None:
            return summary

        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                refs = bam.references
                if self.contig_name not in refs:
                    raise SynthesisOrchestrationError(
                        f"Expected mitochondrial contig '{self.contig_name}' is absent from BAM "
                        f"header @SQ references: {refs[:5]}. Refusing fallback."
                    )

                ref_len = bam.get_reference_length(self.contig_name)

                # Compute coverage with consensus official Phred threshold (Block 6)
                a_counts, c_counts, g_counts, t_counts = bam.count_coverage(
                    self.contig_name, 0, ref_len, quality_threshold=self.min_base_quality
                )
                depths = [
                    a_counts[i] + c_counts[i] + g_counts[i] + t_counts[i]
                    for i in range(ref_len)
                ]

                def _eval_intervals(intervals: List[Tuple[int, int]]) -> Tuple[float, float]:
                    region_depths: List[int] = []
                    for start, end in intervals:
                        for p in range(start, end + 1):
                            if 1 <= p <= ref_len:
                                region_depths.append(depths[p - 1])
                    if not region_depths:
                        return 0.0, 0.0
                    mean_d = sum(region_depths) / len(region_depths)
                    b1x = sum(1 for d in region_depths if d >= 1) / len(region_depths)
                    return round(mean_d, 2), round(b1x, 4)

                summary.cov_dloop_mean, summary.cov_dloop_1x = _eval_intervals(self.DEFAULT_INTERVALS["dloop"])
                summary.cov_rrna_mean, summary.cov_rrna_1x = _eval_intervals(self.DEFAULT_INTERVALS["rrna"])
                summary.cov_cox_mean, summary.cov_cox_1x = _eval_intervals(self.DEFAULT_INTERVALS["cox"])
                summary.cov_nd_mean, summary.cov_nd_1x = _eval_intervals(self.DEFAULT_INTERVALS["nd"])
                summary.cov_cyb_mean, summary.cov_cyb_1x = _eval_intervals(self.DEFAULT_INTERVALS["cyb"])

        except Exception as err:
            logger.error(f"Failed calculating macro-complex coverages on '{bam_path}': {err}")
            raise

        return summary

    def _extract_sequencing_qc(
        self, sample_id: str, contracts: ContractFileSet
    ) -> Tuple[SequencingQCSummary, SequencingQCStatus]:
        """Extract sequencing quality metrics tracking data completeness explicitly."""
        summary = SequencingQCSummary()
        d_b2: Dict[str, Any] = {}
        if contracts.block2 and contracts.block2.is_file():
            try:
                with open(contracts.block2, "r", encoding="utf-8") as fp:
                    d_b2 = json.load(fp)
            except Exception:
                return summary, SequencingQCStatus.QC_INVALID
        else:
            return summary, SequencingQCStatus.QC_MISSING

        summary.raw_reads_total = d_b2.get("raw_reads_total")
        summary.clean_reads_total = d_b2.get("clean_reads_total") or d_b2.get("retained_reads_total")
        summary.retention_fraction = d_b2.get("retention_fraction")
        summary.collapsed_fraction = d_b2.get("collapsed_fraction")
        summary.adapter_trimmed_reads_pct = d_b2.get("adapter_trimmed_reads_pct")
        summary.q20_rate_post = d_b2.get("q20_rate_post") or d_b2.get("q20_rate")
        summary.q30_rate_post = d_b2.get("q30_rate_post") or d_b2.get("q30_rate")
        summary.gc_content_pct = d_b2.get("gc_content_pct")
        summary.mean_read_length_pre = d_b2.get("mean_read_length_pre")
        summary.mean_read_length_post = d_b2.get("mean_read_length_post")

        qc_json = self.operational_files_dir / "preprocessing" / sample_id / f"{sample_id}.qc.json"
        if qc_json.is_file():
            try:
                with open(qc_json, "r", encoding="utf-8") as fp:
                    qc_data = json.load(fp)
                    if "per_base_quality" in qc_data:
                        summary.per_base_quality_cycles = qc_data["per_base_quality"]
                    if "length_distribution" in qc_data:
                        summary.read_length_distribution = {
                            int(k): int(v) for k, v in qc_data["length_distribution"].items()
                        }
            except Exception as err:
                logger.debug(f"Companion QC JSON unreadable for '{sample_id}': {err}")

        is_complete = all(
            v is not None
            for v in (summary.raw_reads_total, summary.clean_reads_total, summary.retention_fraction)
        )
        status = SequencingQCStatus.QC_PRESENT_VALID if is_complete else SequencingQCStatus.QC_PARTIAL
        return summary, status

    # --------------------------------------------------------------------------
    # Sample-Level Orchestration (Verified Inputs First)
    # --------------------------------------------------------------------------
    def process_sample(
        self, sample_id: str
    ) -> Tuple[SampleMasterRecord, FinalVerificationReport, Dict[str, Path]]:
        """
        Execute full synthesis and verification pipeline under the Verified Inputs paradigm:
        1. Ingest contracts from operational_files/ and evaluate logical invariants early.
        2. Execute 4-Tier Provenance Audit to verify physical artifact integrity.
        3. Dispatch downstream computation strictly using verified inputs.
        4. Reconcile mutation catalog and check releaseability.
        5. Package deliverables atomically into results/<category>/<sample_id>/.
        """
        logger.info(f"[{sample_id}] === Starting Block 8 Synthesis and Verification ===")

        # 1. Ingest formal upstream contracts from operational_files_dir
        contracts = self.aggregator.locate_contracts(sample_id)
        record = self.aggregator.aggregate_sample(sample_id, contracts=contracts)

        # 2. Early validation of logical invariants
        violations = record.validate_invariants()
        for v in violations:
            if not any(w.code == v.code and v.message in w.message for w in record.warnings):
                record.warnings.append(
                    CategorizedWarning(
                        category=WarningCategory.PIPELINE_INTEGRITY,
                        code=v.code,
                        severity=v.severity,
                        source_block="SYNTHESIS",
                        message=f"Cross-block invariant anomaly: {v.message}",
                        context=v.context,
                    )
                )
        if any(w.severity == WarningSeverity.CRITICAL for w in record.warnings):
            record.overall_status = "FAIL"

        # 3. 4-TIER PROVENANCE AUDIT (Executed BEFORE file consumption)
        verification_report = self.checker.audit_sample_provenance(record, contracts)
        record.provenance = verification_report.artifact_checks

        tier1_status = "FAIL" if any(c.status == ProvenanceStatus.MISMATCH for c in verification_report.contract_checks) else "PASS"
        tier2_status = "FAIL" if any(a.status in {ProvenanceStatus.MISMATCH, ProvenanceStatus.FILE_MISSING} for a in verification_report.artifact_checks if a.status != ProvenanceStatus.GATED_OMITTED) else "PASS"
        tier3_status = "FAIL" if any(not s.is_valid for s in verification_report.semantic_checks) else "PASS"
        tier4_status = "FAIL" if any(not l.is_consistent for l in verification_report.lineage_checks) else "PASS"

        logger.info(
            f"[{sample_id}] 4-Tier Provenance Audit: "
            f"T1(Contracts)={tier1_status} | T2(Artifacts)={tier2_status} | "
            f"T3(Semantics)={tier3_status} | T4(Lineage)={tier4_status} | "
            f"Overall={verification_report.overall_status}"
        )

        # 4. RESOLUTION OF VERIFIED INPUTS (Verified Inputs Only)
        is_gated = self.is_sample_gated(record)
        prov_map = {p.artifact_name: p for p in verification_report.artifact_checks}

        # Resolve verified BAM with protocol awareness: rescaled if ancient and rescaling was executed, else dedup BAM
        if record.authenticity.rescaling_performed:
            bam_item = prov_map.get("rescaled_ancient_bam")
            if (not bam_item or bam_item.status not in {ProvenanceStatus.VERIFIED, ProvenanceStatus.UNANCHORED}) and not self.enforce_strict_provenance:
                bam_item = prov_map.get("filtered_dedup_bam")
        else:
            bam_item = prov_map.get("filtered_dedup_bam")

        verified_bam_path: Optional[Path] = None
        if bam_item and bam_item.status in {ProvenanceStatus.VERIFIED, ProvenanceStatus.UNANCHORED}:
            verified_bam_path = Path(bam_item.file_path)
            record.macro_complexes = self._compute_macro_complex_coverages(verified_bam_path)
        else:
            logger.warning(
                f"[{sample_id}] Verified BAM unavailable (status={getattr(bam_item, 'status', 'MISSING')}). "
                f"Macro-complex coverage calculations skipped."
            )
            record.macro_complexes = MacroComplexCoverageSummary()

        # Ingest Sequencing QC metrics
        record.sequencing_qc, qc_status = self._extract_sequencing_qc(sample_id, contracts)

        # Resolve verified VCF and reconcile mutation catalog
        vcf_item = prov_map.get("variants_vcf")
        vcf_verified = vcf_item and vcf_item.status in {ProvenanceStatus.VERIFIED, ProvenanceStatus.UNANCHORED}

        if is_gated:
            record.mutations = []
            mutation_catalog_status = MutationCatalogStatus.NOT_APPLICABLE
            logger.info(f"[{sample_id}] Upstream gated sample: mutation catalog N/A.")
        elif not vcf_verified and self.enforce_strict_provenance:
            record.mutations = []
            mutation_catalog_status = MutationCatalogStatus.FAILED
            record.warnings.append(
                CategorizedWarning(
                    category=WarningCategory.PIPELINE_INTEGRITY,
                    code="UNTRUSTED_VCF_BLOCKED",
                    severity=WarningSeverity.CRITICAL,
                    source_block="SYNTHESIS",
                    message="Mutation cataloging halted: variants VCF failed cryptographic provenance audit.",
                    context={"vcf_status": getattr(vcf_item, "status", "MISSING")},
                )
            )
            record.overall_status = "FAIL"
        else:
            try:
                record.mutations = self.cataloguer.build_catalog(
                    sample_id=sample_id,
                    require_vcf=True,
                    min_vaf_reportable=self.min_vaf_reportable,
                )
                mutation_catalog_status = (
                    MutationCatalogStatus.COMPLETE if record.mutations else MutationCatalogStatus.EMPTY_VALID
                )
            except Exception as err:
                logger.error(f"[{sample_id}] Error compiling mutation catalog: {err}")
                record.mutations = []
                mutation_catalog_status = MutationCatalogStatus.FAILED
                record.overall_status = "FAIL"

        # 5. RELEASEABILITY GATE & PROCEDURAL PACKAGING
        is_releasable = (record.overall_status != "FAIL")
        if not is_releasable:
            logger.warning(
                f"[{sample_id}] Releaseability gate: sample evaluated as FAIL. "
                f"Deliverables will contain explicit non-compliance warnings."
            )

        emitted_artifacts = self.packager.package_sample(record)

        logger.info(
            f"[{sample_id}] Synthesis complete: Status={record.overall_status} "
            f"(Releasable={is_releasable}), Catalog={mutation_catalog_status.value}, "
            f"QC={qc_status.value}."
        )

        return record, verification_report, emitted_artifacts

    # --------------------------------------------------------------------------
    # Cohort-Level Aggregation Pipeline (Stale-Checking & Provenance Anchored)
    # --------------------------------------------------------------------------
    def _is_sample_stale(self, sample_id: str, contracts: ContractFileSet) -> bool:
        """Check if an on-disk sample report is outdated relative to its source contracts."""
        record_json = self.dir_master / f"{sample_id}.master.json"
        record_tsv = self.dir_master / f"{sample_id}.summary.tsv"
        if not (record_json.is_file() and record_tsv.is_file()):
            return True

        try:
            mtime_json = record_json.stat().st_mtime
            for b_name in ("block2", "block3", "block4", "block5", "block6", "block7"):
                p: Optional[Path] = getattr(contracts, b_name)
                if p and p.is_file() and p.stat().st_mtime > mtime_json:
                    logger.info(f"[{sample_id}] Contract '{b_name}' has been updated. Invalidating cache.")
                    return True
            return False
        except OSError:
            return True

    def run_cohort(self, sample_ids: Optional[List[str]] = None) -> Dict[str, Path]:
        """
        Aggregate all cohort deliverables across sample-partitioned folders:
        - results/master_summary/cohort_qc_summary.tsv (48 columns)
        - results/master_summary/cohort_master.json
        - results/fasta/cohort_canonical.fasta
        - results/fasta/cohort_heteroplasmy.fasta
        - results/haplogroups/cohort_haplogroups_summary.tsv
        - results/mutations/cohort_mutation_catalog.tsv (Tidy, 1 row per mutation)
        - results/mutations/cohort_mutations_summary.tsv (Horizontal: sample_id | hg | count | list)
        """
        if not sample_ids:
            cohort_cfg_samples = self.config.get("cohort", {}).get("samples")
            if cohort_cfg_samples and isinstance(cohort_cfg_samples, list):
                sample_ids = sorted(list(set(str(s) for s in cohort_cfg_samples)))
            else:
                discovered: Set[str] = set()
                for cand_dir in [self.operational_files_dir / "consensus", self.dir_master]:
                    if cand_dir.is_dir():
                        for item in cand_dir.iterdir():
                            if item.is_dir() and (self.operational_files_dir / "consensus" / item.name / f"{item.name}.block6_contract.json").is_file():
                                discovered.add(item.name)
                            elif item.is_file() and item.name.endswith(".summary.tsv"):
                                discovered.add(item.name.replace(".summary.tsv", ""))
                sample_ids = sorted(list(discovered))

        if not sample_ids:
            logger.warning("No valid samples identified for cohort aggregation.")
            return {}

        logger.info(f"=== Starting Cohort Aggregation: {len(sample_ids)} samples ===")
        cohort_deliverables: Dict[str, Path] = {}

        master_rows: List[List[str]] = []
        master_json_records: Dict[str, Any] = {}
        all_mutations_rows: List[List[str]] = []
        haplo_summary_rows: List[List[str]] = []

        fasta_canonical_out = self.dir_fasta / "cohort_canonical.fasta"
        fasta_het_out = self.dir_fasta / "cohort_heteroplasmy.fasta"

        tmp_can_fasta = fasta_canonical_out.with_name(f".tmp_{fasta_canonical_out.name}_{os.getpid()}")
        tmp_het_fasta = fasta_het_out.with_name(f".tmp_{fasta_het_out.name}_{os.getpid()}")

        try:
            with open(tmp_can_fasta, "w", encoding="utf-8") as fp_can, \
                 open(tmp_het_fasta, "w", encoding="utf-8") as fp_het:

                for sid in sample_ids:
                    contracts = self.aggregator.locate_contracts(sid)
                    record_json_path = self.dir_master / f"{sid}.master.json"
                    record_summary_path = self.dir_master / f"{sid}.summary.tsv"

                    if self._is_sample_stale(sid, contracts):
                        self.process_sample(sid)

                    with open(record_json_path, "r", encoding="utf-8") as jfp:
                        rec_dict = json.load(jfp)

                    with open(record_summary_path, "r", encoding="utf-8") as sfp:
                        reader = csv.reader(sfp, delimiter="\t")
                        _ = next(reader, None)
                        row = next(reader, None)
                        if row:
                            master_rows.append(row)

                    master_json_records[sid] = rec_dict

                    # Read from sample-partitioned mutations folder: results/mutations/<sample_id>/
                    sample_mut_tsv = self.dir_muts / sid / f"{sid}.mutation_catalog.tsv"
                    if sample_mut_tsv.is_file():
                        with open(sample_mut_tsv, "r", encoding="utf-8") as mfp:
                            m_reader = csv.reader(mfp, delimiter="\t")
                            _ = next(m_reader, None)
                            for m_row in m_reader:
                                all_mutations_rows.append(m_row)

                    warnings_raw = rec_dict.get("warnings", [])
                    active_warnings: List[str] = []
                    for w in warnings_raw:
                        if isinstance(w, dict):
                            active_warnings.append(str(w.get("code", "UNKNOWN")))
                        elif isinstance(w, str):
                            active_warnings.append(w)
                        else:
                            active_warnings.append(str(getattr(w, "code", w)))
                    warn_str = "; ".join(active_warnings) if active_warnings else "CLEAN"

                    hg_info = rec_dict.get("haplogroup", {})
                    haplo_summary_rows.append([
                        sid,
                        str(hg_info.get("primary_haplogroup", "UNKNOWN")),
                        format_na(hg_info.get("quality_score"), "{:.4f}"),
                        format_na(hg_info.get("quality_margin"), "{:.4f}"),
                        str(hg_info.get("haplogroup_status", "UNRESOLVED")),
                        str(rec_dict.get("overall_status", "UNKNOWN")),
                        str(rec_dict.get("mutations_count", 0)),
                        warn_str,
                    ])

                    # Read from sample-partitioned fasta folder: results/fasta/<sample_id>/
                    can_fasta = self.dir_fasta / sid / f"{sid}.canonical.fasta"
                    if can_fasta.is_file():
                        with open(can_fasta, "r", encoding="utf-8") as sfp:
                            fp_can.write(sfp.read().strip() + "\n")

                    het_fasta = self.dir_fasta / sid / f"{sid}.heteroplasmy.fasta"
                    if het_fasta.is_file():
                        with open(het_fasta, "r", encoding="utf-8") as sfp:
                            fp_het.write(sfp.read().strip() + "\n")

                fp_can.flush()
                os.fsync(fp_can.fileno())
                fp_het.flush()
                os.fsync(fp_het.fileno())

            tmp_can_fasta.replace(fasta_canonical_out)
            tmp_het_fasta.replace(fasta_het_out)
            cohort_deliverables["cohort_canonical_fasta"] = fasta_canonical_out
            cohort_deliverables["cohort_heteroplasmy_fasta"] = fasta_het_out

        except Exception as err:
            for p in (tmp_can_fasta, tmp_het_fasta):
                if p.is_file():
                    p.unlink(missing_ok=True)
            raise SynthesisOrchestrationError(f"Failed generating cohort multi-FASTA files: {err}") from err

        # 2. Emit 48-column Master QC Summary TSV
        cohort_summary_tsv = self.dir_master / "cohort_qc_summary.tsv"
        self._write_atomic_tsv(cohort_summary_tsv, SAMPLE_MASTER_TSV_HEADERS, master_rows)
        cohort_deliverables["cohort_qc_summary"] = cohort_summary_tsv

        # 3. Emit Cohort Master JSON
        cohort_master_json = self.dir_master / "cohort_master.json"
        cohort_dump = {
            "project_name": self.config.get("project", {}).get("name", "ancient-modern-mtdna-pipeline"),
            "pipeline_version": self.pipeline_version,
            "reference_contig": self.contig_name,
            "reference_length": self.reference_length,
            "cohort_size": len(sample_ids),
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "samples": master_json_records,
        }
        self._write_atomic_text(cohort_master_json, json.dumps(cohort_dump, indent=2) + "\n")
        cohort_deliverables["cohort_master_json"] = cohort_master_json

        # 4. Emit Cohort Haplogroup Summary Table
        haplo_headers = [
            "sample_id",
            "primary_haplogroup",
            "quality_score",
            "quality_margin",
            "haplogroup_status",
            "overall_status",
            "mutations_count",
            "active_warnings",
        ]
        cohort_haplo_tsv = self.dir_haplo / "cohort_haplogroups_summary.tsv"
        self._write_atomic_tsv(cohort_haplo_tsv, haplo_headers, haplo_summary_rows)
        cohort_deliverables["cohort_haplogroups_summary"] = cohort_haplo_tsv

        # 5. Emit Cohort Mutation Catalog (Tidy vertical TSV)
        cohort_muts_tsv = self.dir_muts / "cohort_mutation_catalog.tsv"
        self._write_atomic_tsv(cohort_muts_tsv, MUTATION_TSV_HEADERS, all_mutations_rows)
        cohort_deliverables["cohort_mutation_catalog"] = cohort_muts_tsv

        # 6. Emit Cohort Horizontal Mutations Summary (sample_id | hg | count | list)
        cohort_mut_summary_tsv = self.dir_muts / "cohort_mutations_summary.tsv"
        mut_summary_headers = ["sample_id", "primary_haplogroup", "mutations_count", "mutations"]
        mut_summary_rows: List[List[str]] = []

        for sid in sample_ids:
            rec_dict = master_json_records.get(sid, {})
            muts = rec_dict.get("mutations", [])
            hg = rec_dict.get("haplogroup", {}).get("primary_haplogroup", "UNKNOWN")
            mut_string = " ".join(m.get("mutation_notation", "") for m in muts if m.get("mutation_notation"))
            mut_summary_rows.append([sid, hg, str(len(muts)), mut_string])

        self._write_atomic_tsv(cohort_mut_summary_tsv, mut_summary_headers, mut_summary_rows)
        cohort_deliverables["cohort_mutations_summary"] = cohort_mut_summary_tsv

        logger.info(
            f"✅ Cohort aggregation complete: {len(cohort_deliverables)} deliverables emitted "
            f"for {len(sample_ids)} samples in '{self.results_dir}'."
        )

        return cohort_deliverables

    # --------------------------------------------------------------------------
    # Hardware fsync Atomic Writing Primitives
    # --------------------------------------------------------------------------
    @staticmethod
    def _write_atomic_tsv(target_path: Path, headers: List[str], rows: List[List[str]]) -> None:
        """Atomically write TSV payload with row length validation and fsync flush."""
        target_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = target_path.with_name(f".tmp_{target_path.name}_{os.getpid()}")
        try:
            with open(tmp_path, "w", encoding="utf-8", newline="") as fp:
                writer = csv.writer(fp, delimiter="\t", lineterminator="\n")
                writer.writerow(headers)
                for r_idx, row in enumerate(rows, start=1):
                    if len(row) != len(headers):
                        raise SynthesisOrchestrationError(
                            f"Row {r_idx} in '{target_path.name}' has {len(row)} values; expected {len(headers)}."
                        )
                    writer.writerow(row)
                fp.flush()
                os.fsync(fp.fileno())
            tmp_path.replace(target_path)
        except Exception as err:
            if tmp_path.is_file():
                tmp_path.unlink(missing_ok=True)
            raise SynthesisOrchestrationError(f"Failed atomic TSV write to '{target_path}': {err}") from err

    @staticmethod
    def _write_atomic_text(target_path: Path, content: str) -> None:
        """Atomically write text payload with fsync flush."""
        target_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = target_path.with_name(f".tmp_{target_path.name}_{os.getpid()}")
        try:
            with open(tmp_path, "w", encoding="utf-8") as fp:
                fp.write(content)
                fp.flush()
                os.fsync(fp.fileno())
            tmp_path.replace(target_path)
        except Exception as err:
            if tmp_path.is_file():
                tmp_path.unlink(missing_ok=True)
            raise SynthesisOrchestrationError(f"Failed atomic text write to '{target_path}': {err}") from err