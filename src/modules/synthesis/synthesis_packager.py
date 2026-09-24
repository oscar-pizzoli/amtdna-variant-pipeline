"""
Block 8 Synthesis Packager: Procedure-oriented transactional deliverable serializer.
Dispatches audited artifacts across structured functional directories partitioned by sample_id:
results/
├── master_summary/
├── fasta/<sample_id>/
├── vcf/<sample_id>/
├── haplogroups/<sample_id>/
├── contamination/<sample_id>/
│   ├── schmutzi/
│   └── dissent/
├── mutations/<sample_id>/
├── qc_sequencing/<sample_id>/
└── reports_markdown/<sample_id>/
"""

import csv
from datetime import datetime, timezone
import json
import logging
import math
import os
from pathlib import Path
import shutil
from typing import Any, Dict, List, Optional, Tuple

from src.modules.synthesis.models import (
    MUTATION_TSV_HEADERS,
    SAMPLE_MASTER_TSV_HEADERS,
    CategorizedWarning,
    PackagingConfig,
    ProvenanceItem,
    ProvenanceStatus,
    SampleMasterRecord,
    WarningCategory,
    WarningSeverity,
    format_na,
    format_pct_100,
    format_pct_fraction,
)

logger = logging.getLogger("mtDNA_Pipeline.SynthesisPackager")


class PackagingExecutionError(Exception):
    """Raised when an unrecoverable I/O failure or contract violation occurs during packaging."""
    pass


class SynthesisPackager:
    """
    Exports single-sample and cohort deliverables into structured procedure directories
    partitioned by sample_id, using atomic transactional writes.
    """

    PMD_ANNOTATED_HEADERS = [
        "read_index",
        "p_damage",
        "p_null",
        "likelihood_ratio",
        "pmd_score",
        "is_damage_significant",
    ]

    MUTATIONS_LIST_HEADERS = [
        "sample_id",
        "mutations_count",
        "mutations",
    ]

    def __init__(
        self,
        results_dir: Path,
        operational_files_dir: Optional[Path] = None,
        config: Optional[PackagingConfig] = None,
    ) -> None:
        self.results_dir = Path(results_dir)
        self.operational_files_dir = Path(operational_files_dir) if operational_files_dir else self.results_dir
        self.config = config or PackagingConfig(
            reports_dir=self.results_dir / "master_summary",
            emit_tsv_summary=True,
            emit_mutations_catalog=True,
            emit_json_master=True,
            emit_markdown_sheets=True,
            enforce_strict_provenance=True,
        )

        # Initialize base procedural directories
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

    # --------------------------------------------------------------------------
    # Transactional I/O Helpers
    # --------------------------------------------------------------------------
    @staticmethod
    def _atomic_write_text(target_path: Path, content: str) -> None:
        target_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = target_path.with_name(f".tmp_{target_path.name}_{os.getpid()}")
        try:
            with open(tmp_path, "w", encoding="utf-8") as fp:
                fp.write(content)
            tmp_path.replace(target_path)
        except Exception as err:
            if tmp_path.is_file():
                tmp_path.unlink(missing_ok=True)
            raise PackagingExecutionError(f"Failed atomic write to '{target_path}': {err}") from err

    @staticmethod
    def _atomic_write_tsv(target_path: Path, headers: List[str], rows: List[List[str]]) -> None:
        target_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = target_path.with_name(f".tmp_{target_path.name}_{os.getpid()}")
        try:
            with open(tmp_path, "w", encoding="utf-8", newline="") as fp:
                writer = csv.writer(fp, delimiter="\t", lineterminator="\n")
                writer.writerow(headers)
                for r_idx, row in enumerate(rows, start=1):
                    if len(row) != len(headers):
                        raise PackagingExecutionError(
                            f"Row {r_idx} in '{target_path.name}' has {len(row)} values; expected {len(headers)}."
                        )
                    writer.writerow(row)
            tmp_path.replace(target_path)
        except Exception as err:
            if tmp_path.is_file():
                tmp_path.unlink(missing_ok=True)
            raise PackagingExecutionError(f"Failed atomic TSV export for '{target_path}': {err}") from err

    @staticmethod
    def _atomic_copy_file(src_path: Path, dst_path: Path) -> None:
        dst_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = dst_path.with_name(f".tmp_{dst_path.name}_{os.getpid()}")
        try:
            shutil.copy2(src_path, tmp_path)
            tmp_path.replace(dst_path)
        except Exception as err:
            if tmp_path.is_file():
                tmp_path.unlink(missing_ok=True)
            raise PackagingExecutionError(f"Failed atomic copy from '{src_path}' to '{dst_path}': {err}") from err

    # --------------------------------------------------------------------------
    # 1. Pipeline Artifact Copy & Relocation
    # --------------------------------------------------------------------------
    def relocate_sample_artifacts(self, record: SampleMasterRecord) -> Dict[str, Path]:
        sample_id = record.sample_id
        relocated: Dict[str, Path] = {}

        prov_status_by_path: Dict[str, ProvenanceStatus] = {
            str(Path(p.file_path).resolve()): p.status
            for p in record.provenance
            if p.file_path not in {"UNSPECIFIED", "NOT_PRODUCED"}
        }

        b5_dir = self.operational_files_dir / "contamination" / sample_id
        b6_dir = self.operational_files_dir / "consensus" / sample_id
        b7_dir = self.operational_files_dir / "haplogroup" / sample_id

        # Target directories partitioned by sample_id and analysis method
        sample_fasta_dir = self.dir_fasta / sample_id
        sample_vcf_dir = self.dir_vcf / sample_id
        sample_haplo_dir = self.dir_haplo / sample_id
        sample_contam_dir = self.dir_contam / sample_id
        sample_schmutzi_dir = sample_contam_dir / "schmutzi"
        sample_dissent_dir = sample_contam_dir / "dissent"

        for d in (
            sample_fasta_dir,
            sample_vcf_dir,
            sample_haplo_dir,
            sample_contam_dir,
            sample_schmutzi_dir,
            sample_dissent_dir,
        ):
            d.mkdir(parents=True, exist_ok=True)

        relocation_plan: List[Tuple[str, List[Path], Path, bool]] = [
            (
                "canonical_fasta",
                [b6_dir / f"{sample_id}.consensus.fasta", b6_dir / f"{sample_id}.canonical.fasta"],
                sample_fasta_dir / f"{sample_id}.canonical.fasta",
                True,
            ),
            (
                "heteroplasmy_fasta",
                [b6_dir / f"{sample_id}.consensus_iupac.fasta", b6_dir / f"{sample_id}.consensus.iupac.fasta"],
                sample_fasta_dir / f"{sample_id}.heteroplasmy.fasta",
                False,
            ),
            (
                "indel_fasta",
                [b6_dir / f"{sample_id}.consensus_indel.fasta", b6_dir / f"{sample_id}.consensus.indel.fasta"],
                sample_fasta_dir / f"{sample_id}.indel_aware.fasta",
                False,
            ),
            (
                "variants_vcf",
                [b6_dir / f"{sample_id}.variants.vcf.gz"],
                sample_vcf_dir / f"{sample_id}.variants.vcf.gz",
                True,
            ),
            (
                "variants_vcf_tbi",
                [b6_dir / f"{sample_id}.variants.vcf.gz.tbi", b6_dir / f"{sample_id}.variants.vcf.tbi"],
                sample_vcf_dir / f"{sample_id}.variants.vcf.gz.tbi",
                False,
            ),
            (
                "subconsensus_vcf",
                [b6_dir / f"{sample_id}.subconsensus_heteroplasmies.vcf.gz", b6_dir / f"{sample_id}.subconsensus.vcf.gz"],
                sample_vcf_dir / f"{sample_id}.subconsensus.vcf.gz",
                False,
            ),
            (
                "haplogroup_report",
                [
                    b7_dir / f"{sample_id}.haplogrep.fasta.raw.tsv",
                    b7_dir / f"{sample_id}.haplogrep.vcf.raw.tsv",
                    b7_dir / f"{sample_id}_fasta.txt",
                    b7_dir / f"{sample_id}.haplogrep.txt",
                ],
                sample_haplo_dir / f"{sample_id}.haplogroup_report.tsv",
                False,
            ),
            # --- Dedicated Schmutzi Deliverables ---
            (
                "contamination_deam_pdf",
                [b5_dir / f"{sample_id}.cont.pdf", b5_dir / f"{sample_id}_endo.cont.pdf"],
                sample_schmutzi_dir / f"{sample_id}.cont.pdf",
                False,
            ),
            (
                "contamination_positions_pdf",
                [b5_dir / f"{sample_id}_final.cont.pdf", b5_dir / f"{sample_id}.final.cont.pdf"],
                sample_schmutzi_dir / f"{sample_id}_final.cont.pdf",
                False,
            ),
            (
                "contamination_estimates",
                [b5_dir / f"{sample_id}_final.cont.est", b5_dir / f"{sample_id}.cont.est"],
                sample_schmutzi_dir / f"{sample_id}_final.cont.est",
                False,
            ),
            (
                "contamination_likelihood_grid",
                [b5_dir / f"{sample_id}_final.cont", b5_dir / f"{sample_id}.cont"],
                sample_schmutzi_dir / f"{sample_id}_final.cont",
                False,
            ),
            # --- Dedicated Dissent Deliverables ---
            (
                "dissent_tsv",
                [
                    b5_dir / sample_id / f"{sample_id}.dissent.tsv",
                    b5_dir / f"{sample_id}.dissent.tsv",
                    b5_dir / "dissent" / f"{sample_id}.dissent.tsv",
                    b5_dir / sample_id / sample_id / f"{sample_id}.dissent.tsv",
                ],
                sample_dissent_dir / f"{sample_id}.dissent.tsv",
                False,
            ),
            (
                "fast_screen_json",
                [
                    b5_dir / sample_id / f"{sample_id}.fast_screen.json",
                    b5_dir / f"{sample_id}.fast_screen.json",
                    b5_dir / "dissent" / f"{sample_id}.fast_screen.json",
                    b5_dir / sample_id / sample_id / f"{sample_id}.fast_screen.json",
                ],
                sample_dissent_dir / f"{sample_id}.fast_screen.json",
                False,
            ),
        ]

        is_gated = (
            record.contamination.recommended_action in {"REJECT_CONTAMINATED", "EXCLUDE_OR_MANUAL_REVIEW"}
            or record.authenticity.authentication_status in {"INCONCLUSIVE", "NOISY_OR_INSUFFICIENT"}
        )

        for label, candidates, dst_path, is_crit in relocation_plan:
            src_file = next((p for p in candidates if p.is_file()), None)

            if src_file:
                resolved_key = str(src_file.resolve())
                prov_st = prov_status_by_path.get(resolved_key)

                if self.config.enforce_strict_provenance:
                    if is_crit and prov_st != ProvenanceStatus.VERIFIED:
                        st_name = prov_st.value if prov_st else "UNINSPECTED"
                        raise PackagingExecutionError(
                            f"Refusing to relocate critical artifact '{label}' from '{src_file}': "
                            f"provenance status is '{st_name}' (strict provenance requires VERIFIED)."
                        )
                    elif not is_crit and prov_st in {
                        ProvenanceStatus.MISMATCH,
                        ProvenanceStatus.SEMANTIC_ERROR,
                        ProvenanceStatus.IO_ERROR,
                    }:
                        raise PackagingExecutionError(
                            f"Refusing to relocate non-critical artifact '{label}' from '{src_file}': "
                            f"provenance check failed with status '{prov_st.value}'."
                        )
                else:
                    if prov_st == ProvenanceStatus.MISMATCH:
                        raise PackagingExecutionError(
                            f"Refusing to relocate corrupted artifact '{label}' from '{src_file}': "
                            f"ProvenanceStatus is MISMATCH."
                        )

                self._atomic_copy_file(src_file, dst_path)
                relocated[dst_path.name] = dst_path
            else:
                if is_crit and not is_gated:
                    logger.warning(
                        f"[{sample_id}] Critical artifact '{label}' not found among candidates: "
                        f"{[p.name for p in candidates]}"
                    )

        return relocated

    # --------------------------------------------------------------------------
    # 2. Sequencing QC & Trimming Detail Export
    # --------------------------------------------------------------------------
    def export_sequencing_qc(self, record: SampleMasterRecord) -> Dict[str, Path]:
        artifacts: Dict[str, Path] = {}
        sid = record.sample_id
        sample_qc_dir = self.dir_qc / sid
        sample_qc_dir.mkdir(parents=True, exist_ok=True)

        seq_qc = getattr(record, "sequencing_qc", None)

        if seq_qc and seq_qc.per_base_quality_cycles:
            q_path = sample_qc_dir / f"{sid}.per_base_quality.tsv"
            headers = ["cycle", "mean_quality", "q20_bases_count", "q30_bases_count", "mean_a", "mean_c", "mean_g", "mean_t"]
            rows: List[List[str]] = []
            for item in seq_qc.per_base_quality_cycles:
                rows.append([
                    str(item.get("cycle", "")),
                    format_na(item.get("mean_quality"), "{:.2f}"),
                    str(item.get("q20_bases_count", 0)),
                    str(item.get("q30_bases_count", 0)),
                    format_na(item.get("mean_a"), "{:.4f}"),
                    format_na(item.get("mean_c"), "{:.4f}"),
                    format_na(item.get("mean_g"), "{:.4f}"),
                    format_na(item.get("mean_t"), "{:.4f}"),
                ])
            self._atomic_write_tsv(q_path, headers, rows)
            artifacts["per_base_quality"] = q_path

        if seq_qc and seq_qc.read_length_distribution:
            l_path = sample_qc_dir / f"{sid}.read_length_dist.tsv"
            headers = ["length_bp", "reads_count"]
            rows = [[str(k), str(v)] for k, v in sorted(seq_qc.read_length_distribution.items())]
            self._atomic_write_tsv(l_path, headers, rows)
            artifacts["read_length_dist"] = l_path

        return artifacts

    # --------------------------------------------------------------------------
    # 3. Authenticity & PMD Scores Annotated Export
    # --------------------------------------------------------------------------
    def export_annotated_pmd_scores(self, record: SampleMasterRecord) -> Optional[Path]:
        sid = record.sample_id
        sample_damage_dir = self.operational_files_dir / "damage" / sid
        if not sample_damage_dir.is_dir():
            sample_damage_dir.mkdir(parents=True, exist_ok=True)

        raw_candidates = [
            sample_damage_dir / f"{sid}.pmd_scores.tsv",
            sample_damage_dir / f"{sid}.pmdscores.tsv",
            sample_damage_dir / "pmdscores.tsv",
            sample_damage_dir / "tables" / f"{sid}.pmd_scores.tsv",
        ]
        raw_pmd_path = next((p for p in raw_candidates if p.is_file()), None)
        if not raw_pmd_path:
            return None

        out_path = sample_damage_dir / f"{sid}.pmd_scores.annotated.tsv"
        rows: List[List[str]] = []

        try:
            with open(raw_pmd_path, "r", encoding="utf-8") as fp:
                for idx, line in enumerate(fp, start=1):
                    line_clean = line.strip()
                    if not line_clean:
                        continue
                    parts = line_clean.split()
                    if len(parts) < 4:
                        continue

                    p_dmg = parts[0]
                    p_null = parts[1]
                    lr = parts[2]
                    pmd_sc = parts[3]

                    try:
                        score_float = float(pmd_sc)
                        is_sig = "TRUE" if score_float >= 3.0 else "FALSE"
                    except ValueError:
                        is_sig = "NA"

                    rows.append([str(idx), p_dmg, p_null, lr, pmd_sc, is_sig])

            self._atomic_write_tsv(out_path, self.PMD_ANNOTATED_HEADERS, rows)
            logger.info(f"[{sid}] Exported annotated PMD scores -> {out_path.name}")
            return out_path
        except Exception as err:
            logger.warning(f"[{sid}] Could not generate annotated PMD table: {err}")
            return None

    # --------------------------------------------------------------------------
    # 4. Mutations Catalogs (Vertical Tidy & Horizontal List)
    # --------------------------------------------------------------------------
    def export_mutation_catalog_tsv(self, record: SampleMasterRecord) -> Path:
        sample_muts_dir = self.dir_muts / record.sample_id
        sample_muts_dir.mkdir(parents=True, exist_ok=True)
        out_path = sample_muts_dir / f"{record.sample_id}.mutation_catalog.tsv"
        rows: List[List[str]] = [m.to_tsv_row() for m in record.mutations]
        self._atomic_write_tsv(out_path, MUTATION_TSV_HEADERS, rows)
        return out_path

    def export_mutation_list_tsv(self, record: SampleMasterRecord) -> Path:
        sample_muts_dir = self.dir_muts / record.sample_id
        sample_muts_dir.mkdir(parents=True, exist_ok=True)
        out_path = sample_muts_dir / f"{record.sample_id}.mutations_list.tsv"
        mut_str = " ".join(m.mutation_notation for m in record.mutations)
        row = [record.sample_id, str(len(record.mutations)), mut_str]
        self._atomic_write_tsv(out_path, self.MUTATIONS_LIST_HEADERS, [row])
        return out_path

    def export_master_json(self, record: SampleMasterRecord) -> Path:
        out_path = self.dir_master / f"{record.sample_id}.master.json"
        content = json.dumps(record.to_dict(), indent=2, ensure_ascii=False)
        self._atomic_write_text(out_path, content + "\n")
        return out_path

    def export_sample_summary_tsv(self, record: SampleMasterRecord) -> Path:
        out_path = self.dir_master / f"{record.sample_id}.summary.tsv"
        row = record.to_summary_tsv_row()
        self._atomic_write_tsv(out_path, SAMPLE_MASTER_TSV_HEADERS, [row])
        return out_path

    # --------------------------------------------------------------------------
    # 5. Markdown Forensic Report Render
    # --------------------------------------------------------------------------
    def render_markdown_report(self, record: SampleMasterRecord) -> Path:
        """Render forensic report with detailed missing/confirmed diagnostic mutations."""
        sid = record.sample_id
        sample_reports_dir = self.dir_reports / sid
        sample_reports_dir.mkdir(parents=True, exist_ok=True)
        out_path = sample_reports_dir / f"{sid}.report.md"

        now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

        status_badges = {
            "PASS": "🟢 **PASS**",
            "WARN": "🟡 **WARN (Cautionary Notices)**",
            "FAIL": "🔴 **FAIL (Validation Rejected)**",
        }
        overall_badge = status_badges.get(record.overall_status, f"⚪ **{record.overall_status}**")

        mc = record.macro_complexes
        sq = record.sequencing_qc

        depth_val = record.consensus.mean_consensus_depth
        depth_eval = "⚪ N/A (Missing)" if depth_val is None else ("✅ Optimal" if depth_val >= 10.0 else "⚠️ Low")

        breadth_val = record.alignment.coverage_breadth_1x
        breadth_eval = "⚪ N/A (Missing)" if breadth_val is None else ("✅ High" if breadth_val >= 0.95 else "⚠️ Incomplete")

        n_val = record.consensus.n_fraction_canonical
        n_eval = "⚪ N/A (Missing)" if n_val is None else ("✅ Clean" if n_val <= 0.15 else "⚠️ High N-Content")

        contam_val = record.contamination.contamination_fraction
        contam_eval = "⚪ N/A (Missing)" if contam_val is None else ("✅ Safe" if contam_val <= 0.02 else "⚠️ Contaminated")

        # Formatting confirmed vs missing diagnostic mutations
        missing_muts = getattr(record.haplogroup, "missing_diagnostic_mutations", [])
        found_muts = getattr(record.haplogroup, "found_diagnostic_mutations", [])

        missing_display = (
            f"⚠️ **`{len(missing_muts)} missing`**: `{', '.join(missing_muts)}`"
            if missing_muts
            else "✅ *None (All expected diagnostic mutations confirmed)*"
        )
        found_display = (
            f"`{', '.join(found_muts)}` ({len(found_muts)})"
            if found_muts
            else "⚪ *None detected*"
        )

        lines: List[str] = [
            f"# Forensic Mitochondrial Synthesis Report: `{sid}`",
            f"*Generated: {now_utc} | Pipeline Synthesis Engine v1.2.0*",
            "",
            "---",
            "",
            "## 1. Executive Summary & Verdict",
            "",
            f"* **Final Pipeline Verdict**: {overall_badge}",
            f"* **Biological Authentication**: `{record.authenticity.biological_type.value}`",
            f"* **Resolved Haplogroup**: **`{record.haplogroup.primary_haplogroup}`** "
            f"(Score: `{format_na(record.haplogroup.quality_score, '{:.4f}')}`, Margin: `{format_na(record.haplogroup.quality_margin, '{:.4f}')}`)",
            f"* **Diagnostic Mutations Confirmed**: {found_display}",
            f"* **Missing Diagnostic Mutations**: {missing_display}",
            f"* **Diagnostic Narrative**: {record.get_warning_narrative()}",
            "",
            "### Core Quality & Authenticity Dashboard",
            "",
            "| Metric | Value | Reference / Standard | Evaluation |",
            "| :--- | :--- | :--- | :--- |",
            f"| **Mean Consensus Depth** | `{format_na(depth_val, '{:.1f}')}x` | $\\ge 10.0x$ | {depth_eval} |",
            f"| **Genome Breadth (1x)** | `{format_pct_fraction(breadth_val)}` | $\\ge 95.0\\%$ | {breadth_eval} |",
            f"| **N-Fraction Canonical** | `{format_pct_fraction(n_val)}` | $\\le 15.0\\%$ (SequenceValidator limit) | {n_eval} |",
            f"| **Contamination Estimate** | `{format_pct_fraction(contam_val)}` | $\\le 2.0\\%$ | {contam_eval} |",
            f"| **5'-C>T Damage (Pos 1)** | `{format_pct_fraction(record.authenticity.c_to_t_5p_rate_pos1)}` | aDNA degradation marker | Metric: `{record.authenticity.authentication_status}` |",
            f"| **Validated Mutations** | `{len(record.mutations)}` | Diagnostic & Private catalog | `{record.consensus.heteroplasmy_sites_pass}` pass heteroplasmies |",
            "",
            "---",
            "",
            "## 2. Sequencing Dynamics, Quality & Trimming (Block 2)",
            "",
            f"* **Raw Reads Ingested**: `{format_na(sq.raw_reads_total if sq.raw_reads_total is not None else record.preprocessing.raw_reads_total)}`",
            f"* **Clean Retained Reads**: `{format_na(sq.clean_reads_total if sq.clean_reads_total is not None else record.preprocessing.retained_reads_total)}` ({format_pct_fraction(record.preprocessing.retention_fraction)})",
            f"* **Merged / Collapsed Pairs**: `{format_pct_fraction(record.preprocessing.collapsed_fraction)}`",
            f"* **Adapter Trimmed Reads**: `{format_pct_100(sq.adapter_trimmed_reads_pct)}`",
            f"* **Post-Trimming $\\ge Q30$ Rate**: `{format_pct_fraction(sq.q30_rate_post)}`",
            f"* **GC Content Fraction**: `{format_pct_100(sq.gc_content_pct)}`",
            f"* **Mean Length (Pre $\\to$ Post)**: `{format_na(sq.mean_read_length_pre, '{:.1f}')} bp \\to {format_na(sq.mean_read_length_post, '{:.1f}')} bp`",
            "",
            "---",
            "",
            "## 3. Mitochondrial Functional Macro-Complexes Coverage (Block 6)",
            "",
            "| Functional Macro-Complex | rCRS Coordinates | Length | Mean Depth | Breadth $\\ge 1\\text{x}$ | Evaluation |",
            "| :--- | :--- | :--- | :--- | :--- | :--- |",
        ]

        def _eval_complex(cov1x: Optional[float]) -> str:
            if cov1x is None:
                return "⚪ N/A"
            return "✅ Sufficient" if cov1x >= 0.90 else "⚠️ Drop"

        lines.extend([
            f"| **D-Loop (Control Region)** | 16024-576 | 1,122 bp | `{format_na(mc.cov_dloop_mean, '{:.1f}')}x` | `{format_pct_fraction(mc.cov_dloop_1x)}` | {_eval_complex(mc.cov_dloop_1x)} |",
            f"| **rRNA (12S / 16S)** | 648-3229 | 2,582 bp | `{format_na(mc.cov_rrna_mean, '{:.1f}')}x` | `{format_pct_fraction(mc.cov_rrna_1x)}` | {_eval_complex(mc.cov_rrna_1x)} |",
            f"| **Complex IV (COX1-3)** | Intercalated | 3,038 bp | `{format_na(mc.cov_cox_mean, '{:.1f}')}x` | `{format_pct_fraction(mc.cov_cox_1x)}` | {_eval_complex(mc.cov_cox_1x)} |",
            f"| **Complex I (ND1-6, ND4L)** | Intercalated | 6,361 bp | `{format_na(mc.cov_nd_mean, '{:.1f}')}x` | `{format_pct_fraction(mc.cov_nd_1x)}` | {_eval_complex(mc.cov_nd_1x)} |",
            f"| **Complex III (CYB)** | 14747-15887 | 1,141 bp | `{format_na(mc.cov_cyb_mean, '{:.1f}')}x` | `{format_pct_fraction(mc.cov_cyb_1x)}` | {_eval_complex(mc.cov_cyb_1x)} |",
            "",
            "---",
            "",
            "## 4. Authenticity, Degradation & Contamination (Blocks 4 & 5)",
            "",
            f"* **Enzymatic Library Treatment**: `{record.authenticity.library_treatment}`",
            f"* **Median Insert Length**: `{format_na(record.authenticity.median_fragment_length, '{:.1f}')} bp`",
            f"* **Fraction PMD Score $\\ge 3$**: `{format_na(record.authenticity.pmd_fraction_ge_3, '{:.4f}')}`",
            f"* **MapDamage Quality Rescaling**: `{'YES (rescaled BAM deployed)' if record.authenticity.rescaling_performed else 'NO'}`",
            f"* **Contamination Point Estimate**: `{format_pct_fraction(record.contamination.contamination_fraction)}`",
            f"* **Contamination 95% CI**: `[{format_na(record.contamination.ci_95_lower, '{:.4f}')} - {format_na(record.contamination.ci_95_upper, '{:.4f}')}]`",
            f"* **Gating Decision Applied**: `{record.contamination.recommended_action}`",
            "",
            "---",
            "",
            "## 5. Diagnostic Warnings & Anomalies",
            "",
        ])

        if not record.warnings:
            lines.append("> 🟢 **No diagnostic warnings or technical anomalies reported.**")
        else:
            lines.append("| Category | Severity | Code | Origin | Message |")
            lines.append("| :--- | :--- | :--- | :--- | :--- |")
            for w in record.warnings:
                sev_symbol = "🔴" if w.severity == WarningSeverity.CRITICAL else ("🟡" if w.severity == WarningSeverity.WARNING else "ℹ️")
                lines.append(f"| `{w.category.value}` | {sev_symbol} `{w.severity.value}` | `{w.code}` | `{w.source_block}` | {w.message} |")

        lines.extend([
            "",
            "---",
            "",
            "## 6. Reconciled Mutation Catalog Highlights",
            "",
            f"* **Total Catalogued Mutations**: `{len(record.mutations)}`",
            "",
            "| Position | Notation | Type | Source | Role | VAF | Depth | Spatial Damage? | Filter |",
            "| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |",
        ])

        for m in record.mutations[:40]:
            vaf_s = f"{m.vaf:.2f}" if m.vaf is not None else "NA"
            dp_s = str(m.depth) if m.depth is not None else "NA"
            dmg_s = "True" if m.spatial_damage_suspect is True else ("False" if m.spatial_damage_suspect is False else "NA")
            lines.append(f"| `{m.position}` | **`{m.mutation_notation}`** | `{m.variant_type.value}` | `{m.stream_source}` | `{m.role.value}` | `{vaf_s}` | `{dp_s}` | `{dmg_s}` | `{m.filter_status}` |")

        if len(record.mutations) > 40:
            lines.append(f"| ... | *({len(record.mutations) - 40} additional mutations truncated for brevity)* | ... | ... | ... | ... | ... | ... | ... |")

        lines.extend([
            "",
            "---",
            "",
            "## 7. Chain-of-Custody & Cryptographic Provenance",
            "",
            f"* **Composite Upstream Hash**: `{record.upstream_composite_hash or 'STANDALONE_RUN'}`",
            "",
            "| Artifact | Status | Expected SHA-256 | Actual SHA-256 | Size |",
            "| :--- | :--- | :--- | :--- | :--- |",
        ])

        for p in record.provenance:
            icon = "✅" if p.status == ProvenanceStatus.VERIFIED else ("❌" if p.status == ProvenanceStatus.MISMATCH else "⚠️")
            exp_s = f"`{p.expected_sha256[:12]}...`" if len(p.expected_sha256) > 12 else f"`{p.expected_sha256}`"
            act_s = f"`{p.actual_sha256[:12]}...`" if p.actual_sha256 and len(p.actual_sha256) > 12 else (f"`{p.actual_sha256}`" if p.actual_sha256 else "*None*")
            sz_s = f"{p.file_size_bytes:,} B" if p.file_size_bytes is not None else "NA"
            lines.append(f"| `{p.artifact_name}` | {icon} `{p.status.value}` | {exp_s} | {act_s} | {sz_s} |")

        lines.append("")
        content = "\n".join(lines)
        self._atomic_write_text(out_path, content)
        return out_path

    # --------------------------------------------------------------------------
    # 6. Technical Explanatory Documentation
    # --------------------------------------------------------------------------
    def export_contamination_readme(self) -> Path:
        """Generates an exhaustive technical README in results/contamination/ explaining deliverable artifacts."""
        readme_path = self.dir_contam / "README.md"
        content = """# Forensic Mitochondrial Contamination Assessment Guide (Schmutzi & Dissent)

This directory contains human-interpretable deliverable artifacts produced by the Bayesian contamination estimation framework (Block 5).

---

## Directory Organization

For each evaluated specimen, deliverables are partitioned into method-specific subdirectories:
* **`results/contamination/<sample_id>/schmutzi/`**: Contains Bayesian joint mixture model outputs and diagnostic posterior probability density curves.
* **`results/contamination/<sample_id>/dissent/`**: Contains fast screening triage metrics and high-sensitivity dissenting read fraction tables.

---

## 1. Diagnostic PDF Visualizations (Schmutzi)

* **`schmutzi/{sample_id}.cont.pdf`**
  * **Model:** *Posterior probability for contamination using deamination patterns*
  * **Biological Principle:** Authenticates ancient DNA by contrasting molecules exhibiting post-mortem degradation kinetics (terminal 5'-C>T and 3'-G>A transitions) against undamaged modern molecules.
  * **How to Interpret:**
    * A sharp, narrow posterior distribution centered near 0.00 (0.0% – 0.5%) confirms that undamaged molecules do not represent an exogenous modern contaminant pool.
    * A shifted peak away from zero indicates that undamaged molecules deviate from the expected endogenous background, confirming modern library contamination.

* **`schmutzi/{sample_id}_final.cont.pdf`**
  * **Model:** *Posterior probability for contamination using mitochondrial positions*
  * **Genomic Principle:** Jointly models variant allele frequencies across the mitochondrial genome against an empirical reference database of 197 putative modern human contaminant mitogenomes (`schmutzi_db/197`).
  * **How to Interpret:**
    * **Peak of the Curve:** Corresponds to the maximum a posteriori (MAP) point estimate of contamination.
    * **Curve Envelope:** Represents the 95% Bayesian Credible Interval (CI). A narrow distribution centered below 0.02 (2.0%) confirms high statistical confidence and minimal exogenous contamination.

---

## 2. Numerical Estimates & Dissent Artifacts

* **`schmutzi/{sample_id}_final.cont.est`**
  * Contains three primary numerical values:
    `[95%_CI_lower    point_estimate    95%_CI_upper]`
  * *Example:* `0.00  0.01  0.02` corresponds to a 1.0% contamination point estimate bounded within a [0.0%, 2.0%] credible envelope.

* **`schmutzi/{sample_id}_final.cont`**
  * Exhaustive posterior log-likelihood matrix computed across the discrete contamination grid ($c \\in [0.00, 0.99]$). The maximum log-likelihood value corresponds to the reported point estimate.

* **`dissent/{sample_id}.dissent.tsv`**
  * Position-by-position breakdown of non-consensus (dissenting) bases observed across phylogenetically informative and diagnostic mitochondrial loci.

* **`dissent/{sample_id}.fast_screen.json`**
  * High-concurrency screening metrics evaluating coverage depth and preliminary dissenting fraction to assess bypass eligibility.

---

## 3. Decision Matrix & Action Thresholds

* **LOW CONTAMINATION / HIGH CONFIDENCE (≤ 2.0%):**
  * The specimen is deemed clean and authentic. Standard consensus base calling and heteroplasmy auditing proceed unconstrained (`STANDARD_CALLING`).
* **MODERATE CONTAMINATION (2.0% – 5.0%):**
  * Mild exogenous contamination detected. The pipeline automatically enforces `RESTRICT_HETEROPLASMY`: candidate heteroplasmy thresholds are elevated to Variant Allele Fraction (VAF) ≥ 25% and minimum depth ≥ 20x to suppress contaminant allele misidentification.
* **ELEVATED CONTAMINATION (> 5.0% – 10.0%):**
  * Significant exogenous contamination. Requires read-level PMDtools damage filtering (`PMD_FILTERING_REQUIRED`) or manual forensic review.
* **CRITICAL CONTAMINATION (> 10.0%):**
  * Exogenous sequence burden compromises data integrity. The specimen is gated and excluded from final release (`REJECT_CONTAMINATED`).
* **UNRELIABLE / INSUFFICIENT DATA:**
  * Assigned when statistical power is insufficient (95% CI width > 8.0%, mean coverage depth < 5.0x, or < 10 diagnostic informative positions). Gating defaults to restricted heteroplasmy calling.
"""
        self._atomic_write_text(readme_path, content.strip() + "\n")
        return readme_path

    # --------------------------------------------------------------------------
    # Master Execution Method
    # --------------------------------------------------------------------------
    def package_sample(self, record: SampleMasterRecord) -> Dict[str, Path]:
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

        emitted: Dict[str, Path] = {}

        # 1. Relocate artifacts to sample-partitioned results directories
        relocated = self.relocate_sample_artifacts(record)
        emitted.update(relocated)

        # 2. Export sequencing and trimming QC metrics
        seq_qc = self.export_sequencing_qc(record)
        emitted.update(seq_qc)

        # 3. Export annotated PMD scores inside operational_files/damage/<sample_id>/
        pmd_annotated = self.export_annotated_pmd_scores(record)
        if pmd_annotated:
            emitted["pmd_scores_annotated_tsv"] = pmd_annotated

        # 4. Export tidy mutation catalog and horizontal space-separated list
        if self.config.emit_mutations_catalog:
            emitted["mutation_catalog_tsv"] = self.export_mutation_catalog_tsv(record)
            emitted["mutations_list_tsv"] = self.export_mutation_list_tsv(record)

        # 5. Export structured master JSON
        if self.config.emit_json_master:
            emitted["master_json"] = self.export_master_json(record)

        # 6. Export positional summary TSV row
        if self.config.emit_tsv_summary:
            emitted["summary_tsv"] = self.export_sample_summary_tsv(record)

        # 7. Render human-readable forensic Markdown report
        if self.config.emit_markdown_sheets:
            emitted["report_markdown"] = self.render_markdown_report(record)

        # 8. Export technical explanatory README for contamination deliverables
        readme_path = self.export_contamination_readme()
        emitted["contamination_readme"] = readme_path

        logger.info(
            f"[{record.sample_id}] Packaging completed: {len(emitted)} artifacts emitted "
            f"across sample-partitioned procedural directories (status={record.overall_status})."
        )
        return emitted