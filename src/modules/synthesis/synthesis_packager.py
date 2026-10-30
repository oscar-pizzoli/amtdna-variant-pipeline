"""
Block 8 Synthesis Packager: Procedure-oriented transactional deliverable serializer.
Dispatches audited artifacts across structured functional directories partitioned by sample_id:
results/
├── master_summary/
│   ├── cohort_qc_summary.tsv
│   ├── cohort_master.json
│   ├── WARNINGS_GLOSSARY.md
│   └── <sample_id>.summary.tsv
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
    SAMPLE_MASTER_BAM_TSV_HEADERS,
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
            # Dedicated Schmutzi Deliverables
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
            # Dedicated Dissent Deliverables
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
        """Export detailed quality metrics tables; bypass completely in BAM mode."""
        if record.is_bam_mode:
            return {}

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
        headers = SAMPLE_MASTER_BAM_TSV_HEADERS if record.is_bam_mode else SAMPLE_MASTER_TSV_HEADERS
        self._atomic_write_tsv(out_path, headers, [row])
        return out_path

    # --------------------------------------------------------------------------
    # 5. Markdown Forensic Report Render
    # --------------------------------------------------------------------------
    def render_markdown_report(self, record: SampleMasterRecord) -> Path:
        """Render forensic report with adaptive handling for BAM direct ingestion."""
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
        ]

        if record.is_bam_mode:
            lines.extend([
                "---",
                "",
                "## 2. Input Mode & Sequencing Specification (Block 2)",
                "",
                "> ℹ️ **Direct BAM Passthrough**: The sample was imported directly from an aligned BAM file (`BAM_PASSTHROUGH`). Adapter trimming, pair-collapsing, and per-cycle FASTQ metrics were bypassed as not applicable.",
                "",
            ])
        else:
            lines.extend([
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
            ])

        lines.extend([
            "---",
            "",
            "## 3. Mitochondrial Functional Macro-Complexes Coverage (Block 6)",
            "",
            "| Functional Macro-Complex | rCRS Coordinates | Length | Mean Depth | Breadth $\\ge 1\\text{x}$ | Evaluation |",
            "| :--- | :--- | :--- | :--- | :--- | :--- |",
        ])

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
    # 6. Technical Explanatory Documentation (Contamination)
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
    # 7. Diagnostic Warnings & Error Codes Comprehensive Glossary Export
    # --------------------------------------------------------------------------
    def export_warnings_glossary(self) -> Path:
        """
        Exports an exhaustive, audit-ready reference glossary into results/master_summary/
        systematically documenting all diagnostic warnings, fatal error codes, gating
        actions, and VCF filters across Blocks 3 through 8.
        """
        glossary_path = self.dir_master / "WARNINGS_GLOSSARY.md"
        content = r"""# Forensic Mitochondrial Pipeline: Diagnostic Warnings, Errors & VCF Filter Glossary
*Pipeline Release: v1.2.0 | Master Reference Documentation*

---

## 1. Architectural Overview & Severity Scale

The analytical framework decouples routine biological variation from critical data integrity failures using a formal three-tier severity hierarchy:

* **CRITICAL (🔴)**: Severe data compromise, high exogenous contamination, extreme missing data, or cryptographic lineage discordance. The specimen fails automated release gates (`overall_status = "FAIL"`) and downstream sequence export is blocked.
* **WARNING (🟡)**: Measurable biological anomalies, elevated private polymorphisms, ambiguous phylogenetic placements, or borderline quality scores that warrant specialist review but do not invalidate the specimen (`overall_status = "WARN"`).
* **INFO (ℹ️)**: Benign observational annotations (e.g. sub-clade heteroplasmy refinement or standalone execution notes) providing deeper biological resolution without quality penalties (`overall_status = "PASS"`).

---

## 2. Block 3: Alignment, Deduplication & Coverage Codes

### Diagnostic Warnings (`AlignmentWarningCode`)

* **`LOW_MAPPING_RATE`** (`WARNING`):
  * **Trigger**: The fraction of total input reads successfully mapped to the rCRS mitochondrial contig is below expected library thresholds.
  * **Implication**: Low endogenous mitochondrial capture efficiency or high environmental background DNA burden.
  * **Action**: Review extraction and capture protocols. Inspect nuclear/environmental metagenomic fractions.
* **`HIGH_DUPLICATION_RATE`** (`WARNING`):
  * **Trigger**: Optical and PCR duplicate read pairs exceed `max_duplication_rate_fraction` (default: $\\ge 85.0\\%$).
  * **Implication**: Severe PCR over-amplification resulting from ultra-low starting endogenous molecular complexity.
  * **Action**: Ensure two-coordinate deduplication (`dedup_mode: both`) is active. Check post-deduplication callable depth.
* **`LOW_BREADTH_COVERAGE`** (`WARNING`):
  * **Trigger**: Percentage of mitochondrial bases covered at $\\ge 3\\text{x}$ depth is below threshold ($< 80.0\\%$).
  * **Implication**: Extensive coverage dropouts across the circular genome, compromising full mitogenome reconstruction.
  * **Action**: Cross-reference macro-complex coverage metrics to identify dropout regions (e.g. D-loop or GC-rich complexes).
* **`LOW_DEPTH_COVERAGE`** (`WARNING`):
  * **Trigger**: Mean sequencing depth across the mitochondrial contig is below the minimum threshold ($< 3.0\\text{x}$).
  * **Implication**: Statistical power for majority base consensus calling is marginal; heteroplasmy calling is disabled.
  * **Action**: Sequence at higher depth if additional extract is available, or restrict analysis to major consensus calls.
* **`UNRESOLVED_READ_GROUPS`** (`WARNING`):
  * **Trigger**: Input BAM lacks formal `@RG` header tags or contains malformed identifier fields.
  * **Implication**: Read provenance and lane-specific error tracking cannot be independently verified.
  * **Action**: Pipeline auto-injects default read groups (`ID:rg1`, `PL:ILLUMINA`, `LB:lib1`). Review upstream BAM headers.
* **`ASYMMETRIC_MATE_MAPPING`** (`WARNING`):
  * **Trigger**: Severe imbalance between forward (R1) and reverse (R2) mapping counts in paired-end libraries.
  * **Implication**: Technical sequencing artifact, asymmetric adapter ligation, or chimeric read formation.
  * **Action**: Inspect AdapterRemoval read-merging metrics and adapter contamination profiles.

### Fatal Errors (`AlignmentErrorCode`)

* **`NO_MAPPED_MT_READS`** (`FAIL`): Zero reads aligned to the target mitochondrial contig.
* **`REFERENCE_FASTA_NOT_FOUND`** (`FAIL`): Specified reference FASTA file does not exist on disk.
* **`REFERENCE_CONTIG_MISMATCH`** (`FAIL`): Contig name in alignment stream diverges from reference specification.
* **`ALIGNMENT_BINARY_NOT_FOUND`** (`FAIL`): BWA or Samtools executables missing from PATH or Conda environment.
* **`ALIGNMENT_EXECUTION_FAILED`** (`FAIL`): BWA process terminated with a non-zero exit code.
* **`SAMTOOLS_PIPELINE_FAILED`** (`FAIL`): Samtools pipe (`collate | fixmate -m | sort`) encountered an unrecoverable failure.
* **`CORRUPTED_OUTPUT_BAM`** (`FAIL`): Resulting BAM failed BGZF EOF block validation or index generation.
* **`INDEXING_FAILED`** (`FAIL`): Samtools failed to generate standard `.bai` or `.csi` index.

---

## 3. Block 4: aDNA Damage Assessment & Authenticity Codes

### Diagnostic Warnings (`DamageWarningCode`)

* **`LOW_DAMAGE_ON_NON_UDG`** (`WARNING`):
  * **Trigger**: Untreated library (`library_treatment: NONE`) exhibits terminal 5'-C>T deamination rates $< 3.0\\%$ at position 1.
  * **Implication**: Lack of characteristic post-mortem degradation kinetics; specimen behaves like modern DNA.
  * **Action**: Audit specimen archaeological context. Verify whether undocumented enzymatic repair occurred.
* **`UNEXPECTED_HIGH_DAMAGE_ON_UDG`** (`WARNING`):
  * **Trigger**: UDG-treated library (`library_treatment: UDG_FULL`) displays substantial terminal damage ($> 3.0\\%$).
  * **Implication**: Incomplete enzymatic uracil-DNA glycosylase cleavage or severe library preparation failure.
  * **Action**: Enable read-level PMD score filtering or quality score rescaling.
* **`LOW_READ_COUNT_FOR_DAMAGE_MODEL`** (`WARNING`):
  * **Trigger**: Fewer than `min_reads_for_damage_model` (default: $< 100$) aligned reads available for modeling.
  * **Implication**: Insufficient sample size to fit Bayesian damage models with statistical confidence.
  * **Action**: Statistical rescaling is skipped; pipeline uses raw quality scores.
* **`INSUFFICIENT_TERMINAL_COVERAGE`** (`WARNING`):
  * **Trigger**: Total nucleotide observations at terminal positions 1 and 2 are below threshold ($< 50$).
  * **Implication**: Terminal substitution frequencies are stochastic and susceptible to read-end noise.
  * **Action**: Verdict transitions to `INCONCLUSIVE`. Check fragment length distributions.
* **`UNKNOWN_LIBRARY_TREATMENT`** (`WARNING`):
  * **Trigger**: Library chemistry is unspecified (`library_treatment: UNKNOWN`).
  * **Implication**: Pipeline must assume conservative Non-UDG defaults with spatial distance penalties.
  * **Action**: Explicitly configure `library_treatment` in `default_config.yaml`.
* **`BAYESIAN_CONVERGENCE_FAILURE`** (`WARNING`):
  * **Trigger**: MapDamage2 MCMC optimization failed convergence criteria.
  * **Implication**: Bayesian parameter estimates ($\lambda$, $\delta_s$, $\delta_d$) are unverified.
  * **Action**: MapDamage falls back to empirical mismatch tables for quality rescaling.
* **`MISSING_MD_TAGS`** (`WARNING`):
  * **Trigger**: Input BAM lacks standard SAM `MD:Z` mismatch tags.
  * **Implication**: Alignment mismatch identification requires on-the-fly reference comparison.
  * **Action**: Pipeline auto-injects MD tags using `samtools calmd`.
* **`LOW_PMD_FRACTION`** (`WARNING`):
  * **Trigger**: Fraction of reads with PMD score $\\ge 3.0$ is below expectation ($< 10.0\\%$).
  * **Implication**: Weak single-read ancient DNA authentication support.
  * **Action**: Compare PMD fraction with terminal substitution rates.
* **`ASYMMETRIC_TERMINAL_DAMAGE`** (`WARNING`):
  * **Trigger**: Absolute difference between 5'-C>T and 3'-G>A terminal frequencies exceeds $15.0\\%$.
  * **Implication**: Library preparation bias, asymmetric end-repair, or single-stranded overhang anomalies.
  * **Action**: Inspect library protocol (e.g. single-stranded vs double-stranded library preparation).

### Fatal Errors (`DamageErrorCode`)

* **`MAPDAMAGE_EXECUTION_FAILED`** (`FAIL`): MapDamage2 Python/R scripts failed or timed out.
* **`PMDTOOLS_EXECUTION_FAILED`** (`FAIL`): PMDtools scored zero reads or crashed during execution.
* **`RESCALING_FAILED`** (`FAIL`): Quality score recalibration process failed to output a valid BAM.
* **`CORRUPTED_DAMAGE_OUTPUT`** (`FAIL`): MapDamage `misincorporation.txt` or length tables unparseable.
* **`INPUT_BAM_NOT_FOUND`** (`FAIL`): Deduplicated BAM required for Block 4 does not exist.

### Biological Authentication Verdicts (`AuthenticationVerdict`)

* **`AUTHENTIC_ANCIENT_SUPPORTED`**: Strong terminal deamination and fragmentation confirming ancient origin.
* **`AUTHENTICITY_UNSUPPORTED`**: Absence of deamination and long insert sizes; classified as modern contaminant.
* **`UDG_DEPLETED_EXPECTED`**: Expected flat deamination profile on certified UDG-full library.
* **`INCONCLUSIVE`**: Depth or terminal read observations insufficient to determine authenticity.

---

## 4. Block 5: Contamination Screening & Bayesian Estimation Codes

### Diagnostic Warnings (`ContaminationWarningCode`)

* **`LOW_COVERAGE_FOR_ESTIMATION`** (`WARNING`):
  * **Trigger**: Mean coverage depth is below the threshold required for robust estimation ($< 5.0\\text{x}$).
  * **Implication**: Sparse data increases uncertainty in both deamination-based and position-based Bayesian models.
  * **Action**: Estimation confidence downgraded to `LOW` or `UNRELIABLE`. Heteroplasmy calling restricted.
* **`BROAD_CREDIBLE_INTERVAL`** (`WARNING`):
  * **Trigger**: 95% Bayesian Credible Interval span exceeds threshold ($[CI_{upper} - CI_{lower}] > 4.0\\%$).
  * **Implication**: Wide posterior density indicates substantial uncertainty in the true contaminant burden.
  * **Action**: Confidence downgraded to `MODERATE` or `LOW`. Manual inspection of likelihood curves required.
* **`MCMC_CONVERGENCE_FAILURE`** (`WARNING`):
  * **Trigger**: Schmutzi EM/MCMC optimization failed to reach stationary contamination estimates.
  * **Implication**: Iterative likelihood surface was flat or multimodal.
  * **Action**: Contamination estimates flagged as `UNRELIABLE`.
* **`CONTAMINANT_DEAMINATION_DETECTED`** (`WARNING`):
  * **Trigger**: Schmutzi diagnostic log indicates detectable cytosine deamination within the contaminant fraction.
  * **Implication**: Cross-individual ancient sample contamination or historical archaeological handling admixture.
  * **Action**: Severe forensic alert. Inspect excavation and handling logs.
* **`MULTIPLE_CONTAMINANTS_SUSPECTED`** (`WARNING`):
  * **Trigger**: Evidence of multiple competing modern human haplogroups in the contaminant pool.
  * **Implication**: Sample was handled by multiple individuals without cleanroom containment.
  * **Action**: Review posterior contamination likelihood grid (`.cont`).
* **`DATABASE_NOT_CONFIGURED`** / **`DATABASE_NOT_FOUND`** (`WARNING`):
  * **Trigger**: 197-genome contaminant allele frequency database missing or misconfigured.
  * **Implication**: Position-based Bayesian joint mixture modeling cannot proceed.
  * **Action**: Configure valid database path or trigger automated fallback to DissentEngine.
* **`FALLBACK_TO_DISSENT`** (`WARNING`):
  * **Trigger**: Schmutzi crashed or encountered modern zero-damage `-nan` mathematical instability; auto-recovered.
  * **Implication**: Non-Bayesian allelic dissent screening was used to estimate contamination.
  * **Action**: Review DissentEngine TSV tables for empirical minor allele frequencies.
* **`LOW_INFORMATIVE_SITES`** (`WARNING`):
  * **Trigger**: Fewer than 10 informative divergent loci between endogenous consensus and contaminant database.
  * **Implication**: Statistical resolution to distinguish contaminant and endogenous lineages is limited.
  * **Action**: Contamination estimate confidence reduced.
* **`FAST_SCREEN_ESCALATION`** (`INFO`):
  * **Trigger**: Screening detected $\\ge 3$ dissenting loci or minor allele VAF $\\ge 2.0\\%$.
  * **Implication**: Sample triggered mandatory Phase 2 modeling. Normal workflow progression.
* **`FAST_SCREEN_INSUFFICIENT_DATA`** (`INFO`):
  * **Trigger**: Diagnostic panel coverage insufficient for pre-screening triage.
  * **Implication**: Conservative escalation to Phase 2 modeling enforced to guarantee safety.

### Fatal Errors (`ContaminationErrorCode`)

* **`SCHMUTZI_EXECUTION_FAILED`** (`FAIL`): Schmutzi binary or Perl wrapper terminated abnormally.
* **`DISSENT_EXECUTION_FAILED`** (`FAIL`): DissentEngine pileup iteration failed.
* **`FAST_SCREEN_FAILED`** (`FAIL`): FastContaminationScreen pileup parsing failed.
* **`PREPARER_FAILED`** (`FAIL`): BAM preparation (`calmd` or subsampling) failed.
* **`CORRUPTED_OUTPUT`** (`FAIL`): Contamination output files missing or empty.
* **`INPUT_BAM_NOT_FOUND`** (`FAIL`): BAM required for Block 5 is missing.
* **`REFERENCE_INCOMPATIBLE`** (`FAIL`): Reference contig mismatch between BAM and database.
* **`DATABASE_MISSING`** (`FAIL`): Critical database files unavailable and fallback disabled.

### Consensus Gating Directives (`ConsensusGatingAction`)

* **`STANDARD_CALLING`**: Contamination $\\le 2.0\\%$. Standard variant and heteroplasmy calling permitted.
* **`RESTRICT_HETEROPLASMY`**: Contamination $2.0\\% - 5.0\\%$. Elevated thresholds enforced (VAF $\\ge 25\\%$, depth $\\ge 20\\text{x}$).
* **`USE_ENDOGENOUS_SCHMUTZI_CONSENSUS`**: Contamination elevated; endogenous FASTA generated by Schmutzi deployed.
* **`PMD_FILTERING_REQUIRED`**: Contamination $5.0\\% - 10.0\\%$. Partition reads by PMD score prior to calling.
* **`EXCLUDE_OR_MANUAL_REVIEW`**: Contamination $> 10.0\\%$ or unassessed. Sample gated from release.

---

## 5. Block 6: Consensus, Heteroplasmy & VCF Filter Codes

The following filter codes appear in `results/vcf/<sample_id>/*.vcf.gz`:

| Filter Code | Destination Stream | Definition & Biological Impact | Action Required |
| :--- | :--- | :--- | :--- |
| **`PASS`** | `.variants.vcf.gz` | High-confidence fixed consensus variant or authentic heteroplasmy passing all filters. | Accepted for downstream phylogenetic and population genetic analysis. |
| **`FILTER_EXCLUDED_INTERVAL`** | `.subconsensus.vcf.gz` | Variant falls inside blacklisted artifact intervals (e.g. rCRS 3107N gap region 3105–3111 or HVS Poly-C tracts). | Automatically suppressed from final variant catalog to eliminate known reference artifacts. |
| **`FILTER_LOW_DEPTH`** | `.subconsensus.vcf.gz` | Callable locus depth below threshold (SNV $< 15\\text{x}$, InDel $< 20\\text{x}$). | Reject minor allele call; insufficient depth to rule out sequencing error. |
| **`FILTER_LOW_ALT_READS`** | `.subconsensus.vcf.gz` | Supporting alternative reads below threshold (SNV $< 3$, InDel $< 4$). | Reject candidate variant; stochastic single-molecule error cannot be ruled out. |
| **`FILTER_VAF_RANGE`** | `.subconsensus.vcf.gz` | Variant Allele Frequency falls outside valid heteroplasmy calling window ($< 10\\%$ or $> 80\\%$). | Sub-10% variants filtered as background noise; $>80\\%$ variants promoted to fixed consensus. |
| **`FILTER_STRAND_BIAS`** | `.subconsensus.vcf.gz` | Severe forward/reverse strand asymmetry (Fisher's exact test $p < 0.05$ and $< 1$ read per strand). | Reject call; classic marker of platform-specific sequencing or alignment artifact. |
| **`FILTER_DAMAGE_SUSPECT`** | `.subconsensus.vcf.gz` | Variant carries deamination chemistry (C>T / G>A) enriched at read ends and lacks internal molecule support. | Reject candidate call; post-mortem deamination artifact on damaged ancient DNA fragment. |
| **`FILTER_CONTAMINATION`** | `.subconsensus.vcf.gz` | Minor allele suppressed due to upstream contamination gating constraints (`RESTRICT_HETEROPLASMY`). | Suppressed to prevent reporting exogenous modern contaminant alleles as endogenous heteroplasmies. |
| **`FILTER_AMBIGUOUS`** | `.subconsensus.vcf.gz` | Variant located inside homopolymer slippage tract or contested multi-allelic site. | Excluded; polymerase stutter or reference alignment ambiguity. |
| **`FILTER_REFERENCE_MISMATCH`** | `.subconsensus.vcf.gz` | Genomic reference base at coordinate is non-canonical or legacy placeholder ('N'). | Masked; invalid locus for comparative base calling. |
| **`FILTER_NO_DATA`** | `.subconsensus.vcf.gz` | Zero callable reads observed across coordinate. | Coordinate masked to 'N' in canonical consensus FASTA. |
| **`FILTER_INSUFFICIENT_DATA`** | `.subconsensus.vcf.gz` | Upstream contamination estimation failed to provide a reliable gating threshold. | Strict conservative safety filtering enforced. |

---

## 6. Block 7: Phylogenetics, Node Resolution & Lineage Audit Codes

### Diagnostic Audit Codes (`AuditCode`)

* **`MISSING_DIAGNOSTIC_MUTATIONS`** (`WARNING`):
  * **Trigger**: Specimen lacks $\\ge 2$ expected diagnostic polymorphisms defining the assigned haplogroup node.
  * **Implication**: Potential lineage discordance, local coverage dropout at diagnostic loci, or back-mutation.
  * **Action**: Inspect `missing_diagnostic_mutations` in master JSON; check whether coverage dropped to zero at those exact sites.
* **`EXCESS_PRIVATE_MUTATIONS`** (`WARNING`):
  * **Trigger**: Specimen carries more private mutations than permitted by threshold ($> 5$).
  * **Implication**: Potential unrepaired post-mortem deaminations fixed into consensus, or rare uncharacterized lineage.
  * **Action**: Cross-reference private variants with spatial damage models and fragment-end distance distributions.
* **`DAMAGE_SUSPECT_PRIVATE_MUTATION`** (`WARNING`):
  * **Trigger**: One or more private mutations carry deamination chemistry (5'-C>T or 3'-G>A) and fail spatial internal molecule support.
  * **Implication**: Post-mortem damage artifact misidentified as private evolutionary variation.
  * **Action**: Perform read-level spatial inspection in BAM. Exclude from population genetic networks.
* **`MULTIPLE_TOP_HITS_AMBIGUOUS`** (`WARNING`):
  * **Trigger**: HaploGrep3 reports co-optimal candidate haplogroups with identical Kulczynski distance metrics.
  * **Implication**: Sequence lacks sufficient diagnostic coverage to discriminate between sister clades.
  * **Action**: Review co-optimal clade list in Master JSON. Check whether sub-consensus heteroplasmies provide resolution.
* **`HAPLOGROUP_INPUT_DISCORDANCE`** (`WARNING`):
  * **Trigger**: Consensus FASTA and VCF calls yield discordant haplogroup classifications.
  * **Implication**: IUPAC ambiguities or sub-consensus heteroplasmies altered tree traversal relative to major consensus.
  * **Action**: Inspect VCF variant calls for competing minor alleles.
* **`LOW_HAPLOGROUP_CONFIDENCE`** (`WARNING`):
  * **Trigger**: HaploGrep assignment Kulczynski quality score is below threshold ($< 0.80$).
  * **Implication**: Poor phylogenetic alignment; substantial proportion of expected diagnostic markers are missing or mutated.
  * **Action**: Verify whether consensus has high 'N' fraction or non-mitochondrial read contamination.
* **`SUBCLADE_HETEROPLASMY_REFINEMENT`** (`INFO`):
  * **Trigger**: Consensus FASTA is assigned to a basal node, while VCF calls refine classification to a derived sub-clade.
  * **Implication**: Sub-consensus heteroplasmies capture ongoing intra-individual micro-evolution.
  * **Action**: Informational notice. Confirm that heteroplasmic variants satisfy strand bias and damage filters.
* **`SUBCLADE_FASTA_SPECIFIC`** (`INFO`): Subclade-defining mutations present exclusively in fixed major consensus.
* **`HETEROPLASMY_SUBCLADE_SUPPORTED`** (`INFO`): Subclade classification corroborated by authenticated heteroplasmies.
* **`SUSPECT_HETEROPLASMY_SUBCLADE_REJECTED`** (`WARNING`):
  * **Trigger**: Subclade refinement driven by heteroplasmies that failed damage or contamination audits.
  * **Implication**: Putative derived branch assignment is an artifact of deamination or exogenous contamination.
  * **Action**: Refinement rejected; sample retained at basal consensus node.
* **`CONTAMINANT_AFFINITY_WARNING`** (`CRITICAL`):
  * **Trigger**: Private or minor alleles in sample share diagnostic motifs with the putative contaminant haplogroup.
  * **Implication**: Direct evidence of exogenous contamination polluting variant calls.
  * **Action**: Exclude shared minor alleles; review Bayesian contamination estimates.
* **`TERMINAL_DEL_SUSPECT`** (`WARNING`): Deletion located within terminal fragment end distance; suspected alignment artifact.
* **`INDEL_SOFT_CLIP_SUSPECT`** (`WARNING`): InDel call supported predominantly by soft-clipped alignments.
* **`EXCESS_MISSING_DATA`** (`CRITICAL`): Consensus sequence contains $> 15.0\\%$ masked 'N' bases; haplogroup rejected.

---

## 7. Block 8: Forensic Synthesis, Cryptographic Provenance & Lineage Codes

### Categories & Verification States

* **`PROVENANCE_HASH_MISMATCH`** (`CRITICAL`): Physical file SHA-256 differs from sealed cryptographic contract hash. Indicates file corruption or out-of-order execution.
* **`PROVENANCE_UNANCHORED`** (`WARNING` / `INFO`): Artifact exists but lacks cryptographic root-of-trust hash in execution manifest.
* **`MISSING_UPSTREAM_CONTRACT`** (`CRITICAL`): Upstream block contract JSON missing on disk. Synthesis pipeline cannot guarantee provenance.
* **`SEMANTIC_VALIDATION_ERROR`** (`CRITICAL`): BAM, FASTA, or VCF file failed deep bioinformatic format inspection (e.g. BAM coordinate sort tag missing or malformed VCF).
* **`LINEAGE_DISCORDANCE`** (`CRITICAL`): Entity hash mismatch across block boundaries (e.g. Block 5 analyzed a BAM distinct from Block 4 dedup BAM).
* **`UNTRUSTED_VCF_BLOCKED`** (`CRITICAL`): Mutation cataloging halted because variants VCF failed cryptographic provenance audit.
* **`UPSTREAM_CONTRACT_ANOMALY`** (`CRITICAL` / `WARNING`): Inverted confidence interval bounds or mathematical probability boundary violation ($> 1.0$) in contract payload.
* **`MISSING_CONTAMINATION_METRIC`** (`CRITICAL`): Contamination fraction reported as non-numeric or missing.
* **`SCHEMA_CORRUPTION`** (`CRITICAL`): Contract JSON file truncated or failing JSON parser syntax checks.
"""
        self._atomic_write_text(glossary_path, content.strip() + "\n")
        return glossary_path

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

        # 2. Export sequencing and trimming QC metrics (bypassed in BAM mode)
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

        # 6. Export positional summary TSV row (adaptive 42 or 52 columns)
        if self.config.emit_tsv_summary:
            emitted["summary_tsv"] = self.export_sample_summary_tsv(record)

        # 7. Render human-readable forensic Markdown report
        if self.config.emit_markdown_sheets:
            emitted["report_markdown"] = self.render_markdown_report(record)

        # 8. Export technical explanatory README for contamination deliverables
        readme_path = self.export_contamination_readme()
        emitted["contamination_readme"] = readme_path

        # 9. Export comprehensive diagnostic warnings, errors, and VCF filter glossary
        glossary_path = self.export_warnings_glossary()
        emitted["warnings_glossary"] = glossary_path

        logger.info(
            f"[{record.sample_id}] Packaging completed: {len(emitted)} artifacts emitted "
            f"across sample-partitioned procedural directories (status={record.overall_status})."
        )
        return emitted