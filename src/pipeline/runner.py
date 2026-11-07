"""
End-to-End Modular Ancient and Modern Mitochondrial DNA Analysis Pipeline Runner.

Orchestrates sequential execution of all 8 independent analytical steps:
1. Reference Genome Validation (src.steps.reference)
2. Preprocessing & Discovery (src.steps.discovery)
3. Alignment & Deduplication (src.steps.alignment)
4. aDNA Damage Quantification & Rescaling (src.steps.damage)
5. Contamination Estimation & Gating (src.steps.contamination)
6. Dual-Stream Consensus & Heteroplasmy Calling (src.steps.consensus)
7. Phylogenetic Haplogrouping & Lineage Audit (src.steps.haplogroup)
8. Forensic Synthesis & Root-of-Trust Packaging (src.steps.synthesis)
"""

from datetime import datetime, timezone
import json
import logging
import os
from pathlib import Path
import sys
import time
from typing import Any, Dict, List, Optional, Tuple, Union
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]

from src.steps.reference import ReferenceInfo, ReferenceStep
from src.steps.discovery import DataType, DiscoveryStep, SampleRecord
from src.steps.alignment import AlignmentResult, AlignmentStatus, AlignmentStep
from src.steps.damage import DamageResult, DamageStep
from src.steps.contamination import ContaminationBatchRunner, ContaminationResult, ContaminationStep, ResourceAllocation
from src.steps.consensus import ConsensusQC, ConsensusStep, GatingAction
from src.steps.haplogroup import Block7Contract, HaplogroupStep
from src.steps.synthesis import SampleMasterRecord, SynthesisOrchestrator, SynthesisStep
from src.utils import (
    calculate_sha256,
    compute_composite_hash,
    extract_input_paths,
    load_pipeline_config,
    PipelineConfig,
    resolve_path,
)
from src.pipeline.stage_tracker import StageTracker

logger = logging.getLogger("mtDNA_Pipeline.Runner")


class PipelineExecutionError(Exception):
    """Raised when an unrecoverable failure occurs in pipeline execution."""
    pass


class EndToEndPipelineRunner:
    """Master orchestrator executing the 8 modular pipeline steps."""

    def __init__(
        self,
        config_path: Optional[Path] = None,
        cpu_fraction_override: Optional[float] = None,
        host_cores_override: Optional[int] = None,
        data_type_override: Optional[str] = None,
        cache_mode_override: Optional[str] = None,
    ) -> None:
        self.config_path = config_path or (PROJECT_ROOT / "config" / "default_config.yaml")
        with open(self.config_path, "r", encoding="utf-8") as f:
            self.config = yaml.safe_load(f) or {}

        self.pipeline_config = PipelineConfig.from_yaml(self.config_path)

        # Configuration sections
        self.ref_cfg = self.config.get("reference", {})
        self.input_cfg = self.config.get("input", {})
        self.prep_cfg = self.config.get("preprocessing", {})
        self.aln_cfg = self.config.get("alignment", {})
        self.dedup_cfg = self.config.get("deduplication", {})
        self.filter_cfg = self.config.get("filtering", {})
        self.damage_cfg = self.config.get("damage", {})
        self.contam_cfg = self.config.get("contamination", {})
        self.cons_cfg = self.config.get("consensus", {})
        self.synth_cfg = self.config.get("synthesis", {})
        self.res_cfg = self.config.get("resources", {})
        self.out_cfg = self.config.get("output", {})

        # Override data format from CLI or fallback to YAML
        raw_type = data_type_override or self.input_cfg.get("data_type", "FASTQ")
        self.data_type = str(raw_type).upper()

        # Cache handling policy from YAML ('prompt', 'reuse', or 'overwrite')
        raw_cache_mode = cache_mode_override or self.out_cfg.get("cache_mode", "prompt")
        self.cache_mode = str(raw_cache_mode).lower()

        # Operational workspace vs final release deliverables
        self.results_dir = resolve_path(self.out_cfg.get("results_dir", "results"), PROJECT_ROOT / "results")
        self.operational_files_dir = resolve_path(
            self.out_cfg.get("operational_files_dir", "operational_files"),
            PROJECT_ROOT / "operational_files",
        )
        self.logs_dir = self.operational_files_dir / "logs"
        self.cache_dir = self.operational_files_dir / "cache"

        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.operational_files_dir.mkdir(parents=True, exist_ok=True)
        self.logs_dir.mkdir(parents=True, exist_ok=True)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # CPU resource allocation
        self.cpu_fraction = cpu_fraction_override or float(self.res_cfg.get("cpu_budget_fraction", 0.50))
        self.host_cores = host_cores_override or self.res_cfg.get("host_cores_override")

        self.resource_alloc = ContaminationBatchRunner.calculate_resource_allocation(
            threads_per_sample=int(self.contam_cfg.get("schmutzi", {}).get("threads", 4)),
            cpu_budget_fraction=self.cpu_fraction,
            host_cores_override=self.host_cores,
        )

        self.alignment_threads = int(self.aln_cfg.get("threads", self.resource_alloc.threads_per_schmutzi_instance))

        # Reference setup
        self.ref_fasta_path = resolve_path(self.ref_cfg.get("fasta_path"), PROJECT_ROOT / "data" / "reference" / "rCRS.fasta")
        self.expected_contig = self.ref_cfg.get("expected_contig", "NC_012920.1")
        self.expected_length = int(self.ref_cfg.get("expected_length", 16569))
        self.reference_info: Optional[ReferenceInfo] = None
        self.ref_sha256: str = ""

    def _resolve_path(self, path_str: Optional[str], default: Path) -> Path:
        return resolve_path(path_str, default)

    # --------------------------------------------------------------------------
    # Cache Policy Handling
    # --------------------------------------------------------------------------
    def handle_cache_policy(
        self,
        cache_mode: str = "prompt",
        target_sample_ids: Optional[List[str]] = None,
    ) -> str:
        """Enforces cache policy ('prompt', 'reuse', or 'overwrite')."""
        detected_cache = StageTracker.inspect_operational_cache(
            self.operational_files_dir,
            target_sample_ids=target_sample_ids,
        )

        if not detected_cache:
            logger.info("No pre-existing cache artifacts detected in operational_files/.")
            return "reuse"

        mode = str(cache_mode or "prompt").lower()

        if mode == "overwrite":
            logger.info("Cache policy set to OVERWRITE: purging existing cached artifacts...")
            StageTracker.purge_operational_cache(
                self.operational_files_dir,
                results_dir=self.results_dir,
                target_sample_ids=target_sample_ids,
            )
            return "overwrite"

        if mode == "reuse":
            logger.info(f"Cache policy set to REUSE: retaining cached artifacts for {len(detected_cache)} sample(s).")
            return "reuse"

        if not sys.stdin.isatty():
            logger.info(f"Non-interactive terminal detected with {len(detected_cache)} cached sample(s); defaulting to REUSE.")
            return "reuse"

        print("\n" + "!" * 115)
        print(f"⚠️  EXISTING CACHE DETECTED IN: {self.operational_files_dir}")
        print("-" * 115)
        for s_id, stages in detected_cache.items():
            stages_str = ", ".join(stages)
            print(f"   • Sample '{s_id}': [{stages_str}]")
        print("-" * 115)
        print("How would you like to proceed?")
        print("   [1] Reuse cache     - Resume execution using valid cached stages (Default)")
        print("   [2] Overwrite cache - Purge cached artifacts and recompute everything from scratch")
        print("   [3] Stop / Abort    - Cancel execution and exit immediately")
        print("!" * 115)

        while True:
            try:
                user_choice = input("\nSelect option [1/2/3] (press Enter for 1): ").strip().lower()
            except EOFError:
                user_choice = "1"

            if user_choice in ("", "1", "reuse", "r", "u"):
                logger.info("User selected: REUSE existing cache.")
                return "reuse"
            elif user_choice in ("2", "overwrite", "o", "purge", "p"):
                logger.info("User selected: OVERWRITE cache. Purging operational cache...")
                StageTracker.purge_operational_cache(
                    self.operational_files_dir,
                    results_dir=self.results_dir,
                    target_sample_ids=target_sample_ids,
                )
                return "overwrite"
            elif user_choice in ("3", "stop", "abort", "exit", "q"):
                logger.info("User selected: STOP/ABORT. Exiting pipeline.")
                sys.exit(0)
            else:
                print("Invalid selection. Please enter 1 (Reuse), 2 (Overwrite), or 3 (Stop).")

    # --------------------------------------------------------------------------
    # Step Invocations
    # --------------------------------------------------------------------------
    def validate_reference(self) -> ReferenceInfo:
        """Step 1: Reference Genome Validation."""
        self.reference_info = ReferenceStep.run(
            fasta_path=self.ref_fasta_path,
            expected_contig=self.expected_contig,
            expected_length=self.expected_length,
            auto_index=self.ref_cfg.get("auto_index", True),
            log_dir=self.logs_dir / "reference",
        )
        self.ref_sha256 = self.reference_info.file_sha256 or calculate_sha256(self.ref_fasta_path)
        return self.reference_info

    def run_discovery_and_preprocessing(self, input_dir: Path) -> List[SampleRecord]:
        """Step 2: FASTQ Discovery & Preprocessing."""
        return DiscoveryStep.run(
            input_dir=input_dir,
            config=self.pipeline_config,
            operational_files_dir=self.operational_files_dir,
            cache_dir=self.cache_dir,
        )

    def run_alignment(self, sample: SampleRecord) -> AlignmentResult:
        """Step 3: BWA Alignment & Deduplication."""
        return AlignmentStep.run(
            sample=sample,
            reference_info=self.reference_info,
            config=self.pipeline_config,
            operational_files_dir=self.operational_files_dir,
            cache_dir=self.cache_dir,
            logs_dir=self.logs_dir,
            threads=self.alignment_threads,
        )

    def run_damage(self, sample_id: str, aln_res: AlignmentResult) -> Tuple[DamageResult, Path]:
        """Step 4: Damage Quantification & Rescaling."""
        return DamageStep.run(
            sample_id=sample_id,
            aln_res=aln_res,
            reference_info=self.reference_info,
            config=self.pipeline_config,
            operational_files_dir=self.operational_files_dir,
            cache_dir=self.cache_dir,
            logs_dir=self.logs_dir,
            threads=self.alignment_threads,
        )

    def run_contamination_and_report(self, unrescaled_bams: List[Tuple[str, Path, bool]]) -> Tuple[Path, Path]:
        """Step 5: Two-Tier Contamination Estimation."""
        return ContaminationStep.run(
            unrescaled_bams=unrescaled_bams,
            reference_info=self.reference_info,
            config=self.pipeline_config,
            operational_files_dir=self.operational_files_dir,
            cache_dir=self.cache_dir,
            resource_alloc=self.resource_alloc,
        )

    def run_consensus(
        self,
        sample_id: str,
        b4_contract_path: Path,
        b5_contract_path: Path,
    ) -> Tuple[ConsensusQC, Path, Path, Path]:
        """Step 6: Consensus & Heteroplasmy Calling."""
        return ConsensusStep.run(
            sample_id=sample_id,
            b4_contract_path=b4_contract_path,
            b5_contract_path=b5_contract_path,
            reference_info=self.reference_info,
            config=self.pipeline_config,
            operational_files_dir=self.operational_files_dir,
            cache_dir=self.cache_dir,
        )

    def run_haplogroup(
        self,
        sample_id: str,
        canonical_fasta: Path,
        variants_vcf: Path,
        b4_contract_path: Path,
        b5_contract_path: Path,
        b6_contract_path: Path,
    ) -> Block7Contract:
        """Step 7: Haplogroup Assignment & Cross-Block Audit."""
        return HaplogroupStep.run(
            sample_id=sample_id,
            canonical_fasta=canonical_fasta,
            variants_vcf=variants_vcf,
            b4_contract_path=b4_contract_path,
            b5_contract_path=b5_contract_path,
            b6_contract_path=b6_contract_path,
            reference_info=self.reference_info,
            config=self.pipeline_config,
            operational_files_dir=self.operational_files_dir,
            cache_dir=self.cache_dir,
        )

    # --------------------------------------------------------------------------
    # Main Cohort Execution Loop
    # --------------------------------------------------------------------------
    def execute_cohort(self, input_dir: Path, cache_mode: Optional[str] = None) -> None:
        """Executes complete cohort analysis across all 8 pipeline steps."""
        t_global0 = time.perf_counter()
        effective_cache_mode = (cache_mode or self.cache_mode).lower()

        print("=" * 115)
        print("🚀 EXECUTING COMPLETE MODULAR mtDNA PIPELINE (BLOCKS 1 -> 8)")
        print(f"[*] Total CPU Allocation : {self.resource_alloc.total_allocated_threads}/{self.resource_alloc.total_host_cores} cores ({self.resource_alloc.cpu_budget_fraction:.0%})")
        print(f"[*] Upstream Threads/Job : {self.alignment_threads} threads")
        print(f"[*] Schmutzi Workers     : {self.resource_alloc.max_schmutzi_workers} slots ({self.resource_alloc.threads_per_schmutzi_instance} threads/slot)")
        print(f"[*] Input Format (YAML)  : {self.data_type}")
        print(f"[*] Cache Policy (YAML)  : {effective_cache_mode.upper()}")
        print("=" * 115)

        # Pre-flight Cache Inspection
        self.handle_cache_policy(cache_mode=effective_cache_mode)

        # Step 1: Reference Genome Validation
        self.validate_reference()

        # Step 2: Discovery & Preprocessing
        samples = self.run_discovery_and_preprocessing(input_dir)
        cohort_sample_ids = [s.sample_id for s in samples]

        # Step 3 & 4: Alignment & Damage Profiling
        logger.info("=" * 80)
        logger.info("▶ BLOCKS 3 & 4: ALIGNMENT (bwa-aln) & aDNA DAMAGE PROFILING")
        logger.info("=" * 80)

        unrescaled_bams_for_block5: List[Tuple[str, Path, bool]] = []
        b4_contracts: Dict[str, Path] = {}
        sample_types: Dict[str, str] = {}
        failed_samples: Dict[str, str] = {}

        for sample in samples:
            s_id = sample.sample_id
            logger.info(f"[*] Processing Sample: {s_id}")

            try:
                aln_res = self.run_alignment(sample)
                final_bam = (
                    (aln_res.outputs.dedup_bam if aln_res.outputs and aln_res.outputs.dedup_bam and aln_res.outputs.dedup_bam.is_file() else None)
                    or (aln_res.outputs.filtered_bam if aln_res.outputs and aln_res.outputs.filtered_bam and aln_res.outputs.filtered_bam.is_file() else None)
                    or (self.operational_files_dir / "alignment" / s_id / f"{s_id}.filtered.dedup.bam")
                )

                dam_res, b4_contract = self.run_damage(s_id, aln_res)
                b4_contracts[s_id] = b4_contract

                verdict = DamageStep.extract_verdict_with_contract_fallback(dam_res, b4_contract)
                is_ancient = "AUTHENTIC_ANCIENT" in verdict.upper()
                sample_types[s_id] = "ANCIENT" if is_ancient else "MODERN"
                unrescaled_bams_for_block5.append((s_id, final_bam, is_ancient))

                logger.info(f"[{s_id}] Block 4 Verdict: {verdict} (Classified: {sample_types[s_id]})")
            except Exception as err:
                logger.critical(f"[{s_id}] Blocks 3/4 failure: {err}", exc_info=True)
                failed_samples[s_id] = f"Blocks 3/4: {err}"

        active_sample_ids = [sid for sid in cohort_sample_ids if sid not in failed_samples]
        if not active_sample_ids:
            raise PipelineExecutionError("All samples failed upstream pipeline stages.")

        # Step 5: Contamination Estimation
        tsv_contam, json_contam = self.run_contamination_and_report(unrescaled_bams_for_block5)

        # Step 6: Consensus & Heteroplasmy Calling
        logger.info("=" * 80)
        logger.info("▶ BLOCK 6: DUAL-STREAM CONSENSUS & HETEROPLASMY CALLING")
        logger.info("=" * 80)

        cohort_qc: Dict[str, ConsensusQC] = {}
        b6_contracts: Dict[str, Path] = {}
        consensus_fastas: Dict[str, Path] = {}
        variants_vcfs: Dict[str, Path] = {}

        for s_id in active_sample_ids:
            b4_c = b4_contracts[s_id]
            b5_c = self.operational_files_dir / "contamination" / s_id / f"{s_id}.block5_contract.json"

            logger.info(f"[*] Block 6 Calling: {s_id}")
            try:
                qc_res, can_fa, vcf_gz, b6_c = self.run_consensus(s_id, b4_c, b5_c)
                cohort_qc[s_id] = qc_res
                b6_contracts[s_id] = b6_c
                consensus_fastas[s_id] = can_fa
                variants_vcfs[s_id] = vcf_gz

                logger.info(
                    f"[{s_id}] Block 6 Complete -> Depth: {qc_res.mean_depth:.1f}x | "
                    f"Covered: {qc_res.covered_positions_fraction:.1%} | "
                    f"Het Sites (PASS/Total): {qc_res.heteroplasmy_sites_pass}/{qc_res.heteroplasmy_sites_total}"
                )
            except Exception as err:
                logger.critical(f"[{s_id}] Block 6 failed: {err}", exc_info=True)
                failed_samples[s_id] = f"Block 6: {err}"

        ConsensusStep.export_consensus_cohort_summary(cohort_qc, self.operational_files_dir / "consensus")

        # Step 7: Haplogrouping & Lineage Validation
        logger.info("=" * 80)
        logger.info("▶ BLOCK 7: HAPLOGROUPING & PHYLOGENETIC VALIDATION")
        logger.info("=" * 80)

        cohort_b7_contracts: Dict[str, Block7Contract] = {}
        b7_active_ids = [sid for sid in active_sample_ids if sid not in failed_samples]

        for s_id in b7_active_ids:
            logger.info(f"[*] Block 7 Classification: {s_id}")
            try:
                b7_contract = self.run_haplogroup(
                    sample_id=s_id,
                    canonical_fasta=consensus_fastas[s_id],
                    variants_vcf=variants_vcfs[s_id],
                    b4_contract_path=b4_contracts[s_id],
                    b5_contract_path=self.operational_files_dir / "contamination" / s_id / f"{s_id}.block5_contract.json",
                    b6_contract_path=b6_contracts[s_id],
                )
                cohort_b7_contracts[s_id] = b7_contract

                logger.info(
                    f"[{s_id}] Block 7 Assigned: Primary={b7_contract.primary_haplogroup} | "
                    f"Resolved={b7_contract.resolved_haplogroup} | "
                    f"Score={b7_contract.quality_score:.2f} | "
                    f"Status={b7_contract.haplogroup_status.value}"
                )
            except Exception as err:
                logger.critical(f"[{s_id}] Block 7 failed: {err}", exc_info=True)
                failed_samples[s_id] = f"Block 7: {err}"

        HaplogroupStep.export_haplogroup_cohort_summary(cohort_b7_contracts, self.operational_files_dir / "haplogroup")

        # Step 8: Forensic Synthesis & Cohort Reporting
        cohort_deliverables = SynthesisStep.run(
            sample_ids=cohort_sample_ids,
            reference_info=self.reference_info,
            config=self.pipeline_config,
            operational_files_dir=self.operational_files_dir,
            results_dir=self.results_dir,
        )

        elapsed = time.perf_counter() - t_global0

        print("\n" + "=" * 125)
        print(f"🎉 COMPLETE PIPELINE (BLOCKS 1 -> 8) FINISHED IN {elapsed:.2f}s")
        print("=" * 125)
        print(f"{'Sample ID':<16} | {'Type':<8} | {'Resolved HG':<12} | {'Score':<6} | {'Status':<6} | {'Diagnostic Warnings'}")
        print("-" * 125)

        for s_id in cohort_sample_ids:
            c7 = cohort_b7_contracts.get(s_id)
            if c7:
                b_type = sample_types.get(s_id, "ANCIENT")
                rhg = c7.resolved_haplogroup
                qscore = f"{c7.quality_score:.2f}"
                ov_stat = c7.overall_validation_status.value
                w_str = "Clean"
                print(f"{s_id:<16} | {b_type:<8} | {rhg:<12} | {qscore:<6} | {ov_stat:<6} | {w_str}")
            else:
                print(f"{s_id:<16} | {'FAILED':<8} | {'UNKNOWN':<12} | {'NA':<6} | {'FAIL':<6} | Pipeline processing error")

        if failed_samples:
            print("\n" + "!" * 125)
            print("⚠️ SAMPLES WITH UNRESOLVED FAILURES:")
            for f_id, f_err in failed_samples.items():
                print(f"    - {f_id}: {f_err}")
            print("!" * 125)
            raise PipelineExecutionError(
                f"Cohort execution completed with {len(failed_samples)} failed sample(s): {list(failed_samples.keys())}"
            )

        print("=" * 125)
        print("\n[*] Exported Master Deliverables (Block 8):")
        for label, path in cohort_deliverables.items():
            print(f"    - {label:<28}: {path}")
        print("=" * 125 + "\n")
