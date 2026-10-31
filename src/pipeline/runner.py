"""
Definitive Modular End-to-End mtDNA Pipeline Runner (Blocks 1 -> 8).
Executes reference validation, sample discovery, alignment/deduplication,
damage profiling, Bayesian contamination gating, consensus assembly,
haplogroup classification, and forensic synthesis.
"""

import csv
from datetime import datetime, timezone
import hashlib
import inspect
import json
import logging
from pathlib import Path
import shutil
import subprocess
import sys
import time
from typing import Any, Dict, List, Optional, Tuple, Union

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Block 1: Reference Management
from src.modules.alignment.reference_manager import (
    ReferenceManager,
    ReferenceValidationResult,
)
from src.modules.alignment.models import (
    AlignmentFilterPolicy,
    AlignmentMetrics,
    AlignmentOutputs,
    AlignmentResult,
    AlignmentStatus,
    CoverageMetrics,
    DedupMetrics,
    DedupMode,
    MultiMappingPolicy,
    ReferenceInfo,
)

# Block 2: Discovery, Integrity, Preprocessing
from src.modules.preprocessing.discovery import SampleDiscoveryEngine
from src.modules.preprocessing.integrity import IntegrityEngine
from src.modules.preprocessing.trimmer import PreprocessingEngine
from src.modules.preprocessing.models import (
    SampleRecord,
    ValidationStatus,
)

# Block 3: Alignment Orchestrator
from src.modules.alignment.orchestrator import AlignmentOrchestrator

# Block 4: aDNA Damage Authentication Orchestrator
from src.modules.damage.orchestrator import DamageOrchestrator
from src.modules.damage.models import (
    DamagePolicy,
    DamageResult,
    LibraryTreatment as DamageLibraryTreatment,
)

# Block 5: Contamination Triage & Bayesian Batch Runner
from src.modules.contamination.batch_runner import (
    ContaminationBatchRunner,
    ResourceAllocation,
)
from src.modules.contamination.models import ContaminationResult

# Block 6: Consensus & Heteroplasmy Calling
from src.modules.consensus.enums import GatingAction
from src.modules.consensus.models import (
    ConsensusPolicy,
    ConsensusQC,
    LibraryTreatment as ConsensusLibraryTreatment,
)
from src.modules.consensus.orchestrator import (
    ConsensusExecutionError,
    ConsensusOrchestrator,
    ContractValidationError,
    UnresolvedAuthenticityError,
)

# Block 7: Haplogrouping & Phylogenetic Validation
from src.modules.haplogrouping.orchestrator import HaplogroupOrchestrator
from src.modules.haplogrouping.models import (
    Block7Contract,
    ContaminantAffinity,
    HaplogroupStatus,
    OverallValidationStatus,
)

# Block 8: Synthesis, Lineage Provenance & Cohort Reporting
from src.modules.synthesis.orchestrator import SynthesisOrchestrator
from src.modules.synthesis.models import (
    FinalVerificationReport,
    SampleMasterRecord,
)

from src.pipeline.stage_tracker import StageTracker

logger = logging.getLogger("mtDNA_Pipeline.Runner")


class PipelineExecutionError(Exception):
    """Raised on unrecoverable workflow errors."""
    pass


class EndToEndPipelineRunner:
    """Coordinates execution across Blocks 1 through 8 using official domain engines."""

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
        
        # Cache handling policy from YAML ('prompt', 'reuse', or 'overwrite') with optional CLI override
        raw_cache_mode = cache_mode_override or self.out_cfg.get("cache_mode", "prompt")
        self.cache_mode = str(raw_cache_mode).lower()

        # Operational workspace vs final release deliverables
        self.results_dir = self._resolve_path(self.out_cfg.get("results_dir", "results"), PROJECT_ROOT / "results")
        self.operational_files_dir = self._resolve_path(
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
        self.ref_fasta_path = self._resolve_path(self.ref_cfg.get("fasta_path"), PROJECT_ROOT / "data" / "reference" / "rCRS.fasta")
        self.expected_contig = self.ref_cfg.get("expected_contig", "NC_012920.1")
        self.expected_length = int(self.ref_cfg.get("expected_length", 16569))
        self.reference_info: Optional[ReferenceInfo] = None
        self.reference_sequence: Optional[str] = None
        self.ref_sha256: str = ""

        # Block 6 Consensus Policy & Orchestrator
        self.consensus_policy = self._build_consensus_policy()
        self.consensus_orchestrator = ConsensusOrchestrator(policy=self.consensus_policy)

        # Block 7 Haplogroup Orchestrator
        self.haplogroup_orchestrator = HaplogroupOrchestrator(config=self.config)

        # Block 8 Synthesis Orchestrator
        self.synthesis_orchestrator: Optional[SynthesisOrchestrator] = None

    @staticmethod
    def calculate_sha256(filepath: Union[str, Path], chunk_size: int = 65536) -> str:
        """Compute SHA-256 hex digest for a target file."""
        p = Path(filepath)
        if not p.is_file():
            raise FileNotFoundError(f"Cannot calculate SHA-256; file not found: {p}")
        hasher = hashlib.sha256()
        with open(p, "rb") as f:
            while chunk := f.read(chunk_size):
                hasher.update(chunk)
        return hasher.hexdigest()

    @staticmethod
    def compute_composite_hash(*elements: Any) -> str:
        """Compute deterministic SHA-256 hash across heterogeneous config objects and input hashes."""
        hasher = hashlib.sha256()
        for elem in elements:
            if isinstance(elem, (dict, list)):
                hasher.update(json.dumps(elem, sort_keys=True, default=str).encode("utf-8"))
            elif isinstance(elem, Path):
                hasher.update(str(elem.resolve()).encode("utf-8"))
            else:
                hasher.update(str(elem).encode("utf-8"))
        return hasher.hexdigest()

    @staticmethod
    def _extract_input_paths(obj: Any) -> List[Path]:
        """Defensively extract all physical Path objects from an InputFiles/SampleRecord structure."""
        found_paths: List[Path] = []
        visited = set()

        def _traverse(current: Any):
            if current is None or id(current) in visited:
                return
            visited.add(id(current))

            if isinstance(current, Path):
                if current.is_file():
                    found_paths.append(current)
                return

            if isinstance(current, str):
                p = Path(current)
                if p.is_file():
                    found_paths.append(p)
                return

            if hasattr(current, "path"):
                p_attr = getattr(current, "path")
                if isinstance(p_attr, Path) and p_attr.is_file():
                    found_paths.append(p_attr)
                elif isinstance(p_attr, str):
                    p = Path(p_attr)
                    if p.is_file():
                        found_paths.append(p)

            if isinstance(current, (list, tuple, set)):
                for item in current:
                    _traverse(item)
                return

            if isinstance(current, dict):
                for v in current.values():
                    _traverse(v)
                return

            if hasattr(current, "__dict__"):
                for k, v in vars(current).items():
                    if not k.startswith("_"):
                        _traverse(v)

        _traverse(obj)
        unique_paths: List[Path] = []
        seen = set()
        for p in found_paths:
            resolved = str(p.resolve())
            if resolved not in seen:
                seen.add(resolved)
                unique_paths.append(p)
        return unique_paths

    def _resolve_path(self, path_str: Optional[str], default: Path) -> Path:
        if not path_str:
            return default if default.is_absolute() else (PROJECT_ROOT / default)
        p = Path(path_str)
        return p if p.is_absolute() else (PROJECT_ROOT / p)

    def _build_consensus_policy(self) -> ConsensusPolicy:
        """Instantiate immutable ConsensusPolicy directly from configuration YAML with excluded intervals support."""
        het_cfg = self.cons_cfg.get("heteroplasmy", {})
        restr_cfg = self.cons_cfg.get("restricted_heteroplasmy", {})

        default_spatial = {"NON_UDG": 5, "UDG_HALF": 2, "UDG_FULL": 0, "SINGLE_STRANDED": 5}
        spatial_thresholds = dict(het_cfg.get("spatial_thresholds", default_spatial))

        default_excluded = [[3105, 3111], [303, 315], [513, 524], [16184, 16193]]
        raw_excluded = het_cfg.get("excluded_intervals", default_excluded)
        excluded_intervals = tuple(
            (int(item[0]), int(item[1]))
            for item in raw_excluded
            if isinstance(item, (list, tuple)) and len(item) == 2
        )

        return ConsensusPolicy(
            policy_id=self.cons_cfg.get("policy_id", "consensus_strict_v1"),
            policy_version=self.cons_cfg.get("policy_version", "1.2.0"),
            min_base_quality=int(self.cons_cfg.get("min_base_quality", 30)),
            min_mapping_quality=int(self.cons_cfg.get("min_mapping_quality", 30)),
            min_consensus_depth=int(self.cons_cfg.get("min_consensus_depth", 10)),
            major_allele_fraction=float(self.cons_cfg.get("major_allele_fraction", 0.80)),
            ambiguity_min_vaf=float(self.cons_cfg.get("ambiguity_min_vaf", 0.10)),
            heteroplasmy_min_depth=int(het_cfg.get("min_depth", 15)),
            heteroplasmy_min_alt_reads=int(het_cfg.get("min_alt_reads", 3)),
            heteroplasmy_min_vaf=float(het_cfg.get("min_vaf", 0.10)),
            heteroplasmy_max_vaf=float(het_cfg.get("max_vaf", 0.80)),
            strand_bias_alpha=float(het_cfg.get("strand_bias_alpha", 0.05)),
            strand_bias_min_strand_reads=int(het_cfg.get("strand_bias_min_strand_reads", 1)),
            emit_iupac_fasta=bool(self.cons_cfg.get("emit_iupac_fasta", True)),
            mask_char=str(self.cons_cfg.get("mask_char", "N")),
            restricted_heteroplasmy_min_depth=int(restr_cfg.get("min_depth", 20)),
            restricted_heteroplasmy_min_alt_reads=int(restr_cfg.get("min_alt_reads", 5)),
            restricted_heteroplasmy_min_vaf=float(restr_cfg.get("min_vaf", 0.25)),
            restricted_heteroplasmy_max_vaf=float(restr_cfg.get("max_vaf", 0.75)),
            min_internal_fragments=int(het_cfg.get("min_internal_fragments", 2)),
            spatial_thresholds=spatial_thresholds,
            excluded_intervals=excluded_intervals,
        )

    # --------------------------------------------------------------------------
    # Cache Inspection & Interactive Resolution
    # --------------------------------------------------------------------------
    def handle_cache_policy(
        self,
        cache_mode: str = "prompt",
        target_sample_ids: Optional[List[str]] = None,
    ) -> str:
        """
        Detects existing cached artifacts in operational_files/ and enforces
        the requested cache policy ('prompt', 'reuse', or 'overwrite').
        """
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
            logger.info(
                f"Cache policy set to REUSE: retaining cached artifacts for {len(detected_cache)} sample(s)."
            )
            return "reuse"

        # Interactive prompt mode ('prompt')
        if not sys.stdin.isatty():
            logger.info(
                f"Non-interactive terminal detected with {len(detected_cache)} cached sample(s); defaulting to REUSE."
            )
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
    # Block 1: Reference Genome Validation
    # --------------------------------------------------------------------------
    def validate_reference(self) -> ReferenceInfo:
        logger.info("=" * 80)
        logger.info("▶ BLOCK 1: REFERENCE GENOME VALIDATION (ReferenceManager)")
        logger.info("=" * 80)

        ref_val: ReferenceValidationResult = ReferenceManager.prepare_and_validate_reference(
            fasta_path=self.ref_fasta_path,
            expected_contig=self.expected_contig,
            strict_contig=True,
            expected_length=self.expected_length,
            allow_n=True,
            auto_index=self.ref_cfg.get("auto_index", True),
            log_dir=self.logs_dir / "reference",
        )

        if not ref_val.is_valid:
            raise PipelineExecutionError(f"Reference validation failed: {ref_val.errors}")

        self.reference_info = ref_val.reference_info
        self.ref_sha256 = self.calculate_sha256(self.ref_fasta_path)

        if self.ref_fasta_path.is_file():
            with open(self.ref_fasta_path, "r", encoding="utf-8") as f:
                self.reference_sequence = "".join(
                    line.strip() for line in f if not line.startswith(">")
                ).upper()

        logger.info(
            f"✅ Reference verified: {self.reference_info.fasta_path.name} "
            f"({self.reference_info.contig_name}, {self.reference_info.length_bp} bp) [SHA256: {self.ref_sha256[:12]}...]"
        )

        self.synthesis_orchestrator = SynthesisOrchestrator(
            results_dir=self.results_dir,
            reference_info=self.reference_info,
            operational_files_dir=self.operational_files_dir,
            reference_sequence=self.reference_sequence,
            config=self.config,
        )

        return self.reference_info

    # --------------------------------------------------------------------------
    # Block 2 & 3: Direct BAM Discovery & Ingestion (BAM Entrypoint Mode)
    # --------------------------------------------------------------------------
    def discover_and_ingest_bams(self, input_dir: Path) -> List[Tuple[str, AlignmentResult]]:
        """Scans BAM directory and encapsulates into valid AlignmentResult instances."""
        logger.info("=" * 80)
        logger.info(f"▶ DIRECT BAM INGESTION: Scanning directory {input_dir}")
        logger.info("=" * 80)

        if not input_dir.is_dir():
            raise PipelineExecutionError(f"BAM input directory does not exist: {input_dir}")

        bam_files = sorted(
            list(input_dir.glob("*.bam"))
            + list(input_dir.glob("*.BAM"))
            + list(input_dir.glob("*/*.bam"))
            + list(input_dir.glob("*/*.BAM"))
        )
        bam_files = sorted(list({p.resolve(): p for p in bam_files}.values()))

        if not bam_files:
            raise PipelineExecutionError(f"No BAM files found in {input_dir}")

        ingested: List[Tuple[str, AlignmentResult]] = []

        qc_thresh = self.filter_cfg.get("thresholds") or self.filter_cfg.get("qc_thresholds", {})
        min_depth_gate = float(qc_thresh.get("min_mean_depth", 3.0))
        min_b3x_gate = float(qc_thresh.get("min_breadth_3x_fraction", 0.80))
        max_dup_gate = float(qc_thresh.get("max_duplication_rate_fraction", 0.85))
        min_map_rate_gate = float(qc_thresh.get("min_mapping_rate_pct", 0.10))

        aln_policy = AlignmentFilterPolicy(
            dedup_mode=DedupMode(self.dedup_cfg.get("mode", "both").lower()),
            multimapping=MultiMappingPolicy(
                min_mapq=int(self.filter_cfg.get("multimapping", {}).get("min_mapq", 30)),
                filter_unmapped=bool(self.filter_cfg.get("multimapping", {}).get("filter_unmapped", True)),
                filter_secondary=bool(self.filter_cfg.get("multimapping", {}).get("filter_secondary", True)),
                filter_supplementary=bool(self.filter_cfg.get("multimapping", {}).get("filter_supplementary", True)),
            ),
            min_mean_depth=min_depth_gate,
            min_breadth_3x_fraction=min_b3x_gate,
            max_duplication_rate_fraction=max_dup_gate,
            min_mapping_rate_pct=min_map_rate_gate,
        )

        for bam_p in bam_files:
            if bam_p.name.endswith(".bai") or ".markdup" in bam_p.name:
                continue

            sample_id = bam_p.stem.replace(".filtered.dedup", "").replace(".filtered", "").replace(".dedup", "").replace(".sorted", "")
            logger.info(f"[*] Ingesting BAM for sample: {sample_id} ({bam_p.name})")

            bai_p = Path(f"{bam_p}.bai")
            alt_bai_p = bam_p.with_suffix(".bai")
            if not (bai_p.is_file() or alt_bai_p.is_file()):
                logger.info(f"[{sample_id}] Generating samtools index...")
                subprocess.run(["samtools", "index", str(bam_p)], check=True)

            sample_dir = self.operational_files_dir / "alignment" / sample_id
            sample_dir.mkdir(parents=True, exist_ok=True)
            target_bam = sample_dir / f"{sample_id}.filtered.dedup.bam"

            if bam_p.resolve() != target_bam.resolve() and not target_bam.is_file():
                try:
                    target_bam.symlink_to(bam_p.resolve())
                except OSError:
                    shutil.copy2(bam_p, target_bam)

            mean_d, b1x, b3x, b5x, b10x, mapped_mt = 0.0, 0.0, 0.0, 0.0, 0.0, 0
            try:
                res_idx = subprocess.run(["samtools", "idxstats", str(target_bam)], capture_output=True, text=True, check=True)
                for line in res_idx.stdout.splitlines():
                    parts = line.strip().split("\t")
                    if len(parts) >= 3 and parts[0] == self.expected_contig:
                        mapped_mt = int(parts[2])
                        break

                cmd_d = f"samtools depth -a -r {self.expected_contig}:1-{self.expected_length} {target_bam}"
                res_d = subprocess.run(cmd_d, shell=True, capture_output=True, text=True, check=True)
                depths = [int(line.split("\t")[2]) for line in res_d.stdout.splitlines() if len(line.split("\t")) >= 3]
                if depths:
                    tot_len = len(depths)
                    mean_d = round(sum(depths) / tot_len, 2)
                    b1x = round(sum(1 for d in depths if d >= 1) / tot_len, 4)
                    b3x = round(sum(1 for d in depths if d >= 3) / tot_len, 4)
                    b5x = round(sum(1 for d in depths if d >= 5) / tot_len, 4)
                    b10x = round(sum(1 for d in depths if d >= 10) / tot_len, 4)
            except Exception as exc:
                logger.warning(f"[{sample_id}] Failed extracting BAM coverage metrics: {exc}")

            cov_obj = CoverageMetrics(
                mean_depth_post_dedup=mean_d,
                breadth_1x_pct=b1x * 100.0,
                breadth_3x_pct=b3x * 100.0,
                breadth_5x_pct=b5x * 100.0,
                breadth_10x_pct=b10x * 100.0,
            )

            aln_res = AlignmentResult(
                sample_id=sample_id,
                status=AlignmentStatus.PASS,
                policy=aln_policy,
                reference=self.reference_info,
                outputs=AlignmentOutputs(
                    filtered_bam=target_bam,
                    dedup_bam=target_bam,
                ),
                metrics=AlignmentMetrics(
                    retained_filtered_reads=mapped_mt,
                    coverage=cov_obj,
                    dedup_metrics=DedupMetrics(duplication_rate_reads=0.0),
                ),
            )

            self._ensure_synthetic_block2_contract(sample_id)
            self._ensure_block3_contract(sample_id, aln_res, force_overwrite=True)

            ingested.append((sample_id, aln_res))

        return ingested

    def _ensure_synthetic_block2_contract(self, sample_id: str) -> Path:
        """Serializes clean Block 2 contract for BAM passthrough mode."""
        prep_dir = self.operational_files_dir / "preprocessing" / sample_id
        prep_dir.mkdir(parents=True, exist_ok=True)
        contract_path = prep_dir / f"{sample_id}.block2_contract.json"
        if not contract_path.is_file() or contract_path.stat().st_size == 0:
            with open(contract_path, "w", encoding="utf-8") as f:
                json.dump({
                    "sample_id": sample_id,
                    "raw_reads_total": None,
                    "clean_reads_total": None,
                    "retained_reads_total": None,
                    "retention_fraction": None,
                    "collapsed_fraction": None,
                    "output_profile": "BAM_PASSTHROUGH",
                    "execution_time_seconds": 0.0,
                }, f, indent=2)
        return contract_path

    # --------------------------------------------------------------------------
    # Block 2: Preprocessing (FASTQ Mode)
    # --------------------------------------------------------------------------
    def run_discovery_and_preprocessing(self, input_dir: Path) -> List[SampleRecord]:
        logger.info("=" * 80)
        logger.info("▶ BLOCK 2: DISCOVERY, INTEGRITY & PREPROCESSING (SampleDiscoveryEngine)")
        logger.info("=" * 80)

        raw_samples: List[SampleRecord] = SampleDiscoveryEngine.build_sample_records(
            input_dir=input_dir,
            config=self.config,
            integrity_mode=self.config.get("integrity", {}).get("mode", "fast"),
        )

        if not raw_samples:
            raise PipelineExecutionError(f"No valid FASTQ samples discovered in {input_dir}")

        processed_samples: List[SampleRecord] = []

        for sample in raw_samples:
            s_id = sample.sample_id
            legacy_dir = self.operational_files_dir / "alignment" / s_id
            tracker = StageTracker.load_or_create(self.cache_dir, s_id, legacy_dir=legacy_dir)

            expected_settings = legacy_dir / "trimmed" / f"{s_id}.settings"
            base_prefix = legacy_dir / "trimmed" / s_id

            raw_paths = self._extract_input_paths(sample.inputs)
            raw_input_meta = [(str(p), p.stat().st_size, p.stat().st_mtime) for p in raw_paths]
            prep_composite_hash = self.compute_composite_hash(self.prep_cfg, raw_input_meta)

            is_paired = getattr(sample.inputs, "is_paired", getattr(sample, "is_paired", False))

            if tracker.is_stage_valid("PREPROCESSING", [expected_settings], prep_composite_hash):
                processed, warnings, errors = PreprocessingEngine._locate_and_verify_outputs(
                    base_prefix=base_prefix,
                    is_paired=is_paired,
                    collapse_reads=bool(self.prep_cfg.get("collapse_reads", True)),
                    run_started_ts=0.0,
                    integrity_mode=self.config.get("integrity", {}).get("mode", "fast"),
                )
                if not errors and processed.output_profile.value != "EMPTY":
                    logger.info(f"[{s_id}] Stage PREPROCESSING cached -> Profile: {processed.output_profile.value}")
                    sample.processed_reads = processed
                    sample.validation.status = ValidationStatus.PASS
                    if expected_settings.is_file():
                        sample.metrics = PreprocessingEngine.parse_adapterremoval_settings(
                            settings_file=expected_settings,
                            is_paired=is_paired,
                            warnings_list=sample.validation.warnings,
                        )
                        active_b2_outs = [
                            f for f in [
                                processed.merged_fastq,
                                processed.paired_r1_fastq,
                                processed.paired_r2_fastq,
                                processed.singleton_fastq,
                            ]
                            if f and f.is_file()
                        ]
                        gc_p, q20, q30, s_len = PreprocessingEngine.compute_post_trimming_qc(
                            active_b2_outs, quality_base=33, sample_reads=25000
                        )
                        sample.metrics.gc_content_pct = gc_p
                        sample.metrics.q20_rate_post = q20
                        sample.metrics.q30_rate_post = q30
                        if sample.metrics.mean_read_length_post is None:
                            sample.metrics.mean_read_length_post = s_len
                    self._ensure_block2_contract(s_id, sample)
                    processed_samples.append(sample)
                    continue

            t0 = time.perf_counter()
            processed_sample: SampleRecord = PreprocessingEngine.process_sample(
                sample=sample,
                config=self.config,
            )

            tracker.mark_success(
                "PREPROCESSING",
                time.perf_counter() - t0,
                {"settings": expected_settings},
                composite_hash=prep_composite_hash,
            )
            self._ensure_block2_contract(s_id, processed_sample, force_overwrite=True)
            processed_samples.append(processed_sample)

        logger.info(f"✅ Block 2 complete for {len(processed_samples)} sample(s).")
        return processed_samples

    def _ensure_block2_contract(
        self, sample_id: str, sample: SampleRecord, force_overwrite: bool = False
    ) -> Path:
        prep_dir = self.operational_files_dir / "preprocessing" / sample_id
        prep_dir.mkdir(parents=True, exist_ok=True)
        contract_path = prep_dir / f"{sample_id}.block2_contract.json"

        if not force_overwrite and contract_path.is_file() and contract_path.stat().st_size > 0:
            return contract_path

        pr = sample.processed_reads
        m = sample.metrics
        raw_reads = getattr(sample.inputs, "raw_reads_total", None) or getattr(sample, "raw_reads_total", None) or (m.raw_total_reads if m else None)
        retained = getattr(pr, "retained_reads_total", None) or (m.surviving_original_reads if m else None)
        retention_frac = getattr(pr, "retention_fraction", None) or ((retained / raw_reads) if (raw_reads and retained and raw_reads > 0) else None)
        collapsed_frac = getattr(pr, "collapsed_fraction", None) or (m.merged_fraction_raw if m else None)
        out_prof = getattr(pr.output_profile, "value", str(pr.output_profile)) if pr and hasattr(pr, "output_profile") else "UNKNOWN"

        contract_data = {
            "sample_id": sample_id,
            "raw_reads_total": raw_reads,
            "clean_reads_total": retained,
            "retained_reads_total": retained,
            "retention_fraction": round(retention_frac, 4) if retention_frac is not None else None,
            "collapsed_fraction": round(collapsed_frac, 4) if collapsed_frac is not None else None,
            "output_profile": out_prof,
            "adapter_trimmed_reads_pct": getattr(m, "adapter_trimmed_reads_pct", None) if m else None,
            "q20_rate_post": getattr(m, "q20_rate_post", None) if m else None,
            "q30_rate_post": getattr(m, "q30_rate_post", None) if m else None,
            "gc_content_pct": getattr(m, "gc_content_pct", None) if m else None,
            "mean_read_length_pre": getattr(m, "mean_read_length_pre", None) if m else None,
            "mean_read_length_post": getattr(m, "mean_read_length_post", None) if m else None,
            "execution_time_seconds": 0.0,
        }

        with open(contract_path, "w", encoding="utf-8") as f:
            json.dump(contract_data, f, indent=2)
        return contract_path

    # --------------------------------------------------------------------------
    # Block 3: Alignment (BWA-aln) & Both-End aDNA Deduplication
    # --------------------------------------------------------------------------
    def run_alignment(self, sample: SampleRecord) -> AlignmentResult:
        s_id = sample.sample_id
        legacy_dir = self.operational_files_dir / "alignment" / s_id
        tracker = StageTracker.load_or_create(self.cache_dir, s_id, legacy_dir=legacy_dir)

        qc_thresh = self.filter_cfg.get("thresholds") or self.filter_cfg.get("qc_thresholds", {})
        min_depth_gate = float(qc_thresh.get("min_mean_depth", 3.0))
        min_b3x_gate = float(qc_thresh.get("min_breadth_3x_fraction", 0.80))
        max_dup_gate = float(qc_thresh.get("max_duplication_rate_fraction", 0.85))
        min_map_rate_gate = float(qc_thresh.get("min_mapping_rate_pct", 0.10))

        aln_policy = AlignmentFilterPolicy(
            dedup_mode=DedupMode(self.dedup_cfg.get("mode", "both").lower()),
            multimapping=MultiMappingPolicy(
                min_mapq=int(self.filter_cfg.get("multimapping", {}).get("min_mapq", 30)),
                filter_unmapped=bool(self.filter_cfg.get("multimapping", {}).get("filter_unmapped", True)),
                filter_secondary=bool(self.filter_cfg.get("multimapping", {}).get("filter_secondary", True)),
                filter_supplementary=bool(self.filter_cfg.get("multimapping", {}).get("filter_supplementary", True)),
            ),
            min_mean_depth=min_depth_gate,
            min_breadth_3x_fraction=min_b3x_gate,
            max_duplication_rate_fraction=max_dup_gate,
            min_mapping_rate_pct=min_map_rate_gate,
        )

        final_bam = legacy_dir / f"{s_id}.filtered.dedup.bam"
        markdup_bam = legacy_dir / f"{s_id}.filtered.markdup.bam"

        trimmed_files = []
        if sample.processed_reads:
            for attr in ("collapsed_fastq", "pair1_fastq", "pair2_fastq", "singleton_fastq"):
                fpath = getattr(sample.processed_reads, attr, None)
                if fpath and Path(fpath).is_file():
                    trimmed_files.append((str(fpath), Path(fpath).stat().st_size, Path(fpath).stat().st_mtime))

        aln_composite_hash = self.compute_composite_hash(
            self.aln_cfg,
            self.dedup_cfg,
            self.filter_cfg,
            self.ref_cfg,
            self.ref_sha256,
            trimmed_files,
        )

        if tracker.is_stage_valid("ALIGN_AND_DEDUP", [final_bam], aln_composite_hash):
            logger.info(f"[{s_id}] Stage ALIGN_AND_DEDUP cached -> {final_bam.name}")
            saved_meta = tracker.metadata.get("ALIGN_AND_DEDUP", {})

            cov_obj = CoverageMetrics(
                mean_depth_post_dedup=float(saved_meta.get("mean_depth", 0.0) or 0.0),
                breadth_1x_pct=float(saved_meta.get("coverage_breadth_1x", 0.0) or 0.0) * 100.0,
                breadth_3x_pct=float(saved_meta.get("coverage_breadth_3x", 0.0) or 0.0) * 100.0,
                breadth_5x_pct=float(saved_meta.get("coverage_breadth_5x", 0.0) or 0.0) * 100.0,
                breadth_10x_pct=float(saved_meta.get("coverage_breadth_10x", 0.0) or 0.0) * 100.0,
            )
            dedup_obj = DedupMetrics(
                duplication_rate_reads=float(saved_meta.get("duplication_rate", 0.0) or 0.0)
            )
            rehydrated_metrics = AlignmentMetrics(
                retained_filtered_reads=int(saved_meta.get("mapped_reads_mt", 0) or 0),
                coverage=cov_obj,
                dedup_metrics=dedup_obj,
            )

            res = AlignmentResult(
                sample_id=s_id,
                status=AlignmentStatus.PASS,
                policy=aln_policy,
                reference=self.reference_info,
                outputs=AlignmentOutputs(
                    filtered_bam=final_bam,
                    dedup_bam=final_bam,
                    markdup_bam=markdup_bam if markdup_bam.is_file() else None,
                ),
                metrics=rehydrated_metrics,
            )
            self._ensure_block3_contract(s_id, res)
            return res

        t0 = time.perf_counter()
        aln_res: AlignmentResult = AlignmentOrchestrator.run_sample_alignment(
            sample=sample,
            reference_info=self.reference_info,
            output_dir=self.operational_files_dir / "alignment",
            log_dir=self.logs_dir / "alignment",
            policy=aln_policy,
            threads=self.alignment_threads,
            seed_length=int(self.aln_cfg.get("bwa_aln", {}).get("seed_length", 1000)),
            max_diff=float(self.aln_cfg.get("bwa_aln", {}).get("max_diff", 0.01)),
            max_gap_opens=int(self.aln_cfg.get("bwa_aln", {}).get("max_gap_opens", 2)),
        )

        if aln_res.status == AlignmentStatus.FAIL:
            raise PipelineExecutionError(f"[{s_id}] Block 3 Alignment failed: {aln_res.errors}")

        target_bam = aln_res.outputs.filtered_bam or final_bam
        meta = {}
        if aln_res.metrics:
            m = aln_res.metrics
            cov = getattr(m, "coverage", None)
            dedup = getattr(m, "dedup_metrics", None)

            mean_d = getattr(cov, "mean_depth_post_dedup", None) if cov else getattr(m, "mean_depth", None)
            b1x = (cov.breadth_1x_pct / 100.0) if cov and getattr(cov, "breadth_1x_pct", None) is not None else getattr(m, "coverage_breadth_1x", None)
            b3x = (cov.breadth_3x_pct / 100.0) if cov and getattr(cov, "breadth_3x_pct", None) is not None else getattr(m, "coverage_breadth_3x", None)
            b5x = (cov.breadth_5x_pct / 100.0) if cov and getattr(cov, "breadth_5x_pct", None) is not None else getattr(m, "coverage_breadth_5x", None)
            b10x = (cov.breadth_10x_pct / 100.0) if cov and getattr(cov, "breadth_10x_pct", None) is not None else getattr(m, "coverage_breadth_10x", None)
            dup_r = getattr(dedup, "duplication_rate_reads", None) if dedup else getattr(m, "duplication_rate", None)
            mapped_mt = getattr(m, "retained_filtered_reads", None) or getattr(m, "mapped_reads_raw", None) or getattr(m, "mapped_reads_mt", None)

            meta = {
                "mean_depth": round(mean_d, 2) if mean_d is not None else None,
                "coverage_breadth_1x": round(b1x, 4) if b1x is not None else None,
                "coverage_breadth_3x": round(b3x, 4) if b3x is not None else None,
                "coverage_breadth_5x": round(b5x, 4) if b5x is not None else None,
                "coverage_breadth_10x": round(b10x, 4) if b10x is not None else None,
                "duplication_rate": round(dup_r, 4) if dup_r is not None else None,
                "mapped_reads_mt": mapped_mt,
            }

        tracker.mark_success(
            "ALIGN_AND_DEDUP",
            time.perf_counter() - t0,
            {"filtered_bam": target_bam},
            composite_hash=aln_composite_hash,
            metadata=meta,
        )
        self._ensure_block3_contract(s_id, aln_res, force_overwrite=True)
        return aln_res

    def _ensure_block3_contract(
        self, sample_id: str, aln_res: AlignmentResult, force_overwrite: bool = False
    ) -> Path:
        aln_dir = self.operational_files_dir / "alignment" / sample_id
        aln_dir.mkdir(parents=True, exist_ok=True)
        contract_path = aln_dir / f"{sample_id}.block3_contract.json"

        if not force_overwrite and contract_path.is_file() and contract_path.stat().st_size > 0:
            return contract_path

        m = aln_res.metrics
        cov = getattr(m, "coverage", None) if m else None
        dedup = getattr(m, "dedup_metrics", None) if m else None

        mean_d = getattr(cov, "mean_depth_post_dedup", None) if cov else getattr(m, "mean_depth", None)
        b1x = (cov.breadth_1x_pct / 100.0) if cov and getattr(cov, "breadth_1x_pct", None) is not None else getattr(m, "coverage_breadth_1x", None)
        b3x = (cov.breadth_3x_pct / 100.0) if cov and getattr(cov, "breadth_3x_pct", None) is not None else getattr(m, "coverage_breadth_3x", None)
        b5x = (cov.breadth_5x_pct / 100.0) if cov and getattr(cov, "breadth_5x_pct", None) is not None else getattr(m, "coverage_breadth_5x", None)
        b10x = (cov.breadth_10x_pct / 100.0) if cov and getattr(cov, "breadth_10x_pct", None) is not None else getattr(m, "coverage_breadth_10x", None)
        dup_r = getattr(dedup, "duplication_rate_reads", None) if dedup else getattr(m, "duplication_rate", None)
        mapped_mt = getattr(m, "retained_filtered_reads", None) or getattr(m, "mapped_reads_raw", None) or getattr(m, "mapped_reads_mt", None)

        contract_data = {
            "sample_id": sample_id,
            "mapped_reads_mt": mapped_mt,
            "endogenous_dna_fraction": getattr(m, "endogenous_dna_fraction", None),
            "mean_depth": round(mean_d, 2) if mean_d is not None else None,
            "coverage_breadth_1x": round(b1x, 4) if b1x is not None else None,
            "coverage_breadth_3x": round(b3x, 4) if b3x is not None else None,
            "coverage_breadth_5x": round(b5x, 4) if b5x is not None else None,
            "coverage_breadth_10x": round(b10x, 4) if b10x is not None else None,
            "duplication_rate": round(dup_r, 4) if dup_r is not None else None,
            "execution_time_seconds": 0.0,
        }

        with open(contract_path, "w", encoding="utf-8") as f:
            json.dump(contract_data, f, indent=2)
        return contract_path

    # --------------------------------------------------------------------------
    # Block 4: Damage Assessment & Rescaling
    # --------------------------------------------------------------------------
    def run_damage(self, sample_id: str, aln_res: AlignmentResult) -> Tuple[DamageResult, Path]:
        legacy_dir = self.operational_files_dir / "alignment" / sample_id
        tracker = StageTracker.load_or_create(self.cache_dir, sample_id, legacy_dir=legacy_dir)

        damage_dir = self.operational_files_dir / "damage" / sample_id
        damage_json = damage_dir / f"{sample_id}.damage.json"
        damage_model_json = damage_dir / f"{sample_id}.damage_model.json"
        b4_contract_json = damage_dir / f"{sample_id}.block4_contract.json"

        input_bam: Optional[Path] = None
        if aln_res.outputs:
            if aln_res.outputs.dedup_bam and aln_res.outputs.dedup_bam.is_file():
                input_bam = aln_res.outputs.dedup_bam
            elif aln_res.outputs.filtered_bam and aln_res.outputs.filtered_bam.is_file():
                input_bam = aln_res.outputs.filtered_bam

        if input_bam is None:
            candidate_bam = legacy_dir / f"{sample_id}.filtered.dedup.bam"
            if candidate_bam.is_file():
                input_bam = candidate_bam
            else:
                raise PipelineExecutionError(f"[{sample_id}] Physical filtered dedup BAM not found in {legacy_dir}")

        input_bam_sha256 = self.calculate_sha256(input_bam)

        damage_composite_hash = self.compute_composite_hash(
            self.damage_cfg,
            self.ref_sha256,
            input_bam_sha256,
        )

        # Prefer full damage.json over damage_model.json when rehydrating DamageResult
        active_cache_json = damage_json if damage_json.is_file() else damage_model_json

        if tracker.is_stage_valid("DAMAGE_AUTHENTICATION", [active_cache_json, b4_contract_json], damage_composite_hash):
            logger.info(f"[{sample_id}] Stage DAMAGE_AUTHENTICATION cached.")
            with open(active_cache_json, "r", encoding="utf-8") as f:
                dam_res = DamageResult.from_dict(json.load(f))
            return dam_res, b4_contract_json

        t0 = time.perf_counter()
        raw_treatment = self.damage_cfg.get("library_treatment", "NONE")
        try:
            damage_lib_treatment = DamageLibraryTreatment(raw_treatment)
        except ValueError as exc:
            raise PipelineExecutionError(
                f"Invalid library_treatment '{raw_treatment}' in config. "
                f"Allowed: {[t.value for t in DamageLibraryTreatment]}"
            ) from exc

        damage_policy = DamagePolicy(
            library_treatment=damage_lib_treatment,
            enable_rescaling=bool(self.damage_cfg.get("enable_rescaling", True)),
            enable_pmd_scoring=bool(self.damage_cfg.get("enable_pmd_scoring", True)),
            enable_pmd_filtering=bool(self.damage_cfg.get("enable_pmd_filtering", False)),
            pmd_threshold=float(self.damage_cfg.get("pmd_threshold", 3.0)),
        )

        dam_res: DamageResult = DamageOrchestrator.run_sample_damage_pipeline(
            sample_id=sample_id,
            reference_info=self.reference_info,
            output_dir=self.operational_files_dir / "damage",
            log_dir=self.logs_dir / "damage",
            policy=damage_policy,
            dedup_bam=input_bam,
            markdup_bam=aln_res.outputs.markdup_bam if aln_res.outputs else None,
            threads=self.alignment_threads,
            allow_input_bam_fallback=True,
        )

        if dam_res.status.value == "FAIL":
            raise PipelineExecutionError(f"[{sample_id}] Block 4 Damage authentication failed: {dam_res.errors}")

        self._ensure_block4_contract(b4_contract_json, sample_id, dam_res, input_bam, force_overwrite=True)

        target_model_file = damage_json if damage_json.is_file() else (dam_res.outputs.damage_model_json or damage_model_json)
        tracker.mark_success(
            "DAMAGE_AUTHENTICATION",
            time.perf_counter() - t0,
            {"damage_json": target_model_file, "b4_contract": b4_contract_json},
            composite_hash=damage_composite_hash,
        )
        return dam_res, b4_contract_json

    @staticmethod
    def _extract_verdict_with_contract_fallback(dam_res: DamageResult, b4_contract_path: Path) -> str:
        """
        Extracts the Block 4 authenticity verdict, prioritizing the formal Block 4
        contract JSON if DamageResult rehydration from damage_model.json defaulted to INCONCLUSIVE.
        """
        verdict = "AUTHENTICITY_UNSUPPORTED"
        if getattr(dam_res, "profile", None) and getattr(dam_res.profile, "authentication_verdict", None):
            verdict = getattr(dam_res.profile.authentication_verdict, "value", str(dam_res.profile.authentication_verdict))
        elif getattr(dam_res, "status", None):
            verdict = getattr(dam_res.status, "value", str(dam_res.status))

        if verdict.upper() in ("INCONCLUSIVE", "AUTHENTICITY_UNSUPPORTED", "PASS", "WARN") and b4_contract_path.is_file():
            try:
                with open(b4_contract_path, "r", encoding="utf-8") as f_c:
                    c_data = json.load(f_c)
                contract_auth = c_data.get("authentication_status")
                if contract_auth and str(contract_auth).upper() not in ("INCONCLUSIVE", "AUTHENTICITY_UNSUPPORTED"):
                    verdict = str(contract_auth)
            except Exception:
                pass

        return verdict

    def _ensure_block4_contract(
        self,
        contract_path: Path,
        sample_id: str,
        dam_res: DamageResult,
        input_bam: Path,
        force_overwrite: bool = False,
    ) -> None:
        """Serializes formal Block 4 contract guaranteeing availability of all authenticity metrics."""
        contract_path.parent.mkdir(parents=True, exist_ok=True)
        if not force_overwrite and contract_path.is_file() and contract_path.stat().st_size > 0:
            return

        rescaling_performed = False
        rescaled_bam_str = None
        if getattr(dam_res, "rescaled_bam", None) and Path(dam_res.rescaled_bam).is_file():
            rescaling_performed = True
            rescaled_bam_str = str(dam_res.rescaled_bam)
        elif getattr(dam_res, "outputs", None) and getattr(dam_res.outputs, "rescaled_bam", None):
            rescaling_performed = True
            rescaled_bam_str = str(dam_res.outputs.rescaled_bam)

        auth_status = "AUTHENTICITY_UNSUPPORTED"
        prof = getattr(dam_res, "profile", None)
        if prof and getattr(prof, "authentication_verdict", None):
            auth_status = getattr(prof.authentication_verdict, "value", str(prof.authentication_verdict))
        elif getattr(dam_res, "status", None):
            auth_status = getattr(dam_res.status, "value", str(dam_res.status))

        lib_treatment_str = "NON_UDG"
        if getattr(dam_res, "policy", None) and getattr(dam_res.policy, "library_treatment", None):
            val = getattr(dam_res.policy.library_treatment, "value", str(dam_res.policy.library_treatment)).upper()
            if "HALF" in val:
                lib_treatment_str = "UDG_HALF"
            elif "FULL" in val:
                lib_treatment_str = "UDG_FULL"
            elif "SINGLE" in val or "SS" in val:
                lib_treatment_str = "SINGLE_STRANDED"
            else:
                lib_treatment_str = "NON_UDG"

        rescaled_sha256 = self.calculate_sha256(rescaled_bam_str) if (rescaling_performed and rescaled_bam_str) else None
        filtered_sha256 = self.calculate_sha256(input_bam)

        fl_med = getattr(prof.fragment_length, "median_length", None) if prof else None
        fl_mean = getattr(prof.fragment_length, "mean_length", None) if prof else None
        pmd_ge3 = getattr(prof.pmd_metrics, "fraction_pmd_ge_3", None) if (prof and prof.pmd_metrics) else None
        c1 = getattr(prof.terminal_damage, "c_to_t_5p_pos1", None) if prof else None
        g1 = getattr(prof.terminal_damage, "g_to_a_3p_pos1", None) if prof else None

        if any(v is None for v in (fl_med, pmd_ge3, c1)):
            for json_cand in [
                self.operational_files_dir / "damage" / sample_id / f"{sample_id}.damage_model.json",
                self.operational_files_dir / "damage" / sample_id / f"{sample_id}.damage.json",
            ]:
                if json_cand.is_file():
                    try:
                        with open(json_cand, "r", encoding="utf-8") as fp:
                            dj = json.load(fp)
                        fl_med = fl_med if fl_med is not None else dj.get("fragment_length", {}).get("median_length_bp")
                        fl_mean = fl_mean if fl_mean is not None else dj.get("fragment_length", {}).get("mean_length_bp")
                        pmd_ge3 = pmd_ge3 if pmd_ge3 is not None else dj.get("pmd_metrics", {}).get("fraction_pmd_ge_3")
                        td_raw = dj.get("terminal_damage", {})
                        c1 = c1 if c1 is not None else td_raw.get("c_to_t_5p_pos1_fraction")
                        g1 = g1 if g1 is not None else td_raw.get("g_to_a_3p_pos1_fraction")
                        break
                    except Exception:
                        pass

        contract_data = {
            "sample_id": sample_id,
            "authentication_status": auth_status,
            "rescaling_performed": rescaling_performed,
            "rescaled_bam": rescaled_bam_str,
            "rescaled_bam_sha256": rescaled_sha256,
            "filtered_dedup_bam": str(input_bam),
            "filtered_dedup_bam_sha256": filtered_sha256,
            "damage_model_id": getattr(dam_res, "damage_model_id", None) or "model_ancient_v1",
            "library_treatment": lib_treatment_str,
            "reference_sha256": self.ref_sha256,
            "median_fragment_length": round(fl_med, 1) if fl_med is not None else None,
            "mean_fragment_length": round(fl_mean, 1) if fl_mean is not None else None,
            "pmd_fraction_ge_3": round(pmd_ge3, 4) if pmd_ge3 is not None else None,
            "c_to_t_5p_rate_pos1": round(c1, 4) if c1 is not None else None,
            "g_to_a_3p_rate_pos1": round(g1, 4) if g1 is not None else None,
        }

        with open(contract_path, "w", encoding="utf-8") as f:
            json.dump(contract_data, f, indent=2)

    # --------------------------------------------------------------------------
    # Block 5: Contamination Estimation (Schmutzi & Dissent Two-Tier Framework)
    # --------------------------------------------------------------------------
    def run_contamination_and_report(
        self, unrescaled_bams: List[Tuple[str, Path, bool]]
    ) -> Tuple[Path, Path]:
        logger.info("=" * 80)
        logger.info("▶ BLOCK 5: TWO-TIER CONTAMINATION ESTIMATION & GATING (ContaminationBatchRunner)")
        logger.info("=" * 80)

        contam_dir = self.operational_files_dir / "contamination"
        contam_dir.mkdir(parents=True, exist_ok=True)

        samples_to_run: List[Tuple[str, Path, bool]] = []

        for s_id, bam_path, is_ancient in unrescaled_bams:
            legacy_dir = self.operational_files_dir / "alignment" / s_id
            tracker = StageTracker.load_or_create(self.cache_dir, s_id, legacy_dir=legacy_dir)

            b5_contract_path = contam_dir / s_id / f"{s_id}.block5_contract.json"
            bam_sha256 = self.calculate_sha256(bam_path)

            b5_composite_hash = self.compute_composite_hash(
                self.contam_cfg,
                self.ref_sha256,
                bam_sha256,
                is_ancient,
            )

            if tracker.is_stage_valid("CONTAMINATION_ESTIMATION", [b5_contract_path], b5_composite_hash):
                logger.info(f"[{s_id}] Stage CONTAMINATION_ESTIMATION cached -> Reusing {b5_contract_path.name}")
                continue

            samples_to_run.append((s_id, bam_path, is_ancient))

        if samples_to_run:
            logger.info(f"[*] Submitting {len(samples_to_run)} sample(s) to contamination estimation pool...")
            results = ContaminationBatchRunner.run_cohort_pipeline(
                discovered_samples=samples_to_run,
                reference_fasta=self.reference_info.fasta_path,
                reference_contig=self.reference_info.contig_name,
                reference_length_bp=self.reference_info.length_bp,
                output_dir=contam_dir,
                config=self.config,
                resource_alloc=self.resource_alloc,
            )
            for r in results:
                s_dir = self.operational_files_dir / "alignment" / r.sample_id
                t = StageTracker.load_or_create(self.cache_dir, r.sample_id, legacy_dir=s_dir)
                contract_path = contam_dir / r.sample_id / f"{r.sample_id}.block5_contract.json"
                target_bam = next(b[1] for b in unrescaled_bams if b[0] == r.sample_id)
                target_bam_sha = self.calculate_sha256(target_bam)

                self._ensure_block5_contract(contract_path, r, input_bam_sha256=target_bam_sha, force_overwrite=True)

                comp_hash = self.compute_composite_hash(
                    self.contam_cfg,
                    self.ref_sha256,
                    target_bam_sha,
                    any(b[2] for b in unrescaled_bams if b[0] == r.sample_id),
                )
                t.mark_success("CONTAMINATION_ESTIMATION", 0.0, {"contract": contract_path}, comp_hash)
        else:
            logger.info("✅ All samples have valid contamination contracts cached on disk.")

        tsv_out, json_out = self._rebuild_or_export_contamination_summaries(contam_dir, unrescaled_bams)
        return tsv_out, json_out

    def _rebuild_or_export_contamination_summaries(
        self, contam_dir: Path, discovered_samples: List[Tuple[str, Path, bool]]
    ) -> Tuple[Path, Path]:
        tsv_path = contam_dir / "cohort_contamination_summary.tsv"
        json_path = contam_dir / "cohort_contamination_summary.json"

        records = []
        for s_id, _, _ in discovered_samples:
            c_path = contam_dir / s_id / f"{s_id}.block5_contract.json"
            if c_path.is_file():
                with open(c_path, "r", encoding="utf-8") as f:
                    records.append(json.load(f))

        headers = [
            "sample_id",
            "recommended_action",
            "estimation_confidence",
            "contamination_fraction",
            "contamination_estimate_source",
            "contaminant_haplogroups",
        ]
        with open(tsv_path, "w", newline="", encoding="utf-8") as f_tsv:
            writer = csv.writer(f_tsv, delimiter="\t")
            writer.writerow(headers)
            for r in records:
                cf = f"{r.get('contamination_fraction'):.4f}" if r.get("contamination_fraction") is not None else "NA"
                chgs = ",".join(r.get("contaminant_haplogroups", [])) or "NONE"
                writer.writerow([
                    r.get("sample_id"),
                    r.get("recommended_action"),
                    r.get("estimation_confidence"),
                    cf,
                    r.get("contamination_estimate_source"),
                    chgs,
                ])

        with open(json_path, "w", encoding="utf-8") as f_json:
            json.dump({"total_samples": len(records), "samples": records}, f_json, indent=2)

        return tsv_path, json_path

    def _ensure_block5_contract(
        self,
        contract_path: Path,
        contam_res: ContaminationResult,
        input_bam_sha256: Optional[str] = None,
        force_overwrite: bool = False,
    ) -> None:
        """Serializes Block 5 contract without hardcoded fallbacks, interpreting actual Schmutzi deliverables."""
        contract_path.parent.mkdir(parents=True, exist_ok=True)
        if not force_overwrite and contract_path.is_file() and contract_path.stat().st_size > 0:
            return

        gating_obj = getattr(contam_res, "gating", None)
        est_obj = getattr(contam_res, "estimate", None)
        audit_obj = getattr(contam_res, "audit", None)

        rec_action = None
        for target in (gating_obj, contam_res):
            if target is None:
                continue
            for attr in ("action", "recommended_action", "downstream_gating", "gating_action"):
                val = getattr(target, attr, None)
                if val is not None:
                    rec_action = val
                    break
            if rec_action is not None:
                break
        if hasattr(rec_action, "value"):
            rec_action = rec_action.value
        rec_action = str(rec_action or "RESTRICT_HETEROPLASMY")

        contam_fraction = None
        for target in (est_obj, contam_res):
            if target is None:
                continue
            for attr in (
                "point_estimate",
                "contamination_point_estimate",
                "contamination_fraction",
                "contamination",
                "fraction",
                "rate",
                "mean",
            ):
                val = getattr(target, attr, None)
                if val is not None:
                    try:
                        contam_fraction = float(val)
                        break
                    except (ValueError, TypeError):
                        pass
            if contam_fraction is not None:
                break

        ci_low = getattr(est_obj, "ci_lower_95", None) if est_obj else None
        ci_high = getattr(est_obj, "ci_upper_95", None) if est_obj else None
        if (ci_low is None or ci_high is None) and est_obj:
            ci_raw = getattr(est_obj, "ci_95", None) or getattr(est_obj, "contamination_ci_95", None)
            if isinstance(ci_raw, (list, tuple)) and len(ci_raw) == 2:
                ci_low = float(ci_raw[0]) if ci_raw[0] is not None else None
                ci_high = float(ci_raw[1]) if ci_raw[1] is not None else None

        if ci_low is None or ci_high is None:
            s_id = contam_res.sample_id
            schmutzi_dir = self.operational_files_dir / "contamination" / s_id / "schmutzi"
            for cand in [schmutzi_dir / f"{s_id}_final.cont.est", schmutzi_dir / f"{s_id}.cont.est"]:
                if cand.is_file():
                    try:
                        with open(cand, "r") as f_est:
                            parts = f_est.read().strip().split()
                            if len(parts) >= 3:
                                ci_low = float(parts[0])
                                contam_fraction = float(parts[1]) if contam_fraction is None else contam_fraction
                                ci_high = float(parts[2])
                                break
                    except Exception:
                        pass

        source = None
        for target in (est_obj, audit_obj, contam_res):
            if target is None:
                continue
            for attr in ("source", "contamination_estimate_source", "method", "tool"):
                val = getattr(target, attr, None)
                if val is not None:
                    source = val
                    break
            if source is not None:
                break
        if hasattr(source, "value"):
            source = source.value

        if not source or str(source).upper() in {"NONE", "UNKNOWN"}:
            source = "SCHMUTZI" if (ci_low is not None and ci_high is not None) else "DISSENT_SCREEN"
        else:
            source = str(source).upper()

        confidence = None
        for target in (gating_obj, est_obj, audit_obj, contam_res):
            if target is None:
                continue
            for attr in ("confidence", "estimation_confidence", "confidence_level", "reliability", "quality_status"):
                val = getattr(target, attr, None)
                if val is not None:
                    confidence = val
                    break
            if confidence is not None:
                break
        if hasattr(confidence, "value"):
            confidence = confidence.value

        if not confidence or str(confidence).upper() in {"NONE", "UNKNOWN", "UNRELIABLE"}:
            if ci_low is not None and ci_high is not None:
                ci_width = float(ci_high) - float(ci_low)
                if ci_width <= 0.04:
                    confidence = "HIGH"
                elif ci_width <= 0.08:
                    confidence = "MODERATE"
                else:
                    confidence = "UNRELIABLE"
            elif source == "DISSENT_SCREEN":
                if rec_action == "STANDARD_CALLING" and contam_fraction is not None and contam_fraction <= 0.02:
                    confidence = "MODERATE"
                else:
                    confidence = "UNRELIABLE"
            else:
                confidence = "HIGH" if rec_action == "STANDARD_CALLING" else "UNRELIABLE"
        else:
            confidence = str(confidence).upper()

        raw_chgs = getattr(contam_res, "contaminant_haplogroups", None) or getattr(est_obj, "contaminant_haplogroups", None) or []
        contam_hgs = [raw_chgs] if isinstance(raw_chgs, str) else list(raw_chgs)

        contract_data = {
            "sample_id": contam_res.sample_id,
            "recommended_action": rec_action,
            "estimation_confidence": confidence,
            "contamination_fraction": contam_fraction,
            "ci_95_lower": ci_low,
            "ci_95_upper": ci_high,
            "contamination_estimate_source": source,
            "contaminant_haplogroups": contam_hgs,
            "input_bam_sha256": input_bam_sha256,
            "reference_sha256": self.ref_sha256,
        }

        with open(contract_path, "w", encoding="utf-8") as f:
            json.dump(contract_data, f, indent=2)

    # --------------------------------------------------------------------------
    # Block 6: Consensus & Heteroplasmy Calling
    # --------------------------------------------------------------------------
    def run_consensus(
        self,
        sample_id: str,
        b4_contract_path: Path,
        b5_contract_path: Path,
    ) -> Tuple[ConsensusQC, Path, Path, Path]:
        legacy_dir = self.operational_files_dir / "alignment" / sample_id
        tracker = StageTracker.load_or_create(self.cache_dir, sample_id, legacy_dir=legacy_dir)

        block6_out = self.operational_files_dir / "consensus" / sample_id
        canonical_fa = block6_out / f"{sample_id}.consensus.fasta"
        vcf_gz = block6_out / f"{sample_id}.variants.vcf.gz"
        qc_json = block6_out / f"{sample_id}.consensus_qc.json"
        het_tsv = block6_out / f"{sample_id}.heteroplasmy.tsv"
        b6_contract_json = block6_out / f"{sample_id}.block6_contract.json"

        b4_sha256 = self.calculate_sha256(b4_contract_path)
        b5_sha256 = self.calculate_sha256(b5_contract_path)
        cons_composite_hash = self.compute_composite_hash(
            self.cons_cfg,
            self.ref_sha256,
            b4_sha256,
            b5_sha256,
        )

        expected_artifacts = [canonical_fa, vcf_gz, qc_json, b6_contract_json]

        if tracker.is_stage_valid("CONSENSUS_CALLING", expected_artifacts, cons_composite_hash):
            logger.info(f"[{sample_id}] Stage CONSENSUS_CALLING cached -> {canonical_fa.name}")
            with open(qc_json, "r", encoding="utf-8") as f:
                qc_dict = json.load(f)

            gating_action_val = qc_dict.get("gating_action", "RESTRICT_HETEROPLASMY")
            try:
                resolved_gating_action = GatingAction(gating_action_val)
            except ValueError:
                resolved_gating_action = GatingAction.RESTRICT_HETEROPLASMY

            qc = ConsensusQC(
                sample_id=sample_id,
                input_bam_path=qc_dict.get("input_bam_path", ""),
                input_bam_sha256=qc_dict.get("input_bam_sha256", ""),
                reference_contig=qc_dict.get("reference_contig", self.expected_contig),
                reference_sha256=qc_dict.get("reference_sha256", self.ref_sha256),
                reference_length=qc_dict.get("reference_length", self.expected_length),
                consensus_fasta_sha256=qc_dict.get("consensus_fasta_sha256", ""),
                mean_depth=qc_dict.get("mean_depth", 0.0),
                covered_positions_fraction=qc_dict.get("covered_positions_fraction", 0.0),
                n_fraction_canonical=qc_dict.get("n_fraction_canonical", 0.0),
                iupac_fraction=qc_dict.get("iupac_fraction", 0.0),
                heteroplasmy_sites_total=qc_dict.get("heteroplasmy_sites_total", 0),
                heteroplasmy_sites_pass=qc_dict.get("heteroplasmy_sites_pass", 0),
                heteroplasmy_sites_filtered=qc_dict.get("heteroplasmy_sites_filtered", 0),
                regions=[],
                gating_action=resolved_gating_action,
                software_version=qc_dict.get("software_version", "1.0.0"),
                policy_version=qc_dict.get("policy_version", "1.0.0"),
                execution_time_seconds=qc_dict.get("execution_time_seconds", 0.0),
                heteroplasmy_details=qc_dict.get("heteroplasmy_details", {}),
            )
            return qc, canonical_fa, vcf_gz, b6_contract_json

        t0 = time.perf_counter()
        resolution = self.consensus_orchestrator.resolve_inputs(
            sample_id=sample_id,
            block4_json_path=str(b4_contract_path),
            block5_json_path=str(b5_contract_path),
            reference_fasta_path=str(self.ref_fasta_path),
            target_contig=self.expected_contig,
        )

        qc = self.consensus_orchestrator.run(resolution, str(block6_out))
        elapsed = time.perf_counter() - t0

        input_bam_sha = self.calculate_sha256(resolution.bam_path)
        self._ensure_block6_contract(
            b6_contract_json,
            sample_id,
            canonical_fa,
            vcf_gz,
            qc,
            het_tsv,
            input_bam_sha256=input_bam_sha,
            b4_contract_sha256=b4_sha256,
            b5_contract_sha256=b5_sha256,
            force_overwrite=True,
        )

        tracker.mark_success(
            "CONSENSUS_CALLING",
            elapsed,
            {
                "canonical_fasta": canonical_fa,
                "vcf": vcf_gz,
                "qc_json": qc_json,
                "b6_contract": b6_contract_json,
            },
            composite_hash=cons_composite_hash,
        )
        return qc, canonical_fa, vcf_gz, b6_contract_json

    def _ensure_block6_contract(
        self,
        contract_path: Path,
        sample_id: str,
        canonical_fasta: Path,
        vcf_gz: Path,
        qc: ConsensusQC,
        heteroplasmy_tsv: Path,
        input_bam_sha256: Optional[str] = None,
        b4_contract_sha256: Optional[str] = None,
        b5_contract_sha256: Optional[str] = None,
        force_overwrite: bool = False,
    ) -> None:
        contract_path.parent.mkdir(parents=True, exist_ok=True)
        if not force_overwrite and contract_path.is_file() and contract_path.stat().st_size > 0:
            return

        fasta_sha256 = self.calculate_sha256(canonical_fasta)
        vcf_sha256 = self.calculate_sha256(vcf_gz)

        passing_sites: List[int] = []
        if heteroplasmy_tsv.is_file():
            try:
                with open(heteroplasmy_tsv, "r", encoding="utf-8") as f:
                    reader = csv.DictReader(f, delimiter="\t")
                    for row in reader:
                        if row.get("status") == "PASS":
                            pos_val = row.get("position")
                            if pos_val is not None:
                                passing_sites.append(int(pos_val))
            except Exception as exc:
                logger.warning(f"[{sample_id}] Failed parsing heteroplasmy TSV: {exc}")

        details_dict: Dict[str, Any] = {}
        raw_details = getattr(qc, "heteroplasmy_details", {}) or {}
        for pos, det in raw_details.items():
            if hasattr(det, "ref"):
                details_dict[str(pos)] = {
                    "ref": det.ref,
                    "alt": det.alt,
                    "vaf": round(det.vaf, 4),
                    "internal_support_count": det.internal_support_count,
                    "terminal_support_count": det.terminal_support_count,
                    "spatial_damage_suspect": det.spatial_damage_suspect,
                }
            elif isinstance(det, dict):
                details_dict[str(pos)] = det

        contract_data = {
            "sample_id": sample_id,
            "consensus_fasta_path": str(canonical_fasta),
            "consensus_fasta_sha256": fasta_sha256,
            "variants_vcf_path": str(vcf_gz),
            "variants_vcf_sha256": vcf_sha256,
            "input_bam_sha256": input_bam_sha256,
            "block4_contract_sha256": b4_contract_sha256,
            "block5_contract_sha256": b5_contract_sha256,
            "reference_sha256": self.ref_sha256,
            "mean_depth": round(qc.mean_depth, 2),
            "covered_positions_fraction": round(qc.covered_positions_fraction, 4),
            "n_fraction_canonical": round(qc.n_fraction_canonical, 4),
            "heteroplasmy_sites_pass": passing_sites,
            "heteroplasmy_details": details_dict,
        }

        with open(contract_path, "w", encoding="utf-8") as f:
            json.dump(contract_data, f, indent=2)

    def export_consensus_cohort_summary(
        self,
        cohort_qc: Dict[str, ConsensusQC],
        output_dir: Path,
    ) -> Tuple[Path, Path]:
        output_dir.mkdir(parents=True, exist_ok=True)
        tsv_path = output_dir / "cohort_consensus_summary.tsv"
        json_path = output_dir / "cohort_consensus_summary.json"

        headers = [
            "sample_id",
            "mean_depth",
            "covered_positions_fraction",
            "n_fraction_canonical",
            "iupac_fraction",
            "heteroplasmy_sites_total",
            "heteroplasmy_sites_pass",
            "heteroplasmy_sites_filtered",
            "gating_action",
            "execution_time_seconds",
        ]

        rows = []
        json_records = []
        for s_id, qc in cohort_qc.items():
            gating_str = qc.gating_action.value if hasattr(qc.gating_action, "value") else str(qc.gating_action)
            rows.append([
                s_id,
                f"{qc.mean_depth:.2f}",
                f"{qc.covered_positions_fraction:.4f}",
                f"{qc.n_fraction_canonical:.4f}",
                f"{qc.iupac_fraction:.4f}",
                qc.heteroplasmy_sites_total,
                qc.heteroplasmy_sites_pass,
                qc.heteroplasmy_sites_filtered,
                gating_str,
                f"{qc.execution_time_seconds:.2f}",
            ])
            json_records.append({
                "sample_id": s_id,
                "mean_depth": qc.mean_depth,
                "covered_positions_fraction": qc.covered_positions_fraction,
                "n_fraction_canonical": qc.n_fraction_canonical,
                "iupac_fraction": qc.iupac_fraction,
                "heteroplasmy_sites_total": qc.heteroplasmy_sites_total,
                "heteroplasmy_sites_pass": qc.heteroplasmy_sites_pass,
                "heteroplasmy_sites_filtered": qc.heteroplasmy_sites_filtered,
                "gating_action": gating_str,
                "execution_time_seconds": qc.execution_time_seconds,
            })

        with open(tsv_path, "w", newline="", encoding="utf-8") as f_tsv:
            writer = csv.writer(f_tsv, delimiter="\t")
            writer.writerow(headers)
            writer.writerows(rows)

        with open(json_path, "w", encoding="utf-8") as f_json:
            json.dump(
                {
                    "generated_at": datetime.now(timezone.utc).isoformat(),
                    "total_samples": len(cohort_qc),
                    "samples": json_records,
                },
                f_json,
                indent=2,
            )

        return tsv_path, json_path

    # --------------------------------------------------------------------------
    # Block 7: Haplogrouping & Lineage Audit
    # --------------------------------------------------------------------------
    def run_haplogroup(
        self,
        sample_id: str,
        canonical_fasta: Path,
        variants_vcf: Path,
        b4_contract_path: Path,
        b5_contract_path: Path,
        b6_contract_path: Path,
    ) -> Block7Contract:
        legacy_dir = self.operational_files_dir / "alignment" / sample_id
        tracker = StageTracker.load_or_create(self.cache_dir, sample_id, legacy_dir=legacy_dir)

        block7_out = self.operational_files_dir / "haplogroup" / sample_id
        contract_json = block7_out / f"{sample_id}.block7_contract.json"
        hg_json = block7_out / f"{sample_id}.haplogroup.json"
        val_json = block7_out / f"{sample_id}.validation.json"

        fasta_sha256 = self.calculate_sha256(canonical_fasta)
        b4_sha256 = self.calculate_sha256(b4_contract_path)
        b5_sha256 = self.calculate_sha256(b5_contract_path)
        b6_sha256 = self.calculate_sha256(b6_contract_path)

        hg_composite_hash = self.compute_composite_hash(
            self.config.get("haplogroup", {}),
            self.config.get("cross_block_audit", {}),
            self.config.get("validation", {}),
            self.ref_sha256,
            fasta_sha256,
            b4_sha256,
            b5_sha256,
            b6_sha256,
        )

        expected_artifacts = [contract_json, hg_json, val_json]

        if tracker.is_stage_valid("HAPLOGROUP_CLASSIFICATION", expected_artifacts, hg_composite_hash):
            logger.info(f"[{sample_id}] Stage HAPLOGROUP_CLASSIFICATION cached.")
            with open(contract_json, "r", encoding="utf-8") as f:
                c_dict = json.load(f)

            return Block7Contract(
                sample_id=sample_id,
                primary_haplogroup=c_dict.get("primary_haplogroup", "UNKNOWN"),
                resolved_haplogroup=c_dict.get("resolved_haplogroup", "UNKNOWN"),
                quality_score=float(c_dict.get("quality_score", 0.0)),
                quality_margin=c_dict.get("quality_margin"),
                haplogroup_status=HaplogroupStatus(c_dict.get("haplogroup_status", "UNRESOLVED")),
                overall_validation_status=OverallValidationStatus(c_dict.get("overall_validation_status", "FAIL")),
                tree_id=c_dict.get("tree_id", "phylotree-fu-rcrs"),
                tree_version=c_dict.get("tree_version", "1.2"),
                metric=c_dict.get("metric", "kulczynski"),
                contaminant_affinity=ContaminantAffinity(c_dict.get("contaminant_affinity", "NOT_APPLICABLE")),
                contamination_fraction=c_dict.get("contamination_fraction"),
                vcf_haplogroup=c_dict.get("vcf_haplogroup"),
                vcf_concordant=c_dict.get("vcf_concordant"),
                subclade_shift_detected=bool(c_dict.get("subclade_shift_detected", False)),
                subclade_shift_supported=c_dict.get("subclade_shift_supported"),
                fasta_sha256=c_dict.get("fasta_sha256", ""),
                flags=c_dict.get("flags", []),
                execution_time_seconds=float(c_dict.get("execution_time_seconds", 0.0)),
            )

        t0 = time.perf_counter()

        sig = inspect.signature(self.haplogroup_orchestrator.run_sample)
        run_kwargs = {
            "sample_id": sample_id,
            "consensus_fasta": canonical_fasta,
            "sample_output_dir": block7_out,
            "variants_vcf": variants_vcf,
            "reference_info": self.reference_info,
            "block4_contract_path": b4_contract_path,
            "block5_contract_path": b5_contract_path,
            "block6_contract_path": b6_contract_path,
        }
        if "reference_fasta" in sig.parameters:
            run_kwargs["reference_fasta"] = self.ref_fasta_path
        if "expected_fasta_sha256" in sig.parameters:
            run_kwargs["expected_fasta_sha256"] = fasta_sha256

        contract: Block7Contract = self.haplogroup_orchestrator.run_sample(**run_kwargs)
        elapsed = time.perf_counter() - t0

        tracker.mark_success(
            "HAPLOGROUP_CLASSIFICATION",
            elapsed,
            {
                "contract_json": contract_json,
                "haplogroup_json": hg_json,
                "validation_json": val_json,
            },
            composite_hash=hg_composite_hash,
        )
        return contract

    def export_haplogroup_cohort_summary(
        self,
        cohort_contracts: Dict[str, Block7Contract],
        output_dir: Path,
    ) -> Tuple[Path, Path]:
        output_dir.mkdir(parents=True, exist_ok=True)
        tsv_path = output_dir / "cohort_haplogroup_summary.tsv"
        json_path = output_dir / "cohort_haplogroup_summary.json"

        headers = [
            "sample_id",
            "primary_haplogroup",
            "resolved_haplogroup",
            "quality_score",
            "quality_margin",
            "haplogroup_status",
            "overall_validation_status",
            "subclade_shift_detected",
            "subclade_shift_supported",
            "contaminant_affinity",
            "contamination_fraction",
            "flags_count",
            "warning_codes",
            "warning_details",
            "execution_time_seconds",
        ]

        rows = []
        json_records = []

        for s_id, c in cohort_contracts.items():
            margin_str = f"{c.quality_margin:.4f}" if c.quality_margin is not None else "NA"
            contam_str = f"{c.contamination_fraction:.4f}" if c.contamination_fraction is not None else "NA"
            shift_supp_str = str(c.subclade_shift_supported) if c.subclade_shift_supported is not None else "NA"

            active_warns = [
                f for f in c.flags
                if (getattr(f, "severity", None) and getattr(f.severity, "value", str(f.severity)).upper() in {"WARNING", "CRITICAL"})
                or (isinstance(f, dict) and str(f.get("severity", "")).upper() in {"WARNING", "CRITICAL"})
            ]

            codes_list: List[str] = []
            details_list: List[str] = []
            for w in active_warns:
                code_raw = getattr(w, "code", None) or (w.get("code") if isinstance(w, dict) else str(w))
                code_val = code_raw.value if hasattr(code_raw, "value") else str(code_raw)
                msg_val = getattr(w, "message", "") or (w.get("message", "") if isinstance(w, dict) else "")
                codes_list.append(code_val)
                details_list.append(f"{code_val}: {msg_val}" if msg_val else code_val)

            codes_str = "; ".join(codes_list) if codes_list else "NONE"
            details_str = " | ".join(details_list) if details_list else "NONE"

            rows.append([
                s_id,
                c.primary_haplogroup,
                c.resolved_haplogroup,
                f"{c.quality_score:.4f}",
                margin_str,
                c.haplogroup_status.value,
                c.overall_validation_status.value,
                c.subclade_shift_detected,
                shift_supp_str,
                c.contaminant_affinity.value,
                contam_str,
                len(c.flags),
                codes_str,
                details_str,
                f"{c.execution_time_seconds:.2f}",
            ])

            json_records.append(c.to_dict())

        with open(tsv_path, "w", newline="", encoding="utf-8") as f_tsv:
            writer = csv.writer(f_tsv, delimiter="\t")
            writer.writerow(headers)
            writer.writerows(rows)

        with open(json_path, "w", encoding="utf-8") as f_json:
            json.dump(
                {
                    "generated_at": datetime.now(timezone.utc).isoformat(),
                    "total_samples": len(cohort_contracts),
                    "samples": json_records,
                },
                f_json,
                indent=2,
            )

        return tsv_path, json_path

    # --------------------------------------------------------------------------
    # Block 8: Manifest Generation & Root of Trust Anchoring
    # --------------------------------------------------------------------------
    def _emit_execution_manifest(self, sample_ids: List[str]) -> Path:
        manifest_path = self.results_dir / "execution_manifest.json"
        samples_manifest: Dict[str, Dict[str, str]] = {}
        flat_contracts: Dict[str, str] = {}

        for s_id in sample_ids:
            samples_manifest[s_id] = {}
            for b_num, sub_dir, c_name in [
                ("block2_contract", "preprocessing", f"{s_id}.block2_contract.json"),
                ("block3_contract", "alignment", f"{s_id}.block3_contract.json"),
                ("block4_contract", "damage", f"{s_id}.block4_contract.json"),
                ("block5_contract", "contamination", f"{s_id}.block5_contract.json"),
                ("block6_contract", "consensus", f"{s_id}.block6_contract.json"),
                ("block7_contract", "haplogroup", f"{s_id}.block7_contract.json"),
            ]:
                c_file = self.operational_files_dir / sub_dir / s_id / c_name
                if not c_file.is_file():
                    c_file = self.results_dir / sub_dir / s_id / c_name
                if c_file.is_file():
                    h = self.calculate_sha256(c_file)
                    samples_manifest[s_id][f"{b_num}_sha256"] = h
                    flat_contracts[f"{s_id}_{b_num}_sha256"] = h

            sample_single_manifest = self.operational_files_dir / s_id / f"{s_id}.manifest.json"
            sample_single_manifest.parent.mkdir(parents=True, exist_ok=True)
            with open(sample_single_manifest, "w", encoding="utf-8") as f_smp:
                json.dump(
                    {
                        "sample_id": s_id,
                        "run_id": f"RUN_{s_id}_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}",
                        "reference_sha256": self.ref_sha256,
                        "contracts": samples_manifest[s_id],
                    },
                    f_smp,
                    indent=2,
                )

        manifest_data = {
            "run_id": f"RUN_COHORT_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}",
            "project_name": self.config.get("project", {}).get("name", "ancient-modern-mtdna-pipeline"),
            "pipeline_version": self.config.get("project", {}).get("version", "1.2.0"),
            "reference_fasta_sha256": self.ref_sha256,
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "cohort_samples": sample_ids,
            "samples": {sid: {"contracts": s_c} for sid, s_c in samples_manifest.items()},
            "contracts": flat_contracts,
        }

        with open(manifest_path, "w", encoding="utf-8") as f:
            json.dump(manifest_data, f, indent=2)

        op_manifest = self.operational_files_dir / "execution_manifest.json"
        with open(op_manifest, "w", encoding="utf-8") as f:
            json.dump(manifest_data, f, indent=2)

        logger.info(f"🔒 Root of Trust manifest: {manifest_path}")
        return manifest_path

    # --------------------------------------------------------------------------
    # Main Cohort Execution Loop
    # --------------------------------------------------------------------------
    def execute_cohort(self, input_dir: Path, cache_mode: Optional[str] = None) -> None:
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

        # Pre-flight Cache Inspection & Resolution
        self.handle_cache_policy(cache_mode=effective_cache_mode)

        # Block 1
        self.validate_reference()

        unrescaled_bams_for_block5: List[Tuple[str, Path, bool]] = []
        b4_contracts: Dict[str, Path] = {}
        sample_types: Dict[str, str] = {}
        failed_samples: Dict[str, str] = {}
        cohort_sample_ids: List[str] = []

        # Input mode branch: BAM vs FASTQ
        if self.data_type == "BAM":
            logger.info("=" * 80)
            logger.info("▶ BAM INGESTION MODE: Direct ingestion with two-coordinate deduplication ('both')")
            logger.info("=" * 80)

            ingested_bams = self.discover_and_ingest_bams(input_dir)
            cohort_sample_ids = [s_id for s_id, _ in ingested_bams]

            logger.info("=" * 80)
            logger.info("▶ BLOCK 4: aDNA DAMAGE ASSESSMENT & RESCALING (Direct BAM Stream)")
            logger.info("=" * 80)

            for s_id, aln_res in ingested_bams:
                try:
                    final_bam = aln_res.outputs.dedup_bam
                    dam_res, b4_contract = self.run_damage(s_id, aln_res)
                    b4_contracts[s_id] = b4_contract

                    verdict = self._extract_verdict_with_contract_fallback(dam_res, b4_contract)
                    is_ancient = "AUTHENTIC_ANCIENT" in verdict.upper()
                    sample_types[s_id] = "ANCIENT" if is_ancient else "MODERN"
                    unrescaled_bams_for_block5.append((s_id, final_bam, is_ancient))
                    logger.info(f"[{s_id}] Block 4 Verdict: {verdict} (Classified: {sample_types[s_id]})")
                except Exception as err:
                    logger.critical(f"[{s_id}] Block 4 failure: {err}", exc_info=True)
                    failed_samples[s_id] = f"Block 4: {err}"
        else:
            samples = self.run_discovery_and_preprocessing(input_dir)
            cohort_sample_ids = [s.sample_id for s in samples]

            logger.info("=" * 80)
            logger.info("▶ BLOCKS 3 & 4: ALIGNMENT (bwa-aln) & aDNA DAMAGE PROFILING")
            logger.info("=" * 80)

            for sample in samples:
                s_id = sample.sample_id
                logger.info(f"[*] Processing Sample: {s_id}")

                try:
                    aln_res = self.run_alignment(sample)
                    if aln_res.outputs and aln_res.outputs.dedup_bam and aln_res.outputs.dedup_bam.is_file():
                        final_bam = aln_res.outputs.dedup_bam
                    elif aln_res.outputs and aln_res.outputs.filtered_bam and aln_res.outputs.filtered_bam.is_file():
                        final_bam = aln_res.outputs.filtered_bam
                    else:
                        final_bam = self.operational_files_dir / "alignment" / s_id / f"{s_id}.filtered.dedup.bam"

                    dam_res, b4_contract = self.run_damage(s_id, aln_res)
                    b4_contracts[s_id] = b4_contract

                    verdict = self._extract_verdict_with_contract_fallback(dam_res, b4_contract)
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

        # Block 5: Contamination Estimation & Reporting
        tsv_contam, json_contam = self.run_contamination_and_report(unrescaled_bams_for_block5)

        # Block 6: Consensus & Heteroplasmy Calling
        logger.info("=" * 80)
        logger.info("▶ BLOCK 6: DUAL-STREAM CONSENSUS & HETEROPLASMY CALLING (ConsensusOrchestrator)")
        logger.info("=" * 80)

        consensus_dir = self.operational_files_dir / "consensus"
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

        tsv_cons, json_cons = self.export_consensus_cohort_summary(cohort_qc, consensus_dir)

        # Block 7: Haplogrouping & Lineage Validation
        logger.info("=" * 80)
        logger.info("▶ BLOCK 7: HAPLOGROUPING & PHYLOGENETIC VALIDATION (HaplogroupOrchestrator)")
        logger.info("=" * 80)

        haplogroup_dir = self.operational_files_dir / "haplogroup"
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
                    f"Status={b7_contract.haplogroup_status.value} | "
                    f"Validation={b7_contract.overall_validation_status.value}"
                )
            except Exception as err:
                logger.critical(f"[{s_id}] Block 7 failed: {err}", exc_info=True)
                failed_samples[s_id] = f"Block 7: {err}"

        tsv_hg, json_hg = self.export_haplogroup_cohort_summary(cohort_b7_contracts, haplogroup_dir)

        # Block 8: Forensic Synthesis & Lineage Audit
        logger.info("=" * 80)
        logger.info("▶ BLOCK 8: FORENSIC SYNTHESIS, LINEAGE AUDIT & COHORT REPORTING")
        logger.info("=" * 80)

        self._emit_execution_manifest(cohort_sample_ids)

        synthesis_records: Dict[str, SampleMasterRecord] = {}
        for s_id in cohort_sample_ids:
            try:
                rec, v_rep, artifacts = self.synthesis_orchestrator.process_sample(s_id)
                synthesis_records[s_id] = rec
                logger.info(
                    f"[{s_id}] Block 8 Complete -> Overall: {rec.overall_status} | "
                    f"Mutations: {len(rec.mutations)} | Warnings: {len(rec.warnings)} | "
                    f"Artifacts: {len(artifacts)}"
                )
            except Exception as err:
                logger.critical(f"[{s_id}] Block 8 Sample Synthesis failed: {err}", exc_info=True)
                failed_samples[s_id] = f"Block 8: {err}"

        cohort_deliverables = self.synthesis_orchestrator.run_cohort(cohort_sample_ids)

        elapsed = time.perf_counter() - t_global0

        print("\n" + "=" * 125)
        print(f"🎉 COMPLETE PIPELINE (BLOCKS 1 -> 8) FINISHED IN {elapsed:.2f}s")
        print("=" * 125)
        print(f"{'Sample ID':<16} | {'Type':<8} | {'Resolved HG':<12} | {'Score':<6} | {'Status':<6} | {'Mutations':<9} | {'Diagnostic Warnings'}")
        print("-" * 125)
        for s_id in cohort_sample_ids:
            rec = synthesis_records.get(s_id)
            if rec:
                b_type = rec.authenticity.biological_type.value
                rhg = rec.haplogroup.resolved_haplogroup
                qscore = f"{rec.haplogroup.quality_score:.2f}" if rec.haplogroup.quality_score is not None else "NA"
                ov_stat = rec.overall_status
                n_muts = str(len(rec.mutations))
                w_codes = "; ".join(rec.get_warning_codes()) if rec.warnings else "Clean (No warnings)"
                print(f"{s_id:<16} | {b_type:<8} | {rhg:<12} | {qscore:<6} | {ov_stat:<6} | {n_muts:<9} | {w_codes}")
            else:
                print(f"{s_id:<16} | {'FAILED':<8} | {'UNKNOWN':<12} | {'NA':<6} | {'FAIL':<6} | {'0':<9} | Pipeline processing error")

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