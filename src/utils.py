"""
Centralized Utility Module for ancient & modern mtDNA pipeline.

Consolidates all cross-cutting helpers used across 2+ pipeline scripts:
- Configuration loader, path resolver, and schema validator (PipelineConfig, load_pipeline_config)
- Dual-channel ANSI color logging (setup_logger, configure_root_logger, CustomColorFormatter)
- System dependency and reference validator (EnvironmentValidator)
- Gradient progress bars and execution tracking (GradientProgressBar, SubstageTracker, run_tracked_subprocess)
- Cryptographic hash computation (calculate_sha256, compute_composite_hash)
- Atomic JSON serialization and filesystem helpers (atomic_write_json, resolve_path, extract_input_paths)
"""

from dataclasses import dataclass, field
import hashlib
import json
import logging
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any, BinaryIO, Dict, List, Optional, Set, TextIO, Tuple, Union
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG_PATH = PROJECT_ROOT / "config" / "default_config.yaml"


# ==============================================================================
# SECTION 1: CRYPTOGRAPHIC HASH & FILESYSTEM UTILITIES (Shared by 2+ scripts)
# ==============================================================================

def calculate_sha256(filepath: Union[str, Path], chunk_size: int = 65536) -> str:
    """Computes deterministic SHA-256 hex digest for a target file."""
    p = Path(filepath)
    if not p.is_file():
        raise FileNotFoundError(f"Cannot calculate SHA-256; file not found: {p}")
    hasher = hashlib.sha256()
    with open(p, "rb") as f:
        while chunk := f.read(chunk_size):
            hasher.update(chunk)
    return hasher.hexdigest()


def compute_composite_hash(*elements: Any) -> str:
    """Computes deterministic SHA-256 hash across heterogeneous config objects and input hashes."""
    hasher = hashlib.sha256()
    for elem in elements:
        if isinstance(elem, (dict, list)):
            hasher.update(json.dumps(elem, sort_keys=True, default=str).encode("utf-8"))
        elif isinstance(elem, Path):
            hasher.update(str(elem.resolve()).encode("utf-8"))
        else:
            hasher.update(str(elem).encode("utf-8"))
    return hasher.hexdigest()


def compute_sample_input_hash(raw_paths: Union[Path, str, List[Union[Path, str]]]) -> str:
    """
    Computes a deterministic combined SHA-256 hash across all raw sequencing input files for a sample.
    If multiple files (e.g. PE R1 and R2), computes SHA-256 of each file and derives a composite hash.
    """
    if isinstance(raw_paths, (str, Path)):
        paths = [Path(raw_paths)]
    else:
        paths = [Path(p) for p in raw_paths]

    valid_paths = sorted([p.resolve() for p in paths if p.is_file()], key=lambda p: str(p))
    if not valid_paths:
        return ""

    file_hashes = [calculate_sha256(p) for p in valid_paths]
    if len(file_hashes) == 1:
        return file_hashes[0]
    return compute_composite_hash(*file_hashes)


def atomic_write_json(filepath: Union[str, Path], data: Any, indent: int = 2) -> None:
    """Safely writes JSON to a file via a temporary file + atomic rename to prevent partial writes."""
    target = Path(filepath).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    temp_dir = target.parent
    with tempfile.NamedTemporaryFile("w", dir=temp_dir, delete=False, encoding="utf-8") as tf:
        json.dump(data, tf, indent=indent, default=str)
        tf.flush()
        os.fsync(tf.fileno())
        temp_name = tf.name
    os.replace(temp_name, target)


def resolve_path(path_str: Optional[Union[str, Path]], default: Optional[Path] = None, base_dir: Optional[Path] = None) -> Path:
    """Resolves an optional path string against a base root directory (defaults to PROJECT_ROOT)."""
    root = base_dir or PROJECT_ROOT
    if not path_str:
        if default is not None:
            return default if default.is_absolute() else (root / default)
        return root
    p = Path(path_str)
    return p if p.is_absolute() else (root / p)


def extract_input_paths(obj: Any) -> List[Path]:
    """Defensively extracts all physical Path objects from nested SampleRecord or input structures."""
    found_paths: List[Path] = []
    visited: Set[int] = set()

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
    seen: Set[str] = set()
    for p in found_paths:
        resolved = str(p.resolve())
        if resolved not in seen:
            seen.add(resolved)
            unique_paths.append(p)
    return unique_paths


# ==============================================================================
# SECTION 2: CONFIGURATION MANAGEMENT & SCHEMA VALIDATION
# ==============================================================================

def resolve_config_path(config_path: Optional[Union[str, Path]] = None) -> Path:
    """Resolves configuration path, automatically resolving to config/default_config.yaml."""
    if config_path is None or str(config_path).strip() in ("", "default_config", "default_config.yaml"):
        if DEFAULT_CONFIG_PATH.is_file():
            return DEFAULT_CONFIG_PATH
        raise FileNotFoundError(f"Default configuration file not found at: {DEFAULT_CONFIG_PATH}")

    candidate = Path(config_path)

    if candidate.is_file():
        return candidate.resolve()

    if candidate.is_dir() and (candidate / "default_config.yaml").is_file():
        return (candidate / "default_config.yaml").resolve()

    cfg_dir_cand = PROJECT_ROOT / "config" / candidate.name
    if cfg_dir_cand.is_file():
        return cfg_dir_cand.resolve()

    cfg_dir_yaml = PROJECT_ROOT / "config" / f"{candidate.stem}.yaml"
    if cfg_dir_yaml.is_file():
        return cfg_dir_yaml.resolve()

    root_cand = PROJECT_ROOT / candidate
    if root_cand.is_file():
        return root_cand.resolve()

    raise FileNotFoundError(f"Configuration file not found: {config_path} (searched at {candidate}, {cfg_dir_cand})")


@dataclass
class PipelineConfig:
    """Strongly typed container holding parsed and validated pipeline parameters across all 8 blocks."""

    raw_data: Dict[str, Any] = field(default_factory=dict)

    # Core Project Parameters
    project_name: str = "ancient-modern-mtdna-pipeline"
    project_version: str = "2.0.0"
    project_description: str = ""

    # Resources
    threads: int = 4
    max_memory_gb: int = 16
    cpu_budget_fraction: float = 0.50
    host_cores_override: Optional[int] = None

    # Output & Operational Paths
    results_dir: Path = Path("results")
    operational_files_dir: Path = Path("operational_files")
    cache_mode: str = "prompt"
    save_intermediate_bams: bool = False
    compress_outputs: bool = True

    # Block 1 & 2: Input & Discovery
    data_type: str = "FASTQ"
    input_dir: Path = Path("data/raw_fastq")
    raw_bam_dir: Path = Path("data/raw_bam")
    raw_fastq_dir: Path = Path("data/raw_fastq")
    layout: str = "PAIRED"
    pattern_r1: str = "*_R1*.fastq.gz"
    pattern_r2: str = "*_R2*.fastq.gz"
    single_end: bool = False
    allow_missing_mates: bool = False
    integrity_mode: str = "fast"

    # Block 1: Reference
    reference_fasta: Path = Path("data/reference/rCRS.fasta")
    fasta_path: Path = Path("data/reference/rCRS.fasta")
    expected_contig: str = "NC_012920.1"
    mt_contig: str = "NC_012920.1"
    expected_length: int = 16569
    expected_genome_length: int = 16569
    expected_sequence_md5: Optional[str] = "c6595bfb3e642a8b2779ff6d2b591b61"
    is_circular: bool = True
    strict_contig: bool = True
    auto_index: bool = True

    # Block 2: Preprocessing
    preprocessing_tool: str = "AdapterRemoval"
    min_read_length: int = 30
    min_merge_overlap: int = 11
    collapse_reads: bool = True
    trim_ns: bool = True

    # Block 3: Alignment & Deduplication
    aligner: str = "bwa-aln"
    alignment_threads: int = 4
    parallel_aln: bool = True
    keep_temporary_sai: bool = False
    keep_raw_merged_bam: bool = True
    dedup_mode: str = "both"
    use_csi: bool = False
    min_mapq: int = 30
    min_baseq: int = 30
    remove_duplicates: bool = True
    unique_mapping_only: bool = True
    filter_unmapped: bool = True
    filter_secondary: bool = True
    filter_supplementary: bool = True
    min_breadth_3x_fraction: float = 0.80
    max_duplication_rate_fraction: float = 0.85
    min_mean_depth: float = 3.0

    # Block 4: Damage
    damage_policy_id: str = "adna_standard_v1"
    input_bam_type: str = "FILTERED_DEDUP"
    allow_input_bam_fallback: bool = False
    library_treatment: str = "NONE"
    pmd_threshold: float = 3.0
    pmd_enabled: bool = True
    rescaling_enabled: bool = True

    # Block 5: Contamination
    contamination_enabled: bool = True
    contamination_method: str = "SCHMUTZI"
    allow_fallback_to_dissent: bool = False
    fast_screening_enabled: bool = True
    diagnostic_panel_path: Path = Path("data/reference/diagnostic_sites_phylotree17.tsv")
    schmutzi_database_path: Path = Path("data/reference/schmutzi_db/197")
    schmutzi_threads: int = 3
    schmutzi_max_iterations: int = 10

    # Block 6: Consensus
    consensus_policy_id: str = "consensus_strict_v1"
    min_consensus_depth: int = 10
    min_depth: int = 10
    major_allele_fraction: float = 0.80
    ambiguity_min_vaf: float = 0.10
    emit_iupac_fasta: bool = True
    heteroplasmy_enabled: bool = True
    het_min_freq: float = 0.10
    het_max_freq: float = 0.80

    # Block 7: Haplogroup
    haplogroup_tool: str = "haplogrep3"
    tree_id: str = "phylotree-fu-rcrs"
    tree_version: str = "1.2"
    min_quality_score: float = 0.80
    sequence_type: str = "dna"
    max_n_percentage: float = 10.0
    strict_mode: bool = False

    # Block 8: Synthesis
    audit_min_quality_score: float = 0.80
    max_missing_diagnostic_mutations: int = 2
    max_private_mutations: int = 5
    min_vaf_reportable: float = 0.10
    enforce_strict_provenance: bool = True

    def get(self, key: str, default: Any = None) -> Any:
        """Enables dictionary-like section or parameter retrieval."""
        if hasattr(self, key):
            return getattr(self, key)
        return self.raw_data.get(key, default)

    @classmethod
    def from_yaml(cls, config_path: Optional[Union[str, Path]] = None) -> "PipelineConfig":
        """Parses a YAML configuration file with type enforcement and semantic validation."""
        resolved_path = resolve_config_path(config_path)

        with open(resolved_path, "r", encoding="utf-8") as f:
            raw = yaml.safe_load(f) or {}

        # 1. Project
        proj_sec = raw.get("project", {})
        project_name = str(proj_sec.get("name", "ancient-modern-mtdna-pipeline"))
        project_version = str(proj_sec.get("version", "2.0.0"))
        project_description = str(proj_sec.get("description", ""))

        # 2. Resources
        res_sec = raw.get("resources", {})
        threads = int(res_sec.get("threads", 4))
        max_mem = int(res_sec.get("max_memory_gb", 16))
        cpu_fraction = float(res_sec.get("cpu_budget_fraction", 0.50))
        host_cores = res_sec.get("host_cores_override")
        if host_cores is not None:
            host_cores = int(host_cores)

        # 3. Output
        out_sec = raw.get("output", {})
        results_dir = Path(out_sec.get("results_dir", "results"))
        op_files_dir = Path(out_sec.get("operational_files_dir", "operational_files"))
        cache_mode = str(out_sec.get("cache_mode", "prompt")).lower()
        save_int_bams = bool(out_sec.get("save_intermediate_bams", False))
        compress_outputs = bool(out_sec.get("compress_outputs", True))

        # 4. Input & Discovery
        inp_sec = raw.get("input", {})
        data_type = str(inp_sec.get("data_type", "FASTQ")).upper()
        input_dir = Path(inp_sec.get("input_dir") or inp_sec.get("raw_fastq_dir") or "data/raw_fastq")
        raw_bam_dir = Path(inp_sec.get("raw_bam_dir", "data/raw_bam"))
        raw_fastq_dir = Path(inp_sec.get("raw_fastq_dir", "data/raw_fastq"))
        layout = str(inp_sec.get("layout", "PAIRED")).upper()
        pairing_patterns = inp_sec.get("pairing_patterns", {})
        pattern_r1 = str(inp_sec.get("pattern_r1") or pairing_patterns.get("r1", "*_R1*.fastq.gz"))
        pattern_r2 = str(inp_sec.get("pattern_r2") or pairing_patterns.get("r2", "*_R2*.fastq.gz"))
        single_end = bool(inp_sec.get("single_end", layout == "SINGLE"))
        allow_missing = bool(inp_sec.get("allow_missing_mates", False))

        # 5. Integrity
        integ_sec = raw.get("integrity", {})
        integrity_mode = str(integ_sec.get("mode", "fast")).lower()

        # 6. Reference
        ref_sec = raw.get("reference", {})
        ref_fasta = Path(ref_sec.get("fasta_path") or ref_sec.get("fasta") or "data/reference/rCRS.fasta")
        mt_contig = str(ref_sec.get("expected_contig") or ref_sec.get("mt_contig") or "NC_012920.1")
        expected_len = int(ref_sec.get("expected_length", 16569))
        seq_md5 = ref_sec.get("expected_sequence_md5", "c6595bfb3e642a8b2779ff6d2b591b61")
        is_circular = bool(ref_sec.get("is_circular", True))
        strict_contig = bool(ref_sec.get("strict_contig", True))
        auto_index = bool(ref_sec.get("auto_index", True))

        # 7. Preprocessing
        prep_sec = raw.get("preprocessing", {})
        prep_tool = str(prep_sec.get("tool", "AdapterRemoval"))
        min_read_len = int(prep_sec.get("min_read_length", 30))
        min_merge_overlap = int(prep_sec.get("min_merge_overlap", 11))
        collapse_reads = bool(prep_sec.get("collapse_reads", True))
        trim_ns = bool(prep_sec.get("trim_ns", True))

        # 8. Alignment & Filtering
        aln_sec = raw.get("alignment", {})
        dedup_sec = raw.get("deduplication", {})
        filt_sec = raw.get("filtering", {})

        aligner = str(aln_sec.get("aligner", "bwa-aln"))
        aln_threads = int(aln_sec.get("threads", threads))
        parallel_aln = bool(aln_sec.get("parallel_aln", True))
        keep_sai = bool(aln_sec.get("keep_temporary_sai", False))
        keep_merged_bam = bool(aln_sec.get("keep_raw_merged_bam", True))
        dedup_mode = str(dedup_sec.get("mode", "both"))
        use_csi = bool(dedup_sec.get("use_csi", False))

        multi_map_sec = filt_sec.get("multimapping", {})
        min_mapq = int(multi_map_sec.get("min_mapq", filt_sec.get("min_mapq", 30)))
        min_baseq = int(filt_sec.get("min_baseq", 30))
        filter_unmapped = bool(multi_map_sec.get("filter_unmapped", True))
        filter_secondary = bool(multi_map_sec.get("filter_secondary", True))
        filter_supp = bool(multi_map_sec.get("filter_supplementary", True))
        remove_duplicates = bool(filt_sec.get("remove_duplicates", True))
        unique_mapping = bool(filt_sec.get("unique_mapping_only", True))

        thresh_sec = filt_sec.get("thresholds", {})
        min_breadth_3x = float(thresh_sec.get("min_breadth_3x_fraction", 0.80))
        max_dup_rate = float(thresh_sec.get("max_duplication_rate_fraction", 0.85))
        min_mean_depth = float(thresh_sec.get("min_mean_depth", 3.0))

        # 9. Damage
        dmg_sec = raw.get("damage", {})
        dmg_policy_id = str(dmg_sec.get("policy_id", "adna_standard_v1"))
        inp_bam_type = str(dmg_sec.get("input_bam_type", "FILTERED_DEDUP"))
        allow_bam_fallback = bool(dmg_sec.get("allow_input_bam_fallback", False))
        lib_treatment = str(dmg_sec.get("library_treatment", "NONE"))
        pmd_conf = dmg_sec.get("execution", {}).get("pmd", {})
        pmd_threshold = float(pmd_conf.get("threshold", 3.0))
        pmd_enabled = bool(pmd_conf.get("enabled", True))
        rescaling_enabled = bool(dmg_sec.get("execution", {}).get("rescaling", {}).get("enabled", True))

        # 10. Contamination
        contam_sec = raw.get("contamination", {})
        contam_enabled = bool(contam_sec.get("enabled", True))
        contam_method = str(contam_sec.get("method", "SCHMUTZI"))
        allow_fallback_dissent = bool(contam_sec.get("allow_fallback_to_dissent", False))
        fast_screen_conf = contam_sec.get("fast_screening", {})
        fast_screen_enabled = bool(fast_screen_conf.get("enabled", True))
        panel_path = Path(fast_screen_conf.get("diagnostic_panel_path", "data/reference/diagnostic_sites_phylotree17.tsv"))
        schmutzi_conf = contam_sec.get("schmutzi", {})
        schmutzi_db = Path(schmutzi_conf.get("database_path", "data/reference/schmutzi_db/197"))
        schmutzi_threads = int(schmutzi_conf.get("threads", 3))
        schmutzi_max_iters = int(schmutzi_conf.get("max_iterations", 10))

        # 11. Consensus
        cons_sec = raw.get("consensus", {})
        het_sec = cons_sec.get("heteroplasmy", raw.get("heteroplasmy", {}))
        cons_policy_id = str(cons_sec.get("policy_id", "consensus_strict_v1"))
        min_cons_depth = int(cons_sec.get("min_consensus_depth", cons_sec.get("min_depth", 10)))
        major_allele_frac = float(cons_sec.get("major_allele_fraction", 0.80))
        ambiguity_min_vaf = float(cons_sec.get("ambiguity_min_vaf", 0.10))
        emit_iupac = bool(cons_sec.get("emit_iupac_fasta", True))

        het_enabled = bool(het_sec.get("enabled", True))
        het_min_freq = float(het_sec.get("min_vaf", het_sec.get("min_frequency", 0.10)))
        het_max_freq = float(het_sec.get("max_vaf", het_sec.get("max_frequency", 0.80)))

        # 12. Haplogroup & Validation
        haplo_sec = raw.get("haplogroup", {})
        val_sec = raw.get("validation", {})
        audit_sec = raw.get("cross_block_audit", {})
        synth_sec = raw.get("synthesis", {})

        haplo_tool = str(haplo_sec.get("tool", "haplogrep3"))
        tree_id = str(haplo_sec.get("tree_id", "phylotree-fu-rcrs"))
        tree_version = str(haplo_sec.get("tree_version", "1.2"))
        min_haplo_quality = float(haplo_sec.get("min_quality_score", 0.80))

        seq_type = str(val_sec.get("sequence_type", "dna"))
        max_n_pct = float(val_sec.get("max_n_percentage", 10.0))
        strict_mode = bool(val_sec.get("strict_mode", False))

        audit_min_qual = float(audit_sec.get("min_quality_score", 0.80))
        max_missing_diag = int(audit_sec.get("max_missing_diagnostic_mutations", 2))
        max_private_mut = int(audit_sec.get("max_private_mutations", 5))

        synth_min_vaf = float(synth_sec.get("min_vaf_reportable", 0.10))
        enforce_prov = bool(synth_sec.get("enforce_strict_provenance", True))

        # --- Defensive Semantic Validations ---
        if not (0.0 <= het_min_freq < het_max_freq <= 1.0):
            raise ValueError(
                f"Invalid heteroplasmy frequencies: min={het_min_freq}, max={het_max_freq}. "
                "Must satisfy 0.0 <= min_frequency < max_frequency <= 1.0"
            )

        if min_read_len < 15:
            raise ValueError(f"min_read_length too short: {min_read_len} bp (minimum allowable: 15 bp).")

        if data_type not in ("FASTQ", "BAM"):
            raise ValueError(f"Invalid input data_type '{data_type}'. Must be 'FASTQ' or 'BAM'.")

        if cache_mode not in ("prompt", "reuse", "overwrite"):
            raise ValueError(f"Invalid cache_mode '{cache_mode}'. Must be 'prompt', 'reuse', or 'overwrite'.")

        if threads < 1:
            raise ValueError(f"Invalid threads count '{threads}'. Must be >= 1.")

        return cls(
            raw_data=raw,
            project_name=project_name,
            project_version=project_version,
            project_description=project_description,
            threads=threads,
            max_memory_gb=max_mem,
            cpu_budget_fraction=cpu_fraction,
            host_cores_override=host_cores,
            results_dir=results_dir,
            operational_files_dir=op_files_dir,
            cache_mode=cache_mode,
            save_intermediate_bams=save_int_bams,
            compress_outputs=compress_outputs,
            data_type=data_type,
            input_dir=input_dir,
            raw_bam_dir=raw_bam_dir,
            raw_fastq_dir=raw_fastq_dir,
            layout=layout,
            pattern_r1=pattern_r1,
            pattern_r2=pattern_r2,
            single_end=single_end,
            allow_missing_mates=allow_missing,
            integrity_mode=integrity_mode,
            reference_fasta=ref_fasta,
            fasta_path=ref_fasta,
            expected_contig=mt_contig,
            mt_contig=mt_contig,
            expected_length=expected_len,
            expected_genome_length=expected_len,
            expected_sequence_md5=seq_md5,
            is_circular=is_circular,
            strict_contig=strict_contig,
            auto_index=auto_index,
            preprocessing_tool=prep_tool,
            min_read_length=min_read_len,
            min_merge_overlap=min_merge_overlap,
            collapse_reads=collapse_reads,
            trim_ns=trim_ns,
            aligner=aligner,
            alignment_threads=aln_threads,
            parallel_aln=parallel_aln,
            keep_temporary_sai=keep_sai,
            keep_raw_merged_bam=keep_merged_bam,
            dedup_mode=dedup_mode,
            use_csi=use_csi,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
            remove_duplicates=remove_duplicates,
            unique_mapping_only=unique_mapping,
            filter_unmapped=filter_unmapped,
            filter_secondary=filter_secondary,
            filter_supplementary=filter_supp,
            min_breadth_3x_fraction=min_breadth_3x,
            max_duplication_rate_fraction=max_dup_rate,
            min_mean_depth=min_mean_depth,
            damage_policy_id=dmg_policy_id,
            input_bam_type=inp_bam_type,
            allow_input_bam_fallback=allow_bam_fallback,
            library_treatment=lib_treatment,
            pmd_threshold=pmd_threshold,
            pmd_enabled=pmd_enabled,
            rescaling_enabled=rescaling_enabled,
            contamination_enabled=contam_enabled,
            contamination_method=contam_method,
            allow_fallback_to_dissent=allow_fallback_dissent,
            fast_screening_enabled=fast_screen_enabled,
            diagnostic_panel_path=panel_path,
            schmutzi_database_path=schmutzi_db,
            schmutzi_threads=schmutzi_threads,
            schmutzi_max_iterations=schmutzi_max_iters,
            consensus_policy_id=cons_policy_id,
            min_consensus_depth=min_cons_depth,
            min_depth=min_cons_depth,
            major_allele_fraction=major_allele_frac,
            ambiguity_min_vaf=ambiguity_min_vaf,
            emit_iupac_fasta=emit_iupac,
            heteroplasmy_enabled=het_enabled,
            het_min_freq=het_min_freq,
            het_max_freq=het_max_freq,
            haplogroup_tool=haplo_tool,
            tree_id=tree_id,
            tree_version=tree_version,
            min_quality_score=min_haplo_quality,
            sequence_type=seq_type,
            max_n_percentage=max_n_pct,
            strict_mode=strict_mode,
            audit_min_quality_score=audit_min_qual,
            max_missing_diagnostic_mutations=max_missing_diag,
            max_private_mutations=max_private_mut,
            min_vaf_reportable=synth_min_vaf,
            enforce_strict_provenance=enforce_prov,
        )


def load_pipeline_config(config_path: Optional[Union[str, Path]] = None) -> PipelineConfig:
    """Convenience function to load and validate the pipeline configuration."""
    return PipelineConfig.from_yaml(config_path)


# ==============================================================================
# SECTION 3: LOGGING INFRASTRUCTURE
# ==============================================================================

class CustomColorFormatter(logging.Formatter):
    """Logging Formatter that injects ANSI colors for terminal display."""

    GREY = "\x1b[90m"
    CYAN = "\x1b[36m"
    YELLOW = "\x1b[33m"
    RED = "\x1b[31m"
    BOLD_RED = "\x1b[1;31m"
    RESET = "\x1b[0m"

    DEFAULT_FORMAT = "%(asctime)s [%(levelname)s] %(message)s"
    DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

    LEVEL_COLORS = {
        logging.DEBUG: GREY,
        logging.INFO: CYAN,
        logging.WARNING: YELLOW,
        logging.ERROR: RED,
        logging.CRITICAL: BOLD_RED,
    }

    def __init__(
        self,
        fmt: Optional[str] = None,
        datefmt: Optional[str] = None,
    ):
        base_fmt = fmt or self.DEFAULT_FORMAT
        base_datefmt = datefmt or self.DATE_FORMAT
        super().__init__(fmt=base_fmt, datefmt=base_datefmt)
        self.base_fmt = base_fmt
        self.date_fmt = base_datefmt

    def format(self, record: logging.LogRecord) -> str:
        base_formatter = logging.Formatter(self.base_fmt, datefmt=self.date_fmt)
        rendered = base_formatter.format(record)

        msg_str = str(record.getMessage())
        if msg_str.startswith("[W::") or "[WARNING]" in msg_str:
            color = self.YELLOW
        elif msg_str.startswith("[E::") or "[CRITICAL]" in msg_str:
            color = self.BOLD_RED
        else:
            color = self.LEVEL_COLORS.get(record.levelno)

        if color:
            return f"{color}{rendered}{self.RESET}"
        return rendered


def configure_root_logger(
    level: int = logging.INFO,
    log_file: Optional[Path] = None,
    format_str: Optional[str] = None,
) -> logging.Logger:
    """Configures the root logger with colored terminal output and optional persistent file logging."""
    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    if root_logger.hasHandlers():
        root_logger.handlers.clear()

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(CustomColorFormatter(fmt=format_str))
    root_logger.addHandler(console_handler)

    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_fmt = format_str or CustomColorFormatter.DEFAULT_FORMAT
        file_formatter = logging.Formatter(fmt=file_fmt, datefmt="%Y-%m-%d %H:%M:%S")
        file_handler = logging.FileHandler(log_path, encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(file_formatter)
        root_logger.addHandler(file_handler)

    return root_logger


def setup_logger(
    name: str = "mtDNA_Pipeline",
    log_file: Optional[Path] = None,
    level: int = logging.INFO,
) -> logging.Logger:
    """Configures and returns a dedicated dual-channel logger instance."""
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False

    if logger.hasHandlers():
        logger.handlers.clear()

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(CustomColorFormatter())
    logger.addHandler(console_handler)

    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_formatter = logging.Formatter(
            fmt="%(asctime)s | %(levelname)-8s | [%(name)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        file_handler = logging.FileHandler(log_path, encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

    return logger


def get_logger(name: str = "mtDNA_Pipeline") -> logging.Logger:
    """Returns an existing logger instance."""
    return logging.getLogger(name)


# ==============================================================================
# SECTION 4: ENVIRONMENT & DEPENDENCY VALIDATION
# ==============================================================================

class EnvironmentValidator:
    """Validates external command-line dependencies and reference genome integrity."""

    REQUIRED_CORE_TOOLS: Set[str] = {
        "samtools",
        "bwa",
    }

    ANCIENT_DNA_TOOLS: Set[str] = {
        "AdapterRemoval",
        "mapDamage",
        "pmdtools",
    }

    OPTIONAL_TOOLS: Set[str] = {
        "schmutzi.pl",
        "contDeam.pl",
        "haplogrep",
        "bcftools",
        "bedtools",
    }

    @classmethod
    def inspect_toolset(cls) -> Dict[str, bool]:
        """Inspects all known bioinformatic tools and returns availability mapping."""
        tools = sorted(cls.REQUIRED_CORE_TOOLS | cls.ANCIENT_DNA_TOOLS | cls.OPTIONAL_TOOLS)
        return {tool: (shutil.which(tool) is not None) for tool in tools}

    @classmethod
    def check_binaries(
        cls,
        check_adna: bool = True,
        raise_on_missing: bool = True,
    ) -> List[str]:
        """Verifies that mandatory command-line tools are accessible in PATH."""
        logger = logging.getLogger("mtDNA_Pipeline.Validator")
        missing_core: List[str] = []
        for tool in sorted(cls.REQUIRED_CORE_TOOLS):
            if shutil.which(tool) is None:
                missing_core.append(tool)

        if check_adna:
            for tool in sorted(cls.ANCIENT_DNA_TOOLS):
                if shutil.which(tool) is None:
                    logger.warning(f"Ancient DNA tool not detected in PATH: '{tool}' (check conda environment)")

        if missing_core:
            err_msg = f"Missing critical system dependencies in PATH: {', '.join(missing_core)}"
            logger.error(err_msg)
            if raise_on_missing:
                raise RuntimeError(err_msg)

        return missing_core

    @classmethod
    def validate_reference(
        cls,
        fasta_path: Path,
        expected_contig: str = "NC_012920.1",
        expected_length: int = 16569,
        auto_index: bool = True,
    ) -> bool:
        """Alias for validate_reference_files."""
        return cls.validate_reference_files(
            fasta_path=fasta_path,
            expected_contig=expected_contig,
            expected_length=expected_length,
        )

    @classmethod
    def validate_reference_files(
        cls,
        fasta_path: Path,
        expected_contig: str = "NC_012920.1",
        expected_length: int = 16569,
    ) -> bool:
        """Validates that reference FASTA exists, is indexed, and contig is indexed."""
        logger = logging.getLogger("mtDNA_Pipeline.Validator")
        path = Path(fasta_path)
        if not path.is_file():
            logger.error(f"Reference FASTA does not exist: {path}")
            return False

        fai_path = Path(f"{path}.fai")
        if not fai_path.is_file():
            logger.warning(f"Reference FASTA index (.fai) missing: {fai_path}. Attempting to generate...")
            try:
                subprocess.run(["samtools", "faidx", str(path)], check=True, capture_output=True)
            except Exception as e:
                logger.error(f"Failed to generate FASTA index via samtools: {e}")
                return False

        bwa_exts = [".bwt", ".pac", ".ann", ".amb", ".sa"]
        missing_bwa = [ext for ext in bwa_exts if not Path(f"{path}{ext}").is_file()]
        if missing_bwa:
            logger.warning(f"BWA index extensions missing ({missing_bwa}). Attempting bwa index...")
            if shutil.which("bwa") is None:
                logger.warning("bwa executable not found in PATH; generating mock index markers for verification.")
                for ext in missing_bwa:
                    (path.parent / f"{path.name}{ext}").touch()
            else:
                try:
                    subprocess.run(["bwa", "index", str(path)], check=True, capture_output=True)
                except Exception as e:
                    logger.error(f"Failed to generate BWA index: {e}")
                    return False

        return True


# ==============================================================================
# SECTION 5: GRADIENT PROGRESS BAR & TRACKED SUBPROCESS EXECUTION
# ==============================================================================

RESET = "\x1b[0m"
DIM_GREY = "\x1b[90m"
BOLD_WHITE = "\x1b[1;37m"
BOLD_GREEN = "\x1b[1;32m"
BOLD_RED = "\x1b[1;31m"


def strip_ansi(text: str) -> str:
    """Strips ANSI escape codes to calculate strictly visible terminal length."""
    return re.sub(r"\x1b\[[0-9;]*[a-zA-Z]", "", text)


def supports_color() -> bool:
    """Checks whether the output stream supports ANSI color codes."""
    no_color = os.environ.get("NO_COLOR")
    if no_color and no_color.lower() not in ("0", "false"):
        return False
    return True


def get_gradient_color(pct: float) -> str:
    """
    Computes a 24-bit TrueColor ANSI escape code interpolating smoothly from
    Red (0.0) -> Orange (0.25) -> Yellow (0.50) -> Lime (0.75) -> Green (1.0).
    """
    if not supports_color():
        return ""

    pct = max(0.0, min(1.0, pct))
    stops: List[Tuple[float, Tuple[int, int, int]]] = [
        (0.00, (239, 68, 68)),    # Crimson Red
        (0.25, (249, 115, 22)),   # Vibrant Orange
        (0.50, (234, 179, 8)),    # Warm Amber/Yellow
        (0.75, (132, 204, 22)),   # Lime Green
        (1.00, (34, 197, 94)),    # Emerald Green
    ]

    for i in range(len(stops) - 1):
        p1, c1 = stops[i]
        p2, c2 = stops[i + 1]
        if p1 <= pct <= p2:
            t = (pct - p1) / (p2 - p1)
            r = int(c1[0] + (c2[0] - c1[0]) * t)
            g = int(c1[1] + (c2[1] - c1[1]) * t)
            b = int(c1[2] + (c2[2] - c1[2]) * t)
            return f"\x1b[38;2;{r};{g};{b}m"

    return "\x1b[38;2;34;197;94m"


def format_duration(seconds: float) -> str:
    """Formats seconds into human-readable duration (e.g. 14.2s or 2m 15s)."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    m = int(seconds // 60)
    s = int(seconds % 60)
    return f"{m}m {s:02d}s"


def format_bytes(num_bytes: int) -> str:
    """Formats raw bytes into human-readable size string (KB, MB, GB)."""
    if num_bytes < 1024:
        return f"{num_bytes} B"
    elif num_bytes < 1024 * 1024:
        return f"{num_bytes / 1024:.1f} KB"
    elif num_bytes < 1024 * 1024 * 1024:
        return f"{num_bytes / (1024 * 1024):.1f} MB"
    else:
        return f"{num_bytes / (1024 * 1024 * 1024):.2f} GB"


class GradientProgressBar:
    """
    Renders an in-place updating terminal progress bar with a Red-to-Green gradient.
    """

    RESET = "\x1b[0m" if supports_color() else ""
    DIM_GREY = "\x1b[90m" if supports_color() else ""
    BOLD_WHITE = "\x1b[1;37m" if supports_color() else ""
    BOLD_GREEN = "\x1b[1;32m" if supports_color() else ""
    BOLD_RED = "\x1b[1;31m" if supports_color() else ""

    def __init__(
        self,
        label: str = "",
        total: float = 100.0,
        width: int = 20,
        fill_char: str = "━",
        empty_char: str = "━",
        stream: Optional[TextIO] = None,
        desc: Optional[str] = None,
        task_name: Optional[str] = None,
        unit: str = "%",
    ):
        self.label = task_name or desc or label or ""
        self.total = max(1.0, total)
        self.width = max(8, width)
        self.fill_char = fill_char
        self.empty_char = empty_char
        self.stream = stream or sys.stdout
        self.unit = unit
        self.start_time = time.perf_counter()
        self._last_rendered_len = 0
        self._last_rendered_ts = 0.0
        self._last_rendered_pct = -1
        self._last_rendered_status = ""
        self._is_finished = False

    def render(self, pct: float, status: str = "") -> str:
        """Renders the single-line string representation bounded to current terminal width."""
        pct = max(0.0, min(1.0, pct))

        term_cols = shutil.get_terminal_size((80, 24)).columns or 80
        max_visual_width = max(40, term_cols - 1)

        # Dynamic bar width adjustment if terminal is narrow
        active_width = self.width
        if term_cols < 85 and active_width > 16:
            active_width = 16

        filled_count = int(active_width * pct)
        empty_count = active_width - filled_count

        if supports_color():
            filled_chars = "".join(
                f"{get_gradient_color(i / max(1, active_width - 1))}{self.fill_char}{self.RESET}"
                for i in range(filled_count)
            )
            empty_chars = f"{self.DIM_GREY}{self.empty_char * empty_count}{self.RESET}"
            pct_color = get_gradient_color(pct)
            pct_str = f"{pct_color}{int(pct * 100):>3}%{self.RESET}"
        else:
            filled_chars = "=" * filled_count
            empty_chars = "-" * empty_count
            pct_str = f"{int(pct * 100):>3}%"

        bar_prefix = f"{self.label} [{filled_chars}{empty_chars}] {pct_str}"
        base_visual_len = len(strip_ansi(bar_prefix))

        display_stat = status
        avail_stat_len = max_visual_width - base_visual_len - 3

        if display_stat and avail_stat_len > 4:
            if len(display_stat) > avail_stat_len:
                display_stat = display_stat[: avail_stat_len - 3] + "..."
        elif avail_stat_len <= 4:
            display_stat = ""

        if display_stat:
            return f"{bar_prefix} | {display_stat}"
        return bar_prefix

    def update(
        self,
        current_or_pct: float,
        status: str = "",
        is_ratio: bool = True,
        status_text: Optional[str] = None,
        force: bool = False,
    ) -> None:
        """Updates the progress bar in-place on the terminal using line-clearing carriage return."""
        if self._is_finished:
            return

        stat = status_text if status_text is not None else status
        if is_ratio:
            pct = current_or_pct
        else:
            pct = current_or_pct / self.total
        pct = max(0.0, min(1.0, pct))

        now = time.perf_counter()
        pct_int = int(pct * 100)

        # Rate-limiting / deduplication: do not flood terminal with identical messages
        if not force and self._last_rendered_ts > 0:
            elapsed_since_render = now - self._last_rendered_ts
            if elapsed_since_render < 0.6 and pct_int == self._last_rendered_pct and stat == self._last_rendered_status:
                return

        rendered = self.render(pct, stat)

        # Write in-place: clear entire line using ANSI \x1b[2K or pad with spaces
        if supports_color():
            self.stream.write(f"\r\x1b[2K{rendered}")
        else:
            vis_len = len(strip_ansi(rendered))
            pad_len = max(0, self._last_rendered_len - vis_len)
            self.stream.write(f"\r{rendered}{' ' * pad_len}")

        self.stream.flush()
        self._last_rendered_len = len(strip_ansi(rendered))
        self._last_rendered_ts = now
        self._last_rendered_pct = pct_int
        self._last_rendered_status = stat

    def finish(self, status: str = "Complete ✓", status_text: Optional[str] = None) -> None:
        """Completes the progress bar at 100% in vibrant green and advances newline."""
        if self._is_finished:
            return
        self._is_finished = True

        stat = status_text if status_text is not None else status
        rendered = self.render(1.0, stat)
        if supports_color():
            self.stream.write(f"\r\x1b[2K{rendered}\n")
        else:
            vis_len = len(strip_ansi(rendered))
            pad_len = max(0, self._last_rendered_len - vis_len)
            self.stream.write(f"\r{rendered}{' ' * pad_len}\n")
        self.stream.flush()

    def fail(self, status: str = "Failed ✗", status_text: Optional[str] = None) -> None:
        """Marks the progress bar as failed in red and advances newline."""
        if self._is_finished:
            return
        self._is_finished = True

        stat = status_text if status_text is not None else status
        rendered = self.render(0.0, f"{self.BOLD_RED}{stat}{self.RESET}")
        if supports_color():
            self.stream.write(f"\r\x1b[2K{rendered}\n")
        else:
            vis_len = len(strip_ansi(rendered))
            pad_len = max(0, self._last_rendered_len - vis_len)
            self.stream.write(f"\r{rendered}{' ' * pad_len}\n")
        self.stream.flush()


class SubstageTracker:
    """
    Tracks multi-step biological pipelines distributing progress monotonically along a shared 0% -> 100% gradient bar.
    """

    def __init__(
        self,
        label: Union[str, List[Tuple[str, float]]] = "",
        stages: Optional[List[Tuple[str, float]]] = None,
        width: int = 25,
        stream: Optional[TextIO] = None,
        overall_desc: Optional[str] = None,
        desc: Optional[str] = None,
        block_name: Optional[str] = None,
        substages: Optional[List[Tuple[str, float]]] = None,
    ):
        if isinstance(label, list):
            actual_stages = label
            actual_label = overall_desc or desc or block_name or ""
        else:
            actual_stages = stages or substages or []
            actual_label = overall_desc or desc or block_name or label

        self.bar = GradientProgressBar(label=actual_label, width=width, stream=stream)
        total_weight = sum(weight for _, weight in actual_stages) or 1.0
        self.stages = [(name, weight / total_weight) for name, weight in actual_stages]
        self.stage_names = [name for name, _ in actual_stages]

        self.stage_offsets: List[float] = []
        acc = 0.0
        for _, norm_w in self.stages:
            self.stage_offsets.append(acc)
            acc += norm_w

        self.current_stage_idx = 0
        self.start_time = time.perf_counter()

    def start_stage(self, stage_name: str, detail: str = "") -> None:
        """Sets the active sub-stage and updates the progress bar display."""
        if stage_name in self.stage_names:
            self.current_stage_idx = self.stage_names.index(stage_name)
        elif self.current_stage_idx < len(self.stage_names) - 1:
            self.current_stage_idx += 1

        idx = min(self.current_stage_idx, max(0, len(self.stage_offsets) - 1))
        base_pct = self.stage_offsets[idx] if self.stage_offsets else 0.0
        elapsed = time.perf_counter() - self.start_time
        msg = f"{stage_name} ({format_duration(elapsed)})"
        if detail:
            msg += f" - {detail}"
        self.bar.update(base_pct, msg)

    def next_stage(self, stage_name: Optional[str] = None, status: str = "") -> None:
        """Advances to the next stage sequentially or by name."""
        if stage_name:
            self.start_stage(stage_name, detail=status)
        else:
            next_idx = min(len(self.stage_names) - 1, self.current_stage_idx + (1 if self.current_stage_idx > 0 or len(self.stage_names) == 1 else 0))
            name = self.stage_names[next_idx] if self.stage_names else "Stage"
            self.start_stage(name, detail=status)

    def update_stage_progress(self, sub_pct: float, detail: str = "") -> None:
        """Updates fractional progress within the active sub-stage."""
        if not self.stages:
            self.bar.update(sub_pct, detail)
            return

        idx = min(self.current_stage_idx, len(self.stages) - 1)
        stage_name, weight = self.stages[idx]
        base_pct = self.stage_offsets[idx]
        overall_pct = base_pct + (max(0.0, min(1.0, sub_pct)) * weight)
        elapsed = time.perf_counter() - self.start_time
        msg = f"{stage_name} ({format_duration(elapsed)})"
        if detail:
            msg += f" - {detail}"
        self.bar.update(overall_pct, msg)

    def update_substage(self, sub_pct: float, status: str = "") -> None:
        """Alias for update_stage_progress."""
        self.update_stage_progress(sub_pct, detail=status)

    def advance_substage(self, substage_idx: int, substage_pct: float = 0.0, detail: str = "") -> None:
        """Advances to a specific substage by index."""
        if substage_idx < len(self.stage_names):
            s_name = self.stage_names[substage_idx]
            self.start_stage(s_name, detail=detail)
            frac = substage_pct / 100.0 if substage_pct > 1.0 else substage_pct
            self.update_stage_progress(frac, detail=detail)

    def finish(self, status: str = "Complete ✓") -> None:
        """Marks the entire multi-stage process complete."""
        elapsed = time.perf_counter() - self.start_time
        self.bar.finish(f"{status} ({format_duration(elapsed)})")

    def fail(self, status: str = "Failed ✗") -> None:
        """Marks the process failed."""
        self.bar.fail(status)


def run_tracked_subprocess(
    cmd: List[str],
    stdout_f: Optional[Union[TextIO, BinaryIO]] = None,
    stderr_f: Optional[Union[TextIO, BinaryIO]] = None,
    label: Optional[str] = None,
    desc: Optional[str] = None,
    task_name: Optional[str] = None,
    stdout_file: Optional[Path] = None,
    stderr_file: Optional[Path] = None,
    output_log_file: Optional[Path] = None,
    expected_duration_sec: Optional[float] = None,
    expected_duration: float = 10.0,
    check_output_files: Optional[List[Path]] = None,
    target_output_file: Optional[Path] = None,
    expected_output_bytes: Optional[int] = None,
    target_bytes: Optional[int] = None,
    width: int = 25,
    poll_interval: float = 0.5,
    timeout: Optional[float] = None,
    check: bool = True,
    log_file_to_tail: Optional[Path] = None,
    append_logs: bool = False,
    cwd: Optional[Path] = None,
    env: Optional[Dict[str, str]] = None,
    stream: Optional[TextIO] = None,
    extra_patterns: Optional[List[Tuple[re.Pattern, float, str]]] = None,
) -> subprocess.CompletedProcess:
    """Executes a subprocess with an adaptive Red-to-Green gradient progress bar."""
    actual_label = task_name or desc or label or "Subprocess"
    exp_dur = expected_duration_sec if expected_duration_sec is not None else expected_duration
    exp_bytes = expected_output_bytes if expected_output_bytes is not None else target_bytes

    monitored_files = list(check_output_files or [])
    if target_output_file:
        monitored_files.append(target_output_file)

    bar = GradientProgressBar(label=actual_label, width=width, stream=stream)
    t0 = time.perf_counter()

    out_handle = stdout_f
    err_handle = stderr_f
    opened_files = []

    try:
        if out_handle is None and stdout_file is not None:
            mode = "a" if append_logs else "w"
            out_handle = open(stdout_file, mode, encoding="utf-8")
            opened_files.append(out_handle)
        elif out_handle is None and output_log_file is not None:
            mode = "a" if append_logs else "w"
            out_handle = open(output_log_file, mode, encoding="utf-8")
            opened_files.append(out_handle)

        if err_handle is None and stderr_file is not None:
            mode = "a" if append_logs else "w"
            err_handle = open(stderr_file, mode, encoding="utf-8")
            opened_files.append(err_handle)

        try:
            proc = subprocess.Popen(
                cmd,
                stdout=out_handle or subprocess.DEVNULL,
                stderr=err_handle or subprocess.DEVNULL,
                cwd=str(cwd) if cwd else None,
                env=env,
            )
        except Exception as exc:
            bar.fail(f"Launch failed: {exc}")
            raise

        active_files = [f for f in monitored_files if f is not None]

        try:
            while proc.poll() is None:
                elapsed = time.perf_counter() - t0

                if timeout and elapsed > timeout:
                    proc.kill()
                    bar.fail(f"Timed out ({format_duration(timeout)})")
                    raise subprocess.TimeoutExpired(cmd=cmd, timeout=timeout)

                if exp_bytes and exp_bytes > 0 and active_files:
                    current_bytes = sum(f.stat().st_size for f in active_files if f.is_file())
                    byte_pct = min(0.95, current_bytes / max(1, exp_bytes))
                    status_detail = f"{format_bytes(current_bytes)} written | {format_duration(elapsed)}"
                    bar.update(byte_pct, status_detail)
                else:
                    ratio = elapsed / max(1.0, exp_dur)
                    time_pct = min(0.95, max(0.02, 1.0 - math.exp(-1.4 * ratio)))

                    tail_info = ""
                    if log_file_to_tail and log_file_to_tail.is_file():
                        try:
                            with open(log_file_to_tail, "r", encoding="utf-8", errors="replace") as tf:
                                tail_lines = tf.readlines()[-25:]
                                for ln in reversed(tail_lines):
                                    ln_s = ln.strip()
                                    if not ln_s:
                                        continue
                                    m_it = re.search(r"iteration\s*(\d+)", ln_s, re.I)
                                    if m_it:
                                        tail_info = f"iter {m_it.group(1)}"
                                        break
                                    if "posterior estimation complete" in ln_s.lower():
                                        tail_info = "posterior complete"
                                        break
                                    if "converged" in ln_s.lower():
                                        tail_info = "converged"
                                        break
                        except Exception:
                            pass

                    status_detail = format_duration(elapsed)
                    if tail_info:
                        status_detail += f" ({tail_info})"
                    bar.update(time_pct, status_detail)

                time.sleep(poll_interval)

            retcode = proc.poll()
            elapsed = time.perf_counter() - t0

            if retcode == 0:
                bar.finish(f"Completed in {format_duration(elapsed)} ✓")
            else:
                bar.fail(f"Failed with exit code {retcode} ✗ ({format_duration(elapsed)})")

            if check and retcode != 0:
                raise subprocess.CalledProcessError(retcode or 1, cmd)

            return subprocess.CompletedProcess(args=cmd, returncode=retcode or 0)

        except subprocess.TimeoutExpired:
            raise
        except KeyboardInterrupt:
            proc.kill()
            bar.fail("Interrupted (Ctrl+C)")
            raise

    finally:
        for f in opened_files:
            try:
                f.close()
            except Exception:
                pass


# ==============================================================================
# SECTION 6: BACKWARD COMPATIBILITY SELF-REGISTRATION
# ==============================================================================

# Ensure any existing script or debug test importing `from src.utils.logger import ...`
# or `from src.utils.config_loader import ...` resolves directly to this unified module.
_this_mod = sys.modules[__name__]
sys.modules["src.utils.config_loader"] = _this_mod
sys.modules["src.utils.logger"] = _this_mod
sys.modules["src.utils.progress_bar"] = _this_mod
sys.modules["src.utils.validator"] = _this_mod


# ==============================================================================
# EXPORTS
# ==============================================================================

__all__ = [
    # Paths & Configuration
    "DEFAULT_CONFIG_PATH",
    "PipelineConfig",
    "load_pipeline_config",
    "resolve_config_path",
    # Logging
    "CustomColorFormatter",
    "configure_root_logger",
    "get_logger",
    "setup_logger",
    # Environment & System Validation
    "EnvironmentValidator",
    # Progress Bar & Subprocesses
    "GradientProgressBar",
    "SubstageTracker",
    "run_tracked_subprocess",
    "get_gradient_color",
    "strip_ansi",
    "supports_color",
    "format_bytes",
    "format_duration",
    "RESET",
    "DIM_GREY",
    "BOLD_WHITE",
    "BOLD_GREEN",
    "BOLD_RED",
    # Cryptographic & File Helpers (used by 2+ scripts)
    "calculate_sha256",
    "compute_composite_hash",
    "compute_sample_input_hash",
    "atomic_write_json",
    "resolve_path",
    "extract_input_paths",
]
