"""
Configuration loader, path resolver, and schema validator for mtDNA pipeline.

Supports automatic discovery of 'default_config.yaml' (located at config/default_config.yaml),
strongly-typed parameter access across all 8 pipeline blocks, and defensive validation.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_PATH = PROJECT_ROOT / "config" / "default_config.yaml"


def resolve_config_path(config_path: Optional[Union[str, Path]] = None) -> Path:
    """Resolves configuration path, automatically resolving to config/default_config.yaml.

    Args:
        config_path: Optional path, filename, or alias (e.g. 'default_config').

    Returns:
        Resolved absolute Path to configuration YAML file.

    Raises:
        FileNotFoundError: If the target configuration file does not exist.
    """
    if config_path is None or str(config_path).strip() in ("", "default_config", "default_config.yaml"):
        if DEFAULT_CONFIG_PATH.is_file():
            return DEFAULT_CONFIG_PATH
        raise FileNotFoundError(f"Default configuration file not found at: {DEFAULT_CONFIG_PATH}")

    candidate = Path(config_path)

    # 1. Direct path check
    if candidate.is_file():
        return candidate.resolve()

    # 2. Check if a directory was provided (containing default_config.yaml)
    if candidate.is_dir() and (candidate / "default_config.yaml").is_file():
        return (candidate / "default_config.yaml").resolve()

    # 3. Check inside PROJECT_ROOT / config
    cfg_dir_cand = PROJECT_ROOT / "config" / candidate.name
    if cfg_dir_cand.is_file():
        return cfg_dir_cand.resolve()

    cfg_dir_yaml = PROJECT_ROOT / "config" / f"{candidate.stem}.yaml"
    if cfg_dir_yaml.is_file():
        return cfg_dir_yaml.resolve()

    # 4. Check relative to PROJECT_ROOT
    root_cand = PROJECT_ROOT / candidate
    if root_cand.is_file():
        return root_cand.resolve()

    raise FileNotFoundError(f"Configuration file not found: {config_path} (searched at {candidate}, {cfg_dir_cand})")


@dataclass
class PipelineConfig:
    """Strongly typed container holding parsed and validated pipeline parameters across all 8 blocks."""

    # Raw dictionary holding the entire YAML structure
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

    # Block 1 & Reference Management
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

    # Block 2: Preprocessing & Adapter Trimming
    preprocessing_tool: str = "AdapterRemoval"
    min_read_length: int = 30
    min_merge_overlap: int = 11
    collapse_reads: bool = True
    trim_ns: bool = True

    # Block 3: Alignment, Deduplication & Filtering
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

    # Block 4: Damage Authentication
    damage_policy_id: str = "adna_standard_v1"
    input_bam_type: str = "FILTERED_DEDUP"
    allow_input_bam_fallback: bool = False
    library_treatment: str = "NONE"
    pmd_threshold: float = 3.0
    pmd_enabled: bool = True
    rescaling_enabled: bool = True

    # Block 5: Contamination Estimation & Gating
    contamination_enabled: bool = True
    contamination_method: str = "SCHMUTZI"
    allow_fallback_to_dissent: bool = False
    fast_screening_enabled: bool = True
    diagnostic_panel_path: Path = Path("data/reference/diagnostic_sites_phylotree17.tsv")
    schmutzi_database_path: Path = Path("data/reference/schmutzi_db/197")
    schmutzi_threads: int = 3
    schmutzi_max_iterations: int = 10

    # Block 6: Consensus Assembly & Heteroplasmy Calling
    consensus_policy_id: str = "consensus_strict_v1"
    min_consensus_depth: int = 10
    min_depth: int = 10
    major_allele_fraction: float = 0.80
    ambiguity_min_vaf: float = 0.10
    emit_iupac_fasta: bool = True
    heteroplasmy_enabled: bool = True
    het_min_freq: float = 0.10
    het_max_freq: float = 0.80

    # Block 7: Haplogroup Classification & Validation
    haplogroup_tool: str = "haplogrep3"
    tree_id: str = "phylotree-fu-rcrs"
    tree_version: str = "1.2"
    min_quality_score: float = 0.80
    sequence_type: str = "dna"
    max_n_percentage: float = 10.0
    strict_mode: bool = False

    # Cross-Block & Synthesis
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
        """Parses a YAML configuration file with type enforcement and semantic validation.

        If config_path is omitted, None, or 'default_config', resolves to config/default_config.yaml.

        Args:
            config_path: Path, alias, or None for the default configuration.

        Returns:
            Validated PipelineConfig instance.

        Raises:
            FileNotFoundError: If the configuration file cannot be found.
            ValueError: If semantic validation fails.
        """
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

        # 8. Alignment & Deduplication
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

        # 11. Consensus & Heteroplasmy
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

        # --- Defensive Semantic Validations (Fail-Early) ---
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
    """Convenience function to load and validate the pipeline configuration.

    If config_path is omitted, loads config/default_config.yaml.
    """
    return PipelineConfig.from_yaml(config_path)
