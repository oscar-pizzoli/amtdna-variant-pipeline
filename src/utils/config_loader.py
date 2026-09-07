"""Configuration loader and schema validator for mtDNA pipeline."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict
import yaml


@dataclass
class PipelineConfig:
    """Strongly typed container holding parsed and validated pipeline parameters."""

    # Dati grezzi per moduli downstream avanzati
    raw_data: Dict[str, Any] = field(default_factory=dict)

    # Parametri Core & Ambiente (Blocco 1)
    project_name: str = "ancient-modern-mtdna-pipeline"
    threads: int = 4
    results_dir: Path = Path("results")
    
    # Input & Discovery (Blocco 2)
    data_type: str = "FASTQ"
    input_dir: Path = Path("data/raw_fastq")
    pattern_r1: str = "*_R1*.fastq.gz"
    pattern_r2: str = "*_R2*.fastq.gz"
    single_end: bool = False
    allow_missing_mates: bool = False
    integrity_mode: str = "fast"

    # Reference Genome (Blocco 1 & 3)
    reference_fasta: Path = Path("reference/rCRS.fasta")
    mt_contig: str = "NC_012920.1"
    auto_index: bool = True

    # Preprocessing & Trimming (Blocco 2)
    min_read_length: int = 30
    min_merge_overlap: int = 11
    collapse_reads: bool = True
    trim_ns: bool = True

    # Alignment & Filtering (Blocco 3)
    min_mapq: int = 30
    min_baseq: int = 30
    remove_duplicates: bool = True
    unique_mapping_only: bool = True

    # Consensus & Eteroplasmie (Blocco 6)
    min_depth: int = 10
    het_min_freq: float = 0.20
    het_max_freq: float = 0.80

    # Validazione & Controllo Qualità (Blocco 7)
    strict_mode: bool = False
    expected_genome_length: int = 16569

    @classmethod
    def from_yaml(cls, config_path: Path) -> "PipelineConfig":
        """Parses a YAML configuration file with type enforcement and semantic validation.

        Args:
            config_path: Path to configuration YAML file.

        Returns:
            Validated PipelineConfig instance.

        Raises:
            FileNotFoundError: If the YAML file does not exist.
            ValueError: If semantic validation fails.
        """
        path = Path(config_path).resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Configuration file not found: {path}")

        with open(path, "r", encoding="utf-8") as f:
            raw = yaml.safe_load(f) or {}

        # 1. Parsing Sezione Input & Discovery
        inp_sec = raw.get("input", {})
        data_type = str(inp_sec.get("data_type", "FASTQ")).upper()
        input_dir = Path(inp_sec.get("input_dir", "data/raw_fastq"))
        pattern_r1 = str(inp_sec.get("pattern_r1", "*_R1*.fastq.gz"))
        pattern_r2 = str(inp_sec.get("pattern_r2", "*_R2*.fastq.gz"))
        single_end = bool(inp_sec.get("single_end", False))
        allow_missing = bool(inp_sec.get("allow_missing_mates", False))

        # 2. Parsing Sezione Integrità
        integ_sec = raw.get("integrity", {})
        integrity_mode = str(integ_sec.get("mode", "fast")).lower()

        # 3. Parsing Sezione Reference
        ref_sec = raw.get("reference", {})
        reference_fasta = Path(ref_sec.get("fasta", "reference/rCRS.fasta"))
        mt_contig = str(ref_sec.get("mt_contig", "NC_012920.1"))
        auto_index = bool(ref_sec.get("auto_index", True))

        # 4. Parsing Sezione Preprocessing (con fallback su adapter_removal)
        prep_sec = raw.get("preprocessing", raw.get("adapter_removal", {}))
        filt_sec = raw.get("filtering", {})
        
        min_read_length = int(prep_sec.get("min_read_length", filt_sec.get("min_read_length", 30)))
        min_merge_overlap = int(prep_sec.get("min_merge_overlap", 11))
        collapse_reads = bool(prep_sec.get("collapse_reads", True))
        trim_ns = bool(prep_sec.get("trim_ns", True))

        # 5. Parsing Sezione Filtering & Mapping
        min_mapq = int(filt_sec.get("min_mapq", 30))
        min_baseq = int(filt_sec.get("min_baseq", 30))
        remove_duplicates = bool(filt_sec.get("remove_duplicates", True))
        unique_mapping = bool(filt_sec.get("unique_mapping_only", True))

        # 6. Parsing Sezione Consensus & Heteroplasmy
        cons_sec = raw.get("consensus", {})
        het_sec = raw.get("heteroplasmy", {})
        min_depth = int(cons_sec.get("min_depth", 10))
        het_min_freq = float(het_sec.get("min_frequency", 0.20))
        het_max_freq = float(het_sec.get("max_frequency", 0.80))

        # 7. Parsing Sezione Resources, Validation & Output
        res_sec = raw.get("resources", {})
        val_sec = raw.get("validation", {})
        out_sec = raw.get("output", {})
        
        threads = int(res_sec.get("threads", 4))
        strict_mode = bool(val_sec.get("strict_mode", False))
        expected_len = int(val_sec.get("expected_genome_length", 16569))
        results_dir = Path(out_sec.get("results_dir", "results"))
        project_name = str(raw.get("project", {}).get("name", "ancient-modern-mtdna-pipeline"))

        # Validazione Semantica (Fail-Early)
        if not (0.0 < het_min_freq < het_max_freq <= 1.0):
            raise ValueError(
                f"Invalid heteroplasmy frequencies: min={het_min_freq}, max={het_max_freq}. "
                "Must satisfy 0.0 < min_frequency < max_frequency <= 1.0"
            )

        if min_read_length < 15:
            raise ValueError(f"min_read_length too short: {min_read_length} bp (minimum allowable: 15 bp).")

        return cls(
            raw_data=raw,
            project_name=project_name,
            threads=threads,
            results_dir=results_dir,
            data_type=data_type,
            input_dir=input_dir,
            pattern_r1=pattern_r1,
            pattern_r2=pattern_r2,
            single_end=single_end,
            allow_missing_mates=allow_missing,
            integrity_mode=integrity_mode,
            reference_fasta=reference_fasta,
            mt_contig=mt_contig,
            auto_index=auto_index,
            min_read_length=min_read_length,
            min_merge_overlap=min_merge_overlap,
            collapse_reads=collapse_reads,
            trim_ns=trim_ns,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
            remove_duplicates=remove_duplicates,
            unique_mapping_only=unique_mapping,
            min_depth=min_depth,
            het_min_freq=het_min_freq,
            het_max_freq=het_max_freq,
            strict_mode=strict_mode,
            expected_genome_length=expected_len,
        )
