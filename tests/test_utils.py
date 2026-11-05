"""
Comprehensive unit test suite for the modernized src.utils package.
"""

import logging
from pathlib import Path
import tempfile
import pytest

from src.utils import (
    DEFAULT_CONFIG_PATH,
    PipelineConfig,
    load_pipeline_config,
    resolve_config_path,
    EnvironmentValidator,
    setup_logger,
    get_logger,
    configure_root_logger,
    CustomColorFormatter,
    GradientProgressBar,
    SubstageTracker,
    run_tracked_subprocess,
    get_gradient_color,
    strip_ansi,
    supports_color,
)


def test_resolve_default_config_path():
    """Verifies that resolve_config_path resolves default aliases to config/default_config.yaml."""
    resolved_none = resolve_config_path(None)
    assert resolved_none == DEFAULT_CONFIG_PATH
    assert resolved_none.is_file()

    resolved_alias = resolve_config_path("default_config")
    assert resolved_alias == DEFAULT_CONFIG_PATH

    resolved_file = resolve_config_path("default_config.yaml")
    assert resolved_file == DEFAULT_CONFIG_PATH

    with pytest.raises(FileNotFoundError):
        resolve_config_path("non_existent_config_file_xyz_123.yaml")


def test_pipeline_config_loads_all_blocks():
    """Verifies that PipelineConfig parses default_config.yaml with modern block mappings."""
    cfg = load_pipeline_config()

    # Core project & resources
    assert cfg.project_name == "ancient-modern-mtdna-pipeline"
    assert cfg.threads >= 1
    assert cfg.cpu_budget_fraction > 0.0

    # Input & reference
    assert cfg.data_type in ("FASTQ", "BAM")
    assert cfg.expected_contig == "NC_012920.1"
    assert cfg.expected_length == 16569
    assert cfg.reference_fasta.is_file() or str(cfg.reference_fasta).endswith("rCRS.fasta")

    # Block 2: Preprocessing
    assert cfg.preprocessing_tool == "AdapterRemoval"
    assert cfg.min_read_length >= 15
    assert cfg.min_merge_overlap >= 10

    # Block 3: Alignment & Filtering
    assert cfg.aligner == "bwa-aln"
    assert cfg.min_mapq >= 20
    assert cfg.min_baseq >= 20

    # Block 4: Damage
    assert cfg.damage_policy_id == "adna_standard_v1"
    assert cfg.library_treatment in ("NONE", "HALF", "FULL")

    # Block 5: Contamination
    assert cfg.contamination_method == "SCHMUTZI"
    assert cfg.schmutzi_max_iterations >= 1

    # Block 6: Consensus & Heteroplasmy
    assert cfg.min_consensus_depth >= 5
    assert 0.0 <= cfg.het_min_freq < cfg.het_max_freq <= 1.0

    # Block 7: Haplogroup
    assert cfg.haplogroup_tool == "haplogrep3"
    assert cfg.min_quality_score > 0.0


def test_pipeline_config_semantic_validations():
    """Verifies fail-early checks for erroneous parameters."""
    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as tmp:
        # Invalid heteroplasmy frequency (min > max)
        tmp.write(
            """
            input: { data_type: "FASTQ" }
            consensus: { heteroplasmy: { min_vaf: 0.90, max_vaf: 0.20 } }
            """
        )
        p1 = Path(tmp.name)

    try:
        with pytest.raises(ValueError, match="Invalid heteroplasmy frequencies"):
            PipelineConfig.from_yaml(p1)
    finally:
        p1.unlink(missing_ok=True)

    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as tmp:
        # Invalid min_read_length (< 15)
        tmp.write(
            """
            input: { data_type: "FASTQ" }
            preprocessing: { min_read_length: 8 }
            consensus: { heteroplasmy: { min_vaf: 0.10, max_vaf: 0.80 } }
            """
        )
        p2 = Path(tmp.name)

    try:
        with pytest.raises(ValueError, match="min_read_length too short"):
            PipelineConfig.from_yaml(p2)
    finally:
        p2.unlink(missing_ok=True)

    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as tmp:
        # Invalid data_type
        tmp.write(
            """
            input: { data_type: "CRAM_INVALID" }
            consensus: { heteroplasmy: { min_vaf: 0.10, max_vaf: 0.80 } }
            """
        )
        p3 = Path(tmp.name)

    try:
        with pytest.raises(ValueError, match="Invalid input data_type"):
            PipelineConfig.from_yaml(p3)
    finally:
        p3.unlink(missing_ok=True)


def test_environment_validator_tools():
    """Verifies EnvironmentValidator toolsets and inspection API."""
    assert "fastp" not in EnvironmentValidator.REQUIRED_CORE_TOOLS
    assert "SeqPrep" not in EnvironmentValidator.REQUIRED_CORE_TOOLS
    assert "samtools" in EnvironmentValidator.REQUIRED_CORE_TOOLS
    assert "bwa" in EnvironmentValidator.REQUIRED_CORE_TOOLS
    assert "AdapterRemoval" in EnvironmentValidator.ANCIENT_DNA_TOOLS

    toolset_map = EnvironmentValidator.inspect_toolset()
    assert isinstance(toolset_map, dict)
    assert "samtools" in toolset_map
    assert "bwa" in toolset_map


def test_logger_formatting_and_convenience():
    """Verifies CustomColorFormatter and logger creation helpers."""
    formatter = CustomColorFormatter()
    record_info = logging.LogRecord("test", logging.INFO, "test.py", 1, "Informational message", (), None)
    rendered_info = formatter.format(record_info)
    assert "Informational message" in rendered_info

    # Test warning highlighting on [W:: prefix
    record_w = logging.LogRecord("test", logging.INFO, "test.py", 1, "[W::hts_idx_load3] Old index", (), None)
    rendered_w = formatter.format(record_w)
    assert CustomColorFormatter.YELLOW in rendered_w

    # Test error highlighting on [E:: prefix
    record_e = logging.LogRecord("test", logging.INFO, "test.py", 1, "[E::bwa] Index corrupted", (), None)
    rendered_e = formatter.format(record_e)
    assert CustomColorFormatter.BOLD_RED in rendered_e

    log = get_logger("TestLoggerDomain")
    assert isinstance(log, logging.Logger)
