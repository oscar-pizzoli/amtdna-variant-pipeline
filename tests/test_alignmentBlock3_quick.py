"""Quick Integration & Verification Test for Block 3 (Alignment & Filtering)."""

import sys
from pathlib import Path
import yaml

# Inietta la project root in sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.modules.alignment.models import (
    AlignmentFilterPolicy,
    AlignerType,
    DedupMode,
    MultiMappingPolicy,
    ReadGroupInfo,
    ReferenceInfo,
)
from src.modules.alignment.orchestrator import AlignmentOrchestrator


def test_alignment_block3_quick():
    """Verifies YAML parsing, reference resolution, and Block 3 models instantiation."""
    config_path = PROJECT_ROOT / "config" / "default_config.yaml"
    assert config_path.is_file(), f"File di configurazione non trovato: {config_path}"

    with open(config_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    # 1. Risoluzione Reference FASTA
    ref_rel = cfg["reference"]["fasta_path"]
    ref_path = (PROJECT_ROOT / ref_rel).resolve() if not Path(ref_rel).is_absolute() else Path(ref_rel)
    assert ref_path.is_file(), f"Reference FASTA non trovata: {ref_path}"

    ref_info = ReferenceInfo(
        fasta_path=ref_path,
        contig_name=cfg["reference"]["expected_contig"],
        length_bp=cfg["reference"]["expected_length"],
        is_circular=cfg["reference"].get("is_circular", True),
    )
    assert ref_info.length_bp == 16569, f"Lunghezza rCRS errata: {ref_info.length_bp}"

    # 2. Istanziazione Policy di Allineamento e Filtraggio
    policy = AlignmentFilterPolicy(
        multimapping=MultiMappingPolicy(
            min_mapq=cfg["filtering"]["multimapping"]["min_mapq"],
            filter_unmapped=cfg["filtering"]["multimapping"]["filter_unmapped"],
            filter_secondary=cfg["filtering"]["multimapping"]["filter_secondary"],
            filter_supplementary=cfg["filtering"]["multimapping"]["filter_supplementary"],
        ),
        dedup_mode=DedupMode(cfg["deduplication"]["mode"]),
        min_breadth_3x_fraction=cfg["filtering"]["thresholds"]["min_breadth_3x_fraction"],
        max_duplication_rate_fraction=cfg["filtering"]["thresholds"]["max_duplication_rate_fraction"],
        min_mean_depth=cfg["filtering"]["thresholds"]["min_mean_depth"],
    )
    assert policy.dedup_mode == DedupMode.BOTH

    # 3. Validazione Read Group Header
    rg_info = ReadGroupInfo(
        read_group_id="TEST_RG",
        sample_id="OTZI",
        library_id="LIB1",
        platform=cfg["alignment"]["read_group"]["platform"],
        sequencing_center=cfg["alignment"]["read_group"]["center"],
    )
    rg_str = rg_info.to_bwa_rg_string()
    assert "@RG" in rg_str
    assert "SM:OTZI" in rg_str


if __name__ == "__main__":
    test_alignment_block3_quick()
    print("✅ Block 3 Alignment Models & Configuration Verified (OK)")