"""
Block 8 Smoke Test Suite: Validates chain-of-custody, fail-closed contract ingestion,
biological gating bypasses, and relational mutation catalog integrity (noise gating).
"""

import gzip
import hashlib
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from src.modules.synthesis.contract_aggregator import ContractIngestionError
from src.modules.synthesis.models import WarningSeverity
from src.modules.synthesis.orchestrator import SynthesisOrchestrator


def _compute_sha256_file(filepath: Path) -> str:
    """Calcola l'hash crittografico reale leggendo i byte fisici dal disco."""
    hasher = hashlib.sha256()
    with open(filepath, "rb") as f:
        hasher.update(f.read())
    return hasher.hexdigest()


class DummyRefInfo:
    """Mock unificato del contratto ReferenceInfo per mantenere coerenza sul contig."""
    def __init__(self, length_bp=16569, contig="chrM"):
        self.length_bp = length_bp
        self.contig = contig


@pytest.fixture
def workspace(tmp_path: Path):
    """Generates a structurally sound Block 1-7 artifacts filesystem reflecting TRUE production contracts."""
    sample_id = "SMOKE_01"
    results_dir = tmp_path / "results"

    for d in ["preprocessing", "alignment", "damage", "contamination", "consensus", "haplogroup", "master_summary"]:
        (results_dir / d / sample_id).mkdir(parents=True, exist_ok=True)

    bam_path = results_dir / "damage" / sample_id / f"{sample_id}.rescaled.bam"
    bam_path.write_text("DUMMY_BAM_BINARY", encoding="utf-8")

    bai_path = results_dir / "damage" / sample_id / f"{sample_id}.rescaled.bai"
    bai_path.write_text("DUMMY_BAI_BINARY", encoding="utf-8")

    fasta_path = results_dir / "consensus" / sample_id / f"{sample_id}.consensus.fasta"
    fasta_content = ">SMOKE_01\n" + ("A" * 16569) + "\n"
    fasta_path.write_text(fasta_content, encoding="utf-8")

    indel_fasta_path = results_dir / "consensus" / sample_id / f"{sample_id}.consensus_indel.fasta"
    indel_fasta_path.write_text(fasta_content, encoding="utf-8")

    vcf_path = results_dir / "consensus" / sample_id / f"{sample_id}.variants.vcf.gz"
    vcf_content = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chrM,length=16569>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSMOKE_01\n"
        "chrM\t73\t.\tA\tG\t.\tPASS\tDP=100;AF=0.99\tGT:DP:AD\t1/1:100:0,100\n"
    )
    with gzip.open(vcf_path, "wt", encoding="utf-8") as f:
        f.write(vcf_content)

    tbi_path = Path(f"{vcf_path}.tbi")
    tbi_path.write_text("DUMMY_TABIX_BINARY", encoding="utf-8")

    qc_json_path = results_dir / "consensus" / sample_id / f"{sample_id}.consensus_qc.json"
    qc_json_content = {
        "heteroplasmy_sites_pass": [73],
        "heteroplasmy_details": {
            "73": {"ref": "A", "alt": "G", "vaf": 0.99, "depth": 100, "filter_status": "PASS"}
        }
    }
    qc_json_path.write_text(json.dumps(qc_json_content), encoding="utf-8")

    b2_contract = {"sample_id": sample_id, "raw_reads_total": 10000}
    (results_dir / "preprocessing" / sample_id / f"{sample_id}.block2_contract.json").write_text(json.dumps(b2_contract))

    b3_contract = {"sample_id": sample_id, "mapped_reads_mt": 5000}
    (results_dir / "alignment" / sample_id / f"{sample_id}.block3_contract.json").write_text(json.dumps(b3_contract))

    b4_contract = {
        "sample_id": sample_id,
        "authentication_status": "AUTHENTIC_ANCIENT_SUPPORTED",
        "rescaling_performed": True,
        "rescaled_bam": str(bam_path),
        "rescaled_bam_sha256": _compute_sha256_file(bam_path)
    }
    (results_dir / "damage" / sample_id / f"{sample_id}.block4_contract.json").write_text(json.dumps(b4_contract))

    b5_contract = {
        "sample_id": sample_id,
        "contamination_fraction": 0.015,
        "recommended_action": "STANDARD_CALLING"
    }
    (results_dir / "contamination" / sample_id / f"{sample_id}.block5_contract.json").write_text(json.dumps(b5_contract))

    b6_path = results_dir / "consensus" / sample_id / f"{sample_id}.block6_contract.json"
    b6_contract = {
        "sample_id": sample_id,
        "variants_vcf_path": str(vcf_path),
        "consensus_fasta_sha256": _compute_sha256_file(fasta_path),
        "variants_vcf_sha256": _compute_sha256_file(vcf_path)
    }
    b6_content = json.dumps(b6_contract)
    b6_path.write_text(b6_content)

    b7_contract = {
        "sample_id": sample_id,
        "primary_haplogroup": "H1",
        "fasta_sha256": _compute_sha256_file(fasta_path),
        "variants_vcf_sha256": _compute_sha256_file(vcf_path),
        "block6_contract_sha256": _compute_sha256_file(b6_path),
        "found_mutations": ["73G", "16519C"],
        "overall_validation_status": "PASS",
        "haplogroup_status": "ASSIGNED_CONFIDENT"
    }
    (results_dir / "haplogroup" / sample_id / f"{sample_id}.block7_contract.json").write_text(json.dumps(b7_contract))

    return results_dir, sample_id


@patch("src.modules.synthesis.provenance_checker.pysam", None)
@patch("src.modules.synthesis.mutation_cataloguer.pysam", None)
@patch("src.modules.synthesis.orchestrator.pysam", None)
def test_nominal_synthesis_execution(workspace):
    """Verifies flawless chain-of-custody and tidy catalog generation."""
    results_dir, sample_id = workspace

    orchestrator = SynthesisOrchestrator(
        results_dir=results_dir,
        reference_info=DummyRefInfo(),
    )

    record, report, artifacts = orchestrator.process_sample(sample_id)

    assert record.overall_status == "PASS"
    assert record.authenticity.authentication_status == "AUTHENTIC_ANCIENT_SUPPORTED"
    assert record.haplogroup.primary_haplogroup == "H1"

    assert report.overall_status == "PASS"
    for check in report.artifact_checks:
        if check.artifact_name in {"rescaled_ancient_bam", "canonical_consensus_fasta", "variants_vcf"}:
            assert check.status.value == "VERIFIED"

    assert (results_dir / "master_summary" / f"{sample_id}.master.json").exists()
    assert (results_dir / "mutations" / "single_samples" / f"{sample_id}.mutation_catalog.tsv").exists()


@patch("src.modules.synthesis.provenance_checker.pysam", None)
@patch("src.modules.synthesis.mutation_cataloguer.pysam", None)
@patch("src.modules.synthesis.orchestrator.pysam", None)
def test_cryptographic_mismatch_escalation(workspace):
    """Injects a physical file modification to test chain-of-custody collapse."""
    results_dir, sample_id = workspace

    fasta_path = results_dir / "consensus" / sample_id / f"{sample_id}.consensus.fasta"
    fasta_path.write_text(">SMOKE_01\n" + ("C" * 16569) + "\n", encoding="utf-8")

    orchestrator = SynthesisOrchestrator(results_dir=results_dir, reference_info=DummyRefInfo())

    contracts = orchestrator.aggregator.locate_contracts(sample_id)
    record = orchestrator.aggregator.aggregate_sample(sample_id, contracts=contracts)
    report = orchestrator.checker.audit_sample_provenance(record, contracts)

    assert record.overall_status == "FAIL"
    assert report.overall_status == "FAIL"

    mismatch_flags = [w for w in record.warnings if w.code == "PROVENANCE_HASH_MISMATCH"]
    assert len(mismatch_flags) > 0
    assert mismatch_flags[0].severity == WarningSeverity.CRITICAL


def test_contract_aggregator_corrupted_json(workspace):
    """Simulates I/O truncation or disk write failure on a mandatory contract."""
    results_dir, sample_id = workspace

    b5_path = results_dir / "contamination" / sample_id / f"{sample_id}.block5_contract.json"
    b5_path.write_text("{ \"sample_id\": \"SMOKE_01\", \"contamination_fract", encoding="utf-8")

    orchestrator = SynthesisOrchestrator(results_dir=results_dir, reference_info=DummyRefInfo())

    with pytest.raises(ContractIngestionError) as exc_info:
        orchestrator.process_sample(sample_id)
    assert "Corrupted or truncated JSON" in str(exc_info.value)


@patch("src.modules.synthesis.provenance_checker.pysam", None)
@patch("src.modules.synthesis.mutation_cataloguer.pysam", None)
@patch("src.modules.synthesis.orchestrator.pysam", None)
def test_biological_gating_bypass(workspace):
    """
    Tests pessimistic gating: heavily contaminated sample legitimately skips processing.
    """
    results_dir, sample_id = workspace

    b5_path = results_dir / "contamination" / sample_id / f"{sample_id}.block5_contract.json"
    b5_contract = {
        "sample_id": sample_id,
        "contamination_fraction": 0.45,
        "recommended_action": "REJECT_CONTAMINATED"
    }
    b5_path.write_text(json.dumps(b5_contract))

    b7_path = results_dir / "haplogroup" / sample_id / f"{sample_id}.block7_contract.json"
    b7_data = json.loads(b7_path.read_text(encoding="utf-8"))
    b7_data["overall_validation_status"] = "FAIL"
    b7_data["haplogroup_status"] = "REJECTED"
    b7_path.write_text(json.dumps(b7_data), encoding="utf-8")

    (results_dir / "consensus" / sample_id / f"{sample_id}.consensus.fasta").unlink()
    (results_dir / "consensus" / sample_id / f"{sample_id}.consensus_indel.fasta").unlink()
    vcf_path = results_dir / "consensus" / sample_id / f"{sample_id}.variants.vcf.gz"
    vcf_path.unlink()
    Path(f"{str(vcf_path)}.tbi").unlink()

    orchestrator = SynthesisOrchestrator(results_dir=results_dir, reference_info=DummyRefInfo())
    record, report, _ = orchestrator.process_sample(sample_id)

    assert record.overall_status == "FAIL"

    vcf_check = next((c for c in report.artifact_checks if c.artifact_name == "variants_vcf"), None)
    assert vcf_check is not None
    assert vcf_check.status.value == "GATED_OMITTED"

    assert len(record.mutations) == 0


@patch("src.modules.synthesis.provenance_checker.pysam", None)
@patch("src.modules.synthesis.mutation_cataloguer.pysam", None)
@patch("src.modules.synthesis.orchestrator.pysam", None)
def test_noise_gating_sub_threshold_variants(workspace):
    """
    Validates the MutationCataloguer's noise gate: variants < 10% VAF must be purged.
    """
    results_dir, sample_id = workspace

    vcf_path = results_dir / "consensus" / sample_id / f"{sample_id}.variants.vcf.gz"

    with gzip.open(vcf_path, "rt", encoding="utf-8") as f:
        vcf_content = f.read()
    vcf_content += "chrM\t1000\t.\tC\tT\t.\tFILTER_VAF_RANGE\tDP=100;AF=0.04\tGT:DP:AD\t0/1:100:96,4\n"
    with gzip.open(vcf_path, "wt", encoding="utf-8") as f:
        f.write(vcf_content)

    qc_json_path = results_dir / "consensus" / sample_id / f"{sample_id}.consensus_qc.json"
    qc_data = json.loads(qc_json_path.read_text(encoding="utf-8"))
    qc_data["heteroplasmy_details"]["1000"] = {
        "ref": "C", "alt": "T", "vaf": 0.04, "depth": 100, "filter_status": "FILTER_VAF_RANGE"
    }
    qc_json_path.write_text(json.dumps(qc_data), encoding="utf-8")

    b6_path = results_dir / "consensus" / sample_id / f"{sample_id}.block6_contract.json"
    b6_data = json.loads(b6_path.read_text(encoding="utf-8"))
    b6_data["variants_vcf_sha256"] = _compute_sha256_file(vcf_path)
    b6_content = json.dumps(b6_data)
    b6_path.write_text(b6_content, encoding="utf-8")

    b7_path = results_dir / "haplogroup" / sample_id / f"{sample_id}.block7_contract.json"
    b7_data = json.loads(b7_path.read_text(encoding="utf-8"))
    b7_data["variants_vcf_sha256"] = _compute_sha256_file(vcf_path)
    b7_data["block6_contract_sha256"] = _compute_sha256_file(b6_path)
    b7_path.write_text(json.dumps(b7_data), encoding="utf-8")

    orchestrator = SynthesisOrchestrator(
        results_dir=results_dir,
        reference_info=DummyRefInfo(),
        config={"synthesis": {"min_vaf_reportable": 0.10}}
    )
    record, _, _ = orchestrator.process_sample(sample_id)

    assert record.overall_status == "PASS"

    mut_positions = [m.position for m in record.mutations]
    assert 73 in mut_positions
    assert 1000 not in mut_positions