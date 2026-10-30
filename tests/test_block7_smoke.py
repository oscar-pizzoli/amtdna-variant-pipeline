"""
Block 7 Smoke Test: Complete integration and validation test suite.
"""

import hashlib
import inspect
import json
import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from src.modules.alignment.models import ReferenceInfo
from src.modules.haplogrouping.cross_block_auditor import (
    CrossBlockAuditor,
    PhylotreeGraph,
)
from src.modules.haplogrouping.enums import (
    AuditCode,
    HaplogroupStatus,
    OverallValidationStatus,
    PrivateMutationDamageClass,
)
from src.modules.haplogrouping.haplogrep_engine import (
    HaplogrepEngine,
    HaplogrepExecutionError,
)
from src.modules.haplogrouping.models import HaplogrepRunConfig
from src.modules.haplogrouping.orchestrator import HaplogroupOrchestrator
from src.modules.haplogrouping.sequence_validator import SequenceValidator


def _create_mock_rcrs_sequence(length: int = 16569) -> str:
    pattern = "ATCGTACGATCGTACGATCG"
    repeats = (length // len(pattern)) + 1
    bases = list((pattern * repeats)[:length])

    bases[73 - 1] = "A"       # rCRS A73
    bases[263 - 1] = "A"      # rCRS A263
    bases[315 - 1] = "C"      # rCRS C315
    bases[524 - 1] = "C"      # rCRS C524
    bases[3010 - 1] = "G"     # rCRS G3010 (H1 diagnostic transition: G3010A)
    bases[16069 - 1] = "C"    # rCRS C16069 (C16069T transition)

    return "".join(bases)


def _write_mock_fasta(path: Path, sample_id: str, sequence: str) -> str:
    content = f">{sample_id}\n{sequence}\n"
    path.write_text(content, encoding="utf-8")
    return hashlib.sha256(content.encode("utf-8")).hexdigest()


def _mock_subprocess_run(cmd, *args, **kwargs):
    mock_res = MagicMock()
    mock_res.returncode = 0

    if "--version" in cmd:
        mock_res.stdout = "haplogrep3 3.2.1\n"
        mock_res.stderr = ""
        return mock_res

    if "classify" in cmd:
        out_idx = cmd.index("--out") + 1
        out_path = Path(cmd[out_idx])
        in_path = Path(cmd[cmd.index("--in") + 1])
        in_name = in_path.name.lower()

        is_fasta = any(in_name.endswith(ext) for ext in [".fasta", ".fa", ".fasta.gz", ".fa.gz"])
        is_vcf = any(in_name.endswith(ext) for ext in [".vcf", ".vcf.gz", ".bcf"])

        if is_fasta:
            payload = {
                "samples": [
                    {
                        "sample": "SAMPLE_TEST_01",
                        "haplogroup": "H",
                        "quality": 0.95,
                        "hits": [
                            {
                                "haplogroup": "H",
                                "quality": 0.95,
                                "expected": ["73A", "263A"],
                                "found": ["73A", "263A"],
                                "missing": [],
                                "private": ["315.1C", "524D", "73G", "C16069T", "A263C"],
                            },
                            {
                                "haplogroup": "H1",
                                "quality": 0.88,
                            },
                        ],
                    }
                ]
            }
        elif is_vcf:
            # Subclade H1A differentiated naturally by 3010A
            payload = {
                "samples": [
                    {
                        "sample": "SAMPLE_TEST_01",
                        "haplogroup": "H1A",
                        "quality": 0.94,
                        "hits": [
                            {
                                "haplogroup": "H1A",
                                "quality": 0.94,
                                "expected": ["73A", "263A", "3010A"],
                                "found": ["73A", "263A", "3010A"],
                                "missing": [],
                                "private": [],
                            }
                        ],
                    }
                ]
            }
        else:
            raise ValueError(f"Unrecognized paleogenomic file extension: {in_name}")

        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(payload, f)

        mock_res.stdout = "Classification completed successfully."
        mock_res.stderr = ""
        return mock_res

    return mock_res


@pytest.fixture
def minimal_config():
    return {
        "haplogroup": {
            "tool": "haplogrep3",
            "tree_id": "phylotree-fu-rcrs",
            "tree_version": "1.2",
            "metric": "kulczynski",
            "min_quality_score": 0.80,
            "cross_check_vcf": True,
            "timeout_seconds": 30,
            "retain_raw_output": True,
        },
        "validation": {
            "sequence_type": "CANONICAL",
            "max_n_percentage": 15.0,
        },
        "cross_block_audit": {
            "min_quality_score": 0.80,
            "max_missing_diagnostic_mutations": 2,
            "max_private_mutations": 10,
            "contam_affinity_threshold": 0.015,
            "subclade_contam_max_fraction": 0.03,
            "evaluate_contaminant_closeness": True,
            "flag_damage_in_private_mutations": True,
        },
        "reference": {
            "fasta_path": "data/reference/rCRS.fasta",
            "expected_length": 16569,
        },
    }


def test_phylotree_graph_topological_resolution():
    assert PhylotreeGraph.is_descendant("H", "HV") is True
    assert PhylotreeGraph.is_descendant("HV", "H") is False
    assert PhylotreeGraph.is_descendant("H20", "H2") is False
    assert PhylotreeGraph.is_descendant("H2A", "H2") is True
    assert PhylotreeGraph.is_descendant("H2A1", "H2") is True
    assert PhylotreeGraph.is_descendant("K", "U") is True
    assert PhylotreeGraph.is_descendant("U", "K") is False


def test_mutation_parsing_and_damage_triage():
    auditor = CrossBlockAuditor()
    rcrs_seq = _create_mock_rcrs_sequence()

    parsed_ins = auditor._parse_single_mutation("315.1C", rcrs_seq)
    assert parsed_ins.variant_type == "INS"
    assert parsed_ins.position == 315

    parsed_del = auditor._parse_single_mutation("524D", rcrs_seq)
    assert parsed_del.variant_type == "DEL"
    assert parsed_del.ref_base == "C"

    parsed_snv = auditor._parse_single_mutation("73G", rcrs_seq)
    assert parsed_snv.variant_type == "SNV"
    assert parsed_snv.ref_base == "A"
    assert parsed_snv.alt_base == "G"

    assessments = auditor._audit_private_mutations_chemistry(
        ["315.1C", "524D", "73G", "C16069T", "A263C"],
        {"authentication_status": "AUTHENTIC_ANCIENT_SUPPORTED"},
        rcrs_seq,
    )
    res_map = {a.mutation: a.classification for a in assessments}

    assert res_map["315.1C"] == PrivateMutationDamageClass.UNASSESSED
    assert res_map["524D"] == PrivateMutationDamageClass.UNASSESSED
    assert res_map["73G"] == PrivateMutationDamageClass.TRANSITION_NON_DAMAGE
    assert res_map["C16069T"] == PrivateMutationDamageClass.TRANSITION_REQUIRES_SPATIAL_AUDIT
    assert res_map["A263C"] == PrivateMutationDamageClass.TRANSVERSION_NON_DAMAGE


def test_mutation_parsing_out_of_bounds_safety():
    auditor = CrossBlockAuditor()
    short_ref = "ACGT" * 25

    parsed_ins_oob = auditor._parse_single_mutation("315.1C", short_ref)
    assert parsed_ins_oob.variant_type == "INS"
    assert parsed_ins_oob.position == 315

    parsed_snv_oob = auditor._parse_single_mutation("16069T", short_ref)
    assert parsed_snv_oob.variant_type == "SNV"
    assert parsed_snv_oob.ref_base == "N"

    assessments = auditor._audit_private_mutations_chemistry(["16069T"], {}, short_ref)
    assert assessments[0].classification == PrivateMutationDamageClass.UNASSESSED


def test_reference_concatenation_bomb_circuit_breaker(tmp_path, minimal_config):
    multi_ref = tmp_path / "fake_hg38.fasta"
    multi_ref.write_text(">chr1\nACGT\n>chr2\nACGT\n", encoding="utf-8")
    seq, headers, triggered = SequenceValidator.stream_fasta_sequence(multi_ref)
    assert seq is None
    assert triggered is True


@patch("shutil.which", return_value="/usr/local/bin/haplogrep3")
@patch("subprocess.run", side_effect=_mock_subprocess_run)
def test_block7_successful_orchestration_and_cli_contracts(mock_subp, mock_which, tmp_path, minimal_config):
    """
    Nominal test verifying that transition 3010A with spatial_damage_suspect=False
    successfully promotes the subclade to H1A.
    """
    sample_id = "SAMPLE_TEST_01"
    sample_dir = tmp_path / sample_id
    fasta_path = tmp_path / f"{sample_id}.consensus.fasta"
    vcf_path = tmp_path / f"{sample_id}.variants.vcf.gz"
    vcf_path.touch()

    rcrs_sequence = _create_mock_rcrs_sequence()
    fasta_sha256 = _write_mock_fasta(fasta_path, sample_id=sample_id, sequence=rcrs_sequence)

    ref_fasta = tmp_path / "rCRS.fasta"
    _write_mock_fasta(ref_fasta, sample_id="NC_012920.1", sequence=rcrs_sequence)

    b4_path = tmp_path / "block4_contract.json"
    b4_path.write_text(
        json.dumps({
            "sample_id": sample_id,
            "authentication_status": "AUTHENTIC_ANCIENT_SUPPORTED",
            "rescaling_performed": True,
        }),
        encoding="utf-8",
    )

    b5_path = tmp_path / "block5_contract.json"
    b5_path.write_text(
        json.dumps({
            "sample_id": sample_id,
            "contamination_fraction": 0.01,
            "recommended_action": "STANDARD_CALLING",
            "contaminant_haplogroups": ["T2B"],
        }),
        encoding="utf-8",
    )

    # Block 6 exports heteroplasmy_details verifying 3010A has spatial_damage_suspect=False
    b6_path = tmp_path / "block6_contract.json"
    b6_path.write_text(
        json.dumps({
            "sample_id": sample_id,
            "consensus_fasta_sha256": fasta_sha256,
            "heteroplasmy_sites_pass": [3010],
            "heteroplasmy_details": {
                "3010": {
                    "ref": "G",
                    "alt": "A",
                    "vaf": 0.35,
                    "internal_support_count": 4,
                    "terminal_support_count": 0,
                    "spatial_damage_suspect": False,
                }
            },
        }),
        encoding="utf-8",
    )

    sig = inspect.signature(ReferenceInfo)
    ref_kwargs = {}
    if "contig_name" in sig.parameters:
        ref_kwargs["contig_name"] = "NC_012920.1"
    if "length_bp" in sig.parameters:
        ref_kwargs["length_bp"] = 16569
    if "fasta_path" in sig.parameters:
        ref_kwargs["fasta_path"] = ref_fasta
    if "is_circular" in sig.parameters:
        ref_kwargs["is_circular"] = True
    ref_info = ReferenceInfo(**ref_kwargs)

    orchestrator = HaplogroupOrchestrator(config=minimal_config)
    contract = orchestrator.run_sample(
        sample_id=sample_id,
        consensus_fasta=fasta_path,
        sample_output_dir=sample_dir,
        variants_vcf=vcf_path,
        reference_info=ref_info,
        reference_fasta=ref_fasta,
        expected_fasta_sha256=fasta_sha256,
        block4_contract_path=b4_path,
        block5_contract_path=b5_path,
        block6_contract_path=b6_path,
    )

    classify_calls = [call for call in mock_subp.call_args_list if "classify" in call[0][0]]
    assert len(classify_calls) == 2

    for call in classify_calls:
        cmd = call[0][0]
        kwargs = call[1]
        assert cmd[0] == "/usr/local/bin/haplogrep3"
        assert cmd[1] == "classify"
        assert "--tree" in cmd
        assert "--metric" in cmd
        assert kwargs.get("timeout") == 30

    assert contract.sample_id == sample_id
    assert contract.primary_haplogroup == "H"
    assert contract.quality_score == 0.95
    assert contract.subclade_shift_detected is True
    assert contract.subclade_shift_supported is True
    assert contract.resolved_haplogroup == "H1A"
    assert contract.overall_validation_status in {OverallValidationStatus.PASS, OverallValidationStatus.WARN}


@patch("shutil.which", return_value="/usr/local/bin/haplogrep3")
def test_block7_subclade_rejection_pure_deamination(mock_which, tmp_path, minimal_config):
    """
    Verify that if 3010A has spatial_damage_suspect=True, the subclade shift is rejected.
    """
    def _mock_subp_pure_damage(cmd, *args, **kwargs):
        mock_res = MagicMock(returncode=0)
        if "--version" in cmd:
            mock_res.stdout = "haplogrep3 3.2.1\n"
            return mock_res
        if "classify" in cmd:
            out_idx = cmd.index("--out") + 1
            out_path = Path(cmd[out_idx])
            in_path = Path(cmd[cmd.index("--in") + 1])
            in_name = in_path.name.lower()

            if in_name.endswith((".fasta", ".fa")):
                payload = {
                    "samples": [{
                        "sample": "SAMPLE_DAMAGE_ONLY",
                        "haplogroup": "H",
                        "quality": 0.95,
                        "hits": [{"haplogroup": "H", "quality": 0.95, "found": ["73A", "263A"]}],
                    }]
                }
            else:
                payload = {
                    "samples": [{
                        "sample": "SAMPLE_DAMAGE_ONLY",
                        "haplogroup": "H1A",
                        "quality": 0.94,
                        "hits": [{"haplogroup": "H1A", "quality": 0.94, "found": ["73A", "263A", "3010A"]}],
                    }]
                }
            out_path.parent.mkdir(parents=True, exist_ok=True)
            with open(out_path, "w", encoding="utf-8") as f:
                json.dump(payload, f)
            return mock_res
        return mock_res

    with patch("subprocess.run", side_effect=_mock_subp_pure_damage):
        sample_id = "SAMPLE_DAMAGE_ONLY"
        sample_dir = tmp_path / sample_id
        fasta_path = tmp_path / f"{sample_id}.consensus.fasta"
        vcf_path = tmp_path / f"{sample_id}.variants.vcf.gz"
        vcf_path.touch()

        rcrs_sequence = _create_mock_rcrs_sequence()
        fasta_sha256 = _write_mock_fasta(fasta_path, sample_id=sample_id, sequence=rcrs_sequence)

        ref_fasta = tmp_path / "rCRS.fasta"
        _write_mock_fasta(ref_fasta, sample_id="NC_012920.1", sequence=rcrs_sequence)

        b4_path = tmp_path / "b4.json"
        b4_path.write_text(json.dumps({"sample_id": sample_id, "authentication_status": "AUTHENTIC_ANCIENT_SUPPORTED"}))

        b5_path = tmp_path / "b5.json"
        b5_path.write_text(json.dumps({"sample_id": sample_id, "contamination_fraction": 0.01, "recommended_action": "STANDARD_CALLING"}))

        # Block 6 explicitly flags 3010 as spatial_damage_suspect=True
        b6_path = tmp_path / "b6.json"
        b6_path.write_text(json.dumps({
            "sample_id": sample_id,
            "consensus_fasta_sha256": fasta_sha256,
            "heteroplasmy_sites_pass": [3010],
            "heteroplasmy_details": {
                "3010": {
                    "ref": "G",
                    "alt": "A",
                    "vaf": 0.35,
                    "internal_support_count": 0,
                    "terminal_support_count": 4,
                    "spatial_damage_suspect": True,
                }
            },
        }))

        orchestrator = HaplogroupOrchestrator(config=minimal_config)
        contract = orchestrator.run_sample(
            sample_id=sample_id,
            consensus_fasta=fasta_path,
            sample_output_dir=sample_dir,
            variants_vcf=vcf_path,
            reference_fasta=ref_fasta,
            expected_fasta_sha256=fasta_sha256,
            block4_contract_path=b4_path,
            block5_contract_path=b5_path,
            block6_contract_path=b6_path,
        )

        assert contract.subclade_shift_detected is True
        assert contract.subclade_shift_supported is False
        assert contract.primary_haplogroup == "H"
        assert contract.resolved_haplogroup == "H"
        assert any(f["code"] == AuditCode.SUSPECT_HETEROPLASMY_SUBCLADE_REJECTED.value for f in contract.flags)


@patch("shutil.which", return_value="/usr/local/bin/haplogrep3")
@patch("subprocess.run")
def test_block7_haplogrep_timeout_resilience(mock_subp, mock_which, tmp_path, minimal_config):
    mock_version = MagicMock(returncode=0, stdout="haplogrep3 3.2.1\n", stderr="")

    def _timeout_side_effect(cmd, *args, **kwargs):
        if "--version" in cmd:
            return mock_version
        raise subprocess.TimeoutExpired(cmd=cmd, timeout=30)

    mock_subp.side_effect = _timeout_side_effect

    sample_id = "SAMPLE_TIMEOUT"
    sample_dir = tmp_path / sample_id
    fasta_path = tmp_path / f"{sample_id}.consensus.fasta"
    _write_mock_fasta(fasta_path, sample_id=sample_id, sequence=_create_mock_rcrs_sequence())

    b4 = tmp_path / "b4.json"
    b5 = tmp_path / "b5.json"
    b6 = tmp_path / "b6.json"
    b4.write_text("{}", encoding="utf-8")
    b5.write_text("{}", encoding="utf-8")
    b6.write_text("{}", encoding="utf-8")

    orchestrator = HaplogroupOrchestrator(config=minimal_config)
    contract = orchestrator.run_sample(
        sample_id=sample_id,
        consensus_fasta=fasta_path,
        sample_output_dir=sample_dir,
        block4_contract_path=b4,
        block5_contract_path=b5,
        block6_contract_path=b6,
    )

    assert contract.primary_haplogroup == "REJECTED"
    assert contract.overall_validation_status == OverallValidationStatus.FAIL


@patch("shutil.which", return_value="/usr/local/bin/haplogrep3")
@patch("subprocess.run")
def test_block7_haplogrep_corrupt_json_resilience(mock_subp, mock_which, tmp_path, minimal_config):
    mock_version = MagicMock(returncode=0, stdout="haplogrep3 3.2.1\n", stderr="")

    def _corrupt_json_side_effect(cmd, *args, **kwargs):
        if "--version" in cmd:
            return mock_version
        out_idx = cmd.index("--out") + 1
        out_path = Path(cmd[out_idx])
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text('{"samples": [ { "incomplete": true ', encoding="utf-8")
        return MagicMock(returncode=0, stdout="", stderr="")

    mock_subp.side_effect = _corrupt_json_side_effect

    sample_id = "SAMPLE_CORRUPT"
    sample_dir = tmp_path / sample_id
    fasta_path = tmp_path / f"{sample_id}.consensus.fasta"
    _write_mock_fasta(fasta_path, sample_id=sample_id, sequence=_create_mock_rcrs_sequence())

    b4 = tmp_path / "b4.json"
    b5 = tmp_path / "b5.json"
    b6 = tmp_path / "b6.json"
    b4.write_text("{}", encoding="utf-8")
    b5.write_text("{}", encoding="utf-8")
    b6.write_text("{}", encoding="utf-8")

    orchestrator = HaplogroupOrchestrator(config=minimal_config)
    contract = orchestrator.run_sample(
        sample_id=sample_id,
        consensus_fasta=fasta_path,
        sample_output_dir=sample_dir,
        block4_contract_path=b4,
        block5_contract_path=b5,
        block6_contract_path=b6,
    )

    assert contract.primary_haplogroup == "REJECTED"
    assert contract.overall_validation_status == OverallValidationStatus.FAIL


@patch("shutil.which", return_value="/usr/local/bin/haplogrep3")
def test_missing_upstream_contracts_hard_stop(mock_which, tmp_path, minimal_config):
    sample_id = "SAMPLE_NO_CONTRACT"
    sample_dir = tmp_path / sample_id
    fasta_path = tmp_path / f"{sample_id}.consensus.fasta"
    _write_mock_fasta(fasta_path, sample_id=sample_id, sequence=_create_mock_rcrs_sequence())

    orchestrator = HaplogroupOrchestrator(config=minimal_config)
    contract = orchestrator.run_sample(
        sample_id=sample_id,
        consensus_fasta=fasta_path,
        sample_output_dir=sample_dir,
        block4_contract_path=tmp_path / "non_existent_b4.json",
        block5_contract_path=tmp_path / "non_existent_b5.json",
        block6_contract_path=tmp_path / "non_existent_b6.json",
    )

    assert contract.primary_haplogroup == "REJECTED"
    assert contract.overall_validation_status == OverallValidationStatus.FAIL
    assert any(f["code"] == AuditCode.MISSING_UPSTREAM_CONTRACT.value for f in contract.flags)