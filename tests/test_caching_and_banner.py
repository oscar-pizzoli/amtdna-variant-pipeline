"""
Unit and integration tests for hash-aware caching, sample disambiguation,
compact CLI banner, and prompt option pruning.
"""

from pathlib import Path
import sys
import tempfile
import pytest

from src.pipeline.runner import EndToEndPipelineRunner
from src.pipeline.stage_tracker import StageTracker
from src.utils import (
    calculate_sha256,
    compute_sample_input_hash,
    load_pipeline_config,
)


def test_compute_sample_input_hash(tmp_path: Path):
    """Verifies deterministic SHA-256 calculation for single-end and paired inputs."""
    f1 = tmp_path / "sample1_R1.fastq.gz"
    f2 = tmp_path / "sample1_R2.fastq.gz"
    f1.write_bytes(b"@SEQ1\nACGT\n+\nIIII\n")
    f2.write_bytes(b"@SEQ1\nTGCA\n+\nIIII\n")

    h_single = compute_sample_input_hash(f1)
    assert h_single == calculate_sha256(f1)

    h_paired = compute_sample_input_hash([f1, f2])
    assert len(h_paired) == 64
    assert h_paired != h_single

    # Determinism: order in input list doesn't affect result
    assert compute_sample_input_hash([f2, f1]) == h_paired


def test_stage_tracker_hash_isolation_for_same_sample_name(tmp_path: Path):
    """
    Verifies that two different samples sharing the identical sample_id
    maintain isolated cache status and never overwrite each other's cache validity.
    """
    cache_dir = tmp_path / "cache"

    sample_id = "EGIZIANOANTICO1"
    hash_v1 = "1111111111111111111111111111111111111111111111111111111111111111"
    hash_v2 = "2222222222222222222222222222222222222222222222222222222222222222"

    dummy_artifact = tmp_path / "trimmed.settings"
    dummy_artifact.write_text("Dummy settings")

    # 1. First run for sample V1
    tracker_v1 = StageTracker.load_or_create(cache_dir, sample_id=sample_id, sample_hash=hash_v1)
    assert tracker_v1.sample_hash == hash_v1
    assert "PREPROCESSING" not in tracker_v1.stages

    tracker_v1.mark_success(
        stage_name="PREPROCESSING",
        elapsed_sec=1.5,
        artifacts={"settings": dummy_artifact},
        composite_hash="comp_hash_v1",
    )
    assert tracker_v1.is_stage_valid("PREPROCESSING", [dummy_artifact], "comp_hash_v1") is True

    # 2. Second run with a DIFFERENT sample that shares the same sample_id (hash_v2)
    tracker_v2 = StageTracker.load_or_create(cache_dir, sample_id=sample_id, sample_hash=hash_v2)
    assert tracker_v2.sample_hash == hash_v2
    # Crucial assertion: tracker_v2 MUST NOT contain the PREPROCESSING stage from V1!
    assert "PREPROCESSING" not in tracker_v2.stages

    # 3. Mark success for sample V2
    tracker_v2.mark_success(
        stage_name="PREPROCESSING",
        elapsed_sec=2.0,
        artifacts={"settings": dummy_artifact},
        composite_hash="comp_hash_v2",
    )

    # 4. Re-load sample V1: its cache is still preserved!
    reloaded_v1 = StageTracker.load_or_create(cache_dir, sample_id=sample_id, sample_hash=hash_v1)
    assert reloaded_v1.is_stage_valid("PREPROCESSING", [dummy_artifact], "comp_hash_v1") is True

    # 5. Re-load sample V2: its cache is also preserved!
    reloaded_v2 = StageTracker.load_or_create(cache_dir, sample_id=sample_id, sample_hash=hash_v2)
    assert reloaded_v2.is_stage_valid("PREPROCESSING", [dummy_artifact], "comp_hash_v2") is True


def test_inspect_operational_cache_filters_by_sample_hash(tmp_path: Path):
    """
    Verifies that inspect_operational_cache skips cached stages if the target sample's
    current file hash does not match the cached hash.
    """
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir(parents=True)

    sample_id = "SAMPLE_A"
    hash_old = "aaaa1111aaaa1111aaaa1111aaaa1111aaaa1111aaaa1111aaaa1111aaaa1111"
    hash_new = "bbbb2222bbbb2222bbbb2222bbbb2222bbbb2222bbbb2222bbbb2222bbbb2222"

    tracker = StageTracker.load_or_create(cache_dir, sample_id=sample_id, sample_hash=hash_old)
    tracker.mark_success("PREPROCESSING", 1.0, {}, "comp_hash_old")

    # When inspecting without hash filter, cache is visible
    raw_cache = StageTracker.inspect_operational_cache(tmp_path, target_sample_ids=[sample_id])
    assert sample_id in raw_cache
    assert "PREPROCESSING" in raw_cache[sample_id]

    # When inspecting with current sample hash matching hash_old: cache is valid
    matching_cache = StageTracker.inspect_operational_cache(
        tmp_path,
        target_sample_ids=[sample_id],
        target_sample_hashes={sample_id: hash_old},
    )
    assert sample_id in matching_cache

    # When inspecting with new sample hash (different sample with same name): stale cache is ignored!
    mismatched_cache = StageTracker.inspect_operational_cache(
        tmp_path,
        target_sample_ids=[sample_id],
        target_sample_hashes={sample_id: hash_new},
    )
    assert sample_id not in mismatched_cache


def test_handle_cache_policy_options_and_banner(tmp_path: Path, monkeypatch, capsys):
    """Verifies that the prompt offers only [1] Reuse and [2] Overwrite, and Ctrl+C aborts."""
    op_dir = tmp_path / "operational_files"
    cache_dir = op_dir / "cache"
    cache_dir.mkdir(parents=True)

    tracker = StageTracker.load_or_create(cache_dir, sample_id="EGIZIANOANTICO1", sample_hash="1234567890abcdef")
    tracker.mark_success("ALIGNMENT", 1.0, {}, "hash1")

    runner = EndToEndPipelineRunner()
    runner.operational_files_dir = op_dir

    # Simulate user pressing Enter (default = 1: reuse)
    monkeypatch.setattr("sys.stdin.isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda prompt: "")
    decision = runner.handle_cache_policy(
        cache_mode="prompt",
        target_sample_ids=["EGIZIANOANTICO1"],
        target_sample_hashes={"EGIZIANOANTICO1": "1234567890abcdef"},
    )
    assert decision == "reuse"

    out = capsys.readouterr().out
    # Check that banner is compact and does NOT contain option 3
    assert "Cached Stages Detected" in out
    assert "[1] Reuse cache" in out
    assert "[2] Overwrite cache" in out
    assert "[3]" not in out

    # Simulate user choosing option 2 (overwrite)
    monkeypatch.setattr("builtins.input", lambda prompt: "2")
    decision2 = runner.handle_cache_policy(
        cache_mode="prompt",
        target_sample_ids=["EGIZIANOANTICO1"],
    )
    assert decision2 == "overwrite"
    assert runner.force_recompute is True
    # Crucial assertion: status file and cache artifacts MUST NOT be deleted upfront!
    assert tracker.status_file.is_file()

    # Simulate KeyboardInterrupt (Ctrl+C)
    def raise_keyboard_interrupt(prompt):
        raise KeyboardInterrupt()

    monkeypatch.setattr("builtins.input", raise_keyboard_interrupt)
    with pytest.raises(SystemExit) as exc_info:
        runner.handle_cache_policy(cache_mode="prompt")
    assert exc_info.value.code == 0


def test_overwrite_cache_policy_does_not_purge_files(tmp_path: Path):
    """
    Verifies that choosing overwrite policy sets force_recompute=True
    without purging or unlinking any cached deliverables on disk.
    """
    op_dir = tmp_path / "operational_files"
    cache_dir = op_dir / "cache"
    align_dir = op_dir / "alignment" / "SAMPLE_Z"
    align_dir.mkdir(parents=True)
    cache_dir.mkdir(parents=True)

    dummy_bam = align_dir / "SAMPLE_Z.filtered.dedup.bam"
    dummy_bam.write_bytes(b"BAM_DUMMY_CONTENT")

    tracker = StageTracker.load_or_create(cache_dir, sample_id="SAMPLE_Z")
    tracker.mark_success("ALIGN_AND_DEDUP", 1.0, {"bam": dummy_bam}, "composite_hash_z")

    assert tracker.status_file.is_file()
    assert dummy_bam.is_file()

    runner = EndToEndPipelineRunner()
    runner.operational_files_dir = op_dir

    policy = runner.handle_cache_policy(
        cache_mode="overwrite",
        target_sample_ids=["SAMPLE_Z"],
    )
    assert policy == "overwrite"
    assert runner.force_recompute is True

    # Files MUST still exist on disk so failed operations do not lose prior valid artifacts!
    assert tracker.status_file.is_file()
    assert dummy_bam.is_file()

    # StageTracker with force_recompute=True should now bypass cache
    reloaded = StageTracker.load_or_create(
        cache_dir,
        sample_id="SAMPLE_Z",
        force_recompute=runner.force_recompute,
    )
    assert reloaded.is_stage_valid("ALIGN_AND_DEDUP", [dummy_bam], "composite_hash_z") is False


def test_cli_cache_flags(monkeypatch):
    """Verifies -u / --use-cache and -n / --no-cache CLI options."""
    from src.pipeline.cli import parse_arguments

    # 1. Test -u / --use-cache
    monkeypatch.setattr(sys, "argv", ["pipeline", "-u"])
    args_u = parse_arguments()
    assert args_u.use_cache is True
    assert args_u.no_cache is False
    assert args_u.cache_mode == "reuse"

    monkeypatch.setattr(sys, "argv", ["pipeline", "--use-cache"])
    args_use = parse_arguments()
    assert args_use.cache_mode == "reuse"

    # 2. Test -n / --no-cache
    monkeypatch.setattr(sys, "argv", ["pipeline", "-n"])
    args_n = parse_arguments()
    assert args_n.no_cache is True
    assert args_n.use_cache is False
    assert args_n.cache_mode == "overwrite"

    monkeypatch.setattr(sys, "argv", ["pipeline", "--no-cache"])
    args_no = parse_arguments()
    assert args_no.cache_mode == "overwrite"

    # 3. Test mutual exclusivity (-u and -n together should fail)
    monkeypatch.setattr(sys, "argv", ["pipeline", "-u", "-n"])
    with pytest.raises(SystemExit):
        parse_arguments()

