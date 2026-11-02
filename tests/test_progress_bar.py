"""
Unit tests for the Red-to-Green Gradient Progress Bar and SubstageTracker.
"""

import io
from pathlib import Path
import subprocess
import sys
import tempfile
import pytest

from src.utils.progress_bar import (
    get_gradient_color,
    supports_color,
    format_duration,
    format_bytes,
    GradientProgressBar,
    SubstageTracker,
    run_tracked_subprocess,
    RESET,
    BOLD_RED,
    BOLD_GREEN,
)


def test_gradient_color_stops():
    """Verifies that RGB colors interpolate smoothly from Red to Green."""
    c_start = get_gradient_color(0.0)
    assert "\x1b[38;2;239;68;68m" in c_start  # Red stop

    c_mid = get_gradient_color(0.5)
    assert "\x1b[38;2;234;179;8m" in c_mid   # Amber/Yellow stop

    c_end = get_gradient_color(1.0)
    assert "\x1b[38;2;34;197;94m" in c_end   # Green stop


def test_format_helpers():
    """Verifies duration and byte size formatting."""
    assert format_duration(5.2) == "5.2s"
    assert format_duration(125.0) == "2m 05s"

    assert format_bytes(500) == "500 B"
    assert format_bytes(2048) == "2.0 KB"
    assert format_bytes(5 * 1024 * 1024) == "5.0 MB"


def test_progress_bar_render_and_update():
    """Verifies in-place line rendering and completion."""
    buffer = io.StringIO()
    bar = GradientProgressBar(label="TestStage", width=20, stream=buffer)

    bar.update(0.2, status="initializing")
    bar.update(0.8, status="almost done")
    bar.finish(status="Complete ✓")

    output = buffer.getvalue()
    assert "TestStage" in output
    assert "20%" in output
    assert "80%" in output
    assert "100%" in output
    assert "Complete ✓" in output


def test_substage_tracker_weights_and_advancement():
    """Verifies that SubstageTracker distributes weighted percentages from 0% to 100%."""
    buffer = io.StringIO()
    stages = [
        ("Alignment", 40.0),
        ("Sorting", 30.0),
        ("Deduplication", 30.0),
    ]
    tracker = SubstageTracker(label="PipelineTest", stages=stages, stream=buffer)

    tracker.start_stage("Alignment", detail="running")
    tracker.update_substage(0.5, status="50% reads")
    # Half of 40% is 20%
    line_mid = tracker.bar.render(0.20)
    assert " 20%" in line_mid

    tracker.start_stage("Sorting")
    tracker.start_stage("Deduplication")
    tracker.finish("All Done")

    output = buffer.getvalue()
    assert "PipelineTest" in output
    assert "All Done" in output


def test_run_tracked_subprocess_lifecycle():
    """Verifies that run_tracked_subprocess executes a child process and logs progress."""
    with tempfile.TemporaryDirectory() as tmpdir:
        out_f = Path(tmpdir) / "test.out"
        err_f = Path(tmpdir) / "test.err"

        cmd = [
            sys.executable,
            "-c",
            "import sys; sys.stdout.write('hello stdout\\n'); sys.stderr.write('hello stderr\\n')",
        ]

        res = run_tracked_subprocess(
            cmd=cmd,
            desc="Lifecycle Test",
            stdout_file=out_f,
            stderr_file=err_f,
            expected_duration_sec=1.0,
            check=True,
        )

        assert res.returncode == 0
        assert out_f.read_text(encoding="utf-8").strip() == "hello stdout"
        assert err_f.read_text(encoding="utf-8").strip() == "hello stderr"


def test_log_highlighter_extract_progress():
    """Verifies parsing of sub-stage and block progress for Streamlit HUD."""
    from src.gui.log_highlighter import LogHighlighter

    # Test empty lines
    val, text = LogHighlighter.extract_current_progress([])
    assert val == 0.0

    # Test active progress bar
    lines_substage = [
        "2026-10-31 16:39:38 [INFO] ▶ BLOCK 2: DISCOVERY, INTEGRITY & PREPROCESSING (SampleDiscoveryEngine)",
        "[SAMPLE1] [Trimming Ends] [━━━━━━━━━━━━━            ]  45% | 14.2s",
    ]
    val, text = LogHighlighter.extract_current_progress(lines_substage)
    assert val == 0.45
    assert "45%" in text
    assert "[SAMPLE1] [Trimming Ends]" in text

    # Test block header match
    lines_block = [
        "2026-10-31 16:44:19 [INFO] ▶ BLOCK 5: TWO-TIER CONTAMINATION ESTIMATION & GATING (ContaminationBatchRunner)",
        "2026-10-31 16:44:20 [INFO] Starting contamination batch analysis...",
    ]
    val, text = LogHighlighter.extract_current_progress(lines_block)
    assert val == 0.75
    assert "Block 5/8" in text

    # Test pipeline completion match
    lines_complete = [
        "🎉 COMPLETE PIPELINE (BLOCKS 1 -> 8) FINISHED IN 3839.59s",
        "=== JOB COMPLETED AT 2026-10-31T17:41:24.268555+00:00 (EXIT CODE: 0) ===",
    ]
    val, text = LogHighlighter.extract_current_progress(lines_complete)
    assert val == 1.0
    assert "Complete Pipeline Finished" in text


def test_progress_bar_deduplication_and_bounding():
    """Verifies that multiple consecutive progress updates collapse into a single line."""
    from src.gui.log_highlighter import LogHighlighter
    from src.utils.progress_bar import strip_ansi

    # 1. Test log deduplication: 3 consecutive updates collapse to 1
    log_text = (
        "2026-10-31 16:39:38 [INFO] ▶ BLOCK 5: CONTAMINATION\n"
        "Schmutzi: MCMC [EgizianoAntico2] [━━━━━━━━━]  90% | 3m 40s\n"
        "Schmutzi: MCMC [EgizianoAntico2] [━━━━━━━━━]  92% | 3m 45s\n"
        "Schmutzi: MCMC [EgizianoAntico2] [━━━━━━━━━]  95% | 3m 54s\n"
        "2026-10-31 16:44:19 [INFO] ▶ BLOCK 6: CONSENSUS"
    )
    resolved = LogHighlighter.normalize_carriage_returns(log_text)
    assert len(resolved) == 3
    assert "95%" in resolved[1]

    # 2. Test terminal line length bounding
    bar = GradientProgressBar(label="Schmutzi: MCMC [EgizianoAntico2]", width=20)
    line = bar.render(0.95, status="3m 54s (iter 2) - running extra long status string that would wrap")
    vis_len = len(strip_ansi(line))
    # Line must fit within 80 columns without wrapping
    assert vis_len <= 79
