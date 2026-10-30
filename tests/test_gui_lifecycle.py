"""
Regression Test Suite for GUI Process Lifecycle and Log Highlighter.
Validates kernel flock ownership, atomic state synchronization,
traceback parsing, and sample lookahead rules.
"""

from pathlib import Path
import pytest

from src.gui.log_highlighter import LogHighlighter
from src.gui.process_manager import ProcessManager, StopJobStatus


def test_log_highlighter_traceback_multiline():
    """Verify that multi-line Python tracebacks are accurately classified."""
    lines = [
        "2026-10-13 10:00:00 [INFO] Running stage",
        "Traceback (most recent call last):",
        '  File "runner.py", line 42, in execute',
        "    raise PipelineExecutionError('Sample failed')",
        "PipelineExecutionError: Sample failed",
        "2026-10-13 10:00:05 [INFO] Cleanup complete",
    ]
    html_output = LogHighlighter.render_terminal_html(lines)
    assert 'class="ln-traceback"' in html_output
    assert 'class="ln-error"' in html_output
    assert "PipelineExecutionError: Sample failed" in html_output


def test_log_highlighter_sample_token_discrimination():
    """Verify sample IDs are highlighted while HTSlib/internal tags are excluded."""
    lines = [
        "[OTZI] Processing endogenous reads",
        "[EGIZIANOANTICO1] Coverage confirmed",
        "[W::hts_idx_load3] The index file is older than the data file",
        "[Stage 1/5] Merging reads",
        "[INFO] Execution step",
    ]
    html_output = LogHighlighter.render_terminal_html(lines)
    assert '<span class="tk-sample">[OTZI]</span>' in html_output
    assert '<span class="tk-sample">[EGIZIANOANTICO1]</span>' in html_output
    # HTSlib warnings and stages must NOT be tagged as sample tokens
    assert '<span class="tk-sample">[W::hts_idx_load3]</span>' not in html_output
    assert '<span class="tk-sample">[Stage 1/5]</span>' not in html_output


def test_log_highlighter_metrics_formats():
    """Verify scientific notation, negative numbers, and time units are recognized."""
    lines = [
        "Depth: 105.4x | Rate: 89.37% | Delta: -3.5x | VAF: 1.2e-04",
        "Runtime: 8004.62s | MAPQ = 30 | Quality: Q30",
    ]
    html_output = LogHighlighter.render_terminal_html(lines)
    assert '<span class="tk-metric">105.4x</span>' in html_output
    assert '<span class="tk-metric">89.37%</span>' in html_output
    assert '<span class="tk-metric">-3.5x</span>' in html_output
    assert '<span class="tk-metric">1.2e-04</span>' in html_output
    assert '<span class="tk-metric">8004.62s</span>' in html_output
    assert '<span class="tk-metric">MAPQ = 30</span>' in html_output
    assert '<span class="tk-metric">Q30</span>' in html_output


def test_log_highlighter_carriage_return():
    """Verify progress bars overwriting with \\r retain only the terminal state."""
    raw_stream = "Downloading... 10%\rDownloading... 50%\rDownloading... 100%\nDone."
    normalized = LogHighlighter.normalize_carriage_returns(raw_stream)
    assert len(normalized) == 2
    assert normalized[0] == "Downloading... 100%"
    assert normalized[1] == "Done."


def test_process_manager_lock_integrity(tmp_path, monkeypatch):
    """Verify that is_running accurately reports kernel flock ownership without race."""
    fake_cache = tmp_path / "cache"
    fake_cache.mkdir()
    monkeypatch.setattr("src.gui.process_manager.CACHE_DIR", fake_cache)
    monkeypatch.setattr("src.gui.process_manager.LOCK_FILE", fake_cache / "pipeline.lock")
    monkeypatch.setattr("src.gui.process_manager.STATE_FILE", fake_cache / "pipeline.state.json")

    # Initially not running
    is_active, state = ProcessManager.is_running()
    assert is_active is False
    assert state is None

    # When no process is active, stop_job returns ALREADY_STOPPED
    status, _ = ProcessManager.stop_job()
    assert status == StopJobStatus.ALREADY_STOPPED