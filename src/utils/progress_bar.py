"""
Gradient Progress Bar and Subprocess Execution Tracker for mtDNA Pipeline.

Provides ANSI TrueColor (Red -> Orange -> Yellow -> Lime -> Green) progress bars
for long-running analytical stages (Adapter trimming, BWA alignment, Schmutzi MCMC).
Designed as an independent, modular utility to avoid cluttering core pipeline runners.
"""

import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from typing import BinaryIO, List, Optional, TextIO, Tuple, Union


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
        unit: str = "%",
    ):
        self.label = desc if (desc is not None and not label) else label
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
    Tracks multi-step biological pipelines (e.g. Alignment branches, collation, deduplication, filtering)
    distributing progress monotonically along a shared 0% -> 100% Red-Green gradient bar.
    """

    def __init__(
        self,
        label: Union[str, List[Tuple[str, float]]] = "",
        stages: Optional[List[Tuple[str, float]]] = None,
        width: int = 25,
        stream: Optional[TextIO] = None,
        overall_desc: Optional[str] = None,
        desc: Optional[str] = None,
    ):
        """
        Args:
            label: Descriptive header (e.g. '[sample_id] [Alignment & Dedup]').
            stages: List of (stage_name, relative_weight). Weights are normalized to 1.0.
            width: Width of the progress bar in characters.
            stream: Target output text stream (default: sys.stdout).
        """
        if isinstance(label, list):
            actual_stages = label
            actual_label = overall_desc or desc or ""
        else:
            actual_stages = stages or []
            actual_label = overall_desc or desc or label

        self.bar = GradientProgressBar(label=actual_label, width=width, stream=stream)
        total_weight = sum(weight for _, weight in actual_stages) or 1.0
        self.stages = [(name, weight / total_weight) for name, weight in actual_stages]
        self.stage_names = [name for name, _ in actual_stages]

        # Precompute cumulative offsets for each stage
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
    stdout_file: Optional[Union[str, Path]] = None,
    stderr_file: Optional[Union[str, Path]] = None,
    timeout: Optional[float] = None,
    expected_duration: float = 30.0,
    expected_duration_sec: Optional[float] = None,
    check_output_files: Optional[List[Optional[Path]]] = None,
    target_output_file: Optional[Path] = None,
    target_bytes: Optional[int] = None,
    expected_output_bytes: Optional[int] = None,
    log_file_to_tail: Optional[Path] = None,
    poll_interval: float = 0.5,
    width: int = 18,
    append_logs: bool = False,
    check: bool = False,
) -> subprocess.CompletedProcess:
    """
    Executes a subprocess while displaying a live, updating Red-to-Green gradient progress bar.
    Polls subprocess status, elapsed time, and output growth without blocking blindly.
    """
    actual_label = desc or label or "Subprocess"
    exp_dur = expected_duration_sec if expected_duration_sec is not None else expected_duration
    exp_bytes = expected_output_bytes if expected_output_bytes is not None else target_bytes

    monitored_files = list(check_output_files or [])
    if target_output_file:
        monitored_files.append(target_output_file)

    bar = GradientProgressBar(label=actual_label, width=width)
    t0 = time.perf_counter()

    # Determine files to manage if stdout_f / stderr_f are not already open handles
    out_handle = stdout_f
    err_handle = stderr_f
    opened_files = []

    try:
        if out_handle is None and stdout_file is not None:
            mode = "a" if append_logs else "w"
            out_handle = open(stdout_file, mode, encoding="utf-8")
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

                # 1. Progress Estimation
                if exp_bytes and exp_bytes > 0 and active_files:
                    current_bytes = sum(f.stat().st_size for f in active_files if f.is_file())
                    byte_pct = min(0.95, current_bytes / max(1, exp_bytes))
                    status_detail = f"{format_bytes(current_bytes)} written | {format_duration(elapsed)}"
                    bar.update(byte_pct, status_detail)
                else:
                    # Smooth asymptotic curve: ~50% at 0.5x, ~75% at 1.0x, asymptotic to 95%
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
