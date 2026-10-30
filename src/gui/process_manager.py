"""
Pipeline Background Process Lifecycle Engine.
Guarantees verified PID identity, boot-time tick validation against PID recycling,
cross-platform path translation, and deterministic process group termination.
"""

from dataclasses import dataclass
from datetime import datetime, timezone
from enum import Enum
import fcntl
import json
import logging
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time
from typing import Any, Dict, Optional, Tuple

logger = logging.getLogger("mtDNA_Pipeline.ProcessManager")

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CACHE_DIR = PROJECT_ROOT / "operational_files" / "cache"
LOGS_DIR = PROJECT_ROOT / "operational_files" / "logs"

LOCK_FILE = CACHE_DIR / "pipeline.lock"
STATE_FILE = CACHE_DIR / "pipeline.state.json"
EXIT_CODE_FILE = CACHE_DIR / "pipeline.exit_code"
RUN_LOG_FILE = LOGS_DIR / "pipeline_run.log"


class StopJobStatus(str, Enum):
    """Formal termination verdicts returned by stop_job."""
    STOPPED_CONFIRMED = "STOPPED_CONFIRMED"
    ALREADY_STOPPED = "ALREADY_STOPPED"
    SECURITY_REJECTED = "SECURITY_REJECTED"
    STOP_FAILED = "STOP_FAILED"


@dataclass
class JobState:
    """Active background job metadata descriptor."""
    run_id: str
    pid: int
    pgid: int
    start_time_utc: str
    start_time_ticks: Optional[int]
    data_type: str
    input_dir: str
    cache_mode: str = "reuse"


class PathValidator:
    """Handles bidirectional Windows/WSL path translation and runtime existence checks."""

    @staticmethod
    def is_windows_host() -> bool:
        return sys.platform == "win32"

    @classmethod
    def resolve_paths(cls, raw_path_str: str) -> Tuple[bool, Path, str]:
        """
        Resolves input path returning:
        - exists: True if path physically exists on the current host
        - runtime_path: Path instance valid for the current Python interpreter
        - pipeline_path_arg: POSIX path string formatted for the WSL/Linux pipeline
        """
        cleaned = raw_path_str.strip().strip("'\"")
        if not cleaned:
            return False, PROJECT_ROOT, str(PROJECT_ROOT)

        m_drive = re.match(r"^([a-zA-Z]):[\\/](.*)", cleaned)
        m_unc = cleaned.startswith(r"\\wsl") or cleaned.startswith(r"//wsl")

        if cls.is_windows_host():
            if m_drive:
                win_path = Path(cleaned)
                drive_letter = m_drive.group(1).lower()
                subpath = m_drive.group(2).replace("\\", "/")
                return win_path.is_dir(), win_path, f"/mnt/{drive_letter}/{subpath}"
            elif m_unc:
                win_path = Path(cleaned)
                parts = [p for p in cleaned.replace("\\", "/").split("/") if p]
                wsl_posix = "/" + "/".join(parts[2:]) if len(parts) >= 3 else cleaned
                return win_path.is_dir(), win_path, wsl_posix
            else:
                p = Path(cleaned)
                return p.is_dir(), p, cleaned.replace("\\", "/")
        else:
            if m_drive:
                drive_letter = m_drive.group(1).lower()
                subpath = m_drive.group(2).replace("\\", "/")
                wsl_path = Path(f"/mnt/{drive_letter}/{subpath}")
                return wsl_path.is_dir(), wsl_path, str(wsl_path)
            elif m_unc:
                parts = [p for p in cleaned.replace("\\", "/").split("/") if p]
                wsl_path = Path("/" + "/".join(parts[2:])) if len(parts) >= 3 else Path(cleaned)
                return wsl_path.is_dir(), wsl_path, str(wsl_path)
            else:
                p = Path(cleaned)
                return p.is_dir(), p, str(p)


class ProcessManager:
    """Coordinates asynchronous jobs with kernel locking and PID verification."""

    @staticmethod
    def get_process_starttime_ticks(pid: int) -> Optional[int]:
        """
        Extracts field 22 (starttime jiffies post-boot) from /proc/<pid>/stat.
        Safely isolates the command name to handle parentheses in executable names.
        """
        stat_path = Path(f"/proc/{pid}/stat")
        if not stat_path.is_file():
            return None
        try:
            raw = stat_path.read_text(encoding="utf-8", errors="replace").strip()
            r_paren = raw.rfind(")")
            if r_paren == -1:
                return None
            remainder = raw[r_paren + 1:].split()
            # In /proc/<pid>/stat: field 1 is pid, field 2 is comm (inside parentheses).
            # Field 3 (state) is index 0 of remainder. Field 22 (starttime) is index 19.
            if len(remainder) > 19:
                return int(remainder[19])
            return None
        except (OSError, ValueError):
            return None

    @classmethod
    def _verify_process_identity(cls, pid: int, state: Optional[JobState]) -> bool:
        """
        Airtight protection against PID recycling:
        1. Checks command line tokens.
        2. Checks process creation ticks against state descriptor.
        """
        cmdline_path = Path(f"/proc/{pid}/cmdline")
        if not cmdline_path.is_file():
            return False

        try:
            raw = cmdline_path.read_bytes().decode("utf-8", errors="replace")
            is_pipeline_cmd = "job_runner.py" in raw or "run_complete_pipeline.py" in raw
            if not is_pipeline_cmd:
                return False

            if state and state.start_time_ticks is not None:
                current_ticks = cls.get_process_starttime_ticks(pid)
                if current_ticks is None or current_ticks != state.start_time_ticks:
                    logger.error(
                        f"PID {pid} recycled: recorded ticks={state.start_time_ticks}, "
                        f"active ticks={current_ticks}."
                    )
                    return False

            return True
        except (OSError, PermissionError):
            return False

    @classmethod
    def get_active_state(cls) -> Optional[JobState]:
        """Loads certified JobState if present on disk and valid."""
        if not STATE_FILE.is_file():
            return None
        try:
            with open(STATE_FILE, "r", encoding="utf-8") as fp:
                d = json.load(fp)
            return JobState(
                run_id=d["run_id"],
                pid=int(d["pid"]),
                pgid=int(d["pgid"]),
                start_time_utc=d["start_time_utc"],
                start_time_ticks=d.get("start_time_ticks"),
                data_type=d["data_type"],
                input_dir=d["input_dir"],
                cache_mode=d.get("cache_mode", "reuse"),
            )
        except Exception:
            return None

    @classmethod
    def is_running(cls) -> Tuple[bool, Optional[JobState]]:
        """Determines execution state via advisory flock and PID identity audit."""
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        if not LOCK_FILE.is_file():
            return False, None

        try:
            fd = os.open(LOCK_FILE, os.O_RDWR)
        except OSError:
            return False, None

        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            fcntl.flock(fd, fcntl.LOCK_UN)
            os.close(fd)
            return False, None
        except (BlockingIOError, OSError):
            os.close(fd)
            state = cls.get_active_state()
            if state and cls._verify_process_identity(state.pid, state):
                return True, state
            elif state and not cls._verify_process_identity(state.pid, state):
                logger.warning(f"Lock held but PID {state.pid} failed identity check.")
                return False, None
            return True, state

    @classmethod
    def start_job(
        cls,
        data_type: str,
        input_dir_arg: str,
        cpu_fraction: float,
        log_level: str,
        cache_mode: str = "reuse",
        handshake_timeout: float = 3.5,
    ) -> Tuple[bool, str]:
        """Launches detached supervisor with verified kernel handshake."""
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        LOGS_DIR.mkdir(parents=True, exist_ok=True)

        running, active_state = cls.is_running()
        if running:
            return False, f"Job already active (PID: {active_state.pid if active_state else 'unknown'})."

        if EXIT_CODE_FILE.is_file():
            EXIT_CODE_FILE.unlink(missing_ok=True)
        if STATE_FILE.is_file():
            STATE_FILE.unlink(missing_ok=True)

        run_id = f"run_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}"

        cmd = [
            sys.executable,
            "-u",
            str(PROJECT_ROOT / "src" / "gui" / "job_runner.py"),
            "--run-id", run_id,
            "--data-type", data_type,
            "--input-dir", input_dir_arg,
            "--cpu-budget-fraction", str(cpu_fraction),
            "--cache-mode", cache_mode,
            "--log-level", log_level,
        ]

        env = os.environ.copy()
        env["PIPELINE_RUN_ID"] = run_id

        try:
            subprocess.Popen(
                cmd,
                cwd=str(PROJECT_ROOT),
                start_new_session=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                stdin=subprocess.DEVNULL,
                env=env,
            )
        except Exception as err:
            return False, f"Failed spawning supervisor process: {err}"

        t_deadline = time.time() + handshake_timeout
        while time.time() < t_deadline:
            is_active, state = cls.is_running()
            if is_active and state is not None and state.run_id == run_id:
                return True, f"Job successfully initialized (Run ID: {run_id}, PID: {state.pid}, Cache: {cache_mode.upper()})."
            time.sleep(0.05)

        cls.cleanup_docker_containers(run_id)
        return False, "Supervisor failed to acquire kernel lock within handshake deadline."

    @classmethod
    def cleanup_docker_containers(cls, run_id: Optional[str] = None) -> int:
        """Terminates ONLY docker containers explicitly tagged with pipeline_run_id."""
        if not run_id or not shutil.which("docker"):
            return 0

        terminated_count = 0
        try:
            res = subprocess.run(
                ["docker", "ps", "-q", "--filter", f"label=pipeline_run_id={run_id}"],
                capture_output=True,
                text=True,
                check=False,
            )
            container_ids = [cid.strip() for cid in res.stdout.splitlines() if cid.strip()]
            if container_ids:
                logger.info(f"Stopping {len(container_ids)} Docker containers scoped to run_id={run_id}...")
                subprocess.run(["docker", "stop", "-t", "2"] + container_ids, capture_output=True, check=False)
                subprocess.run(["docker", "rm", "-f"] + container_ids, capture_output=True, check=False)
                terminated_count = len(container_ids)
        except Exception as exc:
            logger.warning(f"Error during run-scoped Docker cleanup: {exc}")

        return terminated_count

    @classmethod
    def stop_job(cls, sigterm_timeout: float = 6.0, sigkill_timeout: float = 2.0) -> Tuple[StopJobStatus, str]:
        """Safely stops active job with dual-layer identity verification before signaling."""
        running, state = cls.is_running()
        if not running or state is None:
            return StopJobStatus.ALREADY_STOPPED, "No active job is currently executing."

        pid = state.pid
        pgid = state.pgid
        run_id = state.run_id

        # Security validation against PID recycling
        if not cls._verify_process_identity(pid, state):
            logger.critical(f"ABORT KILL: PID {pid} failed cmdline/tick validation.")
            return StopJobStatus.SECURITY_REJECTED, f"Process {pid} failed security check. Signal blocked."

        try:
            actual_pgid = os.getpgid(pid)
            if actual_pgid != pgid:
                logger.critical(f"PGID mismatch: expected {pgid}, found {actual_pgid}. Aborting.")
                return StopJobStatus.SECURITY_REJECTED, "Process group mismatch. Termination aborted."
        except OSError as err:
            return StopJobStatus.ALREADY_STOPPED, f"Process {pid} already terminated: {err}"

        # 1. Send SIGTERM to process group
        try:
            os.killpg(pgid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        except OSError as err:
            logger.warning(f"SIGTERM error: {err}")

        # 2. Poll for flock release
        t_deadline = time.time() + sigterm_timeout
        while time.time() < t_deadline:
            is_active, _ = cls.is_running()
            if not is_active:
                cls.cleanup_docker_containers(run_id)
                return StopJobStatus.STOPPED_CONFIRMED, f"Process group {pgid} terminated cleanly via SIGTERM."
            time.sleep(0.15)

        # 3. Escalate to SIGKILL if processes remain active
        logger.warning(f"Process group {pgid} did not yield to SIGTERM. Sending SIGKILL.")
        try:
            os.killpg(pgid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        except OSError as err:
            logger.warning(f"SIGKILL error: {err}")

        # 4. Final verification window
        t_kill_deadline = time.time() + sigkill_timeout
        while time.time() < t_kill_deadline:
            is_active, _ = cls.is_running()
            if not is_active:
                cls.cleanup_docker_containers(run_id)
                return StopJobStatus.STOPPED_CONFIRMED, f"Process group {pgid} terminated forcefully via SIGKILL."
            time.sleep(0.1)

        return StopJobStatus.STOP_FAILED, f"Processes in group {pgid} could not be terminated."

    @staticmethod
    def get_last_exit_code() -> Optional[int]:
        if EXIT_CODE_FILE.is_file():
            try:
                return int(EXIT_CODE_FILE.read_text(encoding="utf-8").strip())
            except (ValueError, OSError):
                return None
        return None