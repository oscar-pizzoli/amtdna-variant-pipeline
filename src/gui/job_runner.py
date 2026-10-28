"""
Detached Job Runner Supervisor.
Holds continuous advisory flock, records start-time ticks for process identity verification,
registers run-scoped environment tags, and writes exit codes atomically via filesystem rename.
"""

import argparse
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CACHE_DIR = PROJECT_ROOT / "operational_files" / "cache"
LOGS_DIR = PROJECT_ROOT / "operational_files" / "logs"

LOCK_FILE = CACHE_DIR / "pipeline.lock"
STATE_FILE = CACHE_DIR / "pipeline.state.json"
EXIT_CODE_FILE = CACHE_DIR / "pipeline.exit_code"
RUN_LOG_FILE = LOGS_DIR / "pipeline_run.log"


def get_self_starttime_ticks() -> int | None:
    """Reads self starttime jiffies from /proc/self/stat to prevent PID recycling."""
    stat_path = Path("/proc/self/stat")
    if not stat_path.is_file():
        return None
    try:
        raw = stat_path.read_text(encoding="utf-8", errors="replace").strip()
        r_paren = raw.rfind(")")
        if r_paren != -1:
            remainder = raw[r_paren + 1:].split()
            if len(remainder) > 19:
                return int(remainder[19])
    except Exception:
        pass
    return None


def atomic_write_text(target_path: Path, content: str) -> None:
    target_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = target_path.with_name(f".tmp_{target_path.name}_{os.getpid()}")
    with open(tmp_path, "w", encoding="utf-8") as fp:
        fp.write(content)
        fp.flush()
        os.fsync(fp.fileno())
    tmp_path.replace(target_path)


def atomic_write_json(target_path: Path, data: dict) -> None:
    target_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = target_path.with_name(f".tmp_{target_path.name}_{os.getpid()}")
    with open(tmp_path, "w", encoding="utf-8") as fp:
        json.dump(data, fp, indent=2)
        fp.flush()
        os.fsync(fp.fileno())
    tmp_path.replace(target_path)


def main() -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    parser = argparse.ArgumentParser(description="Pipeline Job Supervisor Daemon.")
    parser.add_argument("--run-id", required=True, help="Unique run identifier.")
    parser.add_argument("--data-type", default="BAM")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--cpu-budget-fraction", default="0.5")
    parser.add_argument("--cache-mode", default="reuse", choices=["prompt", "reuse", "overwrite"])
    parser.add_argument("--log-level", default="INFO")
    args, unknown = parser.parse_known_args()

    # 1. Acquire advisory lock on LOCK_FILE
    lock_fd = os.open(LOCK_FILE, os.O_CREAT | os.O_RDWR)
    try:
        fcntl.flock(lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except (BlockingIOError, OSError) as err:
        sys.stderr.write(f"Supervisor could not acquire flock: {err}\n")
        os.close(lock_fd)
        sys.exit(2)

    my_pid = os.getpid()
    my_pgid = os.getpgid(0)
    my_ticks = get_self_starttime_ticks()
    start_utc = datetime.now(timezone.utc).isoformat()

    # 2. Write state descriptor including boot start-time ticks and cache mode
    state_payload = {
        "run_id": args.run_id,
        "pid": my_pid,
        "pgid": my_pgid,
        "start_time_utc": start_utc,
        "start_time_ticks": my_ticks,
        "data_type": args.data_type,
        "input_dir": args.input_dir,
        "cache_mode": args.cache_mode,
    }
    atomic_write_json(STATE_FILE, state_payload)

    # 3. Setup output log stream
    archive_log = LOGS_DIR / f"{args.run_id}.log"
    cmd = [
        sys.executable,
        "-u",
        str(PROJECT_ROOT / "scripts" / "run_complete_pipeline.py"),
        "--data-type", args.data_type,
        "--input-dir", args.input_dir,
        "--cpu-budget-fraction", args.cpu_budget_fraction,
        "--cache-mode", args.cache_mode,
        "--log-level", args.log_level,
    ] + unknown

    env = os.environ.copy()
    env["PIPELINE_RUN_ID"] = args.run_id

    ret_code = 1
    with open(RUN_LOG_FILE, "w", encoding="utf-8") as log_fp:
        log_fp.write(f"=== JOB LAUNCHED: {args.run_id} ({start_utc}) ===\n")
        log_fp.write(f"Command: {' '.join(cmd)}\n\n")
        log_fp.flush()

        proc = subprocess.Popen(
            cmd,
            stdout=log_fp,
            stderr=subprocess.STDOUT,
            cwd=str(PROJECT_ROOT),
            env=env,
        )

        ret_code = proc.wait()

        finish_utc = datetime.now(timezone.utc).isoformat()
        log_fp.write(f"\n=== JOB COMPLETED AT {finish_utc} (EXIT CODE: {ret_code}) ===\n")
        log_fp.flush()

    try:
        shutil.copy2(RUN_LOG_FILE, archive_log)
    except OSError:
        pass

    # 4. Atomic terminal status writing
    atomic_write_text(EXIT_CODE_FILE, str(ret_code))
    if STATE_FILE.is_file():
        STATE_FILE.unlink(missing_ok=True)

    fcntl.flock(lock_fd, fcntl.LOCK_UN)
    os.close(lock_fd)

    sys.exit(ret_code)


if __name__ == "__main__":
    main()