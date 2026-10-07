"""
Pipeline Stage Tracker & Caching Engine.
Centralizes disk-backed execution states in operational_files/cache/,
verifies artifact validity, and evaluates cryptographic composite hashes.
"""

from dataclasses import dataclass, field
from datetime import datetime, timezone
import json
import logging
from pathlib import Path
import shutil
from typing import Any, Dict, List, Optional

logger = logging.getLogger("mtDNA_Pipeline.StageTracker")


@dataclass
class StageTracker:
    """Tracks physical artifact validity and cross-stage dependency digests."""

    sample_id: str
    status_file: Path
    stages: Dict[str, str] = field(default_factory=dict)
    runtimes: Dict[str, float] = field(default_factory=dict)
    artifacts: Dict[str, str] = field(default_factory=dict)
    stage_hashes: Dict[str, str] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def load_or_create(
        cls,
        cache_dir: Path,
        sample_id: str,
        legacy_dir: Optional[Path] = None,
    ) -> "StageTracker":
        """
        Load an existing status tracker from the centralized cache directory,
        migrating legacy per-sample alignment state files if present.
        """
        cache_dir.mkdir(parents=True, exist_ok=True)
        status_path = cache_dir / f"{sample_id}.stage_status.json"

        # Migrate from legacy operational_files/alignment/<sample_id>/ if needed
        if not status_path.is_file() and legacy_dir:
            legacy_file = legacy_dir / f"{sample_id}.stage_status.json"
            if legacy_file.is_file():
                try:
                    shutil.copy2(legacy_file, status_path)
                    logger.info(f"[{sample_id}] Migrated status cache to centralized folder: {status_path}")
                except Exception as exc:
                    logger.warning(f"[{sample_id}] Failed migrating legacy cache file: {exc}")

        if status_path.is_file():
            try:
                with open(status_path, "r", encoding="utf-8") as f:
                    data = json.load(f)
                return cls(
                    sample_id=sample_id,
                    status_file=status_path,
                    stages=data.get("stages", {}),
                    runtimes=data.get("runtimes", {}),
                    artifacts=data.get("artifacts", {}),
                    stage_hashes=data.get("stage_hashes", {}),
                    metadata=data.get("metadata", {}),
                )
            except Exception as exc:
                logger.warning(f"[{sample_id}] Corrupted stage tracker on disk, re-initializing: {exc}")

        return cls(sample_id=sample_id, status_file=status_path)

    def is_stage_valid(
        self,
        stage_name: str,
        expected_artifacts: List[Path],
        current_composite_hash: str,
    ) -> bool:
        """Evaluate whether a pipeline stage can be bypassed safely via cached deliverables."""
        if self.stages.get(stage_name) != "SUCCESS":
            return False

        saved_hash = self.stage_hashes.get(stage_name)
        if saved_hash != current_composite_hash:
            logger.info(
                f"[{self.sample_id}] Stage {stage_name} cache invalidated: "
                f"upstream inputs or stage configuration have changed."
            )
            return False

        for art in expected_artifacts:
            if not art.is_file() or art.stat().st_size == 0:
                return False
            if art.suffix == ".bam":
                bai = Path(f"{art}.bai")
                alt_bai = art.with_suffix(".bai")
                if not (bai.is_file() or alt_bai.is_file()):
                    return False
        return True

    def mark_success(
        self,
        stage_name: str,
        elapsed_sec: float,
        artifacts: Dict[str, Path],
        composite_hash: str,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Record stage success, runtime, emitted artifacts, and composite hash to disk."""
        self.stages[stage_name] = "SUCCESS"
        self.runtimes[stage_name] = round(elapsed_sec, 2)
        self.stage_hashes[stage_name] = composite_hash
        for k, path in artifacts.items():
            self.artifacts[f"{stage_name}:{k}"] = str(path)
        if metadata:
            self.metadata[stage_name] = metadata

        try:
            with open(self.status_file, "w", encoding="utf-8") as f:
                json.dump(
                    {
                        "sample_id": self.sample_id,
                        "updated_at": datetime.now(timezone.utc).isoformat(),
                        "stages": self.stages,
                        "runtimes": self.runtimes,
                        "stage_hashes": self.stage_hashes,
                        "artifacts": self.artifacts,
                        "metadata": self.metadata,
                    },
                    f,
                    indent=2,
                )
        except Exception as exc:
            logger.warning(f"[{self.sample_id}] Failed saving status tracker: {exc}")
