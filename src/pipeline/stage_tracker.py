"""
Pipeline Stage Tracker & Caching Engine.
Centralizes disk-backed execution states in operational_files/cache/,
verifies artifact validity, evaluates cryptographic composite hashes,
and provides cache inspection and selective purge utilities.
"""

from dataclasses import dataclass, field
from datetime import datetime, timezone
import json
import logging
from pathlib import Path
import shutil
from typing import Any, Dict, List, Optional, Set

logger = logging.getLogger("mtDNA_Pipeline.StageTracker")

STAGE_SUBDIRS: List[str] = [
    "preprocessing",
    "alignment",
    "damage",
    "contamination",
    "consensus",
    "haplogroup",
]


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
    sample_hash: Optional[str] = None
    force_recompute: bool = False

    @classmethod
    def load_or_create(
        cls,
        cache_dir: Path,
        sample_id: str,
        legacy_dir: Optional[Path] = None,
        sample_hash: Optional[str] = None,
        force_recompute: bool = False,
    ) -> "StageTracker":
        """
        Load an existing status tracker from the centralized cache directory,
        matching sample_id and sample_hash if provided. This ensures that different
        samples sharing the same sample_id never incorrectly reuse stale cache.
        """
        cache_dir.mkdir(parents=True, exist_ok=True)
        target_file: Optional[Path] = None
        data: Dict[str, Any] = {}

        # 1. If sample_hash is provided, check hash-specific status file first
        if sample_hash:
            short_h = sample_hash[:12]
            hash_file = cache_dir / f"{sample_id}_{short_h}.stage_status.json"
            if hash_file.is_file():
                target_file = hash_file

        # 2. Check canonical status file f"{sample_id}.stage_status.json"
        canonical_file = cache_dir / f"{sample_id}.stage_status.json"
        if target_file is None and canonical_file.is_file():
            try:
                with open(canonical_file, "r", encoding="utf-8") as f:
                    cand_data = json.load(f)
                cand_hash = cand_data.get("sample_hash") or cand_data.get("metadata", {}).get("sample_hash")
                if sample_hash and cand_hash and cand_hash != sample_hash:
                    logger.info(
                        f"[{sample_id}] Cached status {canonical_file.name} hash ({cand_hash[:8]}...) "
                        f"does not match current sample hash ({sample_hash[:8]}...). Starting fresh tracker."
                    )
                else:
                    target_file = canonical_file
                    data = cand_data
            except Exception as exc:
                logger.warning(f"[{sample_id}] Failed reading {canonical_file}: {exc}")

        # 3. Check legacy dir migration if needed
        if target_file is None and legacy_dir:
            legacy_file = legacy_dir / f"{sample_id}.stage_status.json"
            if legacy_file.is_file():
                try:
                    with open(legacy_file, "r", encoding="utf-8") as f:
                        cand_data = json.load(f)
                    cand_hash = cand_data.get("sample_hash")
                    if not (sample_hash and cand_hash and cand_hash != sample_hash):
                        dest_file = cache_dir / (
                            f"{sample_id}_{sample_hash[:12]}.stage_status.json"
                            if sample_hash
                            else canonical_file
                        )
                        shutil.copy2(legacy_file, dest_file)
                        target_file = dest_file
                        data = cand_data
                except Exception as exc:
                    logger.warning(f"[{sample_id}] Failed migrating legacy cache: {exc}")

        # If target_file was selected but data not loaded yet
        if target_file and not data:
            try:
                with open(target_file, "r", encoding="utf-8") as f:
                    data = json.load(f)
            except Exception as exc:
                logger.warning(f"[{sample_id}] Corrupted stage tracker {target_file}: {exc}")
                target_file = None

        if target_file and data:
            return cls(
                sample_id=sample_id,
                status_file=target_file,
                stages=data.get("stages", {}),
                runtimes=data.get("runtimes", {}),
                artifacts=data.get("artifacts", {}),
                stage_hashes=data.get("stage_hashes", {}),
                metadata=data.get("metadata", {}),
                sample_hash=data.get("sample_hash") or sample_hash,
                force_recompute=force_recompute,
            )

        # Fresh tracker initialization
        status_path = (
            cache_dir / f"{sample_id}_{sample_hash[:12]}.stage_status.json"
            if sample_hash
            else canonical_file
        )
        return cls(
            sample_id=sample_id,
            status_file=status_path,
            sample_hash=sample_hash,
            force_recompute=force_recompute,
        )

    @classmethod
    def inspect_operational_cache(
        cls,
        operational_files_dir: Path,
        target_sample_ids: Optional[List[str]] = None,
        target_sample_hashes: Optional[Dict[str, str]] = None,
    ) -> Dict[str, List[str]]:
        """
        Scans operational_files/ to detect existing cached stages and sample artifacts.
        If target_sample_hashes is provided, filters or validates that cached stages
        actually correspond to the current sample's content hash.
        Returns a dictionary mapping sample_id -> list of cached stages.
        """
        detected: Dict[str, Set[str]] = {}
        target_set = set(target_sample_ids) if target_sample_ids else None

        # 1. Check centralized stage status JSON files in operational_files/cache/
        cache_dir = operational_files_dir / "cache"
        if cache_dir.is_dir():
            for status_file in sorted(cache_dir.glob("*.stage_status.json")):
                stem = status_file.name.removesuffix(".stage_status.json")
                raw_s_id = stem.split("_")[0] if "_" in stem else stem
                if target_set is not None and raw_s_id not in target_set and stem not in target_set:
                    continue
                try:
                    with open(status_file, "r", encoding="utf-8") as f:
                        data = json.load(f)
                    s_id = data.get("sample_id", raw_s_id)
                    s_hash = data.get("sample_hash") or data.get("metadata", {}).get("sample_hash")

                    # If target_sample_hashes provided, verify hash match
                    if target_sample_hashes and s_id in target_sample_hashes:
                        expected_hash = target_sample_hashes[s_id]
                        if s_hash and expected_hash and s_hash != expected_hash:
                            # Stale cache from a different sample with the same name - skip!
                            continue

                    completed_stages = [
                        st_name for st_name, st_val in data.get("stages", {}).items() if st_val == "SUCCESS"
                    ]
                    if completed_stages:
                        detected.setdefault(s_id, set()).update(completed_stages)
                    else:
                        detected.setdefault(s_id, set()).add("CACHE_TRACKER")
                except Exception:
                    detected.setdefault(raw_s_id, set()).add("CACHE_TRACKER")

        # 2. Check physical module subdirectories in operational_files/
        for sub in STAGE_SUBDIRS:
            mod_dir = operational_files_dir / sub
            if not mod_dir.is_dir():
                continue
            for item in sorted(mod_dir.iterdir()):
                if item.is_dir():
                    s_id = item.name
                    if target_set is not None and s_id not in target_set:
                        continue
                    if any(item.iterdir()):
                        detected.setdefault(s_id, set()).add(sub.upper())

        return {s_id: sorted(list(stages)) for s_id, stages in sorted(detected.items())}

    @classmethod
    def purge_operational_cache(
        cls,
        operational_files_dir: Path,
        results_dir: Optional[Path] = None,
        target_sample_ids: Optional[List[str]] = None,
    ) -> int:
        """
        Purges cached status files, operational sample directories, and stale manifests.
        Preserves operational_files/logs/ to keep active execution logs intact.
        Returns the count of removed paths.
        """
        removed_count = 0
        target_set = set(target_sample_ids) if target_sample_ids else None

        # 1. Purge stage_status.json files in operational_files/cache/
        cache_dir = operational_files_dir / "cache"
        if cache_dir.is_dir():
            for status_file in list(cache_dir.glob("*.stage_status.json")):
                stem = status_file.name.removesuffix(".stage_status.json")
                s_id = stem.split("_")[0] if "_" in stem else stem
                if target_set is None or s_id in target_set or stem in target_set:
                    try:
                        status_file.unlink()
                        removed_count += 1
                    except OSError as exc:
                        logger.warning(f"Failed removing cache status file {status_file}: {exc}")

        # 2. Purge per-sample directories inside module folders
        for sub in STAGE_SUBDIRS:
            mod_dir = operational_files_dir / sub
            if not mod_dir.is_dir():
                continue
            for item in list(mod_dir.iterdir()):
                if item.is_dir():
                    if target_set is None or item.name in target_set:
                        try:
                            shutil.rmtree(item)
                            removed_count += 1
                        except OSError as exc:
                            logger.warning(f"Failed removing operational directory {item}: {exc}")
                elif item.is_file() and target_set is None:
                    try:
                        item.unlink()
                        removed_count += 1
                    except OSError:
                        pass

        # 3. Purge top-level per-sample manifest directories in operational_files/<sample_id>
        reserved_names = set(STAGE_SUBDIRS) | {"cache", "logs", "reference"}
        if operational_files_dir.is_dir():
            for item in list(operational_files_dir.iterdir()):
                if item.is_dir() and item.name not in reserved_names:
                    if target_set is None or item.name in target_set:
                        try:
                            shutil.rmtree(item)
                            removed_count += 1
                        except OSError as exc:
                            logger.warning(f"Failed removing sample manifest folder {item}: {exc}")

            op_manifest = operational_files_dir / "execution_manifest.json"
            if op_manifest.is_file():
                try:
                    op_manifest.unlink()
                    removed_count += 1
                except OSError:
                    pass

        # 4. Clean stale execution_manifest.json in results_dir if present
        if results_dir and results_dir.is_dir():
            res_manifest = results_dir / "execution_manifest.json"
            if res_manifest.is_file():
                try:
                    res_manifest.unlink()
                    removed_count += 1
                except OSError:
                    pass

        logger.info(f"Cache purge completed: {removed_count} cached artifact path(s) removed.")
        return removed_count

    def is_stage_valid(
        self,
        stage_name: str,
        expected_artifacts: List[Path],
        current_composite_hash: str,
    ) -> bool:
        """Evaluate whether a pipeline stage can be bypassed safely via cached deliverables."""
        if self.force_recompute:
            logger.info(
                f"[{self.sample_id}] Stage {stage_name} cache bypassed (force_recompute=True)."
            )
            return False

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

        payload = {
            "sample_id": self.sample_id,
            "sample_hash": self.sample_hash,
            "updated_at": datetime.now(timezone.utc).isoformat(),
            "stages": self.stages,
            "runtimes": self.runtimes,
            "stage_hashes": self.stage_hashes,
            "artifacts": self.artifacts,
            "metadata": self.metadata,
        }

        try:
            with open(self.status_file, "w", encoding="utf-8") as f:
                json.dump(payload, f, indent=2)

            # Also maintain canonical <sample_id>.stage_status.json
            canonical = self.status_file.parent / f"{self.sample_id}.stage_status.json"
            if self.status_file.resolve() != canonical.resolve():
                with open(canonical, "w", encoding="utf-8") as f:
                    json.dump(payload, f, indent=2)
        except Exception as exc:
            logger.warning(f"[{self.sample_id}] Failed saving status tracker: {exc}")