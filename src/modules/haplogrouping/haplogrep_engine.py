"""
Block 7 HaploGrep3 Execution Engine: CLI wrapper handling binary inspection,
versioned phylogeny parameterization, quality margin resolution, and strict metric bounding.
"""

import json
import logging
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from src.modules.haplogrouping.models import (
    HaplogrepHit,
    HaplogrepRunConfig,
    RawHaplogrepResult,
)

logger = logging.getLogger("mtDNA_Pipeline.HaplogrepEngine")


class HaplogrepExecutionError(Exception):
    """Raised when HaploGrep3 execution, schema verification, or output parsing fails."""
    pass


class HaplogrepEngine:
    """Executes HaploGrep3 classification over FASTA and VCF artifacts with strict bounds checking."""

    def __init__(self, run_config: HaplogrepRunConfig):
        self.config = run_config
        self._binary_path = self._resolve_binary(run_config.tool)
        self.detected_version = self.detect_version()

    @staticmethod
    def _resolve_binary(tool_name: str) -> str:
        binary = shutil.which(tool_name)
        if not binary:
            raise HaplogrepExecutionError(f"Binary '{tool_name}' not found in PATH.")
        return binary

    def detect_version(self) -> str:
        try:
            res = subprocess.run(
                [self._binary_path, "--version"],
                capture_output=True,
                text=True,
                check=True,
                timeout=15,
            )
            raw_out = res.stdout.strip() or res.stderr.strip()
            match = re.search(r"(\d+\.\d+(\.\d+)?)", raw_out)
            return match.group(1) if match else raw_out
        except Exception as err:
            logger.warning(f"Could not determine haplogrep3 version: {err}")
            return "unknown"

    def classify_fasta(
        self,
        fasta_path: Path,
        sample_id: str,
        output_dir: Optional[Path] = None,
    ) -> RawHaplogrepResult:
        return self._run_classify(fasta_path, sample_id, "fasta", output_dir)

    def classify_vcf(
        self,
        vcf_path: Path,
        sample_id: str,
        output_dir: Optional[Path] = None,
    ) -> RawHaplogrepResult:
        return self._run_classify(vcf_path, sample_id, "vcf", output_dir)

    def _run_classify(
        self,
        input_path: Path,
        sample_id: str,
        input_type: str,
        output_dir: Optional[Path],
    ) -> RawHaplogrepResult:
        if not input_path.is_file():
            raise FileNotFoundError(f"Input {input_type.upper()} file not found: {input_path}")

        tree_spec = f"{self.config.tree_id}@{self.config.tree_version}"
        persisted_raw_path: Optional[str] = None

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_json = Path(tmp_dir) / f"{sample_id}_{input_type}.json"
            cmd = [
                self._binary_path,
                "classify",
                "--in", str(input_path),
                "--out", str(tmp_json),
                "--tree", tree_spec,
                "--metric", self.config.metric,
            ]

            logger.info(f"[{sample_id}] Executing HaploGrep3 ({input_type}): {' '.join(cmd)}")
            try:
                res = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=True,
                    timeout=self.config.timeout_seconds,
                )
            except subprocess.TimeoutExpired as err:
                raise HaplogrepExecutionError(
                    f"HaploGrep3 {input_type} timed out for '{sample_id}' after {self.config.timeout_seconds}s."
                ) from err
            except subprocess.CalledProcessError as err:
                raise HaplogrepExecutionError(
                    f"HaploGrep3 {input_type} failed for '{sample_id}' (exit {err.returncode}):\n{err.stderr}"
                )

            if not tmp_json.is_file() or tmp_json.stat().st_size == 0:
                raise HaplogrepExecutionError(f"HaploGrep3 did not generate output JSON for '{sample_id}'.")

            with open(tmp_json, "r", encoding="utf-8") as f:
                try:
                    payload = json.load(f)
                except json.JSONDecodeError as err:
                    raise HaplogrepExecutionError(f"Malformed JSON from HaploGrep3 for '{sample_id}': {err}")

            if output_dir is not None and self.config.retain_raw_output:
                output_dir.mkdir(parents=True, exist_ok=True)
                dest_json = output_dir / f"{sample_id}.haplogrep.{input_type}.raw.json"
                shutil.copyfile(tmp_json, dest_json)
                persisted_raw_path = str(dest_json)

        return self._parse_json_payload(payload, sample_id, input_type, persisted_raw_path)

    def _parse_json_payload(
        self,
        data: Any,
        sample_id: str,
        input_type: str,
        raw_path: Optional[str] = None,
    ) -> RawHaplogrepResult:
        records = data if isinstance(data, list) else data.get("samples", [data])
        if not records or not isinstance(records, list):
            raise HaplogrepExecutionError(f"No sample records found in HaploGrep3 payload for '{sample_id}'.")

        sample_entry = None
        for entry in records:
            if not isinstance(entry, dict):
                continue
            entry_id = entry.get("sample", entry.get("id", entry.get("sampleId", "")))
            if entry_id == sample_id or len(records) == 1:
                sample_entry = entry
                break

        if not sample_entry:
            raise HaplogrepExecutionError(f"Sample '{sample_id}' not found in HaploGrep3 output records.")

        raw_candidates: List[Tuple[HaplogrepHit, Dict[str, Any]]] = []
        hits_raw = sample_entry.get("hits", [])
        if isinstance(hits_raw, list) and hits_raw:
            for h in hits_raw:
                if not isinstance(h, dict):
                    continue
                hg = h.get("haplogroup")
                q_val = h.get("quality")
                if hg is None or q_val is None:
                    continue
                q_score = self._validate_quality_score(q_val, sample_id)
                raw_candidates.append(
                    (HaplogrepHit(haplogroup=str(hg).strip(), quality_score=q_score), h)
                )

        if not raw_candidates:
            hg = sample_entry.get("haplogroup")
            q_val = sample_entry.get("quality")
            if hg is None or q_val is None or not str(hg).strip():
                raise HaplogrepExecutionError(f"Sample '{sample_id}' lacks classification fields.")
            q_score = self._validate_quality_score(q_val, sample_id)
            raw_candidates.append(
                (HaplogrepHit(haplogroup=str(hg).strip(), quality_score=q_score), sample_entry)
            )

        # Deterministic sorting: quality descending, then alphabetical tie-breaker
        raw_candidates.sort(key=lambda pair: (-pair[0].quality_score, pair[0].haplogroup))

        best_hit, best_raw_dict = raw_candidates[0]
        max_quality = best_hit.quality_score

        second_quality: Optional[float] = None
        quality_margin: Optional[float] = None
        if len(raw_candidates) > 1:
            second_quality = raw_candidates[1][0].quality_score
            quality_margin = round(max_quality - second_quality, 4)

        top_candidates = [
            hit for hit, _ in raw_candidates
            if abs(hit.quality_score - max_quality) < 1e-6
        ]
        co_optimal = len(top_candidates) > 1

        parse_warnings: List[str] = []
        expected_muts = self._clean_mutations(best_raw_dict.get("expected", []), "expected", parse_warnings)
        found_muts = self._clean_mutations(best_raw_dict.get("found", []), "found", parse_warnings)
        missing_muts = self._clean_mutations(best_raw_dict.get("missing", []), "missing", parse_warnings)
        private_muts = self._clean_mutations(best_raw_dict.get("private", []), "private", parse_warnings)

        return RawHaplogrepResult(
            sample_id=sample_id,
            input_type=input_type,
            haplogroup=best_hit.haplogroup,
            quality_score=max_quality,
            second_best_quality=second_quality,
            quality_margin=quality_margin,
            co_optimal=co_optimal,
            top_hits=top_candidates,
            expected_mutations=expected_muts,
            found_mutations=found_muts,
            missing_mutations=missing_muts,
            private_mutations=private_muts,
            mutation_parse_warnings=parse_warnings,
            haplogrep_version=self.detected_version,
            tree_id=self.config.tree_id,
            tree_version=self.config.tree_version,
            metric=self.config.metric,
            raw_output_path=raw_path,
        )

    def _validate_quality_score(self, raw_value: Any, sample_id: str) -> float:
        """
        Strictly validate quality score based on the configured metric without heuristic scaling.
        Kulczynski metric mathematically spans [0.0, 1.0].
        """
        try:
            val = float(raw_value)
        except (ValueError, TypeError) as exc:
            raise HaplogrepExecutionError(
                f"Non-numeric quality score '{raw_value}' reported for '{sample_id}': {exc}"
            )

        if not (0.0 <= val <= 1.0):
            raise HaplogrepExecutionError(
                f"Quality score {val} for '{sample_id}' violates mathematical bounds "
                f"[0.0, 1.0] for metric '{self.config.metric}'. Tool failure or schema corruption."
            )

        return round(val, 4)

    @staticmethod
    def _clean_mutations(mut_list: Any, field_name: str, warnings: List[str]) -> List[str]:
        cleaned: List[str] = []
        if not isinstance(mut_list, list):
            if mut_list is not None:
                warnings.append(f"Field '{field_name}' is not a list: {type(mut_list).__name__}")
            return cleaned

        for item in mut_list:
            if isinstance(item, str):
                s = item.strip()
                if s:
                    cleaned.append(s)
            elif isinstance(item, dict):
                variant = item.get("mutation") or item.get("variant") or item.get("pos")
                if variant:
                    cleaned.append(str(variant).strip())
                else:
                    warnings.append(f"Unrecognized dict structure in '{field_name}': {item}")
            else:
                warnings.append(f"Unexpected item type in '{field_name}': {type(item).__name__}")

        return cleaned