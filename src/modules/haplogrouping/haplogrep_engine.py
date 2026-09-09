"""
Block 7 HaploGrep3 Execution Engine: CLI wrapper handling binary inspection,
versioned phylogeny parameterization, quality margin resolution, and strict metric bounding.
Supports native HaploGrep3 TSV output as well as legacy JSON payloads.
"""

import csv
import io
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
        metric_str = str(self.config.metric).upper()
        persisted_raw_path: Optional[str] = None

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_out = Path(tmp_dir) / f"{sample_id}_{input_type}.txt"
            cmd = [
                self._binary_path,
                "classify",
                "--in", str(input_path),
                "--out", str(tmp_out),
                "--tree", tree_spec,
                "--metric", metric_str,
                "--extend-report",
                "--hits", "5",
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

            if not tmp_out.is_file() or tmp_out.stat().st_size == 0:
                raise HaplogrepExecutionError(f"HaploGrep3 did not generate output for '{sample_id}'.")

            with open(tmp_out, "r", encoding="utf-8") as f:
                content = f.read()

            if not content.strip():
                raise HaplogrepExecutionError(f"Empty classification output from HaploGrep3 for '{sample_id}'.")

            if output_dir is not None and self.config.retain_raw_output:
                output_dir.mkdir(parents=True, exist_ok=True)
                dest_file = output_dir / f"{sample_id}.haplogrep.{input_type}.raw.txt"
                shutil.copyfile(tmp_out, dest_file)
                persisted_raw_path = str(dest_file)

        if content.strip().startswith(("{", "[")):
            try:
                payload = json.loads(content)
                return self._parse_json_payload(payload, sample_id, input_type, persisted_raw_path)
            except json.JSONDecodeError as err:
                raise HaplogrepExecutionError(f"Malformed JSON payload from HaploGrep3 for '{sample_id}': {err}")
        else:
            return self._parse_tsv_payload(content, sample_id, input_type, persisted_raw_path)

    @staticmethod
    def _clean_tsv_polys(raw_str: Any) -> List[str]:
        """Strip annotations like (globalPrivateMutation) and trailing exclamation marks."""
        if not raw_str or not isinstance(raw_str, str):
            return []
        cleaned_str = re.sub(r"\s*\([^)]*\)", "", raw_str)
        tokens = cleaned_str.strip().split()
        cleaned: List[str] = []
        for t in tokens:
            item = t.strip().rstrip("!")
            if item:
                cleaned.append(item)
        return cleaned

    def _parse_tsv_payload(
        self,
        tsv_content: str,
        sample_id: str,
        input_type: str,
        raw_path: Optional[str] = None,
    ) -> RawHaplogrepResult:
        """Parse native TSV output produced by HaploGrep3."""
        reader = csv.DictReader(io.StringIO(tsv_content), delimiter="\t")
        all_rows = list(reader)

        if not all_rows:
            raise HaplogrepExecutionError(f"No records found in HaploGrep3 TSV output for '{sample_id}'.")

        normalized_rows: List[Dict[str, str]] = []
        for r in all_rows:
            clean_r = {
                str(k).strip().strip('"'): (str(v).strip().strip('"') if v is not None else "")
                for k, v in r.items()
                if k is not None
            }
            normalized_rows.append(clean_r)

        sample_rows: List[Dict[str, str]] = []
        for r in normalized_rows:
            sid = r.get("SampleID") or r.get("Sample ID") or r.get("Sample") or r.get("sample", "")
            if sid == sample_id or len(normalized_rows) == 1:
                sample_rows.append(r)

        if not sample_rows:
            distinct_samples = {
                r.get("SampleID") or r.get("Sample ID") or r.get("Sample") or ""
                for r in normalized_rows
            }
            if len(distinct_samples) == 1:
                sample_rows = normalized_rows

        if not sample_rows:
            raise HaplogrepExecutionError(f"Sample '{sample_id}' not found in HaploGrep3 TSV output.")

        raw_candidates: List[Tuple[HaplogrepHit, Dict[str, str]]] = []
        for r in sample_rows:
            hg = r.get("Haplogroup", "")
            q_raw = r.get("Quality", "")
            if not hg or q_raw == "":
                continue
            q_score = self._validate_quality_score(q_raw, sample_id)
            raw_candidates.append((HaplogrepHit(haplogroup=hg, quality_score=q_score), r))

        if not raw_candidates:
            raise HaplogrepExecutionError(f"Sample '{sample_id}' lacks classification fields in HaploGrep3 TSV.")

        raw_candidates.sort(key=lambda pair: (-pair[0].quality_score, pair[0].haplogroup))

        best_hit, best_row = raw_candidates[0]
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

        found_muts = self._clean_tsv_polys(best_row.get("Found_Polys", ""))
        missing_muts = self._clean_tsv_polys(best_row.get("Not_Found_Polys", ""))
        all_private = self._clean_tsv_polys(best_row.get("Remaining_Polys", ""))

        private_muts = [m for m in all_private if not m.upper().endswith("N")]
        expected_muts = sorted(set(found_muts + missing_muts))

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
            mutation_parse_warnings=[],
            haplogrep_version=self.detected_version,
            tree_id=self.config.tree_id,
            tree_version=self.config.tree_version,
            metric=self.config.metric,
            raw_output_path=raw_path,
        )

    def _parse_json_payload(
        self,
        data: Any,
        sample_id: str,
        input_type: str,
        raw_path: Optional[str] = None,
    ) -> RawHaplogrepResult:
        """Backward-compatible parser for JSON payloads."""
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
        try:
            val = float(raw_value)
        except (ValueError, TypeError) as exc:
            raise HaplogrepExecutionError(
                f"Non-numeric quality score '{raw_value}' reported for '{sample_id}': {exc}"
            )

        if not (0.0 <= val <= 1.0):
            raise HaplogrepExecutionError(
                f"Quality score {val} for '{sample_id}' violates mathematical bounds "
                f"[0.0, 1.0] for metric '{self.config.metric}'."
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
                s = item.strip().rstrip("!")
                if s:
                    cleaned.append(s)
            elif isinstance(item, dict):
                variant = item.get("mutation") or item.get("variant") or item.get("pos")
                if variant:
                    cleaned.append(str(variant).strip().rstrip("!"))
                else:
                    warnings.append(f"Unrecognized dict structure in '{field_name}': {item}")
            else:
                warnings.append(f"Unexpected item type in '{field_name}': {type(item).__name__}")

        return cleaned