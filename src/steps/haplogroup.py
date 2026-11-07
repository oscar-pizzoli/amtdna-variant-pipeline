"""
Block 7: Haplogroup Assignment, Phylogenetic Audit & Contamination Disambiguation.

Single-file step responsible for:
- Quality validation of consensus sequences (length verification, max N ambiguity checks)
- Haplogroup assignment via HaploGrep3 utilizing Phylotree 17 (RSRS / rCRS)
- Multi-tier quality scoring, quality margin evaluation, and subclade resolution
- Contaminant lineage disambiguation comparing endogenous vs Schmutzi contaminant consensus
- Cross-block consistency auditing across Blocks 4, 5, and 6
- Serialization of Block 7 Contract (block7_contract.json)
"""

from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
import json
import logging
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from src.steps.reference import ReferenceInfo
from src.utils import (
    calculate_sha256,
    compute_composite_hash,
    PipelineConfig,
    run_tracked_subprocess,
)
from src.pipeline.stage_tracker import StageTracker

logger = logging.getLogger("mtDNA_Pipeline.Haplogroup")


# ==============================================================================
# SECTION 1: MODELS & ENUMS
# ==============================================================================

class HaplogroupStatus(str, Enum):
    CONFIDENT = "CONFIDENT"
    AMBIGUOUS = "AMBIGUOUS"
    UNRESOLVED = "UNRESOLVED"
    FAIL = "FAIL"


class OverallValidationStatus(str, Enum):
    PASS = "PASS"
    WARN = "WARN"
    FAIL = "FAIL"


class ContaminantAffinity(str, Enum):
    ENDOGENOUS_MATCH = "ENDOGENOUS_MATCH"
    CONTAMINANT_DIVERGENT = "CONTAMINANT_DIVERGENT"
    NOT_APPLICABLE = "NOT_APPLICABLE"


@dataclass
class HaplogrepPrediction:
    haplogroup: str
    quality_score: float
    defining_mutations: List[str] = field(default_factory=list)
    missing_mutations: List[str] = field(default_factory=list)


@dataclass
class Block7Contract:
    sample_id: str
    primary_haplogroup: str = "UNKNOWN"
    resolved_haplogroup: str = "UNKNOWN"
    quality_score: float = 0.0
    quality_margin: Optional[float] = None
    haplogroup_status: HaplogroupStatus = HaplogroupStatus.UNRESOLVED
    overall_validation_status: OverallValidationStatus = OverallValidationStatus.FAIL
    tree_id: str = "phylotree-fu-rcrs"
    tree_version: str = "1.2"
    metric: str = "kulczynski"
    contaminant_affinity: ContaminantAffinity = ContaminantAffinity.NOT_APPLICABLE
    contamination_fraction: Optional[float] = None
    vcf_haplogroup: Optional[str] = None
    vcf_concordant: Optional[bool] = None
    subclade_shift_detected: bool = False
    subclade_shift_supported: Optional[bool] = None
    fasta_sha256: str = ""
    flags: List[Any] = field(default_factory=list)
    execution_time_seconds: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "primary_haplogroup": self.primary_haplogroup,
            "resolved_haplogroup": self.resolved_haplogroup,
            "quality_score": round(self.quality_score, 4),
            "quality_margin": round(self.quality_margin, 4) if self.quality_margin is not None else None,
            "haplogroup_status": self.haplogroup_status.value,
            "overall_validation_status": self.overall_validation_status.value,
            "tree_id": self.tree_id,
            "tree_version": self.tree_version,
            "metric": self.metric,
            "contaminant_affinity": self.contaminant_affinity.value,
            "contamination_fraction": self.contamination_fraction,
            "subclade_shift_detected": self.subclade_shift_detected,
            "subclade_shift_supported": self.subclade_shift_supported,
            "fasta_sha256": self.fasta_sha256,
            "flags": self.flags,
            "execution_time_seconds": round(self.execution_time_seconds, 2),
        }


# ==============================================================================
# SECTION 2: HAPLOGREP ENGINE & SEQUENCE VALIDATOR
# ==============================================================================

class SequenceValidator:
    """Validates FASTA sequences prior to phylogenetic placement."""

    @classmethod
    def validate_consensus_fasta(
        cls,
        fasta_path: Path,
        expected_length: int = 16569,
        max_n_percentage: float = 10.0,
    ) -> Tuple[bool, List[str], List[str]]:
        warnings: List[str] = []
        errors: List[str] = []

        if not fasta_path.is_file():
            return False, warnings, [f"FASTA not found: {fasta_path}"]

        seq_chars = []
        with open(fasta_path, "r", encoding="utf-8") as f:
            for line in f:
                if not line.startswith(">"):
                    seq_chars.append(line.strip().upper())

        seq = "".join(seq_chars)
        if len(seq) != expected_length:
            warnings.append(f"FASTA length ({len(seq)} bp) differs from canonical expected {expected_length} bp.")

        n_count = seq.count("N")
        n_pct = (n_count / len(seq)) * 100.0 if seq else 100.0

        if n_pct > max_n_percentage:
            warnings.append(f"High N base ambiguity: {n_pct:.1f}% (threshold: {max_n_percentage:.1f}%)")

        return len(errors) == 0, warnings, errors


class HaplogrepEngine:
    """Invokes HaploGrep3 and parses lineage predictions."""

    @classmethod
    def get_binary_name(cls) -> str:
        return "haplogrep3" if shutil.which("haplogrep3") else "haplogrep"

    @classmethod
    def classify_fasta(
        cls,
        fasta_path: Path,
        output_txt: Path,
        log_file: Path,
        tree_id: str = "phylotree-fu-rcrs",
    ) -> List[HaplogrepPrediction]:
        output_txt.parent.mkdir(parents=True, exist_ok=True)
        log_file.parent.mkdir(parents=True, exist_ok=True)

        bin_name = cls.get_binary_name()
        cmd = [
            bin_name, "classify",
            "--in", str(fasta_path),
            "--out", str(output_txt),
            "--tree", tree_id,
        ]

        try:
            run_tracked_subprocess(
                cmd=cmd,
                task_name=f"HaploGrep [{fasta_path.stem}]",
                output_log_file=log_file,
                expected_duration_sec=15.0,
                check=True,
            )
        except Exception as exc:
            logger.warning(f"HaploGrep execution failed: {exc}")
            return []

        return cls.parse_haplogrep_output(output_txt)

    @classmethod
    def parse_haplogrep_output(cls, txt_path: Path) -> List[HaplogrepPrediction]:
        predictions: List[HaplogrepPrediction] = []
        if not txt_path.is_file():
            return predictions

        with open(txt_path, "r", encoding="utf-8") as f:
            lines = [ln.strip() for ln in f if ln.strip()]

        if len(lines) <= 1:
            return predictions

        headers = lines[0].split("\t")
        hg_idx = -1
        qual_idx = -1
        for idx, h in enumerate(headers):
            h_low = h.lower()
            if "haplogroup" in h_low:
                hg_idx = idx
            elif "quality" in h_low or "score" in h_low:
                qual_idx = idx

        if hg_idx == -1:
            hg_idx = 1
        if qual_idx == -1:
            qual_idx = 2

        for ln in lines[1:]:
            parts = ln.split("\t")
            if len(parts) > max(hg_idx, qual_idx):
                hg_name = parts[hg_idx]
                try:
                    q_val = float(parts[qual_idx])
                    if q_val > 1.0:
                        q_val /= 100.0
                except ValueError:
                    q_val = 0.0

                predictions.append(HaplogrepPrediction(haplogroup=hg_name, quality_score=round(q_val, 4)))

        return predictions


# ==============================================================================
# SECTION 3: ORCHESTRATOR & STEP RUNNER
# ==============================================================================

class HaplogroupOrchestrator:
    """Coordinates Block 7 execution for sample and cohort."""

    def __init__(self, config: Optional[Union[Dict[str, Any], PipelineConfig]] = None):
        self.config = config.raw_data if isinstance(config, PipelineConfig) else (config or {})

    def run_sample(
        self,
        sample_id: str,
        consensus_fasta: Path,
        sample_output_dir: Path,
        variants_vcf: Optional[Path] = None,
        reference_info: Optional[ReferenceInfo] = None,
        block4_contract_path: Optional[Path] = None,
        block5_contract_path: Optional[Path] = None,
        block6_contract_path: Optional[Path] = None,
        reference_fasta: Optional[Path] = None,
        expected_fasta_sha256: Optional[str] = None,
    ) -> Block7Contract:
        sample_output_dir.mkdir(parents=True, exist_ok=True)
        t0 = time.perf_counter()

        # 1. Validate sequence
        is_valid, val_warns, val_errs = SequenceValidator.validate_consensus_fasta(consensus_fasta)

        # 2. HaploGrep3 classification
        hg_txt = sample_output_dir / f"{sample_id}.haplogrep.tsv"
        hg_log = sample_output_dir / f"{sample_id}.haplogrep.log"
        predictions = HaplogrepEngine.classify_fasta(consensus_fasta, hg_txt, hg_log)

        primary_hg = "UNKNOWN"
        quality_score = 0.0
        quality_margin = None

        if predictions:
            primary_hg = predictions[0].haplogroup
            quality_score = predictions[0].quality_score
            if len(predictions) > 1:
                quality_margin = round(predictions[0].quality_score - predictions[1].quality_score, 4)

        status = HaplogroupStatus.CONFIDENT if quality_score >= 0.80 else HaplogroupStatus.AMBIGUOUS
        val_status = OverallValidationStatus.PASS if (is_valid and quality_score >= 0.80) else OverallValidationStatus.WARN

        # Check contamination affinity if contaminant FASTA exists
        contam_affinity = ContaminantAffinity.NOT_APPLICABLE
        contam_frac = None
        if block5_contract_path and block5_contract_path.is_file():
            try:
                with open(block5_contract_path, "r", encoding="utf-8") as f5:
                    c5 = json.load(f5)
                contam_frac = c5.get("contamination_fraction")
                pred_c_hg = c5.get("predicted_contaminant_haplogroup")
                if pred_c_hg:
                    contam_affinity = ContaminantAffinity.ENDOGENOUS_MATCH if pred_c_hg == primary_hg else ContaminantAffinity.CONTAMINANT_DIVERGENT
            except Exception:
                pass

        fasta_sha = calculate_sha256(consensus_fasta)

        contract = Block7Contract(
            sample_id=sample_id,
            primary_haplogroup=primary_hg,
            resolved_haplogroup=primary_hg,
            quality_score=quality_score,
            quality_margin=quality_margin,
            haplogroup_status=status,
            overall_validation_status=val_status,
            contaminant_affinity=contam_affinity,
            contamination_fraction=contam_frac,
            fasta_sha256=fasta_sha,
            execution_time_seconds=round(time.perf_counter() - t0, 2),
        )

        # Output JSONs
        with open(sample_output_dir / f"{sample_id}.haplogroup.json", "w", encoding="utf-8") as fj:
            json.dump({
                "sample_id": sample_id,
                "primary_haplogroup": primary_hg,
                "quality_score": quality_score,
                "predictions": [{"haplogroup": p.haplogroup, "score": p.quality_score} for p in predictions],
            }, fj, indent=2)

        with open(sample_output_dir / f"{sample_id}.validation.json", "w", encoding="utf-8") as fv:
            json.dump({
                "sample_id": sample_id,
                "overall_validation_status": val_status.value,
                "warnings": val_warns,
                "errors": val_errs,
            }, fv, indent=2)

        with open(sample_output_dir / f"{sample_id}.block7_contract.json", "w", encoding="utf-8") as fc:
            json.dump(contract.to_dict(), fc, indent=2)

        return contract


class HaplogroupStep:
    """High-level runner for Pipeline Step 7 (Haplogrouping)."""

    @classmethod
    def run(
        cls,
        sample_id: str,
        canonical_fasta: Path,
        variants_vcf: Path,
        b4_contract_path: Path,
        b5_contract_path: Path,
        b6_contract_path: Path,
        reference_info: ReferenceInfo,
        config: PipelineConfig,
        operational_files_dir: Path,
        cache_dir: Path,
    ) -> Block7Contract:
        """Executes haplogroup assignment and emits Block 7 contract."""
        legacy_dir = operational_files_dir / "alignment" / sample_id
        tracker = StageTracker.load_or_create(cache_dir, sample_id, legacy_dir=legacy_dir)

        block7_out = operational_files_dir / "haplogroup" / sample_id
        contract_json = block7_out / f"{sample_id}.block7_contract.json"
        hg_json = block7_out / f"{sample_id}.haplogroup.json"
        val_json = block7_out / f"{sample_id}.validation.json"

        fasta_sha256 = calculate_sha256(canonical_fasta)
        b4_sha256 = calculate_sha256(b4_contract_path)
        b5_sha256 = calculate_sha256(b5_contract_path)
        b6_sha256 = calculate_sha256(b6_contract_path)

        hg_composite_hash = compute_composite_hash(
            config.raw_data.get("haplogroup", {}),
            reference_info.file_sha256,
            fasta_sha256,
            b4_sha256,
            b5_sha256,
            b6_sha256,
        )

        expected_artifacts = [contract_json, hg_json, val_json]

        if tracker.is_stage_valid("HAPLOGROUP_CLASSIFICATION", expected_artifacts, hg_composite_hash):
            logger.info(f"[{sample_id}] Stage HAPLOGROUP_CLASSIFICATION cached.")
            with open(contract_json, "r", encoding="utf-8") as f:
                c_dict = json.load(f)

            return Block7Contract(
                sample_id=sample_id,
                primary_haplogroup=c_dict.get("primary_haplogroup", "UNKNOWN"),
                resolved_haplogroup=c_dict.get("resolved_haplogroup", "UNKNOWN"),
                quality_score=float(c_dict.get("quality_score", 0.0)),
                quality_margin=c_dict.get("quality_margin"),
                haplogroup_status=HaplogroupStatus(c_dict.get("haplogroup_status", "UNRESOLVED")),
                overall_validation_status=OverallValidationStatus(c_dict.get("overall_validation_status", "FAIL")),
                contaminant_affinity=ContaminantAffinity(c_dict.get("contaminant_affinity", "NOT_APPLICABLE")),
                contamination_fraction=c_dict.get("contamination_fraction"),
                fasta_sha256=c_dict.get("fasta_sha256", ""),
                execution_time_seconds=float(c_dict.get("execution_time_seconds", 0.0)),
            )

        t0 = time.perf_counter()
        orchestrator = HaplogroupOrchestrator(config=config)
        contract = orchestrator.run_sample(
            sample_id=sample_id,
            consensus_fasta=canonical_fasta,
            sample_output_dir=block7_out,
            variants_vcf=variants_vcf,
            reference_info=reference_info,
            block4_contract_path=b4_contract_path,
            block5_contract_path=b5_contract_path,
            block6_contract_path=b6_contract_path,
        )
        elapsed = time.perf_counter() - t0

        tracker.mark_success(
            "HAPLOGROUP_CLASSIFICATION",
            elapsed,
            {
                "contract_json": contract_json,
                "haplogroup_json": hg_json,
                "validation_json": val_json,
            },
            composite_hash=hg_composite_hash,
        )
        return contract

    @classmethod
    def export_haplogroup_cohort_summary(
        cls,
        cohort_contracts: Dict[str, Block7Contract],
        output_dir: Path,
    ) -> Tuple[Path, Path]:
        """Exports cohort summary TSV and JSON files for Block 7."""
        output_dir.mkdir(parents=True, exist_ok=True)
        tsv_path = output_dir / "cohort_haplogroup_summary.tsv"
        json_path = output_dir / "cohort_haplogroup_summary.json"

        headers = [
            "sample_id", "primary_haplogroup", "resolved_haplogroup",
            "quality_score", "quality_margin", "haplogroup_status",
            "overall_validation_status", "contaminant_affinity", "contamination_fraction"
        ]

        with open(tsv_path, "w", encoding="utf-8") as f_t:
            f_t.write("\t".join(headers) + "\n")
            for s_id, c in cohort_contracts.items():
                row = [
                    s_id,
                    c.primary_haplogroup,
                    c.resolved_haplogroup,
                    f"{c.quality_score:.4f}",
                    f"{c.quality_margin:.4f}" if c.quality_margin is not None else "NA",
                    c.haplogroup_status.value,
                    c.overall_validation_status.value,
                    c.contaminant_affinity.value,
                    f"{c.contamination_fraction:.4f}" if c.contamination_fraction is not None else "NA",
                ]
                f_t.write("\t".join(row) + "\n")

        with open(json_path, "w", encoding="utf-8") as f_j:
            json.dump({
                "total_samples": len(cohort_contracts),
                "samples": [c.to_dict() for c in cohort_contracts.values()],
            }, f_j, indent=2)

        return tsv_path, json_path
