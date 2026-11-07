"""
Block 4: Ancient DNA Damage Quantification, Bayesian Rescaling & Authentication.

Single-file step responsible for:
- Quantification of post-mortem terminal deamination (5' C->T and 3' G->A transitions) via mapDamage2
- Ancient DNA read authentication and damage likelihood scoring via PMDtools (PMD scores >= 3)
- Library preparation treatment awareness (NON_UDG, UDG_HALF, UDG_FULL, SINGLE_STRANDED)
- Bayesian quality score rescaling downweighting damaged transition bases
- Forensic authenticity classification (AUTHENTIC_ANCIENT, MODERN_OR_NO_DAMAGE, INCONCLUSIVE)
- Serialization of Block 4 Contract (block4_contract.json)
"""

from dataclasses import dataclass, field
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
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pysam

from src.steps.reference import ReferenceInfo
from src.steps.alignment import AlignmentResult
from src.utils import (
    calculate_sha256,
    compute_composite_hash,
    PipelineConfig,
    run_tracked_subprocess,
)
from src.pipeline.stage_tracker import StageTracker

logger = logging.getLogger("mtDNA_Pipeline.Damage")


# ==============================================================================
# SECTION 1: MODELS & ENUMS
# ==============================================================================

class LibraryTreatment(str, Enum):
    NONE = "NONE"
    NON_UDG = "NON_UDG"
    UDG_HALF = "UDG_HALF"
    UDG_FULL = "UDG_FULL"
    SINGLE_STRANDED = "SINGLE_STRANDED"


class AuthenticationVerdict(str, Enum):
    AUTHENTIC_ANCIENT = "AUTHENTIC_ANCIENT"
    MODERN_OR_NO_DAMAGE = "MODERN_OR_NO_DAMAGE"
    INCONCLUSIVE = "INCONCLUSIVE"
    AUTHENTICITY_UNSUPPORTED = "AUTHENTICITY_UNSUPPORTED"


class DamageSignalLevel(str, Enum):
    NONE = "NONE"
    LOW = "LOW"
    MODERATE = "MODERATE"
    STRONG = "STRONG"


class DamageStatus(str, Enum):
    PASS = "PASS"
    WARN = "WARN"
    FAIL = "FAIL"


class DamageErrorCode(str, Enum):
    MAPDAMAGE_FAILED = "MAPDAMAGE_FAILED"
    PMDTOOLS_FAILED = "PMDTOOLS_FAILED"
    RESCALING_FAILED = "RESCALING_FAILED"
    BAM_CORRUPTED = "BAM_CORRUPTED"


class DamageWarningCode(str, Enum):
    LOW_DAMAGE_SIGNAL = "LOW_DAMAGE_SIGNAL"
    INSUFFICIENT_READS = "INSUFFICIENT_READS"
    FRAGMENT_LENGTH_ANOMALY = "FRAGMENT_LENGTH_ANOMALY"


@dataclass
class TerminalDamage:
    c_to_t_5p_pos1: Optional[float] = None
    c_to_t_5p_pos2: Optional[float] = None
    g_to_a_3p_pos1: Optional[float] = None
    g_to_a_3p_pos2: Optional[float] = None


@dataclass
class FragmentLengthMetrics:
    mean_length: Optional[float] = None
    median_length: Optional[float] = None


@dataclass
class PMDMetrics:
    total_reads_scored: int = 0
    min_pmd_score: Optional[float] = None
    max_pmd_score: Optional[float] = None
    mean_pmd_score: Optional[float] = None
    median_pmd_score: Optional[float] = None
    p75_pmd_score: Optional[float] = None
    fraction_pmd_ge_0: Optional[float] = None
    fraction_pmd_ge_3: Optional[float] = None
    fraction_pmd_ge_5: Optional[float] = None
    applied_threshold: Optional[float] = None
    reads_retained_post_pmd_filter: Optional[int] = None


@dataclass
class DamageProfile:
    terminal_damage: TerminalDamage = field(default_factory=TerminalDamage)
    fragment_length: FragmentLengthMetrics = field(default_factory=FragmentLengthMetrics)
    pmd_metrics: Optional[PMDMetrics] = None
    signal_level: DamageSignalLevel = DamageSignalLevel.NONE
    authentication_verdict: AuthenticationVerdict = AuthenticationVerdict.INCONCLUSIVE


@dataclass
class DamagePolicy:
    library_treatment: LibraryTreatment = LibraryTreatment.NONE
    enable_rescaling: bool = True
    enable_pmd_scoring: bool = True
    enable_pmd_filtering: bool = False
    pmd_threshold: float = 3.0
    min_c_to_t_5p_pos1: float = 0.03


@dataclass
class DamageOutputFiles:
    rescaled_bam: Optional[Path] = None
    damage_model_json: Optional[Path] = None
    misincorporation_plot: Optional[Path] = None
    length_distribution_plot: Optional[Path] = None


@dataclass
class DamageResult:
    sample_id: str
    status: DamageStatus = DamageStatus.PASS
    policy: DamagePolicy = field(default_factory=DamagePolicy)
    profile: DamageProfile = field(default_factory=DamageProfile)
    outputs: DamageOutputFiles = field(default_factory=DamageOutputFiles)
    damage_model_id: str = "model_ancient_v1"
    rescaled_bam: Optional[Path] = None
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "DamageResult":
        s_id = data.get("sample_id", "UNKNOWN")
        stat_str = data.get("status", "PASS")
        stat = DamageStatus(stat_str) if stat_str in [s.value for s in DamageStatus] else DamageStatus.PASS

        prof_dict = data.get("profile", {})
        td_dict = prof_dict.get("terminal_damage", {})
        fl_dict = prof_dict.get("fragment_length", {})
        pmd_dict = prof_dict.get("pmd_metrics")

        td = TerminalDamage(
            c_to_t_5p_pos1=td_dict.get("c_to_t_5p_pos1"),
            c_to_t_5p_pos2=td_dict.get("c_to_t_5p_pos2"),
            g_to_a_3p_pos1=td_dict.get("g_to_a_3p_pos1"),
            g_to_a_3p_pos2=td_dict.get("g_to_a_3p_pos2"),
        )
        fl = FragmentLengthMetrics(
            mean_length=fl_dict.get("mean_length"),
            median_length=fl_dict.get("median_length"),
        )
        pmd = None
        if pmd_dict:
            pmd = PMDMetrics(
                total_reads_scored=pmd_dict.get("total_reads_scored", 0),
                min_pmd_score=pmd_dict.get("min_pmd_score"),
                max_pmd_score=pmd_dict.get("max_pmd_score"),
                mean_pmd_score=pmd_dict.get("mean_pmd_score"),
                median_pmd_score=pmd_dict.get("median_pmd_score"),
                fraction_pmd_ge_0=pmd_dict.get("fraction_pmd_ge_0"),
                fraction_pmd_ge_3=pmd_dict.get("fraction_pmd_ge_3"),
                fraction_pmd_ge_5=pmd_dict.get("fraction_pmd_ge_5"),
            )

        verdict_str = prof_dict.get("authentication_verdict", "INCONCLUSIVE")
        verdict = AuthenticationVerdict(verdict_str) if verdict_str in [v.value for v in AuthenticationVerdict] else AuthenticationVerdict.INCONCLUSIVE

        profile = DamageProfile(
            terminal_damage=td,
            fragment_length=fl,
            pmd_metrics=pmd,
            authentication_verdict=verdict,
        )

        outs_dict = data.get("outputs", {})
        rescaled_p = Path(outs_dict["rescaled_bam"]) if outs_dict.get("rescaled_bam") else None
        outputs = DamageOutputFiles(rescaled_bam=rescaled_p)

        return cls(
            sample_id=s_id,
            status=stat,
            profile=profile,
            outputs=outputs,
            rescaled_bam=rescaled_p,
            errors=data.get("errors", []),
            warnings=data.get("warnings", []),
        )


# ==============================================================================
# SECTION 2: MAPDAMAGE2 QUANTIFICATION ENGINE
# ==============================================================================

class DamageQuantifierEngine:
    """Runs mapDamage2 to quantify terminal nucleotide misincorporation frequencies."""

    @classmethod
    def run_mapdamage(
        cls,
        bam_path: Path,
        reference_fasta: Path,
        output_dir: Path,
        log_file: Path,
        threads: int = 2,
    ) -> Tuple[TerminalDamage, FragmentLengthMetrics]:
        output_dir.mkdir(parents=True, exist_ok=True)
        log_file.parent.mkdir(parents=True, exist_ok=True)

        cmd = [
            "mapDamage",
            "-i", str(bam_path),
            "-r", str(reference_fasta),
            "-d", str(output_dir),
            "--no-rescale",
        ]

        run_tracked_subprocess(
            cmd=cmd,
            task_name=f"mapDamage [{bam_path.stem}]",
            output_log_file=log_file,
            expected_duration_sec=35.0,
            check=True,
        )

        misinc_txt = output_dir / "misincorporation.txt"
        td = cls.parse_misincorporation_table(misinc_txt)
        fl = cls.parse_length_distribution(output_dir / "lgdistribution.txt")
        return td, fl

    @classmethod
    def parse_misincorporation_table(cls, table_path: Path) -> TerminalDamage:
        td = TerminalDamage(c_to_t_5p_pos1=0.0, c_to_t_5p_pos2=0.0, g_to_a_3p_pos1=0.0, g_to_a_3p_pos2=0.0)
        if not table_path.is_file():
            return td

        with open(table_path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 6:
                    end, pos, ref, base, count, total = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]
                    try:
                        p_int = int(pos)
                        tot_fl = float(total)
                        frac = float(count) / tot_fl if tot_fl > 0 else 0.0
                    except ValueError:
                        continue

                    if end == "5p" and ref == "C" and base == "T":
                        if p_int == 1:
                            td.c_to_t_5p_pos1 = round(frac, 4)
                        elif p_int == 2:
                            td.c_to_t_5p_pos2 = round(frac, 4)
                    elif end == "3p" and ref == "G" and base == "A":
                        if p_int == 1:
                            td.g_to_a_3p_pos1 = round(frac, 4)
                        elif p_int == 2:
                            td.g_to_a_3p_pos2 = round(frac, 4)

        return td

    @classmethod
    def parse_length_distribution(cls, lg_path: Path) -> FragmentLengthMetrics:
        if not lg_path.is_file():
            return FragmentLengthMetrics()

        lengths: List[int] = []
        with open(lg_path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        l_val, count = int(parts[0]), int(parts[1])
                        lengths.extend([l_val] * min(count, 5000))
                    except ValueError:
                        pass

        if not lengths:
            return FragmentLengthMetrics()

        return FragmentLengthMetrics(
            mean_length=round(float(np.mean(lengths)), 1),
            median_length=round(float(np.median(lengths)), 1),
        )


# ==============================================================================
# SECTION 3: PMDTOOLS SCORING ENGINE
# ==============================================================================

class PMDScoringEngine:
    """Runs PMDtools and parses log-likelihood damage scores."""

    @classmethod
    def run_pmdtools(
        cls,
        bam_path: Path,
        reference_fasta: Path,
        output_tsv: Path,
        log_file: Path,
        threshold: float = 3.0,
    ) -> PMDMetrics:
        output_tsv.parent.mkdir(parents=True, exist_ok=True)
        log_file.parent.mkdir(parents=True, exist_ok=True)

        cmd = [
            "pmdtools",
            "--print",
            "--threshold", str(threshold),
            "--reference", str(reference_fasta),
        ]

        with open(bam_path, "rb") as f_in, open(output_tsv, "w", encoding="utf-8") as f_out, open(log_file, "w", encoding="utf-8") as f_err:
            p = subprocess.Popen(cmd, stdin=f_in, stdout=f_out, stderr=f_err)
            p.wait()

        return cls.parse_pmd_scores_table(output_tsv)

    @classmethod
    def parse_pmd_scores_table(cls, scores_tsv: Path) -> PMDMetrics:
        """Parses PMD scoring table, supporting header detection, comments, and column indexing."""
        if not scores_tsv.is_file() or scores_tsv.stat().st_size == 0:
            return PMDMetrics()

        scores: List[float] = []
        score_col_idx: Optional[int] = None

        with open(scores_tsv, "r", encoding="utf-8") as f:
            for line in f:
                line_str = line.strip()
                if not line_str or line_str.startswith("#"):
                    continue

                parts = line_str.split()
                if not parts:
                    continue

                # Header detection: identify if line contains 'score'
                if score_col_idx is None and any("score" in p.lower() for p in parts):
                    indices = [i for i, p in enumerate(parts) if "score" in p.lower()]
                    score_col_idx = indices[0]
                    continue

                val: Optional[float] = None
                if score_col_idx is not None and score_col_idx < len(parts):
                    try:
                        val = float(parts[score_col_idx])
                    except ValueError:
                        pass

                if val is None:
                    # Fallback heuristics: try last column or second column
                    try:
                        val = float(parts[-1])
                    except ValueError:
                        if len(parts) >= 2:
                            try:
                                val = float(parts[1])
                            except ValueError:
                                pass

                if val is not None:
                    scores.append(val)

        if not scores:
            return PMDMetrics()

        arr = np.array(scores, dtype=np.float64)
        total = len(arr)

        return PMDMetrics(
            total_reads_scored=total,
            min_pmd_score=float(np.min(arr)),
            max_pmd_score=float(np.max(arr)),
            mean_pmd_score=float(np.mean(arr)),
            median_pmd_score=float(np.median(arr)),
            fraction_pmd_ge_0=round(float(np.sum(arr >= 0.0) / total), 4),
            fraction_pmd_ge_3=round(float(np.sum(arr >= 3.0) / total), 4),
            fraction_pmd_ge_5=round(float(np.sum(arr >= 5.0) / total), 4),
        )


# ==============================================================================
# SECTION 4: RESCALING ENGINE
# ==============================================================================

class QualityRescaleEngine:
    """Executes mapDamage2 Bayesian quality score rescaling."""

    @classmethod
    def run_rescaling(
        cls,
        bam_path: Path,
        reference_fasta: Path,
        output_dir: Path,
        log_file: Path,
    ) -> Path:
        output_dir.mkdir(parents=True, exist_ok=True)
        log_file.parent.mkdir(parents=True, exist_ok=True)

        cmd = [
            "mapDamage",
            "-i", str(bam_path),
            "-r", str(reference_fasta),
            "-d", str(output_dir),
            "--rescale",
        ]

        run_tracked_subprocess(
            cmd=cmd,
            task_name=f"mapDamage Rescale [{bam_path.stem}]",
            output_log_file=log_file,
            expected_duration_sec=60.0,
            check=True,
        )

        rescaled_bams = list(output_dir.glob("*.rescaled.bam"))
        if not rescaled_bams:
            raise RuntimeError(f"mapDamage completed but no *.rescaled.bam found in {output_dir}")

        rescaled_bam = rescaled_bams[0]
        cmd_idx = ["samtools", "index", str(rescaled_bam)]
        subprocess.run(cmd_idx, check=True, capture_output=True)

        return rescaled_bam


# ==============================================================================
# SECTION 5: QC EVALUATOR & ORCHESTRATOR
# ==============================================================================

class DamageQCValidator:
    """Evaluates terminal damage and PMD scoring to determine biological authenticity."""

    @classmethod
    def evaluate_authenticity(
        cls,
        td: TerminalDamage,
        pmd: Optional[PMDMetrics],
        policy: DamagePolicy,
    ) -> Tuple[AuthenticationVerdict, DamageSignalLevel]:
        c1 = td.c_to_t_5p_pos1 or 0.0
        g1 = td.g_to_a_3p_pos1 or 0.0
        pmd3 = pmd.fraction_pmd_ge_3 if pmd and pmd.fraction_pmd_ge_3 is not None else 0.0

        if policy.library_treatment in (LibraryTreatment.UDG_FULL,):
            verdict = AuthenticationVerdict.AUTHENTIC_ANCIENT
            level = DamageSignalLevel.LOW
        elif c1 >= 0.05 or g1 >= 0.05 or pmd3 >= 0.10:
            verdict = AuthenticationVerdict.AUTHENTIC_ANCIENT
            level = DamageSignalLevel.STRONG if (c1 >= 0.10 or pmd3 >= 0.20) else DamageSignalLevel.MODERATE
        elif c1 >= 0.03 or g1 >= 0.03:
            verdict = AuthenticationVerdict.AUTHENTIC_ANCIENT
            level = DamageSignalLevel.LOW
        else:
            verdict = AuthenticationVerdict.MODERN_OR_NO_DAMAGE
            level = DamageSignalLevel.NONE

        return verdict, level


class DamageOrchestrator:
    """Full coordinator for Block 4 damage analysis."""

    @classmethod
    def run_sample_damage_pipeline(
        cls,
        sample_id: str,
        reference_info: ReferenceInfo,
        output_dir: Path,
        log_dir: Path,
        policy: Optional[DamagePolicy] = None,
        dedup_bam: Optional[Path] = None,
        markdup_bam: Optional[Path] = None,
        threads: int = 2,
        allow_input_bam_fallback: bool = True,
    ) -> DamageResult:
        active_policy = policy or DamagePolicy()
        sample_out = output_dir / sample_id
        sample_log = log_dir / sample_id

        sample_out.mkdir(parents=True, exist_ok=True)
        sample_log.mkdir(parents=True, exist_ok=True)

        target_bam = dedup_bam or markdup_bam
        if target_bam is None or not target_bam.is_file():
            raise RuntimeError(f"Input BAM not found for {sample_id}")

        # 1. mapDamage quantification
        md_dir = sample_out / "mapdamage"
        md_log = sample_log / f"{sample_id}.mapdamage.log"
        td, fl = DamageQuantifierEngine.run_mapdamage(
            target_bam, reference_info.fasta_path, md_dir, md_log, threads=threads
        )

        # 2. PMDtools
        pmd_metrics = None
        if active_policy.enable_pmd_scoring:
            pmd_tsv = sample_out / f"{sample_id}.pmd.tsv"
            pmd_log = sample_log / f"{sample_id}.pmdtools.log"
            try:
                pmd_metrics = PMDScoringEngine.run_pmdtools(
                    target_bam, reference_info.fasta_path, pmd_tsv, pmd_log, threshold=active_policy.pmd_threshold
                )
            except Exception as e:
                logger.warning(f"[{sample_id}] PMDtools failed: {e}")

        # 3. Rescaling
        rescaled_bam = None
        if active_policy.enable_rescaling:
            resc_dir = sample_out / "rescaling"
            resc_log = sample_log / f"{sample_id}.rescale.log"
            try:
                rescaled_bam = QualityRescaleEngine.run_rescaling(
                    target_bam, reference_info.fasta_path, resc_dir, resc_log
                )
            except Exception as e:
                logger.warning(f"[{sample_id}] mapDamage rescaling failed: {e}")

        # 4. Authenticity assessment
        verdict, level = DamageQCValidator.evaluate_authenticity(td, pmd_metrics, active_policy)

        prof = DamageProfile(
            terminal_damage=td,
            fragment_length=fl,
            pmd_metrics=pmd_metrics,
            signal_level=level,
            authentication_verdict=verdict,
        )

        model_json = sample_out / f"{sample_id}.damage_model.json"
        with open(model_json, "w", encoding="utf-8") as f:
            json.dump({
                "sample_id": sample_id,
                "authentication_verdict": verdict.value,
                "terminal_damage": {
                    "c_to_t_5p_pos1": td.c_to_t_5p_pos1,
                    "g_to_a_3p_pos1": td.g_to_a_3p_pos1,
                },
                "fragment_length": {
                    "mean_length_bp": fl.mean_length,
                    "median_length_bp": fl.median_length,
                },
                "pmd_metrics": {
                    "fraction_pmd_ge_3": pmd_metrics.fraction_pmd_ge_3 if pmd_metrics else None,
                },
            }, f, indent=2)

        return DamageResult(
            sample_id=sample_id,
            status=DamageStatus.PASS,
            policy=active_policy,
            profile=prof,
            outputs=DamageOutputFiles(rescaled_bam=rescaled_bam, damage_model_json=model_json),
            rescaled_bam=rescaled_bam,
        )


# ==============================================================================
# SECTION 6: STEP RUNNER & CONTRACT BUILDER (Block 4)
# ==============================================================================

class DamageStep:
    """High-level runner for Pipeline Step 4 (Damage Assessment & Rescaling)."""

    @classmethod
    def run(
        cls,
        sample_id: str,
        aln_res: AlignmentResult,
        reference_info: ReferenceInfo,
        config: PipelineConfig,
        operational_files_dir: Path,
        cache_dir: Path,
        logs_dir: Path,
        threads: int = 2,
    ) -> Tuple[DamageResult, Path]:
        """Runs damage assessment, checks cache validity, and emits Block 4 contract."""
        legacy_dir = operational_files_dir / "alignment" / sample_id
        tracker = StageTracker.load_or_create(cache_dir, sample_id, legacy_dir=legacy_dir)

        damage_dir = operational_files_dir / "damage" / sample_id
        damage_json = damage_dir / f"{sample_id}.damage.json"
        damage_model_json = damage_dir / f"{sample_id}.damage_model.json"
        b4_contract_json = damage_dir / f"{sample_id}.block4_contract.json"

        input_bam = None
        if aln_res.outputs:
            if aln_res.outputs.dedup_bam and aln_res.outputs.dedup_bam.is_file():
                input_bam = aln_res.outputs.dedup_bam
            elif aln_res.outputs.filtered_bam and aln_res.outputs.filtered_bam.is_file():
                input_bam = aln_res.outputs.filtered_bam

        if input_bam is None:
            cand = legacy_dir / f"{sample_id}.filtered.dedup.bam"
            if cand.is_file():
                input_bam = cand
            else:
                raise RuntimeError(f"Physical filtered dedup BAM not found for {sample_id}")

        input_bam_sha256 = calculate_sha256(input_bam)
        damage_cfg = config.raw_data.get("damage", {})
        damage_composite_hash = compute_composite_hash(
            damage_cfg,
            reference_info.file_sha256,
            input_bam_sha256,
        )

        active_cache_json = damage_json if damage_json.is_file() else damage_model_json

        if tracker.is_stage_valid("DAMAGE_AUTHENTICATION", [active_cache_json, b4_contract_json], damage_composite_hash):
            logger.info(f"[{sample_id}] Stage DAMAGE_AUTHENTICATION cached.")
            with open(active_cache_json, "r", encoding="utf-8") as f:
                dam_res = DamageResult.from_dict(json.load(f))
            return dam_res, b4_contract_json

        t0 = time.perf_counter()
        raw_treatment = damage_cfg.get("library_treatment", "NONE")
        try:
            damage_lib_treatment = LibraryTreatment(raw_treatment)
        except ValueError:
            damage_lib_treatment = LibraryTreatment.NONE

        damage_policy = DamagePolicy(
            library_treatment=damage_lib_treatment,
            enable_rescaling=bool(damage_cfg.get("enable_rescaling", True)),
            enable_pmd_scoring=bool(damage_cfg.get("enable_pmd_scoring", True)),
            enable_pmd_filtering=bool(damage_cfg.get("enable_pmd_filtering", False)),
            pmd_threshold=float(damage_cfg.get("pmd_threshold", 3.0)),
        )

        dam_res = DamageOrchestrator.run_sample_damage_pipeline(
            sample_id=sample_id,
            reference_info=reference_info,
            output_dir=operational_files_dir / "damage",
            log_dir=logs_dir / "damage",
            policy=damage_policy,
            dedup_bam=input_bam,
            markdup_bam=aln_res.outputs.markdup_bam if aln_res.outputs else None,
            threads=threads,
        )

        cls.ensure_block4_contract(
            b4_contract_json, sample_id, dam_res, input_bam, reference_info.file_sha256, force_overwrite=True
        )

        target_model_file = damage_json if damage_json.is_file() else (dam_res.outputs.damage_model_json or damage_model_json)
        tracker.mark_success(
            "DAMAGE_AUTHENTICATION",
            time.perf_counter() - t0,
            {"damage_json": target_model_file, "b4_contract": b4_contract_json},
            composite_hash=damage_composite_hash,
        )
        return dam_res, b4_contract_json

    @classmethod
    def ensure_block4_contract(
        cls,
        contract_path: Path,
        sample_id: str,
        dam_res: DamageResult,
        input_bam: Path,
        ref_sha256: str,
        force_overwrite: bool = False,
    ) -> None:
        """Serializes the Block 4 Contract (block4_contract.json)."""
        contract_path.parent.mkdir(parents=True, exist_ok=True)
        if not force_overwrite and contract_path.is_file() and contract_path.stat().st_size > 0:
            return

        rescaling_performed = False
        rescaled_bam_str = None
        if dam_res.rescaled_bam and Path(dam_res.rescaled_bam).is_file():
            rescaling_performed = True
            rescaled_bam_str = str(dam_res.rescaled_bam)
        elif dam_res.outputs and dam_res.outputs.rescaled_bam and Path(dam_res.outputs.rescaled_bam).is_file():
            rescaling_performed = True
            rescaled_bam_str = str(dam_res.outputs.rescaled_bam)

        auth_status = "AUTHENTICITY_UNSUPPORTED"
        prof = dam_res.profile
        if prof and prof.authentication_verdict:
            auth_status = prof.authentication_verdict.value

        rescaled_sha = calculate_sha256(rescaled_bam_str) if (rescaling_performed and rescaled_bam_str) else None
        filtered_sha = calculate_sha256(input_bam)

        fl_med = prof.fragment_length.median_length if prof else None
        fl_mean = prof.fragment_length.mean_length if prof else None
        pmd_ge3 = prof.pmd_metrics.fraction_pmd_ge_3 if (prof and prof.pmd_metrics) else None
        c1 = prof.terminal_damage.c_to_t_5p_pos1 if prof else None
        g1 = prof.terminal_damage.g_to_a_3p_pos1 if prof else None

        contract_data = {
            "sample_id": sample_id,
            "authentication_status": auth_status,
            "rescaling_performed": rescaling_performed,
            "rescaled_bam": rescaled_bam_str,
            "rescaled_bam_sha256": rescaled_sha,
            "filtered_dedup_bam": str(input_bam),
            "filtered_dedup_bam_sha256": filtered_sha,
            "damage_model_id": getattr(dam_res, "damage_model_id", "model_ancient_v1"),
            "library_treatment": dam_res.policy.library_treatment.value if dam_res.policy else "NONE",
            "reference_sha256": ref_sha256,
            "median_fragment_length": round(fl_med, 1) if fl_med is not None else None,
            "mean_fragment_length": round(fl_mean, 1) if fl_mean is not None else None,
            "pmd_fraction_ge_3": round(pmd_ge3, 4) if pmd_ge3 is not None else None,
            "c_to_t_5p_rate_pos1": round(c1, 4) if c1 is not None else None,
            "g_to_a_3p_rate_pos1": round(g1, 4) if g1 is not None else None,
        }

        with open(contract_path, "w", encoding="utf-8") as f:
            json.dump(contract_data, f, indent=2)

    @classmethod
    def extract_verdict_with_contract_fallback(cls, dam_res: DamageResult, b4_contract_path: Path) -> str:
        """Extracts authenticity verdict prioritizing formal contract fallback."""
        verdict = "AUTHENTICITY_UNSUPPORTED"
        if dam_res.profile and dam_res.profile.authentication_verdict:
            verdict = dam_res.profile.authentication_verdict.value

        if verdict.upper() in ("INCONCLUSIVE", "AUTHENTICITY_UNSUPPORTED", "PASS", "WARN") and b4_contract_path.is_file():
            try:
                with open(b4_contract_path, "r", encoding="utf-8") as f_c:
                    c_data = json.load(f_c)
                contract_auth = c_data.get("authentication_status")
                if contract_auth and str(contract_auth).upper() not in ("INCONCLUSIVE", "AUTHENTICITY_UNSUPPORTED"):
                    verdict = str(contract_auth)
            except Exception:
                pass

        return verdict
