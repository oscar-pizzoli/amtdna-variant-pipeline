"""Track A Engine: schmutzi execution via Docker container or local binary.

Bayesian joint mixture modeling, output parsing, diagnostic extraction,
and adaptive zero-damage / modern DNA handling.
"""

from dataclasses import dataclass
import hashlib
import logging
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
from typing import Any, Dict, List, Optional, Tuple
import uuid

import pysam

from src.modules.contamination.models import (
    ContaminationErrorCode,
    ContaminationEstimate,
    ContaminationMethod,
    ContaminationOutputFiles,
    ContaminationWarningCode,
    EstimationConfidence,
    LibraryTreatment,
    ModelQualityDiagnostics,
    SchmutziPolicy,
    WrappingMetadata,
)

logger = logging.getLogger("mtDNA_Pipeline.SchmutziEngine")

NUMBER_PATTERN = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


class SchmutziExecutionError(Exception):
    """Raised when schmutzi binaries fail, timeout, or output artifacts cannot be parsed."""
    pass


@dataclass
class SchmutziExecutionMetadata:
    """Audit metadata tracking schmutzi command execution, versions, and file hashes."""

    tool_version: str = "unknown"
    executed_command: str = ""
    contaminant_database_used: Optional[Path] = None
    contaminant_database_sha256: Optional[str] = None
    runtime_seconds: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "tool_version": self.tool_version,
            "executed_command": self.executed_command,
            "contaminant_database_used": (
                str(self.contaminant_database_used) if self.contaminant_database_used else None
            ),
            "contaminant_database_sha256": self.contaminant_database_sha256,
            "runtime_seconds": round(self.runtime_seconds, 2),
        }


@dataclass
class SchmutziOutput:
    """Standardized output bundle wrapping Schmutzi execution results."""

    estimate: ContaminationEstimate
    diagnostics: ModelQualityDiagnostics
    outputs: ContaminationOutputFiles
    audit: SchmutziExecutionMetadata

    @property
    def cont_est_file(self) -> Optional[Path]:
        return self.outputs.cont_est_file

    @property
    def cont_deam_file(self) -> Optional[Path]:
        return self.outputs.cont_deam_file

    @property
    def endogenous_fasta(self) -> Optional[Path]:
        return self.outputs.endogenous_fasta

    @property
    def contaminant_fasta(self) -> Optional[Path]:
        return self.outputs.contaminant_fasta


class SchmutziEngine:
    """Wrapper and parser for schmutzi Bayesian mtDNA contamination estimation."""

    DOCKER_IMAGE: str = "schmutzi:latest"

    @classmethod
    def is_docker_available(cls) -> bool:
        """Checks if docker CLI is available and the schmutzi image exists."""
        if not shutil.which("docker"):
            return False
        try:
            res = subprocess.run(
                ["docker", "image", "inspect", cls.DOCKER_IMAGE],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=False,
            )
            return res.returncode == 0
        except Exception:
            return False

    @classmethod
    def get_binary_name(cls) -> Optional[str]:
        """Detects if schmutzi is available via Docker image or local binary in PATH."""
        if cls.is_docker_available():
            return f"docker:{cls.DOCKER_IMAGE}"
        for candidate in ["schmutzi.pl", "schmutzi"]:
            if shutil.which(candidate):
                return candidate
        return None

    @classmethod
    def get_schmutzi_version(cls) -> str:
        """Retrieves installed schmutzi version string."""
        if cls.is_docker_available():
            return "docker-container-latest"

        binary = shutil.which("schmutzi.pl") or shutil.which("schmutzi")
        if not binary:
            return "not_found"

        try:
            res = subprocess.run([binary, "--help"], capture_output=True, text=True, check=False)
            raw = f"{res.stdout}\n{res.stderr}".strip()
            match = re.search(r"schmutzi(?:\.pl)?\s+(?:version\s+)?([v\d\.]+)", raw, re.IGNORECASE)
            return match.group(1) if match else "installed_unknown_version"
        except Exception as exc:
            logger.warning(f"Could not determine schmutzi version: {exc}")
            return "unknown"

    @classmethod
    def compute_sha256(cls, target_path: Path) -> Optional[str]:
        """Computes SHA-256 digest of a target file or directory tree."""
        if not target_path.exists():
            return None

        hasher = hashlib.sha256()

        if target_path.is_file():
            with open(target_path, "rb") as f:
                while chunk := f.read(65536):
                    hasher.update(chunk)
            return hasher.hexdigest()

        if target_path.is_dir():
            for sub_file in sorted(target_path.rglob("*")):
                if sub_file.is_file():
                    hasher.update(sub_file.name.encode("utf-8"))
                    with open(sub_file, "rb") as f:
                        while chunk := f.read(65536):
                            hasher.update(chunk)
            return hasher.hexdigest()

        return None

    @classmethod
    def resolve_frequency_database_dir(cls, database_path: Path) -> Path:
        """Resolves the folder containing the .freq allele frequency files."""
        if not database_path.exists():
            raise SchmutziExecutionError(f"Contaminant database path does not exist: {database_path}")

        freqs_sub = database_path / "freqs"
        if freqs_sub.is_dir() and any(freqs_sub.glob("*.freq")):
            return freqs_sub

        if database_path.is_dir() and any(database_path.glob("*.freq")):
            return database_path

        return database_path

    @classmethod
    def compute_bam_mean_depth(
        cls,
        bam_path: Path,
        expected_contig: Optional[str] = None,
    ) -> float:
        """Computes mean sequencing depth across mitochondrial genome, capturing orphan and anomalous mates."""
        if not bam_path.is_file():
            raise SchmutziExecutionError(f"Cannot compute depth: BAM file does not exist: {bam_path}")

        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                target = expected_contig if (expected_contig and expected_contig in bam.references) else bam.references[0]
                target_len = bam.get_reference_length(target)
                if target_len == 0:
                    raise SchmutziExecutionError(f"Target reference contig '{target}' has length 0.")

                total_bases = sum(
                    col.nsegments
                    for col in bam.pileup(
                        contig=target,
                        start=0,
                        stop=target_len,
                        truncate=True,
                        stepper="all",  # Do not discard anomalous pairs or reads with mate on another contig (FLAG 129)
                    )
                )

                return float(total_bases / target_len)
        except Exception as exc:
            raise SchmutziExecutionError(f"Failed to compute depth metrics from {bam_path.name}: {exc}") from exc

    @classmethod
    def parse_cont_est_file(
        cls,
        cont_est_path: Path,
    ) -> Tuple[float, float, float, Optional[float]]:
        """Parses point estimate, 95% CI (lower, upper), and log-likelihood from .cont.est output."""
        if not cont_est_path.is_file() or cont_est_path.stat().st_size == 0:
            raise SchmutziExecutionError(f"Contamination estimates file missing or empty: {cont_est_path}")

        try:
            data_tokens: List[float] = []
            with open(cont_est_path, "r", encoding="utf-8") as f:
                for line in f:
                    line_str = line.strip()
                    if not line_str or line_str.startswith("#"):
                        continue
                    tokens = NUMBER_PATTERN.findall(line_str)
                    for token in tokens:
                        try:
                            data_tokens.append(float(token))
                        except ValueError:
                            continue

            if len(data_tokens) < 3:
                raise SchmutziExecutionError(
                    f"Insufficient numerical fields in {cont_est_path.name}. Found {len(data_tokens)}, expected at least 3."
                )

            point = data_tokens[0] / 100.0 if data_tokens[0] > 1.0 else data_tokens[0]
            low = data_tokens[1] / 100.0 if data_tokens[1] > 1.0 else data_tokens[1]
            high = data_tokens[2] / 100.0 if data_tokens[2] > 1.0 else data_tokens[2]
            log_lik = data_tokens[3] if len(data_tokens) >= 4 else None

            if not (0.0 <= point <= 1.0 and 0.0 <= low <= 1.0 and 0.0 <= high <= 1.0):
                raise SchmutziExecutionError(
                    f"Contamination bounds outside [0.0, 1.0] in {cont_est_path.name}: "
                    f"point={point}, low={low}, high={high}"
                )

            if low > high + 1e-5:
                raise SchmutziExecutionError(
                    f"Malformed CI ordering in {cont_est_path.name}: lower_ci ({low:.4f}) > upper_ci ({high:.4f})"
                )

            ci_span = max(0.0, high - low)
            allowed_tol = max(1e-4, 0.10 * ci_span)

            if point < (low - allowed_tol) or point > (high + allowed_tol):
                raise SchmutziExecutionError(
                    f"Gross statistical inconsistency in {cont_est_path.name}: "
                    f"point estimate ({point:.4f}) falls outside 95% CI [{low:.4f}, {high:.4f}] "
                    f"beyond tolerance threshold ({allowed_tol:.4f})."
                )

            if point < low:
                logger.debug(f"Rounding point estimate from {point:.5f} to lower CI bound {low:.5f}.")
                point = low
            elif point > high:
                logger.debug(f"Rounding point estimate from {point:.5f} to upper CI bound {high:.5f}.")
                point = high

            return point, low, high, log_lik

        except Exception as exc:
            if isinstance(exc, SchmutziExecutionError):
                raise
            raise SchmutziExecutionError(f"Failed parsing .cont.est file {cont_est_path.name}: {exc}") from exc

    @classmethod
    def parse_cont_deam_file(cls, cont_deam_path: Path) -> Optional[float]:
        """Parses preliminary contamination rate estimated from deamination alone (.contDeam.est)."""
        if not cont_deam_path.is_file() or cont_deam_path.stat().st_size == 0:
            return None

        try:
            with open(cont_deam_path, "r", encoding="utf-8") as f:
                for line in f:
                    line_str = line.strip()
                    if not line_str or line_str.startswith("#"):
                        continue
                    tokens = NUMBER_PATTERN.findall(line_str)
                    if tokens:
                        val = float(tokens[0])
                        return val / 100.0 if val > 1.0 else val
        except Exception as exc:
            logger.warning(f"Could not parse contDeam estimate from {cont_deam_path.name}: {exc}")
        return None

    @classmethod
    def parse_schmutzi_log(
        cls,
        log_path: Path,
    ) -> Tuple[bool, int, bool, bool]:
        """Parses execution log for strict MCMC convergence, iterations, and damage warnings."""
        converged = False
        iterations = 0
        warn_damage_pattern = False
        multiple_contaminants = False

        if not log_path.is_file():
            return converged, iterations, warn_damage_pattern, multiple_contaminants

        convergence_pattern = re.compile(
            r"stable\s+contamination\s+rate\s+at\s+iteration\s+#\s*(\d+)",
            re.IGNORECASE,
        )
        iter_pattern = re.compile(r"ITERATION\s+#\s*(\d+)", re.IGNORECASE)

        try:
            with open(log_path, "r", encoding="utf-8", errors="replace") as f:
                for line in f:
                    line_str = line.strip()
                    if not line_str:
                        continue

                    m_conv = convergence_pattern.search(line_str)
                    if m_conv:
                        converged = True
                        try:
                            iterations = max(iterations, int(m_conv.group(1)))
                        except ValueError:
                            pass

                    m_iter = iter_pattern.search(line_str)
                    if m_iter:
                        try:
                            iterations = max(iterations, int(m_iter.group(1)))
                        except ValueError:
                            pass

                    line_lower = line_str.lower()
                    if "deaminated contaminant" in line_lower or "damage on contaminant" in line_lower:
                        warn_damage_pattern = True
                    if "multiple contaminants" in line_lower or "mixture of contaminants" in line_lower:
                        multiple_contaminants = True

        except Exception as exc:
            logger.warning(f"Failed parsing schmutzi log {log_path.name}: {exc}")

        return converged, iterations, warn_damage_pattern, multiple_contaminants

    @classmethod
    def parse_contaminant_haplogroups(
        cls,
        final_mtcont_path: Optional[Path],
        final_cont_path: Optional[Path],
    ) -> List[str]:
        """Extracts the top predicted contaminant database identifiers from mtcont or final.cont."""
        contaminants: List[str] = []

        if final_cont_path and final_cont_path.is_file():
            try:
                with open(final_cont_path, "r", encoding="utf-8") as f:
                    for line in f:
                        line_str = line.strip()
                        if not line_str or line_str.startswith("#"):
                            continue
                        parts = line_str.split()
                        if parts:
                            cont_id = Path(parts[0]).stem
                            if cont_id not in contaminants:
                                contaminants.append(cont_id)
                                break
            except Exception as exc:
                logger.warning(f"Could not parse top contaminant from {final_cont_path.name}: {exc}")

        if not contaminants and final_mtcont_path and final_mtcont_path.is_file():
            try:
                with open(final_mtcont_path, "r", encoding="utf-8") as f:
                    for line in f:
                        line_str = line.strip()
                        if not line_str or line_str.startswith("#"):
                            continue
                        parts = line_str.split()
                        if parts:
                            cont_id = Path(parts[0]).stem
                            if cont_id not in contaminants:
                                contaminants.append(cont_id)
                                if len(contaminants) >= 3:
                                    break
            except Exception as exc:
                logger.warning(f"Could not parse contaminants from {final_mtcont_path.name}: {exc}")

        return contaminants

    @classmethod
    def count_informative_positions_from_logs(
        cls,
        endo_log: Optional[Path],
        cont_log: Optional[Path],
    ) -> int:
        """Counts positions where endogenous and contaminant consensus calls diverge."""
        if not endo_log or not cont_log or not endo_log.is_file() or not cont_log.is_file():
            return 0

        informative_count = 0
        try:
            endo_calls: Dict[int, str] = {}
            with open(endo_log, "r", encoding="utf-8") as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 3 and parts[0].isdigit():
                        endo_calls[int(parts[0])] = parts[2].upper()

            with open(cont_log, "r", encoding="utf-8") as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 3 and parts[0].isdigit():
                        pos = int(parts[0])
                        cont_base = parts[2].upper()
                        if pos in endo_calls:
                            endo_base = endo_calls[pos]
                            if (
                                endo_base != cont_base
                                and endo_base not in ("N", "-")
                                and cont_base not in ("N", "-")
                            ):
                                informative_count += 1
        except Exception as exc:
            logger.warning(f"Could not calculate informative positions from logs: {exc}")

        return informative_count

    @classmethod
    def _find_deliverable(cls, directory: Path, sample_id: str, patterns: List[str]) -> Optional[Path]:
        """Locates generated output files matching multiple common naming patterns."""
        for pattern in patterns:
            direct = directory / f"{sample_id}{pattern}"
            if direct.is_file():
                return direct
            matches = list(directory.glob(f"*{pattern}"))
            if matches:
                return matches[0]
        return None

    @classmethod
    def _wrap_command_for_execution(
        cls,
        cmd: List[str],
        repo_root: Path,
        sample_id: Optional[str] = None,
        stage: str = "run",
    ) -> List[str]:
        """Wraps CLI invocation in Docker container execution with an explicit, traceable container name."""
        if cls.is_docker_available():
            uid = os.getuid()
            gid = os.getgid()
            token = uuid.uuid4().hex[:6]
            safe_id = re.sub(r"[^\w.-]", "_", sample_id) if sample_id else "batch"
            container_name = f"schmutzi_{safe_id}_{stage}_{token}"

            docker_prefix = [
                "docker", "run", "--rm",
                "--name", container_name,
                "-u", f"{uid}:{gid}",
                "-v", f"{repo_root}:{repo_root}",
                "-w", str(repo_root),
                cls.DOCKER_IMAGE,
            ]
            return docker_prefix + cmd
        return cmd

    @classmethod
    def patch_posterior_deam(cls) -> bool:
        """Ensures posteriorDeam.R correctly interprets newline escape sequences in plot titles."""
        r_bin = shutil.which("posteriorDeam.R")
        if not r_bin:
            return False
        try:
            r_path = Path(r_bin)
            if not r_path.is_file():
                return False
            content = r_path.read_text(encoding="utf-8", errors="replace")
            if 'gsub("\\\\n", "\\n"' in content or 'gsub("\\n", "\n"' in content:
                return True
            if "title <- args[3]" in content:
                patched = content.replace(
                    "title <- args[3]",
                    'title <- gsub("\\\\n", "\\n", args[3])'
                )
                r_path.write_text(patched, encoding="utf-8")
                logger.info(f"Auto-patched {r_path} to resolve newline escapes in plot titles.")
                return True
        except Exception as exc:
            logger.debug(f"Could not auto-patch posteriorDeam.R: {exc}")
        return False

    @classmethod
    def run_schmutzi_pipeline(
        cls,
        sample_id: str,
        schmutzi_ready_bam: Path,
        reference_fasta: Path,
        contaminant_database: Path,
        output_dir: Path,
        log_dir: Path,
        policy: SchmutziPolicy,
        use_deam_prior: bool = True,
        wrapping_metadata: Optional[WrappingMetadata] = None,
        expected_contig: Optional[str] = None,
        timeout_seconds: int = 7200,
    ) -> Tuple[ContaminationEstimate, ModelQualityDiagnostics, ContaminationOutputFiles, SchmutziExecutionMetadata]:
        """Executes full schmutzi workflow with adaptive deamination prior and auto-recovery on -nan crash."""
        # Preemptive auto-patch to resolve literal \n escape sequences in PDF titles
        cls.patch_posterior_deam()

        resolved_freqs_dir = cls.resolve_frequency_database_dir(contaminant_database)

        out_sample_dir = output_dir / sample_id
        out_sample_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)

        prefix = out_sample_dir / sample_id
        stdout_log = log_dir / f"{sample_id}.schmutzi.log"
        stderr_log = log_dir / f"{sample_id}.schmutzi.err"

        repo_root = output_dir.resolve().parent.parent

        t0 = time.perf_counter()
        executed_commands: List[str] = []
        active_use_deam = use_deam_prior

        # Defensive resolution of library configuration to string
        lib_val = (
            policy.library_treatment.value
            if hasattr(policy.library_treatment, "value")
            else str(policy.library_treatment)
        ).upper()

        # PHASE 1: contDeam.pl (executed only if use_deam_prior is True)
        if active_use_deam:
            is_single = "SINGLE" in lib_val
            lib_str = "single" if is_single else "double"
            raw_cmd_deam = [
                "contDeam.pl",
                "--library", lib_str,
                "--out", str(prefix),
                "--ref", str(reference_fasta),
                str(schmutzi_ready_bam),
            ]

            exec_cmd_deam = cls._wrap_command_for_execution(
                raw_cmd_deam, repo_root, sample_id=sample_id, stage="contdeam"
            )

            logger.info(f"[{sample_id}] [Phase 1/2] Executing contDeam.pl (library={lib_str}, container={cls.is_docker_available()})...")
            logger.info(f"[{sample_id}] Invocation: {' '.join(exec_cmd_deam)}")
            executed_commands.append(" ".join(exec_cmd_deam))

            try:
                with open(stdout_log, "w", encoding="utf-8") as out_f, open(stderr_log, "w", encoding="utf-8") as err_f:
                    res_deam = subprocess.run(
                        exec_cmd_deam,
                        stdout=out_f,
                        stderr=err_f,
                        check=False,
                        timeout=timeout_seconds,
                    )
            except subprocess.TimeoutExpired as t_exc:
                raise SchmutziExecutionError(f"contDeam timed out after {timeout_seconds}s.") from t_exc
            except OSError as exc:
                raise SchmutziExecutionError(f"Failed to launch contDeam: {exc}") from exc

            if res_deam.returncode != 0:
                err_content = stderr_log.read_text(encoding="utf-8", errors="replace").strip()
                # Intercept -nan crash for zero-damage or modern data
                if "destringify() Unable to convert string=\"-nan\"" in err_content or "string=\"-nan\"" in err_content:
                    logger.warning(
                        f"[{sample_id}] contDeam crashed on zero-damage data (-nan). "
                        f"Auto-recovering: bypassing deamination prior and falling back to schmutzi --notusepred."
                    )
                    active_use_deam = False
                else:
                    raise SchmutziExecutionError(f"contDeam failed (exit code {res_deam.returncode}): {err_content}")
        else:
            logger.info(f"[{sample_id}] [Phase 1/2] Skipping contDeam (modern/zero-damage mode configured).")

        # PHASE 2: schmutzi.pl
        raw_cmd_schmutzi = [
            "schmutzi.pl",
            "--t", str(policy.threads),
            "--ref", str(reference_fasta),
        ]

        if not active_use_deam:
            raw_cmd_schmutzi.append("--notusepred")
            # Defensive injection of placeholder .deam.config to prevent crash at line 578 of schmutzi.pl
            deam_config_file = Path(f"{prefix}.deam.config")
            if not deam_config_file.is_file():
                try:
                    deam_config_file.touch()
                except Exception as exc:
                    logger.warning(f"[{sample_id}] Could not create placeholder .deam.config: {exc}")

        if any(tag in lib_val for tag in ("HALF", "PARTIAL", "USER_HALF", "UDG_HALF")):
            raw_cmd_schmutzi.append("--uselength")

        if policy.max_iterations and policy.max_iterations > 0:
            raw_cmd_schmutzi.extend(["--iterations", str(policy.max_iterations)])

        raw_cmd_schmutzi.extend([
            str(prefix),
            str(resolved_freqs_dir),
            str(schmutzi_ready_bam),
        ])

        exec_cmd_schmutzi = cls._wrap_command_for_execution(
            raw_cmd_schmutzi, repo_root, sample_id=sample_id, stage="schmutzi"
        )

        logger.info(
            f"[{sample_id}] [Phase 2/2] Executing schmutzi.pl "
            f"(threads={policy.threads}, use_deam_prior={active_use_deam})..."
        )
        logger.info(f"[{sample_id}] Invocation: {' '.join(exec_cmd_schmutzi)}")
        executed_commands.append(" ".join(exec_cmd_schmutzi))

        open_mode = "a" if active_use_deam else "w"
        try:
            with open(stdout_log, open_mode, encoding="utf-8") as out_f, open(stderr_log, open_mode, encoding="utf-8") as err_f:
                res_schmutzi = subprocess.run(
                    exec_cmd_schmutzi,
                    stdout=out_f,
                    stderr=err_f,
                    check=False,
                    timeout=timeout_seconds,
                )
        except subprocess.TimeoutExpired as t_exc:
            raise SchmutziExecutionError(f"schmutzi.pl timed out after {timeout_seconds}s.") from t_exc
        except OSError as exc:
            raise SchmutziExecutionError(f"Failed to launch schmutzi.pl: {exc}") from exc

        elapsed = time.perf_counter() - t0

        if res_schmutzi.returncode != 0:
            err_content = stderr_log.read_text(encoding="utf-8", errors="replace").strip()
            raise SchmutziExecutionError(f"schmutzi.pl terminated with exit code {res_schmutzi.returncode}: {err_content}")

        # Deliverable collection
        cont_est_file = cls._find_deliverable(
            out_sample_dir, sample_id, ["_final.cont.est", ".cont.est", "_cont.est"]
        )
        cont_deam_file = cls._find_deliverable(
            out_sample_dir, sample_id, [".contDeam.est", "_final.contDeam.est", "_contDeam.est"]
        ) if active_use_deam else None

        endo_fa = cls._find_deliverable(
            out_sample_dir, sample_id, ["_final_endo.fa", "_final.fasta", "_final_endo.fasta", "_endo.fa"]
        )
        cont_fa = cls._find_deliverable(
            out_sample_dir, sample_id, ["_final_cont.fa", "_final_cont.fasta", "_cont.fa"]
        )
        endo_log = cls._find_deliverable(
            out_sample_dir, sample_id, ["_final_endo.log", "_endo.log"]
        )
        cont_log = cls._find_deliverable(
            out_sample_dir, sample_id, ["_final_cont.log", "_cont.log"]
        )
        final_cont = cls._find_deliverable(
            out_sample_dir, sample_id, ["_final.cont", ".cont"]
        )
        final_mtcont = cls._find_deliverable(
            out_sample_dir, sample_id, ["_final_mtcont.out", ".mtcont.out"]
        )
        ready_bai = Path(f"{schmutzi_ready_bam}.bai")

        outputs = ContaminationOutputFiles(
            schmutzi_ready_bam=schmutzi_ready_bam if schmutzi_ready_bam.is_file() else None,
            schmutzi_ready_bai=ready_bai if ready_bai.is_file() else None,
            cont_est_file=cont_est_file,
            cont_deam_file=cont_deam_file,
            endogenous_fasta=endo_fa,
            contaminant_fasta=cont_fa,
        )

        if not cont_est_file or not cont_est_file.is_file():
            raise SchmutziExecutionError("schmutzi completed but contamination estimates file (.cont.est) was not found.")

        point_est, ci_low, ci_high, log_lik = cls.parse_cont_est_file(cont_est_file)
        prior_deam = cls.parse_cont_deam_file(cont_deam_file) if cont_deam_file and cont_deam_file.is_file() else None
        converged, iterations, warn_damage, multi_cont = cls.parse_schmutzi_log(stdout_log)
        contaminant_hgs = cls.parse_contaminant_haplogroups(final_mtcont, final_cont)
        info_sites = cls.count_informative_positions_from_logs(endo_log, cont_log)

        mean_depth = cls.compute_bam_mean_depth(
            bam_path=schmutzi_ready_bam,
            expected_contig=expected_contig,
        )

        estimate = ContaminationEstimate(
            point_estimate=point_est,
            ci_lower_95=ci_low,
            ci_upper_95=ci_high,
            method_used=ContaminationMethod.SCHMUTZI,
            cont_deam_prior_estimate=prior_deam,
            log_likelihood=log_lik,
        )

        ci_width = max(0.0, ci_high - ci_low)
        diagnostics = ModelQualityDiagnostics(
            mean_depth=mean_depth,
            informative_positions_evaluated=info_sites,
            ci_width=ci_width,
            model_converged=converged,
            iterations_run=iterations,
            warn_contaminant_damage_pattern=warn_damage,
            possible_multiple_contaminants=multi_cont,
            contaminant_haplogroups=contaminant_hgs,
            wrapping=wrapping_metadata or WrappingMetadata(),
        )

        audit = SchmutziExecutionMetadata(
            tool_version=cls.get_schmutzi_version(),
            executed_command="; ".join(executed_commands),
            contaminant_database_used=contaminant_database,
            contaminant_database_sha256=cls.compute_sha256(contaminant_database),
            runtime_seconds=elapsed,
        )

        logger.info(
            f"[{sample_id}] schmutzi completed in {elapsed:.2f}s | "
            f"Contamination: {point_est:.2%} (95% CI: [{ci_low:.2%}, {ci_high:.2%}], width={ci_width:.2%}, "
            f"mean depth={mean_depth:.1f}x, informative sites={info_sites}, converged={converged})"
        )

        return estimate, diagnostics, outputs, audit

    @classmethod
    def run(
        cls,
        sample_id: str,
        input_bam: Path,
        reference_fasta: Path,
        database_dir: Path,
        output_dir: Path,
        threads: int = 4,
        max_iterations: int = 10,
        use_wrapped: bool = False,
        use_deam_prior: bool = True,
        library_treatment: Any = LibraryTreatment.NONE,
        log_dir: Optional[Path] = None,
    ) -> SchmutziOutput:
        """High-level execution adapter matching ContaminationOrchestrator requirements."""
        policy = SchmutziPolicy(
            contaminant_database=database_dir,
            threads=threads,
            max_iterations=max_iterations,
            use_wrapped_reference=use_wrapped,
            library_treatment=library_treatment,
        )
        target_log_dir = log_dir or (output_dir / sample_id / "logs")

        estimate, diagnostics, outputs, audit = cls.run_schmutzi_pipeline(
            sample_id=sample_id,
            schmutzi_ready_bam=input_bam,
            reference_fasta=reference_fasta,
            contaminant_database=database_dir,
            output_dir=output_dir,
            log_dir=target_log_dir,
            policy=policy,
            use_deam_prior=use_deam_prior,
        )

        return SchmutziOutput(
            estimate=estimate,
            diagnostics=diagnostics,
            outputs=outputs,
            audit=audit,
        )