"""Track A Engine: mapDamage2 execution, terminal misincorporation profiling, and Bayesian damage parameter estimation."""

import csv
from dataclasses import dataclass
import logging
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pysam

from src.modules.damage.models import (
    BayesianDamageParameters,
    DamageErrorCode,
    DamageOutputFiles,
    DamagePolicy,
    DamageWarningCode,
    FragmentLengthMetrics,
    LibraryTreatment,
    TerminalDamageFrequencies,
)

logger = logging.getLogger("mtDNA_Pipeline.DamageQuantifier")


class DamageQuantifierError(Exception):
    """Raised when mapDamage execution, pre-flight validation, or output parsing encounters an unrecoverable failure."""
    pass


class DamageQuantifierEngine:
    """Orchestrates mapDamage2 quantification and parses empirical mismatch matrices and Bayesian stats."""

    @classmethod
    def get_binary_name(cls) -> Optional[str]:
        """Detects if mapDamage is registered in system PATH."""
        if shutil.which("mapDamage"):
            return "mapDamage"
        if shutil.which("mapdamage"):
            return "mapdamage"
        return None

    @classmethod
    def get_mapdamage_version(cls) -> Tuple[str, str]:
        """Retrieves installed mapDamage version string."""
        binary = cls.get_binary_name()
        if not binary:
            return "not_found", "mapDamage binary not found in system PATH."

        try:
            res = subprocess.run([binary, "--version"], capture_output=True, text=True, check=False)
            raw = f"{res.stdout}\n{res.stderr}".strip()
            match = re.search(r"mapDamage\s+(?:version\s+)?([\d\.]+)", raw, re.IGNORECASE)
            version_str = match.group(1) if match else raw.split("\n")[0]
            return version_str, raw
        except Exception as exc:
            logger.warning(f"Could not determine mapDamage version: {exc}")
            return "unknown", str(exc)

    @classmethod
    def validate_bam_reference_compatibility(cls, bam_path: Path, ref_fasta: Path) -> None:
        """Validates that contigs and lengths declared in BAM header strictly match the reference FASTA."""
        if not bam_path.is_file():
            raise DamageQuantifierError(f"Input BAM file not found on disk: {bam_path}")
        if not ref_fasta.is_file():
            raise DamageQuantifierError(f"Reference FASTA file not found on disk: {ref_fasta}")

        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam, pysam.FastaFile(str(ref_fasta)) as fa:
                bam_sq = dict(zip(bam.references, bam.lengths))
                fa_sq = dict(zip(fa.references, fa.lengths))

                # Verify that all contigs in BAM exist in the reference FASTA with identical lengths
                for contig, bam_len in bam_sq.items():
                    if contig not in fa_sq:
                        raise DamageQuantifierError(
                            f"BAM contig '{contig}' does not exist in reference FASTA '{ref_fasta.name}'."
                        )
                    if fa_sq[contig] != bam_len:
                        raise DamageQuantifierError(
                            f"Length mismatch for contig '{contig}': BAM={bam_len} bp, FASTA={fa_sq[contig]} bp."
                        )
        except Exception as exc:
            if isinstance(exc, DamageQuantifierError):
                raise
            raise DamageQuantifierError(f"Pre-flight BAM/Reference compatibility check failed: {exc}") from exc

    @classmethod
    def parse_misincorporation_file(
        cls,
        misinc_path: Path,
        window_length: int = 25,
    ) -> Tuple[TerminalDamageFrequencies, List[str]]:
        """Parses mapDamage2 misincorporation.txt using explicit C>T and G>A misincorporation columns.

        Args:
            misinc_path: Path to misincorporation.txt file.
            window_length: Maximum terminal positions to evaluate.

        Returns:
            Tuple of populated TerminalDamageFrequencies and parsing warnings.
        """
        warnings: List[str] = []
        if not misinc_path.is_file() or misinc_path.stat().st_size == 0:
            raise DamageQuantifierError(f"misincorporation.txt missing or empty: {misinc_path}")

        # Initialize accumulation arrays (positions 1..window_length)
        c_ref_5p = [0] * (window_length + 1)
        ct_mut_5p = [0] * (window_length + 1)
        g_ref_3p = [0] * (window_length + 1)
        ga_mut_3p = [0] * (window_length + 1)

        try:
            with open(misinc_path, "r", encoding="utf-8") as f:
                reader = csv.reader(f, delimiter="\t")
                header = None

                for row in reader:
                    if not row or not row[0].strip():
                        continue

                    # Header identification and field normalization
                    if row[0].startswith("#") or "End" in row or "Pos" in row or "Chr" in row:
                        header = [col.lstrip("#").strip().upper() for col in row]
                        continue

                    if not header:
                        continue

                    row_dict = dict(zip(header, [r.strip() for r in row]))

                    end_val = row_dict.get("END", "").lower()
                    try:
                        pos_val = int(row_dict.get("POS", -1))
                    except ValueError:
                        continue

                    if pos_val < 1 or pos_val > window_length:
                        continue

                    # Extract C->T and G->A counts from mapDamage mismatch columns
                    # Standard columns: C>T, G>A, C, G, Total, or individual base counts A, C, G, T
                    try:
                        ct_count = int(row_dict.get("C>T", 0))
                        ga_count = int(row_dict.get("G>A", 0))
                        
                        # Denominators: check whether line specifies reference base context
                        base_ref = row_dict.get("BASE", "").upper()
                        total_col = int(row_dict.get("TOTAL", 0)) if "TOTAL" in row_dict else 0
                        
                        cnt_c = int(row_dict.get("C", 0)) if "C" in row_dict else 0
                        cnt_g = int(row_dict.get("G", 0)) if "G" in row_dict else 0
                        cnt_a = int(row_dict.get("A", 0)) if "A" in row_dict else 0
                        cnt_t = int(row_dict.get("T", 0)) if "T" in row_dict else 0
                    except (ValueError, KeyError):
                        continue

                    # 5' End Accumulation (C -> T)
                    if "5p" in end_val:
                        ct_mut_5p[pos_val] += ct_count
                        # Compute reference C denominator
                        if base_ref == "C" and total_col > 0:
                            c_ref_5p[pos_val] += total_col
                        elif base_ref == "C":
                            c_ref_5p[pos_val] += (cnt_a + cnt_c + cnt_g + cnt_t)
                        else:
                            # If lines are not split by base context, total C sites = C_unmutated + C>T
                            c_ref_5p[pos_val] += cnt_c + ct_count

                    # 3' End Accumulation (G -> A)
                    elif "3p" in end_val:
                        ga_mut_3p[pos_val] += ga_count
                        # Compute reference G denominator
                        if base_ref == "G" and total_col > 0:
                            g_ref_3p[pos_val] += total_col
                        elif base_ref == "G":
                            g_ref_3p[pos_val] += (cnt_a + cnt_c + cnt_g + cnt_t)
                        else:
                            g_ref_3p[pos_val] += cnt_g + ga_count

        except Exception as exc:
            raise DamageQuantifierError(f"Error parsing misincorporation file {misinc_path.name}: {exc}") from exc

        # Construct profiles, counts, and observed frequency vectors
        c_to_t_profile: List[float] = []
        g_to_a_profile: List[float] = []
        c_to_t_obs: List[int] = []
        c_to_t_den: List[int] = []
        g_to_a_obs: List[int] = []
        g_to_a_den: List[int] = []

        for p in range(1, window_length + 1):
            den_5p = c_ref_5p[p]
            obs_5p = ct_mut_5p[p]
            frac_5p = (obs_5p / den_5p) if den_5p > 0 else 0.0
            c_to_t_profile.append(frac_5p)
            c_to_t_obs.append(obs_5p)
            c_to_t_den.append(den_5p)

            den_3p = g_ref_3p[p]
            obs_3p = ga_mut_3p[p]
            frac_3p = (obs_3p / den_3p) if den_3p > 0 else 0.0
            g_to_a_profile.append(frac_3p)
            g_to_a_obs.append(obs_3p)
            g_to_a_den.append(den_3p)

        c1 = c_to_t_profile[0] if c_to_t_profile else 0.0
        c2 = c_to_t_profile[1] if len(c_to_t_profile) > 1 else 0.0
        g1 = g_to_a_profile[0] if g_to_a_profile else 0.0
        g2 = g_to_a_profile[1] if len(g_to_a_profile) > 1 else 0.0

        freqs = TerminalDamageFrequencies(
            c_to_t_5p_pos1=c1,
            c_to_t_5p_pos2=c2,
            c_to_t_5p_profile=c_to_t_profile,
            max_5p_c_to_t=max(c_to_t_profile) if c_to_t_profile else 0.0,
            g_to_a_3p_pos1=g1,
            g_to_a_3p_pos2=g2,
            g_to_a_3p_profile=g_to_a_profile,
            max_3p_g_to_a=max(g_to_a_profile) if g_to_a_profile else 0.0,
            n_c_to_t_5p_pos1_obs=c_to_t_obs[0] if c_to_t_obs else 0,
            n_c_ref_5p_pos1_total=c_to_t_den[0] if c_to_t_den else 0,
            n_g_to_a_3p_pos1_obs=g_to_a_obs[0] if g_to_a_obs else 0,
            n_g_ref_3p_pos1_total=g_to_a_den[0] if g_to_a_den else 0,
            c_to_t_5p_observations=c_to_t_obs,
            c_to_t_5p_denominators=c_to_t_den,
            g_to_a_3p_observations=g_to_a_obs,
            g_to_a_3p_denominators=g_to_a_den,
        )

        return freqs, warnings

    @classmethod
    def parse_length_distribution_file(
        cls,
        lg_dist_path: Path,
    ) -> FragmentLengthMetrics:
        """Parses lgdistribution.txt supporting both 2-column and 3-column (Std/Length/Occurences) formats."""
        if not lg_dist_path.is_file() or lg_dist_path.stat().st_size == 0:
            return FragmentLengthMetrics()

        lengths: List[int] = []
        counts: List[int] = []

        try:
            with open(lg_dist_path, "r", encoding="utf-8") as f:
                header = None
                for line in f:
                    line_str = line.strip()
                    if not line_str or line_str.startswith("#"):
                        continue

                    parts = re.split(r"\s+", line_str)
                    if not parts:
                        continue

                    # Header detection
                    if "Length" in parts or "Occurences" in parts or "Std" in parts:
                        header = [p.strip().lower() for p in parts]
                        continue

                    try:
                        if header and "length" in header and "occurences" in header:
                            ln = int(parts[header.index("length")])
                            cnt = int(parts[header.index("occurences")])
                        elif len(parts) >= 3:
                            # Standard format mapDamage 2.2+: [Std, Length, Occurences]
                            ln = int(parts[1])
                            cnt = int(parts[2])
                        elif len(parts) == 2:
                            # Compact format: [Length, Count]
                            ln = int(parts[0])
                            cnt = int(parts[1])
                        else:
                            continue

                        if cnt > 0:
                            lengths.append(ln)
                            counts.append(cnt)
                    except (ValueError, IndexError):
                        continue

            if not lengths:
                return FragmentLengthMetrics()

            total_fragments = sum(counts)
            if total_fragments == 0:
                return FragmentLengthMetrics()

            lengths_arr = np.array(lengths)
            counts_arr = np.array(counts)

            mean_len = float(np.sum(lengths_arr * counts_arr) / total_fragments)

            # Calculate weighted quantiles
            sort_idx = np.argsort(lengths_arr)
            sorted_lengths = lengths_arr[sort_idx]
            sorted_counts = counts_arr[sort_idx]
            cum_counts = np.cumsum(sorted_counts)
            cum_fractions = cum_counts / total_fragments

            p25 = float(sorted_lengths[np.searchsorted(cum_fractions, 0.25)])
            median = float(sorted_lengths[np.searchsorted(cum_fractions, 0.50)])
            p75 = float(sorted_lengths[np.searchsorted(cum_fractions, 0.75)])
            p90 = float(sorted_lengths[np.searchsorted(cum_fractions, 0.90)])

            return FragmentLengthMetrics(
                total_fragments_sampled=total_fragments,
                mean_length=mean_len,
                median_length=median,
                p25_length=p25,
                p75_length=p75,
                p90_length=p90,
                min_length=int(np.min(sorted_lengths)),
                max_length=int(np.max(sorted_lengths)),
            )

        except Exception as exc:
            logger.warning(f"Failed to calculate fragment length metrics from {lg_dist_path.name}: {exc}")
            return FragmentLengthMetrics()

    @classmethod
    def parse_mcmc_summary_file(
        cls,
        summ_stat_path: Path,
    ) -> BayesianDamageParameters:
        """Parses MCMC summary statistics from Stats_out_MCMC_iter_summ_stat.csv with convergence validation.

        Args:
            summ_stat_path: Path to Stats_out_MCMC_iter_summ_stat.csv file.

        Returns:
            Populated BayesianDamageParameters dataclass instance.
        """
        params = BayesianDamageParameters()
        if not summ_stat_path.is_file() or summ_stat_path.stat().st_size == 0:
            return params

        try:
            with open(summ_stat_path, "r", encoding="utf-8") as f:
                reader = csv.reader(f)
                header = None
                
                for row in reader:
                    if not row:
                        continue
                    if header is None:
                        header = [c.strip().lower() for c in row]
                        continue

                    # mapDamage row structure: Parameter / Stat / Mean / Median / 2.5% / 97.5% / Std
                    row_dict = dict(zip(header, [c.strip() for c in row]))
                    param_name = row_dict.get("param", row_dict.get("parameter", row[0])).strip().lower()

                    def get_stat_val(d: Dict[str, str], row_list: List[str]) -> Optional[float]:
                        for key in ["mean", "median", "estimate", "value"]:
                            if key in d:
                                try:
                                    return float(d[key])
                                except ValueError:
                                    pass
                        # Fallback to second column if header fields are unlabelled
                        if len(row_list) > 1:
                            try:
                                return float(row_list[1])
                            except ValueError:
                                pass
                        return None

                    val = get_stat_val(row_dict, row)
                    if val is None:
                        continue

                    if "lambda" in param_name:
                        params.lambda_param = val
                    elif "delta_s" in param_name or "deltas" in param_name:
                        params.delta_s = val
                    elif "delta_d" in param_name or "deltad" in param_name:
                        params.delta_d = val
                    elif "theta" in param_name:
                        params.theta = val
                    elif "loglik" in param_name or "log_likelihood" in param_name or "log-likelihood" in param_name:
                        params.model_log_likelihood = val

            # Verify presence of all cardinal parameters within theoretical probability bounds [0.0, 1.0]
            if (
                params.lambda_param is not None
                and params.delta_s is not None
                and params.delta_d is not None
                and params.theta is not None
                and not np.isnan(params.lambda_param)
                and not np.isnan(params.delta_s)
                and 0.0 <= params.delta_s <= 1.0
                and 0.0 <= params.delta_d <= 1.0
            ):
                params.model_converged = True

        except Exception as exc:
            logger.warning(f"Could not parse MCMC summary file {summ_stat_path.name}: {exc}")

        return params

    @classmethod
    def run_mapdamage_stats(
        cls,
        sample_id: str,
        input_bam: Path,
        reference_fasta: Path,
        output_dir: Path,
        log_dir: Path,
        policy: DamagePolicy,
        threads: int = 2,
    ) -> Tuple[DamageOutputFiles, TerminalDamageFrequencies, FragmentLengthMetrics, BayesianDamageParameters]:
        """Executes mapDamage2 in statistical estimation mode (without rescaling), collecting tables and MCMC models.

        Args:
            sample_id: Biological sample identifier.
            input_bam: Path to input BAM (dedup or markdup).
            reference_fasta: Path to target reference FASTA (rCRS).
            output_dir: Destination folder for mapDamage outputs.
            log_dir: Destination folder for process stdout/stderr logs.
            policy: Validated DamagePolicy container.
            threads: Number of worker threads.

        Returns:
            Tuple of (DamageOutputFiles, TerminalDamageFrequencies, FragmentLengthMetrics, BayesianDamageParameters).

        Raises:
            DamageQuantifierError: If pre-flight validation fails or mapDamage terminates with errors.
        """
        binary = cls.get_binary_name()
        if not binary:
            raise DamageQuantifierError("mapDamage binary not found in system PATH. Ensure it is installed in Conda.")

        # 1. Pre-flight Validation
        cls.validate_bam_reference_compatibility(input_bam, reference_fasta)

        out_sample_dir = output_dir / sample_id
        out_sample_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)

        stdout_log = log_dir / f"{sample_id}.mapdamage.log"
        stderr_log = log_dir / f"{sample_id}.mapdamage.err"

        # 2. Command invocation (statistical quantification and Bayesian fit only, without --rescale)
        cmd = [
            binary,
            "-i", str(input_bam),
            "-r", str(reference_fasta),
            "-d", str(out_sample_dir),
            "-l", str(policy.terminal_positions_to_profile),
        ]

        if policy.library_treatment == LibraryTreatment.SINGLE_STRANDED:
            cmd.append("--single-stranded")

        logger.info(
            f"[{sample_id}] Executing mapDamage2 Statistical Quantification "
            f"(window={policy.terminal_positions_to_profile} bp, library={policy.library_treatment.value})..."
        )

        # 3. Process Execution
        try:
            with open(stdout_log, "w", encoding="utf-8") as out_f, open(stderr_log, "w", encoding="utf-8") as err_f:
                res = subprocess.run(cmd, stdout=out_f, stderr=err_f, check=False)
        except OSError as exc:
            raise DamageQuantifierError(f"Failed to launch mapDamage process for {sample_id}: {exc}") from exc

        if res.returncode != 0:
            raise DamageQuantifierError(
                f"mapDamage failed with exit code {res.returncode}. Inspect error log: {stderr_log}"
            )

        # 4. Output File Mapping
        misinc_file = out_sample_dir / "misincorporation.txt"
        lg_dist_file = out_sample_dir / "lgdistribution.txt"
        plot_pdf = out_sample_dir / "Fragmisincorporation_plot.pdf"
        summ_stat_file = out_sample_dir / "Stats_out_MCMC_iter_summ_stat.csv"
        correct_prob_file = out_sample_dir / "Stats_out_MCMC_correct_prob.csv"

        if not misinc_file.is_file():
            raise DamageQuantifierError(
                f"mapDamage completed with exit code 0 but '{misinc_file.name}' was not produced."
            )

        outputs = DamageOutputFiles(
            misincorporation_txt=misinc_file if misinc_file.is_file() else None,
            length_distribution_txt=lg_dist_file if lg_dist_file.is_file() else None,
            damage_plot_pdf=plot_pdf if plot_pdf.is_file() else None,
        )

        # 5. Extract profiles and empirical parameters
        terminal_freqs, _ = cls.parse_misincorporation_file(
            misinc_path=misinc_file,
            window_length=policy.terminal_positions_to_profile,
        )
        length_metrics = cls.parse_length_distribution_file(
            lg_dist_path=lg_dist_file,
        )
        bayesian_params = cls.parse_mcmc_summary_file(
            summ_stat_path=summ_stat_file,
        )

        logger.info(
            f"[{sample_id}] Quantification Completed | "
            f"5' C->T (pos1): {terminal_freqs.c_to_t_5p_pos1:.2%} ({terminal_freqs.n_c_to_t_5p_pos1_obs}/{terminal_freqs.n_c_ref_5p_pos1_total} obs) | "
            f"3' G->A (pos1): {terminal_freqs.g_to_a_3p_pos1:.2%} ({terminal_freqs.n_g_to_a_3p_pos1_obs}/{terminal_freqs.n_g_ref_3p_pos1_total} obs) | "
            f"Median Length: {length_metrics.median_length:.1f} bp | "
            f"Bayesian Fit Converged: {bayesian_params.model_converged}"
        )

        return outputs, terminal_freqs, length_metrics, bayesian_params