"""Adapter trimming, conservative read merging, physical output verification, and metrics extraction via AdapterRemoval."""

import gzip
import logging
from pathlib import Path
import re
import shutil
import subprocess
import time
from typing import Any, Dict, List, Optional, Tuple

from src.modules.preprocessing.models import (
    DataType,
    OutputProfile,
    PreprocessingMetrics,
    PreprocessingTool,
    ProcessedReads,
    SampleRecord,
    ToolMetadata,
    ValidationStatus,
)
from src.modules.preprocessing.integrity import IntegrityEngine, IntegrityError
from src.modules.preprocessing.fastq_validator import FastqValidator, FastqValidationError

logger = logging.getLogger("mtDNA_Pipeline.Trimmer")

PROJECT_ROOT = Path(__file__).resolve().parents[3]


class PreprocessingEngine:
    """Orchestrates conservative FASTQ adapter trimming and read merging using AdapterRemoval."""

    @classmethod
    def get_binary_name(cls) -> Optional[str]:
        """Detects if AdapterRemoval is registered as uppercase or lowercase in system PATH."""
        if shutil.which("AdapterRemoval"):
            return "AdapterRemoval"
        if shutil.which("adapterremoval"):
            return "adapterremoval"
        return None

    @classmethod
    def get_adapterremoval_version(cls) -> Tuple[str, str]:
        binary = cls.get_binary_name()
        if not binary:
            return "not_found", "AdapterRemoval binary not found in PATH"

        try:
            res = subprocess.run(
                [binary, "--version"],
                capture_output=True,
                text=True,
                check=False,
                timeout=30,
            )
            raw_output = f"{res.stdout}\n{res.stderr}".strip()
            match = re.search(r"AdapterRemoval\s+(?:ver(?:sion|\.)?)\s*([\d\.]+)", raw_output, re.IGNORECASE)
            version_str = match.group(1) if match else raw_output.split("\n")[0]
            return version_str, raw_output
        except Exception as exc:
            logger.warning(f"Could not determine AdapterRemoval version: {exc}")
            return "unknown", str(exc)

    @classmethod
    def _validate_preprocessing_metrics_invariants(
        cls,
        metrics: PreprocessingMetrics,
        is_paired: bool,
        warnings_list: List[str],
    ) -> None:
        """Enforces arithmetic consistency and structural invariants on parsed AdapterRemoval metrics."""
        # 1. Parity check for paired-end runs
        if is_paired and metrics.raw_total_reads is not None and metrics.raw_total_reads % 2 != 0:
            warnings_list.append(
                f"Structural anomaly: raw_total_reads ({metrics.raw_total_reads}) is odd in a Paired-End run."
            )

        # 2. Biological read conservation check (Raw ≈ Surviving + Discarded)
        if (
            metrics.raw_total_reads is not None
            and metrics.surviving_original_reads is not None
            and metrics.discarded_reads is not None
        ):
            accounted_reads = metrics.surviving_original_reads + metrics.discarded_reads
            if accounted_reads != metrics.raw_total_reads:
                warnings_list.append(
                    f"Metric conservation divergence: raw_total_reads ({metrics.raw_total_reads:,}) != "
                    f"surviving ({metrics.surviving_original_reads:,}) + discarded ({metrics.discarded_reads:,}) "
                    f"[Delta: {metrics.raw_total_reads - accounted_reads}]."
                )

    @classmethod
    def compute_post_trimming_qc(
        cls,
        fastq_paths: List[Optional[Path]],
        quality_base: int = 33,
        sample_reads: int = 25000,
    ) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[float]]:
        """Fast byte-level sampling of FASTQ files to compute GC%, Q20, Q30, and mean length."""
        total_bases = 0
        gc_bases = 0
        q20_bases = 0
        q30_bases = 0
        total_len = 0
        inspected_reads = 0
        q_thresh_20 = quality_base + 20
        q_thresh_30 = quality_base + 30

        for p in fastq_paths:
            if not p or not Path(p).is_file() or Path(p).stat().st_size == 0:
                continue
            open_fn = gzip.open if str(p).lower().endswith((".gz", ".bgz")) else open
            try:
                with open_fn(Path(p), "rb") as stream:
                    while inspected_reads < sample_reads:
                        header = stream.readline()
                        if not header:
                            break
                        seq = stream.readline().rstrip(b"\r\n")
                        plus = stream.readline()
                        qual = stream.readline().rstrip(b"\r\n")
                        if not qual:
                            break

                        inspected_reads += 1
                        s_len = len(seq)
                        total_bases += s_len
                        total_len += s_len

                        # Vectorized byte counts in C
                        gc_bases += seq.count(b"G") + seq.count(b"C") + seq.count(b"g") + seq.count(b"c")
                        q20_bases += sum(1 for b in qual if b >= q_thresh_20)
                        q30_bases += sum(1 for b in qual if b >= q_thresh_30)
            except Exception as exc:
                logger.debug(f"Sampling QC interrupted on {p}: {exc}")
                continue

        if total_bases == 0:
            return None, None, None, None

        gc_pct = round((gc_bases / total_bases) * 100.0, 2)
        q20_rate = round(q20_bases / total_bases, 4)
        q30_rate = round(q30_bases / total_bases, 4)
        mean_len = round(total_len / inspected_reads, 2) if inspected_reads > 0 else None

        return gc_pct, q20_rate, q30_rate, mean_len

    @classmethod
    def parse_adapterremoval_settings(
        cls,
        settings_file: Path,
        is_paired: bool,
        warnings_list: List[str],
    ) -> PreprocessingMetrics:
        """Parses AdapterRemoval .settings output file into structured metrics with strict None handling."""
        metrics = PreprocessingMetrics()
        if not settings_file.is_file():
            warnings_list.append(f"AdapterRemoval settings file not found: {settings_file.name}")
            return metrics

        try:
            with open(settings_file, "r", encoding="utf-8") as stream:
                content = stream.read()
        except OSError as exc:
            warnings_list.append(f"I/O error reading .settings file {settings_file.name}: {exc}")
            return metrics

        try:
            # 1. Raw read counts
            raw_pairs_match = re.search(
                r"Total number of (?:read pairs|read-pairs|pairs):\s+(\d+)",
                content,
                re.IGNORECASE,
            )
            raw_reads_match = re.search(
                r"Total number of reads:\s+(\d+)",
                content,
                re.IGNORECASE,
            )

            if raw_pairs_match:
                metrics.raw_read_pairs = int(raw_pairs_match.group(1))
                metrics.raw_total_reads = metrics.raw_read_pairs * 2
            elif raw_reads_match:
                metrics.raw_total_reads = int(raw_reads_match.group(1))
                if is_paired:
                    metrics.raw_read_pairs = metrics.raw_total_reads // 2
            else:
                warnings_list.append("Could not parse raw read counts from .settings file.")

            # 2. Merged / collapsed fragments
            collapsed_full = re.search(
                r"Number of full-length collapsed pairs:\s+(\d+)",
                content,
                re.IGNORECASE,
            )
            collapsed_trunc = re.search(
                r"Number of truncated collapsed pairs:\s+(\d+)",
                content,
                re.IGNORECASE,
            )
            collapsed_generic = re.search(
                r"Number of collapsed pairs:\s+(\d+)",
                content,
                re.IGNORECASE,
            )

            if collapsed_full or collapsed_trunc:
                c_full = int(collapsed_full.group(1)) if collapsed_full else 0
                c_trunc = int(collapsed_trunc.group(1)) if collapsed_trunc else 0
                metrics.merged_fragments = c_full + c_trunc
            elif collapsed_generic:
                metrics.merged_fragments = int(collapsed_generic.group(1))
            else:
                metrics.merged_fragments = None

            # 3. Surviving unmerged PE pairs
            retained_pairs_match = re.search(
                r"Number of (?:retained mate-pairs|retained pairs|retained mate pairs):\s+(\d+)",
                content,
                re.IGNORECASE,
            )
            retained_reads_match = re.search(
                r"Number of retained reads:\s+(\d+)",
                content,
                re.IGNORECASE,
            )

            if retained_pairs_match:
                metrics.surviving_read_pairs = int(retained_pairs_match.group(1))
            elif is_paired:
                metrics.surviving_read_pairs = None

            # 4. Singletons
            singleton_match = re.search(
                r"Number of (?:singletons|singleton reads):\s+(\d+)",
                content,
                re.IGNORECASE,
            )
            if singleton_match:
                metrics.singleton_reads = int(singleton_match.group(1))
            elif not is_paired and retained_reads_match:
                metrics.singleton_reads = int(retained_reads_match.group(1))
            else:
                metrics.singleton_reads = 0

            # 5. Discarded reads
            discarded_match = re.search(
                r"Number of discarded reads:\s+(\d+)",
                content,
                re.IGNORECASE,
            )
            if discarded_match:
                metrics.discarded_reads = int(discarded_match.group(1))
            else:
                metrics.discarded_reads = None

            # 6. Derived metrics computation
            surv_pairs = metrics.surviving_read_pairs or 0
            m_frags = metrics.merged_fragments or 0
            singles = metrics.singleton_reads or 0

            metrics.surviving_original_reads = (surv_pairs * 2) + (m_frags * 2) + singles
            metrics.surviving_sequenced_fragments = surv_pairs + m_frags + singles

            # 7. Merging fractions (Strict None preservation)
            if is_paired:
                if metrics.raw_read_pairs is not None and metrics.merged_fragments is not None:
                    if metrics.raw_read_pairs > 0:
                        metrics.merged_fraction_raw = metrics.merged_fragments / metrics.raw_read_pairs
                    else:
                        metrics.merged_fraction_raw = 0.0
                else:
                    metrics.merged_fraction_raw = None

                if metrics.merged_fragments is not None and metrics.surviving_read_pairs is not None:
                    total_surviving_pairs = metrics.merged_fragments + metrics.surviving_read_pairs
                    if total_surviving_pairs > 0:
                        metrics.merged_fraction_surviving = metrics.merged_fragments / total_surviving_pairs
                    else:
                        metrics.merged_fraction_surviving = 0.0
                else:
                    metrics.merged_fraction_surviving = None
            else:
                metrics.merged_fraction_raw = None
                metrics.merged_fraction_surviving = None

            # 8. Adapter trimming percentage with polymorphic regex and biological fallback
            m_adapt_pairs = re.search(
                r"Number of (?:well-aligned |retained )?(?:read pairs|mate-pairs|read-pairs|pairs) with adapter(?: sequence)?:\s+(\d+)",
                content,
                re.IGNORECASE,
            )
            m_adapt_reads = re.search(
                r"Number of reads with adapter(?: sequence)?:\s+(\d+)",
                content,
                re.IGNORECASE,
            )
            m_adapt_m1 = re.search(
                r"Number of mate 1 reads with adapter(?: sequence)?:\s+(\d+)",
                content,
                re.IGNORECASE,
            )
            m_adapt_m2 = re.search(
                r"Number of mate 2 reads with adapter(?: sequence)?:\s+(\d+)",
                content,
                re.IGNORECASE,
            )

            if m_adapt_pairs and metrics.raw_read_pairs and metrics.raw_read_pairs > 0:
                metrics.adapter_trimmed_reads_pct = round((int(m_adapt_pairs.group(1)) / metrics.raw_read_pairs) * 100.0, 2)
            elif m_adapt_reads and metrics.raw_total_reads and metrics.raw_total_reads > 0:
                metrics.adapter_trimmed_reads_pct = round((int(m_adapt_reads.group(1)) / metrics.raw_total_reads) * 100.0, 2)
            elif (m_adapt_m1 or m_adapt_m2) and metrics.raw_total_reads and metrics.raw_total_reads > 0:
                tot_adapt = (int(m_adapt_m1.group(1)) if m_adapt_m1 else 0) + (int(m_adapt_m2.group(1)) if m_adapt_m2 else 0)
                metrics.adapter_trimmed_reads_pct = round((tot_adapt / metrics.raw_total_reads) * 100.0, 2)

            # Biological fallback: in heavily degraded aDNA libraries where all reads are merged,
            # the fraction of adapter-trimmed sequences corresponds directly to the collapsed fraction
            if metrics.adapter_trimmed_reads_pct is None and metrics.merged_fraction_raw is not None:
                metrics.adapter_trimmed_reads_pct = round(metrics.merged_fraction_raw * 100.0, 2)

            # 9. Mean read lengths (pre and post trimming)
            m_pre1 = re.search(r"Average length of read 1 before trimming:\s+([\d\.]+)", content, re.IGNORECASE)
            m_pre2 = re.search(r"Average length of read 2 before trimming:\s+([\d\.]+)", content, re.IGNORECASE)
            if m_pre1 and m_pre2:
                metrics.mean_read_length_pre = round((float(m_pre1.group(1)) + float(m_pre2.group(1))) / 2.0, 2)
            elif m_pre1:
                metrics.mean_read_length_pre = round(float(m_pre1.group(1)), 2)

            m_post = re.search(r"Average length of retained reads:\s+([\d\.]+)", content, re.IGNORECASE)
            if m_post:
                metrics.mean_read_length_post = round(float(m_post.group(1)), 2)

        except (ValueError, re.error) as exc:
            warnings_list.append(f"Parsing syntax error in .settings file: {exc}")

        cls._validate_preprocessing_metrics_invariants(metrics, is_paired, warnings_list)
        return metrics

    @classmethod
    def _locate_and_verify_outputs(
        cls,
        base_prefix: Path,
        is_paired: bool,
        collapse_reads: bool,
        run_started_ts: float,
        integrity_mode: str,
    ) -> Tuple[ProcessedReads, List[str], List[str]]:
        """Scans and strictly verifies that newly generated outputs match run parameters and pass integrity gates."""
        processed = ProcessedReads()
        warnings: List[str] = []
        errors: List[str] = []

        candidate_merged = [
            Path(f"{base_prefix}.collapsed.gz"),
            Path(f"{base_prefix}.collapsed.truncated.gz"),
            Path(f"{base_prefix}.collapsed"),
        ]
        candidate_pair1 = [
            Path(f"{base_prefix}.pair1.truncated.gz"),
            Path(f"{base_prefix}.pair1.gz"),
            Path(f"{base_prefix}.pair1.truncated"),
            Path(f"{base_prefix}.pair1"),
        ]
        candidate_pair2 = [
            Path(f"{base_prefix}.pair2.truncated.gz"),
            Path(f"{base_prefix}.pair2.gz"),
            Path(f"{base_prefix}.pair2.truncated"),
            Path(f"{base_prefix}.pair2"),
        ]
        candidate_singletons = [
            Path(f"{base_prefix}.singleton.truncated.gz"),
            Path(f"{base_prefix}.singleton.gz"),
            Path(f"{base_prefix}.truncated.gz"),
            Path(f"{base_prefix}.singleton.truncated"),
            Path(f"{base_prefix}.truncated"),
        ]

        def _is_fresh_artifact(p: Path) -> bool:
            return p.is_file() and p.stat().st_size > 0 and (p.stat().st_mtime >= run_started_ts - 2.0)

        # Select valid files verified by timestamp
        if is_paired and collapse_reads:
            for m in candidate_merged:
                if _is_fresh_artifact(m):
                    processed.merged_fastq = m
                    break

        if is_paired:
            for p1 in candidate_pair1:
                if _is_fresh_artifact(p1):
                    processed.paired_r1_fastq = p1
                    break

            for p2 in candidate_pair2:
                if _is_fresh_artifact(p2):
                    processed.paired_r2_fastq = p2
                    break

        for s in candidate_singletons:
            if _is_fresh_artifact(s):
                processed.singleton_fastq = s
                break

        # Verify paired-end consistency
        if is_paired:
            if (processed.paired_r1_fastq and not processed.paired_r2_fastq) or (
                processed.paired_r2_fastq and not processed.paired_r1_fastq
            ):
                errors.append("Asymmetric paired-end output: found only one mate file for retained pairs.")

        active_outputs = [
            f for f in [
                processed.merged_fastq,
                processed.paired_r1_fastq,
                processed.paired_r2_fastq,
                processed.singleton_fastq,
            ]
            if f is not None
        ]

        if not active_outputs:
            errors.append("EMPTY_OUTPUT_AFTER_TRIMMING: AdapterRemoval produced no usable non-empty FASTQ files.")
            processed.output_profile = OutputProfile.EMPTY
            return processed, warnings, errors

        # Structural integrity check across all outputs
        for out_file in active_outputs:
            try:
                FastqValidator.validate_fastq_structure(out_file, max_records_to_inspect=10000)
                fp = IntegrityEngine.generate_fingerprint(out_file, mode=integrity_mode, is_input=False)
                processed.output_fingerprints.append(fp)
            except (FastqValidationError, IntegrityError) as exc:
                errors.append(f"Corrupted output file {out_file.name}: {exc}")

        # Exhaustive synchronization check on surviving PE pairs
        if processed.paired_r1_fastq and processed.paired_r2_fastq:
            try:
                FastqValidator.validate_paired_end_sync(
                    r1_path=processed.paired_r1_fastq,
                    r2_path=processed.paired_r2_fastq,
                    sample_check_records=5000,
                    verify_total_counts=True,
                )
            except FastqValidationError as sync_exc:
                errors.append(f"Generated paired-end output failed mate synchronization: {sync_exc}")

        # Determine output profile
        has_merged = processed.merged_fastq is not None
        has_paired = processed.paired_r1_fastq is not None and processed.paired_r2_fastq is not None
        has_single = processed.singleton_fastq is not None

        if has_merged and not has_paired and not has_single:
            processed.output_profile = OutputProfile.MERGED_ONLY
        elif has_paired and not has_merged and not has_single:
            processed.output_profile = OutputProfile.PAIRED_ONLY
        elif has_single and not has_merged and not has_paired:
            processed.output_profile = OutputProfile.SINGLETON_ONLY
        elif has_merged or has_paired or has_single:
            processed.output_profile = OutputProfile.MIXED
        else:
            processed.output_profile = OutputProfile.EMPTY

        return processed, warnings, errors

    @classmethod
    def run_adapterremoval(
        cls,
        sample: SampleRecord,
        output_dir: Path,
        log_dir: Path,
        min_read_length: int,
        min_merge_overlap: int,
        collapse_reads: bool,
        trim_ns: bool,
        threads: int,
        quality_base: Optional[int] = None,
        quality_max: Optional[int] = None,
        integrity_mode: str = "fast",
        timeout: int = 1800,
    ) -> ProcessedReads:
        """Executes AdapterRemoval with cross-mate encoding checks, pre-run cleanup, and post-run validation."""
        if sample.data_type != DataType.FASTQ:
            raise ValueError(f"PreprocessingEngine can only process FASTQ data, got: {sample.data_type}")

        binary = cls.get_binary_name()
        if not binary:
            raise RuntimeError("AdapterRemoval binary not found in system PATH.")

        r1_path = sample.inputs.r1_path
        r2_path = sample.inputs.r2_path
        is_paired = sample.inputs.is_paired

        if not r1_path or not r1_path.is_file():
            raise FileNotFoundError(f"[{sample.sample_id}] Forward read (R1) not found: {r1_path}")

        if is_paired and (not r2_path or not r2_path.is_file()):
            raise FileNotFoundError(f"[{sample.sample_id}] Reverse read (R2) not found: {r2_path}")

        output_dir = Path(output_dir).resolve()
        log_dir = Path(log_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)

        base_prefix = output_dir / sample.sample_id

        # 1. Automatic quality encoding detection with R1/R2 cross-check
        if quality_base is None or quality_max is None:
            auto_base_r1, auto_max_r1 = FastqValidator.detect_quality_encoding(r1_path)
            if is_paired and r2_path:
                auto_base_r2, auto_max_r2 = FastqValidator.detect_quality_encoding(r2_path)
                if auto_base_r1 != auto_base_r2:
                    raise RuntimeError(
                        f"[{sample.sample_id}] Conflicting quality score encodings between mate files: "
                        f"R1=Phred+{auto_base_r1}, R2=Phred+{auto_base_r2}."
                    )
            active_quality_base = quality_base if quality_base is not None else auto_base_r1
            active_quality_max = quality_max if quality_max is not None else auto_max_r1
        else:
            active_quality_base = quality_base
            active_quality_max = quality_max

        # 2. Preventive cleanup of residual files (prevent stale outputs)
        for stale_file in output_dir.glob(f"{sample.sample_id}.*"):
            if stale_file.is_file():
                try:
                    stale_file.unlink()
                except OSError as exc:
                    logger.warning(f"Could not remove stale output file {stale_file.name}: {exc}")

        settings_file = Path(f"{base_prefix}.settings")
        stdout_log = log_dir / f"{sample.sample_id}.adapterremoval.log"
        stderr_log = log_dir / f"{sample.sample_id}.adapterremoval.err"

        # 3. Command construction
        cmd = [
            binary,
            "--basename", str(base_prefix),
            "--minlength", str(min_read_length),
            "--qualitybase", str(active_quality_base),
            "--qualitymax", str(active_quality_max),
            "--threads", str(threads),
            "--gzip",
        ]
        if is_paired:
            cmd.extend([
                "--file1", str(r1_path),
                "--file2", str(r2_path),
            ])
            if collapse_reads:
                cmd.extend([
                    "--collapse",
                    "--minadapteroverlap", str(min_merge_overlap),
                ])
        else:
            cmd.extend([
                "--file1", str(r1_path),
            ])

        if trim_ns:
            cmd.append("--trimns")

        # 4. Audit trail
        tool_version, raw_version = cls.get_adapterremoval_version()
        sample.tool_meta = ToolMetadata(
            tool_name=PreprocessingTool.ADAPTER_REMOVAL,
            tool_version=tool_version,
            tool_raw_version=raw_version,
            command_line=" ".join(cmd),
            parameters={
                "min_read_length": min_read_length,
                "min_merge_overlap": min_merge_overlap if (is_paired and collapse_reads) else None,
                "collapse_reads": collapse_reads if is_paired else False,
                "trim_ns": trim_ns,
                "quality_base": active_quality_base,
                "quality_max": active_quality_max,
                "threads": threads,
            },
        )

        logger.info(
            f"[{sample.sample_id}] Executing AdapterRemoval (v{tool_version}) | "
            f"min_length={min_read_length}, collapse={collapse_reads}, trim_ns={trim_ns}, "
            f"qualitybase={active_quality_base}, qualitymax={active_quality_max}"
        )

        sample.mark_started()
        run_started_ts = time.time()
        try:
            with open(stdout_log, "w", encoding="utf-8") as out_f, open(stderr_log, "w", encoding="utf-8") as err_f:
                res = subprocess.run(cmd, stdout=out_f, stderr=err_f, check=False, timeout=timeout)
        except subprocess.TimeoutExpired as t_exc:
            sample.validation.add_error(f"AdapterRemoval timed out after {timeout} seconds.")
            sample.validation.status = ValidationStatus.FAIL
            raise RuntimeError(f"[{sample.sample_id}] AdapterRemoval execution timed out ({timeout}s)") from t_exc
        except OSError as exc:
            sample.validation.add_error(f"Failed to launch AdapterRemoval: {exc}")
            sample.validation.status = ValidationStatus.FAIL
            raise RuntimeError(f"Failed to execute AdapterRemoval for {sample.sample_id}: {exc}") from exc
        finally:
            sample.mark_finished()

        # 5. Exit code verification and stderr capture
        if res.returncode != 0:
            err_content = ""
            if stderr_log.is_file():
                try:
                    err_content = stderr_log.read_text(encoding="utf-8", errors="replace").strip()
                except Exception:
                    pass
            err_msg = (
                f"AdapterRemoval failed with exit code {res.returncode}. "
                f"Log details:\n{err_content or 'No stderr details available'}"
            )
            sample.validation.add_error(err_msg)
            sample.validation.status = ValidationStatus.FAIL
            raise RuntimeError(f"[{sample.sample_id}] {err_msg}")

        # 6. Parse metrics from .settings
        sample.metrics = cls.parse_adapterremoval_settings(
            settings_file=settings_file,
            is_paired=is_paired,
            warnings_list=sample.validation.warnings,
        )

        # 7. Validate and map outputs
        processed, out_warnings, out_errors = cls._locate_and_verify_outputs(
            base_prefix=base_prefix,
            is_paired=is_paired,
            collapse_reads=collapse_reads,
            run_started_ts=run_started_ts,
            integrity_mode=integrity_mode,
        )

        for w in out_warnings:
            sample.validation.add_warning(w)

        if out_errors:
            for e in out_errors:
                sample.validation.add_error(e)
            sample.validation.status = ValidationStatus.FAIL
            raise RuntimeError(f"[{sample.sample_id}] Output validation failed: {'; '.join(out_errors)}")

        # 8. Compute post-trimming quality metrics (GC%, Q20, Q30, and post-read length)
        active_outputs = [
            f for f in [
                processed.merged_fastq,
                processed.paired_r1_fastq,
                processed.paired_r2_fastq,
                processed.singleton_fastq,
            ]
            if f is not None and f.is_file()
        ]
        gc_pct, q20_r, q30_r, sampled_post_len = cls.compute_post_trimming_qc(
            fastq_paths=active_outputs,
            quality_base=active_quality_base,
            sample_reads=25000,
        )
        sample.metrics.gc_content_pct = gc_pct
        sample.metrics.q20_rate_post = q20_r
        sample.metrics.q30_rate_post = q30_r
        if sample.metrics.mean_read_length_post is None:
            sample.metrics.mean_read_length_post = sampled_post_len

        # Fallback for mean_read_length_pre if not parsed from .settings
        if sample.metrics.mean_read_length_pre is None and r1_path and r1_path.is_file():
            _, _, _, sampled_pre_len = cls.compute_post_trimming_qc(
                fastq_paths=[r1_path],
                quality_base=active_quality_base,
                sample_reads=2000,
            )
            sample.metrics.mean_read_length_pre = sampled_pre_len

        # Final guarantee for adapter_trimmed_reads_pct on merged reads
        if sample.metrics.adapter_trimmed_reads_pct is None and sample.metrics.merged_fraction_raw is not None:
            sample.metrics.adapter_trimmed_reads_pct = round(sample.metrics.merged_fraction_raw * 100.0, 2)

        sample.processed_reads = processed
        if sample.validation.status != ValidationStatus.FAIL:
            sample.validation.status = ValidationStatus.PASS

        m_raw_pct = f"{sample.metrics.merged_fraction_raw:.2%}" if sample.metrics.merged_fraction_raw is not None else "N/A"
        m_surv_pct = (
            f"{sample.metrics.merged_fraction_surviving:.2%}"
            if sample.metrics.merged_fraction_surviving is not None
            else "N/A"
        )
        adapt_pct_str = (
            f"{sample.metrics.adapter_trimmed_reads_pct:.2f}%"
            if sample.metrics.adapter_trimmed_reads_pct is not None
            else "N/A"
        )

        logger.info(
            f"[{sample.sample_id}] Preprocessing Completed [Profile: {processed.output_profile.value}] | "
            f"Original Reads Retained: {sample.metrics.surviving_original_reads} | "
            f"Alignable Fragments: {sample.metrics.surviving_sequenced_fragments} "
            f"(Merged: {sample.metrics.merged_fragments}, Retained Pairs: {sample.metrics.surviving_read_pairs}, "
            f"Singletons: {sample.metrics.singleton_reads}) | "
            f"Merged Fraction (Raw): {m_raw_pct} | Merged Fraction (Surviving): {m_surv_pct} | "
            f"Adapter Trimmed: {adapt_pct_str} | GC%: {sample.metrics.gc_content_pct}% | "
            f"Q20: {sample.metrics.q20_rate_post} | Q30: {sample.metrics.q30_rate_post}"
        )
        return processed

    @classmethod
    def process_sample(cls, sample: SampleRecord, config: Any) -> SampleRecord:
        """Dispatches sample preprocessing extracting configuration parameters."""
        if sample.data_type != DataType.FASTQ:
            logger.info(f"[{sample.sample_id}] Ingestion Mode: BAM. Skipping FASTQ preprocessing.")
            return sample

        def _resolve_dir(p_str: Optional[str], default_rel: Path) -> Path:
            if not p_str:
                return default_rel if default_rel.is_absolute() else (PROJECT_ROOT / default_rel)
            p = Path(p_str)
            return p if p.is_absolute() else (PROJECT_ROOT / p)

        if isinstance(config, dict):
            prep_cfg = config.get("preprocessing", {})
            out_cfg = config.get("output", {})
            min_length = int(prep_cfg.get("min_read_length", 30))
            min_overlap = int(prep_cfg.get("min_merge_overlap", 11))
            collapse = bool(prep_cfg.get("collapse_reads", True))
            trim_ns = bool(prep_cfg.get("trim_ns", True))
            threads = int(prep_cfg.get("threads", 4))
            quality_base = prep_cfg.get("quality_base", None)
            quality_max = prep_cfg.get("quality_max", None)
            integrity_mode = str(config.get("integrity", {}).get("mode", "fast"))
            results_dir = _resolve_dir(out_cfg.get("results_dir"), PROJECT_ROOT / "results")
        else:
            min_length = int(getattr(config, "min_read_length", 30))
            min_overlap = int(getattr(config, "min_merge_overlap", 11))
            collapse = bool(getattr(config, "collapse_reads", True))
            trim_ns = bool(getattr(config, "trim_ns", True))
            threads = int(getattr(config, "threads", 4))
            quality_base = getattr(config, "quality_base", None)
            quality_max = getattr(config, "quality_max", None)
            integrity_mode = str(getattr(config, "integrity_mode", "fast"))
            results_dir = _resolve_dir(getattr(config, "results_dir", None), PROJECT_ROOT / "results")

        out_prep_dir = results_dir / "alignment" / sample.sample_id / "trimmed"
        log_dir = results_dir / "logs" / "preprocessing"

        cls.run_adapterremoval(
            sample=sample,
            output_dir=out_prep_dir,
            log_dir=log_dir,
            min_read_length=min_length,
            min_merge_overlap=min_overlap,
            collapse_reads=collapse,
            trim_ns=trim_ns,
            threads=threads,
            quality_base=quality_base,
            quality_max=quality_max,
            integrity_mode=integrity_mode,
        )
        return sample