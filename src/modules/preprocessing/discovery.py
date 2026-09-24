"""Robust sample discovery, token-based filename parsing, deterministic pairing, and multi-lane aggregation engine."""

from dataclasses import dataclass, field
import logging
from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Set, Tuple

from src.modules.preprocessing.models import (
    DataType,
    FileFingerprint,
    InputFiles,
    SampleIdentity,
    SampleRecord,
    ValidationResult,
    ValidationStatus,
)
from src.modules.preprocessing.integrity import IntegrityEngine, IntegrityError
from src.modules.preprocessing.fastq_validator import FastqValidator, FastqValidationError

logger = logging.getLogger("mtDNA_Pipeline.Discovery")


@dataclass(frozen=True)
class PairingKey:
    """Hashable unique identifier for paired-end mate matching."""

    sample_id: str
    library_id: str
    lane: str


@dataclass
class DiscoveredFile:
    """Intermediate metadata container derived from tokenized filename inspection."""

    filepath: Path
    filename: str
    sample_id: str
    library_id: str
    lane: str
    read_type: str  # "R1", "R2", "SE", "BAM", "UNKNOWN"
    file_format: str  # "FASTQ", "BAM", "UNKNOWN"
    compression: str  # "gzip", "bgzip", "none"

    @property
    def pairing_key(self) -> PairingKey:
        """Generates the pairing key for mate matching."""
        return PairingKey(
            sample_id=self.sample_id,
            library_id=self.library_id,
            lane=self.lane,
        )


@dataclass
class SequencingUnit:
    """Represents an atomic sequencing unit (e.g., a specific lane or library) for a sample."""

    unit_id: str
    lane: str
    library_id: str
    is_paired: bool
    r1_path: Optional[Path] = None
    r2_path: Optional[Path] = None
    bam_path: Optional[Path] = None
    fingerprints: List[FileFingerprint] = field(default_factory=list)


class FilenameTokenizer:
    """Deconstructs filenames from right-to-left into semantic sequencing metadata tokens."""

    # Default fallback Illumina read markers
    R1_PATTERNS = re.compile(r"^(r1|1|read1)$", re.IGNORECASE)
    R2_PATTERNS = re.compile(r"^(r2|2|read2)$", re.IGNORECASE)
    LANE_PATTERNS = re.compile(r"^(l\d+|lane\d+)$", re.IGNORECASE)
    SAMPLE_NUM_PATTERNS = re.compile(r"^s\d+$", re.IGNORECASE)
    SET_NUM_PATTERNS = re.compile(r"^\d{3}$")  # e.g., '001'

    @classmethod
    def compile_marker_pattern(
        cls, user_patterns: Optional[List[str]], default_pattern: re.Pattern
    ) -> re.Pattern:
        """Compiles user-defined prefix/suffix strings into an isolated token-matching regex."""
        if not user_patterns:
            return default_pattern
        tokens = []
        for p in user_patterns:
            clean = str(p).strip().lstrip("._-").lower()
            if clean:
                tokens.append(re.escape(clean))
        if not tokens:
            return default_pattern
        return re.compile(f"^({'|'.join(tokens)})$", re.IGNORECASE)

    @classmethod
    def parse(
        cls,
        file_path: Path,
        custom_r1_pattern: Optional[re.Pattern] = None,
        custom_r2_pattern: Optional[re.Pattern] = None,
    ) -> DiscoveredFile:
        """Parses a file path into a structured DiscoveredFile descriptor."""
        path = Path(file_path).resolve()
        name = path.name

        # 1. Compression and extension stripping
        compression = "none"
        working_name = name

        if working_name.lower().endswith(".gz"):
            compression = "gzip"
            working_name = working_name[:-3]
        elif working_name.lower().endswith(".bgz"):
            compression = "bgzip"
            working_name = working_name[:-4]

        file_format = "UNKNOWN"
        if working_name.lower().endswith((".fastq", ".fq")):
            file_format = "FASTQ"
            working_name = re.sub(r"\.(fastq|fq)$", "", working_name, flags=re.IGNORECASE)
        elif working_name.lower().endswith(".bam"):
            file_format = "BAM"
            compression = "bgzip"
            working_name = re.sub(r"\.bam$", "", working_name, flags=re.IGNORECASE)
            return DiscoveredFile(
                filepath=path,
                filename=name,
                sample_id=working_name,
                library_id="lib01",
                lane="lane01",
                read_type="BAM",
                file_format="BAM",
                compression=compression,
            )

        # 2. Right-to-left semantic tokenization
        tokens = re.split(r"[._-]", working_name)

        lane = "L001"
        library_id = "lib01"
        read_type = "SE"
        consumed_idx = len(tokens)

        r1_pat = custom_r1_pattern or cls.R1_PATTERNS
        r2_pat = custom_r2_pattern or cls.R2_PATTERNS

        i = len(tokens) - 1
        while i >= 0:
            token = tokens[i]

            # Discard trailing set numbers (e.g., '001')
            if i == len(tokens) - 1 and cls.SET_NUM_PATTERNS.match(token):
                consumed_idx = i
                i -= 1
                continue

            # Identify R1 or R2 mates
            if r1_pat.match(token):
                read_type = "R1"
                consumed_idx = i
                i -= 1
                continue
            elif r2_pat.match(token):
                read_type = "R2"
                consumed_idx = i
                i -= 1
                continue

            # Identify Lane index (e.g., L001, L1, lane1)
            if cls.LANE_PATTERNS.match(token):
                lane = token.upper()
                consumed_idx = i
                i -= 1
                continue

            # Identify Sample / Library number index (e.g., S1, S92)
            if cls.SAMPLE_NUM_PATTERNS.match(token):
                library_id = token.upper()
                consumed_idx = i
                i -= 1
                continue

            break

        sample_id = "_".join(tokens[:consumed_idx]) if consumed_idx > 0 else working_name
        if not sample_id:
            sample_id = working_name

        return DiscoveredFile(
            filepath=path,
            filename=name,
            sample_id=sample_id,
            library_id=library_id,
            lane=lane,
            read_type=read_type,
            file_format=file_format,
            compression=compression,
        )


class SampleDiscoveryEngine:
    """Orchestrates sample discovery, strict pairing validation, and pre-flight reporting."""

    @staticmethod
    def _extract_config_section(config: Any, section_name: str) -> Dict[str, Any]:
        """Robustly extracts configuration dictionaries supporting dict, dataclass, or wrapper types."""
        if isinstance(config, dict):
            if section_name in config:
                val = config[section_name]
                return val if isinstance(val, dict) else {}
            raw_data = config.get("raw_data")
            if isinstance(raw_data, dict) and section_name in raw_data:
                val = raw_data[section_name]
                return val if isinstance(val, dict) else {}
            return {}

        raw_cfg = getattr(config, "raw_data", None)
        if isinstance(raw_cfg, dict) and section_name in raw_cfg:
            val = raw_cfg[section_name]
            return val if isinstance(val, dict) else {}

        if hasattr(config, section_name):
            val = getattr(config, section_name)
            if isinstance(val, dict):
                return val
            if hasattr(val, "__dict__"):
                return val.__dict__

        return {}

    @classmethod
    def discover_files(
        cls,
        input_dir: Path,
        custom_r1_patterns: Optional[List[str]] = None,
        custom_r2_patterns: Optional[List[str]] = None,
    ) -> List[DiscoveredFile]:
        """Scans input directory for all supported genomic files while filtering out OS/system artifacts."""
        dir_path = Path(input_dir).resolve()
        if not dir_path.is_dir():
            raise FileNotFoundError(f"Input directory does not exist: {dir_path}")

        valid_extensions = {".fastq.gz", ".fq.gz", ".fastq", ".fq", ".bam", ".gz", ".bgz"}

        def _is_valid_candidate(p: Path) -> bool:
            if not p.is_file():
                return False
            name = p.name.lower()
            # Discard hidden files, temporary files, and Windows Zone.Identifier metadata
            if name.startswith((".", "~")) or "zone.identifier" in name or name.endswith((".tmp", ".aria2")):
                return False
            return any(name.endswith(ext) for ext in valid_extensions)

        raw_files = [f for f in dir_path.iterdir() if _is_valid_candidate(f)]

        if not raw_files:
            raise FileNotFoundError(f"No valid sequencing files (.fastq, .fq, .bam, .gz) found in: {dir_path}")

        r1_regex = FilenameTokenizer.compile_marker_pattern(custom_r1_patterns, FilenameTokenizer.R1_PATTERNS)
        r2_regex = FilenameTokenizer.compile_marker_pattern(custom_r2_patterns, FilenameTokenizer.R2_PATTERNS)

        discovered = [
            FilenameTokenizer.parse(f, custom_r1_pattern=r1_regex, custom_r2_pattern=r2_regex)
            for f in sorted(raw_files)
        ]
        logger.info(f"Discovered {len(discovered)} raw file(s) in {dir_path.name}")
        return discovered

    @classmethod
    def pair_and_validate_fastq(
        cls,
        discovered_files: List[DiscoveredFile],
        layout: str = "AUTO",
        allow_missing_mates: bool = False,
        single_end_mode: Optional[bool] = None,
    ) -> Tuple[Dict[str, List[SequencingUnit]], List[str], List[str]]:
        """Executes mate pairing, single-end ingestion, or mixed-mode resolution across samples."""
        if single_end_mode is True:
            mode = "SINGLE"
        else:
            mode = layout.upper() if layout else "AUTO"

        if mode not in {"AUTO", "MIXED", "PAIRED", "SINGLE"}:
            mode = "AUTO"

        r1_map: Dict[PairingKey, DiscoveredFile] = {}
        r2_map: Dict[PairingKey, DiscoveredFile] = {}
        se_list: List[DiscoveredFile] = []

        warnings: List[str] = []
        errors: List[str] = []

        for f in discovered_files:
            if f.file_format != "FASTQ":
                continue

            if f.read_type == "R1":
                if f.pairing_key in r1_map:
                    errors.append(f"Duplicate R1 file detected for key {f.pairing_key}: {f.filename} vs {r1_map[f.pairing_key].filename}")
                r1_map[f.pairing_key] = f
            elif f.read_type == "R2":
                if f.pairing_key in r2_map:
                    errors.append(f"Duplicate R2 file detected for key {f.pairing_key}: {f.filename} vs {r2_map[f.pairing_key].filename}")
                r2_map[f.pairing_key] = f
            else:
                se_list.append(f)

        sample_units: Dict[str, List[SequencingUnit]] = {}

        # 1. Explicit Single-End layout override
        if mode == "SINGLE":
            all_se_files = list(r1_map.values()) + list(r2_map.values()) + se_list
            for f in all_se_files:
                unit = SequencingUnit(
                    unit_id=f"{f.library_id}_{f.lane}",
                    lane=f.lane,
                    library_id=f.library_id,
                    is_paired=False,
                    r1_path=f.filepath,
                )
                sample_units.setdefault(f.sample_id, []).append(unit)
            return sample_units, warnings, errors

        # 2. Pairing resolution for PAIRED, AUTO, and MIXED layouts
        all_keys: Set[PairingKey] = set(r1_map.keys()).union(set(r2_map.keys()))

        for key in sorted(all_keys, key=lambda k: (k.sample_id, k.lane, k.library_id)):
            r1_file = r1_map.get(key)
            r2_file = r2_map.get(key)

            # Matched paired-end unit
            if r1_file and r2_file:
                unit = SequencingUnit(
                    unit_id=f"{key.library_id}_{key.lane}",
                    lane=key.lane,
                    library_id=key.library_id,
                    is_paired=True,
                    r1_path=r1_file.filepath,
                    r2_path=r2_file.filepath,
                )
                sample_units.setdefault(key.sample_id, []).append(unit)

            # R1 present without corresponding R2
            elif r1_file and not r2_file:
                msg = f"Missing mate R2 for file: {r1_file.filename} (Sample: {key.sample_id}, Lane: {key.lane})"
                if mode in {"AUTO", "MIXED"}:
                    if allow_missing_mates:
                        warnings.append(f"{msg} -> Fallback to Single-End.")
                        unit = SequencingUnit(
                            unit_id=f"{key.library_id}_{key.lane}",
                            lane=key.lane,
                            library_id=key.library_id,
                            is_paired=False,
                            r1_path=r1_file.filepath,
                        )
                        sample_units.setdefault(key.sample_id, []).append(unit)
                    else:
                        errors.append(f"FAIL: {msg}. Set allow_missing_mates=True if this sample is Single-End despite R1 marker.")
                else:  # PAIRED mode
                    if allow_missing_mates:
                        warnings.append(f"{msg} -> Fallback to Single-End.")
                        unit = SequencingUnit(
                            unit_id=f"{key.library_id}_{key.lane}",
                            lane=key.lane,
                            library_id=key.library_id,
                            is_paired=False,
                            r1_path=r1_file.filepath,
                        )
                        sample_units.setdefault(key.sample_id, []).append(unit)
                    else:
                        errors.append(f"FAIL: {msg}")

            # Orphan R2 mate
            elif r2_file and not r1_file:
                errors.append(f"Orphan Reverse read R2 without matching R1: {r2_file.filename} (Sample: {key.sample_id})")

        # 3. Handle files without explicit R1/R2 markers
        for f in se_list:
            if mode in {"AUTO", "MIXED"}:
                unit = SequencingUnit(
                    unit_id=f"{f.library_id}_{f.lane}",
                    lane=f.lane,
                    library_id=f.library_id,
                    is_paired=False,
                    r1_path=f.filepath,
                )
                sample_units.setdefault(f.sample_id, []).append(unit)
            elif mode == "PAIRED":
                if allow_missing_mates:
                    warnings.append(f"File '{f.filename}' has no mate indicators; treated as Single-End (allow_missing_mates=True).")
                    unit = SequencingUnit(
                        unit_id=f"{f.library_id}_{f.lane}",
                        lane=f.lane,
                        library_id=f.library_id,
                        is_paired=False,
                        r1_path=f.filepath,
                    )
                    sample_units.setdefault(f.sample_id, []).append(unit)
                else:
                    errors.append(
                        f"File '{f.filename}' does not contain R1/R2 mate markers while layout=PAIRED. "
                        f"Set 'layout: AUTO' in config to process mixed Single-End and Paired-End cohorts."
                    )

        return sample_units, warnings, errors

    @classmethod
    def print_preflight_report(cls, sample_records: List[SampleRecord]) -> None:
        """Prints a structured ASCII Pre-Flight Discovery Summary table to stdout and logger."""
        header = f"{'SAMPLE ID':<18} | {'TYPE':<6} | {'LAYOUT':<6} | {'LANES':<6} | {'STATUS':<8} | {'NOTES'}"
        divider = "-" * 85

        lines = [
            "\n" + divider,
            " PRE-FLIGHT SAMPLE DISCOVERY REPORT",
            divider,
            header,
            divider,
        ]

        for s in sample_records:
            layout = "PE" if s.inputs.is_paired else "SE"
            lanes_count = len(s.fingerprints) // 2 if s.inputs.is_paired else len(s.fingerprints)
            if lanes_count == 0:
                lanes_count = 1

            notes = ""
            if s.validation.errors:
                notes = f"ERR: {s.validation.errors[0]}"
            elif s.validation.warnings:
                notes = f"WARN: {s.validation.warnings[0]}"
            else:
                notes = "Ready for preprocessing"

            lines.append(
                f"{s.sample_id:<18} | {s.data_type.value:<6} | {layout:<6} | {lanes_count:<6} | "
                f"{s.status.value:<8} | {notes}"
            )

        lines.append(divider + "\n")
        summary_text = "\n".join(lines)
        print(summary_text)
        logger.info("Pre-Flight Discovery Report generated.")

    @classmethod
    def build_sample_records(
        cls,
        input_dir: Path,
        config: Any,
        integrity_mode: str = "fast",
    ) -> List[SampleRecord]:
        """Main orchestrator for sample discovery, validation, and SampleRecord construction."""
        input_sec = cls._extract_config_section(config, "input")
        val_sec = cls._extract_config_section(config, "validation")

        allow_missing = bool(input_sec.get("allow_missing_mates", False))

        raw_layout = str(input_sec.get("layout", "")).strip().upper()
        if not raw_layout:
            single_end_bool = input_sec.get("single_end")
            if single_end_bool is True:
                layout = "SINGLE"
            elif single_end_bool is False:
                layout = "PAIRED"
            else:
                layout = "AUTO"
        else:
            layout = raw_layout

        pairing_patterns = input_sec.get("pairing_patterns", {})
        custom_r1 = pairing_patterns.get("r1") if isinstance(pairing_patterns, dict) else None
        custom_r2 = pairing_patterns.get("r2") if isinstance(pairing_patterns, dict) else None

        strict_mode = bool(val_sec.get("strict_mode", False))

        discovered = cls.discover_files(
            input_dir=input_dir,
            custom_r1_patterns=custom_r1,
            custom_r2_patterns=custom_r2,
        )

        sample_records: List[SampleRecord] = []

        # --- Branch 1: BAM Files ---
        bam_files = [f for f in discovered if f.file_format == "BAM"]
        if bam_files:
            for b in bam_files:
                sample_rec = SampleRecord(
                    identity=SampleIdentity(sample_id=b.sample_id),
                    inputs=InputFiles(data_type=DataType.BAM, is_paired=False, bam_path=b.filepath),
                )
                try:
                    fp = IntegrityEngine.generate_fingerprint(b.filepath, is_input=True)
                    sample_rec.fingerprints.append(fp)
                except IntegrityError as exc:
                    sample_rec.validation.add_error(str(exc))

                sample_records.append(sample_rec)

            cls.print_preflight_report(sample_records)
            return sample_records

        # --- Branch 2: FASTQ Files ---
        units_map, warnings, errors = cls.pair_and_validate_fastq(
            discovered_files=discovered,
            layout=layout,
            allow_missing_mates=allow_missing,
        )

        for sample_id, units in units_map.items():
            first_unit = units[0]
            val_res = ValidationResult()

            sample_rec = SampleRecord(
                identity=SampleIdentity(sample_id=sample_id, library_id=first_unit.library_id),
                inputs=InputFiles(
                    data_type=DataType.FASTQ,
                    is_paired=first_unit.is_paired,
                    r1_path=first_unit.r1_path,
                    r2_path=first_unit.r2_path,
                ),
                validation=val_res,
            )

            for u in units:
                try:
                    if u.r1_path:
                        fp1 = IntegrityEngine.generate_fingerprint(u.r1_path, is_input=True)
                        sample_rec.fingerprints.append(fp1)
                    if u.r2_path:
                        fp2 = IntegrityEngine.generate_fingerprint(u.r2_path, is_input=True)
                        sample_rec.fingerprints.append(fp2)

                    if integrity_mode == "strict" and u.is_paired and u.r1_path and u.r2_path:
                        FastqValidator.validate_paired_end_sync(u.r1_path, u.r2_path, sample_check_records=5000)

                except (IntegrityError, FastqValidationError) as exc:
                    sample_rec.validation.add_error(str(exc))

            if sample_rec.validation.status == ValidationStatus.UNTESTED:
                sample_rec.validation.status = ValidationStatus.PASS

            sample_records.append(sample_rec)

        if errors:
            for err in errors:
                logger.error(f"[Discovery Error] {err}")
            if strict_mode:
                cls.print_preflight_report(sample_records)
                raise RuntimeError(f"Discovery aborted due to {len(errors)} critical pairing/input errors.")

        for warn in warnings:
            logger.warning(f"[Discovery Warning] {warn}")

        cls.print_preflight_report(sample_records)
        return sample_records