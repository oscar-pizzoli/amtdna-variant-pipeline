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

    # Common Illumina read markers
    R1_PATTERNS = re.compile(r"^(r1|1|read1)$", re.IGNORECASE)
    R2_PATTERNS = re.compile(r"^(r2|2|read2)$", re.IGNORECASE)
    LANE_PATTERNS = re.compile(r"^(l\d+|lane\d+)$", re.IGNORECASE)
    SAMPLE_NUM_PATTERNS = re.compile(r"^s\d+$", re.IGNORECASE)
    SET_NUM_PATTERNS = re.compile(r"^\d{3}$")  # e.g., '001'

    @classmethod
    def parse(cls, file_path: Path) -> DiscoveredFile:
        """Parses a file path into a structured DiscoveredFile descriptor."""
        path = Path(file_path).resolve()
        name = path.name

        # 1. Riconoscimento compressione ed estensione
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

        # 2. Tokenizzazione da destra verso sinistra
        tokens = re.split(r"[._-]", working_name)

        lane = "L001"
        library_id = "lib01"
        read_type = "SE"
        consumed_idx = len(tokens)

        # Scansione tokens dal fondo
        i = len(tokens) - 1
        while i >= 0:
            token = tokens[i]

            # Scarta eventuale '001' terminale Illumina
            if i == len(tokens) - 1 and cls.SET_NUM_PATTERNS.match(token):
                consumed_idx = i
                i -= 1
                continue

            # Riconoscimento R1 / R2 (solo se separati da delimitatori canonici)
            if cls.R1_PATTERNS.match(token):
                read_type = "R1"
                consumed_idx = i
                i -= 1
                continue
            elif cls.R2_PATTERNS.match(token):
                read_type = "R2"
                consumed_idx = i
                i -= 1
                continue

            # Riconoscimento Lane (L001, L1, lane1)
            if cls.LANE_PATTERNS.match(token):
                lane = token.upper()
                consumed_idx = i
                i -= 1
                continue

            # Riconoscimento Sample/Library Index (S1, S92)
            if cls.SAMPLE_NUM_PATTERNS.match(token):
                library_id = token.upper()
                consumed_idx = i
                i -= 1
                continue

            break

        # Tutto ciò che precede i metadati tecnici è l'ID biologico pulito
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
    def discover_files(cls, input_dir: Path) -> List[DiscoveredFile]:
        """Scans input directory for all supported genomic files."""
        dir_path = Path(input_dir).resolve()
        if not dir_path.is_dir():
            raise FileNotFoundError(f"Input directory does not exist: {dir_path}")

        valid_extensions = {".gz", ".bgz", ".fastq", ".fq", ".bam"}
        raw_files = [
            f for f in dir_path.iterdir()
            if f.is_file() and not f.name.startswith((".", "~")) and any(f.name.lower().endswith(ext) for ext in valid_extensions)
        ]

        if not raw_files:
            raise FileNotFoundError(f"No valid sequencing files (.fastq, .fq, .bam, .gz) found in: {dir_path}")

        discovered = [FilenameTokenizer.parse(f) for f in sorted(raw_files)]
        logger.info(f"Discovered {len(discovered)} raw file(s) in {dir_path.name}")
        return discovered

    @classmethod
    def pair_and_validate_fastq(
        cls,
        discovered_files: List[DiscoveredFile],
        allow_missing_mates: bool = False,
        single_end_mode: bool = False,
    ) -> Tuple[Dict[str, List[SequencingUnit]], List[str], List[str]]:
        """Executes mate pairing or single-end ingestion without requiring artificial file renaming."""
        r1_map: Dict[PairingKey, DiscoveredFile] = {}
        r2_map: Dict[PairingKey, DiscoveredFile] = {}
        se_list: List[DiscoveredFile] = []

        warnings: List[str] = []
        errors: List[str] = []

        # 1. Smistamento iniziale dei file FASTQ
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

        # 2. Ingestione Single-End Universale (nessun suffisso obbligatorio)
        if single_end_mode:
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

        # 3. Pairing Deterministico Paired-End
        all_keys: Set[PairingKey] = set(r1_map.keys()).union(set(r2_map.keys()))

        for key in sorted(all_keys, key=lambda k: (k.sample_id, k.lane, k.library_id)):
            r1_file = r1_map.get(key)
            r2_file = r2_map.get(key)

            # Caso Nominale: Coppia R1 + R2 presente
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

            # Caso Anomalo: R1 senza R2
            elif r1_file and not r2_file:
                msg = f"Missing mate R2 for file: {r1_file.filename} (Sample: {key.sample_id}, Lane: {key.lane})"
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

            # Caso Anomalo: R2 Orfano
            elif r2_file and not r1_file:
                errors.append(f"Orphan Reverse read R2 without matching R1: {r2_file.filename} (Sample: {key.sample_id})")

        # 4. Gestione file privi di R1/R2 in modalità Paired-End
        if se_list:
            if allow_missing_mates:
                for f in se_list:
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
                for f in se_list:
                    errors.append(
                        f"File '{f.filename}' does not contain R1/R2 mate markers while single_end=False. "
                        f"Set 'single_end: true' in config to process Single-End files."
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
            notes = ""
            if s.validation.errors:
                notes = f"ERR: {s.validation.errors[0]}"
            elif s.validation.warnings:
                notes = f"WARN: {s.validation.warnings[0]}"
            else:
                notes = "Ready for preprocessing"

            lines.append(
                f"{s.sample_id:<18} | {s.data_type.value:<6} | {layout:<6} | {len(s.fingerprints):<6} | "
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
        discovered = cls.discover_files(input_dir)

        # Risoluzione robusta della configurazione
        input_sec = cls._extract_config_section(config, "input")
        val_sec = cls._extract_config_section(config, "validation")

        allow_missing = bool(input_sec.get("allow_missing_mates", False))
        single_end = bool(input_sec.get("single_end", False))
        strict_mode = bool(val_sec.get("strict_mode", False))

        sample_records: List[SampleRecord] = []

        # --- RAMO 1: File BAM ---
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

        # --- RAMO 2: File FASTQ ---
        units_map, warnings, errors = cls.pair_and_validate_fastq(
            discovered_files=discovered,
            allow_missing_mates=allow_missing,
            single_end_mode=single_end,
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