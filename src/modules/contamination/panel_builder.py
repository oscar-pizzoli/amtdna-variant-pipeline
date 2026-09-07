"""Diagnostic Panel Builder for mtDNA Fast Screening.

Extracts polymorphic and informative diagnostic sites from schmutzi's 197-genome
allele frequency database and rCRS reference FASTA to generate versioned,
allele-aware, self-describing panel TSVs.
"""

import argparse
from dataclasses import dataclass
import hashlib
import json
import logging
import math
from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Set, Tuple

import pysam
import yaml

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger("mtDNA_Pipeline.PanelBuilder")

# Risoluzione dinamica della radice del progetto
PROJECT_ROOT = Path(__file__).resolve().parents[3]


class PanelBuilderValidationError(Exception):
    """Raised when the diagnostic panel fails post-construction invariant validation."""
    pass


@dataclass
class DiagnosticSite:
    """Represents a validated polymorphic diagnostic site."""

    position: int
    ref_base: str
    primary_nonref_allele: str
    contaminant_maf: float
    database_allele_observations: float
    informativeness_score: float
    allele_breakdown: Dict[str, float]


class DiagnosticPanelBuilder:
    """Parses schmutzi frequency tables and builds validated allele-aware screening panels."""

    @classmethod
    def load_yaml_config(cls, config_path: Optional[Path] = None) -> Dict[str, Any]:
        """Loads default parameters from pipeline_config.yaml if present."""
        target_path = config_path or (PROJECT_ROOT / "config" / "pipeline_config.yaml")
        if not target_path.is_file():
            logger.warning(f"Config file not found at {target_path}, using built-in defaults.")
            return {}

        try:
            with open(target_path, "r", encoding="utf-8") as f:
                return yaml.safe_load(f) or {}
        except Exception as exc:
            logger.warning(f"Failed loading YAML config ({exc}), using built-in defaults.")
            return {}

    @classmethod
    def resolve_path(cls, path_str: Optional[str], default_relative: Path) -> Path:
        """Resolves path string relative to PROJECT_ROOT if not already absolute."""
        if not path_str:
            return default_relative if default_relative.is_absolute() else (PROJECT_ROOT / default_relative)
        p = Path(path_str)
        return p if p.is_absolute() else (PROJECT_ROOT / p)

    @classmethod
    def compute_sha256(cls, file_path: Path) -> str:
        """Calculates SHA-256 checksum of a single file."""
        hasher = hashlib.sha256()
        with open(file_path, "rb") as f:
            while chunk := f.read(65536):
                hasher.update(chunk)
        return hasher.hexdigest()

    @classmethod
    def compute_directory_sha256(cls, dir_path: Path) -> str:
        """Calculates combined SHA-256 checksum of all files in a directory."""
        hasher = hashlib.sha256()
        for f in sorted(dir_path.rglob("*")):
            if f.is_file():
                hasher.update(f.name.encode("utf-8"))
                with open(f, "rb") as fp:
                    while chunk := fp.read(65536):
                        hasher.update(chunk)
        return hasher.hexdigest()

    @classmethod
    def _parse_frequency_line_strict(cls, line: str) -> Optional[Tuple[int, Dict[str, float], bool]]:
        """Strictly parses standard schmutzi frequency file lines.

        Expected format: <pos> <count_A> <count_C> <count_G> <count_T>
        Returns:
            Tuple of (pos, base_values_dict, is_normalized_frequency)
        """
        line = line.strip()
        if not line or line.startswith("#"):
            return None

        parts = re.split(r"[\s,\t]+", line)
        if len(parts) < 5 or not parts[0].isdigit():
            return None

        try:
            pos = int(parts[0])
            vals = [float(p) for p in parts[1:5]]
            bases = ["A", "C", "G", "T"]
            total = sum(vals)
            # Determina se i valori sono frequenze normalizzate (somma ~ 1.0) o conteggi interi
            is_freq = math.isclose(total, 1.0, rel_tol=1e-3, abs_tol=1e-3)
            return pos, dict(zip(bases, vals)), is_freq
        except ValueError:
            return None

    @classmethod
    def parse_schmutzi_freq_directory(
        cls,
        freqs_dir: Path,
        strict_mode: bool = True,
    ) -> Tuple[Dict[int, Dict[str, float]], Dict[int, float], int]:
        """Aggregates allele distributions across frequency tables, handling counts vs normalized frequencies."""
        if not freqs_dir.is_dir():
            raise FileNotFoundError(f"Frequencies directory not found: {freqs_dir}")

        freq_files = sorted(
            list(freqs_dir.glob("*.freq"))
            + list(freqs_dir.glob("*.txt"))
            + list(freqs_dir.glob("*.tsv"))
            + [f for f in freqs_dir.iterdir() if f.is_file() and not f.name.startswith(".")]
        )
        unique_files = sorted(list(set(freq_files)))
        if not unique_files:
            raise FileNotFoundError(f"No frequency data files found in: {freqs_dir}")

        logger.info(f"Parsing {len(unique_files)} frequency files from {freqs_dir} (strict_mode={strict_mode})...")

        aggregated_counts: Dict[int, Dict[str, float]] = {}
        total_observations_per_site: Dict[int, float] = {}
        skipped_lines = 0

        for f_path in unique_files:
            with open(f_path, "r", encoding="utf-8", errors="replace") as f:
                for line in f:
                    parsed = cls._parse_frequency_line_strict(line)
                    if not parsed:
                        skipped_lines += 1
                        continue

                    pos, base_vals, is_freq = parsed
                    if pos not in aggregated_counts:
                        aggregated_counts[pos] = {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0}
                        total_observations_per_site[pos] = 0.0

                    weight = 1.0 if not is_freq else 1.0
                    for b, val in base_vals.items():
                        aggregated_counts[pos][b] += val

                    total_observations_per_site[pos] += (sum(base_vals.values()) if not is_freq else 1.0)

        if skipped_lines > 0:
            logger.warning(f"Skipped {skipped_lines} unparseable lines across frequency files.")

        return aggregated_counts, total_observations_per_site, len(unique_files)

    @classmethod
    def calculate_informativeness(cls, ref_base: str, allele_counts: Dict[str, float], total_obs: float) -> Tuple[float, float, str]:
        """Calculates Minor Allele Frequency (MAF) and non-reference informativeness score."""
        if total_obs <= 0:
            return 0.0, 0.0, ""

        alt_candidates = [(b, cnt) for b, cnt in allele_counts.items() if b != ref_base]
        if not alt_candidates:
            return 0.0, 0.0, ""

        alt_candidates.sort(key=lambda x: x[1], reverse=True)
        top_alt_base, top_alt_count = alt_candidates[0]
        total_alt_counts = sum(cnt for _, cnt in alt_candidates)

        contaminant_maf = top_alt_count / float(total_obs)
        # Informativeness: frazione totale di basi non-reference nel database
        informativeness_score = total_alt_counts / float(total_obs)

        return contaminant_maf, informativeness_score, top_alt_base

    @classmethod
    def validate_panel(
        cls,
        diagnostic_sites: List[DiagnosticSite],
        genome_length: int,
        min_panel_size: int = 10,
    ) -> None:
        """Enforces structural, positional, and biological invariants on the generated panel."""
        if len(diagnostic_sites) < min_panel_size:
            raise PanelBuilderValidationError(
                f"Generated panel size ({len(diagnostic_sites)} sites) is below minimum threshold ({min_panel_size})."
            )

        seen_positions: Set[int] = set()
        prev_pos = 0

        for site in diagnostic_sites:
            if site.position in seen_positions:
                raise PanelBuilderValidationError(f"Duplicate position detected in panel: {site.position}")
            seen_positions.add(site.position)

            if site.position <= prev_pos:
                raise PanelBuilderValidationError(f"Panel is not strictly sorted: {site.position} <= {prev_pos}")
            prev_pos = site.position

            if not (1 <= site.position <= genome_length):
                raise PanelBuilderValidationError(f"Position {site.position} out of mitochondrial range [1, {genome_length}].")

            if site.ref_base not in ("A", "C", "G", "T"):
                raise PanelBuilderValidationError(f"Invalid reference base '{site.ref_base}' at position {site.position}")

            if site.primary_nonref_allele not in ("A", "C", "G", "T"):
                raise PanelBuilderValidationError(f"Invalid non-reference allele '{site.primary_nonref_allele}' at position {site.position}")

            if site.primary_nonref_allele == site.ref_base:
                raise PanelBuilderValidationError(f"Non-reference allele matches reference base ({site.ref_base}) at position {site.position}")

            if not (0.0 <= site.contaminant_maf <= 1.0):
                raise PanelBuilderValidationError(f"MAF {site.contaminant_maf} out of bounds at position {site.position}")

            if site.database_allele_observations <= 0:
                raise PanelBuilderValidationError(f"Database observations <= 0 at position {site.position}")

        logger.info(f"Post-condition validation PASSED for {len(diagnostic_sites):,} diagnostic sites.")

    @classmethod
    def build_panel(
        cls,
        freqs_dir: Path,
        reference_fasta: Path,
        expected_contig: str = "NC_012920.1",
        min_maf_threshold: float = 0.02,
        min_total_observations: float = 10.0,
        min_panel_size: int = 10,
        strict_mode: bool = True,
    ) -> Tuple[List[DiagnosticSite], str, str]:
        """Builds and validates the allele-aware candidate diagnostic panel."""
        if not reference_fasta.is_file():
            raise FileNotFoundError(f"Reference FASTA not found: {reference_fasta}")

        ref_sha = cls.compute_sha256(reference_fasta)
        db_sha = cls.compute_directory_sha256(freqs_dir)

        with pysam.FastaFile(str(reference_fasta)) as fa:
            if expected_contig not in fa.references:
                raise ValueError(f"Contig '{expected_contig}' missing from reference FASTA: {fa.references}")
            ref_seq = fa.fetch(expected_contig)

        genome_len = len(ref_seq)
        raw_counts, total_obs_map, n_files = cls.parse_schmutzi_freq_directory(freqs_dir, strict_mode=strict_mode)

        diagnostic_sites: List[DiagnosticSite] = []
        invalid_positions_skipped = 0

        for pos in sorted(raw_counts.keys()):
            if pos < 1 or pos > genome_len:
                invalid_positions_skipped += 1
                continue

            ref_base = ref_seq[pos - 1].upper()
            if ref_base not in ("A", "C", "G", "T"):
                invalid_positions_skipped += 1
                continue

            counts = raw_counts[pos]
            total_obs = total_obs_map.get(pos, 0.0)
            if total_obs < min_total_observations:
                continue

            maf, info_score, primary_nonref = cls.calculate_informativeness(ref_base, counts, total_obs)
            if maf >= min_maf_threshold and primary_nonref:
                diagnostic_sites.append(
                    DiagnosticSite(
                        position=pos,
                        ref_base=ref_base,
                        primary_nonref_allele=primary_nonref,
                        contaminant_maf=maf,
                        database_allele_observations=total_obs,
                        informativeness_score=info_score,
                        allele_breakdown={b: round(c / total_obs, 4) for b, c in counts.items() if total_obs > 0},
                    )
                )

        if invalid_positions_skipped > 0:
            logger.warning(f"Skipped {invalid_positions_skipped} positions out of reference range [1, {genome_len}].")

        cls.validate_panel(diagnostic_sites, genome_len, min_panel_size=min_panel_size)
        return diagnostic_sites, ref_sha, db_sha

    @classmethod
    def export_panel_tsv(
        cls,
        diagnostic_sites: List[DiagnosticSite],
        output_tsv_path: Path,
        reference_sha: str,
        database_sha: str,
        min_maf: float,
        min_total_observations: float,
        version_tag: str = "phylotree17_schmutzi197_v1.0",
    ) -> Path:
        """Serializes diagnostic sites to an audit-ready, self-describing TSV file."""
        output_tsv_path.parent.mkdir(parents=True, exist_ok=True)

        lines: List[str] = [
            f"# diagnostic_panel_version={version_tag}\n",
            f"# builder_version=panel-builder-v2.0\n",
            f"# reference_sha256={reference_sha}\n",
            f"# schmutzi_database_sha256={database_sha}\n",
            f"# min_maf_threshold={min_maf}\n",
            f"# min_total_observations={min_total_observations}\n",
            f"# total_sites={len(diagnostic_sites)}\n",
            "pos\tref_base\tprimary_nonref_allele\tcontaminant_maf\tinformativeness_score\tdatabase_allele_observations\tallele_breakdown\n",
        ]

        for s in diagnostic_sites:
            allele_json = json.dumps(s.allele_breakdown, sort_keys=True)
            lines.append(
                f"{s.position}\t{s.ref_base}\t{s.primary_nonref_allele}\t"
                f"{s.contaminant_maf:.4f}\t{s.informativeness_score:.4f}\t"
                f"{s.database_allele_observations:.1f}\t{allele_json}\n"
            )

        with open(output_tsv_path, "w", encoding="utf-8") as f_out:
            f_out.writelines(lines)

        panel_sha = cls.compute_sha256(output_tsv_path)
        logger.info(f"Exported self-describing panel -> {output_tsv_path} (SHA-256: {panel_sha})")
        return output_tsv_path


def main() -> None:
    """CLI entrypoint resolving configurations from YAML and command-line flags."""
    config_data = DiagnosticPanelBuilder.load_yaml_config()

    ref_cfg = config_data.get("reference", {})
    contam_cfg = config_data.get("contamination", {})
    schmutzi_cfg = contam_cfg.get("schmutzi", {})
    fast_cfg = contam_cfg.get("fast_screening", {})
    builder_cfg = fast_cfg.get("panel_builder", {})

    default_ref_fasta = DiagnosticPanelBuilder.resolve_path(
        ref_cfg.get("fasta_path"),
        PROJECT_ROOT / "data" / "reference" / "rCRS.fasta",
    )
    default_contig = ref_cfg.get("expected_contig", "NC_012920.1")

    db_base = schmutzi_cfg.get("database_path", "data/reference/schmutzi_db/197")
    default_freqs_dir = DiagnosticPanelBuilder.resolve_path(
        f"{db_base}/freqs" if not db_base.endswith("/freqs") else db_base,
        PROJECT_ROOT / "data" / "reference" / "schmutzi_db" / "197" / "freqs",
    )

    default_output_tsv = DiagnosticPanelBuilder.resolve_path(
        fast_cfg.get("diagnostic_panel_path"),
        PROJECT_ROOT / "data" / "reference" / "diagnostic_sites_phylotree17.tsv",
    )
    default_version = fast_cfg.get("diagnostic_panel_version", "phylotree17_schmutzi197_v1.0")

    parser = argparse.ArgumentParser(
        description="Build allele-aware, self-describing diagnostic panel from schmutzi database."
    )
    parser.add_argument("--config", type=Path, default=PROJECT_ROOT / "config" / "pipeline_config.yaml", help="Path to YAML config.")
    parser.add_argument("--freqs-dir", type=Path, default=default_freqs_dir, help=f"Path to database freqs (default: {default_freqs_dir}).")
    parser.add_argument("--reference-fasta", type=Path, default=default_ref_fasta, help=f"Path to reference FASTA (default: {default_ref_fasta}).")
    parser.add_argument("--contig", type=str, default=default_contig, help=f"Target mitochondrial contig (default: {default_contig}).")
    parser.add_argument("--min-maf", type=float, default=builder_cfg.get("min_maf", 0.02), help="Minimum database MAF to retain site (default: 0.02).")
    parser.add_argument("--min-total-observations", type=float, default=builder_cfg.get("min_total_observations", 10.0), help="Minimum database observations per site (default: 10.0).")
    parser.add_argument("--min-panel-size", type=int, default=builder_cfg.get("min_panel_size", 10), help="Minimum acceptable panel size (default: 10).")
    parser.add_argument("--output-tsv", type=Path, default=default_output_tsv, help=f"Destination TSV path (default: {default_output_tsv}).")
    parser.add_argument("--version-tag", type=str, default=default_version, help=f"Version tag for panel metadata (default: {default_version}).")
    parser.add_argument("--strict-mode", action="store_true", default=True, help="Enforce strict frequency line parsing (default: True).")

    args = parser.parse_args()

    sites, ref_sha, db_sha = DiagnosticPanelBuilder.build_panel(
        freqs_dir=args.freqs_dir,
        reference_fasta=args.reference_fasta,
        expected_contig=args.contig,
        min_maf_threshold=args.min_maf,
        min_total_observations=args.min_total_observations,
        min_panel_size=args.min_panel_size,
        strict_mode=args.strict_mode,
    )

    DiagnosticPanelBuilder.export_panel_tsv(
        diagnostic_sites=sites,
        output_tsv_path=args.output_tsv,
        reference_sha=ref_sha,
        database_sha=db_sha,
        min_maf=args.min_maf,
        min_total_observations=args.min_total_observations,
        version_tag=args.version_tag,
    )


if __name__ == "__main__":
    main()
