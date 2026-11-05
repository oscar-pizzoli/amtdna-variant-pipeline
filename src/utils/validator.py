"""
System dependency and reference genome validator for mtDNA analysis pipeline.

Verifies the presence of required bioinformatic binaries (samtools, bwa, AdapterRemoval,
mapDamage, pmdtools, schmutzi, haplogrep) and validates FASTA reference genome indexation.
"""

import logging
from pathlib import Path
import shutil
import subprocess
from typing import Dict, List, Optional, Set

logger = logging.getLogger("mtDNA_Pipeline.Validator")


class EnvironmentValidator:
    """Validates external command-line dependencies and reference genome integrity."""

    # Mandatory tools for core mapping and processing
    REQUIRED_CORE_TOOLS: Set[str] = {
        "samtools",
        "bwa",
    }

    # Modular tools required across ancient DNA, damage profiling, and haplogrouping
    ANCIENT_DNA_TOOLS: Set[str] = {
        "AdapterRemoval",
        "mapDamage",
        "pmdtools",
    }

    OPTIONAL_TOOLS: Set[str] = {
        "schmutzi.pl",
        "contDeam.pl",
        "haplogrep",
        "bcftools",
        "bedtools",
    }

    @classmethod
    def check_binaries(
        cls,
        check_adna: bool = True,
        raise_on_missing: bool = True,
    ) -> List[str]:
        """Verifies that mandatory command-line tools are accessible in PATH.

        Args:
            check_adna: Whether to also check dedicated ancient DNA analytical tools.
            raise_on_missing: Whether to raise RuntimeError if required core tools are missing.

        Returns:
            List of missing core tools.

        Raises:
            RuntimeError: If raise_on_missing is True and any required core tool is absent.
        """
        missing_core: List[str] = []
        for tool in sorted(cls.REQUIRED_CORE_TOOLS):
            if shutil.which(tool) is None:
                missing_core.append(tool)

        if check_adna:
            for tool in sorted(cls.ANCIENT_DNA_TOOLS):
                if shutil.which(tool) is None:
                    logger.warning(f"Ancient DNA tool not detected in PATH: '{tool}' (check conda environment)")

        if missing_core:
            err_msg = f"Missing critical system dependencies in PATH: {', '.join(missing_core)}"
            logger.error(err_msg)
            if raise_on_missing:
                raise RuntimeError(err_msg)
            return missing_core

        logger.info("Core external dependencies (samtools, bwa) verified successfully.")
        return []

    @classmethod
    def inspect_toolset(cls) -> Dict[str, bool]:
        """Returns a boolean availability map for all pipeline-related binaries."""
        all_tools = cls.REQUIRED_CORE_TOOLS | cls.ANCIENT_DNA_TOOLS | cls.OPTIONAL_TOOLS
        return {tool: shutil.which(tool) is not None for tool in sorted(all_tools)}

    @staticmethod
    def validate_reference(
        reference_fasta: Path,
        auto_index: bool = True,
    ) -> None:
        """Validates FASTA reference existence and verifies/generates required indices (.fai, .bwt).

        Args:
            reference_fasta: Path to the reference FASTA file.
            auto_index: Automatically generate missing samtools (.fai) and BWA indices.

        Raises:
            FileNotFoundError: If the reference FASTA does not exist or indices are missing.
            RuntimeError: If required indexing binaries are missing or fail during execution.
        """
        ref_path = Path(reference_fasta).resolve()
        if not ref_path.is_file():
            raise FileNotFoundError(f"Reference FASTA not found at path: {ref_path}")

        fai_file = ref_path.with_suffix(ref_path.suffix + ".fai")
        bwa_bwt = ref_path.with_suffix(ref_path.suffix + ".bwt")

        # 1. FASTA Index (.fai) via samtools faidx
        if not fai_file.is_file():
            if auto_index:
                if shutil.which("samtools") is None:
                    raise RuntimeError("Cannot index reference: 'samtools' binary not found in PATH.")
                logger.info(f"Generating missing FASTA index (.fai) for: {ref_path.name}")
                try:
                    subprocess.run(
                        ["samtools", "faidx", str(ref_path)],
                        capture_output=True,
                        text=True,
                        check=True,
                    )
                except subprocess.CalledProcessError as exc:
                    raise RuntimeError(f"Failed to generate samtools faidx index: {exc.stderr}") from exc
            else:
                raise FileNotFoundError(f"Missing FASTA index (.fai): {fai_file}")

        # 2. BWA Index (.bwt, .pac, .ann, .amb, .sa)
        if not bwa_bwt.is_file():
            if auto_index:
                if shutil.which("bwa") is None:
                    raise RuntimeError("Cannot index reference: 'bwa' binary not found in PATH.")
                logger.info(f"Generating missing BWA index for: {ref_path.name}")
                try:
                    subprocess.run(
                        ["bwa", "index", str(ref_path)],
                        capture_output=True,
                        text=True,
                        check=True,
                    )
                except subprocess.CalledProcessError as exc:
                    raise RuntimeError(f"Failed to generate BWA index: {exc.stderr}") from exc
            else:
                raise FileNotFoundError(f"Missing BWA index (.bwt): {bwa_bwt}")

        logger.info("Reference genome and indices (.fai, .bwt) validated.")
