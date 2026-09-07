"""System dependency and input integrity validator."""

import shutil
import subprocess
from pathlib import Path
from typing import List, Set
import logging

logger = logging.getLogger("mtDNA_Pipeline.Validator")


class EnvironmentValidator:
    """Validates external binaries, reference indices, and dataset inputs."""

    REQUIRED_CORE_TOOLS: Set[str] = {
        "samtools",
        "bcftools",
        "bedtools",
        "bwa",
    }

    OPTIONAL_ADNA_TOOLS: Set[str] = {
        "mapDamage",
        "pmdtools",
        "fastp",
        "SeqPrep",
    }

    @classmethod
    def check_binaries(cls, check_adna: bool = True) -> None:
        """Verifies that mandatory command-line tools are available in PATH.

        Args:
            check_adna: Whether to validate optional aDNA tools.

        Raises:
            RuntimeError: If any mandatory binary is missing from PATH.
        """
        missing_tools: List[str] = []
        tools_to_check = set(cls.REQUIRED_CORE_TOOLS)
        if check_adna:
            tools_to_check.update(cls.OPTIONAL_ADNA_TOOLS)

        for tool in tools_to_check:
            if shutil.which(tool) is None:
                if tool in cls.REQUIRED_CORE_TOOLS:
                    missing_tools.append(tool)
                else:
                    logger.warning(f"Optional tool not found in PATH: '{tool}'")

        if missing_tools:
            err_msg = f"Missing critical system dependencies in PATH: {', '.join(missing_tools)}"
            logger.error(err_msg)
            raise RuntimeError(err_msg)

        logger.info("All mandatory external dependencies verified successfully.")

    @staticmethod
    def validate_reference(reference_fasta: Path, auto_index: bool = True) -> None:
        """Validates FASTA reference existence and verifies/generates required indices.

        Args:
            reference_fasta: Path to reference FASTA file.
            auto_index: Automatically generate missing samtools (.fai) and BWA indices.

        Raises:
            FileNotFoundError: If the reference FASTA does not exist.
            RuntimeError: If binary is missing or index generation fails.
        """
        ref_path = Path(reference_fasta).resolve()
        if not ref_path.is_file():
            raise FileNotFoundError(f"Reference FASTA not found at path: {ref_path}")

        fai_file = ref_path.with_suffix(ref_path.suffix + ".fai")
        bwa_bwt = ref_path.with_suffix(ref_path.suffix + ".bwt")

        # 1. Verifica/Generazione Indice FASTA (.fai)
        if not fai_file.is_file():
            if auto_index:
                if shutil.which("samtools") is None:
                    raise RuntimeError("Cannot index reference: 'samtools' binary not found in PATH.")
                logger.info(f"Generating missing FASTA index (.fai) for: {ref_path.name}")
                try:
                    res = subprocess.run(["samtools", "faidx", str(ref_path)], capture_output=True, text=True, check=True)
                except subprocess.CalledProcessError as exc:
                    raise RuntimeError(f"Failed to generate samtools faidx index: {exc.stderr}") from exc
            else:
                raise FileNotFoundError(f"Missing FASTA index: {fai_file}")

        # 2. Verifica/Generazione Indice BWA (.bwt, .sa, ecc.)
        if not bwa_bwt.is_file():
            if auto_index:
                if shutil.which("bwa") is None:
                    raise RuntimeError("Cannot index reference: 'bwa' binary not found in PATH.")
                logger.info(f"Generating missing BWA index for: {ref_path.name}")
                try:
                    res = subprocess.run(["bwa", "index", str(ref_path)], capture_output=True, text=True, check=True)
                except subprocess.CalledProcessError as exc:
                    raise RuntimeError(f"Failed to generate BWA index: {exc.stderr}") from exc
            else:
                raise FileNotFoundError(f"Missing BWA index: {bwa_bwt}")

        logger.info("Reference genome and indices validated.")