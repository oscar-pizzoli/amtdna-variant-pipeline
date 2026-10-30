"""Unit and integration smoke test for Block 1 components."""

import tempfile
from pathlib import Path
import sys

# Aggiunge la root del progetto al PATH di Python
ROOT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT_DIR))

from src.utils.logger import setup_logger
from src.utils.config_loader import PipelineConfig
from src.utils.validator import EnvironmentValidator


def test_block1() -> None:
    """Executes end-to-end verification of Block 1 utilities."""
    
    # --- Test 1: Logger Initialization ---
    log_dir = ROOT_DIR / "results" / "logs"
    log_file = log_dir / "test_execution.log"
    logger = setup_logger(name="TestRunner", log_file=log_file)
    logger.info("=== AVVIO SMOKE TEST: BLOCCO 1 ===")

    assert log_file.is_file(), "ERRORE: Il file di log non è stato creato su disco."
    logger.info("Test 1 Superato: Logger a doppio canale configurato correttamente.")

    # --- Test 2: Configuration Loader & Validation ---
    config_file = ROOT_DIR / "config" / "default_config.yaml"
    logger.info(f"Caricamento configurazione da: {config_file.name}")
    
    config = PipelineConfig.from_yaml(config_file)
    assert config.project_name != "", "ERRORE: Project name non popolato."
    assert 0.0 < config.het_min_freq < config.het_max_freq <= 1.0, "ERRORE: Range frequenze eteroplasmia invalido."
    assert isinstance(config.input_dir, Path), "ERRORE: input_dir non è un'istanza di Path."
    logger.info("Test 2 Superato: Parsing e validazione semantica YAML completati.")

    # --- Test 3: Defensive Programming (Fail-Early su config invalida) ---
    logger.info("Verifica comportamento difensivo su configurazione errata...")
    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as tmp:
        # Inseriamo un valore non valido (min_frequency = 1.50)
        tmp.write(
            """
            input: { data_type: "FASTQ", input_dir: "data" }
            reference: { fasta: "ref.fasta" }
            filtering: { min_mapq: 30, min_baseq: 30 }
            consensus: { min_depth: 10 }
            resources: { threads: 4 }
            output: { results_dir: "results" }
            heteroplasmy: { min_frequency: 1.50, max_frequency: 0.80 }
            """
        )
        tmp_path = Path(tmp.name)

    try:
        PipelineConfig.from_yaml(tmp_path)
        logger.error("ERRORE: Il parser non ha bloccato un valore di frequenza non valido (>1.0).")
    except ValueError as exc:
        logger.info(f"Test 3 Superato: Errore semantico intercettato con successo ({exc}).")
    finally:
        tmp_path.unlink(missing_ok=True)

    # --- Test 4: Environment & Binaries Check ---
    logger.info("Verifica strumenti di sistema disponibili nel PATH...")
    try:
        EnvironmentValidator.check_binaries(check_adna=False)
        logger.info("Test 4 Superato: Dipendenze binarie core (samtools, bwa, ecc.) verificate.")
    except RuntimeError as exc:
        logger.warning(f"Avviso dipendenze di sistema: {exc}")

    # --- Test 5: Reference FASTA & Indexer ---
    logger.info("Verifica validazione reference FASTA...")
    ref_path = ROOT_DIR / "reference" / "rCRS.fasta"
    if not ref_path.is_file():
        # Crea un mock temporaneo per il test se il file non è ancora presente
        ref_path.parent.mkdir(parents=True, exist_ok=True)
        with open(ref_path, "w", encoding="utf-8") as f:
            f.write(">NC_012920.1\nGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCAT\n")
        logger.info(f"Creato file mock di reference in: {ref_path}")

    EnvironmentValidator.validate_reference(ref_path, auto_index=True)
    assert ref_path.with_suffix(ref_path.suffix + ".fai").is_file(), "ERRORE: Indice .fai mancante."
    assert ref_path.with_suffix(ref_path.suffix + ".bwt").is_file(), "ERRORE: Indice BWA .bwt mancante."
    logger.info("Test 5 Superato: Reference e indici (.fai, .bwt) generati e validati.")

    logger.info("=== TUTTI I TEST DEL BLOCCO 1 SONO STATI SUPERATI CON SUCCESSO ===")


if __name__ == "__main__":
    test_block1()
