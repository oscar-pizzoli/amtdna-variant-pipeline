"""Centralized logging module for the mtDNA pipeline."""

import logging
import sys
from pathlib import Path
from typing import Optional


class CustomColorFormatter(logging.Formatter):
    """Logging Formatter to add ANSI colors for terminal output."""

    GREY = "\x1b[38;20m"
    GREEN = "\x1b[32;20m"
    YELLOW = "\x1b[33;20m"
    RED = "\x1b[31;20m"
    BOLD_RED = "\x1b[31;1m"
    RESET = "\x1b[0m"
    FORMAT = "%(asctime)s | %(levelname)-8s | [%(name)s] %(message)s"

    FORMATS = {
        logging.DEBUG: GREY + FORMAT + RESET,
        logging.INFO: GREEN + FORMAT + RESET,
        logging.WARNING: YELLOW + FORMAT + RESET,
        logging.ERROR: RED + FORMAT + RESET,
        logging.CRITICAL: BOLD_RED + FORMAT + RESET,
    }

    def format(self, record: logging.LogRecord) -> str:
        log_fmt = self.FORMATS.get(record.levelno, self.FORMAT)
        formatter = logging.Formatter(log_fmt, datefmt="%Y-%m-%d %H:%M:%S")
        return formatter.format(record)


def setup_logger(
    name: str = "mtDNA_Pipeline",
    log_file: Optional[Path] = None,
    level: int = logging.INFO,
) -> logging.Logger:
    """Configures and returns a dual-channel logger instance.

    Args:
        name: Name of the logger domain.
        log_file: Optional filesystem path for persistent log storage.
        level: Logging verbosity level (default: logging.INFO).

    Returns:
        Configured logging.Logger instance.
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False

    # Evita la registrazione di duplicati di handler se richiamato più volte
    if logger.hasHandlers():
        logger.handlers.clear()

    # Channel 1: Console (STDOUT) con formattazione colorata
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(CustomColorFormatter())
    logger.addHandler(console_handler)

    # Channel 2: File (Append-mode) con formattazione standard deterministica
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_formatter = logging.Formatter(
            fmt="%(asctime)s | %(levelname)-8s | [%(name)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

    return logger
