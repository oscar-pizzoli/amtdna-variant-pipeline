"""Centralized logging module for the mtDNA pipeline."""

import logging
import sys
from pathlib import Path
from typing import Optional


class CustomColorFormatter(logging.Formatter):
    """Logging Formatter to add ANSI colors for terminal output."""

    GREY = "\x1b[90m"
    CYAN = "\x1b[36m"
    YELLOW = "\x1b[33m"
    RED = "\x1b[31m"
    BOLD_RED = "\x1b[1;31m"
    RESET = "\x1b[0m"

    DEFAULT_FORMAT = "%(asctime)s [%(levelname)s] %(message)s"
    DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

    LEVEL_COLORS = {
        logging.DEBUG: GREY,
        logging.INFO: CYAN,
        logging.WARNING: YELLOW,
        logging.ERROR: RED,
        logging.CRITICAL: BOLD_RED,
    }

    def __init__(
        self,
        fmt: Optional[str] = None,
        datefmt: Optional[str] = None,
    ):
        base_fmt = fmt or self.DEFAULT_FORMAT
        base_datefmt = datefmt or self.DATE_FORMAT
        super().__init__(fmt=base_fmt, datefmt=base_datefmt)
        self.base_fmt = base_fmt
        self.date_fmt = base_datefmt

    def format(self, record: logging.LogRecord) -> str:
        base_formatter = logging.Formatter(self.base_fmt, datefmt=self.date_fmt)
        rendered = base_formatter.format(record)
        color = self.LEVEL_COLORS.get(record.levelno)
        if color:
            return f"{color}{rendered}{self.RESET}"
        return rendered


def configure_root_logger(
    level: int = logging.INFO,
    log_file: Optional[Path] = None,
    format_str: Optional[str] = None,
) -> logging.Logger:
    """Configures the root logger with colored terminal output and optional persistent file log.

    Args:
        level: Logging verbosity level (default: logging.INFO).
        log_file: Optional filesystem path for persistent log storage.
        format_str: Optional format string override.

    Returns:
        Configured root logging.Logger instance.
    """
    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    if root_logger.hasHandlers():
        root_logger.handlers.clear()

    # Channel 1: Console (STDOUT) con formattazione colorata ANSI
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(CustomColorFormatter(fmt=format_str))
    root_logger.addHandler(console_handler)

    # Channel 2: File persistente senza sequenze ANSI
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_fmt = format_str or CustomColorFormatter.DEFAULT_FORMAT
        file_formatter = logging.Formatter(fmt=file_fmt, datefmt="%Y-%m-%d %H:%M:%S")
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(file_formatter)
        root_logger.addHandler(file_handler)

    return root_logger


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
