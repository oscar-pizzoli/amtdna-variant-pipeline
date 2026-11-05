"""
Centralized logging and color formatting module for the mtDNA analysis pipeline.

Provides dual-channel logging (colorized ANSI terminal output and persistent, clean file logs),
proper level-based styling, and automatic classification of HTSlib/samtools messages.
"""

import logging
from pathlib import Path
import sys
from typing import Optional


class CustomColorFormatter(logging.Formatter):
    """Logging Formatter that injects ANSI colors for terminal display."""

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

        msg_str = str(record.getMessage())
        # Highlight HTSlib / samtools warning and error markers
        if msg_str.startswith("[W::") or "[WARNING]" in msg_str:
            color = self.YELLOW
        elif msg_str.startswith("[E::") or "[CRITICAL]" in msg_str:
            color = self.BOLD_RED
        else:
            color = self.LEVEL_COLORS.get(record.levelno)

        if color:
            return f"{color}{rendered}{self.RESET}"
        return rendered


def configure_root_logger(
    level: int = logging.INFO,
    log_file: Optional[Path] = None,
    format_str: Optional[str] = None,
) -> logging.Logger:
    """Configures the root logger with colored terminal output and optional persistent file logging.

    Args:
        level: Logging verbosity level (default: logging.INFO).
        log_file: Optional filesystem path for clean persistent log storage.
        format_str: Optional format string override.

    Returns:
        Configured root logging.Logger instance.
    """
    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    # Clean existing handlers to prevent duplicate lines
    if root_logger.hasHandlers():
        root_logger.handlers.clear()

    # Channel 1: Console (STDOUT) with ANSI colors
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(CustomColorFormatter(fmt=format_str))
    root_logger.addHandler(console_handler)

    # Channel 2: File (Append-mode) without ANSI escapes
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_fmt = format_str or CustomColorFormatter.DEFAULT_FORMAT
        file_formatter = logging.Formatter(fmt=file_fmt, datefmt="%Y-%m-%d %H:%M:%S")
        file_handler = logging.FileHandler(log_path, encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(file_formatter)
        root_logger.addHandler(file_handler)

    return root_logger


def setup_logger(
    name: str = "mtDNA_Pipeline",
    log_file: Optional[Path] = None,
    level: int = logging.INFO,
) -> logging.Logger:
    """Configures and returns a dedicated dual-channel logger instance.

    Args:
        name: Name of the logger domain.
        log_file: Optional filesystem path for clean persistent log storage.
        level: Logging verbosity level (default: logging.INFO).

    Returns:
        Configured logging.Logger instance.
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False

    if logger.hasHandlers():
        logger.handlers.clear()

    # Channel 1: Console (STDOUT) with ANSI colors
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(CustomColorFormatter())
    logger.addHandler(console_handler)

    # Channel 2: Persistent File
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_formatter = logging.Formatter(
            fmt="%(asctime)s | %(levelname)-8s | [%(name)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        file_handler = logging.FileHandler(log_path, encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

    return logger


def get_logger(name: str = "mtDNA_Pipeline") -> logging.Logger:
    """Convenience getter for existing loggers within the mtDNA pipeline domain."""
    return logging.getLogger(name)
