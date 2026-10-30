#!/usr/bin/env python3
"""Thin CLI wrapper for the modular mtDNA analysis pipeline."""

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.pipeline.cli import main

if __name__ == "__main__":
    main()