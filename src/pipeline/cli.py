"""
Command Line Interface (CLI) parser and entrypoint for the mtDNA analysis pipeline.
"""

import argparse
import logging
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from src.pipeline.runner import EndToEndPipelineRunner

from src.utils.logger import configure_root_logger


def parse_arguments() -> argparse.Namespace:
    """Parse command line options."""
    parser = argparse.ArgumentParser(
        description="End-to-End Ancient and Modern mtDNA Consensus, Damage, and Haplogrouping Pipeline."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=PROJECT_ROOT / "config" / "default_config.yaml",
        help="Path to YAML configuration file.",
    )
    parser.add_argument(
        "--data-type",
        type=str,
        default=None,
        choices=["BAM", "FASTQ"],
        help="Override input sequence format (BAM or FASTQ).",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=None,
        help="Input sequence directory (FASTQ or BAM). Overrides config input_dir.",
    )
    parser.add_argument(
        "--cpu-budget-fraction",
        type=float,
        default=None,
        help="Fraction of total host CPU cores allocated to pipeline execution.",
    )
    parser.add_argument(
        "--host-cores",
        type=int,
        default=None,
        help="Manual host core count override.",
    )
    parser.add_argument(
        "--cache-mode",
        type=str,
        default=None,
        choices=["prompt", "reuse", "overwrite"],
        help="Override YAML output.cache_mode ('prompt', 'reuse', or 'overwrite').",
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging verbosity level (default: INFO).",
    )
    return parser.parse_args()


def resolve_input_directory(orchestrator: "EndToEndPipelineRunner", cli_input_dir: Path | None) -> Path:
    """Resolve physical input directory based on input format and configuration."""
    if cli_input_dir:
        return cli_input_dir

    yaml_input_dir = orchestrator.input_cfg.get("input_dir")
    if yaml_input_dir:
        return orchestrator._resolve_path(yaml_input_dir, PROJECT_ROOT / "data" / "raw_bam")

    if orchestrator.data_type == "BAM":
        target_dir = orchestrator.input_cfg.get("raw_bam_dir", "data/raw_bam")
        return orchestrator._resolve_path(target_dir, PROJECT_ROOT / "data" / "raw_bam")
    else:
        target_dir = orchestrator.input_cfg.get("raw_fastq_dir", "data/raw_fastq")
        return orchestrator._resolve_path(target_dir, PROJECT_ROOT / "data" / "raw_fastq")


def main() -> None:
    """Main CLI execution flow."""
    args = parse_arguments()

    numeric_level = getattr(logging, args.log_level.upper(), logging.INFO)
    configure_root_logger(level=numeric_level)

    from src.pipeline.runner import EndToEndPipelineRunner, PipelineExecutionError

    try:
        orchestrator = EndToEndPipelineRunner(
            config_path=args.config,
            cpu_fraction_override=args.cpu_budget_fraction,
            host_cores_override=args.host_cores,
            data_type_override=args.data_type,
            cache_mode_override=args.cache_mode,
        )

        input_directory = resolve_input_directory(orchestrator, args.input_dir)
        orchestrator.execute_cohort(input_dir=input_directory)
    except KeyboardInterrupt:
        logging.warning("Execution interrupted by user (SIGINT).")
        sys.exit(130)
    except PipelineExecutionError as exc:
        logging.critical(f"Execution halted due to pipeline failure: {exc}")
        sys.exit(1)
    except Exception as exc:
        logging.critical(f"Unhandled fatal error: {exc}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()