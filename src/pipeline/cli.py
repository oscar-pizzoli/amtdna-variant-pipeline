"""
Command Line Interface (CLI) parser and entrypoint for the mtDNA analysis pipeline.
"""

import argparse
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.pipeline.runner import EndToEndPipelineRunner


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
    return parser.parse_args()


def resolve_input_directory(orchestrator: EndToEndPipelineRunner, cli_input_dir: Path | None) -> Path:
    """Resolve physical input directory based on input format and configuration."""
    if cli_input_dir:
        return cli_input_dir

    yaml_input_dir = orchestrator.input_cfg.get("input_dir")
    if orchestrator.data_type == "BAM":
        if yaml_input_dir and "fastq" not in str(yaml_input_dir).lower() and Path(yaml_input_dir).is_dir():
            return orchestrator._resolve_path(yaml_input_dir, PROJECT_ROOT / "data" / "raw_bam")
        if (PROJECT_ROOT / "data" / "raw_bam").is_dir():
            return PROJECT_ROOT / "data" / "raw_bam"
        return orchestrator._resolve_path(yaml_input_dir, PROJECT_ROOT / "data" / "bam")

    return orchestrator._resolve_path(yaml_input_dir, PROJECT_ROOT / "data" / "raw_fastq")


def main() -> None:
    """Main CLI execution flow."""
    args = parse_arguments()
    orchestrator = EndToEndPipelineRunner(
        config_path=args.config,
        cpu_fraction_override=args.cpu_budget_fraction,
        host_cores_override=args.host_cores,
    )

    input_directory = resolve_input_directory(orchestrator, args.input_dir)
    orchestrator.execute_cohort(input_dir=input_directory)


if __name__ == "__main__":
    main()
