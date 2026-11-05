"""
Utility components for the ancient & modern mtDNA pipeline.

Exposes configuration management (PipelineConfig, load_pipeline_config),
logging infrastructure (setup_logger, configure_root_logger, CustomColorFormatter),
system/reference validation (EnvironmentValidator),
and gradient progress bars (GradientProgressBar, SubstageTracker, run_tracked_subprocess).
"""

from src.utils.config_loader import (
    DEFAULT_CONFIG_PATH,
    PipelineConfig,
    load_pipeline_config,
    resolve_config_path,
)
from src.utils.logger import (
    CustomColorFormatter,
    configure_root_logger,
    get_logger,
    setup_logger,
)
from src.utils.progress_bar import (
    BOLD_GREEN,
    BOLD_RED,
    BOLD_WHITE,
    DIM_GREY,
    RESET,
    GradientProgressBar,
    SubstageTracker,
    format_bytes,
    format_duration,
    get_gradient_color,
    run_tracked_subprocess,
    strip_ansi,
    supports_color,
)
from src.utils.validator import EnvironmentValidator

__all__ = [
    # Config
    "DEFAULT_CONFIG_PATH",
    "PipelineConfig",
    "load_pipeline_config",
    "resolve_config_path",
    # Logger
    "CustomColorFormatter",
    "configure_root_logger",
    "get_logger",
    "setup_logger",
    # Validator
    "EnvironmentValidator",
    # Progress Bar
    "GradientProgressBar",
    "SubstageTracker",
    "run_tracked_subprocess",
    "get_gradient_color",
    "strip_ansi",
    "supports_color",
    "format_bytes",
    "format_duration",
    "RESET",
    "DIM_GREY",
    "BOLD_WHITE",
    "BOLD_GREEN",
    "BOLD_RED",
]
