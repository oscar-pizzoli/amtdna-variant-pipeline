"""
End-to-End Pipeline Orchestration and Stage Tracking Package.
"""

from src.pipeline.stage_tracker import StageTracker


def __getattr__(name: str):
    if name in ("EndToEndPipelineRunner", "PipelineExecutionError"):
        from src.pipeline import runner
        return getattr(runner, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = ["StageTracker", "EndToEndPipelineRunner", "PipelineExecutionError"]
