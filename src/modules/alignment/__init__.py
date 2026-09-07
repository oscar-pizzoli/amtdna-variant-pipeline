"""Alignment and Filtering Engine package (Block 3)."""

from src.modules.alignment.models import (
    AlignerType,
    AlignmentErrorCode,
    AlignmentFilterPolicy,
    AlignmentMetrics,
    AlignmentOutputs,
    AlignmentResult,
    AlignmentStatus,
    AlignmentToolMetadata,
    AlignmentWarningCode,
    CoverageMetrics,
    DedupMetrics,
    DedupMode,
    MAPQProfile,
    MultiMappingPolicy,
    MultiMappingStrategy,
    ProvenanceMetrics,
    ReadGroupInfo,
    ReadGroupSource,
    ReferenceInfo,
)

__all__ = [
    # Enums & Status Codes
    "AlignerType",
    "DedupMode",
    "AlignmentStatus",
    "AlignmentWarningCode",
    "AlignmentErrorCode",
    "ReadGroupSource",
    "MultiMappingStrategy",
    # Domain Models & Policies
    "MultiMappingPolicy",
    "AlignmentFilterPolicy",
    "ReadGroupInfo",
    "ReferenceInfo",
    "MAPQProfile",
    "CoverageMetrics",
    "DedupMetrics",
    "ProvenanceMetrics",
    "AlignmentOutputs",
    "AlignmentMetrics",
    "AlignmentToolMetadata",
    "AlignmentResult",
]
