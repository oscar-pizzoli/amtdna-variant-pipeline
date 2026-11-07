"""
Mitochondrial DNA (mtDNA) Variant Calling, Authenticity and Haplogrouping Pipeline.

Modular step architecture:
- utils: Unified configuration, logging, progress bars, and shared helpers
- steps.reference: Step 1 Reference genome verification & indexing
- steps.discovery: Step 2 Preprocessing & sample discovery
- steps.alignment: Step 3 BWA mapping & deduplication
- steps.damage: Step 4 Damage quantification & rescaling
- steps.contamination: Step 5 Two-tier contamination estimation
- steps.consensus: Step 6 Consensus assembly & heteroplasmy calling
- steps.haplogroup: Step 7 Haplogrouping & phylogenetic audit
- steps.synthesis: Step 8 Forensic synthesis & cohort deliverables
"""

from src import utils
from src.steps import (
    alignment,
    consensus,
    contamination,
    damage,
    discovery,
    haplogroup,
    reference,
    synthesis,
)

__all__ = [
    "utils",
    "reference",
    "discovery",
    "alignment",
    "damage",
    "contamination",
    "consensus",
    "haplogroup",
    "synthesis",
]
