"""
Enumerations and vocabulary definitions for the consensus and heteroplasmy engine.
"""

from enum import Enum


class InputRole(str, Enum):
    """Execution input categorization resolved upstream."""
    RESCALED_ANCIENT = "RESCALED_ANCIENT"
    FILTERED_DEDUP_MODERN = "FILTERED_DEDUP_MODERN"
    UNKNOWN = "UNKNOWN"


class GatingAction(str, Enum):
    """Exhaustive gating actions bridging Block 5 triage and Block 6 calling."""
    STANDARD_CALLING = "STANDARD_CALLING"
    RESTRICT_HETEROPLASMY = "RESTRICT_HETEROPLASMY"
    PMD_FILTERING_REQUIRED = "PMD_FILTERING_REQUIRED"
    USE_ENDOGENOUS_SCHMUTZI_CONSENSUS = "USE_ENDOGENOUS_SCHMUTZI_CONSENSUS"
    EXCLUDE_OR_MANUAL_REVIEW = "EXCLUDE_OR_MANUAL_REVIEW"
    REJECT_CONTAMINATED = "REJECT_CONTAMINATED"
    INSUFFICIENT_DATA = "INSUFFICIENT_DATA"


class ContaminationEstimateSource(str, Enum):
    """Source of the contamination point estimate and decision."""
    SCHMUTZI = "SCHMUTZI"
    FAST_SCREEN = "FAST_SCREEN"
    DISSENT_SCREENING = "DISSENT_SCREENING"
    NONE = "NONE"


class CallStatus(str, Enum):
    """Operational call status for heteroplasmy sites."""
    PASS = "PASS"
    FILTERED = "FILTERED"


class ConsensusCallStatus(str, Enum):
    """Decision state for consensus base calling at a specific site."""
    CALLED = "CALLED"
    LOW_DEPTH = "LOW_DEPTH"
    AMBIGUOUS = "AMBIGUOUS"
    NO_USABLE_BASES = "NO_USABLE_BASES"


class ConsensusCallReason(str, Enum):
    """Structured diagnostic reason for consensus site classification."""
    CALLED_MAJORITY = "CALLED_MAJORITY"
    NO_USABLE_BASES = "NO_USABLE_BASES"
    LOW_DEPTH = "LOW_DEPTH"
    AMBIGUOUS_IUPAC = "AMBIGUOUS_IUPAC"
    SECOND_ALLELE_TOO_LOW = "SECOND_ALLELE_TOO_LOW"
    REFERENCE_N_PLACEHOLDER = "REFERENCE_N_PLACEHOLDER"


class FilterReason(str, Enum):
    """Audit reasons for site exclusion or variant filtering."""
    PASS = "PASS"
    FILTER_NO_DATA = "FILTER_NO_DATA"
    FILTER_LOW_DEPTH = "FILTER_LOW_DEPTH"
    FILTER_LOW_ALT_READS = "FILTER_LOW_ALT_READS"
    FILTER_VAF_RANGE = "FILTER_VAF_RANGE"
    FILTER_STRAND_BIAS = "FILTER_STRAND_BIAS"
    FILTER_CONTAMINATION = "FILTER_CONTAMINATION"
    FILTER_INSUFFICIENT_DATA = "FILTER_INSUFFICIENT_DATA"
    FILTER_DAMAGE_SUSPECT = "FILTER_DAMAGE_SUSPECT"
    FILTER_AMBIGUOUS = "FILTER_AMBIGUOUS"
    FILTER_REFERENCE_MISMATCH = "FILTER_REFERENCE_MISMATCH"
    FILTER_EXCLUDED_INTERVAL = "FILTER_EXCLUDED_INTERVAL"


class IUPACAmbiguity(str, Enum):
    """Standard IUPAC nucleotide ambiguity codes and deterministic resolution."""
    A = "A"
    C = "C"
    G = "G"
    T = "T"
    R = "R"  # A or G (Purine)
    Y = "Y"  # C or T (Pyrimidine)
    S = "S"  # G or C (Strong)
    W = "W"  # A or T (Weak)
    K = "K"  # G or T (Keto)
    M = "M"  # A or C (Amino)
    N = "N"  # Any / Unspecified

    @classmethod
    def resolve(cls, base1: str, base2: str) -> str:
        """
        Deterministically resolve two bases into a canonical base or IUPAC code.
        Handles homozygous inputs cleanly without collapsing into 'N'.
        """
        b1, b2 = base1.upper(), base2.upper()
        pair = frozenset({b1, b2})

        if len(pair) == 1:
            single_base = next(iter(pair))
            return single_base if single_base in {"A", "C", "G", "T"} else "N"

        lookup = {
            frozenset({"A", "G"}): "R",
            frozenset({"C", "T"}): "Y",
            frozenset({"G", "C"}): "S",
            frozenset({"A", "T"}): "W",
            frozenset({"G", "T"}): "K",
            frozenset({"A", "C"}): "M",
        }
        return lookup.get(pair, "N")