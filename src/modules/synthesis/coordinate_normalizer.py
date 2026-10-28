"""
Block 8 Coordinate Normalizer: Unified positioning utility resolving coordinate
divergences between VCF records, JSON string keys, and HaploGrep notation.
Guarantees deterministic relational joins across heterogeneous bioinformatics dialects.
"""

from dataclasses import dataclass
from enum import Enum
import logging
import re
from typing import Any, List, Optional, Tuple

logger = logging.getLogger("mtDNA_Pipeline.CoordinateNormalizer")

DEFAULT_MT_LENGTH = 16569


class MutationType(str, Enum):
    """Structural classification of a sequence variant."""
    SNV = "SNV"
    INS = "INS"
    DEL = "DEL"
    COMPLEX = "COMPLEX"


@dataclass(frozen=True)
class CanonicalMutation:
    """Standardized representation of an atomic mitochondrial variant across all pipeline dialects."""
    raw_token: str
    position: Optional[int]
    end_position: Optional[int]          # Populated for deletions and multi-base spans
    insertion_index: int                 # 0 for genomic coordinates; 1, 2, ... for insertions (e.g. 315.1C)
    ref: str
    alt: str
    variant_type: MutationType
    notation: str
    contains_n: bool = False             # True if REF or ALT contains 'N'
    is_ambiguous_iupac: bool = False     # True if variant carries IUPAC ambiguity codes
    is_unparseable: bool = False         # True when parsing syntax fails
    parse_error: Optional[str] = None
    tentative_position: Optional[int] = None

    @property
    def canonical_key(self) -> str:
        """
        Source-independent join key matching VCF, Consensus JSON, and HaploGrep identically.
        Guarantees SNV bi-allelic discrimination, semantic deletion matching, and insertion index parity.
        """
        if self.is_unparseable or self.position is None:
            return f"UNPARSED_{self.raw_token}"
        if self.variant_type == MutationType.INS:
            return f"{self.position}.{self.insertion_index}_INS_{self.alt}"
        if self.variant_type == MutationType.DEL:
            if self.end_position and self.end_position != self.position:
                return f"{self.position}-{self.end_position}_DEL"
            return f"{self.position}_DEL"
        # SNV and COMPLEX include REF to preserve multi-allelic and reference-specific resolution
        return f"{self.position}_{self.ref}>{self.alt}"

    def deletion_length(self, ref_length: int = DEFAULT_MT_LENGTH) -> int:
        """Calculate deleted span with strict handling of circular origin wraparound."""
        if self.variant_type != MutationType.DEL or self.position is None:
            return 0
        if self.end_position is not None:
            if self.end_position >= self.position:
                return (self.end_position - self.position) + 1
            # Circular genome wraparound: e.g. 16569 to 1
            return (ref_length - self.position + 1) + self.end_position
        return len(self.ref) if self.ref not in {".", "-", "N"} else 1

    @property
    def is_circular_wraparound(self) -> bool:
        """Flag events spanning across the circular origin boundary."""
        return bool(
            self.variant_type == MutationType.DEL
            and self.position is not None
            and self.end_position is not None
            and self.position > self.end_position
        )


class MutationList(list):
    """
    Polymorphic container inheriting from list.
    Supports iteration over multi-variant unpackings while allowing direct attribute
    access to the primary variant for full backwards-compatibility.
    """
    def __getattr__(self, name: str) -> Any:
        if self:
            return getattr(self[0], name)
        raise AttributeError(f"'MutationList' object has no attribute '{name}'")


class CoordinateNormalizer:
    """Robust parser, prefix/suffix trimmer, and canonicalizer for mitochondrial variants."""

    CLEAN_EXTRAS_REGEX = re.compile(r"\s*\([^)]*\)")
    IUPAC_AMBIGUOUS = set("RYSWKMBDHV")

    # 1. Deletions: Prioritizing ranges, explicit keywords, and anchored single-base markers
    DEL_RANGE_REGEX = re.compile(
        r"^(?:([ACGTN])>?)?([0-9]{1,5})-([0-9]{1,5})(?:D|DEL)?$", re.IGNORECASE
    )
    DEL_KEYWORD_REGEX = re.compile(
        r"^(?:([ACGTN])>?)?([0-9]{1,5})(?:DEL)$", re.IGNORECASE
    )
    DEL_SINGLE_REGEX = re.compile(
        r"^(?:([ACGTN])>?)?([0-9]{1,5})d$"  # Explicit lowercase 'd' from HaploGrep
    )
    DEL_BARE_D_REGEX = re.compile(
        r"^([0-9]{1,5})D$"  # HaploGrep bare D notation without REF base prefix (e.g. 16519D)
    )

    # 2. Insertions: Explicit dot notation (e.g. "315.1C", "309.2CC")
    INS_REGEX = re.compile(
        r"^([0-9]{1,5})\.([0-9]{1,2})([ACGTNRYSWKMBDHV]+)$", re.IGNORECASE
    )

    # 3. Substitutions (SNVs & IUPAC Ambiguities): e.g. "73G", "A73G", "73R", "A73R", "A>73>G", "A73D"
    SNV_REGEX = re.compile(
        r"^(?:([ACGTN])>?)?([0-9]{1,5})>?(?:([ACGTNRYSWKMBDHV]))$", re.IGNORECASE
    )

    def __init__(
        self,
        reference_length: int = DEFAULT_MT_LENGTH,
        reference_sequence: Optional[str] = None,
    ) -> None:
        self.reference_length = reference_length
        self.reference_sequence = reference_sequence

    @classmethod
    def from_reference_info(
        cls,
        reference_info: Any,
        reference_sequence: Optional[str] = None,
    ) -> "CoordinateNormalizer":
        """Factory initializing normalizer with metadata from ReferenceInfo."""
        ref_len = getattr(reference_info, "length_bp", getattr(reference_info, "length", DEFAULT_MT_LENGTH))
        return cls(reference_length=ref_len, reference_sequence=reference_sequence)

    @classmethod
    def clean_raw_token(cls, raw: Any) -> str:
        """Strip parenthetical annotations, trailing exclamation marks, and whitespace."""
        if raw is None:
            return ""
        s = str(raw).strip()
        s = cls.CLEAN_EXTRAS_REGEX.sub("", s)
        return s.rstrip("!").strip()

    def to_canonical_position(self, pos_val: Any) -> Optional[int]:
        """
        Defensively resolve only integer coordinates bounded in [1, reference_length].
        Rejects decimal coordinates (e.g. 315.1) to protect insertion semantics.
        """
        if pos_val is None:
            return None
        if isinstance(pos_val, float):
            if not pos_val.is_integer():
                return None
            val = int(pos_val)
        elif isinstance(pos_val, int):
            val = pos_val
        else:
            s = str(pos_val).strip()
            if not s.isdigit():
                return None
            try:
                val = int(s)
            except (ValueError, TypeError):
                return None

        return val if 1 <= val <= self.reference_length else None

    def _extract_reference_span(
        self,
        start_pos: int,
        end_pos: int,
        ref_seq_override: Optional[str] = None,
    ) -> str:
        """Extract reference bases, correctly concatenating across the circular origin."""
        seq = ref_seq_override or self.reference_sequence
        if not seq:
            return "N"
        seq_len = len(seq)
        if not (1 <= start_pos <= seq_len):
            return "N"

        if start_pos <= end_pos:
            if end_pos <= seq_len:
                return seq[start_pos - 1 : end_pos].upper()
            return "N"
        else:
            # Circular boundary crossover: start_pos > end_pos (e.g. 16569 -> 1)
            part1 = seq[start_pos - 1 : seq_len]
            part2 = seq[0 : min(end_pos, seq_len)]
            return (part1 + part2).upper()

    def _normalize_vcf_record(
        self,
        pos: int,
        ref: str,
        alt: str,
        ref_seq_override: Optional[str] = None,
    ) -> List[CanonicalMutation]:
        """
        Full VCF representation reducer applying prefix/suffix trimming,
        modular origin wraparound on curr_pos, atomic tandem insertion expansion,
        and MNV-to-SNV decomposition.
        """
        ref_u = ref.upper()
        alt_u = alt.upper()

        # Case A: Trivial single-base substitution
        if len(ref_u) == 1 and len(alt_u) == 1 and alt_u not in {".", "-", "DEL"}:
            notation = f"{ref_u}{pos}{alt_u}"
            has_n = (ref_u == "N") or (alt_u == "N")
            is_iupac = bool(set(alt_u).intersection(self.IUPAC_AMBIGUOUS))
            return [
                CanonicalMutation(
                    raw_token=f"{pos}:{ref_u}>{alt_u}",
                    position=pos,
                    end_position=None,
                    insertion_index=0,
                    ref=ref_u,
                    alt=alt_u,
                    variant_type=MutationType.SNV,
                    notation=notation,
                    contains_n=has_n,
                    is_ambiguous_iupac=is_iupac,
                    is_unparseable=False,
                    tentative_position=pos,
                )
            ]

        # Case B: Pure deletion without anchor
        if alt_u in {".", "-", "", "DEL"}:
            del_len = len(ref_u)
            end_pos = pos + del_len - 1
            if end_pos > self.reference_length:
                end_pos = (end_pos - 1) % self.reference_length + 1
            notation = f"{pos}-{end_pos}d" if end_pos != pos else f"{pos}d"
            return [
                CanonicalMutation(
                    raw_token=f"{pos}:{ref_u}>{alt_u}",
                    position=pos,
                    end_position=end_pos,
                    insertion_index=0,
                    ref=ref_u,
                    alt="DEL",
                    variant_type=MutationType.DEL,
                    notation=notation,
                    contains_n=("N" in ref_u),
                    is_ambiguous_iupac=False,
                    is_unparseable=False,
                    tentative_position=pos,
                )
            ]

        # Case C: Pure insertion without anchor -> Unpack into atomic mutations
        if ref_u in {".", "-", ""}:
            atomic_ins: List[CanonicalMutation] = []
            for idx, base in enumerate(alt_u, start=1):
                notation = f"{pos}.{idx}{base}"
                atomic_ins.append(
                    CanonicalMutation(
                        raw_token=f"{pos}:{ref_u}>{alt_u}[{idx}]",
                        position=pos,
                        end_position=None,
                        insertion_index=idx,
                        ref=".",
                        alt=base,
                        variant_type=MutationType.INS,
                        notation=notation,
                        contains_n=(base == "N"),
                        is_ambiguous_iupac=(base in self.IUPAC_AMBIGUOUS),
                        is_unparseable=False,
                        tentative_position=pos,
                    )
                )
            return atomic_ins

        # Case D: Anchored VCF InDel / MNV -> Two-way prefix/suffix trimming
        common_prefix = 0
        min_len = min(len(ref_u), len(alt_u))
        while common_prefix < min_len and ref_u[common_prefix] == alt_u[common_prefix]:
            common_prefix += 1

        # Strict circular modular wraparound on curr_pos
        curr_pos = pos + common_prefix
        if curr_pos > self.reference_length:
            curr_pos = (curr_pos - 1) % self.reference_length + 1

        trimmed_ref = ref_u[common_prefix:]
        trimmed_alt = alt_u[common_prefix:]

        common_suffix = 0
        min_rem_len = min(len(trimmed_ref), len(trimmed_alt))
        while (
            common_suffix < min_rem_len
            and trimmed_ref[-(common_suffix + 1)] == trimmed_alt[-(common_suffix + 1)]
        ):
            common_suffix += 1

        if common_suffix > 0:
            trimmed_ref = trimmed_ref[:-common_suffix]
            trimmed_alt = trimmed_alt[:-common_suffix]

        # Evaluation post-trimming:
        # D1. Single-base SNV (e.g. POS=72, REF=AC, ALT=AG -> C>G at 73)
        if len(trimmed_ref) == 1 and len(trimmed_alt) == 1:
            notation = f"{trimmed_ref}{curr_pos}{trimmed_alt}"
            has_n = (trimmed_ref == "N") or (trimmed_alt == "N")
            is_iupac = bool(set(trimmed_alt).intersection(self.IUPAC_AMBIGUOUS))
            return [
                CanonicalMutation(
                    raw_token=f"{pos}:{ref_u}>{alt_u}",
                    position=curr_pos,
                    end_position=None,
                    insertion_index=0,
                    ref=trimmed_ref,
                    alt=trimmed_alt,
                    variant_type=MutationType.SNV,
                    notation=notation,
                    contains_n=has_n,
                    is_ambiguous_iupac=is_iupac,
                    is_unparseable=False,
                    tentative_position=curr_pos,
                )
            ]

        # D2. Tandem Insertion (e.g. POS=309, REF=C, ALT=CCT -> 309.1C, 309.2T)
        if not trimmed_ref and trimmed_alt:
            ins_anchor = curr_pos - 1
            if ins_anchor == 0:
                ins_anchor = self.reference_length
            atomic_ins = []
            for idx, base in enumerate(trimmed_alt, start=1):
                notation = f"{ins_anchor}.{idx}{base}"
                atomic_ins.append(
                    CanonicalMutation(
                        raw_token=f"{pos}:{ref_u}>{alt_u}[{idx}]",
                        position=ins_anchor,
                        end_position=None,
                        insertion_index=idx,
                        ref=".",
                        alt=base,
                        variant_type=MutationType.INS,
                        notation=notation,
                        contains_n=(base == "N"),
                        is_ambiguous_iupac=(base in self.IUPAC_AMBIGUOUS),
                        is_unparseable=False,
                        tentative_position=ins_anchor,
                    )
                )
            return atomic_ins

        # D3. Multi-base Deletion (e.g. POS=522, REF=CAC, ALT=C -> 523-524d)
        if trimmed_ref and not trimmed_alt:
            del_start = curr_pos
            del_len = len(trimmed_ref)
            del_end = del_start + del_len - 1
            if del_start > self.reference_length:
                del_start = (del_start - 1) % self.reference_length + 1
            if del_end > self.reference_length:
                del_end = (del_end - 1) % self.reference_length + 1

            final_ref = self._extract_reference_span(del_start, del_end, ref_seq_override)
            if final_ref == "N":
                final_ref = trimmed_ref

            notation = f"{del_start}-{del_end}d" if del_end != del_start else f"{del_start}d"
            return [
                CanonicalMutation(
                    raw_token=f"{pos}:{ref_u}>{alt_u}",
                    position=del_start,
                    end_position=del_end,
                    insertion_index=0,
                    ref=final_ref,
                    alt="DEL",
                    variant_type=MutationType.DEL,
                    notation=notation,
                    contains_n=("N" in final_ref),
                    is_ambiguous_iupac=False,
                    is_unparseable=False,
                    tentative_position=del_start,
                )
            ]

        # D4. Block Substitution (MNV Decomposition into independent SNVs)
        if len(trimmed_ref) == len(trimmed_alt) and len(trimmed_ref) > 1:
            atomic_snvs = []
            for k in range(len(trimmed_ref)):
                r_base = trimmed_ref[k]
                a_base = trimmed_alt[k]
                if r_base != a_base:
                    site_pos = curr_pos + k
                    if site_pos > self.reference_length:
                        site_pos = (site_pos - 1) % self.reference_length + 1
                    notation = f"{r_base}{site_pos}{a_base}"
                    has_n = (r_base == "N") or (a_base == "N")
                    is_iupac = bool(set(a_base).intersection(self.IUPAC_AMBIGUOUS))
                    atomic_snvs.append(
                        CanonicalMutation(
                            raw_token=f"{pos}:{ref_u}>{alt_u}[{site_pos}]",
                            position=site_pos,
                            end_position=None,
                            insertion_index=0,
                            ref=r_base,
                            alt=a_base,
                            variant_type=MutationType.SNV,
                            notation=notation,
                            contains_n=has_n,
                            is_ambiguous_iupac=is_iupac,
                            is_unparseable=False,
                            tentative_position=site_pos,
                        )
                    )
            if atomic_snvs:
                return atomic_snvs

        # D5. Complex Non-decomposable InDel
        end_pos = curr_pos + len(trimmed_ref) - 1
        if end_pos > self.reference_length:
            end_pos = (end_pos - 1) % self.reference_length + 1
        notation = f"{curr_pos}_{trimmed_ref}>{trimmed_alt}"
        return [
            CanonicalMutation(
                raw_token=f"{pos}:{ref_u}>{alt_u}",
                position=curr_pos,
                end_position=end_pos,
                insertion_index=0,
                ref=trimmed_ref,
                alt=trimmed_alt,
                variant_type=MutationType.COMPLEX,
                notation=notation,
                contains_n=("N" in trimmed_ref or "N" in trimmed_alt),
                is_ambiguous_iupac=False,
                is_unparseable=False,
                tentative_position=curr_pos,
            )
        ]

    def parse_mutation(
        self,
        raw_mutation: Any,
        reference_sequence: Optional[str] = None,
        default_pos: Optional[int] = None,
        default_ref: Optional[str] = None,
        default_alt: Optional[str] = None,
    ) -> MutationList:
        """
        Normalize any mutation representation into a MutationList of CanonicalMutation objects.
        Dispatches deterministically across VCF calls, InDel range tokens, and single-base markers.
        """
        token = self.clean_raw_token(raw_mutation)
        ref_seq = reference_sequence or self.reference_sequence

        # ----------------------------------------------------------------------
        # Branch 1: Ingestion via explicit VCF / Caller fields
        # ----------------------------------------------------------------------
        if not token and default_pos is not None:
            pos = self.to_canonical_position(default_pos)
            if pos is None:
                return MutationList([
                    CanonicalMutation(
                        raw_token=f"{default_pos}:{default_ref}>{default_alt}",
                        position=None,
                        end_position=None,
                        insertion_index=0,
                        ref=str(default_ref or "N").upper(),
                        alt=str(default_alt or "N").upper(),
                        variant_type=MutationType.COMPLEX,
                        notation=f"INVALID_POS_{default_pos}",
                        contains_n=False,
                        is_ambiguous_iupac=False,
                        is_unparseable=True,
                        parse_error=f"Position {default_pos} is outside reference range [1, {self.reference_length}].",
                        tentative_position=None,
                    )
                ])

            ref_str = str(default_ref or "N").upper()
            alt_str = str(default_alt or "N").upper()

            return MutationList(self._normalize_vcf_record(pos, ref_str, alt_str, ref_seq))

        if not token:
            return MutationList([
                CanonicalMutation(
                    raw_token="",
                    position=None,
                    end_position=None,
                    insertion_index=0,
                    ref="N",
                    alt="N",
                    variant_type=MutationType.COMPLEX,
                    notation="EMPTY_TOKEN",
                    contains_n=False,
                    is_ambiguous_iupac=False,
                    is_unparseable=True,
                    parse_error="Empty or null mutation token provided.",
                    tentative_position=None,
                )
            ])

        # ----------------------------------------------------------------------
        # Branch 2: Deletion Notations (Ranges, DEL keyword, lowercase 'd', bare 'D')
        # ----------------------------------------------------------------------
        # 2a. Range Deletion (e.g. 523-524d, 523-524DEL, 16569-1d)
        del_range = self.DEL_RANGE_REGEX.match(token)
        if del_range:
            start_pos = self.to_canonical_position(del_range.group(2))
            end_pos = self.to_canonical_position(del_range.group(3))
            if start_pos and end_pos:
                ref_b = self._extract_reference_span(start_pos, end_pos, ref_seq)
                notation = f"{start_pos}-{end_pos}d"
                return MutationList([
                    CanonicalMutation(
                        raw_token=token,
                        position=start_pos,
                        end_position=end_pos,
                        insertion_index=0,
                        ref=ref_b,
                        alt="DEL",
                        variant_type=MutationType.DEL,
                        notation=notation,
                        contains_n=("N" in ref_b),
                        is_ambiguous_iupac=False,
                        is_unparseable=False,
                        tentative_position=start_pos,
                    )
                ])

        # 2b. Explicit DEL Keyword, lowercase 'd', or bare 'D' (without REF base prefix)
        del_match = (
            self.DEL_KEYWORD_REGEX.match(token)
            or self.DEL_SINGLE_REGEX.match(token)
            or self.DEL_BARE_D_REGEX.match(token)
        )
        if del_match:
            pos_str = del_match.group(2) if del_match.lastindex and del_match.lastindex >= 2 else del_match.group(1)
            pos = self.to_canonical_position(pos_str)
            if pos is not None:
                ref_b = self._extract_reference_span(pos, pos, ref_seq)
                return MutationList([
                    CanonicalMutation(
                        raw_token=token,
                        position=pos,
                        end_position=pos,
                        insertion_index=0,
                        ref=ref_b,
                        alt="DEL",
                        variant_type=MutationType.DEL,
                        notation=f"{pos}d",
                        contains_n=("N" in ref_b),
                        is_ambiguous_iupac=False,
                        is_unparseable=False,
                        tentative_position=pos,
                    )
                ])

        # ----------------------------------------------------------------------
        # Branch 3: Insertion Notation (e.g. 315.1C, 309.2CC)
        # ----------------------------------------------------------------------
        ins_match = self.INS_REGEX.match(token)
        if ins_match:
            pos = self.to_canonical_position(ins_match.group(1))
            if pos is not None:
                base_ins_idx = int(ins_match.group(2))
                alt_bases = ins_match.group(3).upper()

                # Unpack multi-base insertion token if present (e.g. 309.2CC -> 309.2C, 309.3C)
                atomic_ins = []
                for offset, base in enumerate(alt_bases):
                    current_idx = base_ins_idx + offset
                    is_iupac = bool(set(base).intersection(self.IUPAC_AMBIGUOUS))
                    atomic_ins.append(
                        CanonicalMutation(
                            raw_token=token,
                            position=pos,
                            end_position=None,
                            insertion_index=current_idx,
                            ref=".",
                            alt=base,
                            variant_type=MutationType.INS,
                            notation=f"{pos}.{current_idx}{base}",
                            contains_n=(base == "N"),
                            is_ambiguous_iupac=is_iupac,
                            is_unparseable=False,
                            tentative_position=pos,
                        )
                    )
                return MutationList(atomic_ins)

        # ----------------------------------------------------------------------
        # Branch 4: Substitutions (SNVs & Full IUPAC Alphabet, including A73D)
        # ----------------------------------------------------------------------
        snv_match = self.SNV_REGEX.match(token)
        if snv_match:
            pos = self.to_canonical_position(snv_match.group(2))
            if pos is not None:
                raw_ref = snv_match.group(1)
                alt_b = snv_match.group(3).upper()
                ref_b = raw_ref.upper() if raw_ref else self._extract_reference_span(pos, pos, ref_seq)

                notation = f"{ref_b}{pos}{alt_b}" if ref_b != "N" else f"{pos}{alt_b}"
                has_n = (alt_b == "N") or (ref_b == "N")
                is_iupac = (alt_b in self.IUPAC_AMBIGUOUS)

                return MutationList([
                    CanonicalMutation(
                        raw_token=token,
                        position=pos,
                        end_position=None,
                        insertion_index=0,
                        ref=ref_b,
                        alt=alt_b,
                        variant_type=MutationType.SNV,
                        notation=notation,
                        contains_n=has_n,
                        is_ambiguous_iupac=is_iupac,
                        is_unparseable=False,
                        tentative_position=pos,
                    )
                ])

        # ----------------------------------------------------------------------
        # Branch 5: Unparseable Syntax Fallback
        # ----------------------------------------------------------------------
        digits = re.findall(r"\d+", token)
        tentative = int(digits[0]) if digits and digits[0].isdigit() else None

        return MutationList([
            CanonicalMutation(
                raw_token=token,
                position=None,
                end_position=None,
                insertion_index=0,
                ref="N",
                alt="N",
                variant_type=MutationType.COMPLEX,
                notation=token,
                contains_n=False,
                is_ambiguous_iupac=False,
                is_unparseable=True,
                parse_error=f"Unrecognized mitochondrial mutation syntax: '{token}'",
                tentative_position=tentative,
            )
        ])


# Module-level default instance for backwards-compatible direct invocation
_DEFAULT_NORMALIZER = CoordinateNormalizer()


def parse_mutation(
    raw_mutation: Any,
    reference_sequence: Optional[str] = None,
    default_pos: Optional[int] = None,
    default_ref: Optional[str] = None,
    default_alt: Optional[str] = None,
) -> MutationList:
    """Convenience module function delegating to the default rCRS normalizer instance."""
    return _DEFAULT_NORMALIZER.parse_mutation(
        raw_mutation=raw_mutation,
        reference_sequence=reference_sequence,
        default_pos=default_pos,
        default_ref=default_ref,
        default_alt=default_alt,
    )