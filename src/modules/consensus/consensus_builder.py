"""
Consensus sequence generator building parallel canonical and IUPAC-encoded FASTA files.

ARCHITECTURAL PRINCIPLE:
Contamination gating and minor-allele validation affect intra-individual heteroplasmy
interpretation, not majority consensus calling. ConsensusBuilder evaluates purely
observational majority signals and ambiguity boundaries, remaining agnostic to
heteroplasmy gating policies, strand-bias statistical tests, and contamination triage.
"""

from typing import Iterable, List, Tuple
import os

from src.modules.consensus.enums import (
    ConsensusCallReason,
    ConsensusCallStatus,
    IUPACAmbiguity,
)
from src.modules.consensus.models import (
    AlleleCounts,
    ConsensusPolicy,
    ConsensusSite,
    ReferenceInfo,
    SitePileup,
)


class ConsensusBuilder:
    """
    Constructs haploid canonical and IUPAC ambiguity consensus sequences
    from a stream of SitePileup observational records.
    """

    BASE_PRIORITY = {"A": 0, "C": 1, "G": 2, "T": 3}

    def __init__(
        self,
        policy: ConsensusPolicy,
        sample_id: str,
        reference_info: ReferenceInfo,
    ):
        """
        Initialize the consensus builder with policy thresholds and reference specifications.

        :param policy: ConsensusPolicy governing depth, quality, and majority thresholds.
        :param sample_id: Unique sample identifier for FASTA header labeling.
        :param reference_info: Reference genome dimensions and metadata contract.
        """
        self.policy = policy
        self.sample_id = sample_id
        self.reference_info = reference_info

    @classmethod
    def rank_canonical_bases(
        cls, counts: AlleleCounts, ref_base: str
    ) -> List[Tuple[str, int]]:
        """
        Rank canonical bases (A, C, G, T) with fully deterministic tie-breaking:
        1. Higher observation count.
        2. Concordance with reference base (preserves reference identity on exact ties).
        3. Stable alphabetical priority (A > C > G > T) as fallback.
        """
        ref_upper = ref_base.upper()
        candidates = [
            ("A", counts.A),
            ("C", counts.C),
            ("G", counts.G),
            ("T", counts.T),
        ]
        candidates.sort(
            key=lambda x: (
                x[1],
                x[0] == ref_upper,
                -cls.BASE_PRIORITY.get(x[0], 99),
            ),
            reverse=True,
        )
        return candidates

    def evaluate_site(self, site: SitePileup) -> ConsensusSite:
        """
        Evaluate a single observational site and determine canonical and IUPAC bases.
        Enforces depth floors, evaluates major allele fraction, and resolves secondary
        ambiguity without conflating it with biological heteroplasmy validation.
        """
        ref_upper = site.reference_base.upper()

        # Handle reference N placeholder (e.g. rCRS position 3107 legacy phantom deletion)
        if ref_upper == "N":
            return ConsensusSite(
                position=site.position,
                reference_base=site.reference_base,
                canonical_base=self.policy.mask_char,
                iupac_base=self.policy.mask_char,
                callable_depth=site.depth_callable,
                major_allele=self.policy.mask_char,
                major_vaf=0.0,
                call_status=ConsensusCallStatus.AMBIGUOUS,
                call_reason=ConsensusCallReason.REFERENCE_N_PLACEHOLDER,
            )

        # Baseline callable depth checks for consensus calling
        if site.depth_callable == 0:
            return ConsensusSite(
                position=site.position,
                reference_base=site.reference_base,
                canonical_base=self.policy.mask_char,
                iupac_base=self.policy.mask_char,
                callable_depth=0,
                major_allele=self.policy.mask_char,
                major_vaf=0.0,
                call_status=ConsensusCallStatus.NO_USABLE_BASES,
                call_reason=ConsensusCallReason.NO_USABLE_BASES,
            )

        if site.depth_callable < self.policy.min_consensus_depth:
            return ConsensusSite(
                position=site.position,
                reference_base=site.reference_base,
                canonical_base=self.policy.mask_char,
                iupac_base=self.policy.mask_char,
                callable_depth=site.depth_callable,
                major_allele=self.policy.mask_char,
                major_vaf=0.0,
                call_status=ConsensusCallStatus.LOW_DEPTH,
                call_reason=ConsensusCallReason.LOW_DEPTH,
            )

        # Rank canonical bases using unified deterministic tie-breaking
        ranked_bases = self.rank_canonical_bases(site.counts, ref_upper)
        top1_base, top1_count = ranked_bases[0]
        top2_base, top2_count = ranked_bases[1]

        major_vaf = top1_count / site.depth_callable
        top2_vaf = top2_count / site.depth_callable

        # Clear majority call: unambiguous consensus
        if major_vaf >= self.policy.major_allele_fraction:
            return ConsensusSite(
                position=site.position,
                reference_base=site.reference_base,
                canonical_base=top1_base,
                iupac_base=top1_base,
                callable_depth=site.depth_callable,
                major_allele=top1_base,
                major_vaf=major_vaf,
                call_status=ConsensusCallStatus.CALLED,
                call_reason=ConsensusCallReason.CALLED_MAJORITY,
            )

        # Contested majority: canonical stream must be masked
        canonical_base = self.policy.mask_char

        # Secondary allele above ambiguity floor warrants IUPAC degenerate symbol
        if top2_vaf >= self.policy.ambiguity_min_vaf:
            iupac_resolved = (
                IUPACAmbiguity.resolve(top1_base, top2_base)
                if self.policy.emit_iupac_fasta
                else self.policy.mask_char
            )
            return ConsensusSite(
                position=site.position,
                reference_base=site.reference_base,
                canonical_base=canonical_base,
                iupac_base=iupac_resolved,
                callable_depth=site.depth_callable,
                major_allele=top1_base,
                major_vaf=major_vaf,
                call_status=ConsensusCallStatus.AMBIGUOUS,
                call_reason=ConsensusCallReason.AMBIGUOUS_IUPAC,
            )

        # Contested majority but secondary allele below ambiguity floor
        return ConsensusSite(
            position=site.position,
            reference_base=site.reference_base,
            canonical_base=canonical_base,
            iupac_base=self.policy.mask_char,
            callable_depth=site.depth_callable,
            major_allele=top1_base,
            major_vaf=major_vaf,
            call_status=ConsensusCallStatus.AMBIGUOUS,
            call_reason=ConsensusCallReason.SECOND_ALLELE_TOO_LOW,
        )

    def build_consensus(
        self, pileup_stream: Iterable[SitePileup]
    ) -> Tuple[List[ConsensusSite], str, str]:
        """
        Iterate over the pileup stream, strictly verifying coordinate continuity (1..N),
        contig dimension invariants, and reference base validity (allowing rCRS legacy N).

        :param pileup_stream: Sequential stream of SitePileup objects.
        :return: Tuple of (List of ConsensusSite, canonical_sequence, iupac_sequence).
        """
        consensus_sites: List[ConsensusSite] = []
        canonical_chars: List[str] = []
        iupac_chars: List[str] = []

        expected_position = 1
        allowed_ref_bases = {"A", "C", "G", "T", "N"}

        for site in pileup_stream:
            # Enforce 1-based strictly sequential gapless coordinates
            if site.position != expected_position:
                raise ValueError(
                    f"Coordinate continuity broken for sample '{self.sample_id}': "
                    f"expected position {expected_position}, observed {site.position}."
                )

            if site.position > self.reference_info.length:
                raise ValueError(
                    f"Position {site.position} exceeds expected reference length "
                    f"({self.reference_info.length}) for contig '{self.reference_info.contig}'."
                )

            ref_upper = site.reference_base.upper()
            if ref_upper not in allowed_ref_bases:
                raise ValueError(
                    f"Invalid non-canonical reference base '{site.reference_base}' "
                    f"at position {site.position} for sample '{self.sample_id}'."
                )

            call = self.evaluate_site(site)
            consensus_sites.append(call)
            canonical_chars.append(call.canonical_base)
            iupac_chars.append(call.iupac_base)

            expected_position += 1

        # Enforce exact total length invariant
        observed_sites_count = expected_position - 1
        if observed_sites_count != self.reference_info.length:
            raise ValueError(
                f"Consensus sequence length mismatch for '{self.sample_id}': "
                f"expected {self.reference_info.length} bp, evaluated {observed_sites_count} sites."
            )

        canonical_seq = "".join(canonical_chars)
        iupac_seq = "".join(iupac_chars)

        return consensus_sites, canonical_seq, iupac_seq

    @staticmethod
    def write_fasta(
        output_path: str,
        header_id: str,
        description: str,
        sequence: str,
        line_width: int = 60,
    ) -> None:
        """
        Write a nucleotide sequence to a standardized FASTA file with fixed line wrapping.
        """
        if line_width <= 0:
            raise ValueError(f"FASTA line_width must be positive, got {line_width}")

        clean_header = header_id.strip()
        clean_desc = description.strip()

        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as out:
            header_line = f">{clean_header} {clean_desc}".strip()
            out.write(f"{header_line}\n")
            for i in range(0, len(sequence), line_width):
                out.write(sequence[i : i + line_width] + "\n")
