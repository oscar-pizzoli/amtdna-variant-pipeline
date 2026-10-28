"""
Consensus sequence generator building parallel canonical, IUPAC-encoded,
and structural InDel-aware FASTA files for downstream phylogenetic haplogrouping.

ARCHITECTURAL PRINCIPLE:
Zero silent fallbacks and zero data fabrication. Uncovered loci remain strictly masked
as 'N'. Structural InDel consensus sequences resolve mitochondrial circular wraparound
across arbitrarily large structural deletions and enforce deterministic multi-allelic arbitration.
"""

from typing import Dict, Iterable, List, Optional, Set, Tuple
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
    IndelObservation,
    ReferenceInfo,
    SitePileup,
)
from src.modules.synthesis.coordinate_normalizer import MutationType


class ConsensusBuilder:
    """
    Constructs haploid canonical, IUPAC ambiguity, and circular InDel-aware consensus
    sequences from a stream of SitePileup observational records.
    """

    BASE_PRIORITY = {"A": 0, "C": 1, "G": 2, "T": 3}

    def __init__(
        self,
        policy: ConsensusPolicy,
        sample_id: str,
        reference_info: ReferenceInfo,
    ):
        self.policy = policy
        self.sample_id = sample_id
        self.reference_info = reference_info

    @classmethod
    def rank_canonical_bases(
        cls, counts: AlleleCounts, ref_base: str
    ) -> List[Tuple[str, int]]:
        """Rank canonical bases (A, C, G, T) with fully deterministic tie-breaking."""
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
        """Evaluate observational site and determine canonical and IUPAC bases."""
        ref_upper = site.reference_base.upper()

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

        ranked_bases = self.rank_canonical_bases(site.counts, ref_upper)
        top1_base, top1_count = ranked_bases[0]
        top2_base, top2_count = ranked_bases[1]

        major_vaf = top1_count / site.depth_callable
        top2_vaf = top2_count / site.depth_callable

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

        canonical_base = self.policy.mask_char

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

    def _select_dominant_indel(
        self,
        site: SitePileup,
    ) -> Optional[Tuple[MutationType, str, int]]:
        """
        Identify if a fixed/majority InDel is anchored at this site with multi-allelic arbitration.
        Arbitrates between concurrent insertions and deletions by picking the highest VAF event
        satisfying minimum depth and majority fraction thresholds.
        """
        candidate_indels: List[Tuple[float, MutationType, str, int]] = []

        for del_obs in site.deletions.values():
            if (
                del_obs.anchor_depth >= self.policy.indel_min_depth
                and del_obs.vaf >= self.policy.major_allele_fraction
            ):
                candidate_indels.append((del_obs.vaf, MutationType.DEL, del_obs.sequence, del_obs.length))

        for ins_obs in site.insertions.values():
            if (
                ins_obs.anchor_depth >= self.policy.indel_min_depth
                and ins_obs.vaf >= self.policy.major_allele_fraction
            ):
                candidate_indels.append((ins_obs.vaf, MutationType.INS, ins_obs.sequence, ins_obs.length))

        if not candidate_indels:
            return None

        # Deterministic sorting: VAF descending -> DEL priority -> Length descending
        candidate_indels.sort(
            key=lambda x: (x[0], 1 if x[1] == MutationType.DEL else 0, x[3]),
            reverse=True,
        )
        best = candidate_indels[0]
        return (best[1], best[2], best[3])

    def _precompute_circular_wraparound_deletions(
        self, pileup_list: List[SitePileup]
    ) -> Set[int]:
        """
        Scan all genomic coordinates across the contig to detect any deletion
        (regardless of size) that wraps around the circular origin (16569 -> 1).
        Returns the set of 1-based coordinates that must be deleted from the contig head.
        """
        wrap_deleted_positions: Set[int] = set()
        ref_len = self.reference_info.length

        for site in pileup_list:
            dom_indel = self._select_dominant_indel(site)
            if dom_indel is not None and dom_indel[0] == MutationType.DEL:
                del_len = dom_indel[2]
                end_pos = site.position + del_len
                if end_pos > ref_len:
                    # Deletion extends past the contig boundary: wrap around to position 1
                    overflow_count = end_pos - ref_len
                    # Collect wrapped bases starting from coordinate 1
                    for p in range(1, min(overflow_count + 1, ref_len + 1)):
                        wrap_deleted_positions.add(p)

        return wrap_deleted_positions

    def build_consensus(
        self, pileup_stream: Iterable[SitePileup]
    ) -> Tuple[List[ConsensusSite], str, str, str]:
        """
        Iterate over the pileup stream, strictly verifying coordinate continuity (1..N),
        contig dimension invariants, and generating:
        1. Canonical 1:1 FASTA string (exact reference length, masked Ns strictly preserved).
        2. IUPAC 1:1 FASTA string (exact reference length).
        3. Structural InDel-aware FASTA string (circular-aware, zero data fabrication).
        """
        pileup_list = list(pileup_stream)
        ref_len = self.reference_info.length

        # Exhaustive contig scan for wraparound deletions across the circular origin
        deleted_coordinates: Set[int] = self._precompute_circular_wraparound_deletions(pileup_list)

        consensus_sites: List[ConsensusSite] = []
        canonical_chars: List[str] = []
        iupac_chars: List[str] = []
        indel_aware_chars: List[str] = []

        expected_position = 1
        allowed_ref_bases = {"A", "C", "G", "T", "N"}

        for site in pileup_list:
            if site.position != expected_position:
                raise ValueError(
                    f"Coordinate continuity broken for sample '{self.sample_id}': "
                    f"expected position {expected_position}, observed {site.position}."
                )

            if site.position > ref_len:
                raise ValueError(
                    f"Position {site.position} exceeds expected reference length "
                    f"({ref_len}) for contig '{self.reference_info.contig}'."
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

            # Circular InDel-Aware Reconstruction
            if site.position not in deleted_coordinates:
                # Zero data fabrication: uncovered sites remain 'N'
                indel_aware_chars.append(call.canonical_base)

                # Skip structural evaluation if the locus is an unsequenced dropout
                if call.call_status != ConsensusCallStatus.NO_USABLE_BASES:
                    dom_indel = self._select_dominant_indel(site)
                    if dom_indel is not None:
                        event_type, seq, length = dom_indel
                        if event_type == MutationType.INS:
                            indel_aware_chars.append(seq.upper())
                        elif event_type == MutationType.DEL:
                            # Mark subsequent downstream coordinates using modular arithmetic
                            for step in range(1, length + 1):
                                target_coord = ((site.position - 1 + step) % ref_len) + 1
                                deleted_coordinates.add(target_coord)

            expected_position += 1

        observed_sites_count = expected_position - 1
        if observed_sites_count != ref_len:
            raise ValueError(
                f"Consensus sequence length mismatch for '{self.sample_id}': "
                f"expected {ref_len} bp, evaluated {observed_sites_count} sites."
            )

        canonical_seq = "".join(canonical_chars)
        iupac_seq = "".join(iupac_chars)
        indel_aware_seq = "".join(indel_aware_chars)

        return consensus_sites, canonical_seq, iupac_seq, indel_aware_seq

    @staticmethod
    def write_fasta(
        output_path: str,
        header_id: str,
        description: str,
        sequence: str,
        line_width: int = 60,
    ) -> None:
        """Write a nucleotide sequence to a standardized FASTA file with fixed line wrapping."""
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