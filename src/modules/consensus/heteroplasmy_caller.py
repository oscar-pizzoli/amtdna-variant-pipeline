"""
Heteroplasmy variant caller evaluating minor base observations, strand bias,
and contamination-aware gating policies.

STATISTICAL INDEPENDENCE PRECONDITION:
All counts received in SitePileup MUST originate from physical fragment deduplication
(single-molecule resolution across mate overlaps as performed by PileupEngine).
This ensures observations in the 2x2 contingency table for Fisher's exact test
are strictly independent, preventing artificial inflation of depth and VAF.
"""

import csv
import math
import os
from typing import Iterable, List, Optional, Tuple

from scipy.stats import fisher_exact

from src.modules.consensus.consensus_builder import ConsensusBuilder
from src.modules.consensus.enums import (
    CallStatus,
    FilterReason,
    GatingAction,
)
from src.modules.consensus.models import (
    AlleleCounts,
    ConsensusPolicy,
    HeteroplasmyCall,
    ReferenceInfo,
    SitePileup,
)


class HeteroplasmyCaller:
    """
    Evaluates intra-individual minor base observations against callable coverage,
    Fisher's exact strand bias test, spatial damage geometry, and upstream contamination gating.
    """

    DAMAGE_TRANSITIONS = {("C", "T"), ("G", "A")}

    def __init__(
        self,
        policy: ConsensusPolicy,
        sample_id: str,
        reference_info: ReferenceInfo,
        gating_action: GatingAction = GatingAction.STANDARD_CALLING,
    ):
        self.policy = policy
        self.sample_id = sample_id
        self.reference_info = reference_info
        self.gating_action = gating_action

    @staticmethod
    def _get_count(counts: AlleleCounts, base: str) -> int:
        """Defensively retrieve nucleotide observation count from AlleleCounts."""
        if hasattr(counts, "get") and callable(getattr(counts, "get")):
            return counts.get(base)
        return getattr(counts, base.upper(), 0)

    @staticmethod
    def calculate_strand_bias_pvalue(
        fwd_major: int,
        rev_major: int,
        fwd_minor: int,
        rev_minor: int,
    ) -> float:
        """Compute two-sided Fisher's exact test p-value on 2x2 contingency table."""
        table = [
            [fwd_major, rev_major],
            [fwd_minor, rev_minor],
        ]
        try:
            _, pvalue = fisher_exact(table, alternative="two-sided")
            return float(pvalue) if not math.isnan(pvalue) else 1.0
        except Exception:
            return 1.0

    def evaluate_site(self, site: SitePileup) -> Optional[HeteroplasmyCall]:
        """
        Evaluate a single observational site for heteroplasmy evidence.
        Computes physical fragment support distributions for damage authentication.
        """
        if site.depth_callable == 0:
            return None

        ref_upper = site.reference_base.upper()
        ranked_bases = ConsensusBuilder.rank_canonical_bases(site.counts, ref_upper)
        top1_base, top1_count = ranked_bases[0]
        top2_base, top2_count = ranked_bases[1]

        if top2_count == 0:
            return None

        alt_reads = top2_count
        alt_vaf = alt_reads / site.depth_callable

        fwd_major = self._get_count(site.forward_counts, top1_base)
        rev_major = self._get_count(site.reverse_counts, top1_base)
        fwd_minor = self._get_count(site.forward_counts, top2_base)
        rev_minor = self._get_count(site.reverse_counts, top2_base)

        strand_bias_pvalue = self.calculate_strand_bias_pvalue(
            fwd_major=fwd_major,
            rev_major=rev_major,
            fwd_minor=fwd_minor,
            rev_minor=rev_minor,
        )

        strand_reads_pass = (
            fwd_minor >= self.policy.strand_bias_min_strand_reads
            and rev_minor >= self.policy.strand_bias_min_strand_reads
        )
        strand_pvalue_pass = strand_bias_pvalue >= self.policy.strand_bias_alpha
        strand_bias_pass = strand_reads_pass and strand_pvalue_pass

        # Evaluate Spatial Molecule Geometry (Deamination vs Authentic Polimorphism)
        alt_distances = site.allele_distances.get(top2_base, [])
        threshold = self.policy.spatial_end_distance_threshold

        internal_count = sum(1 for d in alt_distances if d > threshold)
        terminal_count = sum(1 for d in alt_distances if d <= threshold)

        is_damage_chem = (ref_upper, top2_base) in self.DAMAGE_TRANSITIONS
        if is_damage_chem and threshold > 0:
            # Transition with potential damage: suspect if internal support is under threshold
            spatial_damage_suspect = internal_count < self.policy.min_internal_fragments
        else:
            # Transversions, indels, or UDG_FULL (threshold == 0) are not suspect
            spatial_damage_suspect = False

        # Reference placeholder check (e.g. rCRS position 3107 N)
        if ref_upper not in {"A", "C", "G", "T"}:
            return HeteroplasmyCall(
                position=site.position,
                reference_base=site.reference_base,
                major_allele=top1_base,
                minor_allele=top2_base,
                callable_depth=site.depth_callable,
                alt_reads=alt_reads,
                alt_vaf=alt_vaf,
                forward_ref=fwd_major,
                reverse_ref=rev_major,
                forward_alt=fwd_minor,
                reverse_alt=rev_minor,
                mean_bq=site.mean_bq,
                mean_mq=site.mean_mq,
                strand_bias_pvalue=strand_bias_pvalue,
                status=CallStatus.FILTERED,
                filter_reason=FilterReason.FILTER_REFERENCE_MISMATCH,
                gating_action=self.gating_action,
                internal_support_count=internal_count,
                terminal_support_count=terminal_count,
                spatial_damage_suspect=spatial_damage_suspect,
            )

        # Standard Calling Profile
        std_depth_pass = site.depth_callable >= self.policy.heteroplasmy_min_depth
        std_alt_pass = alt_reads >= self.policy.heteroplasmy_min_alt_reads
        std_vaf_pass = (
            self.policy.heteroplasmy_min_vaf <= alt_vaf <= self.policy.heteroplasmy_max_vaf
        )
        std_pass = std_depth_pass and std_alt_pass and std_vaf_pass and strand_bias_pass

        def _eval_restricted_profile() -> Tuple[bool, bool, bool]:
            r_depth = site.depth_callable >= self.policy.restricted_heteroplasmy_min_depth
            r_alt = alt_reads >= self.policy.restricted_heteroplasmy_min_alt_reads
            r_vaf = (
                self.policy.restricted_heteroplasmy_min_vaf
                <= alt_vaf
                <= self.policy.restricted_heteroplasmy_max_vaf
            )
            r_pass = r_depth and r_alt and r_vaf and strand_bias_pass
            return r_pass, r_depth, r_alt

        # Gating Action Resolution
        if self.gating_action == GatingAction.STANDARD_CALLING:
            if std_pass:
                status = CallStatus.PASS
                filter_reason = FilterReason.PASS
            else:
                status = CallStatus.FILTERED
                if not std_depth_pass:
                    filter_reason = FilterReason.FILTER_LOW_DEPTH
                elif not std_alt_pass:
                    filter_reason = FilterReason.FILTER_LOW_ALT_READS
                elif not strand_bias_pass:
                    filter_reason = FilterReason.FILTER_STRAND_BIAS
                else:
                    filter_reason = FilterReason.FILTER_VAF_RANGE

        elif self.gating_action == GatingAction.RESTRICT_HETEROPLASMY:
            restr_pass, r_depth_pass, r_alt_pass = _eval_restricted_profile()
            if restr_pass:
                status = CallStatus.PASS
                filter_reason = FilterReason.PASS
            else:
                status = CallStatus.FILTERED
                if std_pass:
                    filter_reason = FilterReason.FILTER_CONTAMINATION
                elif not r_depth_pass:
                    filter_reason = FilterReason.FILTER_LOW_DEPTH
                elif not r_alt_pass:
                    filter_reason = FilterReason.FILTER_LOW_ALT_READS
                elif not strand_bias_pass:
                    filter_reason = FilterReason.FILTER_STRAND_BIAS
                else:
                    filter_reason = FilterReason.FILTER_VAF_RANGE

        elif self.gating_action == GatingAction.PMD_FILTERING_REQUIRED:
            restr_pass, r_depth_pass, r_alt_pass = _eval_restricted_profile()
            if restr_pass and not spatial_damage_suspect:
                status = CallStatus.PASS
                filter_reason = FilterReason.PASS
            else:
                status = CallStatus.FILTERED
                if spatial_damage_suspect:
                    filter_reason = FilterReason.FILTER_DAMAGE_SUSPECT
                elif std_pass:
                    filter_reason = FilterReason.FILTER_DAMAGE_SUSPECT
                elif not r_depth_pass:
                    filter_reason = FilterReason.FILTER_LOW_DEPTH
                elif not r_alt_pass:
                    filter_reason = FilterReason.FILTER_LOW_ALT_READS
                elif not strand_bias_pass:
                    filter_reason = FilterReason.FILTER_STRAND_BIAS
                else:
                    filter_reason = FilterReason.FILTER_DAMAGE_SUSPECT

        elif self.gating_action == GatingAction.USE_ENDOGENOUS_SCHMUTZI_CONSENSUS:
            restr_pass, r_depth_pass, r_alt_pass = _eval_restricted_profile()
            if restr_pass:
                status = CallStatus.PASS
                filter_reason = FilterReason.PASS
            else:
                status = CallStatus.FILTERED
                filter_reason = FilterReason.FILTER_CONTAMINATION

        elif self.gating_action in {GatingAction.REJECT_CONTAMINATED, GatingAction.EXCLUDE_OR_MANUAL_REVIEW}:
            status = CallStatus.FILTERED
            filter_reason = FilterReason.FILTER_CONTAMINATION

        elif self.gating_action == GatingAction.INSUFFICIENT_DATA:
            status = CallStatus.FILTERED
            filter_reason = FilterReason.FILTER_INSUFFICIENT_DATA

        else:
            raise ValueError(f"Unhandled GatingAction directive: {self.gating_action}")

        return HeteroplasmyCall(
            position=site.position,
            reference_base=site.reference_base,
            major_allele=top1_base,
            minor_allele=top2_base,
            callable_depth=site.depth_callable,
            alt_reads=alt_reads,
            alt_vaf=alt_vaf,
            forward_ref=fwd_major,
            reverse_ref=rev_major,
            forward_alt=fwd_minor,
            reverse_alt=rev_minor,
            mean_bq=site.mean_bq,
            mean_mq=site.mean_mq,
            strand_bias_pvalue=strand_bias_pvalue,
            status=status,
            filter_reason=filter_reason,
            gating_action=self.gating_action,
            internal_support_count=internal_count,
            terminal_support_count=terminal_count,
            spatial_damage_suspect=spatial_damage_suspect,
        )

    def call_heteroplasmies(
        self, pileup_stream: Iterable[SitePileup]
    ) -> List[HeteroplasmyCall]:
        """Stream observational sites and compile candidate heteroplasmy calls."""
        calls: List[HeteroplasmyCall] = []
        for site in pileup_stream:
            call = self.evaluate_site(site)
            if call is not None:
                calls.append(call)
        return calls

    @staticmethod
    def write_tsv(calls: List[HeteroplasmyCall], output_path: str) -> None:
        """Write all evaluated calls to TSV audit log."""
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        headers = [
            "position",
            "reference_base",
            "major_allele",
            "minor_allele",
            "callable_depth",
            "alt_reads",
            "alt_vaf",
            "forward_ref",
            "reverse_ref",
            "forward_alt",
            "reverse_alt",
            "mean_bq",
            "mean_mq",
            "strand_bias_pvalue",
            "status",
            "filter_reason",
            "gating_action",
            "internal_support_count",
            "terminal_support_count",
            "spatial_damage_suspect",
        ]

        with open(output_path, "w", newline="", encoding="utf-8") as tsvfile:
            writer = csv.writer(tsvfile, delimiter="\t")
            writer.writerow(headers)

            for call in calls:
                mean_bq_str = f"{call.mean_bq:.2f}" if call.mean_bq is not None else "NA"
                mean_mq_str = f"{call.mean_mq:.2f}" if call.mean_mq is not None else "NA"
                writer.writerow([
                    call.position,
                    call.reference_base,
                    call.major_allele,
                    call.minor_allele,
                    call.callable_depth,
                    call.alt_reads,
                    f"{call.alt_vaf:.4f}",
                    call.forward_ref,
                    call.reverse_ref,
                    call.forward_alt,
                    call.reverse_alt,
                    mean_bq_str,
                    mean_mq_str,
                    f"{call.strand_bias_pvalue:.4e}",
                    call.status.value,
                    call.filter_reason.value,
                    call.gating_action.value,
                    call.internal_support_count,
                    call.terminal_support_count,
                    call.spatial_damage_suspect,
                ])