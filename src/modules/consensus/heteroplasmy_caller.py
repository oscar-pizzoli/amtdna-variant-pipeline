"""
Heteroplasmy variant caller evaluating minor base observations, structural InDels,
strand bias, spatial damage geometry, homopolymer slippage, and excluded interval blacklist.
Enforces strict contamination gating and spatial fragment authentication parity across SNVs and InDels.
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
    IndelObservation,
    ReferenceInfo,
    SitePileup,
)
from src.modules.synthesis.coordinate_normalizer import MutationType


class HeteroplasmyCaller:
    """
    Evaluates intra-individual minor base observations and structural InDels against
    callable coverage, Fisher's exact strand bias test, spatial damage geometry,
    homopolymer slippage tracts, excluded interval blacklists, and upstream contamination gating.
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

    def _evaluate_snv_heteroplasmy(
        self,
        site: SitePileup,
        dominant_indel_vaf: float = 0.0,
    ) -> Optional[HeteroplasmyCall]:
        """
        Evaluate single nucleotide minor allele observations.
        Suppresses false SNVs if the locus is structurally dominated by an InDel event
        or falls within an excluded artifact/slippage interval.
        """
        if site.depth_callable == 0:
            return None

        effective_depth = site.depth_callable

        ref_upper = site.reference_base.upper()
        ranked_bases = ConsensusBuilder.rank_canonical_bases(site.counts, ref_upper)
        top1_base, top1_count = ranked_bases[0]
        top2_base, top2_count = ranked_bases[1]

        if top2_count == 0:
            return None

        alt_reads = top2_count
        total_allelic_molecules = max(site.depth_callable, site.indel_anchor_depth)
        alt_vaf = alt_reads / total_allelic_molecules if total_allelic_molecules > 0 else 0.0

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

        # Valutazione della geometria molecolare (supporto interno vs terminale)
        alt_distances = site.read_end_distances.get(top2_base, [])
        threshold = self.policy.spatial_end_distance_threshold

        internal_count = sum(1 for d in alt_distances if d > threshold)
        terminal_count = sum(1 for d in alt_distances if d <= threshold)

        is_damage_chem = (ref_upper, top2_base) in self.DAMAGE_TRANSITIONS
        if is_damage_chem and threshold > 0:
            spatial_damage_suspect = internal_count < self.policy.min_internal_fragments
        else:
            spatial_damage_suspect = False

        # Verifica preliminare base di riferimento non canonica
        if ref_upper not in {"A", "C", "G", "T"}:
            filt_reason = (
                FilterReason.FILTER_EXCLUDED_INTERVAL
                if self.policy.is_position_excluded(site.position)
                else FilterReason.FILTER_REFERENCE_MISMATCH
            )
            return HeteroplasmyCall(
                position=site.position,
                reference_base=site.reference_base,
                major_allele=top1_base,
                minor_allele=top2_base,
                callable_depth=effective_depth,
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
                filter_reason=filt_reason,
                gating_action=self.gating_action,
                variant_type=MutationType.SNV,
                internal_support_count=internal_count,
                terminal_support_count=terminal_count,
                spatial_damage_suspect=spatial_damage_suspect,
            )

        # Soglie per profilo standard
        std_depth_pass = effective_depth >= self.policy.heteroplasmy_min_depth
        std_alt_pass = alt_reads >= self.policy.heteroplasmy_min_alt_reads
        std_vaf_pass = (
            self.policy.heteroplasmy_min_vaf <= alt_vaf <= self.policy.heteroplasmy_max_vaf
        )
        std_pass = std_depth_pass and std_alt_pass and std_vaf_pass and strand_bias_pass

        def _eval_restricted_profile() -> Tuple[bool, bool, bool]:
            r_depth = effective_depth >= self.policy.restricted_heteroplasmy_min_depth
            r_alt = alt_reads >= self.policy.restricted_heteroplasmy_min_alt_reads
            r_vaf = (
                self.policy.restricted_heteroplasmy_min_vaf
                <= alt_vaf
                <= self.policy.restricted_heteroplasmy_max_vaf
            )
            r_pass = r_depth and r_alt and r_vaf and strand_bias_pass
            return r_pass, r_depth, r_alt

        # Routing in base all'azione di gating
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
                if spatial_damage_suspect or std_pass:
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
            restr_pass, _, _ = _eval_restricted_profile()
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

        # Triage intervalli esclusi / blacklist per SNV
        if self.policy.is_position_excluded(site.position):
            status = CallStatus.FILTERED
            filter_reason = FilterReason.FILTER_EXCLUDED_INTERVAL

        return HeteroplasmyCall(
            position=site.position,
            reference_base=site.reference_base,
            major_allele=top1_base,
            minor_allele=top2_base,
            callable_depth=effective_depth,
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
            variant_type=MutationType.SNV,
            internal_support_count=internal_count,
            terminal_support_count=terminal_count,
            spatial_damage_suspect=spatial_damage_suspect,
        )

    def _evaluate_indel_observation(
        self,
        site: SitePileup,
        obs: IndelObservation,
    ) -> Optional[HeteroplasmyCall]:
        """
        Evaluate structural InDel observations with rigorous parity to SNVs:
        Applies physical terminal damage authentication, homopolymer slippage,
        and excluded interval blacklist filtering.
        """
        if not self.policy.enable_indel_calling or obs.anchor_depth == 0:
            return None

        alt_reads = obs.total_count
        if alt_reads == 0:
            return None

        alt_vaf = obs.vaf
        ref_base = site.reference_base.upper()

        if obs.event_type == MutationType.INS:
            major_allele = ref_base
            minor_allele = f"{ref_base}{obs.sequence}"
        else:
            major_allele = f"{ref_base}{obs.sequence}"
            minor_allele = ref_base

        fwd_indel = obs.fwd_count
        rev_indel = obs.rev_count
        fwd_anchor = site.forward_counts.canonical_total
        rev_anchor = site.reverse_counts.canonical_total

        fwd_ref_anchor = max(0, fwd_anchor - fwd_indel)
        rev_ref_anchor = max(0, rev_anchor - rev_indel)

        strand_bias_pvalue = self.calculate_strand_bias_pvalue(
            fwd_major=fwd_ref_anchor,
            rev_major=rev_ref_anchor,
            fwd_minor=fwd_indel,
            rev_minor=rev_indel,
        )

        strand_reads_pass = (
            fwd_indel >= self.policy.strand_bias_min_strand_reads
            and rev_indel >= self.policy.strand_bias_min_strand_reads
        )
        strand_pvalue_pass = strand_bias_pvalue >= self.policy.strand_bias_alpha
        strand_bias_pass = strand_reads_pass and strand_pvalue_pass

        threshold = self.policy.spatial_end_distance_threshold
        internal_count = sum(1 for d in obs.read_end_distances if d > threshold)
        terminal_count = sum(1 for d in obs.read_end_distances if d <= threshold)

        if threshold > 0:
            spatial_damage_suspect = internal_count < self.policy.min_internal_fragments
        else:
            spatial_damage_suspect = False

        # Rilevamento tratti di slippage omopolimerico
        slip_region = self.policy.get_slippage_region(site.position)
        is_slippage = slip_region is not None
        slip_name = slip_region.name if slip_region else None

        # Rilevamento intervalli esclusi (l'ancora o le basi dell'evento intersecano l'intervallo escluso)
        is_excluded = self.policy.is_position_excluded(site.position)
        if not is_excluded:
            if obs.event_type == MutationType.DEL:
                is_excluded = any(
                    self.policy.is_position_excluded(p)
                    for p in range(site.position, site.position + obs.length + 1)
                )
            elif obs.event_type == MutationType.INS:
                is_excluded = self.policy.is_position_excluded(site.position + 1)

        std_depth_pass = obs.anchor_depth >= self.policy.indel_min_depth
        std_alt_pass = alt_reads >= self.policy.indel_min_alt_reads
        std_vaf_pass = (self.policy.indel_min_vaf <= alt_vaf <= self.policy.indel_max_vaf)
        std_pass = std_depth_pass and std_alt_pass and std_vaf_pass and strand_bias_pass

        def _eval_restricted_indel_profile() -> Tuple[bool, bool, bool]:
            r_depth = obs.anchor_depth >= self.policy.restricted_heteroplasmy_min_depth
            r_alt = alt_reads >= self.policy.restricted_heteroplasmy_min_alt_reads
            r_vaf = (
                self.policy.restricted_heteroplasmy_min_vaf
                <= alt_vaf
                <= self.policy.restricted_heteroplasmy_max_vaf
            )
            r_pass = r_depth and r_alt and r_vaf and strand_bias_pass
            return r_pass, r_depth, r_alt

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
            restr_pass, r_depth_pass, r_alt_pass = _eval_restricted_indel_profile()
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
                if spatial_damage_suspect or std_pass:
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
            restr_pass, _, _ = _eval_restricted_indel_profile()
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

        # Triage prioritario: Intervallo Escluso vs Slippage omopolimerico
        if is_excluded:
            status = CallStatus.FILTERED
            filter_reason = FilterReason.FILTER_EXCLUDED_INTERVAL
        elif status == CallStatus.PASS and is_slippage and self.policy.filter_homopolymer_slippage:
            status = CallStatus.FILTERED
            filter_reason = FilterReason.FILTER_AMBIGUOUS

        return HeteroplasmyCall(
            position=site.position,
            reference_base=site.reference_base,
            major_allele=major_allele,
            minor_allele=minor_allele,
            callable_depth=obs.anchor_depth,
            alt_reads=alt_reads,
            alt_vaf=alt_vaf,
            forward_ref=fwd_ref_anchor,
            reverse_ref=rev_ref_anchor,
            forward_alt=fwd_indel,
            reverse_alt=rev_indel,
            mean_bq=site.mean_bq,
            mean_mq=site.mean_mq,
            strand_bias_pvalue=strand_bias_pvalue,
            status=status,
            filter_reason=filter_reason,
            gating_action=self.gating_action,
            variant_type=obs.event_type,
            internal_support_count=internal_count,
            terminal_support_count=terminal_count,
            spatial_damage_suspect=spatial_damage_suspect,
            is_slippage_suspect=is_slippage,
            slippage_region=slip_name,
        )

    def evaluate_site(self, site: SitePileup) -> List[HeteroplasmyCall]:
        """
        Evaluate a single observational site breaking silos across SNV and InDel events.
        Passes dominant InDel VAF to prevent phantom SNV calling on structural loci.
        """
        calls: List[HeteroplasmyCall] = []

        dominant_indel_vaf = 0.0
        all_indels = list(site.insertions.values()) + list(site.deletions.values())
        if all_indels:
            dominant_indel_vaf = max(i.vaf for i in all_indels)

        snv_call = self._evaluate_snv_heteroplasmy(site, dominant_indel_vaf=dominant_indel_vaf)
        if snv_call is not None:
            calls.append(snv_call)

        for ins_obs in site.insertions.values():
            ins_call = self._evaluate_indel_observation(site, ins_obs)
            if ins_call is not None:
                calls.append(ins_call)

        for del_obs in site.deletions.values():
            del_call = self._evaluate_indel_observation(site, del_obs)
            if del_call is not None:
                calls.append(del_call)

        return calls

    def call_heteroplasmies(
        self, pileup_stream: Iterable[SitePileup]
    ) -> List[HeteroplasmyCall]:
        """Stream observational sites and compile candidate SNV and InDel calls."""
        calls: List[HeteroplasmyCall] = []
        for site in pileup_stream:
            site_calls = self.evaluate_site(site)
            calls.extend(site_calls)
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
            "variant_type",
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
            "is_slippage_suspect",
            "slippage_region",
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
                    call.variant_type.value,
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
                    call.is_slippage_suspect,
                    call.slippage_region or "NONE",
                ])