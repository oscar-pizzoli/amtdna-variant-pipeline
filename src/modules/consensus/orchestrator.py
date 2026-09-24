"""
Block 6 Orchestrator: Governs input contract resolution, cryptographically verified
execution pipelines, multi-stream consensus and heteroplasmy generation, structural InDel
emission, multi-allelic collapsing, and dual VCF outputs (operational vs subconsensus audit).
"""

from dataclasses import dataclass
import hashlib
import json
import os
import tempfile
import time
from typing import Dict, List, Optional, Tuple

import pysam

from src.modules.consensus.consensus_builder import ConsensusBuilder
from src.modules.consensus.enums import (
    CallStatus,
    ConsensusCallStatus,
    ContaminationEstimateSource,
    GatingAction,
    InputRole,
)
from src.modules.consensus.heteroplasmy_caller import HeteroplasmyCaller
from src.modules.consensus.models import (
    ConsensusPolicy,
    ConsensusQC,
    ConsensusSite,
    HeteroplasmyCall,
    HeteroplasmySiteDetail,
    InputResolution,
    LibraryTreatment,
    ReferenceInfo,
    RegionCoverage,
    SitePileup,
)
from src.modules.consensus.pileup_engine import PileupEngine
from src.modules.synthesis.coordinate_normalizer import MutationType


class ConsensusExecutionError(Exception):
    """Base exception for all Block 6 execution failures."""


class ContractValidationError(ConsensusExecutionError):
    """Raised when upstream JSON contract validation fails."""


class UnresolvedAuthenticityError(ConsensusExecutionError):
    """Raised when Block 4 authenticity status cannot be safely routed."""


@dataclass
class VcfAlleleCandidate:
    """Internal candidate representation for unified multi-allelic collapsing."""
    pos: int
    ref: str
    alt: str
    qual: str
    filter_status: str
    var_type: str
    af: float
    dp: int
    ad_ref: int
    ad_alt: int
    sb_pv: Optional[float] = None
    spatial_suspect: bool = False
    slippage_tract: Optional[str] = None
    is_operational: bool = True
    is_consensus_fixed: bool = False


class ConsensusOrchestrator:
    """
    Orchestrates the entire Block 6 execution pipeline adhering strictly to
    cross-module audit contracts, zero silent fallbacks, and InDel-aware sorted VCF emission.
    """

    HVS1_START = 16024
    HVS1_END = 16383
    HVS2_START = 73
    HVS2_END = 340
    SOFTWARE_VERSION = "1.2.1"

    def __init__(self, policy: Optional[ConsensusPolicy] = None):
        self.policy = policy if policy is not None else ConsensusPolicy()

    @staticmethod
    def calculate_sha256(filepath: str, chunk_size: int = 65536) -> str:
        """Compute SHA-256 hex digest for a target file."""
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Cannot calculate SHA-256; file not found: {filepath}")

        hasher = hashlib.sha256()
        with open(filepath, "rb") as f:
            while chunk := f.read(chunk_size):
                hasher.update(chunk)
        return hasher.hexdigest()

    def resolve_inputs(
        self,
        sample_id: str,
        block4_json_path: str,
        block5_json_path: str,
        reference_fasta_path: str,
        target_contig: Optional[str] = None,
    ) -> InputResolution:
        """
        Validate and resolve upstream contracts from Block 4 and Block 5 JSON outputs.
        Extracts library treatment and ensures exhaustive routing.
        """
        if not os.path.exists(block4_json_path):
            raise FileNotFoundError(f"Block 4 contract JSON not found: {block4_json_path}")
        if not os.path.exists(block5_json_path):
            raise FileNotFoundError(f"Block 5 contract JSON not found: {block5_json_path}")
        if not os.path.exists(reference_fasta_path):
            raise FileNotFoundError(f"Reference FASTA not found: {reference_fasta_path}")

        with open(block4_json_path, "r", encoding="utf-8") as f4:
            b4_data = json.load(f4)

        with open(block5_json_path, "r", encoding="utf-8") as f5:
            b5_data = json.load(f5)

        with pysam.FastaFile(reference_fasta_path) as fasta:
            resolved_contig = target_contig if target_contig else fasta.references[0]
            ref_length = fasta.get_reference_length(resolved_contig)

        ref_sha256 = self.calculate_sha256(reference_fasta_path)
        reference_info = ReferenceInfo(
            contig=resolved_contig,
            length=ref_length,
            sequence_sha256=ref_sha256,
        )

        auth_status = b4_data.get("authentication_status")
        rescaling_performed = b4_data.get("rescaling_performed", False)
        rescaled_bam = b4_data.get("rescaled_bam")
        dedup_bam = b4_data.get("filtered_dedup_bam") or b4_data.get("bam_path")
        damage_model_id = b4_data.get("damage_model_id")

        raw_treatment = str(b4_data.get("library_treatment", "NON_UDG")).upper()
        if raw_treatment in {"HALF", "UDG_HALF", "HALF_UDG"}:
            library_treatment = LibraryTreatment.UDG_HALF
        elif raw_treatment in {"FULL", "UDG_FULL", "FULL_UDG", "UDG"}:
            library_treatment = LibraryTreatment.UDG_FULL
        else:
            library_treatment = LibraryTreatment.NON_UDG

        if auth_status == "AUTHENTIC_ANCIENT_SUPPORTED":
            if not rescaling_performed or not rescaled_bam:
                raise ContractValidationError(
                    f"Sample '{sample_id}' flagged AUTHENTIC_ANCIENT_SUPPORTED but lacks "
                    f"rescaled BAM or rescaling flag in Block 4 contract."
                )
            selected_bam = rescaled_bam
            selected_role = InputRole.RESCALED_ANCIENT

        elif auth_status == "AUTHENTICITY_UNSUPPORTED":
            if not dedup_bam:
                raise ContractValidationError(
                    f"Sample '{sample_id}' flagged AUTHENTICITY_UNSUPPORTED but lacks dedup BAM."
                )
            selected_bam = dedup_bam
            selected_role = InputRole.FILTERED_DEDUP_MODERN

        elif auth_status == "UDG_DEPLETED_EXPECTED":
            selected_bam = rescaled_bam if (rescaling_performed and rescaled_bam) else dedup_bam
            if not selected_bam:
                raise ContractValidationError(
                    f"Sample '{sample_id}' flagged UDG_DEPLETED_EXPECTED but lacks usable BAM."
                )
            selected_role = (
                InputRole.RESCALED_ANCIENT
                if (rescaling_performed and rescaled_bam)
                else InputRole.FILTERED_DEDUP_MODERN
            )

        elif auth_status in {"INCONCLUSIVE", "NOISY_OR_INSUFFICIENT"}:
            raise UnresolvedAuthenticityError(
                f"Sample '{sample_id}' produced unresolved Block 4 authenticity status: '{auth_status}'. "
                f"Consensus calling halted to prevent unvalidated sequence emission."
            )
        else:
            raise ContractValidationError(
                f"Unrecognized or missing authentication_status '{auth_status}' in Block 4 JSON."
            )

        if not os.path.exists(selected_bam):
            raise FileNotFoundError(f"Resolved BAM file does not exist on disk: {selected_bam}")

        bam_sha256 = self.calculate_sha256(selected_bam)

        rec_action_raw = b5_data.get("recommended_action") or b5_data.get("downstream_gating")
        est_confidence = str(b5_data.get("estimation_confidence", "UNKNOWN")).upper()
        contam_fraction = b5_data.get("contamination_fraction")
        contam_source_raw = str(b5_data.get("contamination_estimate_source", "NONE")).upper()

        try:
            contam_source = ContaminationEstimateSource(contam_source_raw)
        except ValueError:
            contam_source = ContaminationEstimateSource.NONE

        if rec_action_raw == "STANDARD_CALLING":
            if est_confidence not in {"HIGH", "MODERATE"}:
                resolved_gating = GatingAction.RESTRICT_HETEROPLASMY
            else:
                resolved_gating = GatingAction.STANDARD_CALLING
        elif rec_action_raw == "RESTRICT_HETEROPLASMY":
            resolved_gating = GatingAction.RESTRICT_HETEROPLASMY
        elif rec_action_raw == "PMD_FILTERING_REQUIRED":
            resolved_gating = GatingAction.PMD_FILTERING_REQUIRED
        elif rec_action_raw == "USE_ENDOGENOUS_SCHMUTZI_CONSENSUS":
            resolved_gating = GatingAction.USE_ENDOGENOUS_SCHMUTZI_CONSENSUS
        elif rec_action_raw in {"EXCLUDE_OR_MANUAL_REVIEW", "REJECT_CONTAMINATED"}:
            resolved_gating = GatingAction.REJECT_CONTAMINATED
        elif rec_action_raw == "INSUFFICIENT_DATA":
            resolved_gating = GatingAction.INSUFFICIENT_DATA
        else:
            raise ContractValidationError(
                f"Unrecognized Block 5 gating directive: '{rec_action_raw}' for sample '{sample_id}'."
            )

        resolution = InputResolution(
            sample_id=sample_id,
            bam_path=selected_bam,
            bam_sha256=bam_sha256,
            input_role=selected_role,
            reference_fasta_path=reference_fasta_path,
            reference_info=reference_info,
            damage_model_id=damage_model_id,
            gating_action=resolved_gating,
            contamination_fraction=contam_fraction,
            contamination_estimate_source=contam_source,
            library_treatment=library_treatment,
        )

        resolution.validate_for_execution()
        return resolution

    def _calculate_region_coverages(
        self, sites: List[ConsensusSite]
    ) -> List[RegionCoverage]:
        """Compute mean callable depth and covered fraction for HVS-I and HVS-II."""
        def _calc(name: str, start: int, end: int) -> RegionCoverage:
            reg_sites = [s for s in sites if start <= s.position <= end]
            if not reg_sites:
                return RegionCoverage(name=name, start=start, end=end, mean_depth=0.0, covered_fraction=0.0)

            total_depth = sum(s.callable_depth for s in reg_sites)
            covered_count = sum(1 for s in reg_sites if s.callable_depth >= self.policy.min_consensus_depth)
            return RegionCoverage(
                name=name,
                start=start,
                end=end,
                mean_depth=total_depth / len(reg_sites),
                covered_fraction=covered_count / len(reg_sites),
            )

        return [
            _calc("HVS-I", self.HVS1_START, self.HVS1_END),
            _calc("HVS-II", self.HVS2_START, self.HVS2_END),
        ]

    def _write_vcf(
        self,
        output_vcf_gz_path: str,
        sample_id: str,
        ref_info: ReferenceInfo,
        consensus_sites: List[ConsensusSite],
        heteroplasmy_calls: List[HeteroplasmyCall],
        pileup_by_pos: Dict[int, SitePileup],
    ) -> None:
        """
        Emit standardized, strictly coordinate-sorted, bgzip-compressed VCF 4.2 files:
        1. Operational VCF (.variants.vcf.gz): Fixed consensus variants (SNVs + InDels) + PASS heteroplasmies.
        2. Subconsensus VCF (.subconsensus_heteroplasmies.vcf.gz): Slippage repeats, damage artifacts, sub-threshold events.
        """
        os.makedirs(os.path.dirname(os.path.abspath(output_vcf_gz_path)), exist_ok=True)
        subconsensus_vcf_gz_path = output_vcf_gz_path.replace(
            ".variants.vcf.gz", ".subconsensus_heteroplasmies.vcf.gz"
        )

        header_lines = [
            "##fileformat=VCFv4.2",
            f"##fileDate={time.strftime('%Y%m%d')}",
            f"##source=Block6_ConsensusOrchestrator_v{self.SOFTWARE_VERSION}",
            f"##reference={ref_info.contig}",
            f"##contig=<ID={ref_info.contig},length={ref_info.length}>",
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total callable depth at locus">',
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency for each ALT allele">',
            '##INFO=<ID=VARTYPE,Number=A,Type=String,Description="Variant classification: MAJOR_CONSENSUS, HETEROPLASMY_SNV, INS, DEL">',
            '##INFO=<ID=SB_PV,Number=A,Type=Float,Description="Fisher exact test strand bias p-value">',
            '##INFO=<ID=SPATIAL_SUSPECT,Number=A,Type=Integer,Description="Flag (1=suspect, 0=authentic) per ALT allele for post-mortem damage">',
            '##INFO=<ID=SLIPPAGE_TRACT,Number=1,Type=String,Description="Homopolymer slippage tract identifier">',
            '##FILTER=<ID=PASS,Description="All quality and evidence thresholds satisfied">',
            '##FILTER=<ID=FILTER_NO_DATA,Description="Zero usable bases observed">',
            '##FILTER=<ID=FILTER_LOW_DEPTH,Description="Callable depth below calling threshold">',
            '##FILTER=<ID=FILTER_LOW_ALT_READS,Description="Alternative read support below threshold">',
            '##FILTER=<ID=FILTER_VAF_RANGE,Description="Variant allele frequency outside target window">',
            '##FILTER=<ID=FILTER_STRAND_BIAS,Description="Severe strand bias detected by Fisher exact test">',
            '##FILTER=<ID=FILTER_CONTAMINATION,Description="Restricted due to contamination gating">',
            '##FILTER=<ID=FILTER_INSUFFICIENT_DATA,Description="Gating restricted due to unassessed contamination">',
            '##FILTER=<ID=FILTER_DAMAGE_SUSPECT,Description="Filtered due to suspected post-mortem damage">',
            '##FILTER=<ID=FILTER_AMBIGUOUS,Description="Site is ambiguous, contested, or inside homopolymer slippage tract">',
            '##FILTER=<ID=FILTER_REFERENCE_MISMATCH,Description="Reference base is non-canonical or legacy placeholder">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Callable depth">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (REF, ALT1, ALT2...)">',
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}",
        ]

        candidates_by_pos: Dict[int, List[VcfAlleleCandidate]] = {}

        # ----------------------------------------------------------------------
        # Step A: Collect Fixed Consensus Alterations (SNVs & Structural InDels)
        # ----------------------------------------------------------------------
        for site in consensus_sites:
            pos = site.position
            ref = site.reference_base.upper()
            canon = site.canonical_base.upper()
            p_site = pileup_by_pos.get(pos)

            if (
                site.call_status == ConsensusCallStatus.CALLED
                and canon in {"A", "C", "G", "T"}
                and canon != ref
            ):
                dp = site.callable_depth
                ref_cnt = p_site.counts.get(ref) if p_site else 0
                alt_cnt = p_site.counts.get(canon) if p_site else dp
                if (ref_cnt + alt_cnt) == 0 and dp > 0:
                    alt_cnt = dp

                cand = VcfAlleleCandidate(
                    pos=pos,
                    ref=ref,
                    alt=canon,
                    qual="60",
                    filter_status="PASS",
                    var_type="MAJOR_CONSENSUS",
                    af=site.major_vaf,
                    dp=dp,
                    ad_ref=ref_cnt,
                    ad_alt=alt_cnt,
                    is_operational=True,
                    is_consensus_fixed=True,
                )
                candidates_by_pos.setdefault(pos, []).append(cand)

        # ----------------------------------------------------------------------
        # Step B: Collect Heteroplasmy and InDel Observations
        # ----------------------------------------------------------------------
        for call in heteroplasmy_calls:
            pos = call.position
            ref_b = call.reference_base.upper()

            if call.variant_type == MutationType.SNV:
                vcf_ref = ref_b
                # Risoluzione dinamica della polarità per eteroplasmie a frequenza invertita
                if call.minor_allele.upper() != ref_b:
                    vcf_alt = call.minor_allele.upper()
                    af_val = call.alt_vaf
                    ad_alt = call.alt_reads
                    ad_ref = call.forward_ref + call.reverse_ref
                elif call.major_allele.upper() != ref_b:
                    vcf_alt = call.major_allele.upper()
                    af_val = round(1.0 - call.alt_vaf, 4) if call.alt_vaf is not None else 0.60
                    ad_alt = max(0, call.callable_depth - call.alt_reads)
                    ad_ref = call.alt_reads
                else:
                    continue
                var_type = "HETEROPLASMY_SNV"

            elif call.variant_type == MutationType.INS:
                vcf_ref = call.major_allele.upper()
                vcf_alt = call.minor_allele.upper()
                af_val = call.alt_vaf
                ad_alt = call.alt_reads
                ad_ref = call.forward_ref + call.reverse_ref
                var_type = "INS"

            elif call.variant_type == MutationType.DEL:
                vcf_ref = call.major_allele.upper()
                vcf_alt = call.minor_allele.upper()
                af_val = call.alt_vaf
                ad_alt = call.alt_reads
                ad_ref = call.forward_ref + call.reverse_ref
                var_type = "DEL"

            else:
                continue

            if vcf_ref == vcf_alt:
                continue

            is_fixed_indel = af_val >= self.policy.major_allele_fraction
            if is_fixed_indel:
                var_type = "MAJOR_CONSENSUS"

            status_str = call.status.value
            filter_str = "PASS" if status_str == "PASS" else call.filter_reason.value

            is_op = (call.status == CallStatus.PASS and not call.is_slippage_suspect) or is_fixed_indel

            cand = VcfAlleleCandidate(
                pos=pos,
                ref=vcf_ref,
                alt=vcf_alt,
                qual=".",
                filter_status=filter_str,
                var_type=var_type,
                af=af_val,
                dp=call.callable_depth,
                ad_ref=ad_ref,
                ad_alt=ad_alt,
                sb_pv=call.strand_bias_pvalue,
                spatial_suspect=call.spatial_damage_suspect,
                slippage_tract=call.slippage_region,
                is_operational=is_op,
                is_consensus_fixed=is_fixed_indel,
            )

            existing = candidates_by_pos.setdefault(pos, [])
            already_present = any(
                e.ref == cand.ref and e.alt == cand.alt for e in existing
            )
            if not already_present:
                existing.append(cand)

        # ----------------------------------------------------------------------
        # Step C: Multi-Allelic Consolidation & Strict Coordinate Sorting
        # ----------------------------------------------------------------------
        op_lines: List[Tuple[int, str]] = []
        sub_lines: List[Tuple[int, str]] = []

        all_positions = sorted(candidates_by_pos.keys())

        for pos in all_positions:
            cands = candidates_by_pos[pos]

            op_cands = [c for c in cands if c.is_operational]
            sub_cands = [c for c in cands if not c.is_operational]

            for stream_cands, target_list in [(op_cands, op_lines), (sub_cands, sub_lines)]:
                if not stream_cands:
                    continue

                longest_ref = max((c.ref for c in stream_cands), key=len)
                norm_alts: List[str] = []
                for c in stream_cands:
                    if len(c.ref) < len(longest_ref):
                        padding = longest_ref[len(c.ref):]
                        norm_alts.append(c.alt + padding)
                    else:
                        norm_alts.append(c.alt)

                ref_field = longest_ref
                alt_field = ",".join(norm_alts)

                all_filters = {c.filter_status for c in stream_cands}
                filter_field = "PASS" if all_filters == {"PASS"} else ";".join(sorted(all_filters))

                dp_val = max(c.dp for c in stream_cands)
                af_field = ",".join(f"{c.af:.4f}" for c in stream_cands)
                vartype_field = ",".join(c.var_type for c in stream_cands)

                info_parts = [
                    f"DP={dp_val}",
                    f"AF={af_field}",
                    f"VARTYPE={vartype_field}",
                ]

                sb_vals = [f"{c.sb_pv:.4e}" if c.sb_pv is not None else "." for c in stream_cands]
                if any(v != "." for v in sb_vals):
                    info_parts.append(f"SB_PV={','.join(sb_vals)}")

                spatial_flags = ["1" if c.spatial_suspect else "0" for c in stream_cands]
                if any(c.spatial_suspect for c in stream_cands):
                    info_parts.append(f"SPATIAL_SUSPECT={','.join(spatial_flags)}")

                slip_tracts = {c.slippage_tract for c in stream_cands if c.slippage_tract}
                if slip_tracts:
                    info_parts.append(f"SLIPPAGE_TRACT={';'.join(sorted(slip_tracts))}")

                info_field = ";".join(info_parts)

                ad_ref_val = min(c.ad_ref for c in stream_cands)
                ad_parts = [str(ad_ref_val)] + [str(c.ad_alt) for c in stream_cands]
                ad_field = ",".join(ad_parts)

                if any(c.is_consensus_fixed for c in stream_cands):
                    gt_field = "1" if len(stream_cands) == 1 else "1/2"
                else:
                    gt_field = "0/1" if len(stream_cands) == 1 else "1/2"

                fmt_field = "GT:DP:AD"
                val_field = f"{gt_field}:{dp_val}:{ad_field}"

                line = f"{ref_info.contig}\t{pos}\t.\t{ref_field}\t{alt_field}\t60\t{filter_field}\t{info_field}\t{fmt_field}\t{val_field}\n"
                target_list.append((pos, line))

        # ----------------------------------------------------------------------
        # Step D: Transacted Sorted Writing and Indexing
        # ----------------------------------------------------------------------
        def _write_compressed_vcf(entries: List[Tuple[int, str]], out_gz_path: str) -> None:
            entries.sort(key=lambda x: x[0])

            with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False, encoding="utf-8") as tmp_vcf:
                tmp_path = tmp_vcf.name
                for h in header_lines:
                    tmp_vcf.write(h + "\n")
                for _, row in entries:
                    tmp_vcf.write(row)

            try:
                pysam.tabix_compress(tmp_path, out_gz_path, force=True)
                pysam.tabix_index(out_gz_path, preset="vcf", force=True)
            finally:
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)

        _write_compressed_vcf(op_lines, output_vcf_gz_path)

        if self.policy.emit_subconsensus_vcf:
            _write_compressed_vcf(sub_lines, subconsensus_vcf_gz_path)

    def _write_qc_json(self, qc: ConsensusQC, output_path: str) -> None:
        """
        Serialize global consensus QC metrics and strict string-keyed heteroplasmy details.
        Guarantees self-contained position, depth, and spatial provenance across blocks.
        """
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

        serialized_details = {
            str(pos): {
                "position": detail.position,
                "ref": detail.ref,
                "alt": detail.alt,
                "depth": detail.depth,
                "vaf": round(detail.vaf, 4),
                "internal_support_count": detail.internal_support_count,
                "terminal_support_count": detail.terminal_support_count,
                "spatial_damage_suspect": detail.spatial_damage_suspect,
            }
            for pos, detail in qc.heteroplasmy_details.items()
        }

        data = {
            "sample_id": qc.sample_id,
            "input_bam_path": qc.input_bam_path,
            "input_bam_sha256": qc.input_bam_sha256,
            "reference_contig": qc.reference_contig,
            "reference_sha256": qc.reference_sha256,
            "reference_length": qc.reference_length,
            "consensus_fasta_sha256": qc.consensus_fasta_sha256,
            "mean_depth": round(qc.mean_depth, 2),
            "covered_positions_fraction": round(qc.covered_positions_fraction, 4),
            "n_fraction_canonical": round(qc.n_fraction_canonical, 4),
            "iupac_fraction": round(qc.iupac_fraction, 4),
            "heteroplasmy_sites_total": qc.heteroplasmy_sites_total,
            "heteroplasmy_sites_pass": qc.heteroplasmy_sites_pass,
            "heteroplasmy_sites_filtered": qc.heteroplasmy_sites_filtered,
            "regions": [
                {
                    "name": r.name,
                    "start": r.start,
                    "end": r.end,
                    "mean_depth": round(r.mean_depth, 2),
                    "covered_fraction": round(r.covered_fraction, 4),
                }
                for r in qc.regions
            ],
            "gating_action": qc.gating_action.value,
            "software_version": qc.software_version,
            "policy_version": qc.policy_version,
            "execution_time_seconds": round(qc.execution_time_seconds, 2),
            "heteroplasmy_details": serialized_details,
        }
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)

    def run(
        self,
        input_resolution: InputResolution,
        output_dir: str,
    ) -> ConsensusQC:
        """
        Execute the primary Block 6 pipeline:
        1. Calibrate Policy according to LibraryTreatment
        2. Observational extraction (PileupEngine)
        3. Majority consensus (ConsensusBuilder)
        4. Spatial Heteroplasmy Calling (HeteroplasmyCaller)
        5. Artifact generation & String-Keyed Heteroplasmy Details Export
        """
        start_time = time.perf_counter()
        os.makedirs(output_dir, exist_ok=True)

        sample_id = input_resolution.sample_id
        ref_info = input_resolution.reference_info

        treatment_enum = input_resolution.library_treatment
        resolved_threshold = self.policy.spatial_thresholds.get(treatment_enum, 5)

        calibrated_policy = ConsensusPolicy(
            policy_id=self.policy.policy_id,
            policy_version=self.policy.policy_version,
            min_base_quality=self.policy.min_base_quality,
            min_mapping_quality=self.policy.min_mapping_quality,
            min_consensus_depth=self.policy.min_consensus_depth,
            major_allele_fraction=self.policy.major_allele_fraction,
            ambiguity_min_vaf=self.policy.ambiguity_min_vaf,
            heteroplasmy_min_depth=self.policy.heteroplasmy_min_depth,
            heteroplasmy_min_alt_reads=self.policy.heteroplasmy_min_alt_reads,
            heteroplasmy_min_vaf=self.policy.heteroplasmy_min_vaf,
            heteroplasmy_max_vaf=self.policy.heteroplasmy_max_vaf,
            strand_bias_alpha=self.policy.strand_bias_alpha,
            strand_bias_min_strand_reads=self.policy.strand_bias_min_strand_reads,
            enable_indel_calling=self.policy.enable_indel_calling,
            indel_min_depth=self.policy.indel_min_depth,
            indel_min_alt_reads=self.policy.indel_min_alt_reads,
            indel_min_vaf=self.policy.indel_min_vaf,
            indel_max_vaf=self.policy.indel_max_vaf,
            filter_homopolymer_slippage=self.policy.filter_homopolymer_slippage,
            emit_iupac_fasta=self.policy.emit_iupac_fasta,
            emit_subconsensus_vcf=self.policy.emit_subconsensus_vcf,
            mask_char=self.policy.mask_char,
            restricted_heteroplasmy_min_depth=self.policy.restricted_heteroplasmy_min_depth,
            restricted_heteroplasmy_min_alt_reads=self.policy.restricted_heteroplasmy_min_alt_reads,
            restricted_heteroplasmy_min_vaf=self.policy.restricted_heteroplasmy_min_vaf,
            restricted_heteroplasmy_max_vaf=self.policy.restricted_heteroplasmy_max_vaf,
            min_internal_fragments=self.policy.min_internal_fragments,
            spatial_end_distance_threshold=resolved_threshold,
            spatial_thresholds=self.policy.spatial_thresholds,
            slippage_regions=self.policy.slippage_regions,
        )

        pileup_engine = PileupEngine(
            bam_path=input_resolution.bam_path,
            reference_fasta_path=input_resolution.reference_fasta_path,
            policy=calibrated_policy,
            target_contig=ref_info.contig,
            expected_reference_length=ref_info.length,
        )

        pileup_sites: List[SitePileup] = list(pileup_engine.iter_sites())
        pileup_by_pos: Dict[int, SitePileup] = {s.position: s for s in pileup_sites}

        consensus_builder = ConsensusBuilder(
            policy=calibrated_policy,
            sample_id=sample_id,
            reference_info=ref_info,
        )
        (
            consensus_sites,
            canonical_seq,
            iupac_seq,
            indel_aware_seq,
        ) = consensus_builder.build_consensus(pileup_sites)

        if len(canonical_seq) != ref_info.length:
            raise ConsensusExecutionError(
                f"Canonical FASTA length mismatch: expected {ref_info.length}, got {len(canonical_seq)}"
            )

        heteroplasmy_caller = HeteroplasmyCaller(
            policy=calibrated_policy,
            sample_id=sample_id,
            reference_info=ref_info,
            gating_action=input_resolution.gating_action,
        )
        heteroplasmy_calls = heteroplasmy_caller.call_heteroplasmies(pileup_sites)

        canonical_fasta_path = os.path.join(output_dir, f"{sample_id}.consensus.fasta")
        iupac_fasta_path = os.path.join(output_dir, f"{sample_id}.consensus_iupac.fasta")
        indel_fasta_path = os.path.join(output_dir, f"{sample_id}.consensus_indel.fasta")
        heteroplasmy_tsv_path = os.path.join(output_dir, f"{sample_id}.heteroplasmy.tsv")
        vcf_path = os.path.join(output_dir, f"{sample_id}.variants.vcf.gz")
        qc_json_path = os.path.join(output_dir, f"{sample_id}.consensus_qc.json")

        consensus_builder.write_fasta(
            output_path=canonical_fasta_path,
            header_id=sample_id,
            description=f"canonical_consensus contig={ref_info.contig} length={ref_info.length}",
            sequence=canonical_seq,
        )

        if calibrated_policy.emit_iupac_fasta:
            consensus_builder.write_fasta(
                output_path=iupac_fasta_path,
                header_id=sample_id,
                description=f"iupac_consensus contig={ref_info.contig} length={ref_info.length}",
                sequence=iupac_seq,
            )

        consensus_builder.write_fasta(
            output_path=indel_fasta_path,
            header_id=sample_id,
            description=f"consensus indel-aware circular contig={ref_info.contig} length={ref_info.length}",
            sequence=indel_aware_seq,
        )

        fasta_sha256 = self.calculate_sha256(canonical_fasta_path)

        heteroplasmy_caller.write_tsv(
            calls=heteroplasmy_calls,
            output_path=heteroplasmy_tsv_path,
        )

        self._write_vcf(
            output_vcf_gz_path=vcf_path,
            sample_id=sample_id,
            ref_info=ref_info,
            consensus_sites=consensus_sites,
            heteroplasmy_calls=heteroplasmy_calls,
            pileup_by_pos=pileup_by_pos,
        )

        het_details_map: Dict[str, HeteroplasmySiteDetail] = {}
        for c in heteroplasmy_calls:
            ref_b = c.reference_base.upper()
            if c.variant_type == MutationType.SNV:
                if c.minor_allele.upper() != ref_b:
                    alt_b = c.minor_allele.upper()
                    vaf_b = c.alt_vaf
                elif c.major_allele.upper() != ref_b:
                    alt_b = c.major_allele.upper()
                    vaf_b = round(1.0 - c.alt_vaf, 4) if c.alt_vaf is not None else 0.60
                else:
                    alt_b = c.minor_allele.upper()
                    vaf_b = c.alt_vaf
            else:
                alt_b = c.minor_allele.upper()
                vaf_b = c.alt_vaf

            het_details_map[str(c.position)] = HeteroplasmySiteDetail(
                position=c.position,
                ref=ref_b,
                alt=alt_b,
                depth=c.callable_depth,
                vaf=vaf_b,
                internal_support_count=c.internal_support_count,
                terminal_support_count=c.terminal_support_count,
                spatial_damage_suspect=c.spatial_damage_suspect,
            )

        elapsed_sec = time.perf_counter() - start_time
        total_depth = sum(s.callable_depth for s in consensus_sites)
        mean_depth = total_depth / ref_info.length if ref_info.length > 0 else 0.0

        covered_positions = sum(1 for s in consensus_sites if s.callable_depth >= 1)
        covered_fraction = covered_positions / ref_info.length if ref_info.length > 0 else 0.0

        n_count = canonical_seq.count("N")
        n_fraction = n_count / ref_info.length if ref_info.length > 0 else 0.0

        iupac_symbols = {"R", "Y", "S", "W", "K", "M"}
        iupac_count = sum(1 for char in iupac_seq if char in iupac_symbols)
        iupac_fraction = iupac_count / ref_info.length if ref_info.length > 0 else 0.0

        pass_het = sum(1 for c in heteroplasmy_calls if c.status == CallStatus.PASS)
        filtered_het = sum(1 for c in heteroplasmy_calls if c.status == CallStatus.FILTERED)

        region_coverages = self._calculate_region_coverages(consensus_sites)

        qc = ConsensusQC(
            sample_id=sample_id,
            input_bam_path=input_resolution.bam_path,
            input_bam_sha256=input_resolution.bam_sha256,
            reference_contig=ref_info.contig,
            reference_sha256=ref_info.sequence_sha256,
            reference_length=ref_info.length,
            consensus_fasta_sha256=fasta_sha256,
            mean_depth=mean_depth,
            covered_positions_fraction=covered_fraction,
            n_fraction_canonical=n_fraction,
            iupac_fraction=iupac_fraction,
            heteroplasmy_sites_total=len(heteroplasmy_calls),
            heteroplasmy_sites_pass=pass_het,
            heteroplasmy_sites_filtered=filtered_het,
            regions=region_coverages,
            gating_action=input_resolution.gating_action,
            software_version=self.SOFTWARE_VERSION,
            policy_version=calibrated_policy.policy_version,
            execution_time_seconds=elapsed_sec,
            heteroplasmy_details=het_details_map,
        )

        self._write_qc_json(qc, qc_json_path)
        return qc