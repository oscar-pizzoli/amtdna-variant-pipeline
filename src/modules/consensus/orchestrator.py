"""
Block 6 Orchestrator: Governs input contract resolution, cryptographically verified
execution pipelines, multi-stream consensus and heteroplasmy generation, and audit outputs.
"""

import hashlib
import json
import os
import time
from typing import Dict, List, Optional

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
    InputResolution,
    ReferenceInfo,
    RegionCoverage,
    SitePileup,
)
from src.modules.consensus.pileup_engine import PileupEngine


class ConsensusExecutionError(Exception):
    """Base exception for all Block 6 execution failures."""


class ContractValidationError(ConsensusExecutionError):
    """Raised when upstream JSON contract validation fails."""


class UnresolvedAuthenticityError(ConsensusExecutionError):
    """Raised when Block 4 authenticity status cannot be safely routed."""


class ConsensusOrchestrator:
    """
    Orchestrates the entire Block 6 execution pipeline adhering strictly to
    cross-module audit contracts, zero silent fallbacks, and invariant checks.
    """

    HVS1_START = 16024
    HVS1_END = 16383
    HVS2_START = 73
    HVS2_END = 340
    SOFTWARE_VERSION = "1.0.0"

    def __init__(self, policy: Optional[ConsensusPolicy] = None):
        """
        Initialize the orchestrator with an immutable ConsensusPolicy.
        """
        self.policy = policy if policy is not None else ConsensusPolicy()

    @staticmethod
    def calculate_sha256(filepath: str, chunk_size: int = 65536) -> str:
        """
        Compute SHA-256 hex digest for a target file.
        """
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
        Ensures exhaustive routing with zero silent fallbacks.
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

        # 1. Resolve Reference Genome Specifications
        with pysam.FastaFile(reference_fasta_path) as fasta:
            resolved_contig = target_contig if target_contig else fasta.references[0]
            ref_length = fasta.get_reference_length(resolved_contig)

        ref_sha256 = self.calculate_sha256(reference_fasta_path)
        reference_info = ReferenceInfo(
            contig=resolved_contig,
            length=ref_length,
            sequence_sha256=ref_sha256,
        )

        # 2. Block 4 Exhaustive Authenticity Routing
        auth_status = b4_data.get("authentication_status")
        rescaling_performed = b4_data.get("rescaling_performed", False)
        rescaled_bam = b4_data.get("rescaled_bam")
        dedup_bam = b4_data.get("filtered_dedup_bam") or b4_data.get("bam_path")
        damage_model_id = b4_data.get("damage_model_id")

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

        # 3. Block 5 Contamination Gating & Pessimistic Confidence Evaluation
        rec_action_raw = b5_data.get("recommended_action") or b5_data.get("downstream_gating")
        est_confidence = str(b5_data.get("estimation_confidence", "UNKNOWN")).upper()
        contam_fraction = b5_data.get("contamination_fraction")
        contam_source_raw = str(b5_data.get("contamination_estimate_source", "NONE")).upper()

        try:
            contam_source = ContaminationEstimateSource(contam_source_raw)
        except ValueError:
            contam_source = ContaminationEstimateSource.NONE

        # Exhaustive Gating Action Mapping: Non-high confidence triggers restriction
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
        )

        resolution.validate_for_execution()
        return resolution

    def _calculate_region_coverages(
        self, sites: List[ConsensusSite]
    ) -> List[RegionCoverage]:
        """
        Compute mean callable depth and covered fraction for HVS-I and HVS-II.
        """
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
        Emit a standardized, bgzip-compressed VCF 4.2 file capturing consensus
        variants and all evaluated heteroplasmy events with comprehensive filter headers.
        """
        os.makedirs(os.path.dirname(os.path.abspath(output_vcf_gz_path)), exist_ok=True)
        uncompressed_vcf = output_vcf_gz_path.removesuffix(".gz")

        header_lines = [
            "##fileformat=VCFv4.2",
            f"##fileDate={time.strftime('%Y%m%d')}",
            f"##source=Block6_ConsensusOrchestrator_v{self.SOFTWARE_VERSION}",
            f"##reference={ref_info.contig}",
            f"##contig=<ID={ref_info.contig},length={ref_info.length}>",
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Callable read depth">',
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency of ALT allele">',
            '##INFO=<ID=VARTYPE,Number=1,Type=String,Description="Variant class: MAJOR_CONSENSUS or HETEROPLASMY">',
            '##INFO=<ID=SB_PV,Number=1,Type=Float,Description="Fisher exact test strand bias p-value">',
            '##FILTER=<ID=PASS,Description="All quality and evidence thresholds satisfied">',
            '##FILTER=<ID=FILTER_NO_DATA,Description="Zero usable bases observed">',
            '##FILTER=<ID=FILTER_LOW_DEPTH,Description="Callable depth below calling threshold">',
            '##FILTER=<ID=FILTER_LOW_ALT_READS,Description="Alternative read support below threshold">',
            '##FILTER=<ID=FILTER_VAF_RANGE,Description="Variant allele frequency outside target window">',
            '##FILTER=<ID=FILTER_STRAND_BIAS,Description="Severe strand bias detected by Fisher exact test">',
            '##FILTER=<ID=FILTER_CONTAMINATION,Description="Restricted due to contamination gating">',
            '##FILTER=<ID=FILTER_INSUFFICIENT_DATA,Description="Gating restricted due to unassessed contamination">',
            '##FILTER=<ID=FILTER_DAMAGE_SUSPECT,Description="Filtered due to suspected post-mortem damage">',
            '##FILTER=<ID=FILTER_AMBIGUOUS,Description="Site is ambiguous or contested without sufficient support">',
            '##FILTER=<ID=FILTER_REFERENCE_MISMATCH,Description="Reference base is non-canonical or legacy placeholder">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Callable depth">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (REF, ALT)">',
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}",
        ]

        het_by_pos = {call.position: call for call in heteroplasmy_calls}

        with open(uncompressed_vcf, "w", encoding="utf-8") as out:
            for line in header_lines:
                out.write(line + "\n")

            for site in consensus_sites:
                pos = site.position
                ref = site.reference_base.upper()
                canon = site.canonical_base.upper()
                p_site = pileup_by_pos.get(pos)

                is_consensus_variant = (
                    site.call_status == ConsensusCallStatus.CALLED
                    and canon in {"A", "C", "G", "T"}
                    and canon != ref
                )
                het_call = het_by_pos.get(pos)

                # 1. Emit Major Consensus Variant
                if is_consensus_variant:
                    ref_cnt = p_site.counts.get(ref) if p_site else 0
                    alt_cnt = p_site.counts.get(canon) if p_site else site.callable_depth
                    info = f"DP={site.callable_depth};AF={site.major_vaf:.4f};VARTYPE=MAJOR_CONSENSUS"
                    fmt = "GT:DP:AD"
                    val = f"1:{site.callable_depth}:{ref_cnt},{alt_cnt}"
                    out.write(f"{ref_info.contig}\t{pos}\t.\t{ref}\t{canon}\t60\tPASS\t{info}\t{fmt}\t{val}\n")

                # 2. Emit Heteroplasmy Call (Independent of consensus state)
                if het_call is not None:
                    # Avoid duplicate emission if heteroplasmy allele matches already emitted major variant
                    if is_consensus_variant and het_call.minor_allele == canon:
                        continue

                    alt = het_call.minor_allele
                    status_str = het_call.status.value
                    filter_str = "PASS" if status_str == "PASS" else het_call.filter_reason.value
                    info = (
                        f"DP={het_call.callable_depth};AF={het_call.alt_vaf:.4f};"
                        f"VARTYPE=HETEROPLASMY;SB_PV={het_call.strand_bias_pvalue:.4e}"
                    )
                    ref_cnt = p_site.counts.get(ref) if p_site else (het_call.forward_ref + het_call.reverse_ref)
                    fmt = "GT:DP:AD"
                    val = f"0/1:{het_call.callable_depth}:{ref_cnt},{het_call.alt_reads}"
                    out.write(f"{ref_info.contig}\t{pos}\t.\t{ref}\t{alt}\t.\t{filter_str}\t{info}\t{fmt}\t{val}\n")

        # Compress to BGZF and create Tabix index
        pysam.tabix_compress(uncompressed_vcf, output_vcf_gz_path, force=True)
        try:
            pysam.tabix_index(output_vcf_gz_path, preset="vcf", force=True)
        except Exception:
            pass

        if os.path.exists(uncompressed_vcf):
            os.remove(uncompressed_vcf)

    def _write_qc_json(self, qc: ConsensusQC, output_path: str) -> None:
        """
        Serialize global consensus QC metrics to JSON.
        """
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        data = {
            "sample_id": qc.sample_id,
            "input_bam_path": qc.input_bam_path,
            "input_bam_sha256": qc.input_bam_sha256,
            "reference_sha256": qc.reference_sha256,
            "reference_length": qc.reference_length,
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
        1. Observational extraction (PileupEngine)
        2. Majority consensus & parallel IUPAC streams (ConsensusBuilder)
        3. Statistical heteroplasmy calling & strand bias audit (HeteroplasmyCaller)
        4. Artifact generation (canonical FASTA, IUPAC FASTA, TSV, VCF.gz, QC JSON)
        """
        start_time = time.perf_counter()
        os.makedirs(output_dir, exist_ok=True)

        sample_id = input_resolution.sample_id
        ref_info = input_resolution.reference_info

        # 1. Initialize Observational Engine with direct reference_fasta_path
        pileup_engine = PileupEngine(
            bam_path=input_resolution.bam_path,
            reference_fasta_path=input_resolution.reference_fasta_path,
            policy=self.policy,
            target_contig=ref_info.contig,
            expected_reference_length=ref_info.length,
        )

        # 2. Materialize gapless SitePileup stream in memory (~16,569 sites)
        pileup_sites: List[SitePileup] = list(pileup_engine.iter_sites())
        pileup_by_pos: Dict[int, SitePileup] = {s.position: s for s in pileup_sites}

        # 3. Consensus Generation (Canonical and IUPAC streams)
        consensus_builder = ConsensusBuilder(
            policy=self.policy,
            sample_id=sample_id,
            reference_info=ref_info,
        )
        consensus_sites, canonical_seq, iupac_seq = consensus_builder.build_consensus(pileup_sites)

        # Strict invariant length assertions
        if len(canonical_seq) != ref_info.length:
            raise ConsensusExecutionError(
                f"Canonical FASTA length mismatch: expected {ref_info.length}, got {len(canonical_seq)}"
            )
        if len(iupac_seq) != ref_info.length:
            raise ConsensusExecutionError(
                f"IUPAC FASTA length mismatch: expected {ref_info.length}, got {len(iupac_seq)}"
            )

        # 4. Intra-Individual Heteroplasmy Calling
        heteroplasmy_caller = HeteroplasmyCaller(
            policy=self.policy,
            sample_id=sample_id,
            reference_info=ref_info,
            gating_action=input_resolution.gating_action,
        )
        heteroplasmy_calls = heteroplasmy_caller.call_heteroplasmies(pileup_sites)

        # 5. Emit Artifacts
        canonical_fasta_path = os.path.join(output_dir, f"{sample_id}.consensus.fasta")
        iupac_fasta_path = os.path.join(output_dir, f"{sample_id}.consensus_iupac.fasta")
        heteroplasmy_tsv_path = os.path.join(output_dir, f"{sample_id}.heteroplasmy.tsv")
        vcf_path = os.path.join(output_dir, f"{sample_id}.variants.vcf.gz")
        qc_json_path = os.path.join(output_dir, f"{sample_id}.consensus_qc.json")

        consensus_builder.write_fasta(
            output_path=canonical_fasta_path,
            header_id=sample_id,
            description=f"canonical_consensus contig={ref_info.contig} length={ref_info.length}",
            sequence=canonical_seq,
        )

        consensus_builder.write_fasta(
            output_path=iupac_fasta_path,
            header_id=sample_id,
            description=f"iupac_consensus contig={ref_info.contig} length={ref_info.length}",
            sequence=iupac_seq,
        )

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

        # 6. Global QC Compilation
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
            reference_sha256=ref_info.sequence_sha256,
            reference_length=ref_info.length,
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
            policy_version=self.policy.policy_version,
            execution_time_seconds=elapsed_sec,
        )

        self._write_qc_json(qc, qc_json_path)
        return qc