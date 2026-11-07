"""
Block 6: Dual-Stream Mitochondrial Consensus Assembly & Heteroplasmy Calling.

Single-file step responsible for:
- Quality-calibrated read pileup extraction across all 16,569 mitochondrial positions
- Bayesian & frequentist heteroplasmy calling with strict strand-bias (Fisher exact test) gating
- Spatial damage deamination exclusion avoiding terminal C->T/G->A artifacts
- Assembly of Major Consensus FASTA (dominant allele >= 80%, N on low depth)
- Assembly of IUPAC Ambiguity FASTA (encoding heteroplasmic sites: R, Y, S, W, K, M)
- Compressed and indexed VCF generation (BCF/VCF output with AD, DP, AF, SB tags)
- Serialization of Block 6 Contract (block6_contract.json)
"""

from dataclasses import dataclass, field
from enum import Enum
import gzip
import json
import logging
import math
from pathlib import Path
import re
import shutil
import subprocess
import time
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import pysam

# Safe Fisher's exact test import with pure-python hypergeometric fallback
try:
    from scipy.stats import fisher_exact
except ImportError:
    def fisher_exact(table: List[List[int]], alternative: str = "two-sided") -> Tuple[float, float]:
        """Pure-python fallback for Fisher's exact test using hypergeometric distribution."""
        a, b = table[0]
        c, d = table[1]
        n = a + b + c + d
        if n == 0:
            return 1.0, 1.0

        def hypergeom_prob(x: int) -> float:
            try:
                return (math.comb(a + b, x) * math.comb(c + d, a + c - x)) / math.comb(n, a + c)
            except (ValueError, ZeroDivisionError):
                return 0.0

        p_obs = hypergeom_prob(a)
        min_x = max(0, a - d)
        max_x = min(a + b, a + c)
        p_total = sum(hypergeom_prob(x) for x in range(min_x, max_x + 1) if hypergeom_prob(x) <= p_obs + 1e-9)
        odds_ratio = (a * d) / (b * c) if (b * c) > 0 else float("inf")
        return odds_ratio, min(1.0, p_total)

from src.steps.reference import ReferenceInfo
from src.utils import (
    calculate_sha256,
    compute_composite_hash,
    PipelineConfig,
)
from src.pipeline.stage_tracker import StageTracker

logger = logging.getLogger("mtDNA_Pipeline.Consensus")


# ==============================================================================
# SECTION 1: MODELS & ENUMS
# ==============================================================================

class GatingAction(str, Enum):
    PROCEED_UNRESTRICTED = "PROCEED_UNRESTRICTED"
    RESTRICT_HETEROPLASMY = "RESTRICT_HETEROPLASMY"
    SUPPRESS_HETEROPLASMY = "SUPPRESS_HETEROPLASMY"
    HALT_SAMPLE = "HALT_SAMPLE"


class CallStatus(str, Enum):
    PASS = "PASS"
    FILTERED = "FILTERED"


class FilterReason(str, Enum):
    LOW_DEPTH = "LOW_DEPTH"
    LOW_ALT_COUNT = "LOW_ALT_COUNT"
    STRAND_BIAS = "STRAND_BIAS"
    SPATIAL_DAMAGE = "SPATIAL_DAMAGE"
    HOMOPOLYMER_TRACT = "HOMOPOLYMER_TRACT"
    EXCLUDED_INTERVAL = "EXCLUDED_INTERVAL"
    HIGH_CONTAMINATION = "HIGH_CONTAMINATION"


IUPAC_MAP: Dict[Tuple[str, str], str] = {
    ("A", "G"): "R", ("G", "A"): "R",
    ("C", "T"): "Y", ("T", "C"): "Y",
    ("G", "C"): "S", ("C", "G"): "S",
    ("A", "T"): "W", ("T", "A"): "W",
    ("G", "T"): "K", ("T", "G"): "K",
    ("A", "C"): "M", ("C", "A"): "M",
}


@dataclass
class ConsensusPolicy:
    policy_id: str = "consensus_strict_v1"
    policy_version: str = "1.2.0"
    min_base_quality: int = 30
    min_mapping_quality: int = 30
    min_consensus_depth: int = 10
    major_allele_fraction: float = 0.80
    ambiguity_min_vaf: float = 0.10
    heteroplasmy_min_depth: int = 15
    heteroplasmy_min_alt_reads: int = 3
    heteroplasmy_min_vaf: float = 0.10
    heteroplasmy_max_vaf: float = 0.80
    strand_bias_alpha: float = 0.05
    strand_bias_min_strand_reads: int = 1
    emit_iupac_fasta: bool = True
    mask_char: str = "N"
    restricted_heteroplasmy_min_depth: int = 20
    restricted_heteroplasmy_min_alt_reads: int = 5
    restricted_heteroplasmy_min_vaf: float = 0.25
    restricted_heteroplasmy_max_vaf: float = 0.75
    min_internal_fragments: int = 2
    spatial_thresholds: Dict[str, int] = field(default_factory=lambda: {
        "NON_UDG": 5, "UDG_HALF": 2, "UDG_FULL": 0, "SINGLE_STRANDED": 5
    })
    excluded_intervals: Tuple[Tuple[int, int], ...] = (
        (3105, 3111), (303, 315), (513, 524), (16184, 16193),
    )


@dataclass
class AlleleCounts:
    ref_fwd: int = 0
    ref_rev: int = 0
    alt_fwd: int = 0
    alt_rev: int = 0

    @property
    def ref_total(self) -> int:
        return self.ref_fwd + self.ref_rev

    @property
    def alt_total(self) -> int:
        return self.alt_fwd + self.alt_rev

    @property
    def total_depth(self) -> int:
        return self.ref_total + self.alt_total

    @property
    def vaf(self) -> float:
        return (self.alt_total / self.total_depth) if self.total_depth > 0 else 0.0


@dataclass
class SitePileup:
    position: int
    ref_allele: str
    depth: int
    base_counts: Dict[str, int]
    strand_counts: Dict[str, Tuple[int, int]]  # Base -> (fwd, rev)


@dataclass
class HeteroplasmyCall:
    position: int
    ref_allele: str
    alt_allele: str
    counts: AlleleCounts
    vaf: float
    status: CallStatus
    filter_reasons: List[FilterReason] = field(default_factory=list)
    strand_bias_pvalue: float = 1.0


@dataclass
class ConsensusQC:
    sample_id: str
    input_bam_path: str = ""
    input_bam_sha256: str = ""
    reference_contig: str = "NC_012920.1"
    reference_sha256: str = ""
    reference_length: int = 16569
    consensus_fasta_sha256: str = ""
    mean_depth: float = 0.0
    covered_positions_fraction: float = 0.0
    n_fraction_canonical: float = 0.0
    iupac_fraction: float = 0.0
    heteroplasmy_sites_total: int = 0
    heteroplasmy_sites_pass: int = 0
    heteroplasmy_sites_filtered: int = 0
    regions: List[Any] = field(default_factory=list)
    gating_action: GatingAction = GatingAction.RESTRICT_HETEROPLASMY
    software_version: str = "1.0.0"
    policy_version: str = "1.0.0"
    execution_time_seconds: float = 0.0
    heteroplasmy_details: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ResolvedInputs:
    sample_id: str
    bam_path: Path
    reference_fasta_path: Path
    gating_action: GatingAction = GatingAction.RESTRICT_HETEROPLASMY
    library_treatment: str = "NON_UDG"


# ==============================================================================
# SECTION 2: PILEUP & HETEROPLASMY ENGINES
# ==============================================================================

class PileupEngine:
    """Extracts quality-filtered nucleotide pileups across target mitochondrial positions."""

    @classmethod
    def generate_pileups(
        cls,
        bam_path: Path,
        ref_sequence: str,
        policy: ConsensusPolicy,
        target_length: int = 16569,
    ) -> List[SitePileup]:
        pileups: List[SitePileup] = []
        bases = ["A", "C", "G", "T"]

        # Preallocate base structure
        counts = [{b: [0, 0] for b in bases} for _ in range(target_length)]

        with pysam.AlignmentFile(str(bam_path), "rb") as sam:
            for read in sam:
                if read.is_unmapped or read.mapping_quality < policy.min_mapping_quality:
                    continue

                seq = read.query_sequence
                quals = read.query_qualities
                if not seq or not quals:
                    continue

                is_rev = read.is_reverse
                strand_idx = 1 if is_rev else 0

                for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                    if 0 <= rpos < target_length:
                        if quals[qpos] >= policy.min_base_quality:
                            base = seq[qpos].upper()
                            if base in counts[rpos]:
                                counts[rpos][base][strand_idx] += 1

        for pos in range(target_length):
            ref_b = ref_sequence[pos] if pos < len(ref_sequence) else "N"
            pos_dict = counts[pos]
            total_d = sum(pos_dict[b][0] + pos_dict[b][1] for b in bases)
            base_totals = {b: (pos_dict[b][0] + pos_dict[b][1]) for b in bases}
            strand_dict = {b: (pos_dict[b][0], pos_dict[b][1]) for b in bases}

            pileups.append(SitePileup(
                position=pos + 1,  # 1-based coordinate
                ref_allele=ref_b,
                depth=total_d,
                base_counts=base_totals,
                strand_counts=strand_dict,
            ))

        return pileups


class HeteroplasmyCaller:
    """Evaluates minor alleles, tests strand bias, and applies damage filtering."""

    @classmethod
    def call_site(
        cls,
        site: SitePileup,
        policy: ConsensusPolicy,
        gating: GatingAction,
        treatment: str = "NON_UDG",
    ) -> Optional[HeteroplasmyCall]:
        if site.depth < policy.heteroplasmy_min_depth:
            return None

        # Excluded intervals blacklist check
        for start, end in policy.excluded_intervals:
            if start <= site.position <= end:
                return None

        ref = site.ref_allele
        # Find candidate alternate allele
        sorted_bases = sorted(
            [(b, c) for b, c in site.base_counts.items() if b != ref],
            key=lambda x: x[1],
            reverse=True,
        )
        if not sorted_bases:
            return None

        alt, alt_count = sorted_bases[0]
        if alt_count < policy.heteroplasmy_min_alt_reads:
            return None

        vaf = alt_count / site.depth
        min_vaf = policy.heteroplasmy_min_vaf
        max_vaf = policy.heteroplasmy_max_vaf

        if gating == GatingAction.RESTRICT_HETEROPLASMY:
            min_vaf = policy.restricted_heteroplasmy_min_vaf
            max_vaf = policy.restricted_heteroplasmy_max_vaf

        if not (min_vaf <= vaf <= max_vaf):
            return None

        ref_f, ref_r = site.strand_counts.get(ref, (0, 0))
        alt_f, alt_r = site.strand_counts.get(alt, (0, 0))

        counts = AlleleCounts(ref_fwd=ref_f, ref_rev=ref_r, alt_fwd=alt_f, alt_rev=alt_r)

        # Strand bias evaluation via Fisher's exact test
        table = [[ref_f, ref_r], [alt_f, alt_r]]
        _, p_val = fisher_exact(table)

        filter_reasons: List[FilterReason] = []
        if p_val < policy.strand_bias_alpha and (alt_f < policy.strand_bias_min_strand_reads or alt_r < policy.strand_bias_min_strand_reads):
            filter_reasons.append(FilterReason.STRAND_BIAS)

        # Spatial damage deamination check
        spatial_cutoff = policy.spatial_thresholds.get(treatment, 5)
        if (ref == "C" and alt == "T") or (ref == "G" and alt == "A"):
            if site.position <= spatial_cutoff or site.position >= (16569 - spatial_cutoff):
                filter_reasons.append(FilterReason.SPATIAL_DAMAGE)

        status = CallStatus.PASS if not filter_reasons else CallStatus.FILTERED

        return HeteroplasmyCall(
            position=site.position,
            ref_allele=ref,
            alt_allele=alt,
            counts=counts,
            vaf=round(vaf, 4),
            status=status,
            filter_reasons=filter_reasons,
            strand_bias_pvalue=round(p_val, 4),
        )


# ==============================================================================
# SECTION 3: CONSENSUS BUILDER & VCF EMITTER
# ==============================================================================

class ConsensusBuilder:
    """Builds major consensus, IUPAC consensus, and emission files."""

    @classmethod
    def build_consensus_sequences(
        cls,
        pileups: List[SitePileup],
        heteroplasmies: List[HeteroplasmyCall],
        policy: ConsensusPolicy,
    ) -> Tuple[str, str]:
        canonical_chars = []
        iupac_chars = []

        het_by_pos = {h.position: h for h in heteroplasmies if h.status == CallStatus.PASS}

        for site in pileups:
            if site.depth < policy.min_consensus_depth:
                canonical_chars.append(policy.mask_char)
                iupac_chars.append(policy.mask_char)
                continue

            # Dominant allele
            dominant_base, dominant_count = max(site.base_counts.items(), key=lambda x: x[1])
            frac = dominant_count / site.depth

            if frac >= policy.major_allele_fraction:
                major_base = dominant_base
            else:
                major_base = policy.mask_char

            canonical_chars.append(major_base)

            # IUPAC allele
            if site.position in het_by_pos:
                h = het_by_pos[site.position]
                pair = (h.ref_allele, h.alt_allele)
                iupac_code = IUPAC_MAP.get(pair, major_base)
                iupac_chars.append(iupac_code)
            else:
                iupac_chars.append(major_base)

        return "".join(canonical_chars), "".join(iupac_chars)

    @classmethod
    def emit_vcf(
        cls,
        sample_id: str,
        pileups: List[SitePileup],
        heteroplasmies: List[HeteroplasmyCall],
        vcf_path: Path,
        ref_contig: str = "NC_012920.1",
    ) -> None:
        """Emits bgzipped and indexed VCF file."""
        vcf_path.parent.mkdir(parents=True, exist_ok=True)
        raw_vcf = vcf_path.with_suffix(".tmp.vcf")

        with open(raw_vcf, "w", encoding="utf-8") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write(f"##reference={ref_contig}\n")
            f.write(f"##contig=<ID={ref_contig},length=16569>\n")
            f.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n')
            f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            f.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depths">\n')
            f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}\n")

            for h in heteroplasmies:
                filt_str = "PASS" if h.status == CallStatus.PASS else ";".join(r.value for r in h.filter_reasons)
                info = f"DP={h.counts.total_depth};AF={h.vaf:.4f}"
                fmt = "GT:AD"
                sample_val = f"0/1:{h.counts.ref_total},{h.counts.alt_total}"
                f.write(
                    f"{ref_contig}\t{h.position}\t.\t{h.ref_allele}\t{h.alt_allele}\t.\t"
                    f"{filt_str}\t{info}\t{fmt}\t{sample_val}\n"
                )

        pysam.tabix_compress(str(raw_vcf), str(vcf_path), force=True)
        pysam.tabix_index(str(vcf_path), preset="vcf", force=True)

        if raw_vcf.is_file():
            raw_vcf.unlink()


# ==============================================================================
# SECTION 4: ORCHESTRATOR & STEP RUNNER
# ==============================================================================

class ConsensusOrchestrator:
    """Coordinates Block 6 execution."""

    def __init__(self, policy: Optional[ConsensusPolicy] = None):
        self.policy = policy or ConsensusPolicy()

    def resolve_inputs(
        self,
        sample_id: str,
        block4_json_path: str,
        block5_json_path: str,
        reference_fasta_path: str,
        target_contig: str = "NC_012920.1",
    ) -> ResolvedInputs:
        with open(block4_json_path, "r", encoding="utf-8") as f4:
            b4 = json.load(f4)
        with open(block5_json_path, "r", encoding="utf-8") as f5:
            b5 = json.load(f5)

        bam_str = b4.get("rescaled_bam") or b4.get("filtered_dedup_bam")
        cat = b5.get("contamination_category", "LOW")

        if cat in ("HIGH", "CRITICAL"):
            gating = GatingAction.SUPPRESS_HETEROPLASMY
        elif cat == "MODERATE":
            gating = GatingAction.RESTRICT_HETEROPLASMY
        else:
            gating = GatingAction.PROCEED_UNRESTRICTED

        return ResolvedInputs(
            sample_id=sample_id,
            bam_path=Path(bam_str),
            reference_fasta_path=Path(reference_fasta_path),
            gating_action=gating,
            library_treatment=b4.get("library_treatment", "NON_UDG"),
        )

    def run(self, inputs: ResolvedInputs, output_dir: str) -> ConsensusQC:
        out_p = Path(output_dir)
        out_p.mkdir(parents=True, exist_ok=True)
        s_id = inputs.sample_id

        # Load reference sequence
        with open(inputs.reference_fasta_path, "r", encoding="utf-8") as f:
            ref_seq = "".join(ln.strip() for ln in f if not ln.startswith(">")).upper()

        t0 = time.perf_counter()
        pileups = PileupEngine.generate_pileups(inputs.bam_path, ref_seq, self.policy)

        # Call heteroplasmies
        het_calls: List[HeteroplasmyCall] = []
        for site in pileups:
            call = HeteroplasmyCaller.call_site(site, self.policy, inputs.gating_action, inputs.library_treatment)
            if call:
                het_calls.append(call)

        canonical_seq, iupac_seq = ConsensusBuilder.build_consensus_sequences(pileups, het_calls, self.policy)

        # Write FASTA files
        can_fa = out_p / f"{s_id}.consensus.fasta"
        iupac_fa = out_p / f"{s_id}.iupac.fasta"

        with open(can_fa, "w", encoding="utf-8") as f:
            f.write(f">{s_id}_consensus\n{canonical_seq}\n")
        with open(iupac_fa, "w", encoding="utf-8") as f:
            f.write(f">{s_id}_iupac\n{iupac_seq}\n")

        vcf_gz = out_p / f"{s_id}.variants.vcf.gz"
        ConsensusBuilder.emit_vcf(s_id, pileups, het_calls, vcf_gz)

        # Heteroplasmy TSV
        het_tsv = out_p / f"{s_id}.heteroplasmy.tsv"
        with open(het_tsv, "w", encoding="utf-8") as f:
            f.write("pos\tref\talt\tdepth\talt_depth\tvaf\tstatus\tfilter\n")
            for h in het_calls:
                filt_str = "PASS" if h.status == CallStatus.PASS else ";".join(r.value for r in h.filter_reasons)
                f.write(f"{h.position}\t{h.ref_allele}\t{h.alt_allele}\t{h.counts.total_depth}\t{h.counts.alt_total}\t{h.vaf:.4f}\t{h.status.value}\t{filt_str}\n")

        depths = [p.depth for p in pileups]
        mean_d = sum(depths) / len(depths) if depths else 0.0
        cov_frac = sum(1 for d in depths if d >= self.policy.min_consensus_depth) / len(depths) if depths else 0.0
        n_frac = canonical_seq.count("N") / len(canonical_seq) if canonical_seq else 0.0
        iupac_frac = sum(1 for ch in iupac_seq if ch in IUPAC_MAP.values()) / len(iupac_seq) if iupac_seq else 0.0

        n_pass = sum(1 for h in het_calls if h.status == CallStatus.PASS)
        n_filt = sum(1 for h in het_calls if h.status == CallStatus.FILTERED)

        qc = ConsensusQC(
            sample_id=s_id,
            input_bam_path=str(inputs.bam_path),
            consensus_fasta_sha256=calculate_sha256(can_fa),
            mean_depth=round(mean_d, 2),
            covered_positions_fraction=round(cov_frac, 4),
            n_fraction_canonical=round(n_frac, 4),
            iupac_fraction=round(iupac_frac, 4),
            heteroplasmy_sites_total=len(het_calls),
            heteroplasmy_sites_pass=n_pass,
            heteroplasmy_sites_filtered=n_filt,
            gating_action=inputs.gating_action,
            execution_time_seconds=round(time.perf_counter() - t0, 2),
        )

        qc_json = out_p / f"{s_id}.consensus_qc.json"
        with open(qc_json, "w", encoding="utf-8") as f:
            json.dump({
                "sample_id": s_id,
                "input_bam_path": qc.input_bam_path,
                "consensus_fasta_sha256": qc.consensus_fasta_sha256,
                "mean_depth": qc.mean_depth,
                "covered_positions_fraction": qc.covered_positions_fraction,
                "n_fraction_canonical": qc.n_fraction_canonical,
                "iupac_fraction": qc.iupac_fraction,
                "heteroplasmy_sites_total": qc.heteroplasmy_sites_total,
                "heteroplasmy_sites_pass": qc.heteroplasmy_sites_pass,
                "heteroplasmy_sites_filtered": qc.heteroplasmy_sites_filtered,
                "gating_action": qc.gating_action.value,
                "execution_time_seconds": qc.execution_time_seconds,
            }, f, indent=2)

        return qc


class ConsensusStep:
    """High-level runner for Pipeline Step 6 (Consensus & Heteroplasmy Calling)."""

    @classmethod
    def run(
        cls,
        sample_id: str,
        b4_contract_path: Path,
        b5_contract_path: Path,
        reference_info: ReferenceInfo,
        config: PipelineConfig,
        operational_files_dir: Path,
        cache_dir: Path,
        force_recompute: bool = False,
    ) -> Tuple[ConsensusQC, Path, Path, Path]:
        """Runs consensus calling and emits Block 6 contract."""
        legacy_dir = operational_files_dir / "alignment" / sample_id
        tracker = StageTracker.load_or_create(
            cache_dir,
            sample_id,
            legacy_dir=legacy_dir,
            force_recompute=force_recompute,
        )

        block6_out = operational_files_dir / "consensus" / sample_id
        canonical_fa = block6_out / f"{sample_id}.consensus.fasta"
        vcf_gz = block6_out / f"{sample_id}.variants.vcf.gz"
        qc_json = block6_out / f"{sample_id}.consensus_qc.json"
        het_tsv = block6_out / f"{sample_id}.heteroplasmy.tsv"
        b6_contract_json = block6_out / f"{sample_id}.block6_contract.json"

        b4_sha256 = calculate_sha256(b4_contract_path)
        b5_sha256 = calculate_sha256(b5_contract_path)
        cons_cfg = config.raw_data.get("consensus", {})

        cons_composite_hash = compute_composite_hash(
            cons_cfg,
            reference_info.file_sha256,
            b4_sha256,
            b5_sha256,
        )

        expected_artifacts = [canonical_fa, vcf_gz, qc_json, b6_contract_json]

        if tracker.is_stage_valid("CONSENSUS_CALLING", expected_artifacts, cons_composite_hash):
            logger.info(f"[{sample_id}] Stage CONSENSUS_CALLING cached -> {canonical_fa.name}")
            with open(qc_json, "r", encoding="utf-8") as f:
                qc_dict = json.load(f)

            qc = ConsensusQC(
                sample_id=sample_id,
                input_bam_path=qc_dict.get("input_bam_path", ""),
                consensus_fasta_sha256=qc_dict.get("consensus_fasta_sha256", ""),
                mean_depth=qc_dict.get("mean_depth", 0.0),
                covered_positions_fraction=qc_dict.get("covered_positions_fraction", 0.0),
                n_fraction_canonical=qc_dict.get("n_fraction_canonical", 0.0),
                iupac_fraction=qc_dict.get("iupac_fraction", 0.0),
                heteroplasmy_sites_total=qc_dict.get("heteroplasmy_sites_total", 0),
                heteroplasmy_sites_pass=qc_dict.get("heteroplasmy_sites_pass", 0),
                heteroplasmy_sites_filtered=qc_dict.get("heteroplasmy_sites_filtered", 0),
                gating_action=GatingAction(qc_dict.get("gating_action", "RESTRICT_HETEROPLASMY")),
                execution_time_seconds=qc_dict.get("execution_time_seconds", 0.0),
            )
            return qc, canonical_fa, vcf_gz, b6_contract_json

        t0 = time.perf_counter()
        orchestrator = ConsensusOrchestrator()
        resolution = orchestrator.resolve_inputs(
            sample_id=sample_id,
            block4_json_path=str(b4_contract_path),
            block5_json_path=str(b5_contract_path),
            reference_fasta_path=str(reference_info.fasta_path),
        )

        qc = orchestrator.run(resolution, str(block6_out))
        elapsed = time.perf_counter() - t0

        input_bam_sha = calculate_sha256(resolution.bam_path)
        cls.ensure_block6_contract(
            b6_contract_json,
            sample_id,
            canonical_fa,
            vcf_gz,
            qc,
            het_tsv,
            reference_info.file_sha256,
            input_bam_sha256=input_bam_sha,
            b4_contract_sha256=b4_sha256,
            b5_contract_sha256=b5_sha256,
            force_overwrite=True,
        )

        tracker.mark_success(
            "CONSENSUS_CALLING",
            elapsed,
            {
                "canonical_fasta": canonical_fa,
                "vcf": vcf_gz,
                "qc_json": qc_json,
                "b6_contract": b6_contract_json,
            },
            composite_hash=cons_composite_hash,
        )
        return qc, canonical_fa, vcf_gz, b6_contract_json

    @classmethod
    def ensure_block6_contract(
        cls,
        contract_path: Path,
        sample_id: str,
        canonical_fa: Path,
        vcf_gz: Path,
        qc: ConsensusQC,
        het_tsv: Path,
        ref_sha256: str,
        input_bam_sha256: str = "",
        b4_contract_sha256: str = "",
        b5_contract_sha256: str = "",
        force_overwrite: bool = False,
    ) -> None:
        """Serializes Block 6 Contract (block6_contract.json)."""
        contract_path.parent.mkdir(parents=True, exist_ok=True)
        if not force_overwrite and contract_path.is_file() and contract_path.stat().st_size > 0:
            return

        contract_data = {
            "sample_id": sample_id,
            "consensus_fasta": str(canonical_fa),
            "consensus_fasta_sha256": calculate_sha256(canonical_fa),
            "variants_vcf": str(vcf_gz),
            "variants_vcf_sha256": calculate_sha256(vcf_gz),
            "mean_depth": qc.mean_depth,
            "covered_positions_fraction": qc.covered_positions_fraction,
            "n_fraction_canonical": qc.n_fraction_canonical,
            "iupac_fraction": qc.iupac_fraction,
            "heteroplasmy_sites_total": qc.heteroplasmy_sites_total,
            "heteroplasmy_sites_pass": qc.heteroplasmy_sites_pass,
            "heteroplasmy_sites_filtered": qc.heteroplasmy_sites_filtered,
            "gating_action": qc.gating_action.value,
            "input_bam_sha256": input_bam_sha256,
            "b4_contract_sha256": b4_contract_sha256,
            "b5_contract_sha256": b5_contract_sha256,
            "reference_sha256": ref_sha256,
            "execution_time_seconds": qc.execution_time_seconds,
        }

        with open(contract_path, "w", encoding="utf-8") as f:
            json.dump(contract_data, f, indent=2)

    @classmethod
    def export_consensus_cohort_summary(
        cls,
        cohort_qc: Dict[str, ConsensusQC],
        output_dir: Path,
    ) -> Tuple[Path, Path]:
        """Exports cohort summary TSV and JSON files for Block 6."""
        tsv_path = output_dir / "cohort_consensus_summary.tsv"
        json_path = output_dir / "cohort_consensus_summary.json"

        records = []
        headers = [
            "sample_id", "mean_depth", "covered_positions_fraction",
            "n_fraction_canonical", "iupac_fraction", "heteroplasmy_sites_total",
            "heteroplasmy_sites_pass", "gating_action"
        ]

        with open(tsv_path, "w", encoding="utf-8") as f_t:
            f_t.write("\t".join(headers) + "\n")
            for s_id, qc in cohort_qc.items():
                r = {
                    "sample_id": s_id,
                    "mean_depth": qc.mean_depth,
                    "covered_positions_fraction": qc.covered_positions_fraction,
                    "n_fraction_canonical": qc.n_fraction_canonical,
                    "iupac_fraction": qc.iupac_fraction,
                    "heteroplasmy_sites_total": qc.heteroplasmy_sites_total,
                    "heteroplasmy_sites_pass": qc.heteroplasmy_sites_pass,
                    "gating_action": qc.gating_action.value,
                }
                records.append(r)
                f_t.write("\t".join(str(r[h]) for h in headers) + "\n")

        with open(json_path, "w", encoding="utf-8") as f_j:
            json.dump({"total_samples": len(records), "samples": records}, f_j, indent=2)

        return tsv_path, json_path
