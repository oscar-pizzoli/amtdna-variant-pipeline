"""
Block 7 Cross-Block Auditor: Multi-block consistency integration, topological
phylotree graph traversal, reference-aware mutation parsing, and spatial damage triage.
"""

import logging
import re
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Set, Tuple

from src.modules.haplogrouping.enums import (
    AuditCode,
    AuditSeverity,
    ContaminantAffinity,
    OverallValidationStatus,
    PrivateMutationDamageClass,
)
from src.modules.haplogrouping.models import (
    AuditFlag,
    CrossBlockAuditResult,
    PrivateMutationAssessment,
    RawHaplogrepResult,
)

logger = logging.getLogger("mtDNA_Pipeline.CrossBlockAuditor")


@dataclass(frozen=True)
class ParsedMutation:
    """Structured decomposition of a mitochondrial variant."""
    raw: str
    variant_type: str  # "SNV", "INS", "DEL", "UNKNOWN"
    position: int
    ref_base: str
    alt_base: str


class PhylotreeGraph:
    """
    Topological and phylogenetic clade resolution engine.
    Traverses true ancestral paths and resolves macro-haplogroup aliases
    without string-prefix fallacies.
    """

    KNOWN_MACRO_PARENTS: Dict[str, str] = {
        "H": "HV",
        "V": "HV0",
        "HV0": "HV",
        "HV": "R0",
        "R0": "R",
        "J": "JT",
        "T": "JT",
        "JT": "R",
        "K": "U8B",
        "U8B": "U8",
        "U8": "U",
        "U": "R",
        "B": "R",
        "F": "R",
        "P": "R",
        "R": "N",
        "I": "N1A1",
        "N1A1": "N1A",
        "N1A": "N1",
        "N1": "N",
        "W": "N2",
        "N2": "N",
        "X": "N",
        "Y": "N9",
        "N9": "N",
        "A": "N",
        "S": "N",
        "O": "N",
        "N": "L3",
        "C": "M8",
        "Z": "M8",
        "CZ": "M8",
        "M8": "M",
        "D": "M",
        "G": "M",
        "E": "M9",
        "M9": "M",
        "Q": "M",
        "M": "L3",
        "L0": "L",
        "L1": "L",
        "L2": "L",
        "L3": "L",
        "L4": "L",
        "L5": "L",
        "L6": "L",
    }

    @classmethod
    def get_parent(cls, hg: str) -> Optional[str]:
        clean = hg.strip().upper()
        clean = re.sub(r"[\+@'].*$", "", clean)

        if not clean or clean == "L":
            return None

        if clean in cls.KNOWN_MACRO_PARENTS:
            return cls.KNOWN_MACRO_PARENTS[clean]

        digit_match = re.search(r"^(.*?)\d+$", clean)
        if digit_match:
            prefix = digit_match.group(1)
            if prefix:
                return prefix

        if re.search(r"\d", clean):
            letter_match = re.search(r"^(.*?\d)[A-Z]+$", clean)
            if letter_match:
                return letter_match.group(1)

        return None

    @classmethod
    def is_descendant(cls, candidate_child: str, candidate_parent: str) -> bool:
        c_child = candidate_child.strip().upper()
        c_parent = candidate_parent.strip().upper()

        if not c_child or not c_parent:
            return False
        if c_child == c_parent:
            return True

        curr: Optional[str] = c_child
        visited: Set[str] = set()
        max_depth = 60
        depth = 0

        while curr is not None and depth < max_depth:
            if curr in visited:
                break
            visited.add(curr)
            parent = cls.get_parent(curr)
            if parent is not None and parent == c_parent:
                return True
            curr = parent
            depth += 1

        return False


class CrossBlockAuditor:
    """Evaluates phylogenetic assignment against upstream authenticity, contamination, and consensus contracts."""

    TRANSITIONS: Set[Tuple[str, str]] = {
        ("C", "T"),
        ("T", "C"),
        ("A", "G"),
        ("G", "A"),
    }
    DAMAGE_TRANSITIONS: Set[Tuple[str, str]] = {
        ("C", "T"),
        ("G", "A"),
    }

    KNOWN_GATING_ACTIONS: Set[str] = {
        "STANDARD_CALLING",
        "SENSITIVE_CALLING",
        "HIGH_CONTAMINATION_CAUTION",
        "CONTAMINATED_REJECT",
        "FAST_CLEAR",
    }

    INS_REGEX = re.compile(r"^([0-9]+)\.([0-9]+)([ACGTN]+)$", re.IGNORECASE)
    DEL_REGEX = re.compile(r"^([ACGTN]?|>)?([0-9]+)(?:-[0-9]+)?(?:D|DEL)$", re.IGNORECASE)
    SNV_REGEX = re.compile(r"^([ACGTN]?|>)?([0-9]+)(?:>)?([ACGTN])$", re.IGNORECASE)

    def __init__(
        self,
        min_quality_score: float = 0.80,
        max_missing_diagnostic_mutations: int = 2,
        max_private_mutations: int = 5,
        contam_affinity_threshold: float = 0.015,
        subclade_contam_max_fraction: float = 0.03,
        evaluate_contaminant_closeness: bool = True,
        flag_damage_in_private_mutations: bool = True,
    ):
        self.min_quality_score = min_quality_score
        self.max_missing_diagnostic_mutations = max_missing_diagnostic_mutations
        self.max_private_mutations = max_private_mutations
        self.contam_affinity_threshold = contam_affinity_threshold
        self.subclade_contam_max_fraction = subclade_contam_max_fraction
        self.evaluate_contaminant_closeness = evaluate_contaminant_closeness
        self.flag_damage_in_private_mutations = flag_damage_in_private_mutations

    def audit(
        self,
        sample_id: str,
        fasta_result: RawHaplogrepResult,
        vcf_result: Optional[RawHaplogrepResult] = None,
        block4_contract: Optional[Dict[str, Any]] = None,
        block5_contract: Optional[Dict[str, Any]] = None,
        block6_contract: Optional[Dict[str, Any]] = None,
        reference_sequence: Optional[str] = None,
    ) -> CrossBlockAuditResult:
        flags: List[AuditFlag] = []

        # 0. Strict Upstream Contract Verification
        if not block4_contract:
            flags.append(
                AuditFlag(
                    code=AuditCode.MISSING_UPSTREAM_CONTRACT,
                    severity=AuditSeverity.CRITICAL,
                    message=f"Block 4 contract is missing or empty for '{sample_id}'. Provenance broken.",
                    context={"missing_contract": "Block 4"},
                )
            )
        if not block5_contract:
            flags.append(
                AuditFlag(
                    code=AuditCode.MISSING_UPSTREAM_CONTRACT,
                    severity=AuditSeverity.CRITICAL,
                    message=f"Block 5 contract is missing or empty for '{sample_id}'. Provenance broken.",
                    context={"missing_contract": "Block 5"},
                )
            )

        b4_data = block4_contract or {}
        b5_data = block5_contract or {}
        b6_data = block6_contract or {}

        # 1. Quality Score & Ambiguity Auditing
        if fasta_result.quality_score < self.min_quality_score:
            flags.append(
                AuditFlag(
                    code=AuditCode.LOW_HAPLOGROUP_CONFIDENCE,
                    severity=AuditSeverity.WARNING,
                    message=f"Quality score {fasta_result.quality_score:.4f} is below threshold {self.min_quality_score:.4f}.",
                    context={"quality_score": fasta_result.quality_score, "threshold": self.min_quality_score},
                )
            )

        if fasta_result.co_optimal:
            co_hits = [h.haplogroup for h in fasta_result.top_hits]
            flags.append(
                AuditFlag(
                    code=AuditCode.MULTIPLE_TOP_HITS_AMBIGUOUS,
                    severity=AuditSeverity.WARNING,
                    message=f"Multiple co-optimal candidate haplogroups identified: {co_hits}.",
                    context={"co_optimal_haplogroups": co_hits},
                )
            )

        # 2. Mutation Counts
        missing_count = len(fasta_result.missing_mutations)
        if missing_count > self.max_missing_diagnostic_mutations:
            flags.append(
                AuditFlag(
                    code=AuditCode.MISSING_DIAGNOSTIC_MUTATIONS,
                    severity=AuditSeverity.WARNING,
                    message=f"Sample missing {missing_count} diagnostic mutations for node '{fasta_result.haplogroup}'.",
                    context={
                        "missing_count": missing_count,
                        "limit": self.max_missing_diagnostic_mutations,
                        "missing_mutations": fasta_result.missing_mutations,
                    },
                )
            )

        private_count = len(fasta_result.private_mutations)
        if private_count > self.max_private_mutations:
            flags.append(
                AuditFlag(
                    code=AuditCode.EXCESS_PRIVATE_MUTATIONS,
                    severity=AuditSeverity.WARNING,
                    message=f"Sample carries {private_count} private mutations (limit: {self.max_private_mutations}).",
                    context={
                        "private_count": private_count,
                        "limit": self.max_private_mutations,
                        "private_mutations": fasta_result.private_mutations,
                    },
                )
            )

        # 3. Upstream Contamination Contract Verification (Bypass Resolution)
        raw_contam = b5_data.get("contamination_fraction")
        contam_frac: Optional[float] = None
        if raw_contam is None:
            flags.append(
                AuditFlag(
                    code=AuditCode.MISSING_CONTAMINATION_METRIC,
                    severity=AuditSeverity.CRITICAL,
                    message=f"Missing contamination_fraction in Block 5 contract for '{sample_id}'. Risk treated as maximum.",
                    context={"sample_id": sample_id},
                )
            )
        else:
            try:
                contam_frac = float(raw_contam)
            except (ValueError, TypeError):
                flags.append(
                    AuditFlag(
                        code=AuditCode.UPSTREAM_CONTRACT_ANOMALY,
                        severity=AuditSeverity.CRITICAL,
                        message=f"Non-numeric contamination_fraction in Block 5: {raw_contam}",
                        context={"raw_contamination_fraction": raw_contam},
                    )
                )

        gating_action = b5_data.get("recommended_action")
        if gating_action is not None and str(gating_action).upper() not in self.KNOWN_GATING_ACTIONS:
            flags.append(
                AuditFlag(
                    code=AuditCode.UPSTREAM_CONTRACT_ANOMALY,
                    severity=AuditSeverity.CRITICAL,
                    message=f"Unrecognized gating recommended_action in Block 5: {gating_action}",
                    context={"recommended_action": gating_action},
                )
            )

        # 4. Topological Dual-Stream Concordance with Proportional Subclade Refinement
        concordant, resolved_hg, subclade_detected, subclade_supported = (
            self._evaluate_dual_stream_concordance(
                sample_id=sample_id,
                fasta_result=fasta_result,
                vcf_result=vcf_result,
                b4_data=b4_data,
                b6_data=b6_data,
                contam_fraction=contam_frac,
                reference_sequence=reference_sequence,
                flags=flags,
            )
        )

        # 5. Topological Contaminant Affinity
        affinity, contam_hgs = self._evaluate_contaminant_affinity(
            endogenous_hg=fasta_result.haplogroup,
            b5_data=b5_data,
        )

        if self.evaluate_contaminant_closeness and affinity in {
            ContaminantAffinity.EXACT_MATCH,
            ContaminantAffinity.SUBCLADE,
            ContaminantAffinity.ANCESTRAL,
        }:
            is_elevated_contam = (contam_frac is None) or (contam_frac >= self.contam_affinity_threshold)
            gating_str = str(gating_action).upper() if gating_action else ""
            is_non_standard = (gating_str not in {"", "STANDARD_CALLING", "FAST_CLEAR"})

            if is_elevated_contam or is_non_standard:
                frac_str = f"{contam_frac:.2%}" if contam_frac is not None else "UNKNOWN (MAX RISK)"
                flags.append(
                    AuditFlag(
                        code=AuditCode.CONTAMINANT_AFFINITY_WARNING,
                        severity=AuditSeverity.WARNING,
                        message=(
                            f"Endogenous haplogroup '{fasta_result.haplogroup}' shows topological "
                            f"affinity ({affinity.value}) to contaminant candidate '{contam_hgs}' "
                            f"under contamination fraction {frac_str}."
                        ),
                        context={
                            "endogenous_haplogroup": fasta_result.haplogroup,
                            "contaminant_haplogroups": contam_hgs,
                            "affinity": affinity.value,
                            "contamination_fraction": contam_frac,
                        },
                    )
                )

        # 6. Reference-Aware Chemical Screening of Private Mutations
        private_assessments = self._audit_private_mutations_chemistry(
            private_mutations=fasta_result.private_mutations,
            b4_data=b4_data,
            reference_sequence=reference_sequence,
        )

        spatial_audit_muts = [
            p.mutation
            for p in private_assessments
            if p.classification == PrivateMutationDamageClass.TRANSITION_REQUIRES_SPATIAL_AUDIT
        ]
        if spatial_audit_muts and self.flag_damage_in_private_mutations:
            flags.append(
                AuditFlag(
                    code=AuditCode.DAMAGE_SUSPECT_PRIVATE_MUTATION,
                    severity=AuditSeverity.WARNING,
                    message=(
                        f"Sample carries {len(spatial_audit_muts)} private transition(s) with deamination chemistry "
                        f"requiring spatial read-end audit in BAM: {spatial_audit_muts}."
                    ),
                    context={
                        "mutations_requiring_spatial_audit": spatial_audit_muts,
                        "count": len(spatial_audit_muts),
                    },
                )
            )

        overall_status = self._derive_overall_status(flags)

        return CrossBlockAuditResult(
            overall_status=overall_status,
            fasta_vcf_concordant=concordant,
            fasta_haplogroup=fasta_result.haplogroup,
            vcf_haplogroup=vcf_result.haplogroup if vcf_result else None,
            resolved_haplogroup=resolved_hg,
            quality_score=fasta_result.quality_score,
            quality_margin=fasta_result.quality_margin,
            contaminant_affinity=affinity,
            contaminant_haplogroups=contam_hgs,
            contamination_fraction=contam_frac,
            subclade_shift_detected=subclade_detected,
            subclade_shift_supported=subclade_supported,
            private_mutations_assessment=private_assessments,
            flags=flags,
        )

    def _derive_overall_status(self, flags: List[AuditFlag]) -> OverallValidationStatus:
        if any(f.severity == AuditSeverity.CRITICAL for f in flags):
            return OverallValidationStatus.FAIL
        if any(f.severity == AuditSeverity.WARNING for f in flags):
            return OverallValidationStatus.WARN
        return OverallValidationStatus.PASS

    def _evaluate_dual_stream_concordance(
        self,
        sample_id: str,
        fasta_result: RawHaplogrepResult,
        vcf_result: Optional[RawHaplogrepResult],
        b4_data: Dict[str, Any],
        b6_data: Dict[str, Any],
        contam_fraction: Optional[float],
        reference_sequence: Optional[str],
        flags: List[AuditFlag],
    ) -> Tuple[bool, str, bool, Optional[bool]]:
        f_hg = fasta_result.haplogroup.strip().upper()
        if vcf_result is None:
            return True, f_hg, False, None

        v_hg = vcf_result.haplogroup.strip().upper()
        if f_hg == v_hg:
            return True, v_hg, False, None

        # Scenario A: VCF topologically descends from FASTA
        if PhylotreeGraph.is_descendant(candidate_child=v_hg, candidate_parent=f_hg):
            subclade_detected = True
            flags.append(
                AuditFlag(
                    code=AuditCode.SUBCLADE_HETEROPLASMY_REFINEMENT,
                    severity=AuditSeverity.INFO,
                    message=f"VCF classification refined primary node '{f_hg}' to descendant '{v_hg}'.",
                    context={"primary_fasta": f_hg, "vcf_subclade": v_hg},
                )
            )

            diff_muts = set(vcf_result.found_mutations) - set(fasta_result.found_mutations)
            diff_positions = self._extract_positions(diff_muts)

            raw_b6_sites = b6_data.get("heteroplasmy_sites_pass") or []
            b6_pass_sites: Set[int] = set()
            for s in raw_b6_sites:
                try:
                    b6_pass_sites.add(int(s))
                except (ValueError, TypeError):
                    pass

            matching_sites = diff_positions.intersection(b6_pass_sites)
            auth_status = str(b4_data.get("authentication_status") or "").upper()

            # Fail-safe: None means maximum risk
            is_elevated_contam = (
                contam_fraction is None or contam_fraction >= self.subclade_contam_max_fraction
            )

            # Counting differentiating mutation categories
            damage_transitions: List[str] = []
            solid_variants: List[str] = []

            for mut in diff_muts:
                parsed = self._parse_single_mutation(mut, reference_sequence)
                if parsed and parsed.variant_type == "SNV":
                    pair = (parsed.ref_base, parsed.alt_base)
                    if pair in self.DAMAGE_TRANSITIONS and auth_status == "AUTHENTIC_ANCIENT_SUPPORTED":
                        damage_transitions.append(mut)
                    else:
                        solid_variants.append(mut)
                elif parsed and parsed.variant_type in {"INS", "DEL"}:
                    solid_variants.append(mut)

            # If contamination is elevated, reject shift outright
            if is_elevated_contam:
                flags.append(
                    AuditFlag(
                        code=AuditCode.SUSPECT_HETEROPLASMY_SUBCLADE_REJECTED,
                        severity=AuditSeverity.WARNING,
                        message=(
                            f"Subclade refinement to '{v_hg}' rejected due to elevated or unverified contamination "
                            f"(frac={contam_fraction}). Retaining primary consensus '{f_hg}'."
                        ),
                        context={
                            "primary_haplogroup": f_hg,
                            "rejected_subclade": v_hg,
                            "contamination_fraction": contam_fraction,
                        },
                    )
                )
                return True, f_hg, subclade_detected, False

            # Proportional Refinement: If supported by solid variants or verified heteroplasmy sites
            if solid_variants or (matching_sites and not damage_transitions):
                if damage_transitions:
                    flags.append(
                        AuditFlag(
                            code=AuditCode.DAMAGE_SUSPECT_PRIVATE_MUTATION,
                            severity=AuditSeverity.WARNING,
                            message=(
                                f"Subclade '{v_hg}' approved based on solid variants {solid_variants}, but "
                                f"carries suspect deamination transitions {damage_transitions}."
                            ),
                            context={"solid_variants": solid_variants, "suspect_damage": damage_transitions},
                        )
                    )
                flags.append(
                    AuditFlag(
                        code=AuditCode.HETEROPLASMY_SUBCLADE_SUPPORTED,
                        severity=AuditSeverity.INFO,
                        message=f"Subclade '{v_hg}' confirmed by differentiating variants.",
                        context={"positions": sorted(list(matching_sites)), "subclade": v_hg},
                    )
                )
                return True, v_hg, subclade_detected, True

            # If purely driven by deamination transitions without solid support, reject shift
            if damage_transitions and not solid_variants:
                flags.append(
                    AuditFlag(
                        code=AuditCode.SUSPECT_HETEROPLASMY_SUBCLADE_REJECTED,
                        severity=AuditSeverity.WARNING,
                        message=(
                            f"Subclade refinement to '{v_hg}' rejected: differentiating variants "
                            f"{damage_transitions} are solely deamination transitions on authentic ancient background."
                        ),
                        context={
                            "primary_haplogroup": f_hg,
                            "rejected_subclade": v_hg,
                            "damage_transitions": damage_transitions,
                        },
                    )
                )
                return True, f_hg, subclade_detected, False

            return True, f_hg, subclade_detected, False

        # Scenario B: FASTA is topologically deeper than VCF
        if PhylotreeGraph.is_descendant(candidate_child=f_hg, candidate_parent=v_hg):
            flags.append(
                AuditFlag(
                    code=AuditCode.SUBCLADE_FASTA_SPECIFIC,
                    severity=AuditSeverity.INFO,
                    message=f"Consensus FASTA resolved deeper node '{f_hg}' than VCF '{v_hg}'.",
                    context={"primary_fasta": f_hg, "vcf_subclade": v_hg},
                )
            )
            return True, f_hg, False, None

        # Scenario C: Incompatible divergent clades
        flags.append(
            AuditFlag(
                code=AuditCode.HAPLOGROUP_INPUT_DISCORDANCE,
                severity=AuditSeverity.WARNING,
                message=f"Topological divergence between inputs: FASTA='{f_hg}', VCF='{v_hg}'.",
                context={
                    "fasta_haplogroup": f_hg,
                    "vcf_haplogroup": v_hg,
                    "fasta_quality": fasta_result.quality_score,
                    "vcf_quality": vcf_result.quality_score,
                },
            )
        )
        return False, f_hg, False, None

    def _evaluate_contaminant_affinity(
        self,
        endogenous_hg: str,
        b5_data: Dict[str, Any],
    ) -> Tuple[ContaminantAffinity, List[str]]:
        raw_contam_hgs = b5_data.get("contaminant_haplogroups") or []
        if isinstance(raw_contam_hgs, str):
            raw_contam_hgs = [raw_contam_hgs]

        clean_contam_hgs: List[str] = [
            str(c).strip() for c in raw_contam_hgs
            if str(c).strip() and str(c).strip().upper() not in {"NONE", "UNKNOWN", "NA"}
        ]

        if not clean_contam_hgs or endogenous_hg.upper() in {"UNKNOWN", "N/A"}:
            return ContaminantAffinity.NOT_APPLICABLE, clean_contam_hgs

        endo_clean = endogenous_hg.strip().upper()
        affinities: List[ContaminantAffinity] = []

        for chg in clean_contam_hgs:
            target = chg.upper()
            if endo_clean == target:
                affinities.append(ContaminantAffinity.EXACT_MATCH)
            elif PhylotreeGraph.is_descendant(candidate_child=target, candidate_parent=endo_clean):
                affinities.append(ContaminantAffinity.SUBCLADE)
            elif PhylotreeGraph.is_descendant(candidate_child=endo_clean, candidate_parent=target):
                affinities.append(ContaminantAffinity.ANCESTRAL)
            else:
                affinities.append(ContaminantAffinity.DISTANT)

        for level in [
            ContaminantAffinity.EXACT_MATCH,
            ContaminantAffinity.SUBCLADE,
            ContaminantAffinity.ANCESTRAL,
        ]:
            if level in affinities:
                return level, clean_contam_hgs

        return ContaminantAffinity.DISTANT, clean_contam_hgs

    def _parse_single_mutation(
        self,
        mut_str: str,
        reference_sequence: Optional[str],
    ) -> Optional[ParsedMutation]:
        clean = mut_str.strip().upper()

        ins_match = self.INS_REGEX.match(clean)
        if ins_match:
            pos = int(ins_match.group(1))
            alt_b = ins_match.group(3).upper()
            return ParsedMutation(raw=clean, variant_type="INS", position=pos, ref_base=".", alt_base=alt_b)

        del_match = self.DEL_REGEX.match(clean)
        if del_match:
            pos = int(del_match.group(2))
            raw_ref = (del_match.group(1) or "").lstrip(">").upper()
            ref_b = raw_ref if raw_ref in {"A", "C", "G", "T"} else "N"
            if ref_b == "N" and reference_sequence is not None:
                if 1 <= pos <= len(reference_sequence):
                    ref_b = reference_sequence[pos - 1].upper()
            return ParsedMutation(raw=clean, variant_type="DEL", position=pos, ref_base=ref_b, alt_base="DEL")

        snv_match = self.SNV_REGEX.match(clean)
        if snv_match:
            pos = int(snv_match.group(2))
            raw_ref = (snv_match.group(1) or "").lstrip(">").upper()
            alt_b = snv_match.group(3).upper()
            ref_b = raw_ref if raw_ref in {"A", "C", "G", "T"} else "N"
            if ref_b == "N" and reference_sequence is not None:
                if 1 <= pos <= len(reference_sequence):
                    ref_b = reference_sequence[pos - 1].upper()
            return ParsedMutation(raw=clean, variant_type="SNV", position=pos, ref_base=ref_b, alt_base=alt_b)

        return ParsedMutation(raw=clean, variant_type="UNKNOWN", position=0, ref_base="N", alt_base="N")

    def _audit_private_mutations_chemistry(
        self,
        private_mutations: List[str],
        b4_data: Dict[str, Any],
        reference_sequence: Optional[str],
    ) -> List[PrivateMutationAssessment]:
        auth_status = str(b4_data.get("authentication_status") or "").upper()
        is_authentic_ancient = (auth_status == "AUTHENTIC_ANCIENT_SUPPORTED")
        assessments: List[PrivateMutationAssessment] = []

        for mut_str in private_mutations:
            clean_mut = mut_str.strip().upper()
            parsed = self._parse_single_mutation(clean_mut, reference_sequence)

            if not parsed or parsed.variant_type == "UNKNOWN":
                assessments.append(
                    PrivateMutationAssessment(
                        mutation=clean_mut,
                        position=0,
                        ref_base="N",
                        alt_base="N",
                        classification=PrivateMutationDamageClass.UNASSESSED,
                        evidence_source="BLOCK4_CHEMISTRY_SCREEN",
                        rationale="Unrecognized mutation notation.",
                    )
                )
                continue

            if parsed.variant_type == "INS":
                assessments.append(
                    PrivateMutationAssessment(
                        mutation=clean_mut,
                        position=parsed.position,
                        ref_base=parsed.ref_base,
                        alt_base=parsed.alt_base,
                        classification=PrivateMutationDamageClass.UNASSESSED,
                        evidence_source="BLOCK4_CHEMISTRY_SCREEN",
                        rationale="Insertion variant; indels bypass point substitution deamination triage.",
                    )
                )
                continue

            if parsed.variant_type == "DEL":
                # Check for suspect terminal deletion / soft-clip artifact in aDNA
                is_terminal_region = (parsed.position <= 500) or (parsed.position >= 16000)
                if is_authentic_ancient and is_terminal_region:
                    cls = PrivateMutationDamageClass.INDEL_SOFT_CLIP_SUSPECT
                    rationale = (
                        "Terminal deletion in authentic ancient sample; suspect mapping or soft-clipping artifact."
                    )
                else:
                    cls = PrivateMutationDamageClass.UNASSESSED
                    rationale = "Internal deletion variant; indels bypass substitution triage."

                assessments.append(
                    PrivateMutationAssessment(
                        mutation=clean_mut,
                        position=parsed.position,
                        ref_base=parsed.ref_base,
                        alt_base=parsed.alt_base,
                        classification=cls,
                        evidence_source="BLOCK4_CHEMISTRY_SCREEN",
                        rationale=rationale,
                    )
                )
                continue

            ref_b = parsed.ref_base
            pos = parsed.position
            alt_b = parsed.alt_base

            if ref_b not in {"A", "C", "G", "T"}:
                assessments.append(
                    PrivateMutationAssessment(
                        mutation=clean_mut,
                        position=pos,
                        ref_base=ref_b,
                        alt_base=alt_b,
                        classification=PrivateMutationDamageClass.UNASSESSED,
                        evidence_source="BLOCK4_CHEMISTRY_SCREEN",
                        rationale="Reference base unknown in notation and reference sequence unavailable.",
                    )
                )
                continue

            pair = (ref_b, alt_b)
            if pair not in self.TRANSITIONS:
                cls = PrivateMutationDamageClass.TRANSVERSION_NON_DAMAGE
                rationale = "Transversion substitution chemically incompatible with hydrolytic deamination."
            elif pair in self.DAMAGE_TRANSITIONS:
                if is_authentic_ancient:
                    cls = PrivateMutationDamageClass.TRANSITION_REQUIRES_SPATIAL_AUDIT
                    rationale = (
                        "Transition (C->T or G->A) on authentic ancient background. "
                        "Requires spatial terminal read audit in BAM to differentiate deamination from true variant."
                    )
                else:
                    cls = PrivateMutationDamageClass.TRANSITION_UNSUPPORTED_DAMAGE
                    rationale = "Transition observed on unauthenticated or modern reference background."
            else:
                cls = PrivateMutationDamageClass.TRANSITION_NON_DAMAGE
                rationale = "Transition (T->C or A->G) incompatible with typical hydrolytic deamination."

            assessments.append(
                PrivateMutationAssessment(
                    mutation=clean_mut,
                    position=pos,
                    ref_base=ref_b,
                    alt_base=alt_b,
                    classification=cls,
                    evidence_source="BLOCK4_CHEMISTRY_SCREEN",
                    rationale=rationale,
                )
            )

        return assessments

    def _extract_positions(self, mutations: Set[str]) -> Set[int]:
        positions: Set[int] = set()
        for m in mutations:
            parsed = self._parse_single_mutation(m, None)
            if parsed and parsed.position > 0:
                positions.add(parsed.position)
        return positions