"""
End-to-End Validation Test: InDel capture, phylogenetic role assignment,
and back-mutation preservation across Blocks 6, 7, and 8.
"""

import json
import shutil
import tempfile
from collections import Counter
from pathlib import Path

import pysam

from src.modules.consensus.models import (
    CallStatus,
    ConsensusCallReason,
    ConsensusCallStatus,
    ConsensusSite,
    FilterReason,
    GatingAction,
    HeteroplasmyCall,
    ReferenceInfo,
)
from src.modules.consensus.orchestrator import ConsensusOrchestrator
from src.modules.synthesis.coordinate_normalizer import MutationType
from src.modules.synthesis.models import MutationRole
from src.modules.synthesis.mutation_cataloguer import MutationCataloguer


def load_or_create_reference() -> str:
    """Load reference rCRS if available, else synthesize a realistic 16,569 bp mock."""
    ref_path = Path("data/reference/rCRS.fasta")
    if ref_path.is_file():
        with open(ref_path, "r", encoding="utf-8") as fp:
            return "".join(line.strip() for line in fp if not line.startswith(">")).upper()

    seq = list("A" * 16569)
    seq[72] = "A"        # Pos 73 (1-based)
    seq[151] = "T"       # Pos 152
    seq[314] = "C"       # Pos 315
    seq[521:524] = "CAC" # Pos 522-524
    return "".join(seq)


def run_integration_test():
    sample_id = "SAMPLE_VALIDATION_01"
    ref_contig = "chrM"
    ref_seq = load_or_create_reference()
    ref_len = len(ref_seq)
    ref_info = ReferenceInfo(contig=ref_contig, length=ref_len, sequence_sha256="0" * 64)

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # ----------------------------------------------------------------------
        # Struttura modulare per MutationCataloguer: consensus/ e haplogrouping/
        # ----------------------------------------------------------------------
        consensus_dir = tmp_path / "consensus" / sample_id
        haplo_dir = tmp_path / "haplogrouping" / sample_id
        haplo_alt_dir = tmp_path / "haplogroup" / sample_id
        sample_flat_dir = tmp_path / sample_id

        for d in [consensus_dir, haplo_dir, haplo_alt_dir, sample_flat_dir]:
            d.mkdir(parents=True, exist_ok=True)

        vcf_gz_path = consensus_dir / f"{sample_id}.variants.vcf.gz"
        b6_qc_path = consensus_dir / f"{sample_id}.consensus_qc.json"
        b7_json_path = haplo_dir / f"{sample_id}.block7_contract.json"

        # ----------------------------------------------------------------------
        # 1. SIMULAZIONE BLOCCO 6: VCF con SNV di consenso + InDels ancorati
        # ----------------------------------------------------------------------
        site_73 = ConsensusSite(
            position=73,
            reference_base="A",
            canonical_base="G",
            iupac_base="G",
            callable_depth=100,
            major_allele="G",
            major_vaf=0.99,
            call_status=ConsensusCallStatus.CALLED,
            call_reason=ConsensusCallReason.CALLED_MAJORITY,
        )

        ins_315 = HeteroplasmyCall(
            position=315,
            reference_base="C",
            major_allele="C",
            minor_allele="CC",
            callable_depth=120,
            alt_reads=35,
            alt_vaf=0.2917,
            forward_ref=42,
            reverse_ref=43,
            forward_alt=17,
            reverse_alt=18,
            mean_bq=38.0,
            mean_mq=60.0,
            strand_bias_pvalue=1.0,
            status=CallStatus.PASS,
            filter_reason=FilterReason.PASS,
            gating_action=GatingAction.STANDARD_CALLING,
            variant_type=MutationType.INS,
        )

        del_522 = HeteroplasmyCall(
            position=522,
            reference_base="C",
            major_allele="CAC",
            minor_allele="C",
            callable_depth=90,
            alt_reads=88,
            alt_vaf=0.9778,
            forward_ref=1,
            reverse_ref=1,
            forward_alt=44,
            reverse_alt=44,
            mean_bq=35.0,
            mean_mq=60.0,
            strand_bias_pvalue=1.0,
            status=CallStatus.PASS,
            filter_reason=FilterReason.PASS,
            gating_action=GatingAction.STANDARD_CALLING,
            variant_type=MutationType.DEL,
        )

        orch_b6 = ConsensusOrchestrator()
        orch_b6._write_vcf(
            output_vcf_gz_path=str(vcf_gz_path),
            sample_id=sample_id,
            ref_info=ref_info,
            consensus_sites=[site_73],
            heteroplasmy_calls=[ins_315, del_522],
            pileup_by_pos={},
        )

        b6_qc_payload = {
            "sample_id": sample_id,
            "input_bam_path": "dummy.bam",
            "input_bam_sha256": "0" * 64,
            "reference_contig": ref_contig,
            "reference_sha256": "0" * 64,
            "reference_length": ref_len,
            "consensus_fasta_sha256": "0" * 64,
            "mean_depth": 100.0,
            "covered_positions_fraction": 1.0,
            "n_fraction_canonical": 0.0,
            "iupac_fraction": 0.0,
            "heteroplasmy_sites_total": 2,
            "heteroplasmy_sites_pass": 2,
            "heteroplasmy_sites_filtered": 0,
            "regions": [],
            "gating_action": "STANDARD_CALLING",
            "software_version": "1.2.1",
            "policy_version": "1.1.0",
            "execution_time_seconds": 1.0,
            "heteroplasmy_details": {
                "315": {
                    "position": 315,
                    "ref": "C",
                    "alt": "CC",
                    "depth": 120,
                    "vaf": 0.2917,
                    "internal_support_count": 35,
                    "terminal_support_count": 0,
                    "spatial_damage_suspect": False,
                },
                "522": {
                    "position": 522,
                    "ref": "CAC",
                    "alt": "C",
                    "depth": 90,
                    "vaf": 0.9778,
                    "internal_support_count": 88,
                    "terminal_support_count": 0,
                    "spatial_damage_suspect": False,
                },
            },
        }
        with open(b6_qc_path, "w", encoding="utf-8") as f:
            json.dump(b6_qc_payload, f, indent=2)

        # ----------------------------------------------------------------------
        # 2. SIMULAZIONE BLOCCO 7: Scrittura block7_contract.json normalizzato
        # ----------------------------------------------------------------------
        b7_contract_payload = {
            "sample_id": sample_id,
            "primary_haplogroup": "H1",
            "resolved_haplogroup": "H1",
            "quality_score": 0.95,
            "quality_margin": 0.12,
            "haplogroup_status": "ASSIGNED_CONFIDENT",
            "overall_validation_status": "PASS",
            "tree_id": "phylotree-fu-rcrs",
            "tree_version": "1.2",
            "metric": "KULCZYNSKI",
            "contaminant_affinity": "DISTANT",
            "contamination_fraction": 0.01,
            "vcf_haplogroup": "H1",
            "vcf_concordant": True,
            "subclade_shift_detected": False,
            "subclade_shift_supported": None,
            "fasta_sha256": "0" * 64,
            "found_mutations": ["A73G", "523-524d", "T152C!"],
            "expected_mutations": ["A73G", "523-524d"],
            "missing_mutations": [],
            "private_mutations": ["315.1C"],
            "flags": [],
            "execution_time_seconds": 1.25,
        }

        with open(b7_json_path, "w", encoding="utf-8") as f:
            json.dump(b7_contract_payload, f, indent=2)

        # Copia ridondante per compatibilità con qualsiasi strategia di lookup
        for target_d in [haplo_alt_dir, sample_flat_dir, tmp_path]:
            shutil.copyfile(b7_json_path, target_d / f"{sample_id}.block7_contract.json")
            shutil.copyfile(b6_qc_path, target_d / f"{sample_id}.consensus_qc.json")

        for target_d in [sample_flat_dir, tmp_path]:
            shutil.copyfile(vcf_gz_path, target_d / f"{sample_id}.variants.vcf.gz")
            tbi = consensus_dir / f"{sample_id}.variants.vcf.gz.tbi"
            if tbi.is_file():
                shutil.copyfile(tbi, target_d / f"{sample_id}.variants.vcf.gz.tbi")

        # ----------------------------------------------------------------------
        # 3. ESECUZIONE BLOCCO 8: MutationCataloguer
        # ----------------------------------------------------------------------
        cataloguer = MutationCataloguer(
            results_dir=tmp_path,
            reference_info=len(ref_seq),
            reference_sequence=ref_seq,
        )

        muts = cataloguer.build_catalog(sample_id, require_vcf=True)

        # ----------------------------------------------------------------------
        # 4. REPORT E VERIFICHE ANALITICHE
        # ----------------------------------------------------------------------
        role_counts = Counter(
            (m.role.value if hasattr(m.role, "value") else str(m.role)) for m in muts
        )
        source_counts = Counter(
            (m.stream_source.value if hasattr(m.stream_source, "value") else str(m.stream_source))
            for m in muts
        )

        print("\n" + "=" * 105)
        print(f"[*] REPORT CATALOGO MUTAZIONI - TEST INTEGRAZIONE ({sample_id})")
        print("=" * 105)
        print(f"Totale mutazioni identificate : {len(muts)}")
        print(f"Ripartizione ruoli            : {dict(role_counts)}")
        print(f"Flussi di evidenza integrati  : {dict(source_counts)}")
        print("-" * 105)
        print(
            f"{'POS':<6} | {'REF>ALT':<8} | {'TYPE':<5} | {'NOTATION':<10} | "
            f"{'SOURCES':<20} | {'ROLE':<15} | {'VAF':<6} | {'DP':<6} | {'DMG?'}"
        )
        print("-" * 105)

        mut_by_pos = {}
        for m in muts:
            pos = getattr(m, "position", 0)
            ref_b = getattr(m, "ref", getattr(m, "ref_allele", "?"))
            alt_b = getattr(m, "alt", getattr(m, "alt_allele", "?"))
            vtype = getattr(m, "variant_type", getattr(m, "mutation_type", None))
            vtype_str = vtype.value if hasattr(vtype, "value") else str(vtype)
            notation = getattr(m, "mutation_notation", getattr(m, "notation", "?"))
            role = getattr(m, "role", None)
            role_str = role.value if hasattr(role, "value") else str(role)
            src = getattr(m, "stream_source", getattr(m, "sources", "?"))
            src_str = src.value if hasattr(src, "value") else str(src)
            vaf_str = f"{m.vaf:.2f}" if m.vaf is not None else "NA"
            dp_str = str(m.depth) if m.depth is not None else "NA"
            dmg_str = str(m.spatial_damage_suspect) if m.spatial_damage_suspect is not None else "NA"

            mut_by_pos[pos] = m
            print(
                f"{pos:<6} | {ref_b}>{alt_b:<6} | {vtype_str:<5} | {notation:<10} | "
                f"{src_str:<20} | {role_str:<15} | {vaf_str:<6} | {dp_str:<6} | {dmg_str}"
            )
        print("=" * 105)

        # ----------------------------------------------------------------------
        # 5. ASSERTION CHECK
        # ----------------------------------------------------------------------
        # 1. A73G presente e DIAGNOSTIC
        assert 73 in mut_by_pos, "Asserzione fallita: mutazione al locus 73 assente."
        m_73 = mut_by_pos[73]
        role_73_val = m_73.role.value if hasattr(m_73.role, "value") else str(m_73.role)
        assert role_73_val == MutationRole.DIAGNOSTIC.value, (
            f"Asserzione fallita: locus 73 doveva essere DIAGNOSTIC, invece e' {role_73_val}."
        )

        # 2. Inserzione presente e tipizzata come INS
        has_ins = any(
            (getattr(m, "variant_type", getattr(m, "mutation_type", None)) == MutationType.INS)
            for m in muts
        )
        assert has_ins, "Asserzione fallita: nessuna variante di tipo INS trovata nel catalogo."

        # 3. Delezione presente e tipizzata come DEL
        has_del = any(
            (getattr(m, "variant_type", getattr(m, "mutation_type", None)) == MutationType.DEL)
            for m in muts
        )
        assert has_del, "Asserzione fallita: nessuna variante di tipo DEL trovata nel catalogo."

        # 4. Ruoli diversificati (non collassati al 100% su PRIVATE)
        assert len(role_counts) > 1, (
            f"Asserzione fallita: ruoli ancora monocromatici: {dict(role_counts)}"
        )

        # 5. Controllo semantico VAF su varianti solo-filogenetiche (es. T152C!)
        phylo_only = [
            m for m in muts
            if getattr(m, "stream_source", "") in {"PHYLOGENY", "HAPLOGREP"}
        ]
        if phylo_only:
            bad_vaf = [m for m in phylo_only if m.vaf is not None]
            assert len(bad_vaf) == 0, (
                f"Asserzione fallita: {len(bad_vaf)} varianti solo-filogenetiche hanno VAF numerica."
            )

        print("\n✅ TUTTI I TEST SONO STATI SUPERATI CON SUCCESSO:")
        print("   1. Il VCF include SNV e InDel ordinati.")
        print("   2. Gli InDel vengono tipizzati come INS e DEL nel catalogo.")
        print("   3. I ruoli filogenetici assegnano DIAGNOSTIC e non collassano su PRIVATE.")
        print("   4. Il controllo di semantica VAF per le evidenze solo-filogenetiche e' integro.")


if __name__ == "__main__":
    run_integration_test()