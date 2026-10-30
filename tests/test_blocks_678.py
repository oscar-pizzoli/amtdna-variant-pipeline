"""
Pipeline Verification Script: Production-grade single sample executor.
Enforces cryptographic provenance, formal contract assembly (Blocks 4, 5, 6, 7),
deep invariant validation, and mutation cataloguing (Block 8).
"""

import csv
import hashlib
import inspect
import json
import logging
import os
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import pysam
import yaml

from src.modules.consensus.models import (
    ConsensusPolicy,
    LibraryTreatment,
    ReferenceInfo,
)
from src.modules.consensus.orchestrator import ConsensusOrchestrator
from src.modules.haplogrouping.orchestrator import HaplogroupOrchestrator
from src.modules.synthesis.coordinate_normalizer import MutationType
from src.modules.synthesis.models import MutationRole
from src.modules.synthesis.mutation_cataloguer import MutationCataloguer

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("PipelineVerification")


class VerificationError(Exception):
    """Raised when an integrity check, provenance assertion, or contract invariant fails."""
    pass


def calculate_sha256(filepath: Path, chunk_size: int = 65536) -> str:
    """Compute hex SHA-256 digest of a target file."""
    if not filepath.is_file():
        raise VerificationError(f"File non trovato per il calcolo SHA-256: {filepath}")
    hasher = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(chunk_size):
            hasher.update(chunk)
    return hasher.hexdigest()


def load_and_validate_reference(
    ref_fasta_path: Path, expected_length: Optional[int] = 16569
) -> Tuple[ReferenceInfo, str]:
    """
    Strict reference validator enforcing single-record integrity, contig dimension
    invariants, index (.fai) presence, and SHA-256 calculation (Block 1 parity).
    """
    if not ref_fasta_path.is_file():
        raise VerificationError(f"Reference FASTA non trovato: {ref_fasta_path}")

    fai_path = ref_fasta_path.with_name(f"{ref_fasta_path.name}.fai")
    if not fai_path.is_file():
        logger.info(f"Generazione indice FASTA mancante: {fai_path}")
        pysam.faidx(str(ref_fasta_path))
        if not fai_path.is_file():
            raise VerificationError(f"Impossibile creare l'indice FASTA: {fai_path}")

    headers: List[str] = []
    fragments: List[str] = []

    with open(ref_fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                headers.append(stripped[1:].split()[0])
                if len(headers) > 1:
                    raise VerificationError(
                        f"Rilevati record multipli nel riferimento ({headers}). "
                        "Il genoma mitocondriale deve essere rigorosamente a contig singolo."
                    )
            else:
                fragments.append(stripped.upper())

    if not headers or not fragments:
        raise VerificationError(f"File di riferimento vuoto o corrotto: {ref_fasta_path}")

    contig_name = headers[0]
    sequence = "".join(fragments)
    actual_length = len(sequence)

    if expected_length is not None and actual_length != expected_length:
        raise VerificationError(
            f"Dimensione riferimento errata per '{contig_name}': "
            f"attesi {expected_length} bp, riscontrati {actual_length} bp."
        )

    ref_sha256 = calculate_sha256(ref_fasta_path)
    ref_info = ReferenceInfo(
        contig=contig_name,
        length=actual_length,
        sequence_sha256=ref_sha256,
    )
    return ref_info, sequence


def resolve_upstream_contracts(
    results_root: Path, sample_id: str
) -> Tuple[Path, Path, Dict[str, Any], Dict[str, Any]]:
    """
    Deterministically locate and cross-validate Block 4 and Block 5 contracts.
    Enforces sample_id integrity and BAM existence with relative index (.bai).
    """
    b4_matches = sorted(results_root.rglob(f"*{sample_id}*block4_contract.json"))
    b5_matches = sorted(results_root.rglob(f"*{sample_id}*block5_contract.json"))

    if not b4_matches:
        raise VerificationError(f"Contratto Blocco 4 non trovato per '{sample_id}' sotto {results_root}")
    if not b5_matches:
        raise VerificationError(f"Contratto Blocco 5 non trovato per '{sample_id}' sotto {results_root}")

    b4_path = b4_matches[0]
    b5_path = b5_matches[0]

    with open(b4_path, "r", encoding="utf-8") as f:
        b4_data = json.load(f)
    with open(b5_path, "r", encoding="utf-8") as f:
        b5_data = json.load(f)

    b4_sample = b4_data.get("sample_id")
    b5_sample = b5_data.get("sample_id")

    if b4_sample != sample_id:
        raise VerificationError(f"Discrepanza sample_id in Blocco 4: atteso '{sample_id}', trovato '{b4_sample}'")
    if b5_sample != sample_id:
        raise VerificationError(f"Discrepanza sample_id in Blocco 5: atteso '{sample_id}', trovato '{b5_sample}'")

    bam_path_str = (
        b4_data.get("rescaled_bam")
        if b4_data.get("rescaling_performed")
        else (b4_data.get("filtered_dedup_bam") or b4_data.get("bam_path"))
    )
    if not bam_path_str or not Path(bam_path_str).is_file():
        raise VerificationError(f"BAM indicato nel Blocco 4 non trovato: {bam_path_str}")

    bam_path = Path(bam_path_str)
    bai_path = bam_path.with_suffix(".bam.bai")
    alt_bai = Path(f"{bam_path}.bai")
    if not bai_path.is_file() and not alt_bai.is_file():
        raise VerificationError(f"Indice BAM (.bai) mancante per {bam_path}")

    return b4_path, b5_path, b4_data, b5_data


def build_consensus_policy_from_config(cfg: Dict[str, Any]) -> ConsensusPolicy:
    """Build ConsensusPolicy by extracting thresholds from config.yaml."""
    c_cfg = cfg.get("consensus", {})
    s_cfg = cfg.get("spatial_damage", {})

    spatial_thresholds = {
        LibraryTreatment.NON_UDG: s_cfg.get("non_udg_threshold", 5),
        LibraryTreatment.UDG_HALF: s_cfg.get("half_udg_threshold", 3),
        LibraryTreatment.UDG_FULL: s_cfg.get("full_udg_threshold", 1),
    }

    return ConsensusPolicy(
        min_base_quality=int(c_cfg.get("min_base_quality", 20)),
        min_mapping_quality=int(c_cfg.get("min_mapping_quality", 30)),
        min_consensus_depth=int(c_cfg.get("min_consensus_depth", 5)),
        major_allele_fraction=float(c_cfg.get("major_allele_fraction", 0.80)),
        ambiguity_min_vaf=float(c_cfg.get("ambiguity_min_vaf", 0.20)),
        heteroplasmy_min_depth=int(c_cfg.get("heteroplasmy_min_depth", 10)),
        heteroplasmy_min_alt_reads=int(c_cfg.get("heteroplasmy_min_alt_reads", 3)),
        heteroplasmy_min_vaf=float(c_cfg.get("heteroplasmy_min_vaf", 0.10)),
        heteroplasmy_max_vaf=float(c_cfg.get("heteroplasmy_max_vaf", 0.90)),
        enable_indel_calling=bool(c_cfg.get("enable_indel_calling", True)),
        indel_min_depth=int(c_cfg.get("indel_min_depth", 10)),
        indel_min_alt_reads=int(c_cfg.get("indel_min_alt_reads", 3)),
        indel_min_vaf=float(c_cfg.get("indel_min_vaf", 0.10)),
        indel_max_vaf=float(c_cfg.get("indel_max_vaf", 0.90)),
        filter_homopolymer_slippage=bool(c_cfg.get("filter_homopolymer_slippage", True)),
        emit_iupac_fasta=bool(c_cfg.get("emit_iupac_fasta", True)),
        emit_subconsensus_vcf=bool(c_cfg.get("emit_subconsensus_vcf", True)),
        spatial_thresholds=spatial_thresholds,
    )


def ensure_and_verify_block6_contract(
    sample_id: str,
    consensus_dir: Path,
    qc_json_path: Path,
    fasta_path: Path,
    vcf_path: Path,
    indel_fasta_path: Path,
    het_tsv_path: Path,
) -> Path:
    """
    Assemble formal Block 6 contract after cryptographically verifying artifact existence,
    SHA-256 match, and strictly isolating positions that passed all quality filters.
    """
    for path, label in [
        (qc_json_path, "Consensus QC JSON"),
        (fasta_path, "Canonical Consensus FASTA"),
        (vcf_path, "Variants VCF"),
        (indel_fasta_path, "InDel-aware FASTA"),
        (het_tsv_path, "Heteroplasmy TSV"),
    ]:
        if not path.is_file():
            raise VerificationError(f"Artefatto Blocco 6 mancante ({label}): {path}")

    with open(qc_json_path, "r", encoding="utf-8") as f:
        qc_data = json.load(f)

    expected_sha = qc_data.get("consensus_fasta_sha256")
    if not expected_sha:
        raise VerificationError(f"Il file QC non contiene consensus_fasta_sha256: {qc_json_path}")

    actual_sha = calculate_sha256(fasta_path)
    if actual_sha.lower() != expected_sha.lower():
        raise VerificationError(
            f"Violazione integrità FASTA Blocco 6: atteso {expected_sha}, calcolato {actual_sha}"
        )

    pass_positions: Set[int] = set()
    with open(het_tsv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("status", "").upper() == "PASS":
                try:
                    pass_positions.add(int(row["position"]))
                except (ValueError, KeyError):
                    continue

    sorted_pass_positions = sorted(pass_positions)
    het_details = qc_data.get("heteroplasmy_details", {})

    contract_payload = {
        "sample_id": sample_id,
        "consensus_fasta_path": str(fasta_path.resolve()),
        "consensus_fasta_sha256": actual_sha,
        "variants_vcf_path": str(vcf_path.resolve()),
        "consensus_indel_fasta_path": str(indel_fasta_path.resolve()),
        "mean_depth": qc_data.get("mean_depth", 0.0),
        "covered_positions_fraction": qc_data.get("covered_positions_fraction", 0.0),
        "n_fraction_canonical": qc_data.get("n_fraction_canonical", 0.0),
        "heteroplasmy_sites_total": qc_data.get("heteroplasmy_sites_total", 0),
        "heteroplasmy_sites_pass": sorted_pass_positions,
        "gating_action": qc_data.get("gating_action", "STANDARD_CALLING"),
        "heteroplasmy_details": het_details,
    }

    contract_path = consensus_dir / f"{sample_id}.block6_contract.json"
    with open(contract_path, "w", encoding="utf-8") as f:
        json.dump(contract_payload, f, indent=2)

    logger.info(f"[{sample_id}] Contratto formale Blocco 6 certificato: {len(sorted_pass_positions)} siti PASS.")
    return contract_path


def validate_final_pipeline_state(
    sample_id: str,
    muts: List[Any],
    ref_seq: str,
    b7_contract: Any,
) -> None:
    """Deep scientific invariant checks on the synthesized catalog."""
    ref_len = len(ref_seq)
    role_counts = Counter(m.role.value for m in muts)
    type_counts = Counter(m.variant_type.value for m in muts)
    source_counts = Counter(m.stream_source for m in muts)

    if not muts:
        raise VerificationError(f"Catalogo mutazioni vuoto per '{sample_id}'.")

    if role_counts.get("DIAGNOSTIC", 0) == 0:
        raise VerificationError("Nessuna mutazione DIAGNOSTIC identificata: disallineamento filogenetico persistente.")

    if role_counts.get("PRIVATE", 0) == len(muts):
        raise VerificationError("Il 100% delle mutazioni è marcato PRIVATE: mancata ingestione dell'albero filogenetico.")

    if not any("PHYLOGENY" in m.stream_source for m in muts):
        raise VerificationError("Flusso PHYLOGENY assente dalle sorgenti mutazionali.")

    total_indels = type_counts.get("INS", 0) + type_counts.get("DEL", 0)
    if total_indels == 0:
        raise VerificationError("Zero InDel rilevati: la cattura strutturale mitocondriale non è attiva.")

    seen_keys: Set[Tuple[int, str, str]] = set()
    for m in muts:
        pos = m.position
        ref = m.ref
        alt = m.alt

        if not (1 <= pos <= ref_len):
            raise VerificationError(f"Coordinata non valida {pos} (fuori dall'intervallo [1, {ref_len}]).")

        if m.variant_type == MutationType.SNV and ref != "N":
            expected_ref_base = ref_seq[pos - 1]
            if ref != expected_ref_base:
                raise VerificationError(
                    f"Incoerenza reference al locus {pos}: variante dichiara {ref}, rCRS ha {expected_ref_base}."
                )

        if m.stream_source in {"PHYLOGENY", "HAPLOGREP"}:
            if m.vaf is not None or m.depth is not None:
                raise VerificationError(
                    f"Fabbricazione dati al locus {pos}: variante puramente filogenetica con DP/VAF definiti."
                )
        elif "VCF" in m.stream_source:
            if m.depth is None or m.depth < 0:
                raise VerificationError(f"Profondità (DP) non valida o mancante per variante da VCF al locus {pos}.")

        canon_key = (pos, ref, alt)
        if canon_key in seen_keys:
            raise VerificationError(f"Mutazione duplicata rilevata nel catalogo finale: {canon_key}")
        seen_keys.add(canon_key)


def main():
    project_root = Path(__file__).resolve().parent
    if project_root.name == "tests":
        project_root = project_root.parent

    results_root = project_root / "results"
    ref_fasta_path = project_root / "data" / "reference" / "rCRS.fasta"
    config_file = project_root / "config.yaml"

    if not config_file.is_file():
        config_file = project_root / "config" / "default_config.yaml"
    if not config_file.is_file():
        raise VerificationError(f"Configurazione non trovata: {config_file}")

    with open(config_file, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    # Lista dei campioni da verificare
    samples_to_verify = ["SRR30059011", "EGIZIANOANTICO1"]

    # 1. Certificazione Riferimento (Parità Blocco 1)
    ref_len_expected = cfg.get("reference", {}).get("expected_length", 16569)
    ref_info, ref_seq = load_and_validate_reference(ref_fasta_path, expected_length=ref_len_expected)
    logger.info(f"Riferimento validato: {ref_info.contig} ({ref_info.length} bp), SHA-256={ref_info.sequence_sha256[:12]}...")

    policy = build_consensus_policy_from_config(cfg)
    b7_orch = HaplogroupOrchestrator(cfg)

    cataloguer = MutationCataloguer(
        results_dir=results_root,
        reference_info=ref_info,
        reference_sequence=ref_seq,
    )

    for sample_id in samples_to_verify:
        logger.info(f"\n{'=' * 80}\n=== AVVIO VERIFICA RIGOROSA PIPELINE: {sample_id} ===\n{'=' * 80}")

        try:
            # 2. Risoluzione deterministica contratti a monte (Blocco 4 e 5)
            b4_path, b5_path, b4_data, b5_data = resolve_upstream_contracts(results_root, sample_id)
            logger.info(f"Contratti a monte validati per '{sample_id}' (Lineage certificata).")

            # Directory di produzione
            consensus_dir = results_root / "consensus" / sample_id
            haplo_dir = results_root / "haplogroup" / sample_id
            consensus_dir.mkdir(parents=True, exist_ok=True)
            haplo_dir.mkdir(parents=True, exist_ok=True)

            # 3. Esecuzione Blocco 6 con policy derivata dal config
            logger.info(f"[{sample_id}] Esecuzione Blocco 6 (Consensus & VCF)...")
            b6_orch = ConsensusOrchestrator(policy=policy)
            input_resolution = b6_orch.resolve_inputs(
                sample_id=sample_id,
                block4_json_path=str(b4_path),
                block5_json_path=str(b5_path),
                reference_fasta_path=str(ref_fasta_path),
                target_contig=ref_info.contig,
            )
            b6_orch.run(input_resolution, str(consensus_dir))

            # 4. Certificazione artefatti e assemblaggio formale Block 6 Contract
            consensus_fasta = consensus_dir / f"{sample_id}.consensus.fasta"
            variants_vcf = consensus_dir / f"{sample_id}.variants.vcf.gz"
            consensus_qc = consensus_dir / f"{sample_id}.consensus_qc.json"
            consensus_indel_fasta = consensus_dir / f"{sample_id}.consensus_indel.fasta"
            het_tsv = consensus_dir / f"{sample_id}.heteroplasmy.tsv"

            b6_contract_path = ensure_and_verify_block6_contract(
                sample_id=sample_id,
                consensus_dir=consensus_dir,
                qc_json_path=consensus_qc,
                fasta_path=consensus_fasta,
                vcf_path=variants_vcf,
                indel_fasta_path=consensus_indel_fasta,
                het_tsv_path=het_tsv,
            )

            # 5. Esecuzione Blocco 7 (HaploGrep3 Dual-Stream)
            logger.info(f"[{sample_id}] Esecuzione Blocco 7 (HaploGrep3 Dual-Stream)...")
            b7_contract = b7_orch.run_sample(
                sample_id=sample_id,
                consensus_fasta=consensus_fasta,
                sample_output_dir=haplo_dir,
                variants_vcf=variants_vcf,
                consensus_indel_fasta=consensus_indel_fasta,
                reference_fasta=ref_fasta_path,
                block4_contract_path=b4_path,
                block5_contract_path=b5_path,
                block6_contract_path=b6_contract_path,
            )

            if b7_contract.haplogroup_status.value == "REJECTED":
                raise VerificationError(f"Blocco 7 ha rigettato il campione: {b7_contract.flags}")
            logger.info(f"[{sample_id}] Aplogruppo: {b7_contract.primary_haplogroup} (Score: {b7_contract.quality_score:.4f})")

            # 6. Esecuzione Blocco 8 (MutationCataloguer)
            logger.info(f"[{sample_id}] Compilazione catalogo mutazioni Blocco 8...")
            muts = cataloguer.build_catalog(sample_id, require_vcf=True)

            # 7. Validazione profonda dello stato finale
            validate_final_pipeline_state(sample_id, muts, ref_seq, b7_contract)

            # 8. Report Finale
            role_counts = Counter(m.role.value for m in muts)
            type_counts = Counter(m.variant_type.value for m in muts)
            source_counts = Counter(m.stream_source for m in muts)

            print("\n" + "=" * 115)
            print(f"[*] REPORT FINALE CATALOGO MUTAZIONI (PROVENANCE CERTIFICATA): {sample_id}")
            print("=" * 115)
            print(f"Totale varianti registrate    : {len(muts)}")
            print(f"Ripartizione ruoli            : {dict(role_counts)}")
            print(f"Ripartizione tipi             : {dict(type_counts)}")
            print(f"Flussi di evidenza integrati  : {dict(source_counts)}")
            print("-" * 115)
            print(
                f"{'POS':<6} | {'REF>ALT':<10} | {'TYPE':<5} | {'NOTATION':<10} | "
                f"{'SOURCES':<26} | {'ROLE':<15} | {'VAF':<6} | {'DP':<6} | {'DMG?'}"
            )
            print("-" * 115)

            for m in muts:
                vaf_str = f"{m.vaf:.2f}" if m.vaf is not None else "NA"
                dp_str = str(m.depth) if m.depth is not None else "NA"
                dmg_str = str(m.spatial_damage_suspect) if m.spatial_damage_suspect is not None else "NA"
                print(
                    f"{m.position:<6} | {m.ref}>{m.alt:<8} | {m.variant_type.value:<5} | "
                    f"{m.mutation_notation:<10} | {m.stream_source:<26} | {m.role.value:<15} | "
                    f"{vaf_str:<6} | {dp_str:<6} | {dmg_str}"
                )
            print("=" * 115)
            logger.info(f"✅ [{sample_id}] TUTTI I CONTROLLI SUPERATI: integrità certificata.")

        except Exception as exc:
            logger.error(f"❌ Fallimento durante l'elaborazione di {sample_id}: {exc}")
            raise
    logger.info("✅ TUTTI I CONTROLLI SUPERATI: integrità biologica e contratti formali certificati.")


if __name__ == "__main__":
    main()