"""
Block 6 Observational Engine: Gapless 1..N coordinate streaming, single-molecule
mate overlap deduplication, physical fragment boundary awareness, and structural InDel capture.
"""

from collections import defaultdict
from typing import Any, Dict, Iterator, List, Optional, Tuple

import pysam

from src.modules.consensus.models import (
    AlleleCounts,
    ConsensusPolicy,
    IndelObservation,
    PileupFilterMetrics,
    SitePileup,
)
from src.modules.synthesis.coordinate_normalizer import MutationType


class PileupEngine:
    """
    Streams gapless 1-based genomic coordinates across a target contig.
    Resolves overlapping mate pairs into single-molecule observations, distinguishes
    physical fragment boundaries on paired vs merged reads, tracks granular filter
    exclusions, captures anchored insertions/deletions with circular wrapping, and
    computes per-allele distance vectors.
    """

    def __init__(
        self,
        bam_path: str,
        reference_fasta_path: str,
        policy: ConsensusPolicy,
        target_contig: str,
        expected_reference_length: int,
    ):
        self.bam_path = bam_path
        self.reference_fasta_path = reference_fasta_path
        self.policy = policy
        self.target_contig = target_contig
        self.expected_reference_length = expected_reference_length

    def _fetch_reference_base(self, fasta: pysam.FastaFile, position: int) -> str:
        """Fetch 1-based single reference base."""
        base = fasta.fetch(self.target_contig, position - 1, position)
        return base.upper() if base else "N"

    def _fetch_reference_span(self, fasta: pysam.FastaFile, start_pos: int, length: int) -> str:
        """
        Fetch 1-based multi-base reference sequence span supporting circular wraparound
        across the origin (16569 -> 1).
        """
        if length <= 0:
            return ""

        end_pos = start_pos + length - 1
        if end_pos <= self.expected_reference_length:
            seq = fasta.fetch(self.target_contig, start_pos - 1, end_pos)
            return seq.upper() if seq else "N" * length
        else:
            # Wraparound event across the origin boundary
            part1_len = self.expected_reference_length - start_pos + 1
            part2_len = length - part1_len

            seq1 = fasta.fetch(self.target_contig, start_pos - 1, self.expected_reference_length) or ("N" * part1_len)
            seq2 = fasta.fetch(self.target_contig, 0, min(part2_len, self.expected_reference_length)) or ("N" * part2_len)
            return (seq1 + seq2).upper()

    @staticmethod
    def _compute_fragment_end_distance(
        aln: pysam.AlignedSegment,
        qpos: int,
        q_start: int,
        q_end: int,
    ) -> int:
        """
        Compute true physical molecule end distance (in cycles/bases).
        For single-end/merged: min distance to both ends of the read (insert boundaries).
        For paired-end: distance to the external physical break point of the insert.
        """
        if not aln.is_paired:
            d5 = qpos - q_start
            d3 = q_end - qpos - 1
            return max(0, min(d5, d3))

        tlen = abs(aln.template_length)
        d_ext = qpos - q_start

        if tlen > 0 and tlen >= (q_end - q_start):
            dist_to_opposite = tlen - (qpos - q_start)
            return max(0, min(d_ext, dist_to_opposite))

        return max(0, d_ext)

    @staticmethod
    def _resolve_molecule_strand(aln1: pysam.AlignedSegment, aln2: Optional[pysam.AlignedSegment] = None) -> bool:
        """
        Resolve the physical fragment strand under the CANONICAL_FRAGMENT_STRAND policy.
        For paired reads, strand orientation is defined by Read 1 (5' anchor of original molecule).
        """
        if aln2 is None or not aln1.is_paired:
            return aln1.is_reverse

        if aln1.is_read1:
            return aln1.is_reverse
        return aln2.is_reverse

    def _resolve_indel_overlap(
        self,
        mate_info: List[Tuple[pysam.AlignedSegment, int, str, int, int, int]],
        fasta: pysam.FastaFile,
        pos: int,
    ) -> Optional[Tuple[MutationType, str, int, bool, int, int, int]]:
        """
        Deterministically resolve InDel observations across overlapping mate pairs.
        Returns: (event_type, sequence, length, is_reverse, bq, mq, dist_from_end) or None.
        """
        indel_mates = [m for m in mate_info if m[5] != 0]
        if not indel_mates:
            return None

        molecule_is_rev = self._resolve_molecule_strand(mate_info[0][0], mate_info[1][0] if len(mate_info) > 1 else None)
        fragment_dist = min(m[4] for m in mate_info)

        # Case 1: Concordant InDel calls across mates
        unique_indels = {m[5] for m in indel_mates}
        if len(unique_indels) == 1:
            chosen_indel = indel_mates[0][5]
            best_m = max(indel_mates, key=lambda x: (x[3], x[0].mapping_quality))
            aln_obj, qp, bq, mq = best_m[0], best_m[1], best_m[3], best_m[0].mapping_quality

            if chosen_indel > 0:
                ins_seq = aln_obj.query_sequence[qp + 1 : qp + 1 + chosen_indel].upper()
                if ins_seq:
                    return (MutationType.INS, ins_seq, chosen_indel, molecule_is_rev, bq, mq, fragment_dist)
            else:
                del_len = abs(chosen_indel)
                del_seq = self._fetch_reference_span(fasta, pos + 1, del_len)
                if del_seq:
                    return (MutationType.DEL, del_seq, del_len, molecule_is_rev, bq, mq, fragment_dist)

        # Case 2: One mate calls an InDel, other calls reference
        if len(indel_mates) == 1 and len(mate_info) > 1:
            m = indel_mates[0]
            if m[3] >= 20 and m[0].mapping_quality >= 20:
                aln_obj, qp, chosen_indel, bq, mq = m[0], m[1], m[5], m[3], m[0].mapping_quality
                if chosen_indel > 0:
                    ins_seq = aln_obj.query_sequence[qp + 1 : qp + 1 + chosen_indel].upper()
                    if ins_seq:
                        return (MutationType.INS, ins_seq, chosen_indel, molecule_is_rev, bq, mq, fragment_dist)
                else:
                    del_len = abs(chosen_indel)
                    del_seq = self._fetch_reference_span(fasta, pos + 1, del_len)
                    if del_seq:
                        return (MutationType.DEL, del_seq, del_len, molecule_is_rev, bq, mq, fragment_dist)
            return None

        # Case 3: Discordant InDels between mates (e.g. +C vs +CC)
        mate_info_sorted = sorted(indel_mates, key=lambda x: (x[3], x[0].mapping_quality), reverse=True)
        best_m = mate_info_sorted[0]
        second_m = mate_info_sorted[1]

        if (best_m[3] - second_m[3]) > 5:
            chosen_indel = best_m[5]
            aln_obj, qp, bq, mq = best_m[0], best_m[1], best_m[3], best_m[0].mapping_quality
            if chosen_indel > 0:
                ins_seq = aln_obj.query_sequence[qp + 1 : qp + 1 + chosen_indel].upper()
                if ins_seq:
                    return (MutationType.INS, ins_seq, chosen_indel, molecule_is_rev, bq, mq, fragment_dist)
            else:
                del_len = abs(chosen_indel)
                del_seq = self._fetch_reference_span(fasta, pos + 1, del_len)
                if del_seq:
                    return (MutationType.DEL, del_seq, del_len, molecule_is_rev, bq, mq, fragment_dist)

        return None

    def _create_empty_site(self, position: int, ref_base: str) -> SitePileup:
        """Construct a zero-depth observational SitePileup for uncovered loci."""
        return SitePileup(
            position=position,
            reference_base=ref_base,
            depth_pileup_entries=0,
            depth_filtered=0,
            depth_callable=0,
            indel_anchor_depth=0,
            counts=AlleleCounts(),
            forward_counts=AlleleCounts(),
            reverse_counts=AlleleCounts(),
            mean_bq=0.0,
            mean_mq=0.0,
            filter_metrics=PileupFilterMetrics(),
            read_end_distances={},
            insertions={},
            deletions={},
        )

    def iter_sites(self) -> Iterator[SitePileup]:
        """
        Stream gapless observational SitePileup objects from 1 to expected_reference_length.
        Applies single-molecule collapsing and quality filtering with zero silent drops.
        """
        with pysam.FastaFile(self.reference_fasta_path) as fasta, pysam.AlignmentFile(
            self.bam_path, "rb"
        ) as bam:
            pileup_iter = bam.pileup(
                contig=self.target_contig,
                start=0,
                stop=self.expected_reference_length,
                truncate=True,
                stepper="nofilter",
                min_base_quality=0,
                max_depth=200000,
            )

            current_pos = 1
            for col in pileup_iter:
                col_pos = col.reference_pos + 1

                while current_pos < col_pos and current_pos <= self.expected_reference_length:
                    ref_b = self._fetch_reference_base(fasta, current_pos)
                    yield self._create_empty_site(current_pos, ref_b)
                    current_pos += 1

                if current_pos > self.expected_reference_length:
                    break

                ref_b = self._fetch_reference_base(fasta, col_pos)
                site_pileup = self._process_column(col, ref_b, fasta)
                yield site_pileup
                current_pos += 1

            while current_pos <= self.expected_reference_length:
                ref_b = self._fetch_reference_base(fasta, current_pos)
                yield self._create_empty_site(current_pos, ref_b)
                current_pos += 1

    def _process_column(
        self,
        col: pysam.PileupColumn,
        ref_base: str,
        fasta: pysam.FastaFile,
    ) -> SitePileup:
        """
        Execute molecule-aware filtering on a single pileup column:
        1. Parse base calls and CIGAR InDels.
        2. Group paired reads by template (query_name) to collapse overlaps.
        3. Determine true physical molecule termini.
        4. Apply MQ and BQ filters, computing independent SNV and InDel anchor depths.
        """
        pos = col.reference_pos + 1
        total_entries = len(col.pileups)

        unmapped = 0
        qc_fail = 0
        duplicates = 0
        secondary = 0
        supplementary = 0
        structural_del_spans = 0
        refskips = 0
        no_query_position = 0
        collapsed_mates_count = 0
        low_mq = 0
        low_bq = 0
        noncanonical_bases = 0
        raw_indel_obs_count = 0

        # Pass 1: Structural defect filters & soft-clip boundary rejection
        # Entry: (alignment, qpos, q_start, q_end, is_del, indel_val)
        template_reads: Dict[str, List[Tuple[pysam.AlignedSegment, Optional[int], int, int, bool, int]]] = defaultdict(list)

        for read in col.pileups:
            aln: pysam.AlignedSegment = read.alignment

            if aln.is_unmapped:
                unmapped += 1
                continue
            if aln.is_qcfail:
                qc_fail += 1
                continue
            if aln.is_duplicate:
                duplicates += 1
                continue
            if aln.is_secondary:
                secondary += 1
                continue
            if aln.is_supplementary:
                supplementary += 1
                continue
            if read.is_refskip:
                refskips += 1
                continue

            # Molecule covers a deletion spanning this coordinate
            if read.is_del:
                structural_del_spans += 1
                template_reads[aln.query_name].append((aln, None, 0, 0, True, 0))
                continue

            if read.query_position is None:
                no_query_position += 1
                continue

            qpos = read.query_position
            q_start = aln.query_alignment_start if aln.query_alignment_start is not None else 0
            q_end = (
                aln.query_alignment_end
                if aln.query_alignment_end is not None
                else (len(aln.query_sequence) if aln.query_sequence else 0)
            )

            if qpos < q_start or qpos >= q_end:
                no_query_position += 1
                continue

            indel_val = getattr(read, "indel", 0)
            if indel_val != 0:
                raw_indel_obs_count += 1

            template_reads[aln.query_name].append((aln, qpos, q_start, q_end, False, indel_val))

        # Pass 2: Single-molecule deduplication & Physical Fragment Boundary Calculation
        # SNV Tuple: (aln, base, is_reverse, bq, mq, dist_from_end)
        resolved_snv_obs: List[Tuple[pysam.AlignedSegment, str, bool, int, int, int]] = []
        # InDel Tuple: (event_type, sequence, length, is_reverse, bq, mq, dist_from_end)
        resolved_indel_obs: List[Tuple[MutationType, str, int, bool, int, int, int]] = []

        # Count total callable molecules at this anchor position (for InDel VAF denominator)
        indel_anchor_depth = 0

        for qname, obs_list in template_reads.items():
            if len(obs_list) == 1:
                aln, qpos, q_start, q_end, is_deletion_span, indel_val = obs_list[0]

                if is_deletion_span or qpos is None:
                    continue

                base = aln.query_sequence[qpos].upper() if aln.query_sequence else "N"
                bq = aln.query_qualities[qpos] if aln.query_qualities is not None else 0
                mq = aln.mapping_quality

                dist_from_end = self._compute_fragment_end_distance(aln, qpos, q_start, q_end)
                molecule_is_rev = aln.is_reverse

                resolved_snv_obs.append((aln, base, molecule_is_rev, bq, mq, dist_from_end))

                # Track qualifying molecule for InDel anchor denominator
                if mq >= self.policy.min_mapping_quality and bq >= self.policy.min_base_quality:
                    indel_anchor_depth += 1

                # Capture single-read InDel
                if indel_val > 0:
                    ins_seq = aln.query_sequence[qpos + 1 : qpos + 1 + indel_val].upper()
                    if ins_seq:
                        resolved_indel_obs.append((MutationType.INS, ins_seq, indel_val, molecule_is_rev, bq, mq, dist_from_end))
                elif indel_val < 0:
                    del_len = abs(indel_val)
                    del_seq = self._fetch_reference_span(fasta, pos + 1, del_len)
                    if del_seq:
                        resolved_indel_obs.append((MutationType.DEL, del_seq, del_len, molecule_is_rev, bq, mq, dist_from_end))

            else:
                # OVERLAPPING MATE PAIR: Represent the same physical molecule
                collapsed_mates_count += len(obs_list) - 1

                usable_mates = [m for m in obs_list if not m[4] and m[1] is not None]
                if not usable_mates:
                    continue

                mate_info = []
                for aln, qpos, q_start, q_end, _, indel_val in usable_mates:
                    b = aln.query_sequence[qpos].upper() if aln.query_sequence else "N"
                    q = aln.query_qualities[qpos] if aln.query_qualities is not None else 0
                    dist = self._compute_fragment_end_distance(aln, qpos, q_start, q_end)
                    mate_info.append((aln, qpos, b, q, dist, indel_val))

                fragment_dist = min(m[4] for m in mate_info)
                molecule_is_rev = self._resolve_molecule_strand(mate_info[0][0], mate_info[1][0] if len(mate_info) > 1 else None)

                # Resolve SNV observation preserving selected_bq
                unique_bases = {m[2] for m in mate_info}
                if len(unique_bases) == 1:
                    best_mate = max(mate_info, key=lambda m: m[3])
                    selected_bq = best_mate[3]
                    selected_mq = best_mate[0].mapping_quality
                    resolved_snv_obs.append(
                        (best_mate[0], best_mate[2], molecule_is_rev, selected_bq, selected_mq, fragment_dist)
                    )
                else:
                    mate_info.sort(key=lambda m: m[3], reverse=True)
                    best_mate = mate_info[0]
                    second_mate = mate_info[1]

                    selected_bq = best_mate[3]
                    selected_mq = best_mate[0].mapping_quality
                    quality_delta = selected_bq - second_mate[3]

                    if quality_delta > 0:
                        resolved_snv_obs.append(
                            (best_mate[0], best_mate[2], molecule_is_rev, selected_bq, selected_mq, fragment_dist)
                        )
                    else:
                        resolved_snv_obs.append(
                            (best_mate[0], "N", molecule_is_rev, 0, selected_mq, fragment_dist)
                        )

                # Evaluate molecule for InDel anchor denominator
                highest_bq = max(m[3] for m in mate_info)
                highest_mq = max(m[0].mapping_quality for m in mate_info)
                if highest_mq >= self.policy.min_mapping_quality and highest_bq >= self.policy.min_base_quality:
                    indel_anchor_depth += 1

                # Resolve InDel observation across overlapping mates
                resolved_indel = self._resolve_indel_overlap(mate_info, fasta, pos)
                if resolved_indel is not None:
                    resolved_indel_obs.append(resolved_indel)

        # Pass 3: Evaluate MQ, BQ thresholds and compile SNV AlleleCounts
        counts_map = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "other": 0}
        fwd_counts_map = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "other": 0}
        rev_counts_map = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "other": 0}
        read_end_distances: Dict[str, List[int]] = defaultdict(list)

        total_bq = 0
        total_mq = 0
        depth_filtered = 0
        callable_depth = 0

        for aln, base, is_rev, bq, mq, dist in resolved_snv_obs:
            if mq < self.policy.min_mapping_quality:
                low_mq += 1
                continue
            if bq < self.policy.min_base_quality:
                low_bq += 1
                continue

            depth_filtered += 1
            total_bq += bq
            total_mq += mq

            if base in {"A", "C", "G", "T"}:
                callable_depth += 1
                counts_map[base] += 1
                read_end_distances[base].append(dist)
                if is_rev:
                    rev_counts_map[base] += 1
                else:
                    fwd_counts_map[base] += 1
            else:
                noncanonical_bases += 1
                if base == "N":
                    counts_map["N"] += 1
                else:
                    counts_map["other"] += 1

        counts = AlleleCounts(**counts_map)
        fwd_counts = AlleleCounts(**fwd_counts_map)
        rev_counts = AlleleCounts(**rev_counts_map)

        mean_bq = (total_bq / depth_filtered) if depth_filtered > 0 else 0.0
        mean_mq = (total_mq / depth_filtered) if depth_filtered > 0 else 0.0

        # Pass 4: Compile InDel Observations (Evidence collection with spatial distance tracking)
        insertions_map: Dict[str, IndelObservation] = {}
        deletions_map: Dict[str, IndelObservation] = {}

        ins_tracker: Dict[str, Dict[str, Any]] = defaultdict(lambda: {"fwd": 0, "rev": 0, "len": 0, "distances": []})
        del_tracker: Dict[str, Dict[str, Any]] = defaultdict(lambda: {"fwd": 0, "rev": 0, "len": 0, "distances": []})

        for ev_type, seq, ev_len, is_rev, bq, mq, dist in resolved_indel_obs:
            if mq < self.policy.min_mapping_quality or bq < self.policy.min_base_quality:
                continue

            if ev_type == MutationType.INS:
                strand_key = "rev" if is_rev else "fwd"
                ins_tracker[seq][strand_key] += 1
                ins_tracker[seq]["len"] = ev_len
                ins_tracker[seq]["distances"].append(dist)
            elif ev_type == MutationType.DEL:
                strand_key = "rev" if is_rev else "fwd"
                del_tracker[seq][strand_key] += 1
                del_tracker[seq]["len"] = ev_len
                del_tracker[seq]["distances"].append(dist)

        for seq, data in ins_tracker.items():
            insertions_map[seq] = IndelObservation(
                event_type=MutationType.INS,
                anchor_position=pos,
                sequence=seq,
                length=data["len"],
                anchor_depth=indel_anchor_depth,
                fwd_count=data["fwd"],
                rev_count=data["rev"],
                read_end_distances=tuple(data["distances"]),
            )

        for seq, data in del_tracker.items():
            deletions_map[seq] = IndelObservation(
                event_type=MutationType.DEL,
                anchor_position=pos,
                sequence=seq,
                length=data["len"],
                anchor_depth=indel_anchor_depth,
                fwd_count=data["fwd"],
                rev_count=data["rev"],
                read_end_distances=tuple(data["distances"]),
            )

        resolved_indel_count = len(insertions_map) + len(deletions_map)

        # Strictly separate filtered sequencing errors from collapsed biological events
        excluded_reads = (
            unmapped
            + qc_fail
            + duplicates
            + secondary
            + supplementary
            + refskips
            + no_query_position
            + low_mq
            + low_bq
            + noncanonical_bases
        )

        filter_metrics = PileupFilterMetrics(
            low_bq=low_bq,
            low_mq=low_mq,
            duplicates=duplicates,
            unmapped=unmapped,
            qc_fail=qc_fail,
            secondary=secondary,
            supplementary=supplementary,
            refskips=refskips,
            no_query_position=no_query_position,
            noncanonical_bases=noncanonical_bases,
            structural_deletion_spans=structural_del_spans,
            collapsed_mate_observations=collapsed_mates_count,
            raw_indel_observations=raw_indel_obs_count,
            resolved_indel_observations=resolved_indel_count,
            excluded_reads=excluded_reads,
        )

        return SitePileup(
            position=pos,
            reference_base=ref_base,
            depth_pileup_entries=total_entries,
            depth_filtered=depth_filtered,
            depth_callable=callable_depth,
            indel_anchor_depth=indel_anchor_depth,
            counts=counts,
            forward_counts=fwd_counts,
            reverse_counts=rev_counts,
            mean_bq=mean_bq,
            mean_mq=mean_mq,
            filter_metrics=filter_metrics,
            read_end_distances=dict(read_end_distances),
            insertions=insertions_map,
            deletions=deletions_map,
        )