"""
Block 6 Observational Engine: Gapless 1..N coordinate streaming, single-molecule
mate overlap deduplication, physical fragment boundary awareness, and spatial tracking.
"""

from collections import defaultdict
from typing import Dict, Iterator, List, Optional, Tuple

import pysam

from src.modules.consensus.models import (
    AlleleCounts,
    ConsensusPolicy,
    PileupFilterMetrics,
    SitePileup,
)


class PileupEngine:
    """
    Streams gapless 1-based genomic coordinates across a target contig.
    Resolves overlapping mate pairs into single-molecule observations, distinguishes
    physical fragment boundaries on paired vs merged reads, tracks granular filter
    exclusions, and computes per-allele distance vectors.
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

    def _create_empty_site(self, position: int, ref_base: str) -> SitePileup:
        """Construct a zero-depth observational SitePileup for uncovered loci."""
        return SitePileup(
            position=position,
            reference_base=ref_base,
            depth_pileup_entries=0,
            depth_filtered=0,
            depth_callable=0,
            counts=AlleleCounts(),
            forward_counts=AlleleCounts(),
            reverse_counts=AlleleCounts(),
            mean_bq=0.0,
            mean_mq=0.0,
            filter_metrics=PileupFilterMetrics(
                low_bq=0,
                low_mq=0,
                duplicates=0,
                unmapped=0,
                qc_fail=0,
                secondary=0,
                supplementary=0,
                deletions=0,
                refskips=0,
                no_query_position=0,
                noncanonical_bases=0,
                overlapping_mate=0,
                excluded_reads=0,
            ),
            allele_distances={},
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
                col_pos = col.reference_pos + 1  # 1-based

                while current_pos < col_pos and current_pos <= self.expected_reference_length:
                    ref_b = self._fetch_reference_base(fasta, current_pos)
                    yield self._create_empty_site(current_pos, ref_b)
                    current_pos += 1

                if current_pos > self.expected_reference_length:
                    break

                ref_b = self._fetch_reference_base(fasta, col_pos)
                site_pileup = self._process_column(col, ref_b)
                yield site_pileup
                current_pos += 1

            while current_pos <= self.expected_reference_length:
                ref_b = self._fetch_reference_base(fasta, current_pos)
                yield self._create_empty_site(current_pos, ref_b)
                current_pos += 1

    def _process_column(self, col: pysam.PileupColumn, ref_base: str) -> SitePileup:
        """
        Execute molecule-aware filtering on a single pileup column:
        1. Exclude structural defects and soft-clipped bases.
        2. Group paired reads by template (query_name) to collapse overlaps and penalize discordance.
        3. Determine true physical molecule termini:
           - Merged / Overlapping: both boundaries are true termini (min(d5, d3)).
           - Paired unmerged: only cycle 1 (qpos - q_start) is the true molecule terminus.
        4. Apply Mapping Quality and Base Quality thresholds, preserving depth_filtered vs depth_callable.
        """
        pos = col.reference_pos + 1
        total_entries = len(col.pileups)

        unmapped = 0
        qc_fail = 0
        duplicates = 0
        secondary = 0
        supplementary = 0
        deletions = 0
        refskips = 0
        no_query_position = 0
        overlapping_mate = 0
        low_mq = 0
        low_bq = 0
        noncanonical_bases = 0

        # Pass 1: Structural defect filters & strict soft-clip boundary rejection
        # Elements: (alignment, qpos, q_align_start, q_align_end)
        template_reads: Dict[str, List[Tuple[pysam.AlignedSegment, int, int, int]]] = defaultdict(list)

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
            if read.is_del:
                deletions += 1
                continue
            if read.is_refskip:
                refskips += 1
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

            # Strict Anti-Soft-Clipping: reject bases outside aligned genomic coordinates
            if qpos < q_start or qpos >= q_end:
                no_query_position += 1
                continue

            template_reads[aln.query_name].append((aln, qpos, q_start, q_end))

        # Pass 2: Single-molecule deduplication & Physical Fragment Boundary Calculation
        # Output: (aln, qpos, base, bq, physical_dist_from_end)
        resolved_observations: List[Tuple[pysam.AlignedSegment, int, str, int, int]] = []

        for qname, obs_list in template_reads.items():
            if len(obs_list) == 1:
                aln, qpos, q_start, q_end = obs_list[0]
                base = aln.query_sequence[qpos].upper() if aln.query_sequence else "N"
                bq = aln.query_qualities[qpos] if aln.query_qualities is not None else 0

                if aln.is_paired:
                    # PAIRED UNMERGED READ:
                    # Cycle 1 (5' of read sequence, qpos - q_start) is the only true physical terminus.
                    # The 3' end is an arbitrary sequencer cycle cutoff inside the insert.
                    dist_from_end = max(0, qpos - q_start)
                else:
                    # MERGED OR SINGLE-END READ:
                    # Both ends represent true physical molecule break points.
                    d5 = qpos - q_start
                    d3 = q_end - qpos - 1
                    dist_from_end = max(0, min(d5, d3))

                resolved_observations.append((aln, qpos, base, bq, dist_from_end))

            else:
                # OVERLAPPING MATE PAIR: Represent the same physical molecule
                overlapping_mate += len(obs_list) - 1

                mate_info = []
                for aln, qpos, q_start, q_end in obs_list:
                    b = aln.query_sequence[qpos].upper() if aln.query_sequence else "N"
                    q = aln.query_qualities[qpos] if aln.query_qualities is not None else 0
                    # For overlapping mates, each read's 5' represents an external physical end of the insert
                    d5 = qpos - q_start
                    d3 = q_end - qpos - 1
                    dist = max(0, min(d5, d3))
                    mate_info.append((aln, qpos, b, q, dist))

                # Physical fragment distance: shortest distance to any true molecule terminus
                fragment_dist = min(m[4] for m in mate_info)

                # Single-Molecule Allelic Agreement Resolution
                unique_bases = {m[2] for m in mate_info}
                if len(unique_bases) == 1:
                    # Concordant: pick mate with highest Base Quality
                    best_mate = max(mate_info, key=lambda m: m[3])
                    resolved_observations.append(
                        (best_mate[0], best_mate[1], best_mate[2], best_mate[3], fragment_dist)
                    )
                else:
                    # Discordant: mates disagree on the same physical molecule
                    mate_info.sort(key=lambda m: m[3], reverse=True)
                    best_mate = mate_info[0]
                    second_mate = mate_info[1]

                    # Penalize Base Quality by Phred subtraction (BQ = BQ1 - BQ2)
                    quality_delta = best_mate[3] - second_mate[3]
                    if quality_delta > 0:
                        resolved_observations.append(
                            (best_mate[0], best_mate[1], best_mate[2], quality_delta, fragment_dist)
                        )
                    else:
                        # Exact tie: nullify call (will be excluded during quality filter)
                        resolved_observations.append(
                            (best_mate[0], best_mate[1], "N", 0, fragment_dist)
                        )

        # Pass 3: Evaluate MQ, BQ and distinguish depth_filtered vs depth_callable
        # Qualifying observations: (base, is_reverse, bq, mq, dist_from_end)
        qualifying_bases: List[Tuple[str, bool, int, int, int]] = []

        for aln, qpos, base, bq, dist_from_end in resolved_observations:
            mq = aln.mapping_quality
            if mq < self.policy.min_mapping_quality:
                low_mq += 1
                continue

            if bq < self.policy.min_base_quality:
                low_bq += 1
                continue

            qualifying_bases.append((base, aln.is_reverse, bq, mq, dist_from_end))

        # depth_filtered includes all observations passing MQ and BQ (even non-canonical)
        depth_filtered = len(qualifying_bases)

        # Pass 4: Compile canonical allele counts and spatial distance vectors
        counts_map = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "other": 0}
        fwd_counts_map = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "other": 0}
        rev_counts_map = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "other": 0}
        allele_distances: Dict[str, List[int]] = defaultdict(list)

        total_bq = 0
        total_mq = 0
        callable_depth = 0

        for base, is_rev, bq, mq, dist in qualifying_bases:
            total_bq += bq
            total_mq += mq

            if base in {"A", "C", "G", "T"}:
                callable_depth += 1
                counts_map[base] += 1
                allele_distances[base].append(dist)
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

        excluded_reads = (
            unmapped
            + qc_fail
            + duplicates
            + secondary
            + supplementary
            + deletions
            + refskips
            + no_query_position
            + overlapping_mate
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
            deletions=deletions,
            refskips=refskips,
            no_query_position=no_query_position,
            noncanonical_bases=noncanonical_bases,
            overlapping_mate=overlapping_mate,
            excluded_reads=excluded_reads,
        )

        return SitePileup(
            position=pos,
            reference_base=ref_base,
            depth_pileup_entries=total_entries,
            depth_filtered=depth_filtered,
            depth_callable=callable_depth,
            counts=counts,
            forward_counts=fwd_counts,
            reverse_counts=rev_counts,
            mean_bq=mean_bq,
            mean_mq=mean_mq,
            filter_metrics=filter_metrics,
            allele_distances=dict(allele_distances),
        )