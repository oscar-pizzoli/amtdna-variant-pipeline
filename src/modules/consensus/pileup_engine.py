"""
Block 6 Observational Engine: Gapless 1..N coordinate streaming, single-molecule
mate overlap deduplication, and rigorous read-level QC metric accounting.
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
    Resolves overlapping mate pairs into single-molecule observations, evaluates
    canonical nucleotide observations, and tracks exhaustive filter exclusions.
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
        )

    def iter_sites(self) -> Iterator[SitePileup]:
        """
        Stream gapless observational SitePileup objects from 1 to expected_reference_length.
        Applies mate-pair collapsing and quality filtering with zero silent drops.
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
                col_pos = col.reference_pos + 1  # Convert to 1-based

                # Emit gapless zero-depth sites for uncovered coordinates
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

            # Emit remaining tail positions if BAM coverage ends early
            while current_pos <= self.expected_reference_length:
                ref_b = self._fetch_reference_base(fasta, current_pos)
                yield self._create_empty_site(current_pos, ref_b)
                current_pos += 1

    def _process_column(self, col: pysam.PileupColumn, ref_base: str) -> SitePileup:
        """
        Execute molecule-aware filtering on a single pileup column:
        1. Exclude structural defects (unmapped, QC-fail, secondary, duplicate, deletions).
        2. Group surviving paired reads by template (query_name) to collapse overlaps.
        3. Apply Mapping Quality and Base Quality thresholds to the resolved observations.
        """
        pos = col.reference_pos + 1
        total_entries = len(col.pileups)

        # Exclusion counters
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

        # Pass 1: Filter structural defects & group candidates by template (query_name)
        template_reads: Dict[str, List[Tuple[pysam.AlignedSegment, int]]] = defaultdict(list)

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

            template_reads[aln.query_name].append((aln, read.query_position))

        # Pass 2: Overlapping mate deduplication per template
        resolved_observations: List[Tuple[pysam.AlignedSegment, int]] = []

        for qname, observations in template_reads.items():
            if len(observations) == 1:
                resolved_observations.append(observations[0])
            else:
                overlapping_mate += len(observations) - 1

                # Select observation with the highest base quality
                best_obs = max(
                    observations,
                    key=lambda obs: obs[0].query_qualities[obs[1]] if obs[0].query_qualities is not None else 0,
                )
                resolved_observations.append(best_obs)

        # Pass 3: Evaluate MQ, BQ and tally allelic observations
        callable_bases: List[Tuple[str, bool, int, int]] = []

        for aln, qpos in resolved_observations:
            mq = aln.mapping_quality
            if mq < self.policy.min_mapping_quality:
                low_mq += 1
                continue

            bq = aln.query_qualities[qpos] if aln.query_qualities is not None else 0
            if bq < self.policy.min_base_quality:
                low_bq += 1
                continue

            base = aln.query_sequence[qpos].upper()
            if base not in {"A", "C", "G", "T"}:
                noncanonical_bases += 1
                continue

            callable_bases.append((base, aln.is_reverse, bq, mq))

        # Pass 4: Compile metrics and allele counts
        counts_map = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "other": 0}
        fwd_counts_map = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "other": 0}
        rev_counts_map = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "other": 0}

        total_bq = 0
        total_mq = 0

        for base, is_rev, bq, mq in callable_bases:
            total_bq += bq
            total_mq += mq
            key = base if base in {"A", "C", "G", "T", "N"} else "other"
            counts_map[key] += 1
            if is_rev:
                rev_counts_map[key] += 1
            else:
                fwd_counts_map[key] += 1

        counts = AlleleCounts(**counts_map)
        fwd_counts = AlleleCounts(**fwd_counts_map)
        rev_counts = AlleleCounts(**rev_counts_map)

        callable_depth = len(callable_bases)
        mean_bq = (total_bq / callable_depth) if callable_depth > 0 else 0.0
        mean_mq = (total_mq / callable_depth) if callable_depth > 0 else 0.0

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
            depth_filtered=callable_depth,
            depth_callable=callable_depth,
            counts=counts,
            forward_counts=fwd_counts,
            reverse_counts=rev_counts,
            mean_bq=mean_bq,
            mean_mq=mean_mq,
            filter_metrics=filter_metrics,
        )