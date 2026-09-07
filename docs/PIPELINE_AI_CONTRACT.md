# COMPLETE PIPELINE ARCHITECTURE, TEST CONTRACT & API SURFACE
> **Generated:** 2026-09-01 15:10:04 UTC | **Root:** `amtdna-variant-pipeline`

## 1. BLOCK-BY-BLOCK CONTRACT & INVARIANTS

### Block 1: Reference Genome Management
- **Engine**: `src.modules.alignment.reference_manager.ReferenceManager`
- **Primary Entrypoint**: `prepare_and_validate_reference(fasta_path, expected_contig, strict_contig, expected_length, allow_n, auto_index, log_dir)`
- **Invariants**:
  1. `allow_n=True` is strictly required for rCRS (NC_012920.1) due to position 3107.
  2. Validates physical FAI and 5 BWA indices (`.bwt`, `.pac`, `.ann`, `.amb`, `.sa`).
  3. Emits validated `ReferenceValidationResult` containing canonical `ReferenceInfo`.

### Block 2: FASTQ Discovery, Integrity & Preprocessing
- **Engines**:
  - Discovery: `src.modules.preprocessing.discovery.SampleDiscoveryEngine`
  - Integrity: `src.modules.preprocessing.integrity.IntegrityEngine`
  - Preprocessing: `src.modules.preprocessing.trimmer.PreprocessingEngine`
- **Primary Entrypoints**:
  - `SampleDiscoveryEngine.discover_files(input_dir)`
  - `SampleDiscoveryEngine.pair_and_validate_fastq(discovered_files, allow_missing_mates, single_end_mode)`
  - `SampleDiscoveryEngine.build_sample_records(paired_units, orphan_files)`
  - `IntegrityEngine.verify_sample_integrity(sample)`
  - `PreprocessingEngine.process_sample(sample, config)`
- **Invariants**:
  1. AdapterRemoval trims adapters, merges overlapping reads, and isolates branches (`collapsed`, `pair1`/`pair2`, `single`).
  2. Verifies SHA-256 fingerprint on raw inputs.

### Block 3: Alignment, Deduplication, Filtering & QC
- **Orchestrator**: `src.modules.alignment.orchestrator.AlignmentOrchestrator`
- **Primary Entrypoint**: `run_sample_alignment(sample, reference_info, output_dir, log_dir, policy, threads, ...)`
- **Invariants**:
  1. Multi-branch alignment with distinct Read Groups (`@RG LB:PAIRED`, `@RG LB:MERGED`, `@RG LB:SINGLE`).
  2. Canonical ordering: Merge branches -> Samtools collate/fixmate -m -> Coordinate sort -> Markdup -> MAPQ filtering.
  3. Outputs `filtered.dedup.bam` (MAPQ >= 30, flags 0xF04 filtered) and indexes.

### Block 4: aDNA Damage Authentication & Profiling
- **Orchestrator**: `src.modules.damage.orchestrator.DamageOrchestrator`
- **Primary Entrypoint**: `run_sample_damage_pipeline(sample_id, reference_info, output_dir, log_dir, policy, dedup_bam, ...)`
- **Invariants**:
  1. Executes mapDamage2 statistical quantification (5' C->T and 3' G->A).
  2. Generates base-quality rescaled BAM (`rescaled.bam`) and PMD score metrics.
  3. Downstream Block 5 MUST receive the unrescaled deduplicated BAM.

### Block 5: Two-Tier Contamination Triage & Bayesian Modeling
- **Batch Runner**: `src.modules.contamination.batch_runner.ContaminationBatchRunner`
- **Primary Entrypoint**: `run_cohort_pipeline(discovered_samples, reference_fasta, reference_contig, reference_length_bp, output_dir, config, resource_alloc)`
- **Invariants**:
  1. FastContaminationScreen evaluates 354 diagnostic sites concurrently.
  2. `FAST_CLEAR` skips Schmutzi Bayesian modeling completely.
  3. `ESCALATE` and `INSUFFICIENT_DATA` queue bounded execution of Schmutzi in isolated containers.
  4. Global host CPU budget fraction is strictly enforced across all steps.

---

## 2. STRICT OPERATIONAL RULES (DO NOT LIST)
1. **DO NOT** instantiate ad-hoc command wrappers; invoke official domain classes directly.
2. **DO NOT** pass `reference_info` to `prepare_and_validate_reference(fasta_path=...)`.
3. **DO NOT** pass rescaled BAMs to Schmutzi in Block 5.
4. **DO NOT** bypass Schmutzi if FastScreen returns `INSUFFICIENT_DATA`.
5. **DO NOT** open `.sai` alignment output as text streams (handle raw binary only).
6. **DO NOT** mark a stage as SUCCESS without verifying physical file existence and non-zero size on disk.

---

## 3. EXACT CODEBASE API SURFACE (CLASSES, METHODS & SIGNATURES)

### Module `src.modules.alignment.bwa_aligner`
**File**: `src/modules/alignment/bwa_aligner.py`

#### `class BWAAlignerError(Exception)`
> Raised when BWA alignment, streaming execution, or BAM validation fails.

#### `class AlignmentBranchResult`
> Outcome and traceability artifacts for an individual alignment branch.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes branch alignment outcome."""
```

#### `class AlignmentBranchResults`
> Container aggregating all executed alignment branches for a sample.

```python
def all_bams(self) -> List[Path]:
    """Returns physical paths to all successfully generated branch BAM files."""
```

```python
def all_read_groups(self) -> List[ReadGroupInfo]:
    """Returns all injected ReadGroupInfo objects."""
```

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes aggregated branch results."""
```

#### `class BWAAligner`
> Executes conservative BWA-aln passes and streams SAM output directly into compressed BAM files.

```python
@classmethod
def get_bwa_version(cls) -> Tuple[str, str]:
    """Inspects and returns the installed BWA version and raw banner."""
```

```python
@classmethod
def get_samtools_version(cls) -> Tuple[str, str]:
    """Inspects and returns the installed samtools version and raw banner."""
```

```python
@classmethod
def align_single_end_branch(cls, fastq_path: Path, source: ReadGroupSource, sample_id: str, library_id: str, reference_info: ReferenceInfo, output_dir: Path, log_dir: Path, cached_count: Optional[int] = None, threads: int = 4, seed_length: int = 1000, max_diff: float = 0.01, max_gap_opens: int = 2, max_gap_extensions: int = -1, trim_quality_threshold: int = 0, platform: str = 'ILLUMINA', platform_unit: Optional[str] = None, sequencing_center: Optional[str] = None, keep_temporary_sai: bool = False) -> AlignmentBranchResult:
    """Executes full alignment lifecycle for a Single-End branch (Merged, Merged-Truncated, or Singleton)."""
```

```python
@classmethod
def align_paired_end_branch(cls, r1_fastq: Path, r2_fastq: Path, sample_id: str, library_id: str, reference_info: ReferenceInfo, output_dir: Path, log_dir: Path, cached_pairs_count: Optional[int] = None, threads: int = 4, seed_length: int = 1000, max_diff: float = 0.01, max_gap_opens: int = 2, max_gap_extensions: int = -1, trim_quality_threshold: int = 0, platform: str = 'ILLUMINA', platform_unit: Optional[str] = None, sequencing_center: Optional[str] = None, keep_temporary_sai: bool = False, parallel_aln: bool = True) -> AlignmentBranchResult:
    """Executes full alignment lifecycle for intact unmerged Paired-End mates using parallel SA indexing."""
```

```python
@classmethod
def align_sample_branches(cls, sample: SampleRecord, reference_info: ReferenceInfo, output_dir: Path, log_dir: Path, threads: int = 4, seed_length: int = 1000, max_diff: float = 0.01, max_gap_opens: int = 2, max_gap_extensions: int = -1, trim_quality_threshold: int = 0, platform: str = 'ILLUMINA', platform_unit: Optional[str] = None, sequencing_center: Optional[str] = None, keep_temporary_sai: bool = False, parallel_aln: bool = True) -> AlignmentBranchResults:
    """Dispatches all active read branches (merged, merged_truncated, singleton, paired) for a SampleRecord."""
```

---

### Module `src.modules.alignment.dedup_engine`
**File**: `src/modules/alignment/dedup_engine.py`

#### `class DedupEngineError(Exception)`
> Raised when deduplication execution, verification, or parsing fails.

#### `class DedupExecutionResult`
> Outcome and artifacts resulting from duplicate marking and physical deduplication.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes deduplication outcome."""
```

#### `class DedupEngine`
> Manages samtools markdup execution, dual-BAM generation (.markdup vs .dedup), indexing, and metrics parsing.

```python
@classmethod
def index_bam(cls, bam_path: Path, threads: int = 2, use_csi: bool = False) -> Path:
    """Generates coordinate index (.bai or .csi) for a target BAM file."""
```

```python
@classmethod
def verify_bam_integrity(cls, bam_path: Path) -> Tuple[bool, int]:
    """Performs quickcheck integrity validation and returns (is_empty, record_count)."""
```

```python
@classmethod
def verify_dedup_bam(cls, dedup_bam: Path, markdup_record_count: int) -> int:
    """Validates that deduplicated BAM has <= records than markdup BAM and contains zero 0x400 duplicate flags."""
```

```python
@classmethod
def parse_markdup_stats(cls, stats_content_or_file: Union[Path, str]) -> Optional[DedupMetrics]:
    """Parses key-value statistical metrics generated by 'samtools markdup -s / -f'."""
```

```python
@classmethod
def execute_dedup(cls, sample_id: str, input_sorted_bam: Path, output_dir: Path, log_dir: Path, dedup_mode: DedupMode = DedupMode.BOTH, threads: int = 4, use_csi: bool = False, cleanup_on_failure: bool = True) -> DedupExecutionResult:
    """Executes 'samtools markdup -s' to flag duplicates and produces dual BAM outputs (.markdup and/or .dedup)."""
```

---

### Module `src.modules.alignment.filter_engine`
**File**: `src/modules/alignment/filter_engine.py`

#### `class FilterEngineError(Exception)`
> Raised when alignment filtering, verification, or MAPQ profiling fails.

#### `class ExclusionReason(str, Enum)`
> Mutually exclusive categories for read filtering.

#### `class FilterExclusionStats`
> Exclusive first-failure read counts summing to total input reads.

```python
def total_removed(self) -> int:
```

```python
def retention_rate(self) -> float:
```

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes exclusion statistics."""
```

#### `class FilterExecutionResult`
> Outcome, artifacts, and stratified metrics resulting from BAM filtering.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes filtering outcome."""
```

#### `class AlignmentFilterEngine`
> Applies standardized SAM bitwise flag filters, MAPQ thresholds, and generates stratified quality profiles.

```python
@classmethod
def index_bam(cls, bam_path: Path, threads: int = 2, use_csi: bool = False) -> Path:
    """Generates coordinate index (.bai or .csi) for a target BAM file."""
```

```python
@classmethod
def verify_bam_readability(cls, bam_path: Path) -> int:
    """Checks BAM readability, header validity, and quick structural integrity via samtools quickcheck."""
```

```python
@classmethod
def verify_filtered_bam_postconditions(cls, filtered_bam: Path, policy: AlignmentFilterPolicy) -> None:
    """Asserts that output filtered BAM contains zero records violating active filter policies."""
```

```python
@classmethod
def compute_mapq_profiles_and_exclusions(cls, input_bam: Path, policy: AlignmentFilterPolicy) -> Tuple[MAPQProfile, MAPQProfile, FilterExclusionStats]:
    """Scans input BAM to compute raw profile, retained profile, and exclusive first-failure statistics."""
```

```python
@classmethod
def filter_bam_file(cls, input_bam: Path, output_bam: Path, policy: AlignmentFilterPolicy, threads: int, log_dir: Path) -> str:
    """Applies SAM bitwise flag and MAPQ filters using samtools view."""
```

```python
@classmethod
def execute_filtering(cls, sample_id: str, input_markdup_bam: Optional[Path], input_dedup_bam: Optional[Path], output_dir: Path, log_dir: Path, policy: Optional[AlignmentFilterPolicy] = None, threads: int = 4, fast_mode: bool = False, use_csi: bool = False, cleanup_on_failure: bool = True) -> FilterExecutionResult:
    """Applies configured quality filters to marked and/or deduplicated BAM files."""
```

---

### Module `src.modules.alignment.metrics_collector`
**File**: `src/modules/alignment/metrics_collector.py`

#### `class MetricsCollectorError(Exception)`
> Raised when QC metrics extraction, depth calculation, or flagstat validation fails.

#### `class FlagstatMetrics`
> Structured container holding parsed samtools flagstat counts and validation invariants.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes flagstat metrics."""
```

#### `class MetricsCollectionResult`
> Top-level container grouping all intermediate flagstat reports, depth profiles, and consolidated metrics.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes complete metrics collection result."""
```

#### `class MetricsCollector`
> Extracts alignment statistics, depth percentiles, and coverage breadth across pipeline stages.

```python
@classmethod
def parse_flagstat_output(cls, flagstat_content_or_file: Union[Path, str]) -> FlagstatMetrics:
    """Parses samtools flagstat standard text output into a strongly typed FlagstatMetrics object."""
```

```python
@classmethod
def run_flagstat(cls, bam_path: Path, output_txt: Path, threads: int = 2) -> Tuple[FlagstatMetrics, str]:
    """Executes 'samtools flagstat' and saves the report file."""
```

```python
@classmethod
def compute_depth_distribution(cls, bam_path: Path, reference_info: ReferenceInfo, output_depth_tsv: Optional[Path] = None, min_base_quality: int = 0, min_mapping_quality: int = 0) -> Tuple[List[int], float, float, int, float, float]:
    """Calculates per-base mitochondrial depth distribution across all positions using 'samtools depth -aa -r contig'."""
```

```python
@classmethod
def calculate_coverage_metrics(cls, final_bam: Path, reference_info: ReferenceInfo, baseline_bam: Optional[Path] = None, output_depth_tsv: Optional[Path] = None, min_base_quality: int = 0, min_mapping_quality: int = 0) -> CoverageMetrics:
    """Computes comprehensive depth percentiles and breadth of coverage percentages for the target reference."""
```

```python
@classmethod
def extract_branch_mapped_counts(cls, branch_results: AlignmentBranchResults) -> Dict[ReadGroupSource, int]:
    """Inspects branch BAMs using pysam to count mapped records per branch."""
```

```python
@classmethod
def validate_alignment_metrics_invariants(cls, metrics: AlignmentMetrics) -> None:
    """Asserts mathematical consistency across all final alignment metrics."""
```

```python
@classmethod
def assemble_alignment_metrics(cls, branch_results: AlignmentBranchResults, filter_result: FilterExecutionResult, dedup_result: DedupExecutionResult, coverage_metrics: CoverageMetrics, raw_flagstat: Optional[FlagstatMetrics] = None, scan_branch_bams: bool = False) -> AlignmentMetrics:
    """Consolidates metrics from all alignment stages into the domain AlignmentMetrics schema."""
```

```python
@classmethod
def collect_all_metrics(cls, sample_id: str, reference_info: ReferenceInfo, branch_results: AlignmentBranchResults, dedup_result: DedupExecutionResult, filter_result: FilterExecutionResult, output_dir: Path, log_dir: Path, prepared_bam: Optional[Path] = None, threads: int = 2, scan_branch_bams: bool = False, min_base_quality: int = 0, min_mapping_quality: int = 0) -> MetricsCollectionResult:
    """Orchestrates comprehensive flagstat extraction across stages, depth calculation, and metric assembly."""
```

---

### Module `src.modules.alignment.models`
**File**: `src/modules/alignment/models.py`

#### `class AlignerType(str, Enum)`
> Supported alignment engines.

#### `class DedupMode(str, Enum)`
> PCR duplicate handling strategies.

#### `class AlignmentStatus(str, Enum)`
> Top-level execution and validation status.

#### `class AlignmentWarningCode(str, Enum)`
> Specific diagnostic warning codes.

#### `class AlignmentErrorCode(str, Enum)`
> Specific diagnostic failure codes.

#### `class ReadGroupSource(str, Enum)`
> Biological and technical classification of input reads within a Read Group.

#### `class MultiMappingStrategy(str, Enum)`
> Strategy used to identify and filter ambiguous multi-mapping reads.

#### `class MultiMappingPolicy`
> Explicit parameters governing multi-mapping and ambiguous read resolution.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes multi-mapping policy."""
```

#### `class AlignmentFilterPolicy`
> Complete filtering policy combining multi-mapping rules, deduplication, and QC thresholds.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes complete filter policy."""
```

#### `class ReadGroupInfo`
> Represents a validated SAM/BAM Read Group header (@RG).

```python
def __post_init__(self) -> None:
    """Defensive validation ensuring compliant SAM header specifications."""
```

```python
def to_bwa_rg_string(self) -> str:
    """Formats the Read Group string for the BWA -r parameter."""
```

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes Read Group information."""
```

#### `class ReferenceInfo`
> Mitochondrial reference genome metadata dynamically resolved against input FASTA.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes reference metadata."""
```

#### `class MAPQProfile`
> Stratified Mapping Quality (MAPQ) distribution.

```python
def total_mapped_evaluated(self) -> int:
```

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes MAPQ distribution."""
```

#### `class CoverageMetrics`
> Mitochondrial depth, coverage distribution, and breadth percentages [0.0 - 100.0%].

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes coverage metrics."""
```

#### `class DedupMetrics`
> PCR duplicate metrics with explicit read, pair, and single-fragment separation.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes deduplication metrics."""
```

#### `class ProvenanceMetrics`
> Tracks origin and class distribution of input units and their mapped counterparts.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes read provenance metrics."""
```

#### `class AlignmentOutputs`
> Container mapping physical BAM, BAI, and metric output files.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes physical file paths."""
```

#### `class AlignmentMetrics`
> Comprehensive metric summary cleanly separating raw, filtered, and deduplicated counts.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes all alignment metrics."""
```

#### `class AlignmentToolMetadata`
> Audit trail recording environment, versions, hashes, and executed commands.

```python
def mark_started(self) -> None:
    """Records UTC start timestamp."""
```

```python
def mark_finished(self) -> None:
    """Records UTC finish timestamp and calculates runtime."""
```

```python
def add_command(self, cmd: str) -> None:
    """Appends an executed shell command."""
```

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes audit metadata."""
```

#### `class AlignmentResult`
> Top-level record summarizing the complete outcome of Block 3 for a sample.

```python
def add_warning(self, msg: str, code: Optional[AlignmentWarningCode] = None) -> None:
    """Records a human-readable warning and optional structured warning code."""
```

```python
def add_error(self, msg: str, code: AlignmentErrorCode) -> None:
    """Records an error message, sets status strictly to FAIL, and logs error code."""
```

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes the complete alignment result."""
```

---

### Module `src.modules.alignment.orchestrator`
**File**: `src/modules/alignment/orchestrator.py`

#### `class AlignmentOrchestratorError(Exception)`
> Raised when an unrecoverable orchestration failure occurs in Block 3.

#### `class AlignmentOrchestrator`
> Orchestrates end-to-end alignment, deduplication, filtering, and metric collection for ancient and modern mtDNA samples.

```python
@classmethod
def evaluate_qc_thresholds(cls, metrics: AlignmentMetrics, policy: AlignmentFilterPolicy, result: AlignmentResult) -> None:
    """Evaluates final alignment and coverage metrics against configured QC thresholds and assigns diagnostic warnings."""
```

```python
@classmethod
def run_sample_alignment(cls, sample: SampleRecord, reference_info: ReferenceInfo, output_dir: Path, log_dir: Path, policy: Optional[AlignmentFilterPolicy] = None, threads: int = 4, parallel_aln: bool = True, keep_temporary_sai: bool = False, keep_raw_merged_bam: bool = True, use_csi: bool = False, fast_mode: bool = False, scan_branch_bams: bool = False, seed_length: int = 1000, max_diff: float = 0.01, max_gap_opens: int = 2, max_gap_extensions: int = -1, trim_quality_threshold: int = 0, platform: str = 'ILLUMINA', platform_unit: Optional[str] = None, sequencing_center: Optional[str] = None, cleanup_on_failure: bool = True) -> AlignmentResult:
    """Executes the complete Block 3 pipeline for a single sample."""
```

---

### Module `src.modules.alignment.reference_manager`
**File**: `src/modules/alignment/reference_manager.py`

#### `class ReferenceManagerError(Exception)`
> Raised strictly on unrecoverable I/O or external binary execution failures.

#### `class ContigAliasResolver`
> Standardizes mitochondrial contig nomenclature across RefSeq, UCSC, and Ensembl conventions.

```python
@classmethod
def are_aliases(cls, contig_a: str, contig_b: str) -> bool:
    """Checks if two contig names represent equivalent naming conventions."""
```

```python
@classmethod
def resolve_canonical(cls, contig_name: str) -> str:
    """Resolves a contig name to its canonical RefSeq identifier if registered."""
```

#### `class ReferenceValidationResult`
> Summary of comprehensive integrity, sequence, and indexing checks for a reference genome.

```python
def to_dict(self) -> Dict[str, object]:
    """Serializes validation summary."""
```

#### `class ReferenceManager`
> Manages validation, cryptographic auditing, cross-consistency, and indexing of the reference FASTA.

```python
@classmethod
def compute_file_sha256(cls, file_path: Path) -> str:
    """Computes the SHA-256 cryptographic hash across the physical reference file on disk."""
```

```python
@classmethod
def parse_and_hash_sequence(cls, file_path: Path, allow_n: bool = False) -> Tuple[str, str, Optional[str], int, List[str], List[str]]:
    """Parses FASTA reference, validates single-contig structure, and computes biological sequence MD5."""
```

```python
@classmethod
def check_bwa_index(cls, fasta_path: Path) -> Tuple[bool, Optional[str]]:
    """Checks if all required BWA index files exist, are non-empty, and are not older than the FASTA."""
```

```python
@classmethod
def parse_and_validate_fai(cls, fasta_path: Path, expected_contig: str, expected_length: int) -> Tuple[bool, Optional[str]]:
    """Validates that samtools faidx index (.fai) exists, has 1 record, and matches contig and length exactly."""
```

```python
@classmethod
def generate_bwa_index(cls, fasta_path: Path, log_dir: Optional[Path] = None) -> None:
    """Executes 'bwa index' to generate Burrows-Wheeler Transform index files."""
```

```python
@classmethod
def generate_fai_index(cls, fasta_path: Path, log_dir: Optional[Path] = None) -> None:
    """Executes 'samtools faidx' to create sequence dictionary index."""
```

```python
@classmethod
def prepare_and_validate_reference(cls, fasta_path: Path, expected_contig: str = 'NC_012920.1', strict_contig: bool = True, expected_length: Optional[int] = 16569, expected_sequence_md5: Optional[str] = None, is_circular: bool = True, allow_n: bool = False, auto_index: bool = True, log_dir: Optional[Path] = None) -> ReferenceValidationResult:
    """Orchestrates reference verification, cryptographic auditing, alias resolution, indexing, and cross-consistency checks."""
```

---

### Module `src.modules.alignment.samtools_pipeline`
**File**: `src/modules/alignment/samtools_pipeline.py`

#### `class SamtoolsPipelineError(Exception)`
> Raised when BAM merging, streaming execution, sorting, indexing, or verification fails.

#### `class PipelinePreparationResult`
> Output artifacts and audit records resulting from structural BAM preparation.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes preparation outcome."""
```

#### `class SamtoolsPipeline`
> Orchestrates BAM branch merging, UNIX pipe collate -> fixmate -m -> sort, indexing, and baseline stats.

```python
@classmethod
def index_bam(cls, bam_path: Path, threads: int = 2) -> Path:
    """Generates coordinate index (.bai) for a target BAM file."""
```

```python
@classmethod
def generate_samtools_stats(cls, bam_path: Path, output_stats_txt: Path, threads: int = 2) -> str:
    """Runs 'samtools stats' on the prepared BAM to generate comprehensive baseline metrics."""
```

```python
@classmethod
def verify_prepared_bam(cls, bam_path: Path, expected_read_groups: List[ReadGroupInfo], expected_contig: Optional[str] = None) -> int:
    """Performs fast integrity/truncation check, validates header consistency (SO:coordinate, @RG, @SQ), and returns total records."""
```

```python
@classmethod
def merge_branch_bams(cls, branch_bams: List[Path], output_merged_bam: Path, threads: int, log_dir: Path) -> str:
    """Merges multiple raw branch BAM files into a single unified BAM preserving all distinct Read Groups."""
```

```python
@classmethod
def run_collate_fixmate_sort_stream(cls, input_bam: Path, output_sorted_bam: Path, threads: int, log_dir: Path) -> List[str]:
    """Executes streamed UNIX pipe: 'collate -u -O' -> 'fixmate -m -u' -> 'sort -@ sort_threads' without intermediate disk files."""
```

```python
@classmethod
def process_sample_branches(cls, branch_results: AlignmentBranchResults, output_dir: Path, log_dir: Path, threads: int = 4, keep_raw_merged_bam: bool = True, cleanup_on_failure: bool = True, expected_contig: Optional[str] = 'NC_012920.1') -> PipelinePreparationResult:
    """Dispatches full BAM merge, streaming fixmate sync, coordinate sort, indexing, stats, and validation."""
```

---

### Module `src.modules.contamination.batch_runner`
**File**: `src/modules/contamination/batch_runner.py`

#### `class ResourceAllocation`
> Quantitative allocation of system threads and multi-phase worker limits.

```python
def to_dict(self) -> Dict[str, Any]:
```

#### `class BatchSampleRecord`
> Summary record for a single processed sample in the cohort.

```python
def to_tsv_row(self) -> str:
```

#### `class ContaminationBatchRunner`
> Coordinates two-phase cohort contamination estimation and strict resource management.

```python
@classmethod
def load_yaml_config(cls, config_path: Optional[Path] = None) -> Dict[str, Any]:
    """Loads configuration from YAML file."""
```

```python
@classmethod
def resolve_path(cls, path_str: Optional[str], default_rel: Path) -> Path:
    """Resolves path string relative to PROJECT_ROOT."""
```

```python
@classmethod
def detect_available_cores(cls, host_cores_override: Optional[int] = None) -> int:
    """Determines host CPU cores respecting affinity, containers, or explicit overrides."""
```

```python
@classmethod
def calculate_resource_allocation(cls, threads_per_sample: int = 4, cpu_budget_fraction: float = 0.5, explicit_max_workers: Optional[int] = None, host_cores_override: Optional[int] = None, allow_oversubscription: bool = False) -> ResourceAllocation:
    """Calculates safe concurrent worker limits strictly bounding host core allocation."""
```

```python
@classmethod
def discover_bam_inputs(cls, input_dir: Path, pattern: str = '*.bam', sample_manifest: Optional[Path] = None) -> List[Tuple[str, Path]]:
    """Discovers input BAM files via structured manifest or defensive directory scanning."""
```

```python
@classmethod
def run_cohort_pipeline(cls, discovered_samples: List[Tuple[str, Path]], reference_fasta: Path, reference_contig: str, reference_length_bp: int, output_dir: Path, config: Dict[str, Any], resource_alloc: ResourceAllocation) -> List[ContaminationResult]:
    """Orchestrates two-phase cohort execution maximizing concurrency and work conservation."""
```

```python
@classmethod
def export_cohort_summary(cls, results: List[ContaminationResult], output_dir: Path, resource_alloc: ResourceAllocation) -> Tuple[Path, Path]:
    """Generates audit-ready cohort-level TSV and JSON summaries with performance analytics."""
```

```python
def main() -> None:
    """CLI interface for running two-phase multi-sample cohort contamination estimation."""
```

---

### Module `src.modules.contamination.dissent_engine`
**File**: `src/modules/contamination/dissent_engine.py`

#### `class DissentExecutionError(Exception)`
> Raised when BAM reading, pileup generation, or dissent parsing encounters an unrecoverable failure.

#### `class DissentExecutionMetadata`
> Audit metadata tracking dissent calculation parameters, filters, and execution runtimes.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes dissent execution metadata to a dictionary."""
```

#### `class DissentEngine`
> Calculates non-Bayesian allelic dissent and heterogeneity screening metrics across mtDNA sites.

```python
@classmethod
def is_masked_position(cls, pos_1based: int, masked_regions: Set[Tuple[int, int]]) -> bool:
    """Checks if a 1-based genomic position falls within masked homopolymer/artifact regions."""
```

```python
@classmethod
def calculate_effective_dissent_threshold(cls, num_evaluated: int, is_targeted: bool, base_min_dissenting: int = 3, genome_wide_dissent_rate_floor: float = 0.001) -> int:
    """Scales the required dissenting site count based on whether the scan is targeted or genome-wide."""
```

```python
@classmethod
def evaluate_model_sufficiency(cls, num_evaluated: int, is_targeted: bool, target_set_size: Optional[int] = None, min_genome_wide_positions: int = 1000, min_targeted_positions: int = 10) -> bool:
    """Evaluates statistical validity and sufficiency of the evaluated position set."""
```

```python
@classmethod
def calculate_robust_point_estimate(cls, vafs_all: np.ndarray, vafs_dissenting: np.ndarray, effective_min_dissenting: int) -> float:
    """Calculates a continuous point estimate without step-function discontinuities."""
```

```python
@classmethod
def execute_dissent_screening(cls, sample_id: str, input_bam: Path, reference_fasta: Path, output_dir: Path, log_dir: Path, expected_contig: Optional[str] = None, min_base_qual: int = 30, min_mapq: int = 30, min_depth_per_site: int = 10, noise_vaf_floor: float = 0.015, base_min_dissenting_sites: int = 3, target_positions: Optional[Set[int]] = None, custom_masked_regions: Optional[Set[Tuple[int, int]]] = None, sam_flag_filter: int = DEFAULT_SAM_FLAG_FILTER) -> Tuple[ContaminationEstimate, ModelQualityDiagnostics, ContaminationOutputFiles, DissentExecutionMetadata]:
    """Scans BAM pileup with strict flag filters and homopolymer masking to quantify allelic heterogeneity."""
```

---

### Module `src.modules.contamination.fast_screen`
**File**: `src/modules/contamination/fast_screen.py`

#### `class FastScreenExecutionError(Exception)`
> Raised when BAM parsing, pileup iteration, or screening metrics calculation fails.

#### `class FastContaminationScreen`
> Performs rapid allelic dissent profiling on diagnostic positions and conservative triage.

```python
@classmethod
def sanitize_sample_id(cls, sample_id: str) -> str:
    """Sanitizes sample identifiers for safe filesystem and filename operations."""
```

```python
@classmethod
def compute_sha256(cls, file_path: Path) -> Optional[str]:
    """Calculates SHA-256 digest of a target file."""
```

```python
@classmethod
def compute_canonical_panel_hash(cls, panel_data: Union[Dict[int, Dict[str, str]], Set[int]]) -> str:
    """Computes a deterministic SHA-256 hash from in-memory positions and panel alleles."""
```

```python
@classmethod
def load_diagnostic_panel(cls, panel_path: Optional[Path]) -> Tuple[Dict[int, Dict[str, str]], Optional[str]]:
    """Loads 1-based diagnostic positions, reference bases, and known alternative alleles from panel TSV."""
```

```python
@classmethod
def compute_genome_mean_depth(cls, bam_path: Path, contig_name: str, contig_length: int) -> float:
    """Calculates true overall sequencing depth across the entire mitochondrial contig."""
```

```python
@classmethod
def load_cached_result(cls, cache_path: Path, current_bam: Path, current_bam_sha: Optional[str], current_ref_sha: Optional[str], current_panel_sha: Optional[str], current_policy_sha: str, algorithm_version: str, cache_mode: str = 'strict') -> Optional[FastScreenResult]:
    """Loads screening result from cache after strictly verifying provenance fingerprints."""
```

```python
@classmethod
def save_cache_result(cls, cache_path: Path, result: FastScreenResult, input_bam: Path) -> None:
    """Serializes screening result to disk cache with filesystem metadata."""
```

```python
@classmethod
def run_screen(cls, sample_id: str, input_bam: Path, reference_fasta: Path, expected_contig: str, policy: FastScreeningPolicy, output_dir: Path, target_positions: Optional[Union[Dict[int, Dict[str, str]], Set[int]]] = None) -> Tuple[FastScreenResult, Optional[Path]]:
    """Executes rapid position-by-position allelic dissent screening against diagnostic sites."""
```

---

### Module `src.modules.contamination.input_preparer`
**File**: `src/modules/contamination/input_preparer.py`

#### `class SchmutziPreparationError(Exception)`
> Raised when BAM validation, subsampling, calmd MD tag injection, or reference wrapping fails.

#### `class SchmutziReferenceAdapter`
> Adapts mitochondrial reference genomes for circularity awareness in schmutzi.

```python
@classmethod
def compute_sha256(cls, file_path: Path) -> str:
    """Computes SHA-256 digest of a target file."""
```

```python
@classmethod
def prepare_wrapped_reference(cls, reference_info: ReferenceInfo, work_dir: Path, padding_bp: int = 200) -> Tuple[Path, WrappingMetadata]:
    """Generates a wrapped (circularly extended) FASTA reference with flanking padding."""
```

#### `class SchmutziInputPreparer`
> Prepares and strictly validates BAM files to satisfy schmutzi's input requirements.

```python
@classmethod
def run_samtools_quickcheck(cls, bam_path: Path) -> bool:
    """Runs 'samtools quickcheck' to verify file completeness and intact EOF blocks."""
```

```python
@classmethod
def ensure_bai_index(cls, bam_path: Path, reference_contig: Optional[str] = None, threads: int = 4) -> Path:
    """Ensures a valid, queryable BAI index exists."""
```

```python
@classmethod
def compute_contig_mean_depth(cls, bam_path: Path, contig_name: str, contig_length: int) -> float:
    """Calculates exact mean sequencing depth across target contig."""
```

```python
@classmethod
def subsample_bam(cls, input_bam: Path, output_bam: Path, target_fraction: float, seed: int = 42, threads: int = 4) -> None:
    """Executes reproducible pseudo-random subsampling using samtools view -s."""
```

```python
@classmethod
def sample_md_tags_present(cls, bam_path: Path, sample_size: int = 100) -> bool:
    """Inspects the first N mapped records to verify if MD tags are present."""
```

```python
@classmethod
def validate_bam_structure(cls, bam_path: Path, reference_info: ReferenceInfo, allow_non_mt_contigs: bool = False) -> None:
    """Validates BAM header, @SQ sequence length, and mtDNA-only contig constraints."""
```

```python
@classmethod
def run_calmd_atomic(cls, input_bam: Path, reference_fasta: Path, output_bam: Path, threads: int = 4) -> Path:
    """Executes 'samtools calmd -b' writing atomically to avoid corrupted partial outputs."""
```

```python
@classmethod
def prepare_bam_for_schmutzi(cls, sample_id: str, input_bam: Path, reference_info: ReferenceInfo, work_dir: Path, subsampling_policy: Optional[SubsamplingPolicy] = None, threads: int = 4, force_calmd: bool = True, allow_non_mt_contigs: bool = False) -> Tuple[Path, Path, bool, float, Optional[float]]:
    """Prepares BAM for schmutzi."""
```

---

### Module `src.modules.contamination.models`
**File**: `src/modules/contamination/models.py`

#### `class ContaminationMethod(str, Enum)`
> Methodological approach selected for contamination estimation or triage.

#### `class ScreeningDecision(str, Enum)`
> Conservative triage decisions emitted by the FastContaminationScreen engine.

#### `class ContaminationStatus(str, Enum)`
> Overall execution status of Block 5.

#### `class ContaminationCategory(str, Enum)`
> Configured operational categorization of the estimated contamination fraction.

#### `class EstimationConfidence(str, Enum)`
> Statistical confidence level of the estimation based on depth, CI width, and convergence.

#### `class ConsensusGatingAction(str, Enum)`
> Recommended action and calling constraints passed downstream to Block 6 (Consensus & Variant Calling).

#### `class ContaminationWarningCode(str, Enum)`
> Diagnostic warning codes emitted during contamination analysis.

#### `class ContaminationErrorCode(str, Enum)`
> Fatal error codes halting downstream processing.

#### `class SubsamplingPolicy`
> Configuration policy for conditional depth subsampling prior to contamination estimation.

```python
def to_dict(self) -> Dict[str, Any]:
```

#### `class FastScreeningPolicy`
> Execution parameters and escalation criteria for the fast pre-screening triage engine.

```python
def compute_policy_sha256(self) -> str:
```

```python
def to_dict(self) -> Dict[str, Any]:
```

#### `class FastScreenResult`
> Quantitative summary and conservative triage decision produced by FastContaminationScreen.

```python
def validate_invariants(self) -> bool:
```

```python
def to_dict(self) -> Dict[str, Any]:
```

#### `class ContaminationInputContext`
> Unified input context container passed to Block 5 orchestration.

```python
def to_dict(self) -> Dict[str, Any]:
```

#### `class WrappingMetadata`
> Metadata tracking circular mitochondrial origin padding and coordinate transformations.

```python
def to_dict(self) -> Dict[str, Any]:
```

#### `class ContaminationEstimate`
> Quantitative contamination metrics and statistical intervals (Level 1).

```python
def __post_init__(self) -> None:
```

```python
def ci_width(self) -> Optional[float]:
```

```python
def to_dict(self) -> Dict[str, Any]:
```

#### `class ModelQualityDiagnostics`
> Model quality, convergence, depth metrics, CI width, and risk indicators (Level 2).

```python
def to_dict(self) -> Dict[str, Any]:
```

#### `class GatingRecommendation`
> Categorical interpretation and formal constraints for downstream Block 6 execution (Level 3).

```python
def to_dict(self) -> Dict[str, Any]:
```

#### `class ContaminationOutputFiles`
> Physical filesystem paths of Block 5 deliverables.

```python
def to_dict(self) -> Dict[str, Any]:
```

#### `class ContaminationAuditMetadata`
> Audit trail tracking input BAM provenance, reference hashes, tool versions, and commands.

```python
def to_dict(self) -> Dict[str, Any]:
```

#### `class SchmutziPolicy`
> Execution parameters specific to schmutzi Bayesian mixture modeling.

```python
def to_dict(self) -> Dict[str, Any]:
```

#### `class DissentPolicy`
> Execution parameters specific to Track B allelic dissent and heterogeneity screening.

```python
def to_dict(self) -> Dict[str, Any]:
```

#### `class ContaminationPolicy`
> Configurable execution parameters, heuristic boundaries, and gating criteria for Block 5.

```python
def __post_init__(self) -> None:
```

```python
def to_dict(self) -> Dict[str, Any]:
```

#### `class ContaminationResult`
> Top-level deliverable container encapsulating Block 5 execution, quality, audit, and gating.

```python
def mark_started(self) -> None:
```

```python
def mark_finished(self) -> None:
```

```python
def add_warning(self, code: ContaminationWarningCode, message: str) -> None:
```

```python
def add_error(self, code: ContaminationErrorCode, message: str) -> None:
```

```python
def to_dict(self) -> Dict[str, Any]:
```

---

### Module `src.modules.contamination.orchestrator`
**File**: `src/modules/contamination/orchestrator.py`

#### `class ContaminationOrchestrationError(Exception)`
> Raised when contamination workflow execution, pre-flight validation, or pipeline stages fail fatally.

#### `class ContaminationOrchestrator`
> Orchestrates mtDNA contamination estimation, pre-screening triage, and downstream gating.

```python
def __init__(self, config: Optional[Dict[str, Any]] = None, config_path: Optional[Path] = None, project_root: Optional[Path] = None) -> None:
```

```python
def process_sample(self, context: ContaminationInputContext, output_dir: Path) -> ContaminationResult:
    """Executes complete contamination analysis adhering to two-tier screening architecture."""
```

---

### Module `src.modules.contamination.panel_builder`
**File**: `src/modules/contamination/panel_builder.py`

#### `class PanelBuilderValidationError(Exception)`
> Raised when the diagnostic panel fails post-construction invariant validation.

#### `class DiagnosticSite`
> Represents a validated polymorphic diagnostic site.

#### `class DiagnosticPanelBuilder`
> Parses schmutzi frequency tables and builds validated allele-aware screening panels.

```python
@classmethod
def load_yaml_config(cls, config_path: Optional[Path] = None) -> Dict[str, Any]:
    """Loads default parameters from pipeline_config.yaml if present."""
```

```python
@classmethod
def resolve_path(cls, path_str: Optional[str], default_relative: Path) -> Path:
    """Resolves path string relative to PROJECT_ROOT if not already absolute."""
```

```python
@classmethod
def compute_sha256(cls, file_path: Path) -> str:
    """Calculates SHA-256 checksum of a single file."""
```

```python
@classmethod
def compute_directory_sha256(cls, dir_path: Path) -> str:
    """Calculates combined SHA-256 checksum of all files in a directory."""
```

```python
@classmethod
def parse_schmutzi_freq_directory(cls, freqs_dir: Path, strict_mode: bool = True) -> Tuple[Dict[int, Dict[str, float]], Dict[int, float], int]:
    """Aggregates allele distributions across frequency tables, handling counts vs normalized frequencies."""
```

```python
@classmethod
def calculate_informativeness(cls, ref_base: str, allele_counts: Dict[str, float], total_obs: float) -> Tuple[float, float, str]:
    """Calculates Minor Allele Frequency (MAF) and non-reference informativeness score."""
```

```python
@classmethod
def validate_panel(cls, diagnostic_sites: List[DiagnosticSite], genome_length: int, min_panel_size: int = 10) -> None:
    """Enforces structural, positional, and biological invariants on the generated panel."""
```

```python
@classmethod
def build_panel(cls, freqs_dir: Path, reference_fasta: Path, expected_contig: str = 'NC_012920.1', min_maf_threshold: float = 0.02, min_total_observations: float = 10.0, min_panel_size: int = 10, strict_mode: bool = True) -> Tuple[List[DiagnosticSite], str, str]:
    """Builds and validates the allele-aware candidate diagnostic panel."""
```

```python
@classmethod
def export_panel_tsv(cls, diagnostic_sites: List[DiagnosticSite], output_tsv_path: Path, reference_sha: str, database_sha: str, min_maf: float, min_total_observations: float, version_tag: str = 'phylotree17_schmutzi197_v1.0') -> Path:
    """Serializes diagnostic sites to an audit-ready, self-describing TSV file."""
```

```python
def main() -> None:
    """CLI entrypoint resolving configurations from YAML and command-line flags."""
```

---

### Module `src.modules.contamination.qc_evaluator`
**File**: `src/modules/contamination/qc_evaluator.py`

#### `class ContaminationQCValidator`
> Evaluates contamination model validity, assigns confidence tiers, and formulates Block 6 gating directives.

```python
@classmethod
def evaluate_model_confidence(cls, estimate: ContaminationEstimate, diagnostics: ModelQualityDiagnostics, policy: ContaminationPolicy) -> Tuple[EstimationConfidence, List[ContaminationWarningCode], List[str]]:
    """Evaluates statistical reliability based on method-specific quality metrics."""
```

```python
@classmethod
def categorize_contamination(cls, estimate: ContaminationEstimate, confidence: EstimationConfidence, diagnostics: ModelQualityDiagnostics, policy: ContaminationPolicy) -> Tuple[ContaminationCategory, List[ContaminationWarningCode], List[str]]:
    """Classifies contamination rate incorporating point estimate and generalized Credible Interval bounds."""
```

```python
@classmethod
def formulate_gating_recommendation(cls, category: ContaminationCategory, confidence: EstimationConfidence, estimate: ContaminationEstimate, outputs: ContaminationOutputFiles, policy: ContaminationPolicy) -> GatingRecommendation:
    """Formulates actionable gating constraints for Block 6 following a strict safety hierarchy."""
```

```python
@classmethod
def determine_execution_status(cls, error_codes: List[ContaminationErrorCode], warning_codes: List[ContaminationWarningCode], confidence: EstimationConfidence, category: ContaminationCategory, gating_action: ConsensusGatingAction) -> ContaminationStatus:
    """Determines the overall execution status enforcing strict coherence with gating actions."""
```

```python
@classmethod
def serialize_deliverable_json(cls, result: ContaminationResult, output_json_dir: Path) -> Optional[Path]:
    """Serializes the full ContaminationResult into a structured JSON deliverable."""
```

```python
@classmethod
def evaluate_sample(cls, sample_id: str, estimate: ContaminationEstimate, diagnostics: ModelQualityDiagnostics, outputs: ContaminationOutputFiles, audit: ContaminationAuditMetadata, policy: ContaminationPolicy, fallback_reason: Optional[str] = None, output_json_dir: Optional[Path] = None) -> ContaminationResult:
    """Orchestrates quality evaluation, gating formulation, and JSON deliverable export."""
```

---

### Module `src.modules.contamination.schmutzi_engine`
**File**: `src/modules/contamination/schmutzi_engine.py`

#### `class SchmutziExecutionError(Exception)`
> Raised when schmutzi binaries fail, timeout, or output artifacts cannot be parsed.

#### `class SchmutziExecutionMetadata`
> Audit metadata tracking schmutzi command execution, versions, and file hashes.

```python
def to_dict(self) -> Dict[str, Any]:
```

#### `class SchmutziOutput`
> Standardized output bundle wrapping Schmutzi execution results.

```python
def cont_est_file(self) -> Optional[Path]:
```

```python
def cont_deam_file(self) -> Optional[Path]:
```

```python
def endogenous_fasta(self) -> Optional[Path]:
```

```python
def contaminant_fasta(self) -> Optional[Path]:
```

#### `class SchmutziEngine`
> Wrapper and parser for schmutzi Bayesian mtDNA contamination estimation.

```python
@classmethod
def is_docker_available(cls) -> bool:
    """Checks if docker CLI is available and the schmutzi image exists."""
```

```python
@classmethod
def get_binary_name(cls) -> Optional[str]:
    """Detects if schmutzi is available via Docker image or local binary in PATH."""
```

```python
@classmethod
def get_schmutzi_version(cls) -> str:
    """Retrieves installed schmutzi version string."""
```

```python
@classmethod
def compute_sha256(cls, target_path: Path) -> Optional[str]:
    """Computes SHA-256 digest of a target file or directory tree."""
```

```python
@classmethod
def resolve_frequency_database_dir(cls, database_path: Path) -> Path:
    """Resolves the folder containing the .freq allele frequency files."""
```

```python
@classmethod
def compute_bam_mean_depth(cls, bam_path: Path, expected_contig: Optional[str] = None) -> float:
    """Computes mean sequencing depth across mitochondrial genome."""
```

```python
@classmethod
def parse_cont_est_file(cls, cont_est_path: Path) -> Tuple[float, float, float, Optional[float]]:
    """Parses point estimate, 95% CI (lower, upper), and log-likelihood from .cont.est output."""
```

```python
@classmethod
def parse_cont_deam_file(cls, cont_deam_path: Path) -> Optional[float]:
    """Parses preliminary contamination rate estimated from deamination alone (.contDeam.est)."""
```

```python
@classmethod
def parse_schmutzi_log(cls, log_path: Path) -> Tuple[bool, int, bool, bool]:
    """Parses execution log for strict MCMC convergence, iterations, and damage warnings."""
```

```python
@classmethod
def parse_contaminant_haplogroups(cls, final_mtcont_path: Optional[Path], final_cont_path: Optional[Path]) -> List[str]:
    """Extracts the top predicted contaminant database identifiers from mtcont or final.cont."""
```

```python
@classmethod
def count_informative_positions_from_logs(cls, endo_log: Optional[Path], cont_log: Optional[Path]) -> int:
    """Counts positions where endogenous and contaminant consensus calls diverge."""
```

```python
@classmethod
def run_schmutzi_pipeline(cls, sample_id: str, schmutzi_ready_bam: Path, reference_fasta: Path, contaminant_database: Path, output_dir: Path, log_dir: Path, policy: SchmutziPolicy, wrapping_metadata: Optional[WrappingMetadata] = None, expected_contig: Optional[str] = None, timeout_seconds: int = 7200) -> Tuple[ContaminationEstimate, ModelQualityDiagnostics, ContaminationOutputFiles, SchmutziExecutionMetadata]:
    """Executes full two-phase schmutzi workflow via Docker container or native binaries."""
```

```python
@classmethod
def run(cls, sample_id: str, input_bam: Path, reference_fasta: Path, database_dir: Path, output_dir: Path, threads: int = 4, max_iterations: int = 10, use_wrapped: bool = False, library_treatment: LibraryTreatment = LibraryTreatment.NONE) -> SchmutziOutput:
    """High-level execution adapter matching ContaminationOrchestrator requirements."""
```

---

### Module `src.modules.damage.damage_quantifier`
**File**: `src/modules/damage/damage_quantifier.py`

#### `class DamageQuantifierError(Exception)`
> Raised when mapDamage execution, pre-flight validation, or output parsing encounters an unrecoverable failure.

#### `class DamageQuantifierEngine`
> Orchestrates mapDamage2 quantification and parses empirical mismatch matrices and Bayesian stats.

```python
@classmethod
def get_binary_name(cls) -> Optional[str]:
    """Detects if mapDamage is registered in system PATH."""
```

```python
@classmethod
def get_mapdamage_version(cls) -> Tuple[str, str]:
    """Retrieves installed mapDamage version string."""
```

```python
@classmethod
def validate_bam_reference_compatibility(cls, bam_path: Path, ref_fasta: Path) -> None:
    """Validates that contigs and lengths declared in BAM header strictly match the reference FASTA."""
```

```python
@classmethod
def parse_misincorporation_file(cls, misinc_path: Path, window_length: int = 25) -> Tuple[TerminalDamageFrequencies, List[str]]:
    """Parses mapDamage2 misincorporation.txt using explicit C>T and G>A misincorporation columns."""
```

```python
@classmethod
def parse_length_distribution_file(cls, lg_dist_path: Path) -> FragmentLengthMetrics:
    """Parses lgdistribution.txt supporting both 2-column and 3-column (Std/Length/Occurences) formats."""
```

```python
@classmethod
def parse_mcmc_summary_file(cls, summ_stat_path: Path) -> BayesianDamageParameters:
    """Parses MCMC summary statistics from Stats_out_MCMC_iter_summ_stat.csv with real convergence checks."""
```

```python
@classmethod
def run_mapdamage_stats(cls, sample_id: str, input_bam: Path, reference_fasta: Path, output_dir: Path, log_dir: Path, policy: DamagePolicy, threads: int = 2) -> Tuple[DamageOutputFiles, TerminalDamageFrequencies, FragmentLengthMetrics, BayesianDamageParameters]:
    """Executes mapDamage2 in statistical estimation mode (without rescaling), collecting tables and MCMC models."""
```

---

### Module `src.modules.damage.models`
**File**: `src/modules/damage/models.py`

#### `class LibraryTreatment(str, Enum)`
> Experimental enzymatic library treatment protocol.

#### `class DamageInputBAMType(str, Enum)`
> Choice of upstream BAM file fed into damage estimation engines.

#### `class DamageSignalLevel(str, Enum)`
> Categorical assessment of observable empirical deamination signal.

#### `class AuthenticationVerdict(str, Enum)`
> Protocol-aware authenticity diagnosis for ancient DNA authentication.

#### `class DamageStatus(str, Enum)`
> Overall execution status of Block 4.

#### `class DamageWarningCode(str, Enum)`
> Diagnostic warning codes for Block 4 quality control.

#### `class DamageErrorCode(str, Enum)`
> Fatal error codes halting downstream processing.

#### `class TerminalDamageFrequencies`
> Empirical nucleotide substitution rates and supporting observation counts at terminal positions.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes terminal mismatch rates and counts to a dictionary."""
```

#### `class FragmentLengthMetrics`
> Summary statistics of aligned fragment length distribution.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes length distribution statistics to a dictionary."""
```

#### `class BayesianDamageParameters`
> Statistical parameters estimated via mapDamage2 statistical model.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes Bayesian model parameters to a dictionary."""
```

#### `class PMDMetrics`
> Read-level Post-Mortem Damage (PMD) score distribution and filtering statistics.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes PMD metrics to a dictionary."""
```

#### `class DamageOutputFiles`
> Physical filesystem paths of Block 4 deliverables.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes output paths to a dictionary."""
```

#### `class DamagePolicy`
> Configurable execution parameters, heuristic threshold limits, and quality criteria for Block 4.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes policy configuration to a dictionary."""
```

#### `class DamageModelProfile`
> Consolidated statistical and biological damage profile for downstream consumption.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes damage model profile to a structured dictionary."""
```

#### `class DamageResult`
> Top-level container encapsulating execution results and QC diagnostics for Block 4.

```python
def __post_init__(self) -> None:
    """Aligns sample_id across internal profile container."""
```

```python
def mark_started(self) -> None:
    """Records ISO start timestamp."""
```

```python
def mark_finished(self) -> None:
    """Records ISO finish timestamp and computes duration."""
```

```python
def add_warning(self, code: DamageWarningCode, message: str) -> None:
    """Appends a non-fatal warning code and message."""
```

```python
def add_error(self, code: DamageErrorCode, message: str) -> None:
    """Appends a fatal error code and transitions status to FAIL."""
```

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes full Block 4 execution result into a structured dictionary."""
```

---

### Module `src.modules.damage.orchestrator`
**File**: `src/modules/damage/orchestrator.py`

#### `class DamageOrchestrationError(Exception)`
> Raised when an unrecoverable failure occurs during Block 4 orchestration.

#### `class DamageBAMValidationError(Exception)`
> Raised when pre-flight checks on input BAM fail.

#### `class DamageAuditTrail`
> Comprehensive execution metadata consolidating tool versions, commands, and sub-engine metrics.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes orchestration audit trail to a dictionary."""
```

#### `class DamageOrchestrator`
> End-to-end coordinator for Block 4 damage profiling, quality rescaling, and authenticity testing.

```python
@classmethod
def validate_input_bam_preflight(cls, bam_path: Path, expected_contig: str) -> int:
    """Executes full pre-flight integrity, header compatibility, and mapped read count checks on input BAM."""
```

```python
@classmethod
def resolve_input_bam(cls, sample_id: str, policy: DamagePolicy, dedup_bam: Optional[Path] = None, markdup_bam: Optional[Path] = None, custom_bam: Optional[Path] = None, allow_fallback: bool = False) -> Path:
    """Strictly resolves target input BAM path according to configured DamageInputBAMType policy."""
```

```python
@classmethod
def run_sample_damage_pipeline(cls, sample_id: str, reference_info: ReferenceInfo, output_dir: Path, log_dir: Path, policy: Optional[DamagePolicy] = None, dedup_bam: Optional[Path] = None, markdup_bam: Optional[Path] = None, custom_input_bam: Optional[Path] = None, threads: int = 2, allow_input_bam_fallback: bool = False, require_pmd_success: bool = False) -> DamageResult:
    """Executes full Block 4 workflow: mapDamage estimation -> Rescaling -> PMD scoring -> QC Evaluation."""
```

---

### Module `src.modules.damage.pmd_engine`
**File**: `src/modules/damage/pmd_engine.py`

#### `class PMDScoringError(Exception)`
> Raised when PMD scoring, MD-tag injection, or threshold filtering encounters a domain-specific failure.

#### `class PMDExecutionMetadata`
> Audit container tracking tool parameters, MD injection, and execution times.

```python
def add_command(self, cmd: str) -> None:
    """Appends an executed shell command to the audit trail."""
```

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes PMD execution metadata to a dictionary."""
```

#### `class PMDScoringEngine`
> Orchestrates per-read ancient damage likelihood scoring and optional filtering via PMDtools.

```python
@classmethod
def get_binary_command(cls) -> List[str]:
    """Resolves whether pmdtools is installed as a standalone binary or Python script in PATH."""
```

```python
@classmethod
def get_version(cls) -> str:
    """Extracts installed pmdtools version string."""
```

```python
@classmethod
def verify_and_ensure_md_tags(cls, bam_path: Path, reference_fasta: Path, work_dir: Path, sample_check_limit: int = 1000) -> Tuple[Path, bool]:
    """Verifies uniform presence of SAM MD tags across sampled mapped reads, injecting them via calmd if incomplete."""
```

```python
@classmethod
def parse_pmd_scores_table(cls, scores_tsv: Path) -> PMDMetrics:
    """Parses pmdtools per-read TSV output and computes empirical distribution metrics with parsing diagnostics."""
```

```python
@classmethod
def execute_pmd_pipeline(cls, sample_id: str, input_bam: Path, reference_fasta: Path, output_dir: Path, log_dir: Path, policy: DamagePolicy, threads: int = 2, timeout: Optional[int] = 3600, cleanup_on_failure: bool = True) -> Tuple[Optional[Path], Optional[Path], Optional[Path], PMDMetrics, PMDExecutionMetadata]:
    """Executes full PMD scoring stream, parses score distribution, and optionally generates filtered BAM."""
```

---

### Module `src.modules.damage.qc_evaluator`
**File**: `src/modules/damage/qc_evaluator.py`

#### `class DamageQCValidator`
> Evaluates multi-criteria damage evidence and assigns protocol-aware authenticity verdicts.

```python
@classmethod
def evaluate_bayesian_convergence(cls, bayesian: BayesianDamageParameters) -> Tuple[bool, List[DamageWarningCode], List[str]]:
    """Audits mapDamage MCMC Bayesian parameter convergence prior to downstream usage."""
```

```python
@classmethod
def evaluate_damage_signal_level(cls, terminal: TerminalDamageFrequencies, bayesian: BayesianDamageParameters, policy: DamagePolicy, bayesian_usable: bool) -> Tuple[DamageSignalLevel, List[DamageWarningCode], List[str]]:
    """Determines empirical deamination signal strength supported by rigorous statistical denominators."""
```

```python
@classmethod
def determine_authenticity_verdict(cls, signal_level: DamageSignalLevel, terminal: TerminalDamageFrequencies, length_metrics: FragmentLengthMetrics, pmd_metrics: Optional[PMDMetrics], policy: DamagePolicy) -> Tuple[AuthenticationVerdict, List[DamageWarningCode], List[str]]:
    """Assigns protocol-aware authenticity verdict across all experimental library types."""
```

```python
@classmethod
def evaluate_sample(cls, sample_id: str, input_bam_path: Path, reference_contig: str, terminal_damage: TerminalDamageFrequencies, fragment_length: FragmentLengthMetrics, bayesian_params: BayesianDamageParameters, pmd_metrics: Optional[PMDMetrics], outputs: DamageOutputFiles, policy: DamagePolicy, rescaling_applied: bool = False, output_json_dir: Optional[Path] = None) -> DamageResult:
    """Synthesizes Block 4 outputs into DamageModelProfile, DamageResult, and exports damage_model.json."""
```

---

### Module `src.modules.damage.rescale_engine`
**File**: `src/modules/damage/rescale_engine.py`

#### `class DamageRescaleError(Exception)`
> Raised when base quality rescaling, model pre-conditions, or output post-conditions fail.

#### `class RescaleExecutionMetadata`
> Audit trail container capturing tool version, command invocation, and validation parameters.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes rescaling execution metadata to a dictionary."""
```

#### `class QualityRescaleEngine`
> Orchestrates base quality score down-scaling on terminal transitions with comprehensive integrity checks.

```python
@classmethod
def get_binary_name(cls) -> Optional[str]:
    """Detects if mapDamage is registered in system PATH."""
```

```python
@classmethod
def get_version(cls) -> str:
    """Retrieves installed mapDamage version string."""
```

```python
@classmethod
def check_samtools_quickcheck(cls, bam_path: Path) -> bool:
    """Executes samtools quickcheck if available in PATH."""
```

```python
@classmethod
def check_bgzf_eof(cls, bam_path: Path) -> bool:
    """Verifies physical presence of the standard 28-byte BGZF EOF block."""
```

```python
@classmethod
def validate_preconditions(cls, input_bam: Path, reference_fasta: Path, damage_model_dir: Path, expected_contig: Optional[str] = 'NC_012920.1') -> None:
    """Verifies presence and validity of input BAM, FASTA, and mapDamage statistical model files."""
```

```python
@classmethod
def validate_postconditions(cls, input_bam: Path, rescaled_bam: Path, rescaled_bai: Path, expected_contig: Optional[str] = 'NC_012920.1') -> int:
    """Executes multi-tier post-rescaling integrity, header preservation, record parity, and BAI checks."""
```

```python
@classmethod
def run_rescaling(cls, sample_id: str, input_bam: Path, reference_fasta: Path, damage_model_dir: Path, output_dir: Path, log_dir: Path, policy: DamagePolicy, expected_contig: Optional[str] = 'NC_012920.1') -> Tuple[Optional[Path], Optional[Path], RescaleExecutionMetadata]:
    """Executes mapDamage --rescale-only, standardizes deliverables, and enforces post-rescaling audit."""
```

---

### Module `src.modules.preprocessing.bam_validator`
**File**: `src/modules/preprocessing/bam_validator.py`

#### `class BAMValidationResult`
> Comprehensive diagnostic evaluation of a BAM alignment file.

```python
def add_failure(self, reason: str) -> None:
    """Registers a critical failure, strictly transitioning status to FAIL."""
```

```python
def add_warning(self, reason: str) -> None:
    """Registers a non-critical anomaly, setting status to WARN if not already in FAIL."""
```

```python
def mark_passed(self) -> None:
    """Sets status to PASS if no failures or warnings were recorded."""
```

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes validation diagnostics to a structured dictionary."""
```

#### `class BAMValidationEngine`
> Multi-tier validator for BAM file integrity, indexation, contig resolution, and alignment metrics.

```python
@classmethod
def check_bgzf_eof(cls, bam_path: Path) -> bool:
    """Verifies presence of the standard 28-byte BGZF EOF block to detect truncation."""
```

```python
@staticmethod
def locate_index(bam_path: Path) -> Tuple[bool, Optional[str], Optional[Path]]:
    """Checks for standard .bai or .csi indices in both naming conventions."""
```

```python
@classmethod
def resolve_mt_contig(cls, bam_references: List[str], expected_contig: Optional[str] = None, auto_detect: bool = True) -> Optional[str]:
    """Resolves target mitochondrial contig in BAM using exact match or authorized alias list."""
```

```python
@classmethod
def verify_sorting_order(cls, bam: pysam.AlignmentFile, target_contig: str, sample_reads: int = 1000) -> bool:
    """Verifies actual coordinate monotonicity for the target contig by sampling records."""
```

```python
@classmethod
def validate_and_inspect_bam(cls, sample: SampleRecord, expected_contig: str = 'NC_012920.1', reference_fasta: Optional[Path] = None, auto_detect_contig: bool = True, auto_index: bool = True) -> BAMValidationResult:
    """Executes full diagnostic validation on input BAM."""
```

---

### Module `src.modules.preprocessing.discovery`
**File**: `src/modules/preprocessing/discovery.py`

#### `class PairingKey`
> Hashable unique identifier for paired-end mate matching.

#### `class DiscoveredFile`
> Intermediate metadata container derived from tokenized filename inspection.

```python
def pairing_key(self) -> PairingKey:
    """Generates the pairing key for mate matching."""
```

#### `class SequencingUnit`
> Represents an atomic sequencing unit (e.g., a specific lane or library) for a sample.

#### `class FilenameTokenizer`
> Deconstructs filenames from right-to-left into semantic sequencing metadata tokens.

```python
@classmethod
def parse(cls, file_path: Path) -> DiscoveredFile:
    """Parses a file path into a structured DiscoveredFile descriptor."""
```

#### `class SampleDiscoveryEngine`
> Orchestrates sample discovery, strict pairing validation, and pre-flight reporting.

```python
@classmethod
def discover_files(cls, input_dir: Path) -> List[DiscoveredFile]:
    """Scans input directory for all supported genomic files."""
```

```python
@classmethod
def pair_and_validate_fastq(cls, discovered_files: List[DiscoveredFile], allow_missing_mates: bool = False, single_end_mode: bool = False) -> Tuple[Dict[str, List[SequencingUnit]], List[str], List[str]]:
    """Executes deterministic mate pairing with orphan and conflict detection."""
```

```python
@classmethod
def print_preflight_report(cls, sample_records: List[SampleRecord]) -> None:
    """Prints a structured ASCII Pre-Flight Discovery Summary table to stdout and logger."""
```

```python
@classmethod
def build_sample_records(cls, input_dir: Path, config, integrity_mode: str = 'fast') -> List[SampleRecord]:
    """Main orchestrator for sample discovery, validation, and SampleRecord construction."""
```

---

### Module `src.modules.preprocessing.fastq_validator`
**File**: `src/modules/preprocessing/fastq_validator.py`

#### `class FastqValidationError(Exception)`
> Raised when FASTQ structural checks or mate synchronization fail.

#### `class FastqValidator`
> Validates 4-line record structure, base-quality length parity, and PE mate synchronization.

```python
@staticmethod
def validate_fastq_structure(file_path: Path, max_records_to_inspect: Optional[int] = 10000) -> Dict[str, Any]:
    """Inspects structural integrity of a FASTQ file without permissive character decoding."""
```

```python
@classmethod
def validate_paired_end_sync(cls, r1_path: Path, r2_path: Path, sample_check_records: int = 5000) -> bool:
    """Validates that R1 and R2 have equal read counts and matched read identifiers."""
```

---

### Module `src.modules.preprocessing.integrity`
**File**: `src/modules/preprocessing/integrity.py`

#### `class IntegrityError(Exception)`
> Raised when file integrity, existence, or hash computation fails.

#### `class IntegrityEngine`
> Computes SHA-256 hashes, detects file formats, and extracts file metadata.

```python
@staticmethod
def compute_sha256(file_path: Path) -> str:
```

```python
@staticmethod
def detect_format_and_compression(file_path: Path) -> Tuple[str, str]:
```

```python
@classmethod
def fast_count_fastq_reads(cls, file_path: Path) -> int:
```

```python
@classmethod
def generate_fingerprint(cls, file_path: Path, mode: str = 'fast', is_input: bool = True) -> FileFingerprint:
    """Generates a complete metadata and provenance fingerprint for a file."""
```

---

### Module `src.modules.preprocessing.models`
**File**: `src/modules/preprocessing/models.py`

#### `class DataType(str, Enum)`
> Supported input sequencing data types.

#### `class ValidationStatus(str, Enum)`
> Sample processing and quality validation status.

#### `class PreprocessingTool(str, Enum)`
> Supported adapter trimming and merging engines.

#### `class OutputProfile(str, Enum)`
> Classifies the composition of reads surviving preprocessing.

#### `class SampleIdentity`
> Represents the biological and experimental provenance of a sample.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes identity fields to dictionary."""
```

#### `class FileFingerprint`
> Stores cryptographic integrity, file size, and read count for input/output files.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes fingerprint metadata to dictionary."""
```

#### `class InputFiles`
> Strongly validated container for input file paths with strict consistency checks.

```python
def __post_init__(self) -> None:
    """Enforces domain consistency on input paths."""
```

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes input paths to dictionary."""
```

#### `class ProcessedReads`
> Granular container for post-trimming and merged sequence file paths.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes processed read paths to dictionary."""
```

#### `class ToolMetadata`
> Captures command-line invocations, raw version strings, and execution parameters.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes tool metadata to dictionary."""
```

#### `class PreprocessingMetrics`
> Detailed read counts and merging statistics with rigorous metric separation.

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes metrics to dictionary."""
```

#### `class ValidationResult`
> Validation diagnostics and audit issues for the sample.

```python
def add_warning(self, msg: str) -> None:
    """Adds a warning and updates status if not in FAIL."""
```

```python
def add_error(self, msg: str) -> None:
    """Adds an error and transitions status strictly to FAIL."""
```

```python
def to_dict(self) -> Dict[str, Any]:
    """Serializes validation result to dictionary."""
```

#### `class SampleRecord`
> Top-level composable entity tracking a sample across all pipeline blocks.

```python
def sample_id(self) -> str:
```

```python
def data_type(self) -> DataType:
```

```python
def status(self) -> ValidationStatus:
```

```python
def mark_started(self) -> None:
    """Records ISO start timestamp."""
```

```python
def mark_finished(self) -> None:
    """Records ISO finish timestamp and calculates duration."""
```

```python
def to_dict(self) -> Dict[str, Any]:
    """Recursively serializes entire sample entity into a structured JSON/TSV dictionary."""
```

---

### Module `src.modules.preprocessing.orchestrator`
**File**: `src/modules/preprocessing/orchestrator.py`

#### `class RawSampleInput`
> Convenience descriptor for paired-end or single-end FASTQ sample input.

#### `class PreprocessingOrchestrator`
> Orchestrates FASTQ ingestion, validation, and adapter trimming/merging via PreprocessingEngine.

```python
@classmethod
def run_sample(cls, sample_input: RawSampleInput, output_dir: Path, log_dir: Path, min_length: int = 30, min_overlap: int = 11, collapse_reads: bool = True, trim_ns: bool = True, threads: int = 4, integrity_mode: str = 'fast') -> SampleRecord:
    """Executes full sample preprocessing, metrics extraction, and output verification."""
```

---

### Module `src.modules.preprocessing.trimmer`
**File**: `src/modules/preprocessing/trimmer.py`

#### `class PreprocessingEngine`
> Orchestrates conservative FASTQ adapter trimming and read merging using AdapterRemoval.

```python
@classmethod
def get_binary_name(cls) -> Optional[str]:
    """Detects if AdapterRemoval is registered as uppercase or lowercase in system PATH."""
```

```python
@classmethod
def get_adapterremoval_version(cls) -> Tuple[str, str]:
```

```python
@classmethod
def parse_adapterremoval_settings(cls, settings_file: Path, is_paired: bool, warnings_list: List[str]) -> PreprocessingMetrics:
    """Parses AdapterRemoval .settings output file into structured metrics with multi-version resilience."""
```

```python
@classmethod
def run_adapterremoval(cls, sample: SampleRecord, output_dir: Path, log_dir: Path, min_read_length: int, min_merge_overlap: int, collapse_reads: bool, trim_ns: bool, threads: int, integrity_mode: str = 'fast') -> ProcessedReads:
    """Executes AdapterRemoval with conservative parameters, dedicated logging, and physical integrity checks."""
```

```python
@classmethod
def process_sample(cls, sample: SampleRecord, config) -> SampleRecord:
    """Dispatches sample preprocessing extracting configuration parameters."""
```

---

### Module `src.utils.config_loader`
**File**: `src/utils/config_loader.py`

#### `class PipelineConfig`
> Strongly typed container holding parsed and validated pipeline parameters.

```python
@classmethod
def from_yaml(cls, config_path: Path) -> 'PipelineConfig':
    """Parses a YAML configuration file with type enforcement and semantic validation."""
```

---

### Module `src.utils.logger`
**File**: `src/utils/logger.py`

#### `class CustomColorFormatter(logging.Formatter)`
> Logging Formatter to add ANSI colors for terminal output.

```python
def format(self, record: logging.LogRecord) -> str:
```

```python
def setup_logger(name: str = 'mtDNA_Pipeline', log_file: Optional[Path] = None, level: int = logging.INFO) -> logging.Logger:
    """Configures and returns a dual-channel logger instance."""
```

---

### Module `src.utils.validator`
**File**: `src/utils/validator.py`

#### `class EnvironmentValidator`
> Validates external binaries, reference indices, and dataset inputs.

```python
@classmethod
def check_binaries(cls, check_adna: bool = True) -> None:
    """Verifies that mandatory command-line tools are available in PATH."""
```

```python
@staticmethod
def validate_reference(reference_fasta: Path, auto_index: bool = True) -> None:
    """Validates FASTA reference existence and verifies/generates required indices."""
```

---
