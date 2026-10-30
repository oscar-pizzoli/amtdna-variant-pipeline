"""Comprehensive unit and integration test suite for Block 2 (Preprocessing & Ingestion Engine)."""

import gzip
import logging
from pathlib import Path
import shutil
import sys
import tempfile
from typing import Any, Dict
import unittest
from dataclasses import dataclass, field

import pysam

# Risoluzione automatica del PYTHONPATH per importare il modulo src
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.modules.preprocessing import (
    BAMValidationEngine,
    BAMValidationResult,
    DataType,
    DiscoveredFile,
    FastqValidationError,
    FastqValidator,
    FilenameTokenizer,
    InputFiles,
    IntegrityEngine,
    IntegrityError,
    OutputProfile,
    PairingKey,
    PreprocessingEngine,
    PreprocessingMetrics,
    PreprocessingTool,
    ProcessedReads,
    SampleDiscoveryEngine,
    SampleIdentity,
    SampleRecord,
    SequencingUnit,
    ValidationStatus,
)

# Silenzia il logging standard durante l'esecuzione dei test
logging.disable(logging.CRITICAL)


@dataclass
class MockPipelineConfig:
    """Mock container mimicking PipelineConfig for testing purposes."""

    results_dir: Path
    threads: int = 2
    min_read_length: int = 20
    min_merge_overlap: int = 11
    collapse_reads: bool = True
    trim_ns: bool = True
    integrity_mode: str = "fast"
    single_end: bool = False
    allow_missing_mates: bool = False
    pattern_r1: str = "*_R1*.fastq.gz"
    pattern_r2: str = "*_R2*.fastq.gz"
    data_type: str = "FASTQ"
    raw_data: Dict[str, Any] = field(default_factory=dict)


class TestBlock2Models(unittest.TestCase):
    """Tests domain models and defensive validation in dataclasses."""

    def test_input_files_valid_pe_fastq(self):
        inputs = InputFiles(
            data_type=DataType.FASTQ,
            is_paired=True,
            r1_path=Path("sample_R1.fq.gz"),
            r2_path=Path("sample_R2.fq.gz"),
        )
        self.assertEqual(inputs.data_type, DataType.FASTQ)
        self.assertTrue(inputs.is_paired)

    def test_input_files_inconsistent_pe_without_r2_raises(self):
        with self.assertRaises(ValueError):
            InputFiles(
                data_type=DataType.FASTQ,
                is_paired=True,
                r1_path=Path("sample_R1.fq.gz"),
                r2_path=None,
            )

    def test_input_files_inconsistent_bam_with_fastq_raises(self):
        with self.assertRaises(ValueError):
            InputFiles(
                data_type=DataType.BAM,
                is_paired=False,
                bam_path=Path("sample.bam"),
                r1_path=Path("sample_R1.fq.gz"),
            )

    def test_sample_record_lifecycle_and_duration(self):
        rec = SampleRecord(
            identity=SampleIdentity(sample_id="BSP71", library_id="L1"),
            inputs=InputFiles(
                data_type=DataType.FASTQ,
                is_paired=False,
                r1_path=Path("BSP71_R1.fq.gz"),
            ),
        )
        rec.mark_started()
        rec.mark_finished()
        self.assertIsNotNone(rec.processing_started)
        self.assertIsNotNone(rec.processing_finished)
        self.assertGreaterEqual(rec.runtime_seconds, 0.0)


class TestIntegrityEngine(unittest.TestCase):
    """Tests binary streaming hashing, format detection, and fast read counting."""

    def setUp(self):
        self.test_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_fast_count_and_sha256_on_valid_fastq_gz(self):
        fq_path = self.test_dir / "sample_R1.fastq.gz"
        records = [
            b"@READ1\nACGTACGT\n+\nIIIIIIII\n",
            b"@READ2\nTGCATGCA\n+\nIIIIIIII\n",
            b"@READ3\nGGCCAAAT\n+\nIIIIIIII\n",
        ]
        with gzip.open(fq_path, "wb") as f:
            f.writelines(records)

        count = IntegrityEngine.fast_count_fastq_reads(fq_path)
        self.assertEqual(count, 3)

        sha = IntegrityEngine.compute_sha256(fq_path)
        self.assertEqual(len(sha), 64)

        fp = IntegrityEngine.generate_fingerprint(fq_path, mode="strict", is_input=True)
        self.assertEqual(fp.read_count, 3)
        self.assertEqual(fp.sha256_checksum, sha)

    def test_truncated_fastq_raises_integrity_error(self):
        bad_fq = self.test_dir / "corrupted.fastq.gz"
        with gzip.open(bad_fq, "wb") as f:
            f.write(b"@READ1\nACGT\n+\n")  # Solo 3 righe invece di 4

        with self.assertRaises(IntegrityError):
            IntegrityEngine.fast_count_fastq_reads(bad_fq)


class TestFastqValidator(unittest.TestCase):
    """Tests structural syntax and paired-end header synchronization checks."""

    def setUp(self):
        self.test_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_valid_paired_end_synchronization(self):
        r1 = self.test_dir / "s1_R1.fastq.gz"
        r2 = self.test_dir / "s1_R2.fastq.gz"

        with gzip.open(r1, "wb") as f1, gzip.open(r2, "wb") as f2:
            f1.write(b"@M001:1:1:1 1:N:0:1\nACGT\n+\nIIII\n@M001:1:1:2 1:N:0:1\nTGCA\n+\nIIII\n")
            f2.write(b"@M001:1:1:1 2:N:0:1\nACGT\n+\nIIII\n@M001:1:1:2 2:N:0:1\nTGCA\n+\nIIII\n")

        self.assertTrue(FastqValidator.validate_paired_end_sync(r1, r2))

    def test_mismatched_sequence_quality_length_raises(self):
        bad_fq = self.test_dir / "length_mismatch.fastq.gz"
        with gzip.open(bad_fq, "wb") as f:
            f.write(b"@READ1\nACGTACGT\n+\nIIII\n")

        with self.assertRaises(FastqValidationError):
            FastqValidator.validate_fastq_structure(bad_fq)

    def test_desynchronized_paired_headers_raises(self):
        r1 = self.test_dir / "s1_R1.fastq.gz"
        r2 = self.test_dir / "s1_R2.fastq.gz"

        with gzip.open(r1, "wb") as f1, gzip.open(r2, "wb") as f2:
            f1.write(b"@SAMPLE_A:1:1\nACGT\n+\nIIII\n")
            f2.write(b"@SAMPLE_B:1:1\nACGT\n+\nIIII\n")

        with self.assertRaises(FastqValidationError):
            FastqValidator.validate_paired_end_sync(r1, r2)


class TestDiscoveryEngine(unittest.TestCase):
    """Tests right-to-left tokenized parsing, multi-lane aggregation, and pairing."""

    def setUp(self):
        self.test_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_illumina_token_extraction(self):
        path = Path("/data/MAR1A_S2_L002_R1_001.fastq.gz")
        disc = FilenameTokenizer.parse(path)
        self.assertEqual(disc.sample_id, "MAR1A")
        self.assertEqual(disc.library_id, "S2")
        self.assertEqual(disc.lane, "L002")
        self.assertEqual(disc.read_type, "R1")

    def test_multi_lane_aggregation(self):
        """Validates that multi-lane files for the same sample produce a single aggregated SampleRecord."""
        files = [
            self.test_dir / "SAMPLE_S1_L001_R1_001.fastq.gz",
            self.test_dir / "SAMPLE_S1_L001_R2_001.fastq.gz",
            self.test_dir / "SAMPLE_S1_L002_R1_001.fastq.gz",
            self.test_dir / "SAMPLE_S1_L002_R2_001.fastq.gz",
        ]
        for f in files:
            with gzip.open(f, "wb") as stream:
                stream.write(b"@READ1\nACGT\n+\nIIII\n")

        cfg = MockPipelineConfig(results_dir=self.test_dir)
        records = SampleDiscoveryEngine.build_sample_records(
            input_dir=self.test_dir,
            config=cfg,
            integrity_mode="fast",
        )

        self.assertEqual(len(records), 1)
        sample = records[0]
        self.assertEqual(sample.sample_id, "SAMPLE")
        self.assertEqual(len(sample.fingerprints), 4)
        self.assertEqual(sample.status, ValidationStatus.PASS)


class TestBAMValidator(unittest.TestCase):
    """Tests BAM integrity, BGZF truncation on real BAMs, and contig verification."""

    def setUp(self):
        self.test_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def _create_synthetic_bam(self, bam_path: Path) -> None:
        """Helper generating a physically valid, indexed coordinate-sorted BAM."""
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"LN": 16569, "SN": "NC_012920.1"}, {"LN": 248956422, "SN": "chr1"}],
        }
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out_bam:
            read = pysam.AlignedSegment()
            read.query_name = "test_read_1"
            read.query_sequence = "ACGT" * 10
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 100
            read.mapping_quality = 30
            read.cigar = ((0, 40),)
            read.query_qualities = pysam.qualitystring_to_array("I" * 40)
            out_bam.write(read)

        pysam.index(str(bam_path))

    def test_real_bam_eof_and_deliberate_truncation(self):
        """Generates a valid BAM, verifies EOF pass, strips the 28-byte BGZF block, and asserts fail."""
        bam_path = self.test_dir / "valid.bam"
        self._create_synthetic_bam(bam_path)

        # 1. BAM integro -> PASS
        self.assertTrue(BAMValidationEngine.check_bgzf_eof(bam_path))

        # 2. Troncamento deliberato degli ultimi 64 byte (rimuove il blocco finale BGZF)
        with open(bam_path, "r+b") as f:
            f.seek(0, 2)
            size = f.tell()
            f.truncate(size - 64)

        # 3. BAM tronco -> FAIL
        self.assertFalse(BAMValidationEngine.check_bgzf_eof(bam_path))

    def test_full_bam_inspection_workflow(self):
        """Tests end-to-end BAMValidationEngine diagnostic inspection."""
        bam_path = self.test_dir / "sample.bam"
        self._create_synthetic_bam(bam_path)

        sample = SampleRecord(
            identity=SampleIdentity(sample_id="BAM_Sample"),
            inputs=InputFiles(data_type=DataType.BAM, is_paired=False, bam_path=bam_path),
        )

        res = BAMValidationEngine.validate_and_inspect_bam(
            sample=sample,
            expected_contig="NC_012920.1",
            auto_detect_contig=True,
            auto_index=True,
        )

        self.assertTrue(res.is_valid)
        self.assertEqual(res.status, ValidationStatus.PASS)
        self.assertEqual(res.detected_mt_contig, "NC_012920.1")
        self.assertEqual(res.mt_mapped_reads, 1)


class TestTrimmerSettingsParser(unittest.TestCase):
    """Unit test for AdapterRemoval .settings metrics parsing."""

    def setUp(self):
        self.test_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_settings_parser_metrics_separation(self):
        settings_path = self.test_dir / "sample.settings"
        content = """
        Total number of read pairs: 1000
        Number of full-length collapsed pairs: 600
        Number of truncated collapsed pairs: 100
        Number of retained mate-pairs: 200
        Number of singletons: 50
        Number of discarded reads: 50
        """
        settings_path.write_text(content, encoding="utf-8")

        warnings = []
        metrics = PreprocessingEngine.parse_adapterremoval_settings(
            settings_file=settings_path,
            is_paired=True,
            warnings_list=warnings,
        )

        self.assertEqual(metrics.raw_read_pairs, 1000)
        self.assertEqual(metrics.raw_total_reads, 2000)
        self.assertEqual(metrics.merged_fragments, 700)
        self.assertEqual(metrics.surviving_read_pairs, 200)
        self.assertEqual(metrics.singleton_reads, 50)
        self.assertEqual(metrics.surviving_sequenced_fragments, 950)
        self.assertAlmostEqual(metrics.merged_fraction_raw, 0.70)
        self.assertAlmostEqual(metrics.merged_fraction_surviving, 700 / 900, places=4)


class TestBlock2IntegrationEndToEnd(unittest.TestCase):
    """End-to-end integration tests executing actual discovery, integrity validation, and AdapterRemoval binary."""

    def setUp(self):
        self.test_dir = Path(tempfile.mkdtemp())
        self.input_dir = self.test_dir / "raw_input"
        self.results_dir = self.test_dir / "results"
        self.input_dir.mkdir()
        self.results_dir.mkdir()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    @unittest.skipUnless(
        shutil.which("AdapterRemoval") or shutil.which("adapterremoval"),
        "AdapterRemoval binary not installed in PATH",
    )
    def test_block2_paired_end_end_to_end(self):
        """Full PE integration: Raw FASTQs -> Discovery -> AdapterRemoval -> ProcessedReads -> Logs -> QC."""
        r1_path = self.input_dir / "AncientPE_S1_L001_R1_001.fastq.gz"
        r2_path = self.input_dir / "AncientPE_S1_L001_R2_001.fastq.gz"

        # Frammento biologico e il suo complementare inverso per R2
        insert_seq = b"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        insert_rc = insert_seq.translate(bytes.maketrans(b"ACGT", b"TGCA"))[::-1]

        illumina_adapter1 = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
        illumina_adapter2 = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

        with gzip.open(r1_path, "wb") as f1, gzip.open(r2_path, "wb") as f2:
            for i in range(25):
                r1_full = insert_seq + illumina_adapter1
                r2_full = insert_rc + illumina_adapter2
                f1.write(f"@READ_{i} 1:N:0:1\n".encode() + r1_full + b"\n+\n" + b"I" * len(r1_full) + b"\n")
                f2.write(f"@READ_{i} 2:N:0:1\n".encode() + r2_full + b"\n+\n" + b"I" * len(r2_full) + b"\n")

        # 1. Configurazione Runtime
        cfg = MockPipelineConfig(
            results_dir=self.results_dir,
            min_read_length=20,
            min_merge_overlap=11,
            collapse_reads=True,
        )

        # 2. Discovery & Ingestion
        records = SampleDiscoveryEngine.build_sample_records(
            input_dir=self.input_dir,
            config=cfg,
            integrity_mode="strict",
        )
        self.assertEqual(len(records), 1)
        sample = records[0]
        self.assertEqual(sample.sample_id, "AncientPE")
        self.assertEqual(len(sample.fingerprints), 2)

        # 3. Preprocessing Execution
        processed_sample = PreprocessingEngine.process_sample(sample=sample, config=cfg)

        # 4. Verifica Strutturale ProcessedReads
        self.assertEqual(processed_sample.status, ValidationStatus.PASS)
        self.assertIsNotNone(processed_sample.processed_reads.merged_fastq)
        self.assertTrue(processed_sample.processed_reads.merged_fastq.is_file())
        self.assertGreater(processed_sample.processed_reads.merged_fastq.stat().st_size, 0)
        self.assertIn(
            processed_sample.processed_reads.output_profile,
            [OutputProfile.MERGED_ONLY, OutputProfile.MIXED],
        )

        # 5. Verifica Log su File Fisici
        log_dir = self.results_dir / "logs" / "preprocessing"
        stdout_log = log_dir / f"{sample.sample_id}.adapterremoval.log"
        stderr_log = log_dir / f"{sample.sample_id}.adapterremoval.err"
        self.assertTrue(stdout_log.is_file())
        self.assertTrue(stderr_log.is_file())

        # 6. Verifica Metriche Preprocessing
        self.assertGreaterEqual(processed_sample.metrics.merged_fragments, 1)
        self.assertGreater(processed_sample.metrics.surviving_sequenced_fragments, 0)

    @unittest.skipUnless(
        shutil.which("AdapterRemoval") or shutil.which("adapterremoval"),
        "AdapterRemoval binary not installed in PATH",
    )
    def test_block2_single_end_end_to_end(self):
        """Full SE integration: Single-End FASTQ -> Discovery -> AdapterRemoval -> ProcessedReads."""
        se_path = self.input_dir / "AncientSE.fastq.gz"
        read_seq = b"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" + b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"

        with gzip.open(se_path, "wb") as f:
            for i in range(20):
                f.write(f"@READ_{i}\n".encode() + read_seq + b"\n+\n" + b"I" * len(read_seq) + b"\n")

        cfg = MockPipelineConfig(
            results_dir=self.results_dir,
            min_read_length=20,
            collapse_reads=False,
            single_end=True,
            raw_data={"input": {"single_end": True}},
        )

        records = SampleDiscoveryEngine.build_sample_records(
            input_dir=self.input_dir,
            config=cfg,
            integrity_mode="fast",
        )
        self.assertEqual(len(records), 1)
        sample = records[0]
        self.assertFalse(sample.inputs.is_paired)

        processed_sample = PreprocessingEngine.process_sample(sample=sample, config=cfg)
        self.assertEqual(processed_sample.status, ValidationStatus.PASS)
        self.assertIsNotNone(processed_sample.processed_reads.singleton_fastq)
        self.assertTrue(processed_sample.processed_reads.singleton_fastq.is_file())


if __name__ == "__main__":
    unittest.main(verbosity=2)
