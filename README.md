

# 🧬 mtDNA Forensic & Ancient DNA Analysis Pipeline

[![Pipeline Version](https://img.shields.io/badge/version-1.2.0-blue.svg)](config/default_config.yaml)
[![Python](https://img.shields.io/badge/python-3.11%20%7C%203.12-brightgreen.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20WSL2%20%7C%20Windows-lightgrey.svg)]()

An end-to-end, forensic-grade bioinformatics pipeline designed for ancient (aDNA) and modern mitochondrial DNA (mtDNA) authentication, contamination quantification, high-resolution consensus calling, and phylogenetic haplogroup assignment.

---

## 📑 Table of Contents
- [Architecture Overview](#-architecture-overview)
- [Key Features](#-key-features)
- [Prerequisites & System Requirements](#-prerequisites--system-requirements)
- [Installation & Environment Setup](#-installation--environment-setup)
- [Reference Genomes & Databases](#-reference-genomes--databases)
- [Quickstart: Windows (WSL2 + Graphical HUD)](#-quickstart-windows-wsl2--graphical-hud)
- [Quickstart: Linux & HPC (CLI & Native GUI)](#-quickstart-linux--hpc-cli--native-gui)
- [Example Public Datasets (Where to Download)](#-example-public-datasets-where-to-download)
- [In-Depth Workflow Mechanics (Blocks 1–8)](#-in-depth-workflow-mechanics-blocks-18)
- [Configuration Guide (`default_config.yaml`)](#-configuration-guide-default_configyaml)
- [Production Deliverables & Forensic Reporting](#-production-deliverables--forensic-reporting)
- [Troubleshooting & FAQ](#-troubleshooting--faq)
- [Citation & References](#-citation--references)

---

## 🏛 Architecture Overview

The pipeline employs a strictly compartmentalized **8-Block Modular Architecture** where analytical stages communicate through tamper-evident, cryptographically sealed JSON contracts:

```text
=============================================================================
              🖥️  INTERFACE LAYER (ENTRYPOINTS)
=============================================================================
  [ GUI: Streamlit HUD ]                   [ CLI: Terminal Interface ]
  (src/gui/app.py, process_manager.py)     (src/pipeline/cli.py, job_runner.py)
               \                                 /
                \                               /
=============================================================================
              ⚙️  CORE ORCHESTRATOR & INFRASTRUCTURE
=============================================================================
 [ EndToEndPipelineRunner (runner.py) ] <--> [ StageTracker & ConfigLoader ]
       (Workflow Orchestration)                 (Atomic Caching & Invariant Audit)
                      |
                      | 🔒 Cryptographic JSON Contracts (SHA-256 Provenance)
                      v
=============================================================================
              🧬  BIOINFORMATICS ENGINE (BLOCKS 1 -> 7)
=============================================================================
                      |
 [ BLOCK 1: Reference Manager ]     --> Auto-indexing, rCRS NC_012920.1 validation
                      |
                      v (FASTQ input)
 [ BLOCK 2: Preprocessing ]         --> Adapter trimming & read merging (AdapterRemoval)
                      |
                      v (Direct BAM mode bypasses B2)
 [ BLOCK 3: Alignment & Dedup ]     --> Seeding-free BWA-aln (-l 1000 -n 0.01), 2-coord dedup
                      |
                      v
 [ BLOCK 4: Damage Authentication ] --> MapDamage2 / PMDtools, Bayesian quality score rescaling
                      |
                      v
 [ BLOCK 5: Contamination Gating ]  --> Two-tier gating: FastScreen Triage & Schmutzi MCMC
                      |                       |
                      |                       |-- [BIOLOGICAL GATING: Automated Safe Fallback]
                      v
 [ BLOCK 6: Consensus Assembly ]    --> Canonical, IUPAC, & Subconsensus VCF/FASTA + Spatial QC
                      |
                      v
 [ BLOCK 7: Haplogrouping ]         --> HaploGrep3 (Phylotree 17), Dual-Stream VCF/FASTA Audit
                      |
=============================================================================
              ⚖️  FORENSIC SYNTHESIS & AUDIT (BLOCK 8)
=============================================================================
                      |
 [ BLOCK 8: Synthesis Orchestrator ]
                      |
  +-------------------+-------------------+-------------------+
  |                   |                   |                   |
[ 🔐 4-Tier Provenance ] [ 📚 Mutation Catalog ] [ 📦 Packaging Engine ] [ 📊 Cohort Aggregator ]
 (Contract Validation)     (Tidy Variant Matrix)   (Atomic TSV/Reports)    (Master JSON Manifest)

```

---

## ⚡ Key Features

* **Cross-Platform Dual Interface**: Run headless on high-performance Linux compute nodes or launch via Windows desktop batch files with an interactive web-based Streamlit HUD.
* **Detached Supervisor Daemon**: GUI processes run detached under kernel-level advisory file locks (`fcntl.flock`) and track process start-time ticks to prevent PID-recycling hazards.
* **Damage-Aware Alignment**: Non-destructive alignment via `bwa-aln` tuned for ultrashort degraded DNA fragments (disabled seeding `-l 1000`, strict mismatch tolerance `-n 0.01`).
* **Two-Tier Contamination Defense**:
* *Tier 1 (FastScreen)*: Millisecond-level allelic triage across diagnostic Phylotree sites (`FAST_CLEAR` bypass).
* *Tier 2 (Bayesian MCMC)*: Full *Schmutzi* mixture modeling against a 197-genome contaminant database with automated subsampling (~40x) and *DissentEngine* fallback for modern/zero-damage specimens.


* **Biological Gating Actions**: Samples with unresolvable or critical contamination are prevented from emitting false-positive consensus sequences (`RESTRICT_HETEROPLASMY`, `PMD_FILTERING_REQUIRED`, or `EXCLUDE_OR_MANUAL_REVIEW`).
* **Spatial Damage Filtering**: Removes spurious deamination artifacts by filtering terminal variant loci against library-specific distance thresholds (Non-UDG, UDG-half, UDG-full).
* **Forensic 4-Tier Provenance Audit**: Verifies cryptographic hash continuity (SHA-256) across all contract boundaries before emitting deliverables.

---

## 💻 Prerequisites & System Requirements

* **Operating System**:
* **Linux**: Ubuntu 20.04+, Debian 11+, Fedora, Arch/Manjaro, RHEL 8+
* **Windows**: Windows 10/11 64-bit with **WSL2** (Windows Subsystem for Linux with Ubuntu 22.04 LTS recommended)


* **Conda Package Manager (Required)**:
* You **must** have **Miniconda** or **Anaconda** installed inside your Linux distribution (or inside WSL2 on Windows) to manage bioinformatics binaries (`bwa`, `samtools`, `AdapterRemoval`, `mapDamage2`, `haplogrep3`) and Python packages.


* **Hardware**:
* **CPU**: Minimum 4 physical cores (8+ recommended for multi-sample Schmutzi parallelization)
* **RAM**: Minimum 8 GB (16 GB+ recommended for large BAM deduplication and MCMC sampling)
* **Disk Space**: ~25 GB (including reference genomes, Schmutzi frequency databases, and Conda packages)



---

## 📦 Installation & Environment Setup

### Step 0: Install & Initialize Miniconda (If Conda is not yet installed)

If you do not have Conda installed on your Linux system or inside your WSL2 Ubuntu terminal, download and install **Miniconda**:

```bash
wget "(https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh)" -O /tmp/miniconda.sh
bash /tmp/miniconda.sh -b -p $HOME/miniconda
$HOME/miniconda/bin/conda init bash
$HOME/miniconda/bin/conda init zsh

```

> ⚠️ **CRITICAL**: After running `conda init`, you **MUST close and restart your terminal window** (or reboot your WSL/system session) so the shell loads Conda properly. If you already had a system-wide Conda installed, run `conda init bash` (or `conda init zsh` if using Zsh) and restart your terminal before proceeding.

### Step 1: Clone the Repository

Clone the repository directly into your Linux home directory (or inside WSL2):

```bash
cd ~
git clone "https://github.com/oscar-pizzoli/amtdna-variant-pipeline.git"
cd amtdna-variant-pipeline

```

> **Important**: Keep the directory path as `~/amtdna-variant-pipeline` if you intend to use the pre-configured Windows startup launchers.

### Step 2: Create and Activate the Conda Environment

Create the `adna_qc` environment directly from the locked `environment.yml` specification file (this automatically installs Python 3.12, `pip`, R packages, all Bioconda tools, and Python dependencies):

```bash
conda env create -f environment.yml
conda activate adna_qc

```

*(Alternative manual setup if you do not use `environment.yml`:)*

```bash
conda create -n adna_qc -c conda-forge -c bioconda python=3.12 pip bwa samtools adapterremoval=2.3.3 mapdamage2 pmdtools haplogrep3 -y
conda activate adna_qc
pip install -r requirements.txt

```

---

## 🧬 Reference Genomes & Databases

The pipeline requires the Revised Cambridge Reference Sequence (rCRS) and the Schmutzi 197 contaminant allele frequency database. Make sure your `adna_qc` environment is active (`conda activate adna_qc`) before running these commands.

### 1. rCRS Mitochondrial Reference

Download the canonical RefSeq rCRS FASTA (`NC_012920.1`) into `data/reference/rCRS.fasta` and build its indexes:

```bash
mkdir -p data/reference
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=fasta" > data/reference/rCRS.fasta
samtools faidx data/reference/rCRS.fasta
bwa index data/reference/rCRS.fasta

```

### 2. Schmutzi 197 Contaminant Allele Database

Extract the official 197 human mitogenome allele frequency database (`share/schmutzi/alleleFreqMT/197`) from the official Schmutzi repository into `data/reference/schmutzi_db/197`:

```bash
mkdir -p data/reference/schmutzi_db
git clone --depth 1 "https://github.com/grenaud/schmutzi.git" /tmp/schmutzi_repo
cp -r /tmp/schmutzi_repo/share/schmutzi/alleleFreqMT/197 data/reference/schmutzi_db/
rm -rf /tmp/schmutzi_repo

```

### 3. Diagnostic Screening Panel

Generate or verify the presence of `data/reference/diagnostic_sites_phylotree17.tsv`:

```bash
python -m src.modules.contamination.panel_builder \
    --freqs-dir data/reference/schmutzi_db/197/freqs \
    --reference-fasta data/reference/rCRS.fasta \
    --output-tsv data/reference/diagnostic_sites_phylotree17.tsv

```

---

## 🪟 Quickstart: Windows (WSL2 + Graphical HUD)

Windows users can operate the entire pipeline without typing terminal commands once the WSL2 `adna_qc` environment is installed.

### Running via Batch Launchers

Navigate to the `Windows_Startup/` directory in Windows File Explorer:

1. **Double-click `launch_hud_gui.bat**`:
* Spawns the background WSL2 service within the `adna_qc` environment.
* Automatically opens your default web browser to `http://localhost:8501`.


2. **Configure Pipeline Settings in the Sidebar**:
* **Sequence Input Mode**: Toggle between `BAM` or `FASTQ`.
* **Input Directory**: Enter standard Windows paths (e.g., `C:\ancient_dna\samples`) or WSL POSIX paths. Paths are automatically translated.
* **Cache Handling Policy**: Choose between `Ask on Start (Prompt)`, `Reuse Existing Cache`, or `Overwrite (Purge & Recompute)`.
* **Theme Palette**: Select your preferred color scheme.


3. **Click `▶ Start Pipeline**`:
* If previous computational cache is detected, an interactive dialog will appear asking whether to **Reuse**, **Overwrite**, or **Cancel**.
* Monitor live color-coded execution logs in the real-time log viewer.
* Once complete, browse the Cohort QC Dashboard, inspect individual specimens, and download publication-ready TSVs and Markdown reports.



To run purely in the Windows Command Prompt without opening a browser, double-click `Windows_Startup/run_pipeline_cli.bat`.

---

## 🐧 Quickstart: Linux & HPC (CLI & Native GUI)

### A. Terminal Execution (Production & Slurm Batch Jobs)

1. **Standard execution** using settings defined in `config/default_config.yaml`:

```bash
conda activate adna_qc
python scripts/run_complete_pipeline.py

```

2. **Custom FASTQ execution** (overriding sequence mode, input folder, and CPU budget):

```bash
python scripts/run_complete_pipeline.py \
    --data-type FASTQ \
    --input-dir data/raw_fastq \
    --cpu-budget-fraction 0.75 \
    --log-level INFO

```

3. **Direct BAM input mode** (bypasses Block 2 Adapter Trimming):

```bash
python scripts/run_complete_pipeline.py \
    --data-type BAM \
    --input-dir data/raw_bam

```

4. **Force clean recalculation** (purges existing cache and recomputes from scratch):

```bash
python scripts/run_complete_pipeline.py --cache-mode overwrite

```

### B. Launching the Graphical HUD on Linux

To run the Streamlit interface locally or over SSH tunneling:

```bash
conda activate adna_qc
streamlit run src/gui/app.py --server.port 8501

```

---

## 🔬 Example Public Datasets (Where to Download)

To test the pipeline on authentic archaeological specimens with established authenticity profiles, download public datasets from the European Nucleotide Archive (ENA) or NCBI Sequence Read Archive (SRA):

### 1. Ötzi the Tyrolean Iceman (~3300 BCE, Haplogroup K1f)

* **Accession**: ENA Project PRJEB37984 / Run ERR4058200 (or SRA SRR10560731)
* **Profile**: Shotgun WGS library, low endogenous mtDNA fraction (~0.11%), typical terminal deamination (~7.5% C>T), narrow Schmutzi contamination (~1.0%).
* **Download command**:

```bash
mkdir -p data/raw_fastq
fastq-dump --split-files --gzip SRR10560731 -O data/raw_fastq/

```

### 2. Ancient Egyptian Mummies (Abusir el-Meleq, Haplogroups W6, U3b)

* **Accession**: ENA Study PRJEB12470 (Schuenemann et al., Nature Communications 2017)
* **Profile**: Hybridization capture enrichment, high coverage (>150x), high endogenous content (>99%), severe terminal damage (8–15% C>T), low fragment lengths (50–80 bp).
* **Download command (Sample JK2888 / Run ERR1830678)**:

```bash
mkdir -p data/raw_fastq
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR183/008/ERR1830678/ERR1830678.fastq.gz" -P data/raw_fastq/

```

### 3. Modern Human Control (Clean Negative Control)

* **Accession**: 1000 Genomes Project (e.g., sample NA12878)
* **Profile**: Zero terminal deamination, modern fragmentation profile (>150 bp), used to test the DissentEngine bypass path.

---

## 🔍 In-Depth Workflow Mechanics (Blocks 1–8)

### Block 1: Reference Genome Validation (`ReferenceManager`)

* Loads and cryptographically checks the rCRS mitochondrial genome (`NC_012920.1`, 16,569 bp).
* Verifies sequence checksums (MD5/SHA-256), builds missing BWA and `.fai` indexes automatically, and checks for circular coordinate wraparound support.

### Block 2: Read Discovery & Preprocessing (`PreprocessingEngine`)

* Validates FASTQ read integrity (detecting truncated blocks or corrupted gzip streams) and infers Phred score encoding (Phred+33 vs Phred+64).
* Executes `AdapterRemoval`: strips sequencing adapters, removes terminal ambiguous `N`s, filters reads under 30 bp, and merges overlapping paired-end reads into reconstructed full-length singleton fragments.

### Block 3: Alignment, Quality Filtering & Deduplication (`AlignmentOrchestrator`)

* Aligns reads to the rCRS using `bwa-aln` with seeding disabled (`-l 1000 -n 0.01 -o 2`) to ensure that post-mortem terminal mismatches do not penalize alignment scores.
* Executes two-coordinate duplicate removal (`mode: both` via `samtools markdup`), retaining independent physical DNA molecules.
* Applies strict filters: removes unmapped reads, secondary alignments (`0x100`), supplementary alignments (`0x800`), and low-quality mapping (`MAPQ < 30`).

### Block 4: aDNA Damage Profiling & Quality Rescaling (`DamageOrchestrator`)

* Quantifies cytosine deamination rates (C to T at 5' ends, G to A at 3' ends) and fragment length distributions via `mapDamage2`.
* Fits statistical models to compute Bayesian damage downscaling for base qualities: bases susceptible to post-mortem deamination have their Phred scores adjusted downwards to prevent false-positive mutation calls.
* Emits individual read degradation metrics via `PMDtools`.
* Compiles the sealed `block4_contract.json` declaring the authenticity classification (`AUTHENTIC_ANCIENT_SUPPORTED`, `INCONCLUSIVE`, or `MODERN`).

### Block 5: Two-Tier Contamination Estimation & Gating (`ContaminationBatchRunner`)

* **Tier 1 (Fast Pre-Screening)**: Evaluates allelic dissent across 354 diagnostic Phylotree loci. Samples showing zero discordant signals receive a `FAST_CLEAR` verdict, bypassing expensive Bayesian MCMC iterations.
* **Tier 2 (Modeling Pool)**:
* Samples certified as ancient enter **Schmutzi**: the sample depth is downsampled conditionally to 40x (preventing MCMC matrix overflow), MD tags are injected via `samtools calmd`, and iterative Bayesian deamination modeling is executed against 197 modern contaminant genomes.
* Samples certified as modern/zero-damage are routed directly to `DissentEngine` to avoid zero-damage numeric singularities (`-nan`).


* **Formulates formal gating constraints for Block 6**:
* **<= 2%**: `STANDARD_CALLING` (Unrestricted calling).
* **2% - 5%**: `RESTRICT_HETEROPLASMY` (Suppresses low-frequency variants; raises VAF floor to 25%).
* **5% - 10%**: `PMD_FILTERING_REQUIRED` / `USE_ENDOGENOUS_SCHMUTZI_CONSENSUS`.
* **> 10%**: `EXCLUDE_OR_MANUAL_REVIEW`.



### Block 6: Dual-Stream Consensus & Heteroplasmy Calling (`ConsensusOrchestrator`)

* **Generates**:
* `sample.consensus.fasta`: Canonical consensus sequence (fixed major alleles >= 80%).
* `sample.consensus_iupac.fasta`: IUPAC ambiguity-coded consensus capturing authentic sub-clonal variations.
* `sample.variants.vcf.gz`: Complete VCF documenting depth, VAF, strand bias, and statistical filters.


* **Spatial Damage Filtration**: Suppresses candidate heteroplasmies located within the terminal ends of reads (e.g., within 5 bp for Non-UDG libraries) to prevent residual deamination errors from masquerading as biological heteroplasmies.
* Blacklists known homopolymer artifact regions (HVS-I/II poly-C tracts and the artificial 3107N site).

### Block 7: Phylogenetic Classification & Lineage Audit (`HaplogroupOrchestrator`)

* Assigns mitochondrial haplogroups using `HaploGrep3` based on the Phylotree 17 database (`phylotree-fu-rcrs@1.2`) using the Kulczynski distance metric.
* Conducts dual-stream concordance cross-checks (classifying FASTA and VCF streams simultaneously).
* Audits private and missing diagnostic mutations against contaminant haplogroups identified by Schmutzi in Block 5 to detect subtle chimeric mixtures.

### Block 8: Forensic Synthesis & Reporting (`SynthesisOrchestrator`)

* **Evaluates the 4-Tier Provenance Audit**:
* *Tier 1 (Contracts)*: Cryptographic completeness of all Block 2–7 contracts.
* *Tier 2 (Artifacts)*: Physical presence, size, and EOF block integrity of BAM, VCF, and FASTA files.
* *Tier 3 (Semantics)*: Cross-stage numerical sanity checks (e.g., consensus length matching reference length).
* *Tier 4 (Lineage)*: SHA-256 digest chain verification anchored to `execution_manifest.json`.


* Emits production-grade deliverables: Markdown forensic case reports, tidy mutation matrices, cohort summaries, and the finalized Root of Trust manifest.

---

## ⚙️ Configuration Guide (`default_config.yaml`)

All parameters are configured in `config/default_config.yaml`. Key sections include:

```yaml
resources:
  cpu_budget_fraction: 0.50
  threads: 8

output:
  results_dir: "results"
  operational_files_dir: "operational_files"
  cache_mode: "prompt"

input:
  data_type: "FASTQ"
  raw_fastq_dir: "data/raw_fastq"
  raw_bam_dir: "data/raw_bam"

alignment:
  bwa_aln:
    seed_length: 1000
    max_diff: 0.01

deduplication:
  mode: "both"

damage:
  library_treatment: "NONE"
  execution:
    rescaling:
      enabled: true

contamination:
  subsampling:
    enabled: true
    target_depth: 40.0
  thresholds:
    low_threshold: 0.02
    moderate_threshold: 0.05
    elevated_threshold: 0.10

```

---

## 📊 Production Deliverables & Forensic Reporting

Upon pipeline completion, final deliverables are written to `results/`:

```text
results/
├── execution_manifest.json          # Root of Trust: cryptographic SHA-256 lineage manifest
├── fasta/
│   ├── cohort_canonical.fasta       # Multi-FASTA of all passing canonical consensus sequences
│   └── cohort_heteroplasmy.fasta    # Multi-FASTA with IUPAC ambiguity codes
├── haplogroups/
│   └── cohort_haplogroups_summary.tsv # Complete haplogroup classification summary table
├── master_summary/
│   ├── cohort_master.json           # Aggregated JSON database of all specimens
│   ├── cohort_qc_summary.tsv        # Master forensic quality control summary table
│   └── WARNINGS_GLOSSARY.md         # Diagnostic warning codes glossary
├── mutations/
│   ├── cohort_mutation_catalog.tsv  # Tidy variant table across all cohort specimens
│   └── cohort_mutations_summary.tsv # Polymorphic site frequencies
└── reports_markdown/
    └── <SAMPLE_ID>/
        └── <SAMPLE_ID>.report.md    # Self-contained forensic case inspection report

```

---

## 🛠 Troubleshooting & FAQ

#### 1. `CondaError: Run 'conda init' before 'conda activate'`

Run `conda init bash` (or `conda init zsh` if your terminal uses Zsh), then **close and reopen your terminal window** before running `conda activate adna_qc`.

#### 2. Why was my ancient sample routed to DissentEngine instead of Schmutzi?

If Block 4 classifies a specimen as `MODERN` or `INCONCLUSIVE` (e.g., terminal deamination < 3%), the pipeline bypasses Schmutzi to prevent mathematical division-by-zero crashes (`-nan` errors in `contDeam.pl`). Check the `block4_contract.json` file in `operational_files/damage/<sample>/` to inspect observed damage rates.

#### 3. How do I force a complete recalculation when input data changes?

If cached results prevent stages from re-running, select **Overwrite (Purge & Recompute)** from the HUD sidebar or pass `--cache-mode overwrite` on the CLI:

```bash
python scripts/run_complete_pipeline.py --cache-mode overwrite

```

#### 4. What does the warning RESTRICT_HETEROPLASMY mean in my report?

When a sample exhibits moderate exogenous contamination (2% - 5%) or low estimation confidence, the pipeline raises the minimum Variant Allele Frequency (VAF) threshold to 25% and minimum depth to 20x. This prevents low-frequency modern contaminant reads from being misidentified as authentic somatic heteroplasmies.

#### 5. WSL Path Resolution Issues on Windows

Ensure that input directories entered into the HUD exist and are accessible from WSL2. Windows paths (e.g., `C:\data`) are automatically translated to `/mnt/c/data`. If using mapped network drives, mount them inside WSL2 first.

---

## 📜 Citation & References

If you use this pipeline in academic or forensic casework, please cite:

* **BWA**: Li H. & Durbin R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics, 25(14):1754-1760.
* **AdapterRemoval**: Schubert M. et al. (2016). AdapterRemoval v2: rapid adapter trimming, identification, and read merging. BMC Research Notes, 9:88.
* **mapDamage2**: Jónsson H. et al. (2013). mapDamage2.0: fast approximate Bayesian estimates of ancient DNA damage parameters. Bioinformatics, 29(13):1682-1684.
* **Schmutzi**: Renaud G. et al. (2015). schmutzi: estimation of contamination and endogenous mitochondrial consensus calling for ancient DNA. Genome Biology, 16:224.
* **HaploGrep3**: Schönherr S. et al. (2023). Haplogrep 3 - an interactive tool for mitochondrial DNA haplogroup classification. Nucleic Acids Research, 51(W1):W207–W212.

```

---
