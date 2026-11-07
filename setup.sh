#!/usr/bin/env bash
# ==============================================================================
# 🧬 Ancient & Modern mtDNA Variant Pipeline - Automated Environment Setup
# ==============================================================================
# This script automates the complete dependency installation and environment
# preparation for Debian/Ubuntu systems based on repository guidelines.
#
# Tasks performed:
#   1. Detect Conda installation and initialize shell environment
#   2. Check and install essential Debian system tools (curl, git, etc.)
#   3. Configure Bioconda and Conda-Forge channels
#   4. Create or update the 'adna_qc' Conda environment (environment.yml & requirements.txt)
#   5. Download and index the rCRS mitochondrial reference genome (NC_012920.1)
#   6. Clone and install the Schmutzi 197 contaminant allele frequency database
#   7. Build the diagnostic Phylotree 17 pre-screening panel
#   8. Initialize project directory structure and verify pre-flight health
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# Visual Styling & Output Helpers
# ------------------------------------------------------------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${CYAN}[INFO]${NC} $1"
}

log_step() {
    echo -e "\n${BOLD}${BLUE}==>${NC} ${BOLD}$1${NC}"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

# Resolve project root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

echo -e "${BOLD}${CYAN}"
echo "================================================================================"
echo "    🧬 mtDNA FORENSIC & ANCIENT DNA ANALYSIS PIPELINE - SETUP & INSTALL"
echo "================================================================================"
echo -e "${NC}"

# ------------------------------------------------------------------------------
# STEP 1: Conda Detection and Shell Initialization
# ------------------------------------------------------------------------------
log_step "Step 1/8: Detecting Conda package manager..."

CONDA_HOOK_PATH=""
CONDA_BIN=""

# 1. Check if conda binary is in PATH
if command -v conda &>/dev/null; then
    CONDA_BIN="$(command -v conda)"
    CONDA_BASE_DIR="$(dirname "$(dirname "${CONDA_BIN}")")"
    if [[ -f "${CONDA_BASE_DIR}/etc/profile.d/conda.sh" ]]; then
        CONDA_HOOK_PATH="${CONDA_BASE_DIR}/etc/profile.d/conda.sh"
    fi
fi

# 2. Check standard installation directories if not found in PATH
if [[ -z "${CONDA_HOOK_PATH}" ]]; then
    POSSIBLE_LOCATIONS=(
        "${HOME}/miniconda3"
        "${HOME}/miniconda"
        "${HOME}/miniforge3"
        "${HOME}/anaconda3"
        "/opt/conda"
        "/usr/local/miniconda"
        "/usr/local/miniconda3"
        "/usr/local/anaconda3"
    )
    for loc in "${POSSIBLE_LOCATIONS[@]}"; do
        if [[ -f "${loc}/etc/profile.d/conda.sh" ]]; then
            CONDA_HOOK_PATH="${loc}/etc/profile.d/conda.sh"
            CONDA_BIN="${loc}/bin/conda"
            export PATH="${loc}/bin:${PATH}"
            break
        fi
    done
fi

if [[ -z "${CONDA_HOOK_PATH}" ]]; then
    log_error "Conda could not be found on this system."
    echo -e "${YELLOW}"
    echo "This pipeline requires Miniconda or Anaconda to install all bioinformatics binaries."
    echo "Please install Miniconda first by executing:"
    echo ""
    echo "  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh"
    echo "  bash /tmp/miniconda.sh -b -p \$HOME/miniconda"
    echo "  \$HOME/miniconda/bin/conda init bash"
    echo "  source \$HOME/miniconda/etc/profile.d/conda.sh"
    echo ""
    echo "Then re-run this setup script: bash setup.sh"
    echo -e "${NC}"
    exit 1
fi

log_info "Sourcing Conda hook from: ${CONDA_HOOK_PATH}"
# shellcheck disable=SC1090
source "${CONDA_HOOK_PATH}"
log_success "Conda detected: $("${CONDA_BIN}" --version)"

# ------------------------------------------------------------------------------
# STEP 2: Debian / Ubuntu System Dependencies
# ------------------------------------------------------------------------------
log_step "Step 2/8: Checking Debian system prerequisites..."

REQUIRED_SYS_TOOLS=("curl" "wget" "git" "bzip2" "unzip" "ca-certificates")
MISSING_SYS_TOOLS=()

for tool in "${REQUIRED_SYS_TOOLS[@]}"; do
    if ! command -v "${tool}" &>/dev/null; then
        MISSING_SYS_TOOLS+=("${tool}")
    fi
done

if [[ ${#MISSING_SYS_TOOLS[@]} -gt 0 ]]; then
    log_info "Missing system tools: ${MISSING_SYS_TOOLS[*]}"
    if [[ $(id -u) -eq 0 ]]; then
        log_info "Running apt-get update and installing missing tools..."
        apt-get update -y && apt-get install -y "${MISSING_SYS_TOOLS[@]}"
    elif command -v sudo &>/dev/null; then
        log_info "Attempting to install missing tools via sudo apt-get..."
        sudo apt-get update -y && sudo apt-get install -y "${MISSING_SYS_TOOLS[@]}" || {
            log_warn "Could not automatically install system tools with sudo. Please run: sudo apt install -y ${MISSING_SYS_TOOLS[*]}"
        }
    else
        log_warn "Neither root nor sudo available. Please ensure ${MISSING_SYS_TOOLS[*]} are installed via apt."
    fi
else
    log_success "All required Debian system tools are installed."
fi

# ------------------------------------------------------------------------------
# STEP 3: Configure Conda Channels
# ------------------------------------------------------------------------------
log_step "Step 3/8: Configuring Conda channel priorities..."

conda config --add channels defaults 2>/dev/null || true
conda config --add channels bioconda 2>/dev/null || true
conda config --add channels conda-forge 2>/dev/null || true
conda config --set channel_priority flexible 2>/dev/null || true
log_success "Channels configured (conda-forge, bioconda, defaults)."

# ------------------------------------------------------------------------------
# STEP 4: Conda Environment Setup ('adna_qc')
# ------------------------------------------------------------------------------
log_step "Step 4/8: Building 'adna_qc' environment from environment.yml..."

ENV_NAME="adna_qc"
ENV_EXISTS=false

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
    ENV_EXISTS=true
fi

if [[ "${ENV_EXISTS}" = true ]]; then
    log_info "Environment '${ENV_NAME}' already exists. Updating packages to match environment.yml..."
    conda env update -n "${ENV_NAME}" -f environment.yml --prune
    log_success "Environment '${ENV_NAME}' updated successfully."
else
    log_info "Creating new environment '${ENV_NAME}' from environment.yml (this may take a few minutes)..."
    conda env create -f environment.yml
    log_success "Environment '${ENV_NAME}' created successfully."
fi

log_info "Activating '${ENV_NAME}' environment..."
conda activate "${ENV_NAME}"

if [[ -f "requirements.txt" ]]; then
    log_info "Ensuring all pip dependencies are satisfied..."
    pip install -r requirements.txt --quiet
    log_success "Pip dependencies verified."
fi

# ------------------------------------------------------------------------------
# STEP 5: Reference Genome Download & Indexing (rCRS)
# ------------------------------------------------------------------------------
log_step "Step 5/8: Validating & preparing rCRS mitochondrial reference..."

REF_DIR="${SCRIPT_DIR}/data/reference"
REF_FASTA="${REF_DIR}/rCRS.fasta"
mkdir -p "${REF_DIR}"

if [[ ! -s "${REF_FASTA}" ]]; then
    log_info "Downloading canonical rCRS FASTA (NC_012920.1) from NCBI Entrez..."
    curl -fsSL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=fasta" -o "${REF_FASTA}"
    log_success "rCRS.fasta downloaded successfully."
else
    log_info "rCRS.fasta already exists."
fi

# Index with samtools faidx
if [[ ! -f "${REF_FASTA}.fai" ]]; then
    log_info "Generating SAMtools FAI index..."
    samtools faidx "${REF_FASTA}"
    log_success "SAMtools index created."
else
    log_info "SAMtools FAI index verified."
fi

# Index with BWA
if [[ ! -f "${REF_FASTA}.bwt" ]]; then
    log_info "Generating BWA index for rCRS..."
    bwa index "${REF_FASTA}"
    log_success "BWA index created."
else
    log_info "BWA index verified."
fi

# ------------------------------------------------------------------------------
# STEP 6: Schmutzi 197 Contaminant Allele Database
# ------------------------------------------------------------------------------
log_step "Step 6/8: Installing Schmutzi 197 contaminant allele frequency database..."

SCHMUTZI_DB_DIR="${REF_DIR}/schmutzi_db/197"
mkdir -p "${REF_DIR}/schmutzi_db"

if [[ ! -d "${SCHMUTZI_DB_DIR}/freqs" ]] || [[ -z "$(ls -A "${SCHMUTZI_DB_DIR}/freqs" 2>/dev/null)" ]]; then
    log_info "Cloning official Schmutzi repository to extract 197 database..."
    TMP_SCHMUTZI_REPO="$(mktemp -d -t schmutzi_repo_XXXXXX)"
    git clone --depth 1 "https://github.com/grenaud/schmutzi.git" "${TMP_SCHMUTZI_REPO}"
    
    if [[ -d "${TMP_SCHMUTZI_REPO}/share/schmutzi/alleleFreqMT/197" ]]; then
        cp -r "${TMP_SCHMUTZI_REPO}/share/schmutzi/alleleFreqMT/197" "${REF_DIR}/schmutzi_db/"
        log_success "Schmutzi 197 database copied to data/reference/schmutzi_db/197."
    else
        log_error "Could not find share/schmutzi/alleleFreqMT/197 in cloned Schmutzi repository."
        rm -rf "${TMP_SCHMUTZI_REPO}"
        exit 1
    fi
    rm -rf "${TMP_SCHMUTZI_REPO}"
else
    log_info "Schmutzi 197 database already present and verified."
fi

# ------------------------------------------------------------------------------
# STEP 7: Diagnostic Screening Panel Construction
# ------------------------------------------------------------------------------
log_step "Step 7/8: Generating Phylotree 17 diagnostic screening panel..."

PANEL_TSV="${REF_DIR}/diagnostic_sites_phylotree17.tsv"

if [[ ! -s "${PANEL_TSV}" ]]; then
    log_info "Building diagnostic screening panel via panel_builder module..."
    python -m src.modules.contamination.panel_builder \
        --freqs-dir "${REF_DIR}/schmutzi_db/197/freqs" \
        --reference-fasta "${REF_FASTA}" \
        --output-tsv "${PANEL_TSV}"
    log_success "Diagnostic panel generated: ${PANEL_TSV}"
else
    log_info "Diagnostic panel already present and verified: ${PANEL_TSV}"
fi

# ------------------------------------------------------------------------------
# STEP 8: Directory Skeleton & Health Verification
# ------------------------------------------------------------------------------
log_step "Step 8/8: Initializing directory skeleton and verifying pipeline..."

mkdir -p "${SCRIPT_DIR}/data/raw_fastq"
mkdir -p "${SCRIPT_DIR}/data/raw_bam"
mkdir -p "${SCRIPT_DIR}/operational_files"
mkdir -p "${SCRIPT_DIR}/results"
mkdir -p "${SCRIPT_DIR}/logs"

chmod +x "${SCRIPT_DIR}/scripts/"*.py 2>/dev/null || true
chmod +x "${SCRIPT_DIR}/setup.sh" 2>/dev/null || true

log_info "Executing pipeline CLI pre-flight check..."
if python scripts/run_complete_pipeline.py --help >/dev/null 2>&1; then
    log_success "Pipeline pre-flight check PASSED."
else
    log_error "Pipeline pre-flight execution check failed."
    exit 1
fi

echo -e "\n${BOLD}${GREEN}================================================================================${NC}"
echo -e "${BOLD}${GREEN}    🎉 SETUP COMPLETED SUCCESSFULLY!${NC}"
echo -e "${BOLD}${GREEN}================================================================================${NC}"
echo -e "Your forensic aDNA/mtDNA analysis environment is completely configured."
echo ""
echo -e "${BOLD}To start using the pipeline:${NC}"
echo ""
echo -e "  1. Activate the Conda environment:"
echo -e "     ${CYAN}conda activate ${ENV_NAME}${NC}"
echo ""
echo -e "  2. Run via Command Line Interface (CLI):"
echo -e "     ${CYAN}python scripts/run_complete_pipeline.py --data-type FASTQ --input-dir data/raw_fastq${NC}"
echo ""
echo -e "  3. Or launch the Graphical Web HUD (Streamlit):"
echo -e "     ${CYAN}streamlit run src/gui/app.py${NC}"
echo ""
echo -e "  (On Windows with WSL2, you can also double-click ${BOLD}Windows_Startup/launch_hud_gui.bat${NC})"
echo -e "${BOLD}${GREEN}================================================================================${NC}\n"
