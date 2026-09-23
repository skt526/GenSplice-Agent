#!/usr/bin/env bash
# ==============================================================================
# GenSplice-Agent Automated System Pre-flight Check & One-Click Installer
# ==============================================================================

set -e

# ANSI Color Definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[1;36m'
RESET='\033[0m'
DIM='\033[2m'
NC='\033[0m' # No Color

ENV_NAME="gensplice-agent"
FORCE_INSTALL=false

for arg in "$@"; do
    if [ "$arg" = "--force" ] || [ "$arg" = "-f" ]; then
        FORCE_INSTALL=true
    fi
done

echo -e "${BLUE}====================================================${NC}"
echo -e "${BLUE}     GenSplice-Agent: Automated One-Click Installer ${NC}"
echo -e "${BLUE}====================================================${NC}"

# 1. Hardware Specification Auto-Detection (Linux & macOS)
if command -v nproc &>/dev/null; then
    CPU_CORES=$(nproc)
else
    CPU_CORES=$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
fi

# RAM Calculation in GB
if [ -f /proc/meminfo ]; then
    RAM_TOTAL_KB=$(grep MemTotal /proc/meminfo | awk '{print $2}')
    RAM_TOTAL_GB=$(awk "BEGIN {printf \"%.0f\", $RAM_TOTAL_KB / 1024 / 1024}")
elif command -v sysctl &>/dev/null; then
    RAM_BYTES=$(sysctl -n hw.memsize 2>/dev/null || echo 34359738368)
    RAM_TOTAL_GB=$(( RAM_BYTES / 1024 / 1024 / 1024 ))
else
    RAM_TOTAL_GB=32
fi

# Disk Free Space Calculation in GB
DISK_FREE_KB=$(df -k . | awk 'NR==2 {print $4}')
DISK_FREE_GB=$(awk "BEGIN {printf \"%.0f\", $DISK_FREE_KB / 1024 / 1024}")

# 2. Display Measured Hardware Specs
echo -e "CPU Threads           : ${YELLOW}${CPU_CORES}${NC}"
echo -e "Total RAM Size        : ${YELLOW}${RAM_TOTAL_GB} GB${NC}"
echo -e "Available Disk Space  : ${YELLOW}${DISK_FREE_GB} GB${NC}"
echo -e "----------------------------------------------------"

# 3. Minimum Requirements Thresholds
MIN_RAM_GB=16
MIN_DISK_GB=20
MIN_CPU=2

IS_SUITABLE=true
FAIL_REASONS=()

if [ "$RAM_TOTAL_GB" -lt "$MIN_RAM_GB" ]; then
    IS_SUITABLE=false
    FAIL_REASONS+=("Low RAM: Current ${RAM_TOTAL_GB} GB (Recommended for STAR Alignment: 30GB+)")
fi

if [ "$DISK_FREE_GB" -lt "$MIN_DISK_GB" ]; then
    IS_SUITABLE=false
    FAIL_REASONS+=("Low Disk Space: Current ${DISK_FREE_GB} GB (Recommended: 50GB+)")
fi

# 4. Gatekeeper Evaluation
if [ "$IS_SUITABLE" = false ]; then
    echo -e "${YELLOW}[Notice: Resource Warning]${NC}"
    for reason in "${FAIL_REASONS[@]}"; do
        echo -e "  - ${YELLOW}${reason}${NC}"
    done
    echo -e "  ${GREEN}✔ Proceeding with setup...${NC}"
else
    echo -e "${GREEN}[SUITABLE: Proceeding with Installation]${NC}"
    echo -e "System specifications are suitable for transcriptomics analysis."
fi

# 5. Calculate Thread & RAM Allocation
if [ "$CPU_CORES" -gt 2 ]; then
    OPTIMAL_THREADS=$((CPU_CORES - 2))
else
    OPTIMAL_THREADS=1
fi

STAR_RAM_GB=$(( RAM_TOTAL_GB * 75 / 100 ))
if [ "$STAR_RAM_GB" -lt 4 ]; then
    STAR_RAM_GB=4
fi
STAR_RAM_BYTES=$(( STAR_RAM_GB * 1024 * 1024 * 1024 ))

echo -e "----------------------------------------------------"
echo -e "System Resource Optimization Settings:"
echo -e "  - Pipeline Threads           : ${GREEN}${OPTIMAL_THREADS}${NC} threads (out of ${CPU_CORES} total)"
echo -e "  - Dynamic STAR BAM Sort RAM  : ${GREEN}${STAR_RAM_GB} GB${NC} (${STAR_RAM_BYTES} bytes)"
echo -e "----------------------------------------------------"

# 6. Auto-generate / Update config.yaml
cat <<EOF > config.yaml
# GenSplice-Agent Auto-generated Configuration
system:
  total_cpu: ${CPU_CORES}
  assigned_threads: ${OPTIMAL_THREADS}
  total_ram_gb: ${RAM_TOTAL_GB}
  star_bam_sort_ram_gb: ${STAR_RAM_GB}
  star_bam_sort_ram_bytes: ${STAR_RAM_BYTES}

reference:
  organism: "human"
  ref_dir: "./human-ref"
  fasta: "./human-ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
  gtf: "./human-ref/Homo_sapiens.GRCh38.113.gtf"
  star_index: "./human-ref/star_index"

threads: ${OPTIMAL_THREADS}

inputs:
  control_dir: "./inputs/control"
  treatment_dir: "./inputs/treatment"

outputs:
  clean_fq: "./outputs/01_clean_fq"
  aligned_bam: "./outputs/02_aligned_bam"
  deg: "./outputs/03_deg"
  rmats: "./outputs/04_rmats"
EOF

echo -e "✔ Updated 'config.yaml' with assigned_threads=${OPTIMAL_THREADS} and STAR RAM=${STAR_RAM_GB}GB."

# 7. Create directory structure
mkdir -p inputs/control inputs/treatment outputs

# 8. Conda/Miniforge Auto-Detection & Auto-Installation
echo -e "\n${BLUE}[1/2] Detecting Conda / Miniforge Environment...${NC}"

CONDA_CMD=""

# A. Check existing command in PATH
if command -v mamba &>/dev/null; then
    CONDA_CMD="mamba"
elif command -v conda &>/dev/null; then
    CONDA_CMD="conda"
elif command -v micromamba &>/dev/null; then
    CONDA_CMD="micromamba"
fi

# B. Check standard installation directories if not in PATH
if [ -z "$CONDA_CMD" ]; then
    for candidate in "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/anaconda3" "/opt/conda"; do
        if [ -f "$candidate/etc/profile.d/conda.sh" ]; then
            echo -e "  Found Conda installation at ${candidate}. Sourcing profile..."
            source "$candidate/etc/profile.d/conda.sh"
            CONDA_CMD="conda"
            break
        elif [ -x "$candidate/bin/conda" ]; then
            export PATH="$candidate/bin:$PATH"
            CONDA_CMD="conda"
            break
        fi
    done
fi

# C. Automatic Download & Install Miniforge if Conda is completely absent
if [ -z "$CONDA_CMD" ]; then
    echo -e "${YELLOW}Notice: Conda/Miniforge is not detected on this system.${NC}"
    echo -e "${GREEN}⚡ Automatically installing Miniforge3 into \$HOME/miniforge3...${NC}"
    
    OS_TYPE=$(uname -s)
    ARCH_TYPE=$(uname -m)
    MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-${OS_TYPE}-${ARCH_TYPE}.sh"
    INSTALLER_SCRIPT="/tmp/miniforge_installer.sh"
    
    echo -e "  Downloading Miniforge3 installer from: ${MINIFORGE_URL}"
    if command -v curl &>/dev/null; then
        curl -L --progress-bar "$MINIFORGE_URL" -o "$INSTALLER_SCRIPT"
    elif command -v wget &>/dev/null; then
        wget -O "$INSTALLER_SCRIPT" "$MINIFORGE_URL"
    else
        echo -e "${RED}Error: Neither curl nor wget was found. Please install curl or wget first.${NC}"
        exit 1
    fi
    
    echo -e "  Running Miniforge3 batch installation..."
    bash "$INSTALLER_SCRIPT" -b -p "$HOME/miniforge3"
    rm -f "$INSTALLER_SCRIPT"
    
    if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
        source "$HOME/miniforge3/etc/profile.d/conda.sh"
        "$HOME/miniforge3/bin/conda" init bash 2>/dev/null || true
        CONDA_CMD="conda"
        echo -e "${GREEN}✔ Successfully installed Miniforge3 at \$HOME/miniforge3!${NC}"
    else
        echo -e "${RED}Error: Miniforge3 installation failed. Please check internet connection or install Conda manually.${NC}"
        exit 1
    fi
fi

# 9. Create or Update Conda Environment
echo -e "\n${BLUE}[2/2] Managing Conda Environment '${ENV_NAME}'...${NC}"

ENV_EXISTS=false
if $CONDA_CMD env list | grep -qE "^${ENV_NAME}\s"; then
    ENV_EXISTS=true
fi

# Check if environment binaries are present
ENV_BIN_DIR="$HOME/miniforge3/envs/${ENV_NAME}/bin"
if [ ! -d "$ENV_BIN_DIR" ]; then
    ENV_BIN_DIR="$HOME/miniconda3/envs/${ENV_NAME}/bin"
fi

IS_ENV_COMPLETE=false
if [ "$ENV_EXISTS" = true ] && [ -x "${ENV_BIN_DIR}/fastp" ] && [ -x "${ENV_BIN_DIR}/STAR" ] && [ -x "${ENV_BIN_DIR}/rmats.py" ]; then
    IS_ENV_COMPLETE=true
fi

if [ "$ENV_EXISTS" = true ]; then
    if [ "$FORCE_INSTALL" = true ] || [ "$IS_ENV_COMPLETE" = false ]; then
        echo -e "Updating existing '${ENV_NAME}' Conda environment from environment.yml..."
        $CONDA_CMD env update -n "$ENV_NAME" -f environment.yml --prune
    else
        echo -e "${GREEN}✔ Environment '${ENV_NAME}' is already installed and verified. Skipping re-creation.${NC}"
        echo -e "${DIM}(Use './install --force' to force re-updating all dependencies)${RESET}"
    fi
else
    echo -e "Creating new '${ENV_NAME}' Conda environment from environment.yml..."
    $CONDA_CMD env create -f environment.yml
fi

echo -e "\n${GREEN}====================================================${NC}"
echo -e "${GREEN}  GenSplice-Agent One-Click Setup Completed! 🎉    ${NC}"
echo -e "${GREEN}====================================================${NC}"
echo -e "Next steps:"
echo -e "  1. Activate Conda Environment:"
echo -e "     ${CYAN}source \$HOME/miniforge3/etc/profile.d/conda.sh && conda activate ${ENV_NAME}${NC}"
echo -e "  2. Test Pipeline:"
echo -e "     ${CYAN}./test${NC}"
echo -e "  3. Download Reference Genome & Run Analysis:"
echo -e "     ${CYAN}./ref human && ./GenSplice${NC}"
echo -e ""
