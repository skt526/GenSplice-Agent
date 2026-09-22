#!/usr/bin/env bash
# ==============================================================================
# GenSplice-Agent System Pre-flight Check & One-Click Conda Installer
# ==============================================================================

set -e

# ANSI Color Definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[1;36m'
NC='\033[0m' # No Color

ENV_NAME="gensplice-agent"
FORCE_INSTALL=false

for arg in "$@"; do
    if [ "$arg" = "--force" ] || [ "$arg" = "-f" ]; then
        FORCE_INSTALL=true
    fi
done

echo -e "${BLUE}====================================================${NC}"
echo -e "${BLUE}     GenSplice-Agent: System Pre-flight Check       ${NC}"
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
echo -e "CPU Number (Threads) : ${YELLOW}${CPU_CORES}${NC}"
echo -e "Total RAM Size       : ${YELLOW}${RAM_TOTAL_GB} GB${NC}"
echo -e "Available Disk Space : ${YELLOW}${DISK_FREE_GB} GB${NC}"
echo -e "----------------------------------------------------"

# 3. Minimum Requirements Thresholds
# - STAR Human/Mouse alignment RAM requirement: >= 30 GB
# - Free Disk requirement for Index & BAM files: >= 50 GB
# - Minimum CPU Cores: >= 4
MIN_RAM_GB=30
MIN_DISK_GB=50
MIN_CPU=4

IS_SUITABLE=true
FAIL_REASONS=()

# RAM Validation
if [ "$RAM_TOTAL_GB" -lt "$MIN_RAM_GB" ]; then
    IS_SUITABLE=false
    FAIL_REASONS+=("Insufficient RAM: Current ${RAM_TOTAL_GB} GB (Minimum required for STAR Human/Mouse alignment: ${MIN_RAM_GB} GB)")
fi

# Disk Space Validation
if [ "$DISK_FREE_GB" -lt "$MIN_DISK_GB" ]; then
    IS_SUITABLE=false
    FAIL_REASONS+=("Insufficient Disk Space: Current ${DISK_FREE_GB} GB (Minimum required for BAM/intermediates: ${MIN_DISK_GB} GB)")
fi

# CPU Cores Validation
if [ "$CPU_CORES" -lt "$MIN_CPU" ]; then
    IS_SUITABLE=false
    FAIL_REASONS+=("Insufficient CPU Cores: Current ${CPU_CORES} threads (Minimum recommended: ${MIN_CPU} cores)")
fi

# 4. Gatekeeper Evaluation
if [ "$IS_SUITABLE" = false ]; then
    echo -e "${RED}[UNSUITABLE: Installation Aborted]${NC}"
    echo -e "${RED}Current system specs do not meet the minimum requirements for GenSplice-Agent:${NC}"
    for reason in "${FAIL_REASONS[@]}"; do
        echo -e "  - ${RED}${reason}${NC}"
    done
    
    if [ "$FORCE_INSTALL" = true ]; then
        echo -e "\n${YELLOW}⚠ Warning: --force flag detected. Proceeding despite system requirement warnings...${NC}"
    else
        echo -e "\nInstallation aborted. Please expand system resources or use '--force' flag to override."
        exit 1
    fi
else
    echo -e "${GREEN}[SUITABLE: Proceeding with Installation]${NC}"
    echo -e "System specifications are suitable for large-scale transcriptomics and splicing analysis."
fi

# 5. Calculate (n - 2) Optimal Thread Allocation
if [ "$CPU_CORES" -gt 2 ]; then
    OPTIMAL_THREADS=$((CPU_CORES - 2))
else
    OPTIMAL_THREADS=1
fi

echo -e "----------------------------------------------------"
echo -e "System Stability Optimization:"
echo -e "  Allocating ${GREEN}${OPTIMAL_THREADS}${NC} threads (n - 2 out of ${CPU_CORES} total)"
echo -e "  to prevent background OS/GUI freezing or screen stuttering."
echo -e "----------------------------------------------------"

# 6. Auto-generate / Update config.yaml with Hardware & Thread Settings
cat <<EOF > config.yaml
# GenSplice-Agent Auto-generated Configuration
system:
  total_cpu: ${CPU_CORES}
  assigned_threads: ${OPTIMAL_THREADS}
  total_ram_gb: ${RAM_TOTAL_GB}

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

echo -e "Configuration file 'config.yaml' created with assigned_threads=${OPTIMAL_THREADS}.\n"

# 7. Package Manager & Conda Environment Setup
echo -e "${BLUE}[1/2] Creating directory structure...${NC}"
mkdir -p inputs/control inputs/treatment outputs

echo -e "${BLUE}[2/2] Detecting Package Manager & Setting up Conda Environment...${NC}"
CONDA_CMD=""
if command -v mamba &>/dev/null; then
    CONDA_CMD="mamba"
elif command -v conda &>/dev/null; then
    CONDA_CMD="conda"
elif command -v micromamba &>/dev/null; then
    CONDA_CMD="micromamba"
fi

if [ -n "$CONDA_CMD" ]; then
    echo -e "Found Conda command: ${CONDA_CMD}"
    
    if $CONDA_CMD env list | grep -qE "^${ENV_NAME}\s"; then
        echo -e "Updating existing '${ENV_NAME}' Conda environment..."
        $CONDA_CMD env update -n "$ENV_NAME" -f environment.yml --prune
    else
        echo -e "Creating new '${ENV_NAME}' Conda environment from environment.yml..."
        $CONDA_CMD env create -f environment.yml
    fi
    
    echo -e "\n${GREEN}====================================================${NC}"
    echo -e "${GREEN}  GenSplice-Agent Installation Completed Successfully! ${NC}"
    echo -e "${GREEN}====================================================${NC}"
    echo -e "Next steps:"
    echo -e "  - For Testing:  ./install ➔ ./test"
    echo -e "  - For Analysis: conda activate ${ENV_NAME} ➔ ./ref human ➔ ./GenSplice"
else
    echo -e "${YELLOW}Warning: Conda was not detected in PATH. Please install Conda/Miniforge to manage bioinformatics binaries.${NC}"
fi
