#!/usr/bin/env bash
# ==============================================================================
# GenSplice-Agent Conda/Mamba One-Click Installer
# ==============================================================================

set -e

# Color definitions
BOLD="\033[1m"
GREEN="\033[1;32m"
CYAN="\033[1;36m"
YELLOW="\033[1;33m"
RED="\033[1;31m"
RESET="\033[0m"
DIM="\033[2m"

ENV_NAME="gensplice-agent"

echo -e "${CYAN}${BOLD}"
echo "======================================================================"
echo "         GenSplice-Agent Conda Environment Installer                  "
echo "======================================================================"
echo -e "${RESET}"

# Step 1: Detect Conda / Mamba / Micromamba
echo -e "${BOLD}[1/5] Detecting package manager (Mamba / Conda)...${RESET}"
CONDA_CMD=""
if command -v mamba &>/dev/null; then
    CONDA_CMD="mamba"
elif command -v conda &>/dev/null; then
    CONDA_CMD="conda"
elif command -v micromamba &>/dev/null; then
    CONDA_CMD="micromamba"
else
    echo -e "  ${YELLOW}⚠ Warning: Neither Conda nor Mamba was found in PATH.${RESET}"
    echo -e "  Falling back to Python virtualenv setup..."
    
    # Fallback to python venv if conda is not installed
    if command -v python3 &>/dev/null; then
        PYTHON_BIN="python3"
    else
        PYTHON_BIN="python"
    fi
    
    VENV_DIR=".venv"
    if [ ! -d "$VENV_DIR" ]; then
        $PYTHON_BIN -m venv "$VENV_DIR"
    fi
    source "$VENV_DIR/bin/activate"
    pip install --upgrade pip setuptools wheel --quiet
    pip install -r requirements.txt --quiet
    echo -e "  ${GREEN}✔ Python venv environment created at ./${VENV_DIR}.${RESET}"
    echo -e "  ${YELLOW}Note: To use bioinformatics tools (fastp, STAR, rmats), please install Conda/Miniforge.${RESET}"
    exit 0
fi

echo -e "  ${GREEN}✔ Found Package Manager: ${CONDA_CMD}${RESET}"

# Step 2: Check if Conda environment exists
echo -e "\n${BOLD}[2/5] Checking Conda environment ('${ENV_NAME}')...${RESET}"
ENV_EXISTS=false
if $CONDA_CMD env list | grep -qE "^${ENV_NAME}\s"; then
    ENV_EXISTS=true
fi

# Step 3: Create or Update Conda Environment
echo -e "\n${BOLD}[3/5] Installing/Updating Conda environment from environment.yml...${RESET}"
if [ "$ENV_EXISTS" = true ]; then
    echo -e "  ${DIM}Updating existing '${ENV_NAME}' environment...${RESET}"
    $CONDA_CMD env update -n "$ENV_NAME" -f environment.yml --prune
else
    echo -e "  ${DIM}Creating new '${ENV_NAME}' environment (this may take a few minutes)...${RESET}"
    $CONDA_CMD env create -f environment.yml
fi
echo -e "  ${GREEN}✔ Conda environment '${ENV_NAME}' is ready.${RESET}"

# Step 4: Verification of CLI & Python dependencies
echo -e "\n${BOLD}[4/5] Verifying installed CLI tools & Python dependencies...${RESET}"

# Get conda env bin path
CONDA_BASE=$($CONDA_CMD info --base 2>/dev/null || echo "$HOME/miniconda3")
ENV_BIN="${CONDA_BASE}/envs/${ENV_NAME}/bin"

# If conda eval available
if [ -d "$ENV_BIN" ]; then
    export PATH="${ENV_BIN}:$PATH"
fi

CLI_TOOLS=("fastp" "STAR" "rmats.py" "featureCounts" "R" "python")
for tool in "${CLI_TOOLS[@]}"; do
    if command -v "$tool" &>/dev/null || [ -x "${ENV_BIN}/${tool}" ]; then
        echo -e "  ${GREEN}✔ CLI Tool '$tool' verified.${RESET}"
    else
        echo -e "  ${YELLOW}⚠ Warning: CLI Tool '$tool' not directly executable (will be available inside conda env).${RESET}"
    fi
done

# Step 5: Verification of Python packages
echo -e "\n${BOLD}[5/5] Verifying Python dashboard dependencies...${RESET}"
PYTHON_ENV_BIN="${ENV_BIN}/python"
if [ ! -x "$PYTHON_ENV_BIN" ]; then
    PYTHON_ENV_BIN="python"
fi

PY_PACKAGES=("streamlit" "polars" "pyarrow" "plotly" "pandas" "google.genai")
for pkg in "${PY_PACKAGES[@]}"; do
    if $PYTHON_ENV_BIN -c "import $pkg" &>/dev/null; then
        echo -e "  ${GREEN}✔ Python Package '$pkg' verified.${RESET}"
    else
        echo -e "  ${RED}✘ Warning: Python Package '$pkg' import failed.${RESET}"
    fi
done

# Completion Notice
echo -e "\n${CYAN}${BOLD}======================================================================"
echo "           GenSplice-Agent Setup Complete! 🎉                         "
echo "======================================================================"
echo -e "${RESET}"
echo -e "To activate the environment and run GenSplice-Agent:"
echo -e "  ${YELLOW}1. Activate Conda Env:${RESET} conda activate ${ENV_NAME}"
echo -e "  ${YELLOW}2. Run Dashboard:${RESET}     streamlit run app.py"
echo -e ""
