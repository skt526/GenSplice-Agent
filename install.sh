#!/usr/bin/env bash
# ==============================================================================
# GenSplice-Agent One-Click Dependency Installer
# ==============================================================================

set -e

# Color definitions for visual progress
BOLD="\033[1m"
GREEN="\033[1;32m"
CYAN="\033[1;36m"
YELLOW="\033[1;33m"
RED="\033[1;31m"
RESET="\033[0m"
DIM="\033[2m"

# Print banner
echo -e "${CYAN}${BOLD}"
echo "======================================================================"
echo "           GenSplice-Agent Environment Installer                      "
echo "======================================================================"
echo -e "${RESET}"

# Step 1: Detect Python 3
echo -e "${BOLD}[1/5] Checking Python 3 environment...${RESET}"
if command -v python3 &>/dev/null; then
    PYTHON_BIN="python3"
elif command -v python &>/dev/null; then
    PYTHON_BIN="python"
else
    echo -e "${RED}✘ Error: Python 3 is not installed or not found in PATH.${RESET}"
    echo "Please install Python 3.9+ and try again."
    exit 1
fi

PY_VERSION=$($PYTHON_BIN -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}")')
echo -e "  ${GREEN}✔ Found Python ${PY_VERSION} (${PYTHON_BIN})${RESET}"

# Step 2: Set up Virtual Environment
VENV_DIR=".venv"
echo -e "\n${BOLD}[2/5] Setting up Virtual Environment (${VENV_DIR})...${RESET}"
if [ ! -d "$VENV_DIR" ]; then
    echo -e "  ${DIM}Creating virtual environment in ./${VENV_DIR}...${RESET}"
    $PYTHON_BIN -m venv "$VENV_DIR"
    echo -e "  ${GREEN}✔ Virtual environment created successfully.${RESET}"
else
    echo -e "  ${GREEN}✔ Virtual environment already exists at ./${VENV_DIR}.${RESET}"
fi

# Activate virtual environment
source "$VENV_DIR/bin/activate"

# Step 3: Upgrade Pip
echo -e "\n${BOLD}[3/5] Upgrading package manager (pip)...${RESET}"
pip install --upgrade pip setuptools wheel --quiet
echo -e "  ${GREEN}✔ pip is up to date.${RESET}"

# Step 4: Install Required Packages
echo -e "\n${BOLD}[4/5] Installing core dependencies from requirements.txt...${RESET}"
if [ -f "requirements.txt" ]; then
    # Install with progress indication
    pip install -r requirements.txt --progress-bar on
    echo -e "  ${GREEN}✔ All packages installed successfully.${RESET}"
else
    echo -e "${RED}✘ Error: requirements.txt file not found.${RESET}"
    exit 1
fi

# Step 5: Verification Check
echo -e "\n${BOLD}[5/5] Verifying installed packages...${RESET}"
PACKAGES=("streamlit" "polars" "pyarrow" "plotly" "pandas" "google.genai")

for pkg in "${PACKAGES[@]}"; do
    if python -c "import $pkg" &>/dev/null; then
        echo -e "  ${GREEN}✔ Module '$pkg' verified.${RESET}"
    else
        echo -e "  ${RED}✘ Warning: Module '$pkg' import failed.${RESET}"
    fi
done

# Completion Notice
echo -e "\n${CYAN}${BOLD}======================================================================"
echo "           GenSplice-Agent Setup Complete! 🎉                         "
echo "======================================================================"
echo -e "${RESET}"
echo -e "To start using GenSplice-Agent:"
echo -e "  ${YELLOW}1. Activate environment:${RESET} source .venv/bin/activate"
echo -e "  ${YELLOW}2. Run dashboard:${RESET}       streamlit run app.py"
echo -e ""
