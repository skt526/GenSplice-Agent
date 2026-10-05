#!/usr/bin/env bash
# ==============================================================================
# GenSplice-Agent: GSE190087 (PRJNA786720 / PRJNA785817) Dataset Downloader
# ==============================================================================
# Biological Study:
#   "SF3B1 mutation–mediated sensitization to H3B-8800 splicing inhibitor in CLL"
#   Life Science Alliance (2023) | PMID: 36723485
#
# Features:
#   - Direct HTTP/FTP download from ENA mirrors (Resumable & Fast)
#   - Supports aria2c, curl, and wget
#   - Interactive menu or direct CLI mode flags
#   - Auto-sorts paired-end reads into inputs/control and inputs/treatment
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

PYTHON_BIN=""
for candidate in \
    "$HOME/miniforge3/envs/gensplice-agent/bin/python3" \
    "$HOME/miniforge3/envs/gensplice-agent/bin/python" \
    "$HOME/miniconda3/envs/gensplice-agent/bin/python3" \
    "$HOME/miniconda3/envs/gensplice-agent/bin/python" \
    "/home/skt526/miniconda3/envs/gensplice-agent/bin/python3" \
    "/home/skt526/miniconda3/envs/gensplice-agent/bin/python" \
    "/root/miniconda3/envs/gensplice-agent/bin/python3" \
    "/root/miniconda3/envs/gensplice-agent/bin/python" \
    "/opt/conda/envs/gensplice-agent/bin/python3" \
    "/opt/conda/envs/gensplice-agent/bin/python" \
    "$CONDA_PREFIX/bin/python3" \
    "$CONDA_PREFIX/bin/python" \
    ".venv/bin/python3" \
    ".venv/bin/python" \
    /opt/homebrew/bin/python3 \
    /usr/bin/python3 \
    python3; do
    if [ -n "$candidate" ] && "$candidate" --version >/dev/null 2>&1; then
        PYTHON_BIN="$candidate"
        if [ "$PYTHON_BIN" != "python3" ] && [[ "$PYTHON_BIN" == */* ]]; then
            ENV_BIN="$(dirname "$PYTHON_BIN")"
            export PATH="$ENV_BIN:$PATH"
            export CONDA_PREFIX="$(dirname "$ENV_BIN")"
            export CONDA_DEFAULT_ENV="gensplice-agent"
        fi
        break
    fi
done

if [ -z "$PYTHON_BIN" ]; then
    PYTHON_BIN="python3"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="$SCRIPT_DIR/download_GSE190087.py"

show_menu() {
    echo -e "${CYAN}${BOLD}"
    echo "======================================================================"
    echo "       GenSplice-Agent: GSE190087 CLL SF3B1 Dataset Downloader        "
    echo "======================================================================"
    echo -e "${RESET}"
    echo -e "${BOLD}Study:${RESET} SF3B1 Mutation & Splicing Modulator H3B-8800 in CLL"
    echo -e "${BOLD}Target Organism:${RESET} Homo sapiens (Human MEC-1 cell line, 150bp PE)"
    echo -e "${BOLD}Required Reference:${RESET} GRCh38 (${CYAN}./ref human${RESET})"
    echo -e "----------------------------------------------------------------------"
    echo -e "Select dataset mode to download:"
    echo -e ""
    echo -e "  ${GREEN}1) Benchmark (Recommended)${RESET}"
    echo -e "     Control: 3x SF3B1 Wild-Type | Treatment: 3x SF3B1 K700E mutant"
    echo -e "     (~9.2 GB, full statistical power for DEG + rMATS splicing)"
    echo -e ""
    echo -e "  ${CYAN}2) Quick Test${RESET}"
    echo -e "     Control: 1x WT rep1 | Treatment: 1x K700E rep1"
    echo -e "     (~3.1 GB, fast end-to-end pipeline sanity test)"
    echo -e ""
    echo -e "  ${YELLOW}3) Drug Response in Wild-Type (WT DMSO vs WT + H3B-8800)${RESET}"
    echo -e "     Control: 3x WT Vehicle | Treatment: 3x WT H3B-8800 inhibitor"
    echo -e "     (~9.0 GB, test splicing modulator effect on wild-type cells)"
    echo -e ""
    echo -e "  ${YELLOW}4) Drug Response in Mutant (K700E DMSO vs K700E + H3B-8800)${RESET}"
    echo -e "     Control: 3x K700E Vehicle | Treatment: 3x K700E H3B-8800 inhibitor"
    echo -e "     (~8.9 GB, test compound synergy on mutant cells)"
    echo -e ""
    echo -e "  ${DIM}5) Dry Run (Inspect download URLs and destination paths)${RESET}"
    echo -e ""
    echo -e "  ${RED}6) Exit${RESET}"
    echo -e "----------------------------------------------------------------------"
    read -rp "Enter choice [1-6] (Default: 1): " choice
    choice=${choice:-1}

    case "$choice" in
        1)
            "$PYTHON_BIN" "$PYTHON_SCRIPT" --mode benchmark
            ;;
        2)
            "$PYTHON_BIN" "$PYTHON_SCRIPT" --mode quick
            ;;
        3)
            "$PYTHON_BIN" "$PYTHON_SCRIPT" --mode drug_wt
            ;;
        4)
            "$PYTHON_BIN" "$PYTHON_SCRIPT" --mode drug_k700e
            ;;
        5)
            "$PYTHON_BIN" "$PYTHON_SCRIPT" --mode benchmark --dry-run
            ;;
        6)
            echo -e "${DIM}Download cancelled.${RESET}"
            exit 0
            ;;
        *)
            echo -e "${RED}Invalid option: $choice${RESET}"
            exit 1
            ;;
    esac
}

# If arguments provided, pass directly to python script
if [ $# -gt 0 ]; then
    case "$1" in
        -h|--help|help)
            echo -e "Usage: ./download_GSE190087.sh [mode|options]"
            echo -e ""
            echo -e "Available modes:"
            echo -e "  ${GREEN}benchmark${RESET}   : 3 WT vs 3 K700E (Default biological benchmark)"
            echo -e "  ${CYAN}quick${RESET}       : 1 WT vs 1 K700E (Fast pipeline test)"
            echo -e "  ${YELLOW}drug_wt${RESET}     : 3 WT vs 3 WT+H3B-8800"
            echo -e "  ${YELLOW}drug_k700e${RESET}  : 3 K700E vs 3 K700E+H3B-8800"
            echo -e "  ${DIM}dry-run${RESET}     : Check URLs and paths without downloading"
            echo -e ""
            echo -e "Examples:"
            echo -e "  ./download_GSE190087.sh benchmark"
            echo -e "  ./download_GSE190087.sh quick"
            echo -e "  ./download_GSE190087.sh --dry-run"
            exit 0
            ;;
        benchmark|quick|drug_wt|drug_k700e|all)
            MODE="$1"
            shift
            "$PYTHON_BIN" "$PYTHON_SCRIPT" --mode "$MODE" "$@"
            ;;
        --dry-run)
            "$PYTHON_BIN" "$PYTHON_SCRIPT" --mode benchmark --dry-run
            ;;
        *)
            "$PYTHON_BIN" "$PYTHON_SCRIPT" "$@"
            ;;
    esac
else
    # Interactive menu
    show_menu
fi
