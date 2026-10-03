#!/usr/bin/env bash
# ==============================================================================
# GenSplice-Agent Public Stress Dataset Downloader (download_data.sh)
# Downloads authentic Human Heat Shock RNA-seq FASTQ files (Control vs Stress)
# Directly from ENA / NCBI SRA high-speed HTTPS endpoints.
# ==============================================================================

set -e

# Terminal Colors
CYAN="\033[1;36m"
GREEN="\033[1;32m"
YELLOW="\033[1;33m"
RED="\033[1;31m"
BOLD="\033[1m"
DIM="\033[2m"
RESET="\033[0m"

echo -e "${CYAN}${BOLD}"
echo "======================================================================"
echo "      GenSplice-Agent RNA-seq Stress Dataset Downloader              "
echo "======================================================================"
echo -e "${RESET}"

show_help() {
    echo -e "Usage: ./download_data.sh [GSE140053|GSE197993]"
    echo -e ""
    echo -e "Available benchmark datasets:"
    echo -e "  ${CYAN}GSE140053${RESET} : Human Heat Shock (Lightweight / Quick Test: ~3.3 GB total)"
    echo -e "               Cell line: HeLa / HEK293 (WT 37°C Control vs 42°C Heat Shock 1h)"
    echo -e "               Replicates: 3 Control vs 3 Treatment (~10M-16M reads/sample)"
    echo -e "  ${CYAN}GSE197993${RESET} : Human Heat Shock (Deep Sequencing: ~18.1 GB total)"
    echo -e "               Cell line: Human cells (siNT 37°C Control vs 42°C Heat Shock)"
    echo -e "               Replicates: 3 Control vs 3 Treatment (~50M-80M reads/sample)"
    echo -e ""
    echo -e "If no argument is given, an interactive selection menu will be prompted."
    exit 0
}

if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
    show_help
fi

DATASET_CHOICE="$1"

# Interactive Selection
if [ -z "$DATASET_CHOICE" ]; then
    echo -e "${BOLD}Select a Human Heat Shock dataset to download into inputs/:${RESET}\n"
    echo -e "  ${CYAN}1)${RESET} GSE140053 (Lightweight / Rapid Benchmark - ~3.3 GB total, 6 samples)"
    echo -e "     ${DIM}Paper: Stress-induced nuclear condensation of NELF (Mol Cell 2021)${RESET}"
    echo -e "  ${CYAN}2)${RESET} GSE197993 (Deep Sequencing Benchmark - ~18.1 GB total, 6 samples)"
    echo -e "     ${DIM}Paper: RPRD1B required for heat shock response (2022)${RESET}"
    echo -e "  ${CYAN}3)${RESET} Exit\n"
    read -p "Enter choice [1-3]: " MENU_CHOICE

    case $MENU_CHOICE in
        1) DATASET_CHOICE="GSE140053" ;;
        2) DATASET_CHOICE="GSE197993" ;;
        3) echo "Exiting..."; exit 0 ;;
        *) echo -e "${RED}Invalid selection.${RESET}"; exit 1 ;;
    esac
fi

DATASET_CHOICE=$(echo "$DATASET_CHOICE" | tr '[:lower:]' '[:upper:]')

mkdir -p inputs/control inputs/treatment

download_file() {
    local url="$1"
    local target="$2"
    local desc="$3"
    local srr="$4"

    echo -e "\n${BOLD}${desc}${RESET} ${DIM}(SRA: ${srr})${RESET}"
    echo -e "  Target: ${GREEN}${target}${RESET}"

    if [ -f "$target" ]; then
        if gzip -t "$target" 2>/dev/null; then
            echo -e "  ${GREEN}✔ File already exists and verified intact. Skipping download.${RESET}"
            return 0
        else
            echo -e "  ${YELLOW}⚠ Existing file is corrupted. Re-downloading...${RESET}"
            rm -f "$target"
        fi
    fi

    echo -e "  Downloading from: ${DIM}${url}${RESET}"
    curl -L -C - --fail --retry 3 --progress-bar "$url" -o "$target"

    echo -e "  Verifying gzip integrity..."
    if gzip -t "$target" 2>/dev/null; then
        echo -e "  ${GREEN}✔ Verified successfully.${RESET}"
    else
        echo -e "  ${RED}✘ Verification failed for ${target}. Please re-run to resume/repair.${RESET}"
        exit 1
    fi
}

case "$DATASET_CHOICE" in
    GSE140053|1)
        echo -e "${BOLD}Target Dataset:${RESET} ${GREEN}GSE140053 (Human Heat Shock WT - ~3.3 GB)${RESET}"
        echo -e "${DIM}Study: Stress-induced nuclear condensation of NELF drives transcriptional downregulation${RESET}"
        echo -e "${DIM}Conditions: Normal 37°C (Control, n=3) vs Heat Shock 42°C (Treatment, n=3)${RESET}\n"

        # Control Samples (WT 37°C Non-Heat Shock)
        download_file \
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/018/SRR12801718/SRR12801718.fastq.gz" \
            "inputs/control/GSE140053_Control_Rep1.fastq.gz" \
            "[1/6] Downloading Control Rep 1" \
            "SRR12801718"

        download_file \
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/019/SRR12801719/SRR12801719.fastq.gz" \
            "inputs/control/GSE140053_Control_Rep2.fastq.gz" \
            "[2/6] Downloading Control Rep 2" \
            "SRR12801719"

        download_file \
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/020/SRR12801720/SRR12801720.fastq.gz" \
            "inputs/control/GSE140053_Control_Rep3.fastq.gz" \
            "[3/6] Downloading Control Rep 3" \
            "SRR12801720"

        # Treatment Samples (WT 42°C Heat Shock 1h)
        download_file \
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/012/SRR12801712/SRR12801712.fastq.gz" \
            "inputs/treatment/GSE140053_HeatShock_Rep1.fastq.gz" \
            "[4/6] Downloading Heat Shock Rep 1" \
            "SRR12801712"

        download_file \
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/013/SRR12801713/SRR12801713.fastq.gz" \
            "inputs/treatment/GSE140053_HeatShock_Rep2.fastq.gz" \
            "[5/6] Downloading Heat Shock Rep 2" \
            "SRR12801713"

        download_file \
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/014/SRR12801714/SRR12801714.fastq.gz" \
            "inputs/treatment/GSE140053_HeatShock_Rep3.fastq.gz" \
            "[6/6] Downloading Heat Shock Rep 3" \
            "SRR12801714"
        ;;

    GSE197993|2)
        echo -e "${BOLD}Target Dataset:${RESET} ${GREEN}GSE197993 (Human Heat Shock Deep - ~18.1 GB)${RESET}"
        echo -e "${DIM}Study: RPRD1B is required for proper heat shock response${RESET}"
        echo -e "${DIM}Conditions: Normal 37°C (Control, n=3) vs Heat Shock 42°C (Treatment, n=3)${RESET}\n"

        # Control Samples (siNT 37°C Non-Heat Shock)
        download_file \
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR182/025/SRR18248925/SRR18248925.fastq.gz" \
            "inputs/control/GSE197993_Control_Rep1.fastq.gz" \
            "[1/6] Downloading Control Rep 1" \
            "SRR18248925"

        download_file \
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR182/021/SRR18248921/SRR18248921.fastq.gz" \
            "inputs/control/GSE197993_Control_Rep2.fastq.gz" \
            "[2/6] Downloading Control Rep 2" \
            "SRR18248921"

        download_file \
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR182/017/SRR18248917/SRR18248917.fastq.gz" \
            "inputs/control/GSE197993_Control_Rep3.fastq.gz" \
            "[3/6] Downloading Control Rep 3" \
            "SRR18248917"

        # Treatment Samples (siNT 42°C Heat Shock)
        download_file \
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR182/024/SRR18248924/SRR18248924.fastq.gz" \
            "inputs/treatment/GSE197993_HeatShock_Rep1.fastq.gz" \
            "[4/6] Downloading Heat Shock Rep 1" \
            "SRR18248924"

        download_file \
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR182/020/SRR18248920/SRR18248920.fastq.gz" \
            "inputs/treatment/GSE197993_HeatShock_Rep2.fastq.gz" \
            "[5/6] Downloading Heat Shock Rep 2" \
            "SRR18248920"

        download_file \
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR182/016/SRR18248916/SRR18248916.fastq.gz" \
            "inputs/treatment/GSE197993_HeatShock_Rep3.fastq.gz" \
            "[6/6] Downloading Heat Shock Rep 3" \
            "SRR18248916"
        ;;

    *)
        echo -e "${RED}Error: Unknown dataset '${DATASET_CHOICE}'.${RESET}"
        echo -e "Available options: GSE140053 or GSE197993"
        exit 1
        ;;
esac

echo -e "\n${CYAN}${BOLD}======================================================================"
echo "           Dataset Download & Verification Complete! 🎉               "
echo "======================================================================"
echo -e "${RESET}"
echo -e "Samples placed in ${GREEN}inputs/${RESET}:"
echo -e "  - ${BOLD}Control:${RESET}   $(ls -1 inputs/control/*.fastq.gz 2>/dev/null | wc -l | tr -d ' ') files in inputs/control/"
echo -e "  - ${BOLD}Treatment:${RESET} $(ls -1 inputs/treatment/*.fastq.gz 2>/dev/null | wc -l | tr -d ' ') files in inputs/treatment/"
echo -e ""
echo -e "Next Steps:"
echo -e "  1. If human reference genome is not built yet:  ${CYAN}./ref human${RESET}"
echo -e "  2. Run the automated GenSplice pipeline:        ${CYAN}./GenSplice.sh${RESET}"
echo -e ""
