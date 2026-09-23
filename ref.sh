#!/usr/bin/env bash
# ==============================================================================
# GenSplice-Agent Reference Genome Downloader
# ==============================================================================

set -e

# Color definitions for terminal output
BOLD="\033[1m"
GREEN="\033[1;32m"
CYAN="\033[1;36m"
YELLOW="\033[1;33m"
RED="\033[1;31m"
RESET="\033[0m"
DIM="\033[2m"

# Print Header Banner
echo -e "${CYAN}${BOLD}"
echo "======================================================================"
echo "           GenSplice-Agent Reference Genome Downloader                "
echo "======================================================================"
echo -e "${RESET}"

# Function to display help
show_help() {
    echo -e "Usage: ./ref [organism_key]"
    echo -e ""
    echo -e "Available organism keys:"
    echo -e "  human       - Human (Homo sapiens - GRCh38 / Ensembl)"
    echo -e "  mouse       - Mouse (Mus musculus - GRCm39 / Ensembl)"
    echo -e "  rice        - Rice (Oryza sativa - IRGSP-1.0 / Ensembl Plants)"
    echo -e "  arabidopsis - Arabidopsis (Arabidopsis thaliana - TAIR10 / Ensembl Plants)"
    echo -e "  maize       - Maize/Corn (Zea mays - Zm-B73-NAM-5.0 / Ensembl Plants)"
    echo -e ""
    echo -e "If no argument is passed, an interactive menu will be displayed."
    exit 0
}

if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
    show_help
fi

ORGANISM_CHOICE="$1"

# Interactive Selection Menu if argument not provided
if [ -z "$ORGANISM_CHOICE" ]; then
    echo -e "${BOLD}Please select a Reference Genome to download:${RESET}"
    echo -e "  ${CYAN}1)${RESET} Human (Homo sapiens - GRCh38)"
    echo -e "  ${CYAN}2)${RESET} Mouse (Mus musculus - GRCm39)"
    echo -e "  ${CYAN}3)${RESET} Rice (Oryza sativa - IRGSP-1.0)"
    echo -e "  ${CYAN}4)${RESET} Arabidopsis thaliana (TAIR10)"
    echo -e "  ${CYAN}5)${RESET} Maize / Corn (Zea mays - B73)"
    echo -e "  ${CYAN}6)${RESET} Exit"
    echo -e ""
    read -p "Enter choice [1-6]: " MENU_SELECTION

    case $MENU_SELECTION in
        1) ORGANISM_CHOICE="human" ;;
        2) ORGANISM_CHOICE="mouse" ;;
        3) ORGANISM_CHOICE="rice" ;;
        4) ORGANISM_CHOICE="arabidopsis" ;;
        5) ORGANISM_CHOICE="maize" ;;
        6) echo "Exiting..."; exit 0 ;;
        *) echo -e "${RED}Invalid selection.${RESET}"; exit 1 ;;
    esac
fi

# Convert choice to lowercase
ORGANISM_CHOICE=$(echo "$ORGANISM_CHOICE" | tr '[:upper:]' '[:lower:]')

# Define Genome Specs based on Organism Choice
case "$ORGANISM_CHOICE" in
    human)
        ORGANISM_NAME="human"
        SPECIES_DESC="Homo sapiens (GRCh38)"
        FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz"
        FASTA_FILE="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        GTF_FILE="Homo_sapiens.GRCh38.113.gtf.gz"
        ;;
    mouse|mice)
        ORGANISM_NAME="mouse"
        SPECIES_DESC="Mus musculus (GRCm39)"
        FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
        GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz"
        FASTA_FILE="Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
        GTF_FILE="Mus_musculus.GRCm39.113.gtf.gz"
        ;;
    rice)
        ORGANISM_NAME="rice"
        SPECIES_DESC="Oryza sativa (IRGSP-1.0)"
        FASTA_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz"
        GTF_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/gtf/oryza_sativa/Oryza_sativa.IRGSP-1.0.60.gtf.gz"
        FASTA_FILE="Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz"
        GTF_FILE="Oryza_sativa.IRGSP-1.0.60.gtf.gz"
        ;;
    arabidopsis)
        ORGANISM_NAME="arabidopsis"
        SPECIES_DESC="Arabidopsis thaliana (TAIR10)"
        FASTA_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
        GTF_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.60.gtf.gz"
        FASTA_FILE="Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
        GTF_FILE="Arabidopsis_thaliana.TAIR10.60.gtf.gz"
        ;;
    maize|corn)
        ORGANISM_NAME="maize"
        SPECIES_DESC="Zea mays (Zm-B73-REFERENCE-NAM-5.0)"
        FASTA_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz"
        GTF_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/gtf/zea_mays/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.60.gtf.gz"
        FASTA_FILE="Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz"
        GTF_FILE="Zea_mays.Zm-B73-REFERENCE-NAM-5.0.60.gtf.gz"
        ;;
    *)
        echo -e "${RED}Error: Unknown organism '${ORGANISM_CHOICE}'.${RESET}"
        echo -e "Available options: human, mouse, rice, arabidopsis, maize"
        exit 1
        ;;
esac

TARGET_DIR="./${ORGANISM_NAME}-ref"

echo -e "\n${BOLD}[1/4] Target Reference Folder:${RESET} ${GREEN}${TARGET_DIR}${RESET}"
echo -e "${DIM}Species: ${SPECIES_DESC}${RESET}"
echo -e "${YELLOW}${BOLD}⚠ Notice: Downloading reference genome FASTA/GTF (5GB+) and building STAR index may take 15 to 30+ minutes depending on your internet connection and hardware environment.${RESET}"

# Step 1: Create target directory
mkdir -p "$TARGET_DIR"

# Step 2: Download FASTA
echo -e "\n${BOLD}[2/4] Downloading Genome FASTA file...${RESET}"
echo -e "${DIM}Source: ${FASTA_URL}${RESET}"
if [ -f "${TARGET_DIR}/${FASTA_FILE}" ]; then
    echo -e "  ${GREEN}✔ FASTA file already exists in ${TARGET_DIR}.${RESET}"
else
    curl -L --progress-bar "$FASTA_URL" -o "${TARGET_DIR}/${FASTA_FILE}"
    echo -e "  ${GREEN}✔ Downloaded ${FASTA_FILE}.${RESET}"
fi

# Step 3: Download GTF
echo -e "\n${BOLD}[3/4] Downloading Gene Annotation GTF file...${RESET}"
echo -e "${DIM}Source: ${GTF_URL}${RESET}"
if [ -f "${TARGET_DIR}/${GTF_FILE}" ]; then
    echo -e "  ${GREEN}✔ GTF file already exists in ${TARGET_DIR}.${RESET}"
else
    curl -L --progress-bar "$GTF_URL" -o "${TARGET_DIR}/${GTF_FILE}"
    echo -e "  ${GREEN}✔ Downloaded ${GTF_FILE}.${RESET}"
fi

# Step 4: Decompressing files for STAR / rMATS readiness
echo -e "\n${BOLD}[4/4] Decompressing FASTA and GTF for STAR/rMATS indexing...${RESET}"
FASTA_UNCOMPRESSED="${FASTA_FILE%.gz}"
GTF_UNCOMPRESSED="${GTF_FILE%.gz}"

if [ ! -f "${TARGET_DIR}/${FASTA_UNCOMPRESSED}" ]; then
    echo -e "  ${DIM}Decompressing FASTA (${FASTA_FILE})...${RESET}"
    gzip -dc "${TARGET_DIR}/${FASTA_FILE}" > "${TARGET_DIR}/${FASTA_UNCOMPRESSED}"
    echo -e "  ${GREEN}✔ Created ${FASTA_UNCOMPRESSED}${RESET}"
else
    echo -e "  ${GREEN}✔ Decompressed FASTA already exists.${RESET}"
fi

if [ ! -f "${TARGET_DIR}/${GTF_UNCOMPRESSED}" ]; then
    echo -e "  ${DIM}Decompressing GTF (${GTF_FILE})...${RESET}"
    gzip -dc "${TARGET_DIR}/${GTF_FILE}" > "${TARGET_DIR}/${GTF_UNCOMPRESSED}"
    echo -e "  ${GREEN}✔ Created ${GTF_UNCOMPRESSED}${RESET}"
else
    echo -e "  ${GREEN}✔ Decompressed GTF already exists.${RESET}"
fi

# Save Metadata JSON
cat <<EOF > "${TARGET_DIR}/metadata.json"
{
  "organism": "${ORGANISM_NAME}",
  "species": "${SPECIES_DESC}",
  "download_timestamp": "$(date -u +"%Y-%m-%dT%H:%M:%SZ")",
  "fasta_url": "${FASTA_URL}",
  "gtf_url": "${GTF_URL}",
  "fasta_file": "${FASTA_UNCOMPRESSED}",
  "gtf_file": "${GTF_UNCOMPRESSED}"
}
EOF

# Auto-update config.yaml reference section
if [ -f "config.yaml" ]; then
    echo -e "\n${BOLD}[Auto-Config] Updating 'config.yaml' reference configuration...${RESET}"
    python3 -c "
import yaml
try:
    with open('config.yaml', 'r') as f:
        cfg = yaml.safe_load(f) or {}
    if 'reference' not in cfg:
        cfg['reference'] = {}
    cfg['reference']['organism'] = '${ORGANISM_NAME}'
    cfg['reference']['ref_dir'] = '${TARGET_DIR}'
    cfg['reference']['fasta'] = '${TARGET_DIR}/${FASTA_UNCOMPRESSED}'
    cfg['reference']['gtf'] = '${TARGET_DIR}/${GTF_UNCOMPRESSED}'
    cfg['reference']['star_index'] = '${TARGET_DIR}/star_index'
    with open('config.yaml', 'w') as f:
        yaml.dump(cfg, f, default_flow_style=False)
    print('  ✔ Automatically updated config.yaml to target reference: ${ORGANISM_NAME}')
except Exception as e:
    print('  ⚠ Notice: Could not update config.yaml:', e)
"
fi

# Completion Notice
echo -e "\n${CYAN}${BOLD}======================================================================"
echo "           Reference Download Complete! 🧬                            "
echo "======================================================================"
echo -e "${RESET}"
echo -e "Files saved in: ${GREEN}${TARGET_DIR}/${RESET}"
echo -e "  - FASTA: ${TARGET_DIR}/${FASTA_UNCOMPRESSED}"
echo -e "  - GTF:   ${TARGET_DIR}/${GTF_UNCOMPRESSED}"
echo -e "  - Meta:  ${TARGET_DIR}/metadata.json"
echo -e "  - Config: config.yaml updated"
echo -e ""
