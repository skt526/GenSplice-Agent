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
    echo -e "  ${BOLD}Mammals:${RESET}"
    echo -e "    human, mouse, rat, pig, cow, dog, macaque, chimpanzee"
    echo -e "  ${BOLD}Model Animals & Birds:${RESET}"
    echo -e "    zebrafish, drosophila, celegans, xenopus, chicken"
    echo -e "  ${BOLD}Plants:${RESET}"
    echo -e "    rice, arabidopsis, maize, wheat, soybean, tomato, potato, barley"
    echo -e "  ${BOLD}Fungi:${RESET}"
    echo -e "    yeast"
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
    echo -e "${BOLD}Please select a Reference Genome to download:${RESET}\n"
    echo -e "${YELLOW}--- Mammals (포유류) ---${RESET}"
    echo -e "  ${CYAN}1)${RESET} Human (Homo sapiens - GRCh38)"
    echo -e "  ${CYAN}2)${RESET} Mouse (Mus musculus - GRCm39)"
    echo -e "  ${CYAN}3)${RESET} Rat (Rattus norvegicus - mRatBN7.2)"
    echo -e "  ${CYAN}4)${RESET} Pig (Sus scrofa - Sscrofa11.1)"
    echo -e "  ${CYAN}5)${RESET} Cow (Bos taurus - ARS-UCD1.3)"
    echo -e "  ${CYAN}6)${RESET} Dog (Canis lupus familiaris - ROS_Cfam_1.0)"
    echo -e "  ${CYAN}7)${RESET} Macaque (Macaca mulatta - Mmul_10)"
    echo -e "  ${CYAN}8)${RESET} Chimpanzee (Pan troglodytes - Pan_tro_3.0)"
    echo -e "\n${YELLOW}--- Model Animals & Birds (대표 모델동물/조류) ---${RESET}"
    echo -e "  ${CYAN}9)${RESET} Zebrafish (Danio rerio - GRCz11)"
    echo -e " ${CYAN}10)${RESET} Drosophila / Fruit Fly (Drosophila melanogaster - BDGP6.46)"
    echo -e " ${CYAN}11)${RESET} C. elegans / Nematode (Caenorhabditis elegans - WBcel235)"
    echo -e " ${CYAN}12)${RESET} Xenopus / Frog (Xenopus tropicalis - UCB_Xtro_10.0)"
    echo -e " ${CYAN}13)${RESET} Chicken (Gallus gallus - GRCg7b)"
    echo -e "\n${YELLOW}--- Plants (식물) ---${RESET}"
    echo -e " ${CYAN}14)${RESET} Rice (Oryza sativa - IRGSP-1.0)"
    echo -e " ${CYAN}15)${RESET} Arabidopsis thaliana (TAIR10)"
    echo -e " ${CYAN}16)${RESET} Maize / Corn (Zea mays - B73)"
    echo -e " ${CYAN}17)${RESET} Wheat (Triticum aestivum - IWGSC)"
    echo -e " ${CYAN}18)${RESET} Soybean (Glycine max - v2.1)"
    echo -e " ${CYAN}19)${RESET} Tomato (Solanum lycopersicum - SL3.0)"
    echo -e " ${CYAN}20)${RESET} Potato (Solanum tuberosum - SolTub_3.0)"
    echo -e " ${CYAN}21)${RESET} Barley (Hordeum vulgare - MorexV3)"
    echo -e "\n${YELLOW}--- Fungi (균류) ---${RESET}"
    echo -e " ${CYAN}22)${RESET} Yeast (Saccharomyces cerevisiae - R64-1-1)"
    echo -e "\n  ${CYAN}23)${RESET} Exit"
    echo -e ""
    read -p "Enter choice [1-23]: " MENU_SELECTION

    case $MENU_SELECTION in
        1) ORGANISM_CHOICE="human" ;;
        2) ORGANISM_CHOICE="mouse" ;;
        3) ORGANISM_CHOICE="rat" ;;
        4) ORGANISM_CHOICE="pig" ;;
        5) ORGANISM_CHOICE="cow" ;;
        6) ORGANISM_CHOICE="dog" ;;
        7) ORGANISM_CHOICE="macaque" ;;
        8) ORGANISM_CHOICE="chimpanzee" ;;
        9) ORGANISM_CHOICE="zebrafish" ;;
        10) ORGANISM_CHOICE="drosophila" ;;
        11) ORGANISM_CHOICE="celegans" ;;
        12) ORGANISM_CHOICE="xenopus" ;;
        13) ORGANISM_CHOICE="chicken" ;;
        14) ORGANISM_CHOICE="rice" ;;
        15) ORGANISM_CHOICE="arabidopsis" ;;
        16) ORGANISM_CHOICE="maize" ;;
        17) ORGANISM_CHOICE="wheat" ;;
        18) ORGANISM_CHOICE="soybean" ;;
        19) ORGANISM_CHOICE="tomato" ;;
        20) ORGANISM_CHOICE="potato" ;;
        21) ORGANISM_CHOICE="barley" ;;
        22) ORGANISM_CHOICE="yeast" ;;
        23) echo "Exiting..."; exit 0 ;;
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
    rat)
        ORGANISM_NAME="rat"
        SPECIES_DESC="Rattus norvegicus (mRatBN7.2)"
        FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz"
        GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.113.gtf.gz"
        FASTA_FILE="Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz"
        GTF_FILE="Rattus_norvegicus.mRatBN7.2.113.gtf.gz"
        ;;
    pig)
        ORGANISM_NAME="pig"
        SPECIES_DESC="Sus scrofa (Sscrofa11.1)"
        FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz"
        GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.113.gtf.gz"
        FASTA_FILE="Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz"
        GTF_FILE="Sus_scrofa.Sscrofa11.1.113.gtf.gz"
        ;;
    cow)
        ORGANISM_NAME="cow"
        SPECIES_DESC="Bos taurus (ARS-UCD1.3)"
        FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.3.dna.toplevel.fa.gz"
        GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/bos_taurus/Bos_taurus.ARS-UCD1.3.113.gtf.gz"
        FASTA_FILE="Bos_taurus.ARS-UCD1.3.dna.toplevel.fa.gz"
        GTF_FILE="Bos_taurus.ARS-UCD1.3.113.gtf.gz"
        ;;
    dog)
        ORGANISM_NAME="dog"
        SPECIES_DESC="Canis lupus familiaris (ROS_Cfam_1.0)"
        FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/canis_lupus_familiaris/dna/Canis_lupus_familiaris.ROS_Cfam_1.0.dna.toplevel.fa.gz"
        GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/canis_lupus_familiaris/Canis_lupus_familiaris.ROS_Cfam_1.0.113.gtf.gz"
        FASTA_FILE="Canis_lupus_familiaris.ROS_Cfam_1.0.dna.toplevel.fa.gz"
        GTF_FILE="Canis_lupus_familiaris.ROS_Cfam_1.0.113.gtf.gz"
        ;;
    macaque)
        ORGANISM_NAME="macaque"
        SPECIES_DESC="Macaca mulatta (Mmul_10)"
        FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz"
        GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.113.gtf.gz"
        FASTA_FILE="Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz"
        GTF_FILE="Macaca_mulatta.Mmul_10.113.gtf.gz"
        ;;
    chimpanzee)
        ORGANISM_NAME="chimpanzee"
        SPECIES_DESC="Pan troglodytes (Pan_tro_3.0)"
        FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/pan_troglodytes/dna/Pan_troglodytes.Pan_tro_3.0.dna.toplevel.fa.gz"
        GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/pan_troglodytes/Pan_troglodytes.Pan_tro_3.0.113.gtf.gz"
        FASTA_FILE="Pan_troglodytes.Pan_tro_3.0.dna.toplevel.fa.gz"
        GTF_FILE="Pan_troglodytes.Pan_tro_3.0.113.gtf.gz"
        ;;
    zebrafish)
        ORGANISM_NAME="zebrafish"
        SPECIES_DESC="Danio rerio (GRCz11)"
        FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz"
        GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/danio_rerio/Danio_rerio.GRCz11.113.gtf.gz"
        FASTA_FILE="Danio_rerio.GRCz11.dna.primary_assembly.fa.gz"
        GTF_FILE="Danio_rerio.GRCz11.113.gtf.gz"
        ;;
    drosophila|fruitfly)
        ORGANISM_NAME="drosophila"
        SPECIES_DESC="Drosophila melanogaster (BDGP6.46)"
        FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz"
        GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.113.gtf.gz"
        FASTA_FILE="Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz"
        GTF_FILE="Drosophila_melanogaster.BDGP6.46.113.gtf.gz"
        ;;
    celegans|worm)
        ORGANISM_NAME="celegans"
        SPECIES_DESC="Caenorhabditis elegans (WBcel235)"
        FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz"
        GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.113.gtf.gz"
        FASTA_FILE="Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz"
        GTF_FILE="Caenorhabditis_elegans.WBcel235.113.gtf.gz"
        ;;
    xenopus|frog)
        ORGANISM_NAME="xenopus"
        SPECIES_DESC="Xenopus tropicalis (UCB_Xtro_10.0)"
        FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/xenopus_tropicalis/dna/Xenopus_tropicalis.UCB_Xtro_10.0.dna.toplevel.fa.gz"
        GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/xenopus_tropicalis/Xenopus_tropicalis.UCB_Xtro_10.0.113.gtf.gz"
        FASTA_FILE="Xenopus_tropicalis.UCB_Xtro_10.0.dna.toplevel.fa.gz"
        GTF_FILE="Xenopus_tropicalis.UCB_Xtro_10.0.113.gtf.gz"
        ;;
    chicken)
        ORGANISM_NAME="chicken"
        SPECIES_DESC="Gallus gallus (GRCg7b)"
        FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz"
        GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.113.gtf.gz"
        FASTA_FILE="Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz"
        GTF_FILE="Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.113.gtf.gz"
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
    wheat)
        ORGANISM_NAME="wheat"
        SPECIES_DESC="Triticum aestivum (IWGSC)"
        FASTA_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.toplevel.fa.gz"
        GTF_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/gtf/triticum_aestivum/Triticum_aestivum.IWGSC.60.gtf.gz"
        FASTA_FILE="Triticum_aestivum.IWGSC.dna.toplevel.fa.gz"
        GTF_FILE="Triticum_aestivum.IWGSC.60.gtf.gz"
        ;;
    soybean)
        ORGANISM_NAME="soybean"
        SPECIES_DESC="Glycine max (Glycine_max_v2.1)"
        FASTA_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/glycine_max/dna/Glycine_max.Glycine_max_v2.1.dna.toplevel.fa.gz"
        GTF_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/gtf/glycine_max/Glycine_max.Glycine_max_v2.1.60.gtf.gz"
        FASTA_FILE="Glycine_max.Glycine_max_v2.1.dna.toplevel.fa.gz"
        GTF_FILE="Glycine_max.Glycine_max_v2.1.60.gtf.gz"
        ;;
    tomato)
        ORGANISM_NAME="tomato"
        SPECIES_DESC="Solanum lycopersicum (SL3.0)"
        FASTA_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/solanum_lycopersicum/dna/Solanum_lycopersicum.SL3.0.dna.toplevel.fa.gz"
        GTF_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/gtf/solanum_lycopersicum/Solanum_lycopersicum.SL3.0.60.gtf.gz"
        FASTA_FILE="Solanum_lycopersicum.SL3.0.dna.toplevel.fa.gz"
        GTF_FILE="Solanum_lycopersicum.SL3.0.60.gtf.gz"
        ;;
    potato)
        ORGANISM_NAME="potato"
        SPECIES_DESC="Solanum tuberosum (SolTub_3.0)"
        FASTA_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/solanum_tuberosum/dna/Solanum_tuberosum.SolTub_3.0.dna.toplevel.fa.gz"
        GTF_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/gtf/solanum_tuberosum/Solanum_tuberosum.SolTub_3.0.60.gtf.gz"
        FASTA_FILE="Solanum_tuberosum.SolTub_3.0.dna.toplevel.fa.gz"
        GTF_FILE="Solanum_tuberosum.SolTub_3.0.60.gtf.gz"
        ;;
    barley)
        ORGANISM_NAME="barley"
        SPECIES_DESC="Hordeum vulgare (MorexV3)"
        FASTA_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/hordeum_vulgare/dna/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa.gz"
        GTF_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/gtf/hordeum_vulgare/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.60.gtf.gz"
        FASTA_FILE="Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa.gz"
        GTF_FILE="Hordeum_vulgare.MorexV3_pseudomolecules_assembly.60.gtf.gz"
        ;;
    yeast)
        ORGANISM_NAME="yeast"
        SPECIES_DESC="Saccharomyces cerevisiae (R64-1-1)"
        FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"
        GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.113.gtf.gz"
        FASTA_FILE="Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"
        GTF_FILE="Saccharomyces_cerevisiae.R64-1-1.113.gtf.gz"
        ;;
    *)
        echo -e "${RED}Error: Unknown organism '${ORGANISM_CHOICE}'.${RESET}"
        echo -e "Use './ref --help' to see all available options."
        exit 1
        ;;
esac

TARGET_DIR="./${ORGANISM_NAME}-ref"

echo -e "\n${BOLD}[1/4] Target Reference Folder:${RESET} ${GREEN}${TARGET_DIR}${RESET}"
echo -e "${DIM}Species: ${SPECIES_DESC}${RESET}"
echo -e "${YELLOW}${BOLD}⚠ Notice: Downloading reference genome FASTA/GTF (up to 5GB+) and building STAR index may take 15 to 30+ minutes depending on your internet connection and hardware environment.${RESET}"

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
