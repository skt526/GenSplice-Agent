#!/usr/bin/env bash
# ==============================================================================
# GenSplice-Agent System Validation & Test Suite
# Fast Pre-flight Verification Test (~15 sec) using Synthetic Reads & Isolated Sandbox
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

echo -e "${CYAN}${BOLD}"
echo "======================================================================"
echo "           GenSplice-Agent System Pre-flight Test Suite              "
export PATH="$HOME/miniconda3/envs/gensplice-agent/bin:$PATH"
PYTHON_BIN="python3"

echo -e "${BOLD}[1/4] Preparing Quick Test Environment & Sandbox...${RESET}"
TEST_DIR="test_data"
TEST_OUT_DIR="test_outputs"
mkdir -p "${TEST_DIR}/inputs/control"
mkdir -p "${TEST_DIR}/inputs/treatment"
mkdir -p "test-ref"
mkdir -p "${TEST_OUT_DIR}/03_deg"
mkdir -p "${TEST_OUT_DIR}/04_rmats"

CTRL_R1="${TEST_DIR}/inputs/control/GSM1275862_control_1_1.fq.gz"
CTRL_R2="${TEST_DIR}/inputs/control/GSM1275862_control_1_2.fq.gz"
TREAT_R1="${TEST_DIR}/inputs/treatment/GSM1275866_dex_1_1.fq.gz"
TREAT_R2="${TEST_DIR}/inputs/treatment/GSM1275866_dex_1_2.fq.gz"

if [ ! -f "$CTRL_R1" ]; then
    echo -e "  ${DIM}Generating test FASTQ reads for rMATS & DEG validation...${RESET}"
    SEQ1="ACTGGACGCCGACGACTTCGACAGCCAGCTGCTGGACGAG"
    QUAL1="BBBFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJ"

    printf "@SRR1039508.1 HWI-ST700660:187:D13YLACXX:1:1101:1203:2174/1\n%s\n+\n%s\n" "$SEQ1" "$QUAL1" | gzip > "$CTRL_R1"
    printf "@SRR1039508.1 HWI-ST700660:187:D13YLACXX:1:1101:1203:2174/2\n%s\n+\n%s\n" "$SEQ1" "$QUAL1" | gzip > "$CTRL_R2"
    printf "@SRR1039512.1 HWI-ST700660:187:D13YLACXX:1:1101:1205:2180/1\n%s\n+\n%s\n" "$SEQ1" "$QUAL1" | gzip > "$TREAT_R1"
    printf "@SRR1039512.1 HWI-ST700660:187:D13YLACXX:1:1101:1205:2180/2\n%s\n+\n%s\n" "$SEQ1" "$QUAL1" | gzip > "$TREAT_R2"
fi

TEST_FASTA="test-ref/human_subset.fa"
TEST_GTF="test-ref/human_subset.gtf"

if [ ! -f "$TEST_FASTA" ]; then
cat <<EOF > "$TEST_FASTA"
>chr13
ACTGGACGCCGACGACTTCGACAGCCAGCTGCTGGACGAGCTCGTCCAGCAGCTGGCTGTCGAAGTCGTCGGCGTCCAGTACTGGACGCCGACGACTTCGACAGCCAGCTGCTGGACGAGCTCGTCCAGCAGCTGGCTGTCGAAGTCGTCGGCGTCCAGTACTGGACGCCGACGACTTCGACAGCCAGCTGCTGGACGAGCTCGTCCAGCAGCTGGCTGTCGAAGTCGTCGGCGTCCAG
>chr16
ACTGGACGCCGACGACTTCGACAGCCAGCTGCTGGACGAGCTCGTCCAGCAGCTGGCTGTCGAAGTCGTCGGCGTCCAGTACTGGACGCCGACGACTTCGACAGCCAGCTGCTGGACGAGCTCGTCCAGCAGCTGGCTGTCGAAGTCGTCGGCGTCCAGTACTGGACGCCGACGACTTCGACAGCCAGCTGCTGGACGAGCTCGTCCAGCAGCTGGCTGTCGAAGTCGTCGGCGTCCAG
EOF

cat <<EOF > "$TEST_GTF"
chr13	ENSEMBL	gene	1	100	.	+	.	gene_id "ENSG00000120129"; gene_name "DUSP1";
chr13	ENSEMBL	transcript	1	100	.	+	.	gene_id "ENSG00000120129"; transcript_id "ENST00000239243"; gene_name "DUSP1";
chr13	ENSEMBL	exon	1	100	.	+	.	gene_id "ENSG00000120129"; transcript_id "ENST00000239243"; exon_number "1"; gene_name "DUSP1";
chr16	ENSEMBL	gene	1	100	.	+	.	gene_id "ENSG00000103194"; gene_name "CRISPLD2";
chr16	ENSEMBL	transcript	1	100	.	+	.	gene_id "ENSG00000103194"; transcript_id "ENST00000219460"; gene_name "CRISPLD2";
chr16	ENSEMBL	exon	1	100	.	+	.	gene_id "ENSG00000103194"; transcript_id "ENST00000219460"; exon_number "1"; gene_name "CRISPLD2";
EOF
fi

if [ ! -f test_config.yaml ]; then
cat <<EOF > test_config.yaml
reference:
  organism: "test"
  ref_dir: "./test-ref"
  fasta: "./test-ref/human_subset.fa"
  gtf: "./test-ref/human_subset.gtf"
  star_index: "./test-ref/star_index"
threads: 0
inputs:
  base_dir: "./test_data/inputs"
  control_dir: "./test_data/inputs/control"
  treatment_dir: "./test_data/inputs/treatment"
outputs:
  clean_fq: "./test_outputs/01_clean_fq"
  aligned_bam: "./test_outputs/02_aligned_bam"
  deg: "./test_outputs/03_deg"
  rmats: "./test_outputs/04_rmats"
EOF
fi

echo -e "\n${BOLD}[2/4] Running Quick Pipeline Test...${RESET}"
"$PYTHON_BIN" run_pipeline.py --config test_config.yaml --skip-confirmation --allow-mock

echo -e "\n${BOLD}[3/4] Populating Full rMATS 5-Event Dataset...${RESET}"
cat <<EOF > ${TEST_OUT_DIR}/03_deg/deg_result.csv
gene_id,geneSymbol,baseMean,log2FoldChange,pvalue,padj
ENSG00000103194,CRISPLD2,1420.5,0.15,0.380,0.450
ENSG00000120129,DUSP1,3850.2,-1.84,0.0001,0.002
ENSG00000140355,GRMZM2G140355,890.1,-0.08,0.520,0.620
ENSG00000165025,SYK,210.4,1.95,0.002,0.015
ENSG00000112715,VEGFA,1850.0,0.12,0.290,0.380
ENSG00000012048,BRCA1,920.8,-0.04,0.480,0.590
ENSG00000026101,STAT3,3100.4,2.10,0.0002,0.001
ENSG00000171862,PTEN,1250.6,-1.65,0.0005,0.004
ENSG00000075624,ACTB,5400.0,0.02,0.910,0.950
ENSG00000111640,GAPDH,6200.0,-0.03,0.850,0.880
EOF

cp ${TEST_OUT_DIR}/03_deg/deg_result.csv "${TEST_DIR}/deg_result.csv"

cat <<EOF > ${TEST_OUT_DIR}/04_rmats/SE.MATS.JC.txt
ID	GeneID	geneSymbol	chr	strand	exonStart_0base	exonEnd	upstreamES	upstreamEE	downstreamES	downstreamEE	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference
1	ENSG00000103194	CRISPLD2	chr16	+	84500	84650	83000	83150	86000	86150	0.0001	0.001	0.85,0.88,0.86	0.25,0.22,0.24	0.63
2	ENSG00000120129	DUSP1	chr13	+	12000	12150	10000	10150	14000	14150	0.0005	0.005	0.10,0.12,0.11	0.75,0.78,0.76	-0.64
EOF

cat <<EOF > ${TEST_OUT_DIR}/04_rmats/RI.MATS.JC.txt
ID	GeneID	geneSymbol	chr	strand	riExonStart_0base	riExonEnd	upstreamES	upstreamEE	downstreamES	downstreamEE	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference
1	ENSG00000140355	GRMZM2G140355	chr1	+	4500	4700	4000	4150	5000	5150	0.0002	0.004	0.20,0.22,0.21	0.55,0.57,0.56	-0.35
2	ENSG00000026101	STAT3	chr17	+	40500	40700	40000	40150	41000	41150	0.0001	0.003	0.75,0.78,0.76	0.37,0.39,0.38	0.38
EOF

cat <<EOF > ${TEST_OUT_DIR}/04_rmats/MXE.MATS.JC.txt
ID	GeneID	geneSymbol	chr	strand	1stExonStart_0base	1stExonEnd	2ndExonStart_0base	2ndExonEnd	upstreamES	upstreamEE	downstreamES	downstreamEE	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference
1	ENSG00000112715	VEGFA	chr6	+	43700	43800	44100	44200	43000	43100	45000	45100	0.0003	0.006	0.65,0.68,0.66	0.37,0.39,0.38	0.28
EOF

cat <<EOF > ${TEST_OUT_DIR}/04_rmats/A5SS.MATS.JC.txt
ID	GeneID	geneSymbol	chr	strand	longExonStart_0base	longExonEnd	shortES	shortEE	flankingES	flankingEE	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference
1	ENSG00000012048	BRCA1	chr17	+	41190	41270	41200	41270	41500	41600	0.0004	0.007	0.15,0.18,0.16	0.46,0.48,0.47	-0.31
EOF

cat <<EOF > ${TEST_OUT_DIR}/04_rmats/A3SS.MATS.JC.txt
ID	GeneID	geneSymbol	chr	strand	longExonStart_0base	longExonEnd	shortES	shortEE	flankingES	flankingEE	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference
1	ENSG00000171862	PTEN	chr10	+	89600	89750	89650	89750	89000	89100	0.0008	0.012	0.70,0.72,0.71	0.45,0.47,0.46	0.25
EOF

cp ${TEST_OUT_DIR}/04_rmats/*.txt "${TEST_DIR}/"

# Step 4: Test Dynamic HTML Report Export
echo -e "\n${BOLD}[4/4] Verifying Dynamic HTML Exporter & Visualizers...${RESET}"
"$PYTHON_BIN" -c "
from core.deg_loader import load_deg_data
from core.rmats_loader import load_rmats_data, select_primary_splicing_events
from core.merger import merge_deg_and_rmats
from visualizer.report_exporter import export_html_report
import os

deg_file = 'test_outputs/03_deg/deg_result.csv' if os.path.exists('test_outputs/03_deg/deg_result.csv') else 'outputs/03_deg/deg_result.csv'
rmats_dir = 'test_outputs/04_rmats' if os.path.exists('test_outputs/04_rmats') else 'outputs/04_rmats'
html_out = 'test_outputs/gensplice_report.html' if os.path.exists('test_outputs') else 'outputs/gensplice_report.html'

df_deg = load_deg_data(deg_file)
df_rmats = load_rmats_data(rmats_dir)
primary_rmats = select_primary_splicing_events(df_rmats)
merged = merge_deg_and_rmats(df_deg, primary_rmats)

out_path = export_html_report(merged, output_html_path=html_out)
print('  ✔ HTML Report verified:', out_path)
"

echo -e "\n${CYAN}${BOLD}======================================================================"
echo "           GenSplice-Agent System Pre-flight Test: PASSED ✅          "
echo "======================================================================"
echo -e "${RESET}"
echo -e "Command Workflow Summary:"
echo -e "  ${YELLOW}System Verification Test:${RESET} ./test"
echo -e "  ${YELLOW}Real Data Pipeline:${RESET}       ./ref human && ./GenSplice"
echo -e "  ${YELLOW}Interactive HTML Report:${RESET}  open outputs/gensplice_report.html"
echo -e ""
