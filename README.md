# GenSplice-Agent 🧬

> **Automated RNA-seq Alternative Splicing Analysis & Interactive Reporting Pipeline**

`GenSplice-Agent` is an automated computational pipeline for differential alternative splicing analysis from raw RNA-seq data. It orchestrates quality control (`fastp`), splice-aware genome alignment (`STAR`), and event-level alternative splicing quantification (`rMATS`), producing a self-contained interactive HTML report (`gensplice_report.html`). The report provides an Alternative Splicing Volcano Plot, dynamic significance threshold filtering, exon-intron Sashimi plots, isoform functional annotations, and RT-qPCR primer designs.

---

## 🚦 Execution Workflow

### 1. Installation & Environment Setup
Install required bioinformatics tools and Python dependencies via Conda/Bioconda:
```bash
./install
```

### 2. Reference Genome Preparation
Download and index reference genome files (FASTA, GTF, STAR index). Supported organisms are **Human** (*Homo sapiens*, GRCh38) and **Mouse** (*Mus musculus*, GRCm39):
```bash
./ref
```
*(Select `1` for Human or `2` for Mouse in the interactive terminal menu)*

### 3. Pipeline Execution
Place raw FASTQ sequencing files in `inputs/control` and `inputs/treatment`, then execute the pipeline:
```bash
./GenSplice
```

---

## 📦 Required Dependencies (Conda/Bioconda)

| Category | Tool / Package | Recommended Version | Description |
| :--- | :--- | :--- | :--- |
| **Upstream Tools** | `fastp` | `>= 0.23.4` | Read quality control and adapter trimming |
| | `STAR` | `>= 2.7.11a` | Splice-aware genome alignment (2-pass mode) |
| | `rmats` | `>= 4.3.0` | Alternative splicing event quantification (SE, RI, A5SS, A3SS, MXE) |
| | `subread` | `>= 2.0.6` | `featureCounts` read summarization |
| **Data Processing & Visualization** | `polars` | `>= 0.20.0` | Columnar table processing |
| | `plotly` | `>= 5.20.0` | Alternative Splicing Volcano Plot & Sashimi exon visualizer |
| | `pyyaml` | `>= 6.0` | Configuration parser |

---

## 📁 Directory Architecture

```text
GenSplice-Agent/
├── install.sh              # Environment installer (Conda/Bioconda)
├── ref.sh                  # Reference genome & annotation downloader (human/mouse)
├── GenSplice               # Pipeline execution wrapper script
├── run_pipeline.py         # Pipeline orchestrator (Python)
├── config.yaml             # Pipeline configuration
├── environment.yml         # Conda environment definition
├── inputs/                 # Input directory for FASTQ files
│   ├── control/            # Control replicate FASTQ files (.fastq / .fq.gz)
│   └── treatment/          # Treatment replicate FASTQ files (.fastq / .fq.gz)
├── {organism}-ref/         # Reference genome (FASTA, GTF) and STAR index
└── outputs/                # Analysis outputs and report
    ├── 01_clean_fq/        # fastp trimmed FASTQ files
    ├── 02_aligned_bam/     # Coordinate-sorted BAM alignment files
    ├── 03_deg/             # Gene-level quantification matrices
    ├── 04_rmats/           # rMATS alternative splicing event outputs
    └── gensplice_report.html # Standalone interactive HTML report dashboard
```

---

## 📥 Input Data Specification (FASTQ)

Place paired-end or single-end sequencing files into the `inputs/` subdirectories:

```text
inputs/
├── control/
│   ├── sampleA_1.fq.gz
│   └── sampleA_2.fq.gz
└── treatment/
    ├── sampleB_1.fq.gz
    └── sampleB_2.fq.gz
```

### 1. Supported File Extensions
- Compressed: `.fq.gz`, `.fastq.gz`
- Uncompressed: `.fq`, `.fastq`

### 2. Paired-End Naming Conventions
- Pattern 1: `*_1.fq.gz` / `*_2.fq.gz` or `*_1.fastq.gz` / `*_2.fastq.gz`
- Pattern 2: `*_R1.fastq.gz` / `*_R2.fastq.gz` or `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz`

### 3. Single-End Reads
Single files are automatically detected as single-end reads, and parameters are adjusted accordingly.

---

## 📄 Documentation
- [BLUEPRINT.md](BLUEPRINT.md): Architecture, methodology, and computational framework.
- [environment.yml](environment.yml): Conda environment definition.
