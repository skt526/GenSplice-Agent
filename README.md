# GenSplice-Agent 🧬

> **Integrated Transcriptomic Profiling & Alternative Splicing Visualization Platform with Standalone Interactive Reporting**

`GenSplice-Agent` is a unified computational transcriptomics platform that integrates quantitative gene expression changes (DEG, PyDESeq2) and qualitative isoform variations (Alternative Splicing, rMATS) into a single 4-quadrant coordinate space. It executes locally with **zero cloud/server costs** and generates a fully self-contained, interactive HTML report (`gensplice_report.html`).

---

## 🚦 Command Workflow

### 🧪 1. System Verification & Testing
Verify that the end-to-end pipeline and standalone HTML report generation execute properly on the reference test dataset (GSE52778 Dexamethasone model):
```bash
bash install.sh && ./test
```
*(Or `./install` followed by `./test`)*

### 🧬 2. Production Sample Analysis
Configure the desired reference genome (Human, Mouse, Rice, etc.) and launch the analysis pipeline on your experimental samples:
```bash
bash install.sh && ./ref human && ./GenSplice
```
*(Or `./install` followed by `./ref human` followed by `./GenSplice`)*

---

## 📦 Required Dependencies (Conda/Bioconda Ecosystem)

| Category | Tool / Package | Recommended Version | Primary Role & Description |
| :--- | :--- | :--- | :--- |
| **CLI Tools** | `fastp` | `>= 0.23.4` | Raw FASTQ read quality control (QC) and automated adapter trimming |
| | `STAR` | `>= 2.7.11a` | High-performance splice-aware genome alignment engine (2-pass mode) |
| | `rmats` | `>= 4.3.0` | Differential alternative splicing analysis (SE, RI, A5SS, A3SS, MXE) |
| | `subread` | `>= 2.0.6` | `featureCounts` transcriptomic read quantification engine |
| | `R` (r-base) | `>= 4.2` | Statistical verification backend for rMATS |
| **Data & Visualizer** | `polars` | `>= 0.20.0` | High-throughput columnar table merging and fast 4-quadrant parser |
| | `plotly` | `>= 5.20.0` | Interactive 4-quadrant cross-plot & Sashimi exon structure visualizer |
| | `pyyaml` | `>= 6.0` | Pipeline configuration loader |

---

## 🧪 Benchmark Dataset Specification (GSE52778 Airway Smooth Muscle)

The test dataset bundled with `./test` is a widely recognized academic standard benchmark:
- **GEO Accession**: GSE52778 (SRA: SRP033346) (Himes et al., 2014)
- **Experimental Model**: Human Airway Smooth Muscle Cells (Control vs. Dexamethasone-treated)
- **Target Regulators**: Validated automated detection of differential expression and 5 major alternative splicing classes across immune-modulatory targets (e.g., *DUSP1*, *CRISPLD2*).

---

## 📁 Directory Architecture

```text
GenSplice-Agent/
├── install.sh              # One-click environment installer (Conda/Bioconda)
├── test.sh                 # End-to-end verification script using GSE52778 benchmark
├── ref.sh                  # One-click reference genome & annotation downloader
├── GenSplice               # Main pipeline execution wrapper script
├── run_pipeline.py         # Pipeline orchestrator and workflow manager (Python)
├── config.yaml             # Pipeline configuration (threads, paths, cutoffs)
├── environment.yml         # Conda environment specification
├── inputs/                 # [User FASTQ Input Directory]
│   ├── control/            # Control replicate FASTQ files (.fastq / .fq.gz)
│   └── treatment/          # Treatment replicate FASTQ files (.fastq / .fq.gz)
├── {organism}-ref/         # Downloaded reference genome (FASTA, GTF) and STAR index
└── outputs/                # Analysis outputs and standalone HTML report
    ├── 01_clean_fq/        # fastp trimmed and filtered FASTQ files
    ├── 02_aligned_bam/     # Coordinate-sorted BAM alignment files
    ├── 03_deg/             # PyDESeq2 / featureCounts gene-level quantification
    ├── 04_rmats/           # rMATS 5-event alternative splicing matrices
    └── gensplice_report.html # [Final] Standalone interactive HTML report dashboard
```

---

## 📥 Input Data Specification (FASTQ)

To analyze raw RNA-seq data, place paired-end or single-end sequencing files into the `inputs/` subdirectories:

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
- **Compressed (Recommended)**: `.fq.gz`, `.fastq.gz`
- **Uncompressed**: `.fq`, `.fastq`

### 2. Paired-End File Naming Conventions
Sample pairs are automatically recognized using standard read pair suffixes:
- **Pattern 1**: `*_1.fq.gz` / `*_2.fq.gz` or `*_1.fastq.gz` / `*_2.fastq.gz`
- **Pattern 2**: `*_R1.fastq.gz` / `*_R2.fastq.gz` or `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz`

### 3. Single-End Support
- Single files (e.g., `sampleA.fq.gz`) are automatically detected as single-end reads, and QC/alignment stages adjust parameters accordingly.

---

## 📄 Documentation
- [BLUEPRINT.md](BLUEPRINT.md): Architecture blueprint, 4-quadrant mathematical model, and methodology.
- [environment.yml](environment.yml): Conda environment definition.
