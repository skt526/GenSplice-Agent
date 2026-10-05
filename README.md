# GenSplice-Agent 🧬

> **Automated RNA-seq Alternative Splicing Analysis, Isoform Functional Annotation & Interactive Reporting Pipeline**

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![rMATS 4.3+](https://img.shields.io/badge/rMATS-4.3.0+-green.svg)](https://rnaseq-mats.sourceforge.net/)
[![STAR 2.7+](https://img.shields.io/badge/STAR-2.7.11a+-orange.svg)](https://github.com/alexdobin/STAR)

`GenSplice-Agent` is a fully automated, end-to-end computational pipeline for differential alternative splicing quantification and functional validation from raw RNA-seq data. It streamlines quality control (`fastp`), splice-aware genome alignment (`STAR`), gene-level differential expression (`featureCounts` / DESeq2 / CPM), and event-level alternative splicing quantification (`rMATS`). 

The pipeline produces a self-contained, interactive Light Mode HTML dashboard (`gensplice_report.html`) that synchronizes an Alternative Splicing Volcano Plot, dynamic threshold sliders, exon-intron Sashimi plots, sample-level replicate profiles, transcript functional impairment evaluation, live NCBI/PubMed literature querying, and RT-qPCR primer design.

---

## 🌟 Key Capabilities

- **Automated Upstream Orchestration**: Single-command execution from raw FASTQ files to an interactive HTML report, with automatic read-type (single-end / paired-end) and average read length detection.
- **Alternative Splicing Volcano Framework**: Real-time 2D Volcano Plot ($\Delta\text{PSI}$ vs. $-\log_{10}\text{FDR}$) with shaded significance zones (Inclusion Favored, Exclusion Favored, Non-Significant Background).
- **Dynamic Threshold Control**: Client-side reactive sliders for $\Delta\text{PSI}$ effect size and rMATS FDR significance, instantly updating event counts, candidate gene sets, and downstream visualizations.
- **Sample-Level Replicate PSI & Junction Depth Profile**: Decomposes biological replicate concordance, mean PSI values, and stacked raw Inclusion Junction Counts (IJC) vs. Skipping Junction Counts (SJC).
- **Heuristic Functional Impairment Priority Index (0–100)**: Evaluates CDS reading frame shifts ($\Delta L \pmod 3 \ne 0$), de novo premature termination codons (PTC), and canonical 50–55 nt Nonsense-Mediated mRNA Decay (NMD) rules.
- **Isoform-Specific RT-qPCR Primer Designer**: Thermodynamic nearest-neighbor design (Primer3 engine) spanning exon-exon junctions (Inclusion vs. Exclusion) using authentic reference FASTA sequences.
- **Synchronized Multi-Card Navigation**: Selecting a gene in the Volcano plot, table, or dropdown instantly updates the Sashimi plot, Replicate profile, Isoform matrix, Primer table, and live NCBI/PubMed literature.
- **Fault-Tolerant Checkpointing**: Granular stage tracking (`pipeline_checkpoint.json`) allows interrupted runs to resume seamlessly. Previous outputs are automatically archived with timestamps (`outputs_{YYMMDD}_{HHMMSS}/`).

---

## 🚦 Quickstart Workflow

### 1. Installation & Environment Setup
Clone the repository and install all required bioinformatics tools and Python dependencies via Conda:
```bash
git clone https://github.com/skt526/GenSplice-Agent.git
cd GenSplice-Agent
./install
```

### 2. Reference Genome Preparation
Download and index Ensembl reference genomes (FASTA, GTF, and STAR 2-pass index). Both **Human** (*Homo sapiens*, GRCh38) and **Mouse** (*Mus musculus*, GRCm39) are supported:
```bash
./ref
```
*(Select `1` for Human or `2` for Mouse in the interactive terminal menu, or pass arguments directly: `./ref human`)*

### 3. Place Input FASTQ Files
Organize your raw sequencing files into the `inputs/control` and `inputs/treatment` directories:
```text
inputs/
├── control/
│   ├── control_rep1_1.fq.gz
│   └── control_rep1_2.fq.gz
└── treatment/
    ├── treatment_rep1_1.fq.gz
    └── treatment_rep1_2.fq.gz
```

> [!TIP]
> **Need test data?** You can automatically download the genuine human CLL splicing benchmark dataset (**GSE190087** / SF3B1 K700E):
> ```bash
> ./download quick       # Fast 1 vs 1 pair test (~6 GB)
> # or
> ./download benchmark   # Full 3 vs 3 replicate benchmark (~17 GB)
> ```

### 4. Execute the Pipeline
Run the main pipeline wrapper:
```bash
./GenSplice
```
*(Alternatively, execute directly via Python: `python run_pipeline.py`)*

Once complete, open `outputs/gensplice_report.html` in any modern web browser.

---

## 📦 Software Stack & Dependencies

All dependencies are automatically managed via the `gensplice-agent` Conda environment:

| Category | Tool / Package | Recommended Version | Purpose |
| :--- | :--- | :--- | :--- |
| **Upstream QC** | `fastp` | `>= 0.23.4` | Adapter trimming, poly-G filtering, QC reporting |
| **Alignment** | `STAR` | `>= 2.7.11a` | 2-pass splice-aware reference genome mapping |
| **Quantification** | `subread` (`featureCounts`) | `>= 2.0.6` | Gene-level read count summarization |
| **Splicing Engine** | `rMATS-turbo` (`rmats.py`) | `>= 4.3.0` | Splicing event detection (SE, RI, MXE, A5SS, A3SS) |
| **Primer Engine** | `primer3-py` / `pyfaidx` | `>= 2.0.0` | Thermodynamic primer design & genomic FASTA extraction |
| **Data Processing** | `polars` / `pandas` | `>= 0.20.0` | High-throughput columnar table processing & filtering |
| **Interactive UI** | `plotly` / `jinja2` | `>= 5.20.0` | Client-side reactive HTML visualization |

---

## 📁 Directory Architecture

```text
GenSplice-Agent/
├── install                 # Shell installer (Conda environment builder)
├── install.sh              # Installer source script
├── ref                     # Reference genome setup wrapper
├── ref.sh                  # Downloader & STAR index generator (human/mouse)
├── download                # Benchmark dataset downloader wrapper
├── download_GSE190087.sh   # ENA direct download script (GSE190087 / PRJNA786720)
├── download_GSE190087.py   # Python dataset downloader engine
├── GenSplice               # Main pipeline execution wrapper
├── GenSplice.sh            # Pipeline execution script
├── run_pipeline.py         # 5-step pipeline orchestrator
├── config.yaml             # User-configurable pipeline settings
├── config.py               # Theme colors, cutoffs, and default parameters
├── environment.yml         # Conda environment specification
│
├── core/                   # Core computational modules
│   ├── fastq_detector.py   # Automatic sample discovery and pair matching
│   ├── fastp_parser.py     # Average read length auto-detection from QC JSON
│   ├── deg_calculator.py   # Gene expression quantification & FDR calculation
│   ├── deg_loader.py       # DEG table loader & Ensembl-to-symbol mapping
│   ├── rmats_loader.py     # Multi-replicate rMATS event parser & noise filter
│   ├── merger.py           # 4-quadrant merger & splicing status classification
│   ├── isoform_annotator.py# CDS frame shift, PTC projection, & NMD evaluator
│   ├── primer_designer.py  # Isoform-specific RT-qPCR primer generator (Primer3)
│   └── enrichment.py       # Genuine GO & KEGG enrichment query engine
│
├── visualizer/             # Interactive plotting & report generator
│   ├── as_volcano_plot.py  # Alternative Splicing Volcano Plot builder
│   ├── exon_structure.py   # 2D Sashimi-style exon-intron structure renderer
│   ├── enrichment_plot.py  # Horizontal bar chart enrichment visualizer
│   └── report_exporter.py  # Standalone Light Mode HTML dashboard exporter
│
├── inputs/                 # Input FASTQ directories
│   ├── control/            # Control group FASTQ files (.fq.gz, .fastq.gz)
│   └── treatment/          # Treatment group FASTQ files (.fq.gz, .fastq.gz)
│
├── {organism}-ref/         # Reference genome (FASTA, GTF) & STAR index
│
└── outputs/                # Analysis outputs & interactive report
    ├── 01_clean_fq/        # fastp trimmed FASTQ files & HTML/JSON QC reports
    ├── 02_aligned_bam/     # STAR coordinate-sorted BAM files
    ├── 03_deg/             # featureCounts matrix & full deg_result.csv
    ├── 04_rmats/           # rMATS JC/JCEC output files for all 5 event types
    ├── pipeline_checkpoint.json # Granular stage execution checkpoints
    ├── pipeline_status.log # Timestamped operational logs
    └── gensplice_report.html # Standalone interactive report dashboard
```

---

## 📥 Input Data Specifications

Place your sequencing files into `inputs/control` and `inputs/treatment`. Supported naming patterns:

### 1. Paired-End Reads (Auto-Paired)
- Illumina standard: `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz`
- Shorthand standard: `*_R1.fastq.gz` / `*_R2.fastq.gz`
- SRA / EBI standard: `*_1.fq.gz` / `*_2.fq.gz` or `*_1.fastq.gz` / `*_2.fastq.gz`

### 2. Single-End Reads
- Any single file (e.g. `sampleA.fastq.gz` or `SRR12345.fq.gz`) is automatically identified as single-end. The pipeline sets STAR and rMATS parameters (`-t single`) accordingly.

### 3. File Integrity Verification
Before executing QC or alignment, the orchestrator validates gzip compression integrity (`gzip -t`) to catch incomplete downloads or truncated files early.

---

## ⚙️ Configuration (`config.yaml`)

Pipeline parameters can be customized via `config.yaml`:

```yaml
system:
  total_ram_gb: 16                  # Total available system memory (GB)
  star_bam_sort_ram_gb: 12          # Memory allocated for STAR BAM sorting (GB)
  assigned_threads: 0               # Thread count (0 = auto-allocate 80% CPU cores)

reference:
  organism: "human"                 # "human" (GRCh38) or "mouse" (GRCm39)
  ref_dir: "./human-ref"
  fasta: "./human-ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
  gtf: "./human-ref/Homo_sapiens.GRCh38.113.gtf"
  star_index: "./human-ref/star_index"

splicing:
  min_junction_reads: 10            # Minimum total junction reads across replicates
  delta_psi_cutoff: 0.10            # Default |ΔPSI| effect size threshold
  fdr_cutoff: 0.05                  # Default rMATS FDR significance threshold

inputs:
  base_dir: "./inputs"
  control_dir: "./inputs/control"
  treatment_dir: "./inputs/treatment"
  lib_type: "fr-unstranded"         # "fr-unstranded", "fr-firststrand", or "fr-secondstrand"
```

---

## 📊 Interactive HTML Report Dashboard

The generated report (`outputs/gensplice_report.html`) is completely standalone (no web server required) and provides six synchronized analysis cards:

1. **Alternative Splicing Volcano Plot & Control Panel**:
   - Interactive Plotly volcano plot displaying $\Delta\text{PSI}$ vs. $-\log_{10}\text{FDR}$.
   - Live sliders for $\Delta\text{PSI}$ and FDR recalculate event classifications and KPI metrics on the fly.
   - One-click export of filtered splicing target tables to CSV.
2. **Visual Exon-Intron Sashimi Plot**:
   - Displays event-specific genomic architecture (SE, RI, MXE, A5SS, A3SS).
   - Shows inclusion/exclusion junction arcs annotated with raw junction read counts.
   - Multi-event selector allows switching between multiple splicing events on the same gene.
3. **Sample-Level Replicate PSI & Junction Depth Profile**:
   - Replicate-level PSI bar charts and stacked IJC vs. SJC read depth charts.
   - Biological replicate concordance metric (e.g. `100% Concordant`).
   - Detailed sample table showing individual replicate values and deviations from group means.
4. **Event-Level Isoform Annotation Matrix**:
   - Classifies transcript alterations: CDS reading frame shifts vs. in-frame variations.
   - Predicts premature termination codons (PTC) and evaluates canonical 50–55 nt NMD degradation vulnerability.
   - Literature-grounded **Heuristic Functional Impairment Priority Index (0–100)** with High, Moderate, and Low priority tiers.
5. **Synchronized NCBI Gene Details & Live NIH PubMed Explorer**:
   - Live querying of NCBI Gene metadata (official gene name, chromosomal location, functional summary).
   - Real-time PubMed literature search via NIH E-utilities API displaying recent peer-reviewed publications.
6. **Isoform-Specific RT-qPCR Primer Designer**:
   - Primer3-calculated primer pairs specifically targeting Inclusion vs. Exclusion isoform junctions.
   - Reports forward/reverse sequences ($5' \rightarrow 3'$), melting temperatures ($T_m$), GC content, and amplicon lengths.

---

## 🛡️ Reliability & Data Integrity

- **Strict Empirical Data**: Synthetic or mock data generation has been removed from all modules. All metrics, p-values, junction counts, and gene models derive strictly from the user's input data and official reference annotations.
- **glibc Memory Protection**: rMATS post-stage includes an automated single-thread retry mechanism to prevent C-heap memory corruption on high-core Linux systems.
- **Polars Order Determinism**: Splicing event prioritization enforces `maintain_order=True` across multithreaded operations to ensure fully reproducible target ranking.
- **Safe Coordinate Parsing**: Robust coordinate parsing protects downstream isoform consequence evaluation from unformatted or novel junction annotations.

---

## 📄 License & Citation

GenSplice-Agent is released under the **MIT License**.

If you use GenSplice-Agent in your research, please cite:
- **STAR**: Dobin et al., *Bioinformatics* (2013). doi:10.1093/bioinformatics/bts635
- **rMATS**: Shen et al., *PNAS* (2014). doi:10.1073/pnas.1413973111
- **fastp**: Chen et al., *Bioinformatics* (2018). doi:10.1093/bioinformatics/bty560
- **Primer3**: Untergasser et al., *Nucleic Acids Res.* (2012). doi:10.1093/nar/gks596
