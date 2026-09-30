# GenSplice-Agent Architecture Blueprint

GenSplice-Agent is an integrated computational biology platform developed to systematically decouple and cross-evaluate quantitative transcriptional abundance changes ($\text{Log}_2\text{FC}$, DEG) against qualitative post-transcriptional alternative splicing variations ($\Delta\text{PSI}$, rMATS) within a unified coordinate system.

---

### 1. Biological Rationale & The 4-Quadrant Framework

Conventional RNA-seq pipelines disproportionately focus on total transcript abundance (Differential Expression Analysis), systematically missing critical regulatory switches where total mRNA remains stable but alternative exon inclusion/skipping introduces functional domain shifts or triggers Nonsense-Mediated Decay (NMD).

GenSplice-Agent combines PyDESeq2 expression estimates with rMATS splicing indices to categorize the transcriptome into four distinct biological quadrants:

```text
                  ▲ ΔPSI (Differential Splicing Inclusion Index)
                  │
     [Q2] Splicing-Driven Only   │   [Q1] Dual Responders
   (Abundance Invariant, AS Switch)│ (Both Abundance and Isoform Altered)
   ★ Hidden Master Regulators ★   │
──────────────────┼──────────────────▶ Log2 Fold Change (Expression Abundance)
     [Q3] Invariant Background   │   [Q4] Expression-Driven Only
       (Homeostatic Baseline)     │  (Quantity Shift, Isoform Unchanged)
                  │
```

- **Quadrant 1 (Q1 - Dual Responders)**: Transcripts exhibiting significant changes in both transcriptional abundance ($|\text{Log}_2\text{FC}| \ge \text{cutoff}$) and alternative splicing ($|\Delta\text{PSI}| \ge \text{cutoff}$). Splicing in Q1 often acts synergistically to amplify or modulate pathway output.
- **Quadrant 2 (Q2 - Splicing-Driven Regulators)**: Transcripts undergoing significant alternative splicing switches without significant alterations in overall expression level. These represent key post-transcriptional regulators invisible to standard DEG screens.
- **Quadrant 3 (Q3 - Invariant Baseline)**: Transcripts maintaining homeostatic stability across both transcriptional and splicing axes.
- **Quadrant 4 (Q4 - Expression-Driven Only)**: Classical differential expression targets altering overall cellular abundance while preserving constitutive transcript architecture.

---

### 2. Methodological & Computational Modules

GenSplice-Agent incorporates four primary deterministic computational modules:

#### 🧬 Module A: Paired-Sample DEG Quantification Engine (`core/deg_calculator.py`)
- Executes Negative Binomial generalized linear model testing via PyDESeq2.
- Automatically handles paired-donor experimental designs to absorb patient/replicate baseline variance, enhancing true biological signal detection.

#### 🔬 Module B: Event-Level Isoform Annotation & Functional Impairment Engine (`core/isoform_annotator.py`)
- Evaluates alternative splicing events across all five canonical types: Skipped Exon (SE), Retained Intron (RI), Alternative 5' Splice Site (A5SS), Alternative 3' Splice Site (A3SS), and Mutually Exclusive Exons (MXE).
- Assesses functional impairment risk based on:
  - **CDS Reading Frame Alteration**: Detects non-triplet indels ($\Delta L \pmod 3 \neq 0$).
  - **Premature Termination Codon (PTC) Formation**: Identifies stop codon introduction upstream of the terminal exon.
  - **Nonsense-Mediated mRNA Decay (NMD)**: Flags degradation-sensitive transcripts according to canonical boundary rules.
  - **Quantitative Loss-of-Function (LoF) Score (0–100%)**: Ranks targets into High, Moderate, and Low Functional Risk tiers.

#### 🧪 Module C: Isoform-Specific RT-qPCR Primer Designer (`core/primer_designer.py`)
- Designs discriminative primer pairs targeting inclusion-specific vs. exclusion-specific splice junctions.
- Automatically computes primer melting temperatures ($T_m$), GC content percentage, amplicon sizes, and quality scores to streamline wet-lab bench validation.

#### 📊 Module D: Self-Contained Interactive HTML Dashboard (`visualizer/report_exporter.py`)
- Eliminates cloud hosting and recurring server dependencies by generating a single lightweight HTML dashboard (`outputs/gensplice_report.html`).
- Features client-side recalculation engines in pure JavaScript/Plotly:
  - Dynamic Log2FC and ΔPSI threshold sliders with real-time KPI re-tallying.
  - Standardized horizontal bar charts for GO Biological Process and KEGG Pathway enrichment.
  - Interactive Exon-Intron Sashimi track visualizer responsive to gene selection.
  - Live NIH NCBI & PubMed REST API literature retrieval via CORS-enabled E-utilities.

---

### 3. Pipeline Software Stack & Compatibility

| Component | Software / Library | Role |
| :--- | :--- | :--- |
| Read QC & Trimming | `fastp` | High-throughput quality control and automated adapter clipping |
| Read Alignment | `STAR` | 2-pass splice-aware genomic alignment |
| Gene Quantification | `featureCounts` (subread) | Exon-level and gene-level read summarization |
| Splicing Detection | `rMATS` | Statistical junction count and PSI calculation |
| Expression Testing | `PyDESeq2` | Differential gene expression modeling |
| Columnar Processing | `Polars` | Sub-millisecond table joins and multi-criteria filtering |
| Visualization | `Plotly.js` | Browser-side interactive cross-plotting and sashimi tracks |
| Reference Metadata | NIH NCBI E-utilities | Direct programmatic gene summary and literature access |

---

### 4. Publication Strategy & Target Venues

GenSplice-Agent is designed as a standalone, reproducible software tool suited for peer-reviewed computational biology and bioinformatics journals:

- **Target Journals**:
  - *Bioinformatics* (Application Notes)
  - *Briefings in Bioinformatics*
  - *BMC Bioinformatics* / *PLOS Computational Biology*
  - *Frontiers in Genetics* / *Frontiers in Bioinformatics*
  - *Scientific Reports*

- **Key Reviewer Selling Points**:
  1. **100% Deterministic Reproducibility**: Pure algorithmic execution without non-deterministic external dependencies.
  2. **End-to-End Automation**: Ingests raw FASTQ files and generates bench-ready RT-qPCR primer pairs and publication-quality figures without requiring multi-language orchestration.
  3. **Zero-Cost Deployment**: Requires no persistent web servers, databases, or third-party paid subscriptions.
