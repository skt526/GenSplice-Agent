# GenSplice-Agent Architecture Blueprint

GenSplice-Agent is an automated computational pipeline for profiling and visualizing alternative splicing events from raw RNA-seq sequencing data. It quantifies five major classes of alternative splicing using rMATS, provides event-level functional consequence annotations, and generates a standalone interactive HTML dashboard.

---

### 1. Biological Context & Alternative Splicing Volcano Framework

Alternative splicing is a primary mechanism diversifying the eukaryotic transcriptome. Changes in splice site selection can alter protein functional domains, shift reading frames, or trigger nonsense-mediated mRNA decay (NMD) without necessarily altering overall gene-level transcriptional abundance.

GenSplice-Agent visualizes differential alternative splicing events using an Alternative Splicing Volcano Plot:

```text
                  ▲ −log10(rMATS FDR) [Significance]
                  │
   Exclusion      │      Inclusion
   Favored        │      Favored
 (ΔPSI ≤ −cutoff) │    (ΔPSI ≥ +cutoff)
                  │
 ─────────────────┼─────────────────▶ ΔPSI (Exon Inclusion Difference)
          Non-Significant Background
```

- **Inclusion Favored**: Events with $\Delta\text{PSI} \ge +\text{cutoff}$ and $\text{FDR} \le \text{cutoff}$ (increased exon/intron inclusion in treatment compared to control).
- **Exclusion Favored**: Events with $\Delta\text{PSI} \le -\text{cutoff}$ and $\text{FDR} \le \text{cutoff}$ (increased exon skipping in treatment compared to control).
- **Non-Significant Background**: Splicing events below the user-defined $\Delta\text{PSI}$ or FDR cutoffs.

The interactive dashboard allows users to adjust both $\Delta\text{PSI}$ and FDR thresholds dynamically, immediately updating event counts, target gene lists, Sashimi exon structures, and downstream primer designs.

---

### 2. Methodological & Computational Modules

GenSplice-Agent comprises four primary computational modules:

#### 🧬 Module A: Automated Upstream Processing & Splicing Quantification
- **Read QC and Trimming (`fastp`)**: Performs automated adapter detection, quality filtering, and read length calculation.
- **Genome Alignment (`STAR`)**: Executes 2-pass splice-aware alignment using Ensembl reference genomes.
- **Alternative Splicing Quantification (`rMATS`)**: Quantifies junction and exon reads across five canonical event classes:
  - Skipped Exon (SE)
  - Retained Intron (RI)
  - Alternative 5' Splice Site (A5SS)
  - Alternative 3' Splice Site (A3SS)
  - Mutually Exclusive Exons (MXE)

#### 🔬 Module B: Isoform Functional Impairment Evaluation (`core/isoform_annotator.py`)
- Evaluates putative functional consequences of detected splicing shifts:
  - **CDS Reading Frame Assessment**: Evaluates whether spliced segment lengths maintain coding triplet frames ($\Delta L \pmod 3 = 0$) or introduce frame shifts.
  - **Premature Termination Codon (PTC) Detection**: Predicts candidate stop codon positions relative to the terminal exon junction.
  - **Nonsense-Mediated mRNA Decay (NMD) Heuristic**: Flags transcripts sensitive to degradation based on the 50–55 nt rule upstream of the last exon-exon junction.
  - **Functional Risk Stratification**: Ranks events into High, Moderate, and Low risk tiers.

#### 🧪 Module C: Isoform-Specific RT-qPCR Primer Design (`core/primer_designer.py`)
- Generates discriminative primer pairs targeting inclusion-specific and exclusion-specific splice junctions.
- Reports primer sequences, melting temperatures ($T_m$), GC content, amplicon lengths, and specificity metrics to assist experimental validation.

#### 📊 Module D: Standalone Interactive HTML Dashboard (`visualizer/report_exporter.py`)
- Generates a self-contained HTML report (`outputs/gensplice_report.html`) without server or database dependencies.
- Interactive features include:
  - Alternative Splicing Volcano Plot with dynamic threshold recalculation in pure JavaScript / Plotly.
  - Interactive Sashimi-style exon-intron structure visualizer responsive to gene selection.
  - Synchronized NCBI Gene summary and PubMed literature query interface.
  - Gene Ontology (GO Biological Process) and KEGG pathway enrichment visualizer.

---

### 3. Pipeline Software Stack & Compatibility

| Component | Software / Library | Version | Role |
| :--- | :--- | :--- | :--- |
| Read QC & Trimming | `fastp` | `>= 0.23.4` | Raw read quality filtering and trimming |
| Read Alignment | `STAR` | `>= 2.7.11a` | Splice-aware reference genome alignment |
| Gene Summarization | `featureCounts` (subread) | `>= 2.0.6` | Exon- and gene-level read summarization |
| Splicing Detection | `rMATS` | `>= 4.3.0` | Junction read counting and statistical testing |
| Table Processing | `Polars` / `Pandas` | `>= 0.20.0` | Columnar data operations and merging |
| Visualization | `Plotly.js` | `>= 5.20.0` | Client-side interactive plotting |

**Supported Organisms**:
- **Human** (*Homo sapiens*): Ensembl GRCh38
- **Mouse** (*Mus musculus*): Ensembl GRCm39

---

### 4. Technical Characteristics

1. **Local Execution**: Runs entirely in user environments without external cloud data transmission.
2. **Deterministic Checkpointing**: Tracks stage completion in `outputs/pipeline_checkpoint.json` to allow resumption without re-running completed steps.
3. **Hardware Resource Management**: Automatically allocates available CPU cores and memory limits to balance throughput and stability.
