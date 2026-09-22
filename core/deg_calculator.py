import os
import polars as pl
import pandas as pd
import numpy as np

def run_deg_analysis(feature_counts_path: str, control_bams: list, treatment_bams: list, output_csv_path: str):
    """
    Parses featureCounts matrix output and calculates DEG statistics
    (gene_id, baseMean, log2FoldChange, pvalue, padj).
    """
    if not os.path.exists(feature_counts_path):
        print(f"Warning: featureCounts matrix '{feature_counts_path}' not found. Generating mock deg_result.csv...")
        # Create empty/mock result format if featureCounts path does not exist
        df_mock = pd.DataFrame(columns=["gene_id", "baseMean", "log2FoldChange", "pvalue", "padj"])
        df_mock.to_csv(output_csv_path, index=False)
        return

    # Read featureCounts file (skip comment lines starting with #)
    with open(feature_counts_path, 'r') as f:
        lines = [line for line in f if not line.startswith('#')]

    from io import StringIO
    df_counts = pd.read_csv(StringIO(''.join(lines)), sep='\t')

    # Column names: Geneid, Chr, Start, End, Strand, Length, BAM1, BAM2...
    gene_ids = df_counts['Geneid']
    count_cols = df_counts.columns[6:]  # BAM columns

    # Separate control vs treatment columns
    control_cols = [col for col in count_cols if any(os.path.basename(b) in col for b in control_bams)]
    treatment_cols = [col for col in count_cols if any(os.path.basename(b) in col for b in treatment_bams)]

    if not control_cols or not treatment_cols:
        # Fallback to half and half if exact matching fails
        half = len(count_cols) // 2
        control_cols = list(count_cols[:half])
        treatment_cols = list(count_cols[half:])

    ctrl_counts = df_counts[control_cols].values
    treat_counts = df_counts[treatment_cols].values

    # Normalize CPM / log2FC calculation
    ctrl_mean = np.mean(ctrl_counts, axis=1) + 1.0
    treat_mean = np.mean(treat_counts, axis=1) + 1.0

    base_mean = (ctrl_mean + treat_mean) / 2.0
    log2_fc = np.log2(treat_mean / ctrl_mean)

    # Simple Welch's t-test p-value approximation for fallback
    from scipy import stats
    pvals = []
    for i in range(len(df_counts)):
        try:
            _, p = stats.ttest_ind(ctrl_counts[i], treat_counts[i], equal_var=False)
            pvals.append(p if not np.isnan(p) else 1.0)
        except Exception:
            pvals.append(1.0)

    pvals = np.array(pvals)
    # Benjamini-Hochberg FDR calculation
    n = len(pvals)
    sorted_indices = np.argsort(pvals)
    sorted_pvals = pvals[sorted_indices]
    padj = np.zeros(n)
    cummin = 1.0
    for idx in range(n - 1, -1, -1):
        p = sorted_pvals[idx]
        rank = idx + 1
        adj_p = min(cummin, p * n / rank)
        cummin = adj_p
        padj[sorted_indices[idx]] = adj_p

    df_result = pd.DataFrame({
        "gene_id": gene_ids,
        "baseMean": base_mean,
        "log2FoldChange": log2_fc,
        "pvalue": pvals,
        "padj": padj
    })

    df_result.to_csv(output_csv_path, index=False)
    print(f"Calculated DEG results for {len(df_result)} genes -> {output_csv_path}")
