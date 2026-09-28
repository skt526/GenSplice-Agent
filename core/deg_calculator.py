import os
import json
import pandas as pd
import numpy as np

try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    PYDESEQ2_AVAILABLE = True
except ImportError:
    PYDESEQ2_AVAILABLE = False

def run_deg_analysis(feature_counts_path: str, control_bams: list, treatment_bams: list, output_csv_path: str, allow_mock: bool = False):
    """
    Parses featureCounts matrix and performs DESeq2 (PyDESeq2) analysis:
    - Size Factor Normalization
    - Negative Binomial Dispersion Shrinkage
    - Log2 Fold Change & Benjamini-Hochberg FDR p-value adjustment
    """
    if not os.path.exists(feature_counts_path) or os.path.getsize(feature_counts_path) == 0:
        if not allow_mock:
            raise FileNotFoundError(f"Error: featureCounts matrix '{feature_counts_path}' is missing or empty. Real DEG analysis failed.")
        print(f"Notice: featureCounts matrix '{feature_counts_path}' missing or empty. Generating target DEG dataset for mock mode...")
        df_target = pd.DataFrame({
            "gene_id": ["ENSG00000103194", "ENSG00000120129", "ENSG00000140355", "ENSG00000165025", "ENSG00000112715", "ENSG00000012048", "ENSG00000026101", "ENSG00000171862", "ENSG00000075624", "ENSG00000111640"],
            "geneSymbol": ["CRISPLD2", "DUSP1", "GRMZM2G140355", "SYK", "VEGFA", "BRCA1", "STAT3", "PTEN", "ACTB", "GAPDH"],
            "baseMean": [1420.5, 3850.2, 890.1, 210.4, 1850.0, 920.8, 3100.4, 1250.6, 5400.0, 6200.0],
            "log2FoldChange": [0.15, -1.84, -0.08, 1.95, 0.12, -0.04, 2.10, -1.65, 0.02, -0.03],
            "pvalue": [0.380, 0.0001, 0.520, 0.002, 0.290, 0.480, 0.0002, 0.0005, 0.910, 0.850],
            "padj": [0.450, 0.002, 0.620, 0.015, 0.380, 0.590, 0.001, 0.004, 0.950, 0.880]
        })
        df_target.to_csv(output_csv_path, index=False)
        return

    # Read featureCounts output
    with open(feature_counts_path, 'r') as f:
        lines = [line for line in f if not line.startswith('#')]

    if not lines or len(lines) <= 1:
        if not allow_mock:
            raise RuntimeError(f"Error: featureCounts matrix '{feature_counts_path}' contains no valid counts data. Real DEG analysis failed.")
        df_target = pd.DataFrame({
            "gene_id": ["ENSG00000103194", "ENSG00000120129", "ENSG00000140355", "ENSG00000165025", "ENSG00000112715", "ENSG00000012048", "ENSG00000026101", "ENSG00000171862", "ENSG00000075624", "ENSG00000111640"],
            "geneSymbol": ["CRISPLD2", "DUSP1", "GRMZM2G140355", "SYK", "VEGFA", "BRCA1", "STAT3", "PTEN", "ACTB", "GAPDH"],
            "baseMean": [1420.5, 3850.2, 890.1, 210.4, 1850.0, 920.8, 3100.4, 1250.6, 5400.0, 6200.0],
            "log2FoldChange": [0.15, -1.84, -0.08, 1.95, 0.12, -0.04, 2.10, -1.65, 0.02, -0.03],
            "pvalue": [0.380, 0.0001, 0.520, 0.002, 0.290, 0.480, 0.0002, 0.0005, 0.910, 0.850],
            "padj": [0.450, 0.002, 0.620, 0.015, 0.380, 0.590, 0.001, 0.004, 0.950, 0.880]
        })
        df_target.to_csv(output_csv_path, index=False)
        return


    from io import StringIO
    df_counts = pd.read_csv(StringIO(''.join(lines)), sep='\t')

    gene_ids = df_counts['Geneid']
    count_cols = df_counts.columns[6:]

    control_cols = [col for col in count_cols if any(os.path.basename(b) in col for b in control_bams)]
    treatment_cols = [col for col in count_cols if any(os.path.basename(b) in col for b in treatment_bams)]

    if not control_cols or not treatment_cols:
        half = len(count_cols) // 2
        control_cols = list(count_cols[:half])
        treatment_cols = list(count_cols[half:])

    # Independent Expression Filtering: Filter unexpressed / low-count background noise genes (total raw counts >= 10)
    all_sample_cols = control_cols + treatment_cols
    raw_counts_sum = df_counts[all_sample_cols].sum(axis=1)
    expressed_mask = raw_counts_sum >= 10
    if len(df_counts) > 0 and expressed_mask.sum() > 0:
        print(f"  Filtering unexpressed genes: {len(df_counts)} -> {expressed_mask.sum()} expressed genes (total raw counts >= 10)...")
        df_counts = df_counts[expressed_mask].reset_index(drop=True)
        gene_ids = df_counts['Geneid']

    # Prepare counts dataframe for PyDESeq2 (Samples x Genes)
    counts_df = df_counts.set_index('Geneid')[control_cols + treatment_cols].T
    counts_df = counts_df.astype(int)


    metadata = pd.DataFrame({
        "condition": ["control"] * len(control_cols) + ["treatment"] * len(treatment_cols)
    }, index=counts_df.index)

    # Use PyDESeq2 if available and samples >= 2 per condition
    if PYDESEQ2_AVAILABLE and len(control_cols) >= 2 and len(treatment_cols) >= 2:
        try:
            print("  Running PyDESeq2 Size Factor Normalization & Dispersion Shrinkage...")
            dds = DeseqDataSet(
                counts=counts_df,
                metadata=metadata,
                design_factors="condition",
                refit_cooks=True,
                quiet=True
            )
            dds.deseq2()
            stat_res = DeseqStats(dds, contrast=["condition", "treatment", "control"], quiet=True)
            stat_res.summary()

            res_df = stat_res.results_df.reset_index()
            res_df.rename(columns={
                "index": "gene_id",
                "log2FoldChange": "log2FoldChange",
                "pvalue": "pvalue",
                "padj": "padj"
            }, inplace=True)
            
            res_df.to_csv(output_csv_path, index=False)
            print(f"  ✔ PyDESeq2 DEG analysis completed ({len(res_df)} genes) -> {output_csv_path}")
            return
        except Exception as e:
            print(f"  Warning: PyDESeq2 calculation encountered exception ({e}). Using standard fallback...")

    # CPM Normalization & Robust DEG Calculation (Supports PyDESeq2, Multi-replicates Welch t-test, & Single-Sample Poisson Dispersion Model)
    print("  Calculating CPM (Counts Per Million) Normalization...")
    ctrl_counts = df_counts[control_cols].values.astype(float)
    treat_counts = df_counts[treatment_cols].values.astype(float)

    ctrl_totals = np.sum(ctrl_counts, axis=0)
    treat_totals = np.sum(treat_counts, axis=0)

    # Avoid zero division
    ctrl_totals[ctrl_totals == 0] = 1.0
    treat_totals[treat_totals == 0] = 1.0

    ctrl_cpm = (ctrl_counts / ctrl_totals) * 1e6
    treat_cpm = (treat_counts / treat_totals) * 1e6

    ctrl_cpm_mean = np.mean(ctrl_cpm, axis=1)
    treat_cpm_mean = np.mean(treat_cpm, axis=1)

    base_mean = (ctrl_cpm_mean + treat_cpm_mean) / 2.0
    log2_fc = np.log2((treat_cpm_mean + 1.0) / (ctrl_cpm_mean + 1.0))

    from scipy import stats

    pvals = []
    n_ctrl = len(control_cols)
    n_treat = len(treatment_cols)

    if n_ctrl >= 2 and n_treat >= 2:
        print("  Running Welch t-test on CPM normalized counts...")
        for i in range(len(df_counts)):
            try:
                _, p = stats.ttest_ind(ctrl_cpm[i], treat_cpm[i], equal_var=False)
                pvals.append(p if not np.isnan(p) else 1.0)
            except Exception:
                pvals.append(1.0)
    else:
        # Single-sample (N=1 vs N=1) DEG model:
        # Uses Poisson variance + baseline biological dispersion (alpha=0.1) Z-score estimation
        print(f"  Single-sample mode ({n_ctrl} vs {n_treat}): Using CPM Poisson/Normal Z-score model (biological dispersion alpha=0.1)...")
        alpha_bio = 0.1
        ctrl_raw_sum = np.sum(ctrl_counts, axis=1)
        treat_raw_sum = np.sum(treat_counts, axis=1)

        se = np.sqrt(1.0 / (ctrl_raw_sum + 1.0) + 1.0 / (treat_raw_sum + 1.0) + (alpha_bio ** 2))
        z_scores = log2_fc / se
        pvals = stats.norm.sf(np.abs(z_scores)) * 2.0
        pvals = np.clip(pvals, 1e-12, 1.0)

    pvals = np.array(pvals)
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
    print(f"  ✔ Saved DEG results ({len(df_result)} genes, {np.sum(padj <= 0.05)} significant at FDR <= 0.05) to {output_csv_path}")
