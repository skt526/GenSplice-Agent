"""
GenSplice-Agent DEG Loader Module (Polars)
"""

import os
import polars as pl

def load_deg_data(filepath: str) -> pl.DataFrame:
    """
    Loads and standardizes DEG (DESeq2/edgeR) CSV/TSV output files using Polars.
    Returns a standardized Polars DataFrame containing:
    ['gene_id', 'geneSymbol', 'log2FoldChange', 'deg_pvalue', 'deg_fdr']
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"DEG file not found: {filepath}")

    # Auto-detect delimiter
    delimiter = "\t" if filepath.endswith((".tsv", ".txt")) else ","
    
    df = pl.read_csv(filepath, separator=delimiter, ignore_errors=True)

    # Standardize column mapping
    col_mapping = {}
    cols = df.columns
    
    for c in cols:
        clower = c.lower()
        if clower in ["gene_id", "geneid", "gene"]:
            col_mapping[c] = "gene_id"
        elif clower in ["genesymbol", "symbol", "gene_name"]:
            col_mapping[c] = "geneSymbol"
        elif clower in ["log2foldchange", "logfc", "log2fc"]:
            col_mapping[c] = "log2FoldChange"
        elif clower in ["padj", "fdr", "adj.p.val", "qvalue"]:
            col_mapping[c] = "deg_fdr"
        elif clower in ["pvalue", "p.val", "pval"]:
            col_mapping[c] = "deg_pvalue"

    df = df.rename(col_mapping)

    # Required columns check & fallback
    if "gene_id" not in df.columns:
        first_col = df.columns[0]
        df = df.rename({first_col: "gene_id"})

    if "geneSymbol" not in df.columns:
        # Strip quote marks or dot extensions from ENSEMBL IDs if present
        df = df.with_columns(pl.col("gene_id").alias("geneSymbol"))

    if "log2FoldChange" not in df.columns:
        df = df.with_columns(pl.lit(0.0).alias("log2FoldChange"))
    else:
        df = df.with_columns(pl.col("log2FoldChange").cast(pl.Float64, strict=False).fill_null(0.0))

    if "deg_fdr" not in df.columns:
        if "deg_pvalue" in df.columns:
            df = df.with_columns(pl.col("deg_pvalue").cast(pl.Float64, strict=False).alias("deg_fdr"))
        else:
            df = df.with_columns(pl.lit(1.0).alias("deg_fdr"))
    else:
        df = df.with_columns(pl.col("deg_fdr").cast(pl.Float64, strict=False).fill_null(1.0))

    if "deg_pvalue" not in df.columns:
        df = df.with_columns(pl.col("deg_fdr").alias("deg_pvalue"))

    # Clean text columns
    df = df.with_columns([
        pl.col("gene_id").cast(pl.Utf8).str.strip_chars('"\''),
        pl.col("geneSymbol").cast(pl.Utf8).str.strip_chars('"\'')
    ])

    return df.select(["gene_id", "geneSymbol", "log2FoldChange", "deg_pvalue", "deg_fdr"])
