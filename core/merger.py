"""
GenSplice-Agent 4-Quadrant Merger Engine (Polars)
"""

import polars as pl
from config import DEFAULT_LOG2FC_CUTOFF, DEFAULT_DELTA_PSI_CUTOFF, DEFAULT_DEG_FDR_CUTOFF, DEFAULT_AS_FDR_CUTOFF

def merge_deg_and_rmats(
    df_deg: pl.DataFrame,
    df_rmats: pl.DataFrame,
    log2fc_cutoff: float = DEFAULT_LOG2FC_CUTOFF,
    delta_psi_cutoff: float = DEFAULT_DELTA_PSI_CUTOFF,
    deg_fdr_cutoff: float = DEFAULT_DEG_FDR_CUTOFF,
    as_fdr_cutoff: float = DEFAULT_AS_FDR_CUTOFF
) -> pl.DataFrame:
    """
    Merges DEG data and primary rMATS data on geneSymbol / gene_id,
    and applies the 4-quadrant statistical classification algorithm.
    """
    if df_deg.height == 0 and df_rmats.height == 0:
        return pl.DataFrame()

    # Outer join to ensure no genes from either DEG or AS are lost
    merged = df_deg.join(df_rmats, on="geneSymbol", how="outer", suffix="_rmats")

    # Resolve duplicate or missing columns from outer join
    gene_id_expr = pl.coalesce([pl.col("gene_id"), pl.col("gene_id_rmats"), pl.col("geneSymbol")])
    
    merged = merged.with_columns([
        gene_id_expr.alias("gene_id"),
        pl.col("log2FoldChange").fill_null(0.0),
        pl.col("deg_fdr").fill_null(1.0),
        pl.col("deg_pvalue").fill_null(1.0),
        pl.col("delta_psi").fill_null(0.0),
        pl.col("as_fdr").fill_null(1.0),
        pl.col("as_pvalue").fill_null(1.0),
        pl.col("event_type").fill_null("None"),
        pl.col("coordinates").fill_null("N/A")
    ])

    # Clean up auxiliary join column if present
    if "gene_id_rmats" in merged.columns:
        merged = merged.drop("gene_id_rmats")

    # Evaluate significance flags
    merged = merged.with_columns([
        ( (pl.col("log2FoldChange").abs() >= log2fc_cutoff) & (pl.col("deg_fdr") <= deg_fdr_cutoff) ).alias("is_deg_sig"),
        ( (pl.col("delta_psi").abs() >= delta_psi_cutoff) & (pl.col("as_fdr") <= as_fdr_cutoff) ).alias("is_as_sig")
    ])

    # Apply 4-Quadrant Classification
    merged = merged.with_columns(
        pl.when(pl.col("is_deg_sig") & pl.col("is_as_sig"))
        .then(pl.lit("Q1"))
        .when((~pl.col("is_deg_sig")) & pl.col("is_as_sig"))
        .then(pl.lit("Q2"))
        .when(pl.col("is_deg_sig") & (~pl.col("is_as_sig")))
        .then(pl.lit("Q4"))
        .otherwise(pl.lit("Q3"))
        .alias("quadrant")
    )

    return merged

def get_quadrant_kpis(df_merged: pl.DataFrame) -> dict:
    """
    Computes summary KPI stats for total genes and count per quadrant.
    """
    if df_merged.height == 0:
        return {"total": 0, "Q1": 0, "Q2": 0, "Q3": 0, "Q4": 0}

    counts = df_merged.group_by("quadrant").len().to_dict(as_series=False)
    quad_dict = dict(zip(counts["quadrant"], counts["len"]))

    return {
        "total": df_merged.height,
        "Q1": quad_dict.get("Q1", 0),
        "Q2": quad_dict.get("Q2", 0),
        "Q3": quad_dict.get("Q3", 0),
        "Q4": quad_dict.get("Q4", 0)
    }
