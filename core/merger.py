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
    as_fdr_cutoff: float = DEFAULT_AS_FDR_CUTOFF,
    deduplicate_genes: bool = True
) -> pl.DataFrame:
    """
    Merges DEG data and rMATS data on normalized gene_id / geneSymbol,
    and applies the 4-quadrant statistical classification algorithm.
    When deduplicate_genes=True, keeps the primary representative event per gene (for 4-quadrant gene plots).
    When deduplicate_genes=False, retains ALL alternative splicing events across all genes (for event-level isoform analysis).
    """
    if (df_deg is None or df_deg.height == 0) and (df_rmats is None or df_rmats.height == 0):
        return pl.DataFrame()

    if df_deg is None or df_deg.height == 0:
        df_deg = pl.DataFrame(schema={
            "gene_id": pl.Utf8, "geneSymbol": pl.Utf8, "log2FoldChange": pl.Float64,
            "deg_pvalue": pl.Float64, "deg_fdr": pl.Float64
        })

    if df_rmats is None or df_rmats.height == 0:
        df_rmats = pl.DataFrame(schema={
            "gene_id": pl.Utf8, "geneSymbol": pl.Utf8, "event_type": pl.Utf8,
            "delta_psi": pl.Float64, "as_pvalue": pl.Float64, "as_fdr": pl.Float64,
            "coordinates": pl.Utf8
        })

    # Prepare normalized join key for df_deg
    df_deg_prepared = df_deg.with_columns([
        pl.col("gene_id").cast(pl.Utf8).str.split(".").list.first().alias("gene_id_clean"),
        pl.col("geneSymbol").cast(pl.Utf8).str.split(".").list.first().alias("symbol_clean")
    ]).with_columns([
        pl.coalesce([pl.col("symbol_clean"), pl.col("gene_id_clean")]).alias("join_key")
    ])

    # Prepare normalized join key for df_rmats
    df_rmats_prepared = df_rmats.with_columns([
        pl.col("gene_id").cast(pl.Utf8).str.split(".").list.first().alias("gene_id_clean"),
        pl.col("geneSymbol").cast(pl.Utf8).alias("symbol_clean_rmats")
    ]).with_columns([
        pl.coalesce([pl.col("gene_id_clean"), pl.col("symbol_clean_rmats")]).alias("join_key")
    ]).drop(["gene_id_clean", "symbol_clean_rmats"])

    # If df_deg lacks real gene symbols (e.g. geneSymbol == gene_id), extract clean gene_id -> geneSymbol map from rMATS
    rmats_symbol_map = (
        df_rmats
        .filter(pl.col("geneSymbol").is_not_null() & ~pl.col("geneSymbol").str.starts_with("ENSG") & ~pl.col("geneSymbol").str.starts_with("AT"))
        .select([
            pl.col("gene_id").cast(pl.Utf8).str.split(".").list.first().alias("join_key"),
            pl.col("geneSymbol").alias("mapped_symbol")
        ])
        .unique(subset=["join_key"])
    )

    if rmats_symbol_map.height > 0:
        df_deg_prepared = df_deg_prepared.join(rmats_symbol_map, on="join_key", how="left")
        df_deg_prepared = df_deg_prepared.with_columns(
            pl.coalesce([pl.col("mapped_symbol"), pl.col("geneSymbol")]).alias("geneSymbol")
        ).drop("mapped_symbol")

    df_deg_prepared = df_deg_prepared.drop(["gene_id_clean", "symbol_clean"])

    # Outer join on normalized join_key
    merged = df_deg_prepared.join(df_rmats_prepared, on="join_key", how="outer", suffix="_rmats")

    # Resolve duplicate or missing columns from outer join
    final_symbol = pl.coalesce([
        pl.col("geneSymbol_rmats") if "geneSymbol_rmats" in merged.columns else pl.col("geneSymbol"),
        pl.col("geneSymbol"),
        pl.col("gene_id"),
        pl.col("join_key")
    ])
    final_gene_id = pl.coalesce([
        pl.col("gene_id"),
        pl.col("gene_id_rmats") if "gene_id_rmats" in merged.columns else pl.col("gene_id"),
        pl.col("join_key")
    ])

    merged = merged.with_columns([
        final_symbol.alias("geneSymbol"),
        final_gene_id.alias("gene_id"),
        pl.col("log2FoldChange").fill_null(0.0),
        pl.col("deg_fdr").fill_null(1.0),
        pl.col("deg_pvalue").fill_null(1.0),
        pl.col("delta_psi").fill_null(0.0),
        pl.col("as_fdr").fill_null(1.0),
        pl.col("as_pvalue").fill_null(1.0),
        pl.col("event_type").fill_null("None"),
        pl.col("coordinates").fill_null("N/A"),
        pl.col("inc_counts").fill_null(0) if "inc_counts" in merged.columns else pl.lit(0).alias("inc_counts"),
        pl.col("exc_counts").fill_null(0) if "exc_counts" in merged.columns else pl.lit(0).alias("exc_counts")
    ])

    # Remove temporary helper columns
    drop_cols = [c for c in ["join_key", "gene_id_rmats", "geneSymbol_rmats"] if c in merged.columns]
    if drop_cols:
        merged = merged.drop(drop_cols)

    # Conditionally deduplicate by geneSymbol keeping the row with lowest FDR / highest significance
    if deduplicate_genes and "geneSymbol" in merged.columns:
        merged = merged.sort(["as_fdr", "deg_fdr"], descending=[False, False]).unique(subset=["geneSymbol"], keep="first")

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
