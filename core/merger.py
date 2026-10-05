"""
GenSplice-Agent Differential Expression & Alternative Splicing Merger Engine (Polars)
"""

import polars as pl
from config import (
    DEFAULT_LOG2FC_CUTOFF, DEFAULT_DELTA_PSI_CUTOFF, DEFAULT_DEG_FDR_CUTOFF, DEFAULT_AS_FDR_CUTOFF,
    NOISE_LOG2FC_CUTOFF, NOISE_DELTA_PSI_CUTOFF, NOISE_FDR_CUTOFF, DEFAULT_FILTER_NOISE
)

def merge_deg_and_rmats(
    df_deg: pl.DataFrame,
    df_rmats: pl.DataFrame,
    log2fc_cutoff: float = DEFAULT_LOG2FC_CUTOFF,
    delta_psi_cutoff: float = DEFAULT_DELTA_PSI_CUTOFF,
    deg_fdr_cutoff: float = DEFAULT_DEG_FDR_CUTOFF,
    as_fdr_cutoff: float = DEFAULT_AS_FDR_CUTOFF,
    deduplicate_genes: bool = True,
    filter_noise: bool = DEFAULT_FILTER_NOISE,
    noise_log2fc_cutoff: float = NOISE_LOG2FC_CUTOFF,
    noise_delta_psi_cutoff: float = NOISE_DELTA_PSI_CUTOFF,
    noise_fdr_cutoff: float = NOISE_FDR_CUTOFF
) -> pl.DataFrame:
    """
    Merges DEG data and rMATS data on normalized gene_id / geneSymbol,
    and applies alternative splicing statistical classification (Inclusion, Exclusion, Non-Significant).
    When deduplicate_genes=True, keeps the primary representative event per gene (sorted by significance).
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

    # Clean up gene_id and geneSymbol in both datasets
    df_deg_prepared = df_deg.with_columns([
        pl.col("gene_id").cast(pl.Utf8).str.split(".").list.first().alias("gene_id_clean"),
        pl.col("geneSymbol").cast(pl.Utf8).str.split(".").list.first().alias("symbol_clean")
    ])

    df_rmats_prepared = df_rmats.with_columns([
        pl.col("gene_id").cast(pl.Utf8).str.split(".").list.first().alias("gene_id_clean"),
        pl.col("geneSymbol").cast(pl.Utf8).str.split(".").list.first().alias("symbol_clean")
    ])

    # If df_deg lacks real gene symbols (e.g. geneSymbol == gene_id), extract clean mapping from rMATS
    rmats_symbol_map = (
        df_rmats_prepared
        .filter(pl.col("symbol_clean").is_not_null() & ~pl.col("symbol_clean").str.starts_with("ENS") & ~pl.col("symbol_clean").str.starts_with("AT"))
        .select([
            pl.col("gene_id_clean"),
            pl.col("symbol_clean").alias("mapped_symbol")
        ])
        .unique(subset=["gene_id_clean"])
    )

    if rmats_symbol_map.height > 0:
        df_deg_prepared = df_deg_prepared.join(rmats_symbol_map, on="gene_id_clean", how="left")
        df_deg_prepared = df_deg_prepared.with_columns(
            pl.coalesce([pl.col("mapped_symbol"), pl.col("symbol_clean")]).alias("symbol_clean"),
            pl.coalesce([pl.col("mapped_symbol"), pl.col("geneSymbol")]).alias("geneSymbol")
        ).drop("mapped_symbol")

    # If df_deg only has symbol and no ENSEMBL ID, map from symbol -> gene_id_clean using rMATS
    rmats_id_map = (
        df_rmats_prepared
        .filter(pl.col("gene_id_clean").is_not_null() & (pl.col("gene_id_clean").str.starts_with("ENS") | pl.col("gene_id_clean").str.starts_with("AT")))
        .select([
            pl.col("symbol_clean"),
            pl.col("gene_id_clean").alias("mapped_gid")
        ])
        .unique(subset=["symbol_clean"])
    )

    if rmats_id_map.height > 0:
        df_deg_prepared = df_deg_prepared.join(rmats_id_map, on="symbol_clean", how="left")
        df_deg_prepared = df_deg_prepared.with_columns(
            pl.coalesce([pl.col("mapped_gid"), pl.col("gene_id_clean")]).alias("gene_id_clean"),
            pl.coalesce([pl.col("mapped_gid"), pl.col("gene_id")]).alias("gene_id")
        ).drop("mapped_gid")

    # Consistent canonical join key: prioritize clean gene_id, then clean symbol
    df_deg_prepared = df_deg_prepared.with_columns(
        pl.coalesce([pl.col("gene_id_clean"), pl.col("symbol_clean")]).alias("join_key")
    ).drop(["gene_id_clean", "symbol_clean"])

    df_rmats_prepared = df_rmats_prepared.with_columns(
        pl.coalesce([pl.col("gene_id_clean"), pl.col("symbol_clean")]).alias("join_key")
    ).drop(["gene_id_clean", "symbol_clean"])

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

    # Exclude completely unperturbed background noise genes (zero DEG shift AND zero AS shift with non-significant FDR)
    if filter_noise:
        is_unperturbed_noise = (
            ( (pl.col("log2FoldChange").abs() <= noise_log2fc_cutoff) & (pl.col("deg_fdr") >= noise_fdr_cutoff) ) &
            ( (pl.col("delta_psi").abs() <= noise_delta_psi_cutoff) & (pl.col("as_fdr") >= noise_fdr_cutoff) )
        )
        noise_cnt = merged.filter(is_unperturbed_noise).height
        if noise_cnt > 0:
            merged = merged.filter(~is_unperturbed_noise)
            print(f"  [4-Quadrant Merger] Excluded {noise_cnt} completely unperturbed noise genes (zero DEG & zero AS shift). Remaining active genes: {merged.height}")

    # Conditionally deduplicate by geneSymbol keeping the row with lowest FDR / highest significance
    if deduplicate_genes and "geneSymbol" in merged.columns:
        merged = merged.sort(["as_fdr", "deg_fdr"], descending=[False, False]).unique(subset=["geneSymbol"], keep="first", maintain_order=True)

    # Classify Alternative Splicing Status: Inclusion Favored, Exclusion Favored, Non-Significant
    is_as_sig = (pl.col("delta_psi").abs() >= delta_psi_cutoff) & (pl.col("as_fdr") <= as_fdr_cutoff)
    merged = merged.with_columns(
        pl.when(is_as_sig & (pl.col("delta_psi") > 0))
        .then(pl.lit("Inclusion"))
        .when(is_as_sig & (pl.col("delta_psi") < 0))
        .then(pl.lit("Exclusion"))
        .otherwise(pl.lit("Non-Significant"))
        .alias("splicing_status")
    )

    return merged

def get_splicing_kpis(
    df_splicing: pl.DataFrame,
    delta_psi_cutoff: float = DEFAULT_DELTA_PSI_CUTOFF,
    fdr_cutoff: float = DEFAULT_AS_FDR_CUTOFF
) -> dict:
    """
    Computes summary KPI stats for alternative splicing events:
    total, significant, inclusion, exclusion, non_significant.
    """
    if df_splicing is None or df_splicing.height == 0:
        return {
            "total": 0, "significant": 0, "inclusion": 0, "exclusion": 0, "non_significant": 0
        }

    # Filter out records without splicing if event_type is present
    if "event_type" in df_splicing.columns:
        valid_df = df_splicing.filter(pl.col("event_type").is_not_null() & (pl.col("event_type") != "None"))
        if valid_df.height > 0:
            df_splicing = valid_df

    total = df_splicing.height
    is_sig = (df_splicing["delta_psi"].abs() >= delta_psi_cutoff) & (df_splicing["as_fdr"] <= fdr_cutoff)
    is_inc = is_sig & (df_splicing["delta_psi"] > 0)
    is_exc = is_sig & (df_splicing["delta_psi"] < 0)

    inc_cnt = int(is_inc.sum())
    exc_cnt = int(is_exc.sum())
    sig_cnt = inc_cnt + exc_cnt
    nonsig_cnt = total - sig_cnt

    return {
        "total": total,
        "significant": sig_cnt,
        "inclusion": inc_cnt,
        "exclusion": exc_cnt,
        "non_significant": nonsig_cnt
    }
