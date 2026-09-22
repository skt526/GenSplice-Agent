"""
GenSplice-Agent rMATS Loader Module (Polars)
"""

import os
import glob
import polars as pl

EVENT_FILES = {
    "SE": "SE.MATS.JC.txt",
    "RI": "RI.MATS.JC.txt",
    "MXE": "MXE.MATS.JC.txt",
    "A5SS": "A5SS.MATS.JC.txt",
    "A3SS": "A3SS.MATS.JC.txt"
}

def load_rmats_data(rmats_dir: str) -> pl.DataFrame:
    """
    Loads rMATS output files (5 event types: SE, RI, MXE, A5SS, A3SS) from rmats_dir using Polars.
    Returns a unified Polars DataFrame containing all splicing events.
    """
    if not os.path.exists(rmats_dir):
        raise FileNotFoundError(f"rMATS directory not found: {rmats_dir}")

    dfs = []
    
    for event_type, filename in EVENT_FILES.items():
        filepath = os.path.join(rmats_dir, filename)
        if not os.path.exists(filepath):
            # Check for alternative naming if necessary
            matches = glob.glob(os.path.join(rmats_dir, f"*{event_type}*JC*.txt"))
            if matches:
                filepath = matches[0]
            else:
                continue

        try:
            df = pl.read_csv(filepath, separator="\t", ignore_errors=True)
            if df.height == 0:
                continue

            # Column standardization
            col_map = {}
            for c in df.columns:
                clower = c.lower()
                if clower in ["geneid", "gene_id"]:
                    col_map[c] = "gene_id"
                elif clower in ["genesymbol", "symbol"]:
                    col_map[c] = "geneSymbol"
                elif clower in ["incleveldifference", "deltapsi", "dpsi"]:
                    col_map[c] = "delta_psi"
                elif clower == "fdr":
                    col_map[c] = "as_fdr"
                elif clower in ["pvalue", "pval"]:
                    col_map[c] = "as_pvalue"
                elif clower == "chr":
                    col_map[c] = "chr"
                elif clower == "strand":
                    col_map[c] = "strand"

            df = df.rename(col_map)

            if "gene_id" not in df.columns:
                first_col = df.columns[0]
                df = df.rename({first_col: "gene_id"})

            if "geneSymbol" not in df.columns:
                df = df.with_columns(pl.col("gene_id").alias("geneSymbol"))

            if "delta_psi" not in df.columns:
                df = df.with_columns(pl.lit(0.0).alias("delta_psi"))
            else:
                df = df.with_columns(pl.col("delta_psi").cast(pl.Float64, strict=False).fill_null(0.0))

            if "as_fdr" not in df.columns:
                df = df.with_columns(pl.lit(1.0).alias("as_fdr"))
            else:
                df = df.with_columns(pl.col("as_fdr").cast(pl.Float64, strict=False).fill_null(1.0))

            if "as_pvalue" not in df.columns:
                df = df.with_columns(pl.col("as_fdr").alias("as_pvalue"))

            # Construct coordinate summary
            coord_cols = [c for c in df.columns if any(k in c.lower() for k in ["exon", "ri", "start", "end"])]
            if coord_cols:
                # Pick up to first 2-3 coordinate columns for clean coordinate representation
                sample_coords = coord_cols[:3]
                df = df.with_columns(
                    pl.concat_str([pl.col(c).cast(pl.Utf8) for c in sample_coords], separator="-").alias("coordinates")
                )
            else:
                df = df.with_columns(pl.lit("N/A").alias("coordinates"))

            df = df.with_columns([
                pl.lit(event_type).alias("event_type"),
                pl.col("gene_id").cast(pl.Utf8).str.strip_chars('"\''),
                pl.col("geneSymbol").cast(pl.Utf8).str.strip_chars('"\'')
            ])

            selected_cols = ["gene_id", "geneSymbol", "event_type", "delta_psi", "as_pvalue", "as_fdr", "coordinates"]
            if "chr" in df.columns:
                selected_cols.append("chr")
            if "strand" in df.columns:
                selected_cols.append("strand")

            dfs.append(df.select(selected_cols))

        except Exception as e:
            print(f"Warning: Failed to parse rMATS file {filepath}: {e}")

    if not dfs:
        # Return empty schema-compliant DataFrame
        return pl.DataFrame({
            "gene_id": pl.Series([], dtype=pl.Utf8),
            "geneSymbol": pl.Series([], dtype=pl.Utf8),
            "event_type": pl.Series([], dtype=pl.Utf8),
            "delta_psi": pl.Series([], dtype=pl.Float64),
            "as_pvalue": pl.Series([], dtype=pl.Float64),
            "as_fdr": pl.Series([], dtype=pl.Float64),
            "coordinates": pl.Series([], dtype=pl.Utf8)
        })

    full_df = pl.concat(dfs, how="diagonal")
    return full_df

def select_primary_splicing_events(df_rmats: pl.DataFrame) -> pl.DataFrame:
    """
    Selects the primary (most significant) splicing event for each gene:
    Ranked by highest |delta_psi| and lowest as_fdr.
    """
    if df_rmats.height == 0:
        return df_rmats

    df_ranked = df_rmats.with_columns(
        pl.col("delta_psi").abs().alias("abs_delta_psi")
    ).sort(["geneSymbol", "as_fdr", "abs_delta_psi"], descending=[False, False, True])

    # Group by geneSymbol and pick the top event
    primary_df = df_ranked.unique(subset=["geneSymbol"], keep="first")
    return primary_df.drop("abs_delta_psi")
