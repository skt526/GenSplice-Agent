"""
GenSplice-Agent rMATS Loader Module (Polars)
"""

import os
import glob
import polars as pl
from config import NOISE_DELTA_PSI_CUTOFF, NOISE_FDR_CUTOFF, DEFAULT_FILTER_NOISE, DEFAULT_MIN_JUNCTION_READS

EVENT_FILES = {
    "SE": "SE.MATS.JC.txt",
    "RI": "RI.MATS.JC.txt",
    "MXE": "MXE.MATS.JC.txt",
    "A5SS": "A5SS.MATS.JC.txt",
    "A3SS": "A3SS.MATS.JC.txt"
}

def load_rmats_data(
    rmats_dir: str,
    filter_noise: bool = DEFAULT_FILTER_NOISE,
    noise_delta_psi_cutoff: float = NOISE_DELTA_PSI_CUTOFF,
    noise_fdr_cutoff: float = NOISE_FDR_CUTOFF,
    min_junction_reads: int = DEFAULT_MIN_JUNCTION_READS
) -> pl.DataFrame:
    """
    Loads rMATS output files (5 event types: SE, RI, MXE, A5SS, A3SS) from rmats_dir using Polars.
    Excludes unperturbed background noise splicing events (|dPSI| <= cutoff & FDR >= cutoff) when filter_noise=True.
    Returns a unified Polars DataFrame containing all active splicing events.
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
                elif clower in ["inclevel1", "inc_level_1"]:
                    col_map[c] = "inc_level_1"
                elif clower in ["inclevel2", "inc_level_2"]:
                    col_map[c] = "inc_level_2"
                elif clower in ["ijc_sample_1", "ijc1"]:
                    col_map[c] = "ijc_sample_1"
                elif clower in ["sjc_sample_1", "sjc1"]:
                    col_map[c] = "sjc_sample_1"
                elif clower in ["ijc_sample_2", "ijc2"]:
                    col_map[c] = "ijc_sample_2"
                elif clower in ["sjc_sample_2", "sjc2"]:
                    col_map[c] = "sjc_sample_2"

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

            # Construct event-specific clean coordinate summary
            chr_col = pl.col("chr").cast(pl.Utf8) if "chr" in df.columns else pl.lit("chr1")
            
            if event_type == "SE" and all(c in df.columns for c in ["upstreamES", "upstreamEE", "exonStart_0base", "exonEnd", "downstreamES", "downstreamEE"]):
                coords_expr = pl.concat_str([
                    chr_col, pl.lit(":"),
                    pl.col("upstreamES").cast(pl.Utf8), pl.lit("-"), pl.col("upstreamEE").cast(pl.Utf8), pl.lit(":"),
                    pl.col("exonStart_0base").cast(pl.Utf8), pl.lit("-"), pl.col("exonEnd").cast(pl.Utf8), pl.lit(":"),
                    pl.col("downstreamES").cast(pl.Utf8), pl.lit("-"), pl.col("downstreamEE").cast(pl.Utf8)
                ])
            elif event_type == "RI" and all(c in df.columns for c in ["upstreamES", "upstreamEE", "riExonStart_0base", "riExonEnd", "downstreamES", "downstreamEE"]):
                coords_expr = pl.concat_str([
                    chr_col, pl.lit(":"),
                    pl.col("upstreamES").cast(pl.Utf8), pl.lit("-"), pl.col("upstreamEE").cast(pl.Utf8), pl.lit(":"),
                    pl.col("riExonStart_0base").cast(pl.Utf8), pl.lit("-"), pl.col("riExonEnd").cast(pl.Utf8), pl.lit(":"),
                    pl.col("downstreamES").cast(pl.Utf8), pl.lit("-"), pl.col("downstreamEE").cast(pl.Utf8)
                ])
            elif event_type == "MXE" and all(c in df.columns for c in ["upstreamES", "upstreamEE", "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "downstreamES", "downstreamEE"]):
                coords_expr = pl.concat_str([
                    chr_col, pl.lit(":"),
                    pl.col("upstreamES").cast(pl.Utf8), pl.lit("-"), pl.col("upstreamEE").cast(pl.Utf8), pl.lit(":"),
                    pl.col("1stExonStart_0base").cast(pl.Utf8), pl.lit("-"), pl.col("1stExonEnd").cast(pl.Utf8), pl.lit(":"),
                    pl.col("2ndExonStart_0base").cast(pl.Utf8), pl.lit("-"), pl.col("2ndExonEnd").cast(pl.Utf8), pl.lit(":"),
                    pl.col("downstreamES").cast(pl.Utf8), pl.lit("-"), pl.col("downstreamEE").cast(pl.Utf8)
                ])
            elif event_type in ["A5SS", "A3SS"] and all(c in df.columns for c in ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"]):
                coords_expr = pl.concat_str([
                    chr_col, pl.lit(":"),
                    pl.col("flankingES").cast(pl.Utf8), pl.lit("-"), pl.col("flankingEE").cast(pl.Utf8), pl.lit(":"),
                    pl.col("longExonStart_0base").cast(pl.Utf8), pl.lit("-"), pl.col("longExonEnd").cast(pl.Utf8), pl.lit(":"),
                    pl.col("shortES").cast(pl.Utf8), pl.lit("-"), pl.col("shortEE").cast(pl.Utf8)
                ])
            else:
                coord_cols = [c for c in df.columns if any(k in c.lower() for k in ["exon", "ri", "start", "end"])]
                if coord_cols:
                    sample_coords = coord_cols[:4]
                    coords_expr = pl.concat_str([chr_col, pl.lit(":")] + [pl.col(c).cast(pl.Utf8) for c in sample_coords], separator="-")
                else:
                    coords_expr = pl.lit("N/A")

            df = df.with_columns(coords_expr.alias("coordinates"))

            df = df.with_columns([
                pl.lit(event_type).alias("event_type"),
                pl.col("gene_id").cast(pl.Utf8).str.strip_chars('"\''),
                pl.col("geneSymbol").cast(pl.Utf8).str.strip_chars('"\'')
            ])

            # Parse real junction read counts if present in rMATS output
            has_ijc = "IJC_SAMPLE_1" in df.columns and "SJC_SAMPLE_1" in df.columns
            if has_ijc:
                def _sum_counts(val_series):
                    sums = []
                    for s in val_series.to_list():
                        if s is None:
                            sums.append(0)
                            continue
                        try:
                            parts = [int(float(x.strip())) for x in str(s).replace('"', '').split(',') if x.strip() and x.strip() != 'None']
                            sums.append(sum(parts))
                        except Exception:
                            sums.append(0)
                    return sums

                ijc1_sum = _sum_counts(df["IJC_SAMPLE_1"])
                sjc1_sum = _sum_counts(df["SJC_SAMPLE_1"])
                ijc2_sum = _sum_counts(df["IJC_SAMPLE_2"]) if "IJC_SAMPLE_2" in df.columns else [0] * df.height
                sjc2_sum = _sum_counts(df["SJC_SAMPLE_2"]) if "SJC_SAMPLE_2" in df.columns else [0] * df.height
                
                df = df.with_columns([
                    pl.Series("inc_counts", [i1 + i2 for i1, i2 in zip(ijc1_sum, ijc2_sum)], dtype=pl.Int64),
                    pl.Series("exc_counts", [s1 + s2 for s1, s2 in zip(sjc1_sum, sjc2_sum)], dtype=pl.Int64)
                ])
            else:
                df = df.with_columns([
                    pl.lit(0).alias("inc_counts"),
                    pl.lit(0).alias("exc_counts")
                ])

            selected_cols = ["gene_id", "geneSymbol", "event_type", "delta_psi", "as_pvalue", "as_fdr", "coordinates", "inc_counts", "exc_counts"]
            if "chr" in df.columns:
                selected_cols.append("chr")
            if "strand" in df.columns:
                selected_cols.append("strand")
            for opt_c in ["inc_level_1", "inc_level_2", "ijc_sample_1", "sjc_sample_1", "ijc_sample_2", "sjc_sample_2"]:
                if opt_c in df.columns:
                    df = df.with_columns(pl.col(opt_c).cast(pl.Utf8).fill_null(""))
                    selected_cols.append(opt_c)

            sub_df = df.select(selected_cols)

            if min_junction_reads > 0:
                # Filter low-coverage junction read events (< min_junction_reads) to eliminate false-positive artifacts
                is_low_coverage = (pl.col("inc_counts") + pl.col("exc_counts")) < min_junction_reads
                sub_df = sub_df.filter(~is_low_coverage)

            if filter_noise:
                # Exclude unperturbed background noise events (|dPSI| <= cutoff and FDR >= cutoff)
                is_as_noise = (pl.col("delta_psi").abs() <= noise_delta_psi_cutoff) & (pl.col("as_fdr") >= noise_fdr_cutoff)
                sub_df = sub_df.filter(~is_as_noise)

            dfs.append(sub_df)

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
            "coordinates": pl.Series([], dtype=pl.Utf8),
            "inc_counts": pl.Series([], dtype=pl.Int64),
            "exc_counts": pl.Series([], dtype=pl.Int64)
        })

    full_df = pl.concat(dfs, how="diagonal")

    if filter_noise:
        # Final safeguard against any unperturbed noise events
        is_as_noise = (pl.col("delta_psi").abs() <= noise_delta_psi_cutoff) & (pl.col("as_fdr") >= noise_fdr_cutoff)
        noise_cnt = full_df.filter(is_as_noise).height
        if noise_cnt > 0:
            full_df = full_df.filter(~is_as_noise)
        print(f"  [rMATS Loader] Loaded {full_df.height} active splicing events (unperturbed noise events with |dPSI| <= {noise_delta_psi_cutoff} & FDR >= {noise_fdr_cutoff} excluded).")

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
