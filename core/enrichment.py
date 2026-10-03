"""
GenSplice-Agent GO Term & KEGG Pathway Enrichment Module
Provides real-time dynamic enrichment analysis using gseapy / Enrichr API.
Enforces strict scientific data integrity: NO synthetic or mock enrichment data.
"""

import math
import os
import pandas as pd
import polars as pl
import numpy as np

# Cache dict to prevent redundant API calls
_ENRICHMENT_CACHE = {}

def fetch_enrichment(gene_list: list[str], gene_sets: list[str], top_n: int = 10) -> pd.DataFrame:
    """
    Fetches real enrichment results for a list of gene symbols from specified gene_sets via gseapy/Enrichr.
    Strictly returns genuine biological enrichment. If the API is offline or query fails,
    returns an empty DataFrame and emits a clear notice without fabricating fake pathway data.
    """
    clean_genes = [str(g).strip().upper() for g in gene_list if g and str(g).strip() and str(g).upper() != "NAN"]
    clean_genes = sorted(list(set(clean_genes)))
    
    empty_df = pd.DataFrame(columns=["Term", "Overlap", "P-value", "Adjusted P-value", "Combined Score", "Genes", "Gene_set", "log_p"])
    
    if not clean_genes:
        return empty_df
        
    cache_key = (tuple(clean_genes), tuple(sorted(gene_sets)), top_n)
    if cache_key in _ENRICHMENT_CACHE:
        return _ENRICHMENT_CACHE[cache_key]
        
    try:
        import gseapy as gp
        enr = gp.enrichr(
            gene_list=clean_genes,
            gene_sets=gene_sets,
            outdir=None,
            no_plot=True
        )
        df_res = enr.results
        if df_res is not None and not df_res.empty:
            df_res["P-value"] = pd.to_numeric(df_res["P-value"], errors="coerce").fillna(1.0)
            df_res["Adjusted P-value"] = pd.to_numeric(df_res["Adjusted P-value"], errors="coerce").fillna(1.0)
            df_res["log_p"] = -np.log10(df_res["P-value"].replace(0, 1e-15))
            df_res = df_res.sort_values(by="P-value", ascending=True).head(top_n)
            _ENRICHMENT_CACHE[cache_key] = df_res
            return df_res
    except ImportError:
        print("  [Warning] gseapy is not installed. Genuine enrichment analysis unavailable.")
    except Exception as e:
        print(f"  [Warning] Enrichr API query failed: {e}. Network offline or server unresponsive.")

    print("  [Notice] No mock enrichment will be generated; strictly returning genuine results.")
    _ENRICHMENT_CACHE[cache_key] = empty_df
    return empty_df
