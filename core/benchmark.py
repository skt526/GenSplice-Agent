"""
GenSplice-Agent High-Performance Benchmarking Engine (Polars vs Pandas)
Measures execution speed and memory footprint across transcriptomic dataset scale.
"""

import time
import tracemalloc
import polars as pl
import pandas as pd
import numpy as np

def benchmark_polars_vs_pandas(num_genes: int = 10000) -> dict:
    """
    Benchmarks inner/outer join and 4-quadrant classification performance
    between Polars and Pandas on mock DEG & rMATS datasets of `num_genes`.
    
    Returns a dictionary with execution times (ms) and peak memory usage (MB).
    """
    # 1. Generate Synthetic Datasets
    np.random.seed(42)
    gene_ids = [f"ENSG00000{i:06d}.12" for i in range(num_genes)]
    gene_symbols = [f"GENE_{i}" for i in range(num_genes)]
    log2fc = np.random.normal(0, 1.2, num_genes)
    deg_fdr = np.random.uniform(0.0001, 1.0, num_genes)
    
    delta_psi = np.random.uniform(-0.5, 0.5, num_genes)
    as_fdr = np.random.uniform(0.0001, 1.0, num_genes)
    event_types = np.random.choice(["SE", "RI", "A5SS", "A3SS", "MXE"], num_genes)
    
    # -------------------------------------------------------------
    # Polars Benchmark
    # -------------------------------------------------------------
    tracemalloc.start()
    t0 = time.perf_counter()
    
    pl_deg = pl.DataFrame({
        "gene_id": gene_ids,
        "geneSymbol": gene_symbols,
        "log2FoldChange": log2fc,
        "deg_fdr": deg_fdr,
        "deg_pvalue": deg_fdr
    })
    
    pl_rmats = pl.DataFrame({
        "gene_id": gene_ids,
        "geneSymbol": gene_symbols,
        "event_type": event_types,
        "delta_psi": delta_psi,
        "as_fdr": as_fdr,
        "as_pvalue": as_fdr,
        "coordinates": ["chr1:100-200"] * num_genes
    })
    
    # Perform Join & Classification using Polars logic
    pl_deg_clean = pl_deg.with_columns(pl.col("gene_id").str.split(".").list.first().alias("join_key"))
    pl_rmats_clean = pl_rmats.with_columns(pl.col("gene_id").str.split(".").list.first().alias("join_key"))
    
    pl_merged = pl_deg_clean.join(pl_rmats_clean, on="join_key", how="outer", suffix="_rmats")
    pl_merged = pl_merged.with_columns([
        ((pl.col("log2FoldChange").abs() >= 0.5) & (pl.col("deg_fdr") <= 0.05)).alias("is_deg_sig"),
        ((pl.col("delta_psi").abs() >= 0.1) & (pl.col("as_fdr") <= 0.05)).alias("is_as_sig")
    ]).with_columns(
        pl.when(pl.col("is_deg_sig") & pl.col("is_as_sig")).then(pl.lit("Q1"))
        .when((~pl.col("is_deg_sig")) & pl.col("is_as_sig")).then(pl.lit("Q2"))
        .when(pl.col("is_deg_sig") & (~pl.col("is_as_sig"))).then(pl.lit("Q4"))
        .otherwise(pl.lit("Q3")).alias("quadrant")
    )
    
    polars_time_ms = (time.perf_counter() - t0) * 1000
    _, polars_peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    polars_mem_mb = polars_peak_bytes / (1024 * 1024)
    
    # -------------------------------------------------------------
    # Pandas Benchmark
    # -------------------------------------------------------------
    tracemalloc.start()
    t0 = time.perf_counter()
    
    pd_deg = pd.DataFrame({
        "gene_id": gene_ids,
        "geneSymbol": gene_symbols,
        "log2FoldChange": log2fc,
        "deg_fdr": deg_fdr,
        "deg_pvalue": deg_fdr
    })
    
    pd_rmats = pd.DataFrame({
        "gene_id": gene_ids,
        "geneSymbol": gene_symbols,
        "event_type": event_types,
        "delta_psi": delta_psi,
        "as_fdr": as_fdr,
        "as_pvalue": as_fdr,
        "coordinates": ["chr1:100-200"] * num_genes
    })
    
    pd_deg["join_key"] = pd_deg["gene_id"].apply(lambda x: x.split(".")[0])
    pd_rmats["join_key"] = pd_rmats["gene_id"].apply(lambda x: x.split(".")[0])
    
    pd_merged = pd.merge(pd_deg, pd_rmats, on="join_key", how="outer", suffixes=("", "_rmats"))
    
    is_deg_sig = (pd_merged["log2FoldChange"].abs() >= 0.5) & (pd_merged["deg_fdr"] <= 0.05)
    is_as_sig = (pd_merged["delta_psi"].abs() >= 0.1) & (pd_merged["as_fdr"] <= 0.05)
    
    conditions = [
        is_deg_sig & is_as_sig,
        (~is_deg_sig) & is_as_sig,
        is_deg_sig & (~is_as_sig)
    ]
    choices = ["Q1", "Q2", "Q4"]
    pd_merged["quadrant"] = np.select(conditions, choices, default="Q3")
    
    pandas_time_ms = (time.perf_counter() - t0) * 1000
    _, pandas_peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    pandas_mem_mb = pandas_peak_bytes / (1024 * 1024)
    
    speedup = pandas_time_ms / max(polars_time_ms, 0.0001)
    
    return {
        "num_genes": num_genes,
        "polars_time_ms": round(polars_time_ms, 2),
        "pandas_time_ms": round(pandas_time_ms, 2),
        "speedup_factor": round(speedup, 1),
        "polars_mem_mb": round(polars_mem_mb, 3),
        "pandas_mem_mb": round(pandas_mem_mb, 3)
    }

def run_full_benchmark_suite() -> list:
    """
    Runs benchmark suite across multiple dataset sizes (1K, 10K, 50K).
    """
    results = []
    for scale in [1000, 10000, 50000]:
        res = benchmark_polars_vs_pandas(num_genes=scale)
        results.append(res)
    return results

if __name__ == "__main__":
    suite = run_full_benchmark_suite()
    print("=== GenSplice-Agent Polars Performance Benchmark Results ===")
    for res in suite:
        print(f"Dataset: {res['num_genes']} genes | Polars: {res['polars_time_ms']}ms ({res['polars_mem_mb']}MB) | Pandas: {res['pandas_time_ms']}ms ({res['pandas_mem_mb']}MB) | Speedup: {res['speedup_factor']}x")
