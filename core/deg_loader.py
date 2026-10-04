import os
import polars as pl
from config import NOISE_LOG2FC_CUTOFF, NOISE_FDR_CUTOFF, DEFAULT_FILTER_NOISE

def parse_gtf_gene_map(gtf_path: str = None) -> dict:
    """
    Parses GTF file to build a mapping from gene_id (with/without version suffix) -> gene_name (geneSymbol).
    """
    gene_map = {}
    candidate_gtfs = [gtf_path] if gtf_path else []
    
    # Auto-discover reference GTF files if not explicitly provided
    for default_dir in ["human-ref", "mouse-ref", "ref"]:
        if os.path.isdir(default_dir):
            for fname in os.listdir(default_dir):
                if fname.endswith(".gtf"):
                    candidate_gtfs.append(os.path.join(default_dir, fname))

    import re
    gene_id_pattern = re.compile(r'gene_id\s+"([^"]+)"')
    gene_name_pattern = re.compile(r'(?:gene_name|gene_symbol)\s+"([^"]+)"')

    for gpath in candidate_gtfs:
        if gpath and os.path.exists(gpath):
            try:
                with open(gpath, 'r', encoding='utf-8', errors='ignore') as f:
                    for line in f:
                        if line.startswith('#'):
                            continue
                        gid_match = gene_id_pattern.search(line)
                        if gid_match:
                            gid = gid_match.group(1).strip()
                            gname_match = gene_name_pattern.search(line)
                            gname = gname_match.group(1).strip() if gname_match else gid
                            gene_map[gid] = gname
                            gene_map[gid.split('.')[0]] = gname
                if gene_map:
                    break
            except Exception:
                pass

    return gene_map

def load_deg_data(
    filepath: str,
    gtf_path: str = None,
    filter_noise: bool = DEFAULT_FILTER_NOISE,
    noise_log2fc_cutoff: float = NOISE_LOG2FC_CUTOFF,
    noise_fdr_cutoff: float = NOISE_FDR_CUTOFF
) -> pl.DataFrame:
    """
    Loads and standardizes DEG (DESeq2/edgeR) CSV/TSV output files using Polars.
    Excludes unperturbed background noise genes (|log2FC| <= cutoff & FDR >= cutoff) when filter_noise=True.
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

    # Clean text columns
    df = df.with_columns([
        pl.col("gene_id").cast(pl.Utf8).str.strip_chars('"\'').alias("gene_id")
    ])

    if "geneSymbol" not in df.columns:
        df = df.with_columns(pl.col("gene_id").alias("geneSymbol"))
    else:
        df = df.with_columns(pl.col("geneSymbol").cast(pl.Utf8).str.strip_chars('"\'').alias("geneSymbol"))

    # Map ENSEMBL IDs to Gene Symbols if GTF is available and geneSymbol equals gene_id
    gtf_map = parse_gtf_gene_map(gtf_path)
    if gtf_map:
        clean_ids = df.select(pl.col("gene_id").str.split(".").list.first()).to_series().to_list()
        mapped_symbols = [gtf_map.get(gid, orig) for gid, orig in zip(clean_ids, df.select("geneSymbol").to_series().to_list())]
        df = df.with_columns(pl.Series("geneSymbol", mapped_symbols))

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

    res = df.select(["gene_id", "geneSymbol", "log2FoldChange", "deg_pvalue", "deg_fdr"])

    if filter_noise:
        # Exclude unperturbed noise genes with zero expression shift and non-significant FDR
        is_noise = (pl.col("log2FoldChange").abs() <= noise_log2fc_cutoff) & (pl.col("deg_fdr") >= noise_fdr_cutoff)
        noise_cnt = res.filter(is_noise).height
        if noise_cnt > 0:
            res = res.filter(~is_noise)
            print(f"  [DEG Loader] Excluded {noise_cnt} unperturbed noise genes (|log2FC| <= {noise_log2fc_cutoff}, FDR >= {noise_fdr_cutoff}). Remaining active genes: {res.height}")

    return res
