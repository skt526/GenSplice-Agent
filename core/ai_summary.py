"""
Legacy alias for core.ncbi_explorer module.
All automated queries use official NCBI & PubMed E-utilities APIs without LLM reliance.
"""
from core.ncbi_explorer import fetch_ncbi_gene_summary, fetch_pubmed_literature

__all__ = ["fetch_ncbi_gene_summary", "fetch_pubmed_literature"]
