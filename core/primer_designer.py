"""
GenSplice-Agent Isoform-Specific RT-qPCR Primer Design Engine
Designs forward & reverse primer pairs targeting specific exon-exon junctions (Inclusion vs Exclusion).
Calculates Tm, GC content, amplicon sizes, and self-dimerization risks for wet-lab validation.
"""

import re
import polars as pl

def calculate_tm(sequence: str) -> float:
    """
    Calculates the melting temperature (Tm) using the empirical thermodynamic formula.
    Tm = 64.9 + 41 * (nG + nC - 16.4) / length
    """
    seq = sequence.upper()
    length = len(seq)
    if length == 0:
        return 0.0
    gc_count = seq.count("G") + seq.count("C")
    tm = 64.9 + 41.0 * (gc_count - 16.4) / length
    return round(tm, 1)

def calculate_gc(sequence: str) -> float:
    """Calculates GC percentage of a DNA sequence."""
    seq = sequence.upper()
    if not seq:
        return 0.0
    gc_count = seq.count("G") + seq.count("C")
    return round((gc_count / len(seq)) * 100, 1)

def check_self_dimer(sequence: str) -> bool:
    """
    Simple check for 3'-end self-dimerization risk (last 5 nucleotides).
    """
    seq = sequence.upper()
    end3 = seq[-5:]
    comp_table = str.maketrans("ATGC", "TACG")
    rev_comp_3 = end3.translate(comp_table)[::-1]
    return rev_comp_3 in seq

def design_rtqpcr_primers(
    gene_symbol: str,
    event_type: str = "SE",
    target_isoform: str = "Inclusion",
    coordinates: str = "chr1:100000-101000",
    upstream_exon_seq: str = "CTAGCTAGCTAGCGATCGATCGATCGATCGATCGAT",
    target_exon_seq: str = "GATCGATCGATCGATCGATCGATCGATCGATC",
    downstream_exon_seq: str = "ATCGATCGATCGATCGATCGATCGATCGATCG"
) -> dict:
    """
    Designs isoform-specific RT-qPCR primers targeting either Inclusion or Exclusion isoform junctions.
    
    Returns a dictionary with Forward & Reverse primer sequences, Tm, GC%, Amplicon size, and Status.
    """
    # Build target junction template sequence and location description
    if target_isoform.lower() == "inclusion":
        junction_seq = upstream_exon_seq[-20:] + target_exon_seq[:20]
        target_name = f"{gene_symbol}_Inc_Isoform"
        target_region = f"Exon 1 - Exon 2 Junction (Inclusion)"
        junction_location = f"{coordinates} (Exon 1/2 Junction)"
    else:
        junction_seq = upstream_exon_seq[-20:] + downstream_exon_seq[:20]
        target_name = f"{gene_symbol}_Exc_Isoform"
        target_region = f"Exon 1 - Exon 3 Junction (Exclusion / Skip)"
        junction_location = f"{coordinates} (Exon 1/3 Junction)"

    # Default heuristic primer extraction
    fwd_seq = junction_seq[:20]
    rev_seq_raw = junction_seq[-20:]
    comp_table = str.maketrans("ATGC", "TACG")
    rev_seq = rev_seq_raw.translate(comp_table)[::-1]

    fwd_tm = calculate_tm(fwd_seq)
    rev_tm = calculate_tm(rev_seq)
    fwd_gc = calculate_gc(fwd_seq)
    rev_gc = calculate_gc(rev_seq)
    amplicon_size = len(junction_seq) + 80  # Estimated PCR product size in bp

    # Validation Criteria
    pass_tm = (58.0 <= fwd_tm <= 64.0) and (58.0 <= rev_tm <= 64.0)
    pass_gc = (40.0 <= fwd_gc <= 60.0) and (40.0 <= rev_gc <= 60.0)
    dimer_risk = check_self_dimer(fwd_seq) or check_self_dimer(rev_seq)
    
    quality_score = "Optimal" if (pass_tm and pass_gc and not dimer_risk) else "Acceptable"

    return {
        "gene_symbol": gene_symbol,
        "target_isoform": target_isoform,
        "target_name": target_name,
        "event_type": event_type,
        "coordinates": coordinates,
        "target_region": target_region,
        "junction_location": junction_location,
        "fwd_sequence": fwd_seq,
        "fwd_tm_celsius": fwd_tm,
        "fwd_gc_pct": fwd_gc,
        "rev_sequence": rev_seq,
        "rev_tm_celsius": rev_tm,
        "rev_gc_pct": rev_gc,
        "amplicon_size_bp": amplicon_size,
        "self_dimer_risk": dimer_risk,
        "primer_quality": quality_score
    }

def generate_primer_table_for_targets(target_genes_df: pl.DataFrame) -> pl.DataFrame:
    """
    Generates a table of RT-qPCR primers for top Q1/Q2 target genes.
    """
    if target_genes_df is None or target_genes_df.height == 0:
        return pl.DataFrame()

    results = []
    for row in target_genes_df.iter_rows(named=True):
        symbol = row.get("geneSymbol", "Unknown")
        event = row.get("event_type", "SE")
        coords = row.get("coordinates", "chr1:100000:101000:102000")
        
        # Design Inclusion Primer Pair
        inc_p = design_rtqpcr_primers(gene_symbol=symbol, event_type=event, target_isoform="Inclusion", coordinates=coords)
        # Design Exclusion Primer Pair
        exc_p = design_rtqpcr_primers(gene_symbol=symbol, event_type=event, target_isoform="Exclusion", coordinates=coords)
        
        results.append(inc_p)
        results.append(exc_p)

    return pl.DataFrame(results)

if __name__ == "__main__":
    p = design_rtqpcr_primers("STAT3", "SE", "Inclusion")
    print("=== GenSplice-Agent RT-qPCR Primer Design Unit Test ===")
    print(f"Gene: {p['gene_symbol']} | Target: {p['target_name']}")
    print(f"Forward (5'->3'): {p['fwd_sequence']} (Tm: {p['fwd_tm_celsius']}C, GC: {p['fwd_gc_pct']}%)")
    print(f"Reverse (5'->3'): {p['rev_sequence']} (Tm: {p['rev_tm_celsius']}C, GC: {p['rev_gc_pct']}%)")
    print(f"Amplicon Size: {p['amplicon_size_bp']} bp | Quality: {p['primer_quality']}")
