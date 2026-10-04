"""
GenSplice-Agent Standalone HTML Exporter Module
Light Mode interactive HTML report with:
- Dedicated Alternative Splicing Volcano Plot & Target Selector (ΔPSI vs -log₁₀ FDR)
- Splicing Direction Classification (Inclusion Favored, Exclusion Favored, Non-Significant Background)
- Real-Time Dynamic Sliders (ΔPSI effect size & rMATS FDR significance)
- Synchronized NCBI Gene Explorer, Isoform Annotation Matrix, Sashimi Plot, and Primer Designer
- GO Term & KEGG Pathway Enrichment Charts & Tables focused on Significant Splicing Targets
"""

import os
import json
import time
import math
import polars as pl
import pandas as pd
import plotly.graph_objects as go
from visualizer.as_volcano_plot import build_as_volcano_plot
from visualizer.quadrant_plot import build_quadrant_plot
from visualizer.enrichment_plot import build_enrichment_chart
from core.merger import get_splicing_kpis, get_quadrant_kpis
from core.enrichment import fetch_enrichment
from core.isoform_annotator import annotate_isoform_events
from visualizer.exon_structure import plot_exon_structure, resolve_gene_exon_coords
from core.primer_designer import generate_primer_table_for_targets

def export_html_report(
    df_merged: pl.DataFrame,
    output_html_path: str = "outputs/gensplice_report.html",
    ai_insights: str = None,
    log2fc_cutoff: float = 1.0,
    delta_psi_cutoff: float = 0.1,
    deg_fdr_cutoff: float = 0.05,
    as_fdr_cutoff: float = 0.05,
    df_all_events: pl.DataFrame = None,
    fasta_path: str = None,
    gtf_path: str = None,
    organism: str = "Homo sapiens"
) -> str:
    """
    Exports a self-contained Light Mode interactive HTML report with an Alternative Splicing Volcano Plot
    and downstream functional validation workflows.
    """
    os.makedirs(os.path.dirname(output_html_path), exist_ok=True)

    if df_all_events is None:
        df_all_events = df_merged

    # 1. Filter for alternative splicing records
    if "event_type" in df_merged.columns:
        df_splicing_events = df_merged.filter(
            pl.col("event_type").is_not_null() & (pl.col("event_type") != "None")
        )
        if df_splicing_events.height == 0:
            df_splicing_events = df_merged
    else:
        df_splicing_events = df_merged

    # 2. Build Alternative Splicing Volcano Plot
    fig_volcano = build_as_volcano_plot(
        df_splicing=df_splicing_events,
        delta_psi_cutoff=delta_psi_cutoff,
        fdr_cutoff=as_fdr_cutoff
    )
    volcano_html = fig_volcano.to_html(full_html=False, include_plotlyjs="cdn", div_id="plotly-as-volcano-div")

    # 3. KPI metrics for Alternative Splicing Volcano Plot
    kpis = get_splicing_kpis(df_splicing_events, delta_psi_cutoff=delta_psi_cutoff, fdr_cutoff=as_fdr_cutoff)

    # 4. Identify initial top significant splicing target for Visual Exon-Intron Sashimi Plot
    top_gene_symbol = ""
    top_coords = "N/A"
    top_event = "N/A"
    top_delta_psi = 0.0
    top_inc_counts = 0
    top_exc_counts = 0

    # Look for significant events first (|dPSI| >= cutoff and as_fdr <= cutoff)
    sig_events = df_splicing_events.filter(
        (pl.col("delta_psi").abs() >= delta_psi_cutoff) &
        (pl.col("as_fdr") <= as_fdr_cutoff)
    ).sort(["as_fdr", "delta_psi"], descending=[False, True])

    if sig_events.height > 0:
        row0 = sig_events.to_dicts()[0]
        top_gene_symbol = row0.get("geneSymbol") or ""
        top_coords = row0.get("coordinates") or "N/A"
        top_event = row0.get("event_type") or "SE"
        top_delta_psi = float(row0.get("delta_psi", 0.0) if row0.get("delta_psi") is not None else 0.0)
        top_inc_counts = int(row0.get("inc_counts", 0) if row0.get("inc_counts") is not None else 0)
        top_exc_counts = int(row0.get("exc_counts", 0) if row0.get("exc_counts") is not None else 0)

        fig_sashimi = plot_exon_structure(
            gene_symbol=top_gene_symbol,
            event_type=top_event,
            coordinates=top_coords,
            inc_counts=top_inc_counts,
            exc_counts=top_exc_counts,
            delta_psi=top_delta_psi
        )
        sashimi_html = fig_sashimi.to_html(full_html=False, include_plotlyjs=False, div_id="plotly-sashimi-div")
        sashimi_title_text = f"<b>{top_gene_symbol} Exon Structure & Splicing Sashimi Plot</b> ({top_event} | ΔPSI = {top_delta_psi:+.2f})"
    else:
        # Informative zero-state sashimi canvas when no genes meet cutoff
        fig_sashimi = go.Figure()
        fig_sashimi.add_annotation(
            text=f"<b>No alternative splicing targets meet current significance cutoffs (|ΔPSI| ≥ {delta_psi_cutoff:.2f}, FDR ≤ {as_fdr_cutoff:.3f})</b><br><span style='color:#64748B;'>Adjust ΔPSI or FDR cutoff sliders in the control panel to identify candidate splicing events</span>",
            xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="#64748B")
        )
        fig_sashimi.update_layout(
            template="plotly_white", paper_bgcolor="#FFFFFF", plot_bgcolor="#F8FAFC",
            xaxis=dict(showgrid=False, showticklabels=False), yaxis=dict(showgrid=False, showticklabels=False),
            height=260
        )
        sashimi_html = fig_sashimi.to_html(full_html=False, include_plotlyjs=False, div_id="plotly-sashimi-div")
        sashimi_title_text = "<b>No Alternative Splicing Targets Meeting Cutoffs</b>"

    # 5. Extract Candidate Splicing Target Genes for downstream analyses
    candidate_symbols = set()
    if sig_events.height > 0 and "geneSymbol" in sig_events.columns:
        candidate_symbols.update(
            g.strip() for g in sig_events["geneSymbol"].drop_nulls().to_list() if g and g.strip()
        )

    # Filter df_all_events: strictly preserve candidate splicing genes or empty state (no synthetic fallbacks)
    if candidate_symbols and "geneSymbol" in df_all_events.columns:
        target_candidates = df_all_events.filter(pl.col("geneSymbol").is_in(list(candidate_symbols)))
    else:
        target_candidates = df_all_events.clear()

    # Exclude non-splicing records (e.g. event_type is None or "None")
    if "event_type" in target_candidates.columns:
        target_candidates_as = target_candidates.filter(
            pl.col("event_type").is_not_null() & (pl.col("event_type") != "None")
        )
        if target_candidates_as.height > 0:
            target_candidates = target_candidates_as

    # 6. Splicing metadata map for client-side interactive Sashimi Plot re-rendering
    sashimi_dict = {}
    sashimi_events_map = {}
    for r in target_candidates.iter_rows(named=True):
        sym = r.get("geneSymbol")
        if sym:
            sym_clean = sym.strip()
            sym_upper = sym_clean.upper()
            ev = r.get("event_type", "SE")
            coords = r.get("coordinates", "")
            dpsi = float(r.get("delta_psi", 0.0) if r.get("delta_psi") is not None else 0.0)
            c_info = resolve_gene_exon_coords(sym_clean, ev, coords)
            inc_c = int(r.get("inc_counts", 0) if r.get("inc_counts") is not None else 0)
            exc_c = int(r.get("exc_counts", 0) if r.get("exc_counts") is not None else 0)
            entry = {
                "gene_symbol": sym_clean,
                "event_type": ev,
                "coordinates": coords if coords and coords != "N/A" else f"{c_info['chrom']}:{c_info['ex1'][0]}-{c_info['ex1'][1]}:{c_info['ex2'][0]}-{c_info['ex2'][1]}:{c_info['ex3'][0]}-{c_info['ex3'][1]}",
                "chrom": c_info["chrom"],
                "delta_psi": dpsi,
                "inc_counts": inc_c,
                "exc_counts": exc_c,
                "ex1": list(c_info["ex1"]),
                "ex2": list(c_info["ex2"]),
                "ex3": list(c_info["ex3"]),
                "nums": c_info.get("nums", [])
            }
            if sym_upper not in sashimi_events_map:
                sashimi_events_map[sym_upper] = []
            entry["event_index"] = len(sashimi_events_map[sym_upper])
            sashimi_events_map[sym_upper].append(entry)

            if sym_upper not in sashimi_dict or abs(dpsi) > abs(sashimi_dict[sym_upper].get("delta_psi", 0.0)):
                sashimi_dict[sym_clean] = entry
                sashimi_dict[sym_upper] = entry

    if top_gene_symbol and top_gene_symbol not in sashimi_dict:
        c_info = resolve_gene_exon_coords(top_gene_symbol, top_event, top_coords)
        entry = {
            "gene_symbol": top_gene_symbol,
            "event_type": top_event,
            "coordinates": top_coords,
            "chrom": c_info["chrom"],
            "delta_psi": float(top_delta_psi),
            "inc_counts": top_inc_counts,
            "exc_counts": top_exc_counts,
            "ex1": list(c_info["ex1"]),
            "ex2": list(c_info["ex2"]),
            "ex3": list(c_info["ex3"]),
            "nums": c_info.get("nums", []),
            "event_index": 0
        }
        sashimi_dict[top_gene_symbol] = entry
        sashimi_dict[top_gene_symbol.upper()] = entry
        sashimi_events_map[top_gene_symbol.upper()] = [entry]

    sashimi_data_json = json.dumps(sashimi_dict)
    sashimi_events_json = json.dumps(sashimi_events_map)

    # 7. RT-qPCR Primer Designer Engine (Synchronized with Targets & Reference FASTA)
    primer_df = generate_primer_table_for_targets(target_candidates, fasta_path=fasta_path)
    primer_records = primer_df.to_dicts() if primer_df.height > 0 else []
    primer_records_json = json.dumps(primer_records)

    first_primer_gene = top_gene_symbol if top_gene_symbol else "-"
    first_coords = top_coords if top_gene_symbol else "-"
    first_event = top_event if top_gene_symbol else "-"

    primer_rows_html = ""
    if primer_df.height > 0 and top_gene_symbol:
        matching_rows = primer_df.filter(pl.col("gene_symbol") == first_primer_gene)
        if matching_rows.height == 0 and primer_records:
            first_primer_gene = primer_records[0]["gene_symbol"]
            first_coords = primer_records[0].get("coordinates", "N/A")
            first_event = primer_records[0].get("event_type", "SE")
            matching_rows = primer_df.filter(pl.col("gene_symbol") == first_primer_gene)

        # Render only the primary event's primer pairs initially
        if "event_index" in matching_rows.columns:
            display_rows = matching_rows.filter(pl.col("event_index") == 0)
            if display_rows.height == 0:
                display_rows = matching_rows.head(2)
        else:
            display_rows = matching_rows.head(2)

        for r in display_rows.iter_rows(named=True):
            isoform_badge = '<span class="badge" style="background:#EFF6FF; color:#2563EB; border:1px solid #3B82F6;">Inclusion Isoform</span>' if r['target_isoform'].lower() == 'inclusion' else '<span class="badge" style="background:#FDF2F8; color:#DB2777; border:1px solid #EC4899;">Exclusion Isoform</span>'
            fwd_tm = f"{r['fwd_tm_celsius']:.1f}°C" if r.get('fwd_tm_celsius') is not None else "N/A"
            fwd_gc = f"{r['fwd_gc_pct']:.1f}%" if r.get('fwd_gc_pct') is not None else "N/A"
            rev_tm = f"{r['rev_tm_celsius']:.1f}°C" if r.get('rev_tm_celsius') is not None else "N/A"
            rev_gc = f"{r['rev_gc_pct']:.1f}%" if r.get('rev_gc_pct') is not None else "N/A"
            amp_size = f"{r['amplicon_size_bp']} bp" if r.get('amplicon_size_bp') is not None else "N/A"
            fwd_seq = f"<code style='font-weight:700; color:#1E40AF;'>{r['fwd_sequence']}</code>" if r.get('fwd_sequence') and r['fwd_sequence'] != 'N/A' else "<span style='color:#64748B;'>N/A</span>"
            rev_seq = f"<code style='font-weight:700; color:#1E40AF;'>{r['rev_sequence']}</code>" if r.get('rev_sequence') and r['rev_sequence'] != 'N/A' else "<span style='color:#64748B;'>N/A</span>"

            primer_rows_html += f"""
            <tr>
                <td><b>{r['gene_symbol']}</b></td>
                <td>{isoform_badge}</td>
                <td><b>{r.get('target_region', r['target_isoform'] + ' Junction')}</b><br><span style="font-size:11px; color:#64748B;">Coords: <code>{r.get('coordinates', 'N/A')}</code></span></td>
                <td>{fwd_seq}<br><span style="font-size:11px; color:#64748B;">Tm: {fwd_tm} | GC: {fwd_gc}</span></td>
                <td>{rev_seq}<br><span style="font-size:11px; color:#64748B;">Tm: {rev_tm} | GC: {rev_gc}</span></td>
                <td><b>{amp_size}</b></td>
                <td><span class="badge" style="background:#F0FDF4; color:#16A34A; border:1px solid #22C55E;">{r['primer_quality']}</span></td>
            </tr>
            """
    else:
        primer_rows_html = '<tr><td colspan="7" style="text-align:center; color:#64748B; padding:16px;">No significant alternative splicing targets meet current threshold cutoffs.</td></tr>'

    # 8. Splicing Target Genes for GO / KEGG Enrichment
    splicing_genes = list(candidate_symbols) if candidate_symbols else []
    if not splicing_genes:
        splicing_genes = df_splicing_events.select("geneSymbol").to_series().to_list()
    all_genes = df_merged.select("geneSymbol").to_series().to_list()

    # Pre-generate GO & KEGG Enrichment for Significant Splicing Targets
    kegg_lib = "KEGG_2019_Mouse" if ("mouse" in str(organism).lower() or "mus" in str(organism).lower()) else "KEGG_2021_Human"
    df_go = fetch_enrichment(splicing_genes, gene_sets=["GO_Biological_Process_2023"], top_n=10)
    df_kegg = fetch_enrichment(splicing_genes, gene_sets=[kegg_lib], top_n=10)

    go_records = []
    if df_go is not None and not df_go.empty:
        for _, r in df_go.iterrows():
            pval = float(r.get('P-value', 1.0))
            adj_pval = float(r.get('Adjusted P-value', pval))
            log_p = round(-math.log10(max(pval, 1e-50)), 2) if pval > 0 else 50.0
            go_records.append({
                "term": str(r.get('Term', '')),
                "overlap": str(r.get('Overlap', '')),
                "pvalue": f"{pval:.2e}",
                "adjPvalue": f"{adj_pval:.2e}",
                "logP": log_p,
                "genes": str(r.get('Genes', ''))
            })

    kegg_records = []
    if df_kegg is not None and not df_kegg.empty:
        for _, r in df_kegg.iterrows():
            pval = float(r.get('P-value', 1.0))
            adj_pval = float(r.get('Adjusted P-value', pval))
            log_p = round(-math.log10(max(pval, 1e-50)), 2) if pval > 0 else 50.0
            kegg_records.append({
                "term": str(r.get('Term', '')),
                "overlap": str(r.get('Overlap', '')),
                "pvalue": f"{pval:.2e}",
                "adjPvalue": f"{adj_pval:.2e}",
                "logP": log_p,
                "genes": str(r.get('Genes', ''))
            })

    go_records_json = json.dumps(go_records)
    kegg_records_json = json.dumps(kegg_records)

    # 9. Compact JSON dataset for client-side JS engine (lightweight ~1.5MB)
    essential_cols = [c for c in [
        "geneSymbol", "gene_id", "delta_psi", "as_fdr", "as_pvalue",
        "event_type", "coordinates", "inc_counts", "exc_counts",
        "log2FoldChange", "deg_fdr"
    ] if c in df_splicing_events.columns]

    # Priority sorting: Significant splicing targets first -> sorted by as_fdr -> delta_psi
    df_sorted = df_splicing_events.with_columns([
        ((pl.col("delta_psi").abs() >= delta_psi_cutoff) & (pl.col("as_fdr") <= as_fdr_cutoff))
        .alias("is_sig_target")
    ]).sort(["is_sig_target", "as_fdr", "delta_psi"], descending=[True, False, True])

    df_compact = df_sorted.select(essential_cols)

    # Round floats to eliminate unnecessary JSON string precision bloat
    round_exprs = []
    if "delta_psi" in df_compact.columns:
        round_exprs.append(pl.col("delta_psi").round(4))
    if "as_fdr" in df_compact.columns:
        round_exprs.append(pl.col("as_fdr").round(6))
    if "as_pvalue" in df_compact.columns:
        round_exprs.append(pl.col("as_pvalue").round(6))
    if "log2FoldChange" in df_compact.columns:
        round_exprs.append(pl.col("log2FoldChange").round(4))
    if "deg_fdr" in df_compact.columns:
        round_exprs.append(pl.col("deg_fdr").round(6))

    if round_exprs:
        df_compact = df_compact.with_columns(round_exprs)

    raw_data_json = df_compact.to_pandas().to_json(orient="records")

    # 10. Dynamically build gene -> pathways map from actual GO and KEGG enrichment results
    gene_pathway_map = {}
    for df_enr, key in [(df_go, "go"), (df_kegg, "kegg")]:
        if df_enr is not None and not df_enr.empty:
            df_enr_pd = df_enr.to_pandas() if isinstance(df_enr, pl.DataFrame) else df_enr
            for _, row in df_enr_pd.iterrows():
                term = str(row.get("Term", ""))
                genes_str = str(row.get("Genes", ""))
                for g in genes_str.split(";"):
                    g_clean = g.strip().upper()
                    if g_clean and g_clean != "N/A":
                        if g_clean not in gene_pathway_map:
                            gene_pathway_map[g_clean] = {"go": [], "kegg": []}
                        if term not in gene_pathway_map[g_clean][key]:
                            gene_pathway_map[g_clean][key].append(term)
    gene_pathway_json = json.dumps(gene_pathway_map)

    # 11. Annotate candidate splicing target events using real GTF CDS structure
    isoform_source_df = target_candidates
    isoform_df = annotate_isoform_events(isoform_source_df, gtf_path=gtf_path)

    isoform_rows_list = []
    isoform_records_list = []
    top_sym_upper = top_gene_symbol.strip().upper() if top_gene_symbol else ""

    for _, r in isoform_df.iterrows():
        gene_sym = str(r['geneSymbol']).strip()
        gene_sym_upper = gene_sym.upper()

        if top_sym_upper and gene_sym_upper == top_sym_upper:
            nmd_badge_style = "background:#FEF2F2; color:#EF4444; border:1px solid #EF4444;" if "NMD Sensitive" in str(r['nmd_prediction']) else "background:#EFF6FF; color:#2563EB; border:1px solid #3B82F6;"
            risk_badge_style = "background:#FEF2F2; color:#DC2626; border:1px solid #EF4444; font-weight:800;" if "High" in str(r['impairment_tier']) else ("background:#FFFBEB; color:#D97706; border:1px solid #F59E0B; font-weight:800;" if "Moderate" in str(r['impairment_tier']) else "background:#F0FDF4; color:#16A34A; border:1px solid #22C55E; font-weight:800;")
            isoform_rows_list.append(f"""
                <tr class="isoform-data-row" data-gene="{gene_sym_upper}">
                    <td><b>{gene_sym}</b></td>
                    <td><span class="badge" style="{risk_badge_style}">{r['impairment_tier']}</span></td>
                    <td><span style="font-size:12px; font-weight:600; color:#334155;">{r['primary_dysfunction_cause']}</span></td>
                    <td><span class="badge" style="background:#F1F5F9; color:#334155; font-weight:700;">{r['event_type']}</span></td>
                    <td><code>{r['transcript_id']}</code></td>
                    <td><span style="font-size:12px; color:#64748B;">{r['coordinates']}</span></td>
                    <td><b style="color:#2563EB;">{r['delta_psi']:.3f}</b></td>
                    <td><b style="color:{'#DC2626' if r['log2FoldChange'] > 0 else '#2563EB'};">{r['log2FoldChange']:.3f}</b></td>
                    <td>{r['cds_frame']}</td>
                    <td><span class="badge" style="{nmd_badge_style}">{r['nmd_prediction']}</span><br><span style="font-size:11px; color:#64748B;">{r['ptc_position']}</span></td>
                    <td><b>{r['protein_domain']}</b></td>
                    <td><span style="font-size:12px; color:#475569;">{r['localization_consequence']}</span></td>
                </tr>
            """)

        isoform_records_list.append({
            "gene_symbol": gene_sym,
            "impairment_tier": str(r.get("impairment_tier", "")),
            "primary_dysfunction_cause": str(r.get("primary_dysfunction_cause", "")),
            "event_type": str(r.get("event_type", "")),
            "transcript_id": str(r.get("transcript_id", "")),
            "coordinates": str(r.get("coordinates", "")),
            "delta_psi": float(r.get("delta_psi", 0.0)),
            "log2FoldChange": float(r.get("log2FoldChange", 0.0)),
            "cds_frame": str(r.get("cds_frame", "")),
            "nmd_prediction": str(r.get("nmd_prediction", "")),
            "ptc_position": str(r.get("ptc_position", "")),
            "protein_domain": str(r.get("protein_domain", "")),
            "localization_consequence": str(r.get("localization_consequence", ""))
        })

    if not isoform_rows_list:
        isoform_rows_list.append(f"""
            <tr><td colspan="12" style="text-align:center; color:#64748B; padding:20px; font-weight:600;">Active splicing target records loaded ({len(isoform_records_list)} events). Select a gene to view isoform details.</td></tr>
        """)

    isoform_table_rows_html = "\n".join(isoform_rows_list)
    isoform_records_json = json.dumps(isoform_records_list)

    # 12. Generate Self-Contained HTML Document
    full_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GenSplice-Agent Interactive Light Mode Report</title>
    <style>
        :root {{
            --bg-body: #F8FAFC;
            --bg-card: #FFFFFF;
            --border-color: #E2E8F0;
            --text-main: #0F172A;
            --text-muted: #64748B;
            --blue-inc: #2563EB;
            --red-exc: #E11D48;
            --purple-target: #8B5CF6;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            background-color: var(--bg-body);
            color: var(--text-main);
            margin: 0;
            padding: 24px;
            line-height: 1.5;
        }}
        .header {{
            background: linear-gradient(135deg, #FFFFFF 0%, #F1F5F9 100%);
            border: 1px solid var(--border-color);
            color: var(--text-main);
            padding: 28px 32px;
            border-radius: 16px;
            margin-bottom: 24px;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.05);
        }}
        .header h1 {{ margin: 0 0 8px 0; font-size: 28px; font-weight: 800; color: #0F172A; }}
        .header p {{ margin: 0; color: var(--text-muted); font-size: 15px; }}

        .kpi-container {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 16px;
            margin-top: 20px;
            margin-bottom: 12px;
        }}
        .kpi-card {{
            background: var(--bg-card);
            border-radius: 12px;
            padding: 18px;
            border: 1px solid var(--border-color);
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.05);
            text-align: center;
        }}
        .kpi-card .value {{ font-size: 32px; font-weight: 800; margin-top: 4px; }}
        .kpi-card .label {{ font-size: 12px; color: var(--text-muted); text-transform: uppercase; letter-spacing: 0.5px; font-weight: 700; }}
        
        .card {{
            background: var(--bg-card);
            border-radius: 16px;
            padding: 24px;
            margin-bottom: 24px;
            border: 1px solid var(--border-color);
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.05);
            transition: all 0.3s ease;
        }}
        .card h2 {{ margin-top: 0; font-size: 20px; color: var(--text-main); border-bottom: 2px solid var(--border-color); padding-bottom: 12px; }}

        @keyframes card-glow {{
            0% {{ box-shadow: 0 0 0 rgba(139, 92, 246, 0); transform: scale(1); }}
            50% {{ box-shadow: 0 0 25px rgba(139, 92, 246, 0.45); transform: scale(1.008); }}
            100% {{ box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.05); transform: scale(1); }}
        }}
        .card-updating {{
            animation: card-glow 0.65s cubic-bezier(0.4, 0, 0.2, 1);
        }}

        .chart-layout-grid {{
            display: grid;
            grid-template-columns: 1fr 340px;
            gap: 24px;
            align-items: start;
        }}

        .right-control-panel {{
            background: #F8FAFC;
            border-radius: 14px;
            padding: 20px;
            border: 2px solid #3B82F6;
            box-shadow: 0 4px 12px rgba(59, 130, 246, 0.1);
            display: flex;
            flex-direction: column;
            gap: 20px;
        }}
        .right-control-panel h3 {{
            margin: 0;
            font-size: 16px;
            color: var(--text-main);
            border-bottom: 1px solid var(--border-color);
            padding-bottom: 8px;
        }}

        .slider-group {{
            display: flex;
            flex-direction: column;
            gap: 6px;
        }}
        .slider-group label {{
            font-weight: 700;
            font-size: 13px;
            display: flex;
            justify-content: space-between;
        }}
        .slider-group input[type=range] {{
            width: 100%;
            height: 6px;
            cursor: pointer;
            border-radius: 3px;
        }}

        .slider-as input[type=range] {{ accent-color: var(--blue-inc); }}
        .slider-fdr input[type=range] {{ accent-color: #4F46E5; }}

        .btn-download-csv {{
            background-color: #3B82F6;
            color: #FFFFFF;
            border: none;
            border-radius: 8px;
            padding: 12px 16px;
            font-size: 14px;
            font-weight: 700;
            cursor: pointer;
            display: flex;
            align-items: center;
            justify-content: center;
            gap: 8px;
            box-shadow: 0 2px 6px rgba(59, 130, 246, 0.3);
            transition: all 0.2s ease;
        }}
        .btn-download-csv:hover {{
            background-color: #2563EB;
            transform: translateY(-1px);
        }}

        .btn-generate-enrichment {{
            background: linear-gradient(135deg, #8B5CF6 0%, #6D28D9 100%);
            color: #FFFFFF;
            border: none;
            border-radius: 10px;
            padding: 14px 20px;
            font-size: 15px;
            font-weight: 800;
            cursor: pointer;
            display: flex;
            align-items: center;
            justify-content: center;
            gap: 10px;
            box-shadow: 0 4px 14px rgba(139, 92, 246, 0.4);
            transition: all 0.25s cubic-bezier(0.4, 0, 0.2, 1);
            width: 100%;
            margin-top: 10px;
            position: relative;
        }}
        .btn-generate-enrichment:hover {{
            background: linear-gradient(135deg, #7C3AED 0%, #5B21B6 100%);
            transform: translateY(-2px) scale(1.02);
            box-shadow: 0 6px 20px rgba(139, 92, 246, 0.6);
        }}
        .btn-generate-enrichment:active {{
            transform: translateY(1px) scale(0.97);
        }}

        .legend-title-desc {{
            font-size: 14px;
            font-weight: 800;
            color: #0F172A;
            margin-bottom: 6px;
        }}

        .legend-guide {{
            background: #FFFFFF;
            border-radius: 8px;
            padding: 12px;
            border: 1px solid var(--border-color);
            font-size: 13px;
            display: flex;
            flex-direction: column;
            gap: 8px;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 8px;
            font-weight: 600;
        }}
        .legend-dot {{
            width: 12px;
            height: 12px;
            border-radius: 50%;
            display: inline-block;
        }}

        .data-table {{ width: 100%; border-collapse: collapse; margin-top: 0; font-size: 14px; color: var(--text-main); }}
        .data-table th, .data-table td {{ padding: 12px 16px; text-align: left; border-bottom: 1px solid var(--border-color); }}
        .data-table th {{ background-color: #F1F5F9; font-weight: 600; color: #475569; position: sticky; top: 0; z-index: 2; box-shadow: 0 1px 2px rgba(0,0,0,0.05); }}
        .badge {{ display: inline-block; padding: 3px 8px; border-radius: 6px; font-size: 12px; font-weight: 600; }}
        .badge-blue {{ background-color: rgba(59, 130, 246, 0.15); color: #3B82F6; border: 1px solid #3B82F6; }}
        .footer {{ text-align: center; color: var(--text-muted); font-size: 13px; margin-top: 40px; padding-top: 20px; border-top: 1px solid var(--border-color); }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 GenSplice-Agent Interactive Light Report</h1>
        <p>Alternative Splicing Volcano Plot & Downstream Functional Validation Engine</p>
    </div>

    <!-- Alternative Splicing Volcano Plot Section -->
    <div class="card">
        <h2>🌋 Alternative Splicing Volcano Plot & Target Selector</h2>
        <div class="chart-layout-grid">
            <!-- Left: Plot Canvas -->
            <div id="plot-wrapper">
                {volcano_html}
            </div>

            <!-- Right: Interactive Controls & CSV Download Panel -->
            <div class="right-control-panel">
                <h3>🎛️ Splicing Threshold Controls</h3>
                
                <div class="slider-group slider-as">
                    <label for="psi-slider">
                        <span style="color: var(--blue-inc);">|ΔPSI| Cutoff (Effect Size):</span>
                        <span id="psi-val">{delta_psi_cutoff:.2f}</span>
                    </label>
                    <input type="range" id="psi-slider" min="0.01" max="0.5" step="0.01" value="{delta_psi_cutoff}">
                </div>

                <div class="slider-group slider-fdr">
                    <label for="fdr-slider">
                        <span style="color: #4F46E5;">rMATS FDR Cutoff (Significance):</span>
                        <span id="fdr-val">{as_fdr_cutoff:.3f}</span>
                    </label>
                    <input type="range" id="fdr-slider" min="0.001" max="0.1" step="0.005" value="{as_fdr_cutoff}">
                </div>

                <button class="btn-download-csv" id="download-csv-btn">
                    📥 Download Filtered Splicing Targets (.csv)
                </button>

                <!-- Color Classification Section (Outside Plot) -->
                <div>
                    <div class="legend-title-desc">Splicing Classification:</div>
                    <div class="legend-guide">
                        <div class="legend-item">
                            <span class="legend-dot" style="background: var(--blue-inc);"></span>
                            <span><b>Inclusion Favored:</b> ΔPSI ≥ +θ & FDR ≤ α</span>
                        </div>
                        <div class="legend-item">
                            <span class="legend-dot" style="background: var(--red-exc);"></span>
                            <span><b>Exclusion Favored:</b> ΔPSI ≤ -θ & FDR ≤ α</span>
                        </div>
                        <div class="legend-item">
                            <span class="legend-dot" style="background: #94A3B8;"></span>
                            <span><b>Non-Significant:</b> Invariant Background</span>
                        </div>
                    </div>
                </div>

                <!-- Generate GO / KEGG Button directly below Color Classification -->
                <button class="btn-generate-enrichment" id="generate-enrichment-btn">
                    🚀 Generate GO / KEGG
                </button>
            </div>
        </div>

        <!-- KPI Cards Positioned DIRECTLY BELOW Volcano Plot -->
        <div class="kpi-container">
            <div class="kpi-card">
                <div class="label">Total Splicing Events</div>
                <div class="value" id="kpi-total">{kpis['total']}</div>
            </div>
            <div class="kpi-card" style="border-top: 4px solid var(--purple-target);">
                <div class="label" style="color: var(--purple-target);">Significant Targets</div>
                <div class="value" style="color: var(--purple-target);" id="kpi-sig">{kpis['significant']}</div>
            </div>
            <div class="kpi-card" style="border-top: 4px solid var(--blue-inc);">
                <div class="label" style="color: var(--blue-inc);">Inclusion Favored</div>
                <div class="value" style="color: var(--blue-inc);" id="kpi-inc">{kpis['inclusion']}</div>
            </div>
            <div class="kpi-card" style="border-top: 4px solid var(--red-exc);">
                <div class="label" style="color: var(--red-exc);">Exclusion Favored</div>
                <div class="value" style="color: var(--red-exc);" id="kpi-exc">{kpis['exclusion']}</div>
            </div>
            <div class="kpi-card" style="border-top: 4px solid #94A3B8;">
                <div class="label" style="color: #64748B;">Non-Significant</div>
                <div class="value" style="color: #64748B;" id="kpi-nonsig">{kpis['non_significant']}</div>
            </div>
        </div>
    </div>

    <!-- GO Term Enrichment Section -->
    <div class="card" id="go-card-section">
        <div style="display: flex; justify-content: space-between; align-items: center; border-bottom: 2px solid var(--border-color); padding-bottom: 12px; margin-bottom: 12px;">
            <h2 id="go-title" style="margin: 0; border: none; padding: 0;">🧬 GO Term Biological Process Enrichment (Significant Splicing Targets)</h2>
            <button class="btn-download-csv" id="download-go-csv-btn">
                📥 Download GO Table (.csv)
            </button>
        </div>
        <div id="go-notice" style="font-size: 13px; color: #2563EB; font-weight: 600; margin-top: 8px;"></div>
        <div id="go-chart-container" style="margin-top: 16px;">
            <div id="go-placeholder" style="padding: 40px 20px; text-align: center; color: #64748B; background: #F8FAFC; border: 1px dashed #CBD5E1; border-radius: 8px;">
                <p style="margin: 0; font-size: 15px; font-weight: 600;">📊 GO Biological Process Enrichment Ready</p>
                <p style="margin: 6px 0 0; font-size: 13px; color: #94A3B8;">Click the <b>'🚀 Generate GO / KEGG'</b> button in the Threshold Controls panel above to compute and visualize enriched terms for active splicing targets.</p>
            </div>
            <div id="plotly-go-div" style="width: 100%; height: 380px; display: none;"></div>
        </div>
    </div>

    <!-- KEGG Pathway Enrichment Section -->
    <div class="card" id="kegg-card-section">
        <div style="display: flex; justify-content: space-between; align-items: center; border-bottom: 2px solid var(--border-color); padding-bottom: 12px; margin-bottom: 12px;">
            <h2 id="kegg-title" style="margin: 0; border: none; padding: 0;">🛤️ KEGG Pathway Enrichment (Significant Splicing Targets)</h2>
            <button class="btn-download-csv" id="download-kegg-csv-btn">
                📥 Download KEGG Table (.csv)
            </button>
        </div>
        <div id="kegg-notice" style="font-size: 13px; color: #7C3AED; font-weight: 600; margin-top: 8px;"></div>
        <div id="kegg-chart-container" style="margin-top: 16px;">
            <div id="kegg-placeholder" style="padding: 40px 20px; text-align: center; color: #64748B; background: #F8FAFC; border: 1px dashed #CBD5E1; border-radius: 8px;">
                <p style="margin: 0; font-size: 15px; font-weight: 600;">🛤️ KEGG Pathway Enrichment Ready</p>
                <p style="margin: 6px 0 0; font-size: 13px; color: #94A3B8;">Click the <b>'🚀 Generate GO / KEGG'</b> button in the Threshold Controls panel above to compute and visualize enriched pathways for active splicing targets.</p>
            </div>
            <div id="plotly-kegg-div" style="width: 100%; height: 380px; display: none;"></div>
        </div>
    </div>

    <!-- SECTION 1: NCBI / PubMed Automated Gene & Literature Explorer Card -->
    <div class="card" id="ncbi-pubmed-card">
        <h2 style="color: #0F172A; border-bottom: 2px solid #E2E8F0; padding-bottom: 10px; margin-top: 0;">📚 NCBI / PubMed Automated Gene & Literature Explorer</h2>
        
        <div style="margin-top: 16px; margin-bottom: 16px;">
            <label for="ncbi-gene-select" style="font-weight: 700; font-size: 14px; margin-right: 10px;">🔍 Select Splicing Target Gene:</label>
            <select id="ncbi-gene-select" style="padding: 8px 14px; border-radius: 8px; border: 1px solid #CBD5E1; font-weight: 700; font-size: 14px; background: #FFFFFF; cursor: pointer;">
            </select>
        </div>

        <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 20px;" id="ncbi-detail-grid">
            <div style="background-color: #FFFFFF; border: 1px solid #CBD5E1; border-radius: 12px; padding: 20px;">
                <h4 style="margin-top: 0; color: #0F172A; border-bottom: 2px solid #E2E8F0; padding-bottom: 8px;" id="ncbi-gene-title">🧬 NCBI Gene Details</h4>
                <p id="ncbi-official-name"><b>Official Name:</b> Loading...</p>
                <p id="ncbi-meta"><b>NCBI Gene ID:</b> <code>-</code> | <b>Map Location:</b> -</p>
                <p id="ncbi-summary-desc" style="font-size: 0.92rem; color: #475569; line-height: 1.55; margin-top: 10px;">Select a gene to inspect description.</p>
            </div>
            <div>
                <h4 style="margin-top: 0; color: #0F172A;" id="pubmed-list-title">📖 Recent PubMed Literature</h4>
                <div id="pubmed-papers-container" style="display: flex; flex-direction: column; gap: 10px;">
                </div>
            </div>
        </div>
    </div>

    <!-- SECTION 2: Event-Level Isoform Annotation & NMD Prediction Card -->
    <div class="card" id="isoform-annotation-card">
        <div style="display: flex; justify-content: space-between; align-items: center; border-bottom: 2px solid #E2E8F0; padding-bottom: 10px; margin-bottom: 12px;">
            <h2 style="color: #0F172A; margin: 0; border: none; padding: 0;">🧬 Event-Level Isoform Annotation & NMD Prediction Matrix</h2>
            <button class="btn-download-csv" id="download-isoform-csv-btn">
                📥 Download Filtered Isoforms (.csv)
            </button>
        </div>

        <div style="display: flex; justify-content: space-between; align-items: center; margin-top: 10px; margin-bottom: 12px; background: #F8FAFC; padding: 10px 16px; border-radius: 8px; border: 1px solid #E2E8F0;">
            <span style="font-size: 13px; color: #334155; font-weight: 600;">
                🔗 Synchronized Target: <b id="isoform-filter-gene" style="color: #1E40AF; font-size: 14px;">-</b> <span id="isoform-filter-count" style="font-size: 12px; color: #64748B; margin-left: 6px;"></span>
            </span>
            <button id="isoform-toggle-view-btn" style="background: #FFFFFF; border: 1px solid #CBD5E1; border-radius: 6px; padding: 5px 12px; font-size: 12px; color: #334155; cursor: pointer; font-weight: 700; box-shadow: 0 1px 2px rgba(0,0,0,0.05);">
                👁️ Show All Active Genes
            </button>
        </div>

        <div style="overflow-x: auto; overflow-y: auto; max-height: 420px; border: 1px solid #CBD5E1; border-radius: 8px;">
            <table class="data-table" id="isoform-table">
                <thead>
                    <tr>
                        <th>Gene</th>
                        <th>Impairment Risk (%)</th>
                        <th>Primary Dysfunction Cause</th>
                        <th>Event Type</th>
                        <th>Transcript ID</th>
                        <th>Coordinates</th>
                        <th>ΔPSI</th>
                        <th>Log₂FC</th>
                        <th>CDS Reading Frame</th>
                        <th>PTC & NMD Status</th>
                        <th>Protein Domain Impact</th>
                        <th>Localization Consequence</th>
                    </tr>
                </thead>
                <tbody id="isoform-table-body">
                    {isoform_table_rows_html}
                </tbody>
            </table>
        </div>
    </div>

    <!-- SECTION 3: Visual Exon-Intron Sashimi Structure Engine Card -->
    <div class="card" id="sashimi-card">
        <div style="display: flex; justify-content: space-between; align-items: center; border-bottom: 2px solid #E2E8F0; padding-bottom: 10px; margin-bottom: 12px; flex-wrap: wrap; gap: 10px;">
            <h2 style="color: #0F172A; margin: 0; border: none; padding: 0;">🧩 Visual Exon-Intron Structure & Sashimi Engine</h2>
            <div style="display: flex; align-items: center; gap: 12px; flex-wrap: wrap;">
                <div style="display: flex; align-items: center; gap: 6px;">
                    <label for="sashimi-gene-select" style="font-size: 13px; font-weight: 700; color: #334155;">Target Gene:</label>
                    <select id="sashimi-gene-select" style="padding: 6px 12px; border-radius: 6px; border: 1px solid #CBD5E1; font-size: 13px; font-weight: 700; background: #FFFFFF; cursor: pointer;"></select>
                </div>
                <div style="display: flex; align-items: center; gap: 6px;">
                    <label for="sashimi-event-select" id="sashimi-event-select-label" style="font-size: 13px; font-weight: 700; color: #334155;">Splicing Event:</label>
                    <select id="sashimi-event-select" style="padding: 6px 12px; border-radius: 6px; border: 1px solid #CBD5E1; font-size: 13px; font-weight: 600; background: #FFFFFF; cursor: pointer;"></select>
                </div>
            </div>
        </div>
        <div id="sashimi-chart-container" style="margin-top: 8px; width: 100%;">
            {sashimi_html}
        </div>
        <div style="text-align: center; margin-top: 12px; padding-top: 10px; border-top: 1px solid #E2E8F0; font-size: 15px; font-weight: 700; color: #1E293B;">
            📊 <span id="sashimi-title-text">{sashimi_title_text}</span>
        </div>
    </div>

    <!-- SECTION 4: Isoform-Specific RT-qPCR Primer Designer Card -->
    <div class="card" id="primer-card">
        <div style="display: flex; justify-content: space-between; align-items: center; border-bottom: 2px solid #E2E8F0; padding-bottom: 10px; margin-bottom: 16px; flex-wrap: wrap; gap: 8px;">
            <h2 style="color: #0F172A; margin: 0; border: none; padding: 0;">🧪 Isoform-Specific RT-qPCR Primer Designer</h2>
            <span style="font-size: 13px; color: #64748B; font-weight: 600;">🔗 Synchronized with Visual Exon-Intron Selection</span>
        </div>
        
        <!-- Target Gene Location & Exon Region Banner -->
        <div id="primer-location-banner" style="background: #F1F5F9; border-left: 4px solid #3B82F6; padding: 12px 16px; border-radius: 8px; margin-bottom: 16px; font-size: 14px; color: #334155;">
            <b>Selected Target:</b> <span id="primer-banner-gene" style="font-weight: 800; color: #1E40AF;">{first_primer_gene}</span> | 
            <b>Target Coordinates:</b> <code id="primer-banner-coords" style="color: #0F172A; font-weight: 700;">{first_coords}</code> | 
            <b>Targeted Splicing Event:</b> <span id="primer-banner-event" class="badge" style="background: #E0E7FF; color: #3730A3;">{first_event}</span>
        </div>

        <div style="overflow-x: auto; border: 1px solid #CBD5E1; border-radius: 8px;">
            <table class="data-table">
                <thead>
                    <tr>
                        <th>Target Gene</th>
                        <th>Target Isoform</th>
                        <th>Targeted Splice Junction / Region</th>
                        <th>Forward Primer (5'->3') & Tm</th>
                        <th>Reverse Primer (5'->3') & Tm</th>
                        <th>Amplicon Size</th>
                        <th>Primer Quality</th>
                    </tr>
                </thead>
                <tbody id="primer-table-body">
                    {primer_rows_html}
                </tbody>
            </table>
        </div>
    </div>

    <div class="footer">
        Generated automatically by GenSplice-Agent Pipeline • Light Mode Theme
    </div>

    <!-- Real-Time Threshold, Point Color Recalculation & Dynamic CSV Script -->
    <script>
        const rawGeneData = {raw_data_json};
        const genePathways = {gene_pathway_json};
        const primerRecords = {primer_records_json};
        const sashimiGeneData = {sashimi_data_json};
        const sashimiEventsMap = {sashimi_events_json};
        const isoformRecords = {isoform_records_json};
        const baseGoRecords = {go_records_json};
        const baseKeggRecords = {kegg_records_json};
        const currentOrganism = "{organism}";
        const liveGeneDataCache = new Map();
        let activeGeneQuery = "";
        let currentGoData = [];
        let currentKeggData = [];
        let isEnrichmentGenerated = false;
        
        const psiSlider = document.getElementById('psi-slider');
        const fdrSlider = document.getElementById('fdr-slider');
        const psiValLabel = document.getElementById('psi-val');
        const fdrValLabel = document.getElementById('fdr-val');
        const downloadBtn = document.getElementById('download-csv-btn');
        const generateBtn = document.getElementById('generate-enrichment-btn');
        const downloadGoBtn = document.getElementById('download-go-csv-btn');
        const downloadKeggBtn = document.getElementById('download-kegg-csv-btn');
        const downloadIsoformBtn = document.getElementById('download-isoform-csv-btn');
        const toggleIsoformBtn = document.getElementById('isoform-toggle-view-btn');
        const ncbiSelect = document.getElementById('ncbi-gene-select');
        const sashimiGeneSelect = document.getElementById('sashimi-gene-select');
        const sashimiEventSelect = document.getElementById('sashimi-event-select');
        const sashimiEventSelectLabel = document.getElementById('sashimi-event-select-label');

        let currentFilteredGenes = [];

        // NCBI E-utilities Rate Limiter & Resilience Layer (NCBI max 3 req/sec per IP)
        let lastNcbiRequestTime = 0;
        const MIN_NCBI_INTERVAL_MS = 380; // Strictly throttled to ~2.6 req/sec to prevent HTTP 429
        const NCBI_TOOL_PARAMS = "tool=GenSpliceAgent&email=gensplice_report%40ncbi.nlm.nih.gov";

        async function safeNcbiFetch(url, maxRetries = 2) {{
            const separator = url.includes("?") ? "&" : "?";
            const fullUrl = url.includes("tool=") ? url : `${{url}}${{separator}}${{NCBI_TOOL_PARAMS}}`;

            for (let attempt = 0; attempt <= maxRetries; attempt++) {{
                const now = Date.now();
                const waitMs = Math.max(0, MIN_NCBI_INTERVAL_MS - (now - lastNcbiRequestTime));
                if (waitMs > 0) {{
                    await new Promise(res => setTimeout(res, waitMs));
                }}
                lastNcbiRequestTime = Date.now();

                try {{
                    const resp = await fetch(fullUrl);
                    if (resp.status === 429) {{
                        console.warn(`NCBI Rate limit (429) on attempt ${{attempt + 1}}. Backing off...`);
                        if (attempt < maxRetries) {{
                            await new Promise(res => setTimeout(res, 800 * (attempt + 1)));
                            continue;
                        }}
                        return {{ ok: false, status: 429, data: null }};
                    }}
                    if (!resp.ok) {{
                        console.warn(`NCBI fetch HTTP ${{resp.status}} on attempt ${{attempt + 1}}`);
                        if (attempt < maxRetries) {{
                            await new Promise(res => setTimeout(res, 500 * (attempt + 1)));
                            continue;
                        }}
                        return {{ ok: false, status: resp.status, data: null }};
                    }}
                    const data = await resp.json();
                    return {{ ok: true, status: resp.status, data: data }};
                }} catch (err) {{
                    console.warn(`NCBI fetch network/parse error on attempt ${{attempt + 1}}:`, err);
                    if (attempt < maxRetries) {{
                        await new Promise(res => setTimeout(res, 600 * (attempt + 1)));
                        continue;
                    }}
                    return {{ ok: false, status: 0, error: err }};
                }}
            }}
            return {{ ok: false, status: 0, data: null }};
        }}

        async function fetchNcbiGene(upper) {{
            try {{
                const orgTerm = encodeURIComponent(currentOrganism || "Homo sapiens");
                let sUrl = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=${{encodeURIComponent(upper)}}[Gene+Name]+AND+${{orgTerm}}[Organism]&retmode=json`;
                let sRes = await safeNcbiFetch(sUrl);
                let idList = (sRes.ok && sRes.data?.esearchresult?.idlist) ? sRes.data.esearchresult.idlist : [];

                if (!idList.length) {{
                    sUrl = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=${{encodeURIComponent(upper)}}[Gene+Name]&retmode=json`;
                    sRes = await safeNcbiFetch(sUrl);
                    idList = (sRes.ok && sRes.data?.esearchresult?.idlist) ? sRes.data.esearchresult.idlist : [];
                }}

                if (idList.length > 0) {{
                    const geneId = idList[0];
                    const sumUrl = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=${{geneId}}&retmode=json`;
                    const sumRes = await safeNcbiFetch(sumUrl);
                    if (sumRes.ok && sumRes.data?.result?.[geneId]) {{
                        const ginfo = sumRes.data.result[geneId];
                        return {{
                            ncbi_id: geneId,
                            official_name: ginfo.description || `${{upper}} (${{currentOrganism || 'Homo sapiens'}})`,
                            chromosome: ginfo.maplocation || "N/A",
                            summary: ginfo.summary || `Official NCBI Gene record for ${{upper}}.`
                        }};
                    }}
                }}
            }} catch (err) {{
                console.warn("NCBI Gene live fetch warning:", err);
            }}
            return {{
                ncbi_id: "N/A",
                official_name: `${{upper}} (${{currentOrganism || 'Homo sapiens'}})`,
                chromosome: "N/A",
                summary: `Official NCBI Gene record for ${{upper}}. Click link above to view details directly on NCBI.`
            }};
        }}

        async function fetchPubMedLiterature(upper) {{
            try {{
                const pSearchUrl = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=${{encodeURIComponent(upper)}}&retmax=3&sort=pub_date&retmode=json`;
                const pRes = await safeNcbiFetch(pSearchUrl);
                if (!pRes.ok) {{
                    return {{ papers: [], totalCount: 0, error: true, rateLimited: pRes.status === 429 }};
                }}
                const pData = pRes.data || {{}};
                const pmidList = pData.esearchresult?.idlist || [];
                const totalCount = parseInt(pData.esearchresult?.count || "0", 10);

                if (pmidList.length > 0) {{
                    const pSumUrl = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=${{pmidList.join(",")}}&retmode=json`;
                    const pSumRes = await safeNcbiFetch(pSumUrl);
                    if (!pSumRes.ok) {{
                        return {{ papers: [], totalCount: totalCount, error: true, rateLimited: pSumRes.status === 429 }};
                    }}
                    const pSumData = pSumRes.data || {{}};
                    const pMap = pSumData.result || {{}};
                    const papers = [];

                    pmidList.forEach(pmid => {{
                        const p = pMap[pmid] || {{}};
                        const authorsList = p.authors || [];
                        let authorStr = "Unknown";
                        if (authorsList.length === 1) {{
                            authorStr = authorsList[0].name || "Unknown";
                        }} else if (authorsList.length > 1) {{
                            authorStr = `${{authorsList[0].name || "Unknown"}} et al.`;
                        }}

                        const pubDate = (p.pubdate || p.sortpubdate || "").split(" ")[0] || "Recent";
                        const journal = p.source || "Journal";
                        const title = p.title || `Research on ${{upper}} gene function`;

                        papers.push({{
                            pmid: pmid,
                            title: title,
                            authors: authorStr,
                            journal: journal,
                            pub_date: pubDate,
                            url: `https://pubmed.ncbi.nlm.nih.gov/${{pmid}}/`
                        }});
                    }});

                    return {{ papers: papers, totalCount: totalCount, error: false }};
                }}
                return {{ papers: [], totalCount: totalCount, error: false }};
            }} catch (err) {{
                console.warn("PubMed live fetch warning:", err);
                return {{ papers: [], totalCount: 0, error: true }};
            }}
        }}

        function displayNcbiData(symbol, data) {{
            const ncbi = (data && data.ncbi) ? data.ncbi : {{}};
            const geneId = (ncbi.ncbi_id && ncbi.ncbi_id !== 'N/A') ? ncbi.ncbi_id : null;
            const ncbiLink = geneId ? `https://www.ncbi.nlm.nih.gov/gene/${{geneId}}` : `https://www.ncbi.nlm.nih.gov/gene/?term=${{encodeURIComponent(symbol)}}`;
            const pubmedSearchUrl = `https://pubmed.ncbi.nlm.nih.gov/?term=${{encodeURIComponent(symbol)}}`;

            const geneTitle = document.getElementById('ncbi-gene-title');
            if (geneTitle) {{
                geneTitle.innerHTML = `🧬 NCBI Gene Details: <code style="color:#2563EB;">${{symbol}}</code> <a href="${{ncbiLink}}" target="_blank" rel="noopener noreferrer" style="font-size:12px; margin-left:8px; color:#2563EB; text-decoration:none; font-weight:700;">🔗 Open in NCBI ↗</a>`;
            }}
            const officialName = document.getElementById('ncbi-official-name');
            if (officialName) {{
                officialName.innerHTML = `<b>Official Name:</b> ${{ncbi.official_name || (symbol + ' (' + (currentOrganism || 'Homo sapiens') + ')')}}`;
            }}
            const meta = document.getElementById('ncbi-meta');
            if (meta) {{
                const idDisplay = geneId ? `<a href="${{ncbiLink}}" target="_blank" rel="noopener noreferrer" style="color:#2563EB; font-weight:700;">${{geneId}}</a>` : `<a href="${{ncbiLink}}" target="_blank" rel="noopener noreferrer" style="color:#2563EB; font-weight:700;">Search NCBI</a>`;
                meta.innerHTML = `<b>NCBI Gene ID:</b> ${{idDisplay}} | <b>Map Location:</b> ${{ncbi.chromosome || 'N/A'}}`;
            }}
            const summaryDesc = document.getElementById('ncbi-summary-desc');
            if (summaryDesc) {{
                summaryDesc.textContent = ncbi.summary || `Official NCBI summary for ${{symbol}}.`;
            }}

            const pubmedInfo = (data && data.pubmed) ? data.pubmed : {{ papers: [], totalCount: 0, error: false }};
            const pubmedList = pubmedInfo.papers || [];
            const pubContainer = document.getElementById('pubmed-papers-container');
            const pubTitle = document.getElementById('pubmed-list-title');

            if (pubTitle) {{
                pubTitle.innerHTML = `📖 Recent PubMed Literature <span style="font-size:12px; font-weight:600; color:#2563EB;">(Live NIH API)</span>`;
            }}

            if (pubContainer) {{
                if (pubmedList.length > 0) {{
                    let pubHtml = '';
                    pubmedList.forEach(paper => {{
                        pubHtml += `
                        <div style="background-color: #F8FAFC; border-left: 4px solid #3B82F6; border-radius: 8px; padding: 12px 14px; border: 1px solid #E2E8F0; margin-bottom: 8px;">
                            <a href="${{paper.url}}" target="_blank" rel="noopener noreferrer" style="text-decoration: none; font-weight: 700; color: #1D4ED8; font-size: 0.95rem; line-height: 1.45; display: block;">🔗 ${{paper.title}}</a>
                            <div style="font-size: 0.82rem; color: #64748B; margin-top: 6px; display: flex; flex-wrap: wrap; gap: 8px; align-items: center;">
                                <span><b>Authors:</b> ${{paper.authors}}</span>
                                <span>•</span>
                                <span><b>Journal:</b> ${{paper.journal}} (${{paper.pub_date}})</span>
                                <span>•</span>
                                <span><b>PMID:</b> <a href="${{paper.url}}" target="_blank" rel="noopener noreferrer" style="color:#2563EB; font-weight:700;">${{paper.pmid}}</a></span>
                            </div>
                        </div>`;
                    }});

                    const totalMsg = pubmedInfo.totalCount > 0 ? `all ${{pubmedInfo.totalCount.toLocaleString()}} PubMed papers` : `all papers`;
                    pubHtml += `
                    <div style="text-align: right; margin-top: 6px;">
                        <a href="${{pubmedSearchUrl}}" target="_blank" rel="noopener noreferrer" style="font-size: 12px; font-weight: 700; color: #2563EB; text-decoration: none;">
                            🔍 View ${{totalMsg}} for "${{symbol}}" ↗
                        </a>
                    </div>`;
                    pubContainer.innerHTML = pubHtml;
                }} else if (pubmedInfo.error) {{
                    const retryBtnHtml = `<button onclick="renderNcbiGeneDetails('${{symbol}}')" style="margin-top:10px; padding:6px 14px; background:#2563EB; color:#FFF; font-weight:700; border:none; border-radius:6px; cursor:pointer; font-size:12px; box-shadow:0 1px 2px rgba(0,0,0,0.1);">🔄 Retry NIH Query</button>`;
                    const rateLimitMsg = pubmedInfo.rateLimited 
                        ? "NCBI API rate limit reached (too many requests in a short period). Please wait a moment and click Retry."
                        : "NIH PubMed connection temporarily unavailable or restricted by network/browser policies.";
                    pubContainer.innerHTML = `
                    <div style="padding:16px; background:#FEF2F2; border-radius:8px; border:1px solid #FECDD3; text-align:center;">
                        <p style="color:#DC2626; font-size:13px; font-weight:700; margin:0 0 6px;">⚠️ ${{rateLimitMsg}}</p>
                        <p style="color:#64748B; font-size:12px; margin:0 0 8px;">You can retry the query or inspect PubMed publications directly in a new tab:</p>
                        <div style="display:flex; justify-content:center; gap:10px; margin-top:8px;">
                            ${{retryBtnHtml}}
                            <a href="${{pubmedSearchUrl}}" target="_blank" rel="noopener noreferrer" style="display:inline-block; margin-top:10px; padding:6px 14px; background:#475569; color:#FFF; font-weight:700; border-radius:6px; text-decoration:none; font-size:12px; line-height:20px;">🔍 Open in PubMed ↗</a>
                        </div>
                    </div>`;
                }} else {{
                    pubContainer.innerHTML = `
                    <div style="padding:20px; background:#F8FAFC; border-radius:8px; border:1px dashed #CBD5E1; text-align:center;">
                        <p style="color:#64748B; margin:0 0 10px; font-size:13px;">No recent indexed PubMed papers found directly matching <b>${{symbol}}</b>.</p>
                        <a href="${{pubmedSearchUrl}}" target="_blank" rel="noopener noreferrer" style="display:inline-block; padding:8px 14px; background:#3B82F6; color:#FFF; font-weight:700; border-radius:6px; text-decoration:none; font-size:13px;">🔍 Search PubMed directly for "${{symbol}}" ↗</a>
                    </div>`;
                }}
            }}
        }}

        async function renderNcbiGeneDetails(symbol) {{
            if (!symbol) {{
                const geneTitle = document.getElementById('ncbi-gene-title');
                if (geneTitle) geneTitle.innerHTML = `🧬 NCBI Gene Details: <span style="color:#64748B;">No target selected</span>`;
                const officialName = document.getElementById('ncbi-official-name');
                if (officialName) officialName.innerHTML = `<b>Official Name:</b> <i>No splicing targets meet current thresholds.</i>`;
                const meta = document.getElementById('ncbi-meta');
                if (meta) meta.innerHTML = `<b>NCBI Gene ID:</b> - | <b>Map Location:</b> -`;
                const summaryDesc = document.getElementById('ncbi-summary-desc');
                if (summaryDesc) summaryDesc.textContent = `No active significant splicing targets found under current threshold cutoffs. Adjust ΔPSI or FDR sliders to view targets.`;
                const pubElem = document.getElementById('pubmed-papers-container');
                if (pubElem) pubElem.innerHTML = `<p style="color:#64748B;"><i>No target gene selected.</i></p>`;
                return;
            }}

            const upper = symbol.trim().toUpperCase();
            activeGeneQuery = upper;

            // 1. Instant display from in-memory session cache if already fetched in this session
            if (liveGeneDataCache.has(upper)) {{
                displayNcbiData(upper, liveGeneDataCache.get(upper));
                return;
            }}

            // 2. Active loading UI state while fetching live
            const geneTitle = document.getElementById('ncbi-gene-title');
            if (geneTitle) {{
                geneTitle.innerHTML = `🧬 NCBI Gene Details: <code style="color:#2563EB;">${{upper}}</code> <span style="font-size:12px; color:#3B82F6; font-weight:600;">(Fetching live from NIH NCBI...)</span>`;
            }}
            const officialName = document.getElementById('ncbi-official-name');
            if (officialName) {{
                officialName.innerHTML = `<b>Official Name:</b> <i>Connecting to NIH NCBI E-utilities API...</i>`;
            }}
            const meta = document.getElementById('ncbi-meta');
            if (meta) {{
                meta.innerHTML = `<b>NCBI Gene ID:</b> <code>Loading...</code> | <b>Map Location:</b> Loading...`;
            }}
            const summaryDesc = document.getElementById('ncbi-summary-desc');
            if (summaryDesc) {{
                summaryDesc.innerHTML = `<span style="color:#64748B; font-style:italic;">🔄 Fetching official NCBI gene records & chromosomal location in real-time...</span>`;
            }}
            const pubElem = document.getElementById('pubmed-papers-container');
            if (pubElem) {{
                pubElem.innerHTML = `
                <div style="padding:22px; background:#F8FAFC; border-radius:8px; border:1px dashed #CBD5E1; text-align:center;">
                    <div style="font-size:14px; font-weight:700; color:#1E40AF; margin-bottom:4px;">🔄 Querying PubMed in real-time for ${{upper}}...</div>
                    <p style="color:#64748B; margin:0; font-size:12px;">Fetching latest peer-reviewed literature via NIH PubMed E-utilities API.</p>
                </div>`;
            }}

            // 3. Sequenced live API query to NCBI Gene and PubMed (no burst collision)
            const ncbiResult = await fetchNcbiGene(upper);
            if (activeGeneQuery !== upper) return; // Prevent race condition if user switched gene
            const pubmedResult = await fetchPubMedLiterature(upper);
            if (activeGeneQuery !== upper) return; // Prevent race condition if user switched gene

            const combinedData = {{
                ncbi: ncbiResult,
                pubmed: pubmedResult
            }};

            // Only cache in session memory if query completed without error!
            if (!pubmedResult.error) {{
                liveGeneDataCache.set(upper, combinedData);
            }}
            displayNcbiData(upper, combinedData);
        }}

        function updateNcbiDropdown(splicingGeneList) {{
            if (!ncbiSelect && !sashimiGeneSelect) return;
            const uniqueGenes = [...new Set((splicingGeneList || []).filter(Boolean))];
            
            const prevVal = ncbiSelect ? ncbiSelect.value : (sashimiGeneSelect ? sashimiGeneSelect.value : "");
            
            [ncbiSelect, sashimiGeneSelect].forEach(sel => {{
                if (!sel) return;
                sel.innerHTML = '';
                if (!uniqueGenes.length) {{
                    const opt = document.createElement('option');
                    opt.value = "";
                    opt.textContent = "No splicing targets meeting thresholds";
                    sel.appendChild(opt);
                }} else {{
                    uniqueGenes.forEach(g => {{
                        const opt = document.createElement('option');
                        opt.value = g;
                        opt.textContent = g;
                        sel.appendChild(opt);
                    }});
                }}
            }});

            if (!uniqueGenes.length) {{
                handleMasterGeneSelection("");
                return;
            }}

            const chosen = (prevVal && uniqueGenes.includes(prevVal)) ? prevVal : uniqueGenes[0];
            if (ncbiSelect) ncbiSelect.value = chosen;
            if (sashimiGeneSelect) sashimiGeneSelect.value = chosen;
            if (chosen !== activeGeneQuery) {{
                handleMasterGeneSelection(chosen);
            }}
        }}

        function renderSashimiPlot(geneSymbol, eventIndex = 0) {{
            const sashimiDiv = document.getElementById('plotly-sashimi-div');
            if (!sashimiDiv || !window.Plotly) return;

            if (!geneSymbol) {{
                Plotly.react(sashimiDiv, [], {{
                    title: "<b>No active splicing target gene selected meeting current thresholds</b>",
                    paper_bgcolor: "#FFFFFF", plot_bgcolor: "#F8FAFC", height: 260
                }});
                const titleElem = document.getElementById('sashimi-title-text');
                if (titleElem) {{
                    titleElem.innerHTML = `<b>No Splicing Target Selected</b>`;
                }}
                if (sashimiEventSelect) sashimiEventSelect.style.display = 'none';
                if (sashimiEventSelectLabel) sashimiEventSelectLabel.style.display = 'none';
                return;
            }}

            const upper = (geneSymbol || '').trim().toUpperCase();
            const events = (sashimiEventsMap && sashimiEventsMap[upper]) || [];

            if (sashimiEventSelect) {{
                sashimiEventSelect.innerHTML = '';
                if (events.length > 0) {{
                    events.forEach((ev, idx) => {{
                        const opt = document.createElement('option');
                        opt.value = idx;
                        const sign = ev.delta_psi >= 0 ? '+' : '';
                        const coordShort = (ev.coordinates || '').split(':')[1] || ev.coordinates || 'N/A';
                        opt.textContent = `Event ${{idx + 1}}: ${{ev.event_type}} (${{coordShort}}) [ΔPSI=${{sign}}${{parseFloat(ev.delta_psi).toFixed(2)}}]`;
                        sashimiEventSelect.appendChild(opt);
                    }});
                    if (events.length > 1) {{
                        const allOpt = document.createElement('option');
                        allOpt.value = -1;
                        allOpt.textContent = `All Events of ${{geneSymbol}} (${{events.length}} events)`;
                        sashimiEventSelect.appendChild(allOpt);
                    }}
                    sashimiEventSelect.value = eventIndex;
                    sashimiEventSelect.disabled = false;
                }} else {{
                    const opt = document.createElement('option');
                    opt.value = 0;
                    opt.textContent = `Event 1: SE (Primary)`;
                    sashimiEventSelect.appendChild(opt);
                    sashimiEventSelect.disabled = true;
                }}
            }}

            const evIdx = parseInt(eventIndex !== undefined && eventIndex !== null ? eventIndex : 0, 10);
            let info = (evIdx >= 0 && events.length > evIdx ? events[evIdx] : (events.length > 0 ? events[0] : null)) || 
                       (sashimiGeneData && (sashimiGeneData[upper] || sashimiGeneData[geneSymbol])) || 
                       (rawGeneData && rawGeneData.find(g => (g.geneSymbol || '').toUpperCase() === upper));

            if (!info) {{
                info = {{ gene_symbol: geneSymbol, event_type: "SE", coordinates: "N/A", delta_psi: 0.0, inc_counts: 0, exc_counts: 0 }};
            }}

            const symbol = info.gene_symbol || geneSymbol;
            const eventType = (info.event_type || "SE").toUpperCase();
            const deltaPsi = parseFloat(info.delta_psi !== undefined && info.delta_psi !== null ? info.delta_psi : 0.0);

            // Update title text below the plot
            const titleElem = document.getElementById('sashimi-title-text');
            if (titleElem) {{
                const sign = deltaPsi >= 0 ? '+' : '';
                const evNumText = evIdx === -1
                    ? ` [All ${{events.length}} Events]` 
                    : (events.length > 1 ? ` [Event ${{evIdx + 1}} of ${{events.length}}]` : '');
                titleElem.innerHTML = `<b>${{symbol}} Exon Structure & Splicing Sashimi Plot${{evNumText}}</b> (${{eventType}} | ΔPSI = ${{sign}}${{deltaPsi.toFixed(2)}})`;
            }}

            const inc_counts = (info.inc_counts !== undefined && info.inc_counts !== null) ? info.inc_counts : 0;
            const exc_counts = (info.exc_counts !== undefined && info.exc_counts !== null) ? info.exc_counts : 0;

            let chrom = info.chrom || "chr1";
            let ex1 = info.ex1, ex2 = info.ex2, ex3 = info.ex3;

            if (!ex1 || !ex2 || !ex3) {{
                const coords = info.coordinates || "";
                const mChr = coords.match(/(chr[0-9XYM]+)/i);
                if (mChr) chrom = mChr[1].toLowerCase();
                
                const nums = (coords.match(/\\d+/g) || []).map(Number);
                if (nums.length >= 6) {{
                    ex1 = [nums[0], nums[1]];
                    ex2 = [nums[2], nums[3]];
                    ex3 = [nums[4], nums[5]];
                }} else if (nums.length >= 4) {{
                    const p1 = nums[0], p2 = nums[1], p3 = nums[2], p4 = nums[3];
                    const span = Math.max(400, p4 - p1);
                    ex1 = [p1, p2];
                    ex2 = [p2 + Math.round(span * 0.25), p2 + Math.round(span * 0.45)];
                    ex3 = [p3, p4];
                }} else if (nums.length >= 2) {{
                    const s = nums[0], e = nums[1], diff = Math.max(150, e - s);
                    ex1 = [s - diff * 3, s - diff * 2];
                    ex2 = [s, e];
                    ex3 = [e + diff * 2, e + diff * 3];
                }} else {{
                    chrom = info.chrom || "chr1";
                    ex1 = [1000, 1200];
                    ex2 = [1800, 2000];
                    ex3 = [2600, 2800];
                }}
            }}

            const color_constitutive = "#2563EB";
            const color_alt = "#EF4444";
            const color_retained = "#F59E0B";
            const color_mxe_b = "#9333EA";
            const color_inc_arc = "#10B981";
            const color_exc_arc = "#8B5CF6";

            const shapes = [];
            const traces = [];

            function addBox(x0, x1, y_c, h, col, name, tip) {{
                shapes.push({{
                    type: "rect",
                    x0: x0, x1: x1,
                    y0: y_c - h/2, y1: y_c + h/2,
                    fillcolor: col,
                    line: {{ color: "#0F172A", width: 1.5 }},
                    layer: "above"
                }});
                traces.push({{
                    x: [(x0 + x1) / 2], y: [y_c],
                    mode: "markers", marker: {{ size: 1, color: "rgba(0,0,0,0)" }},
                    name: name,
                    hovertemplate: (tip || `<b>${{name}}</b><br>Coords: ${{chrom}}:${{x0.toLocaleString()}}-${{x1.toLocaleString()}} (${{(x1-x0).toLocaleString()}} bp)`) + "<extra></extra>"
                }});
            }}

            let x_min = ex1[0] - 250;
            let x_max = ex3[1] + 250;
            let tick_inc = `Inclusion Isoform (${{eventType}})`;
            let tick_exc = `Exclusion Isoform`;

            if (eventType === "RI") {{
                x_min = ex1[0] - 250;
                x_max = ex2[1] + 250;
                tick_inc = "Inclusion (Retained Intron)";
                tick_exc = "Exclusion (Spliced Intron)";

                traces.push({{ x: [ex1[0], ex2[1]], y: [2, 2], mode: "lines", line: {{ color: "#94A3B8", width: 2, dash: "dash" }}, hoverinfo: "skip", showlegend: false }});
                traces.push({{ x: [ex1[0], ex2[1]], y: [-2, -2], mode: "lines", line: {{ color: "#94A3B8", width: 2, dash: "dash" }}, hoverinfo: "skip", showlegend: false }});

                addBox(ex1[0], ex1[1], 2, 0.8, color_constitutive, "Exon 1 (Upstream)");
                addBox(ex1[1], ex2[0], 2, 0.5, color_retained, "Retained Intron Block", `<b>Retained Intron Block</b><br>Span: ${{chrom}}:${{ex1[1].toLocaleString()}}-${{ex2[0].toLocaleString()}}<br>Status: Intron retained in mRNA`);
                addBox(ex2[0], ex2[1], 2, 0.8, color_constitutive, "Exon 2 (Downstream)");

                traces.push({{
                    x: [(ex1[1] + ex2[0]) / 2], y: [2.7],
                    mode: "text", text: [`Retained Intron Reads: ${{inc_counts}}`],
                    textposition: "top center", textfont: {{ color: "#D97706", size: 12, family: "Arial Black" }},
                    name: "Intron Retention Signal"
                }});

                addBox(ex1[0], ex1[1], -2, 0.8, color_constitutive, "Exon 1 (Upstream)");
                addBox(ex2[0], ex2[1], -2, 0.8, color_constitutive, "Exon 2 (Downstream)");

                const mid_ri = (ex1[1] + ex2[0]) / 2;
                traces.push({{
                    x: [ex1[1], mid_ri, ex2[0]], y: [-2.4, -3.5, -2.4],
                    mode: "lines+text", line: {{ color: color_exc_arc, width: 3.5 }},
                    text: ["", `Spliced Intron Reads: ${{exc_counts}}`, ""],
                    textposition: "bottom center", textfont: {{ color: color_exc_arc, size: 11, family: "Arial Bold" }},
                    name: "Spliced Junction Arc"
                }});

            }} else if (eventType === "MXE") {{
                const ex2a = ex2;
                const gap = Math.round((ex3[0] - ex2a[1]) / 3);
                const ex2b_len = ex2a[1] - ex2a[0];
                const ex2b = [ex2a[1] + gap, ex2a[1] + gap + ex2b_len];
                tick_inc = "Isoform 1 (Exon 2A)";
                tick_exc = "Isoform 2 (Exon 2B)";

                traces.push({{ x: [ex1[0], ex3[1]], y: [2, 2], mode: "lines", line: {{ color: "#94A3B8", width: 2, dash: "dash" }}, hoverinfo: "skip", showlegend: false }});
                traces.push({{ x: [ex1[0], ex3[1]], y: [-2, -2], mode: "lines", line: {{ color: "#94A3B8", width: 2, dash: "dash" }}, hoverinfo: "skip", showlegend: false }});

                addBox(ex1[0], ex1[1], 2, 0.8, color_constitutive, "Exon 1");
                addBox(ex2a[0], ex2a[1], 2, 0.8, color_alt, "Exon 2A (Isoform 1)");
                addBox(ex3[0], ex3[1], 2, 0.8, color_constitutive, "Exon 3");

                const mid_1a = (ex1[1] + ex2a[0]) / 2;
                traces.push({{
                    x: [ex1[1], mid_1a, ex2a[0]], y: [2.4, 3.2, 2.4],
                    mode: "lines+text", line: {{ color: color_inc_arc, width: 3 }},
                    text: ["", `Exon 2A Reads: ${{inc_counts}}`, ""],
                    textposition: "top center", name: "Junction 1 -> 2A"
                }});
                const mid_a3 = (ex2a[1] + ex3[0]) / 2;
                traces.push({{
                    x: [ex2a[1], mid_a3, ex3[0]], y: [2.4, 3.2, 2.4],
                    mode: "lines", line: {{ color: color_inc_arc, width: 3 }},
                    name: "Junction 2A -> 3"
                }});

                addBox(ex1[0], ex1[1], -2, 0.8, color_constitutive, "Exon 1");
                addBox(ex2b[0], ex2b[1], -2, 0.8, color_mxe_b, "Exon 2B (Isoform 2)");
                addBox(ex3[0], ex3[1], -2, 0.8, color_constitutive, "Exon 3");

                const mid_1b = (ex1[1] + ex2b[0]) / 2;
                traces.push({{
                    x: [ex1[1], mid_1b, ex2b[0]], y: [-2.4, -3.4, -2.4],
                    mode: "lines+text", line: {{ color: color_exc_arc, width: 3 }},
                    text: ["", `Exon 2B Reads: ${{exc_counts}}`, ""],
                    textposition: "bottom center", name: "Junction 1 -> 2B"
                }});
                const mid_b3 = (ex2b[1] + ex3[0]) / 2;
                traces.push({{
                    x: [ex2b[1], mid_b3, ex3[0]], y: [-2.4, -3.4, -2.4],
                    mode: "lines", line: {{ color: color_exc_arc, width: 3 }},
                    name: "Junction 2B -> 3"
                }});

            }} else if (eventType === "A5SS") {{
                const ext_len = Math.max(80, Math.round((ex1[1] - ex1[0]) * 0.6));
                const ex1_long = [ex1[0], ex1[1] + ext_len];
                const ex2_down = ex3;
                x_max = ex2_down[1] + 250;
                tick_inc = "Long Isoform (Distal 5' Donor)";
                tick_exc = "Short Isoform (Proximal 5' Donor)";

                traces.push({{ x: [ex1[0], ex2_down[1]], y: [2, 2], mode: "lines", line: {{ color: "#94A3B8", width: 2, dash: "dash" }}, hoverinfo: "skip", showlegend: false }});
                traces.push({{ x: [ex1[0], ex2_down[1]], y: [-2, -2], mode: "lines", line: {{ color: "#94A3B8", width: 2, dash: "dash" }}, hoverinfo: "skip", showlegend: false }});

                addBox(ex1[0], ex1[1], 2, 0.8, color_constitutive, "Core Exon 1");
                addBox(ex1[1], ex1_long[1], 2, 0.8, color_alt, "Alt 5' Extension", "<b>Alternative 5' Extension</b><br>Distal splice donor site");
                addBox(ex2_down[0], ex2_down[1], 2, 0.8, color_constitutive, "Exon 2 (Downstream)");

                const mid_long = (ex1_long[1] + ex2_down[0]) / 2;
                traces.push({{
                    x: [ex1_long[1], mid_long, ex2_down[0]], y: [2.4, 3.3, 2.4],
                    mode: "lines+text", line: {{ color: color_inc_arc, width: 3 }},
                    text: ["", `Long 5' Reads: ${{inc_counts}}`, ""],
                    textposition: "top center", name: "Distal 5' Splice Junction"
                }});

                addBox(ex1[0], ex1[1], -2, 0.8, color_constitutive, "Core Exon 1");
                addBox(ex2_down[0], ex2_down[1], -2, 0.8, color_constitutive, "Exon 2 (Downstream)");

                const mid_short = (ex1[1] + ex2_down[0]) / 2;
                traces.push({{
                    x: [ex1[1], mid_short, ex2_down[0]], y: [-2.4, -3.4, -2.4],
                    mode: "lines+text", line: {{ color: color_exc_arc, width: 3 }},
                    text: ["", `Short 5' Reads: ${{exc_counts}}`, ""],
                    textposition: "bottom center", name: "Proximal 5' Splice Junction"
                }});

            }} else if (eventType === "A3SS") {{
                const ext_len = Math.max(80, Math.round((ex3[1] - ex3[0]) * 0.5));
                const ex2_long = [ex3[0] - ext_len, ex3[1]];
                tick_inc = "Long Isoform (Proximal 3' Acceptor)";
                tick_exc = "Short Isoform (Distal 3' Acceptor)";

                traces.push({{ x: [ex1[0], ex3[1]], y: [2, 2], mode: "lines", line: {{ color: "#94A3B8", width: 2, dash: "dash" }}, hoverinfo: "skip", showlegend: false }});
                traces.push({{ x: [ex1[0], ex3[1]], y: [-2, -2], mode: "lines", line: {{ color: "#94A3B8", width: 2, dash: "dash" }}, hoverinfo: "skip", showlegend: false }});

                addBox(ex1[0], ex1[1], 2, 0.8, color_constitutive, "Exon 1 (Upstream)");
                addBox(ex2_long[0], ex3[0], 2, 0.8, color_alt, "Alt 3' Extension", "<b>Alternative 3' Extension</b><br>Proximal splice acceptor site");
                addBox(ex3[0], ex3[1], 2, 0.8, color_constitutive, "Core Exon 2");

                const mid_long3 = (ex1[1] + ex2_long[0]) / 2;
                traces.push({{
                    x: [ex1[1], mid_long3, ex2_long[0]], y: [2.4, 3.3, 2.4],
                    mode: "lines+text", line: {{ color: color_inc_arc, width: 3 }},
                    text: ["", `Long 3' Reads: ${{inc_counts}}`, ""],
                    textposition: "top center", name: "Proximal 3' Splice Junction"
                }});

                addBox(ex1[0], ex1[1], -2, 0.8, color_constitutive, "Exon 1 (Upstream)");
                addBox(ex3[0], ex3[1], -2, 0.8, color_constitutive, "Core Exon 2");

                const mid_short3 = (ex1[1] + ex3[0]) / 2;
                traces.push({{
                    x: [ex1[1], mid_short3, ex3[0]], y: [-2.4, -3.4, -2.4],
                    mode: "lines+text", line: {{ color: color_exc_arc, width: 3 }},
                    text: ["", `Short 3' Reads: ${{exc_counts}}`, ""],
                    textposition: "bottom center", name: "Distal 3' Splice Junction"
                }});

            }} else {{
                // Default: SE (Skipped Exon / Cassette Exon)
                traces.push({{ x: [ex1[0], ex3[1]], y: [2, 2], mode: "lines", line: {{ color: "#94A3B8", width: 2, dash: "dash" }}, hoverinfo: "skip", showlegend: false }});
                traces.push({{ x: [ex1[0], ex3[1]], y: [-2, -2], mode: "lines", line: {{ color: "#94A3B8", width: 2, dash: "dash" }}, hoverinfo: "skip", showlegend: false }});

                addBox(ex1[0], ex1[1], 2, 0.8, color_constitutive, "Upstream Exon 1");
                addBox(ex2[0], ex2[1], 2, 0.8, color_alt, "Alternative Skipped Exon 2", "<b>Alternative Skipped Exon</b><br>Cassette exon included in dominant isoform");
                addBox(ex3[0], ex3[1], 2, 0.8, color_constitutive, "Downstream Exon 3");

                const mid_1 = (ex1[1] + ex2[0]) / 2;
                traces.push({{
                    x: [ex1[1], mid_1, ex2[0]], y: [2.4, 3.2, 2.4],
                    mode: "lines+text", line: {{ color: color_inc_arc, width: 3 }},
                    text: ["", `Inclusion J1: ${{inc_counts}}`, ""],
                    textposition: "top center", name: "Inclusion Junction J1"
                }});
                const mid_2 = (ex2[1] + ex3[0]) / 2;
                traces.push({{
                    x: [ex2[1], mid_2, ex3[0]], y: [2.4, 3.2, 2.4],
                    mode: "lines", line: {{ color: color_inc_arc, width: 3 }},
                    name: "Inclusion Junction J2"
                }});

                addBox(ex1[0], ex1[1], -2, 0.8, color_constitutive, "Upstream Exon 1");
                addBox(ex3[0], ex3[1], -2, 0.8, color_constitutive, "Downstream Exon 3");

                const mid_exc = (ex1[1] + ex3[0]) / 2;
                traces.push({{
                    x: [ex1[1], mid_exc, ex3[0]], y: [-2.4, -3.5, -2.4],
                    mode: "lines+text", line: {{ color: color_exc_arc, width: 3 }},
                    text: ["", `Exclusion J3 (Skipping): ${{exc_counts}}`, ""],
                    textposition: "bottom center", name: "Exclusion Junction J3 (Skipping)"
                }});
            }}

            const layout = {{
                title: null,
                xaxis: {{
                    title: `<b>Genomic Coordinate (${{chrom}})</b>`,
                    range: [x_min, x_max], showgrid: true, gridcolor: "#F1F5F9",
                    tickformat: ",d"
                }},
                yaxis: {{
                    automargin: true,
                    showticklabels: true,
                    tickvals: [2, -2],
                    ticktext: [`<b>${{tick_inc}}</b>`, `<b>${{tick_exc}}</b>`],
                    range: [-4.8, 4.8], showgrid: false
                }},
                plot_bgcolor: "#FFFFFF",
                paper_bgcolor: "#FFFFFF",
                showlegend: true,
                legend: {{ orientation: "h", yanchor: "bottom", y: 1.02, xanchor: "right", x: 1, font: {{ size: 11 }} }},
                margin: {{ l: 180, r: 40, t: 40, b: 40 }},
                height: 380,
                shapes: shapes
            }};

            Plotly.newPlot(sashimiDiv, traces, layout, {{ responsive: true, displayModeBar: false }});
        }}

        function handleMasterGeneSelection(geneSymbol, eventIndex = 0) {{
            if (ncbiSelect && ncbiSelect.value !== geneSymbol) ncbiSelect.value = geneSymbol;
            if (sashimiGeneSelect && sashimiGeneSelect.value !== geneSymbol) sashimiGeneSelect.value = geneSymbol;
            renderNcbiGeneDetails(geneSymbol);
            filterIsoformTable(geneSymbol);
            renderSashimiPlot(geneSymbol, eventIndex);
            renderPrimerDetails(geneSymbol, eventIndex);
        }}

        function getActiveSplicingGeneSet() {{
            if (currentFilteredGenes && currentFilteredGenes.length > 0) {{
                return new Set(currentFilteredGenes
                    .filter(g => g.current_splicing_status === 'Inclusion' || g.current_splicing_status === 'Exclusion')
                    .map(g => (g.geneSymbol || '').toUpperCase())
                    .filter(Boolean));
            }}
            return new Set((rawGeneData || [])
                .filter(g => g.current_splicing_status === 'Inclusion' || g.current_splicing_status === 'Exclusion')
                .map(g => (g.geneSymbol || '').toUpperCase())
                .filter(Boolean));
        }}

        function updateIsoformDownloadButton() {{
            const dlBtn = document.getElementById('download-isoform-csv-btn');
            if (!dlBtn) return;
            const activeSet = getActiveSplicingGeneSet();
            const activeRecords = (isoformRecords || []).filter(r => activeSet.has((r.gene_symbol || '').toUpperCase()));
            dlBtn.textContent = `📥 Download Filtered Isoforms (${{activeRecords.length}}) (.csv)`;
            dlBtn.title = `Export ${{activeRecords.length}} isoform annotations for ${{activeSet.size}} active splicing target genes meeting current thresholds`;
        }}

        let isoformShowAll = false;
        function renderIsoformRow(r) {{
            const isNmd = (r.nmd_prediction || '').includes("NMD Sensitive");
            const nmdStyle = isNmd 
                ? "background:#FEF2F2; color:#EF4444; border:1px solid #EF4444;" 
                : "background:#EFF6FF; color:#2563EB; border:1px solid #3B82F6;";
            const tier = r.impairment_tier || '';
            const riskStyle = tier.includes("High") 
                ? "background:#FEF2F2; color:#DC2626; border:1px solid #EF4444; font-weight:800;" 
                : (tier.includes("Moderate") 
                    ? "background:#FFFBEB; color:#D97706; border:1px solid #F59E0B; font-weight:800;" 
                    : "background:#F0FDF4; color:#16A34A; border:1px solid #22C55E; font-weight:800;");
            const dpsi = (r.delta_psi !== undefined && r.delta_psi !== null) ? Number(r.delta_psi).toFixed(3) : "0.000";
            const fc = (r.log2FoldChange !== undefined && r.log2FoldChange !== null) ? Number(r.log2FoldChange).toFixed(3) : "0.000";
            const fcColor = (r.log2FoldChange || 0) > 0 ? "#DC2626" : "#2563EB";

            return `
                <tr class="isoform-data-row" data-gene="${{(r.gene_symbol || '').toUpperCase()}}" data-coords="${{r.coordinates || ''}}" data-event="${{r.event_type || ''}}" style="cursor: pointer;" title="Click to inspect Exon Structure & Primers">
                    <td><b>${{r.gene_symbol || ''}}</b></td>
                    <td><span class="badge" style="${{riskStyle}}">${{r.impairment_tier || 'N/A'}}</span></td>
                    <td><span style="font-size:12px; font-weight:600; color:#334155;">${{r.primary_dysfunction_cause || 'N/A'}}</span></td>
                    <td><span class="badge" style="background:#F1F5F9; color:#334155; font-weight:700;">${{r.event_type || 'N/A'}}</span></td>
                    <td><code>${{r.transcript_id || 'N/A'}}</code></td>
                    <td><span style="font-size:12px; color:#64748B;">${{r.coordinates || 'N/A'}}</span></td>
                    <td><b style="color:#2563EB;">${{dpsi}}</b></td>
                    <td><b style="color:${{fcColor}};">${{fc}}</b></td>
                    <td>${{r.cds_frame || 'N/A'}}</td>
                    <td><span class="badge" style="${{nmdStyle}}">${{r.nmd_prediction || 'N/A'}}</span><br><span style="font-size:11px; color:#64748B;">${{r.ptc_position || 'N/A'}}</span></td>
                    <td><b>${{r.protein_domain || 'N/A'}}</b></td>
                    <td><span style="font-size:12px; color:#475569;">${{r.localization_consequence || 'N/A'}}</span></td>
                </tr>
            `;
        }}

        function filterIsoformTable(geneSymbol) {{
            const upper = (geneSymbol || '').trim().toUpperCase();
            const filterLabel = document.getElementById('isoform-filter-gene');
            const countLabel = document.getElementById('isoform-filter-count');
            const toggleBtn = document.getElementById('isoform-toggle-view-btn');
            const tbody = document.getElementById('isoform-table-body');

            const activeSet = getActiveSplicingGeneSet();
            const activeGeneCount = activeSet.size;

            updateIsoformDownloadButton();

            if (filterLabel) {{
                if (isoformShowAll) {{
                    filterLabel.textContent = `All Active Genes (${{activeGeneCount}})`;
                }} else {{
                    filterLabel.textContent = upper || 'None Selected';
                }}
            }}

            let matching = [];
            if (isoformShowAll) {{
                matching = (isoformRecords || []).filter(r => activeSet.has((r.gene_symbol || '').toUpperCase()));
            }} else {{
                if (upper && activeSet.has(upper)) {{
                    matching = (isoformRecords || []).filter(r => (r.gene_symbol || '').toUpperCase() === upper);
                }}
            }}

            const visibleCount = matching.length;

            if (tbody) {{
                if (!matching.length) {{
                    const msg = (activeSet.size === 0)
                        ? 'No significant alternative splicing genes meet current threshold cutoffs.'
                        : ((!upper || !activeSet.has(upper))
                            ? `${{upper || 'Selected gene'}} is outside active splicing thresholds.`
                            : `No event annotations recorded for ${{upper}}.`);
                    tbody.innerHTML = `<tr><td colspan="12" style="text-align:center; color:#64748B; padding:20px; font-weight:600;">${{msg}}</td></tr>`;
                }} else {{
                    const maxDomRows = 500;
                    const rowsToRender = matching.slice(0, maxDomRows);
                    let html = rowsToRender.map(renderIsoformRow).join('');
                    if (matching.length > maxDomRows) {{
                        html += `<tr><td colspan="12" style="text-align:center; color:#64748B; padding:12px; font-weight:600; background:#F8FAFC;">(Showing first ${{maxDomRows}} of ${{matching.length}} isoform events. Use 'Download Filtered Isoforms (.csv)' to export all records.)</td></tr>`;
                    }}
                    tbody.innerHTML = html;
                }}
            }}

            if (countLabel) {{
                if (isoformShowAll) {{
                    if (activeGeneCount > 0) {{
                        countLabel.textContent = `(Showing all ${{visibleCount}} isoform event${{visibleCount > 1 ? 's' : ''}} across ${{activeGeneCount}} gene${{activeGeneCount > 1 ? 's' : ''}} meeting thresholds)`;
                    }} else {{
                        countLabel.textContent = `(0 active splicing targets meeting current thresholds)`;
                    }}
                }} else {{
                    if (upper && activeSet.has(upper)) {{
                        countLabel.textContent = visibleCount > 0 
                            ? `(${{visibleCount}} isoform event${{visibleCount > 1 ? 's' : ''}} for ${{upper}})` 
                            : `(No event annotations recorded for ${{upper}})`;
                    }} else if (upper) {{
                        countLabel.textContent = `(${{upper}} is outside active splicing thresholds)`;
                    }} else {{
                        countLabel.textContent = `(0 active targets selected)`;
                    }}
                }}
            }}

            if (toggleBtn) {{
                if (isoformShowAll) {{
                    toggleBtn.textContent = upper ? `🎯 Show Selected Gene (${{upper}}) Only` : `🎯 Show Selected Gene Only`;
                }} else {{
                    toggleBtn.textContent = `👁️ Show All Active Genes (${{activeGeneCount}})`;
                }}
            }}
        }}

        function toggleIsoformView() {{
            isoformShowAll = !isoformShowAll;
            const currentGene = ncbiSelect ? ncbiSelect.value : '';
            filterIsoformTable(currentGene);
        }}

        function downloadIsoformCSV() {{
            if (!isoformRecords || !isoformRecords.length) return;
            const activeSet = getActiveSplicingGeneSet();
            const activeRecords = isoformRecords.filter(r => activeSet.has((r.gene_symbol || '').toUpperCase()));

            if (!activeRecords.length) {{
                alert("No isoform records match the current threshold criteria. Please adjust thresholds to export isoforms.");
                return;
            }}

            const headers = ["Gene", "Impairment Tier", "Primary Dysfunction Cause", "Event Type", "Transcript ID", "Coordinates", "Delta PSI", "Log2FC", "CDS Frame", "PTC & NMD Status", "PTC Position", "Protein Domain Impact", "Localization Consequence"];
            let csvContent = headers.join(",") + "\\n";
            activeRecords.forEach(r => {{
                const row = [
                    `"${{r.gene_symbol || ''}}"`,
                    `"${{r.impairment_tier || ''}}"`,
                    `"${{r.primary_dysfunction_cause || ''}}"`,
                    `"${{r.event_type || ''}}"`,
                    `"${{r.transcript_id || ''}}"`,
                    `"${{r.coordinates || ''}}"`,
                    (r.delta_psi !== undefined ? r.delta_psi : 0).toFixed(4),
                    (r.log2FoldChange !== undefined ? r.log2FoldChange : 0).toFixed(4),
                    `"${{r.cds_frame || ''}}"`,
                    `"${{r.nmd_prediction || ''}}"`,
                    `"${{r.ptc_position || ''}}"`,
                    `"${{r.protein_domain || ''}}"`,
                    `"${{r.localization_consequence || ''}}"`
                ];
                csvContent += row.join(",") + "\\n";
            }});

            const psiCut = psiSlider ? parseFloat(psiSlider.value).toFixed(2) : '0.10';
            const fdrCut = fdrSlider ? parseFloat(fdrSlider.value).toFixed(3) : '0.050';
            const blob = new Blob([csvContent], {{ type: 'text/csv;charset=utf-8;' }});
            const link = document.createElement("a");
            link.href = URL.createObjectURL(blob);
            link.setAttribute("download", `gensplice_isoform_annotations_dPSI${{psiCut}}_FDR${{fdrCut}}.csv`);
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }}

        function renderPrimerDetails(geneSymbol, eventIndex = 0) {{
            const tbody = document.getElementById('primer-table-body');
            const bannerGene = document.getElementById('primer-banner-gene');
            const bannerCoords = document.getElementById('primer-banner-coords');
            const bannerEvent = document.getElementById('primer-banner-event');

            if (!geneSymbol || !primerRecords) {{
                if (tbody) tbody.innerHTML = `<tr><td colspan="7" style="text-align:center; color:#64748B; padding:16px;">No active splicing target selected under current thresholds.</td></tr>`;
                if (bannerGene) bannerGene.textContent = "-";
                if (bannerCoords) bannerCoords.textContent = "-";
                if (bannerEvent) bannerEvent.textContent = "-";
                return;
            }}

            const upper = (geneSymbol || '').trim().toUpperCase();
            const events = (sashimiEventsMap && sashimiEventsMap[upper]) || [];
            const geneRows = primerRecords.filter(r => (r.gene_symbol || '').toUpperCase() === upper);

            const evIdx = parseInt(eventIndex !== undefined && eventIndex !== null ? eventIndex : 0, 10);
            let targetEvent = (evIdx >= 0 && events.length > evIdx) ? events[evIdx] : (events.length > 0 ? events[0] : null);

            let rows = [];
            if (evIdx === -1) {{
                rows = geneRows;
                if (bannerGene) bannerGene.textContent = geneSymbol;
                if (bannerCoords) bannerCoords.textContent = `All Coordinates (${{events.length}} events)`;
                if (bannerEvent) bannerEvent.textContent = `All Events (${{events.length}} total)`;
            }} else {{
                if (targetEvent) {{
                    if (bannerGene) bannerGene.textContent = geneSymbol;
                    if (bannerCoords) bannerCoords.textContent = targetEvent.coordinates || 'N/A';
                    const evNum = events.length > 1 ? `Event ${{evIdx + 1}} of ${{events.length}}` : `Event 1`;
                    if (bannerEvent) bannerEvent.textContent = `${{targetEvent.event_type}} (${{evNum}})`;

                    rows = geneRows.filter(r => r.event_index === evIdx);
                    if (!rows.length && targetEvent.coordinates) {{
                        rows = geneRows.filter(r => r.coordinates === targetEvent.coordinates);
                    }}
                    if (!rows.length) {{
                        rows = geneRows.slice(0, 2);
                    }}
                }} else {{
                    const ginfo = (sashimiGeneData && (sashimiGeneData[upper] || sashimiGeneData[geneSymbol])) || {{}};
                    if (bannerGene) bannerGene.textContent = geneSymbol;
                    if (bannerCoords) bannerCoords.textContent = (geneRows.length && geneRows[0].coordinates) ? geneRows[0].coordinates : (ginfo.coordinates || 'N/A');
                    if (bannerEvent) bannerEvent.textContent = (geneRows.length && geneRows[0].event_type) ? geneRows[0].event_type : (ginfo.event_type || 'SE');
                    rows = geneRows.slice(0, 2);
                }}
            }}

            if (!tbody) return;

            if (!rows.length) {{
                const targetCoordStr = targetEvent ? targetEvent.coordinates : 'N/A';
                tbody.innerHTML = `<tr><td colspan="7" style="text-align:center; color:#64748B; padding:16px;">Specific primer pairs not pre-calculated for <b>${{geneSymbol}}</b> (Target coordinates: <code>${{targetCoordStr}}</code>).</td></tr>`;
                return;
            }}

            let html = '';
            rows.forEach(r => {{
                const isInc = (r.target_isoform || '').toLowerCase() === 'inclusion';
                const badge = isInc ?
                    '<span class="badge" style="background:#EFF6FF; color:#2563EB; border:1px solid #3B82F6;">Inclusion Isoform</span>' :
                    '<span class="badge" style="background:#FDF2F8; color:#DB2777; border:1px solid #EC4899;">Exclusion Isoform</span>';
                
                const fwdTm = (r.fwd_tm_celsius !== null && r.fwd_tm_celsius !== undefined) ? Number(r.fwd_tm_celsius).toFixed(1) + '°C' : 'N/A';
                const fwdGc = (r.fwd_gc_pct !== null && r.fwd_gc_pct !== undefined) ? Number(r.fwd_gc_pct).toFixed(1) + '%' : 'N/A';
                const revTm = (r.rev_tm_celsius !== null && r.rev_tm_celsius !== undefined) ? Number(r.rev_tm_celsius).toFixed(1) + '°C' : 'N/A';
                const revGc = (r.rev_gc_pct !== null && r.rev_gc_pct !== undefined) ? Number(r.rev_gc_pct).toFixed(1) + '%' : 'N/A';
                const ampSize = (r.amplicon_size_bp !== null && r.amplicon_size_bp !== undefined) ? r.amplicon_size_bp + ' bp' : 'N/A';
                const fwdSeqHtml = (r.fwd_sequence && r.fwd_sequence !== 'N/A') ? `<code style="font-weight:700; color:#1E40AF;">${{r.fwd_sequence}}</code>` : `<span style="color:#64748B;">N/A</span>`;
                const revSeqHtml = (r.rev_sequence && r.rev_sequence !== 'N/A') ? `<code style="font-weight:700; color:#1E40AF;">${{r.rev_sequence}}</code>` : `<span style="color:#64748B;">N/A</span>`;

                html += `
                <tr>
                    <td><b>${{r.gene_symbol}}</b></td>
                    <td>${{badge}}</td>
                    <td><b>${{r.target_region || (r.target_isoform + ' Junction')}}</b><br><span style="font-size:11px; color:#64748B;">Coords: <code>${{r.coordinates || 'N/A'}}</code></span></td>
                    <td>${{fwdSeqHtml}}<br><span style="font-size:11px; color:#64748B;">Tm: ${{fwdTm}} | GC: ${{fwdGc}}</span></td>
                    <td>${{revSeqHtml}}<br><span style="font-size:11px; color:#64748B;">Tm: ${{revTm}} | GC: ${{revGc}}</span></td>
                    <td><b>${{ampSize}}</b></td>
                    <td><span class="badge" style="background:#F0FDF4; color:#16A34A; border:1px solid #22C55E;">${{r.primer_quality}}</span></td>
                </tr>
                `;
            }});
            tbody.innerHTML = html;
        }}

        function updateThresholds() {{
            const psiCut = parseFloat(psiSlider.value);
            const fdrCut = parseFloat(fdrSlider.value);

            psiValLabel.textContent = psiCut.toFixed(2);
            fdrValLabel.textContent = fdrCut.toFixed(3);

            let incCount = 0, excCount = 0, nonSigCount = 0;
            currentFilteredGenes = [];

            const volcanoPoints = {{
                Inclusion: {{ x: [], y: [], text: [] }},
                Exclusion: {{ x: [], y: [], text: [] }},
                NonSig: {{ x: [], y: [], text: [] }}
            }};

            rawGeneData.forEach(g => {{
                const dpsi = g.delta_psi !== undefined && g.delta_psi !== null ? Number(g.delta_psi) : 0.0;
                const fdr = g.as_fdr !== undefined && g.as_fdr !== null ? Number(g.as_fdr) : 1.0;
                const logFdr = (fdr <= 0) ? 20.0 : -Math.log10(Math.max(fdr, 1e-20));
                const isSig = (Math.abs(dpsi) >= psiCut) && (fdr <= fdrCut);

                let status = "Non-Significant";
                let statusLabel = "Non-Significant Background";

                if (isSig && dpsi > 0) {{
                    status = "Inclusion";
                    statusLabel = `Significant Inclusion Favored (ΔPSI ≥ +${{psiCut.toFixed(2)}})`;
                    incCount++;
                }} else if (isSig && dpsi < 0) {{
                    status = "Exclusion";
                    statusLabel = `Significant Exclusion Favored (ΔPSI ≤ -${{psiCut.toFixed(2)}})`;
                    excCount++;
                }} else {{
                    nonSigCount++;
                }}

                const sign = dpsi >= 0 ? '+' : '';
                const hoverText = `<b>${{g.geneSymbol || 'Unknown'}}</b> (${{g.gene_id || 'N/A'}})<br>` +
                    `Status: <b>${{statusLabel}}</b><br>` +
                    `Event: ${{g.event_type || 'N/A'}} | Coords: ${{g.coordinates || 'N/A'}}<br>` +
                    `ΔPSI: ${{sign}}${{dpsi.toFixed(3)}}<br>` +
                    `rMATS FDR: ${{fdr.toExponential(2)}} (-log₁₀ FDR: ${{logFdr.toFixed(2)}})<br>` +
                    `rMATS p-value: ${{(g.as_pvalue || fdr).toExponential(2)}}<br>` +
                    `Read Counts: Inc=${{g.inc_counts || 0}} | Exc=${{g.exc_counts || 0}}`;

                const ptKey = (status === "Inclusion") ? "Inclusion" : ((status === "Exclusion") ? "Exclusion" : "NonSig");
                volcanoPoints[ptKey].x.push(dpsi);
                volcanoPoints[ptKey].y.push(logFdr);
                volcanoPoints[ptKey].text.push(hoverText);

                const updatedGene = {{ ...g, current_splicing_status: status }};
                currentFilteredGenes.push(updatedGene);
            }});

            document.getElementById('kpi-total').textContent = rawGeneData.length;
            document.getElementById('kpi-sig').textContent = incCount + excCount;
            document.getElementById('kpi-inc').textContent = incCount;
            document.getElementById('kpi-exc').textContent = excCount;
            document.getElementById('kpi-nonsig').textContent = nonSigCount;

            const activeSplicingGenes = currentFilteredGenes
                .filter(g => g.current_splicing_status === "Inclusion" || g.current_splicing_status === "Exclusion")
                .map(g => (g.geneSymbol || '').toUpperCase());

            updateNcbiDropdown(activeSplicingGenes);

            const currentGene = ncbiSelect ? ncbiSelect.value : '';
            filterIsoformTable(currentGene);
            updateIsoformDownloadButton();

            const volcanoDiv = document.getElementById('plotly-as-volcano-div');
            if (volcanoDiv && window.Plotly) {{
                let maxDataX = 0.8;
                let maxDataY = 3.5;
                rawGeneData.forEach(g => {{
                    const dpsi = Math.abs(g.delta_psi || 0);
                    const fdr = g.as_fdr !== undefined && g.as_fdr !== null ? Number(g.as_fdr) : 1.0;
                    const logFdr = (fdr <= 0) ? 20.0 : -Math.log10(Math.max(fdr, 1e-20));
                    if (isFinite(dpsi) && dpsi > maxDataX) maxDataX = dpsi;
                    if (isFinite(logFdr) && logFdr > maxDataY) maxDataY = logFdr;
                }});

                const fdr_threshold_y = -Math.log10(Math.max(fdrCut, 1e-20));
                const max_x = Math.max(Math.ceil((maxDataX + 0.1) * 100) / 100, psiCut + 0.1, 1.05);
                const max_y = Math.max(Math.ceil((maxDataY + 1.0) * 10) / 10, Math.ceil((fdr_threshold_y + 1.5) * 10) / 10, 5.0);

                const newShapes = [
                    {{ type: 'rect', x0: -max_x, x1: max_x, y0: 0, y1: fdr_threshold_y, fillcolor: 'rgba(241, 245, 249, 0.50)', line: {{ width: 0 }}, layer: 'below' }},
                    {{ type: 'rect', x0: -psiCut, x1: psiCut, y0: fdr_threshold_y, y1: max_y, fillcolor: 'rgba(241, 245, 249, 0.50)', line: {{ width: 0 }}, layer: 'below' }},
                    {{ type: 'rect', x0: psiCut, x1: max_x, y0: fdr_threshold_y, y1: max_y, fillcolor: 'rgba(37, 99, 235, 0.08)', line: {{ width: 0 }}, layer: 'below' }},
                    {{ type: 'rect', x0: -max_x, x1: -psiCut, y0: fdr_threshold_y, y1: max_y, fillcolor: 'rgba(225, 29, 72, 0.08)', line: {{ width: 0 }}, layer: 'below' }},

                    {{ type: 'line', x0: -max_x, x1: max_x, y0: fdr_threshold_y, y1: fdr_threshold_y, line: {{ color: '#6366F1', width: 2, dash: 'dash' }} }},
                    {{ type: 'line', x0: psiCut, x1: psiCut, y0: 0, y1: max_y, line: {{ color: '#2563EB', width: 2, dash: 'dash' }} }},
                    {{ type: 'line', x0: -psiCut, x1: -psiCut, y0: 0, y1: max_y, line: {{ color: '#E11D48', width: 2, dash: 'dash' }} }}
                ];

                const newAnnotations = [
                    {{
                        x: -max_x * 0.65, y: max_y * 0.95,
                        xref: 'x', yref: 'y',
                        text: '<b>Exclusion Favored</b> (ΔPSI < 0)',
                        showarrow: false,
                        font: {{ size: 12, color: '#BE123C', family: 'sans-serif' }},
                        bgcolor: 'rgba(255, 255, 255, 0.8)',
                        bordercolor: '#FECDD3',
                        borderwidth: 1,
                        borderpad: 4
                    }},
                    {{
                        x: max_x * 0.65, y: max_y * 0.95,
                        xref: 'x', yref: 'y',
                        text: '<b>Inclusion Favored</b> (ΔPSI > 0)',
                        showarrow: false,
                        font: {{ size: 12, color: '#1D4ED8', family: 'sans-serif' }},
                        bgcolor: 'rgba(255, 255, 255, 0.8)',
                        bordercolor: '#BFDBFE',
                        borderwidth: 1,
                        borderpad: 4
                    }},
                    {{
                        x: max_x * 0.96, y: fdr_threshold_y + 0.15,
                        xref: 'x', yref: 'y',
                        text: `FDR = ${{fdrCut.toFixed(3)}}`,
                        showarrow: false,
                        font: {{ size: 11, color: '#4F46E5', family: 'sans-serif' }},
                        xanchor: 'right'
                    }}
                ];

                const updatedTraces = [
                    {{
                        x: volcanoPoints.Inclusion.x, y: volcanoPoints.Inclusion.y, text: volcanoPoints.Inclusion.text,
                        mode: 'markers', name: 'Inclusion Favored',
                        marker: {{ color: '#2563EB', size: 10, opacity: 0.9, line: {{ width: 0.8, color: '#FFFFFF' }} }},
                        hoverinfo: 'text', visible: true
                    }},
                    {{
                        x: volcanoPoints.Exclusion.x, y: volcanoPoints.Exclusion.y, text: volcanoPoints.Exclusion.text,
                        mode: 'markers', name: 'Exclusion Favored',
                        marker: {{ color: '#E11D48', size: 10, opacity: 0.9, line: {{ width: 0.8, color: '#FFFFFF' }} }},
                        hoverinfo: 'text', visible: true
                    }},
                    {{
                        x: volcanoPoints.NonSig.x, y: volcanoPoints.NonSig.y, text: volcanoPoints.NonSig.text,
                        mode: 'markers', name: 'Non-Significant Background',
                        marker: {{ color: '#94A3B8', size: 6, opacity: 0.5, line: {{ width: 0.5, color: '#FFFFFF' }} }},
                        hoverinfo: 'text', visible: true
                    }}
                ];

                Plotly.react(volcanoDiv, updatedTraces, volcanoDiv.layout);
                Plotly.relayout(volcanoDiv, {{
                    shapes: newShapes,
                    annotations: newAnnotations,
                    'xaxis.range': [-max_x, max_x],
                    'yaxis.range': [0, max_y]
                }});
            }}

            if (isEnrichmentGenerated) {{
                generateEnrichment();
            }}
        }}

        function generateEnrichment() {{
            const psiCut = parseFloat(psiSlider.value).toFixed(2);
            const fdrCut = parseFloat(fdrSlider.value).toFixed(3);
            
            const activeGenes = currentFilteredGenes
                .filter(g => g.current_splicing_status === "Inclusion" || g.current_splicing_status === "Exclusion")
                .map(g => (g.geneSymbol || '').toUpperCase());
            const activeSet = new Set(activeGenes);

            const goPlaceholder = document.getElementById('go-placeholder');
            const keggPlaceholder = document.getElementById('kegg-placeholder');
            if (goPlaceholder) goPlaceholder.style.display = 'none';
            if (keggPlaceholder) keggPlaceholder.style.display = 'none';

            const goDiv = document.getElementById('plotly-go-div');
            const keggDiv = document.getElementById('plotly-kegg-div');
            if (goDiv) goDiv.style.display = 'block';
            if (keggDiv) keggDiv.style.display = 'block';

            if (!activeGenes.length) {{
                currentGoData = [];
                currentKeggData = [];
                if (goDiv && window.Plotly) {{
                    Plotly.react(goDiv, [], {{
                        title: "<b>No significant splicing target genes found under current thresholds (|ΔPSI| ≥ " + psiCut + ", FDR ≤ " + fdrCut + ")</b>",
                        paper_bgcolor: "#FFFFFF", plot_bgcolor: "#F8FAFC", height: 350
                    }});
                }}
                if (keggDiv && window.Plotly) {{
                    Plotly.react(keggDiv, [], {{
                        title: "<b>No significant splicing target genes found under current thresholds (|ΔPSI| ≥ " + psiCut + ", FDR ≤ " + fdrCut + ")</b>",
                        paper_bgcolor: "#FFFFFF", plot_bgcolor: "#F8FAFC", height: 350
                    }});
                }}
                return;
            }}

            function buildProcessedEnrichment(records) {{
                if (!records || !records.length) return [];
                const matched = [];
                records.forEach(r => {{
                    const rawGenes = (r.genes || '').split(';').map(g => g.trim().toUpperCase()).filter(Boolean);
                    const matching = rawGenes.filter(g => activeSet.has(g));
                    if (matching.length > 0) {{
                        matched.push({{
                            term: r.term,
                            count: matching.length,
                            overlap: `${{matching.length}}/${{activeGenes.length}}`,
                            pvalue: r.pvalue,
                            adjPvalue: r.adjPvalue,
                            logP: parseFloat(r.logP || 0),
                            genes: matching.join('; ')
                        }});
                    }}
                }});

                if (matched.length > 0) {{
                    return matched.sort((a,b) => b.logP - a.logP).slice(0, 10);
                }}

                return records.map(r => ({{
                    term: r.term,
                    count: (r.genes || '').split(';').length,
                    overlap: r.overlap,
                    pvalue: r.pvalue,
                    adjPvalue: r.adjPvalue,
                    logP: parseFloat(r.logP || 0),
                    genes: r.genes
                }})).slice(0, 10);
            }}

            const sortedGo = buildProcessedEnrichment(baseGoRecords);
            const sortedKegg = buildProcessedEnrichment(baseKeggRecords);

            currentGoData = sortedGo;
            currentKeggData = sortedKegg;

            if (goDiv && window.Plotly) {{
                if (sortedGo.length) {{
                    const yVal = sortedGo.map(item => item.term.length > 55 ? item.term.slice(0,55) + '...' : item.term).reverse();
                    const xVal = sortedGo.map(item => parseFloat(item.logP)).reverse();
                    const hoverText = sortedGo.map(item => `<b>${{item.term}}</b><br>-log₁₀(p): ${{item.logP}}<br>p-value: ${{item.pvalue}} | FDR: ${{item.adjPvalue}}<br>Overlap: ${{item.overlap}}<br>Associated Genes: ${{item.genes}}`).reverse();
                    Plotly.react(goDiv, [{{
                        x: xVal,
                        y: yVal,
                        type: 'bar',
                        orientation: 'h',
                        text: xVal.map(v => v.toFixed(1)),
                        textposition: 'auto',
                        marker: {{
                            color: '#3B82F6',
                            line: {{ color: '#1E40AF', width: 1 }}
                        }},
                        hoverinfo: 'text',
                        hovertext: hoverText
                    }}], {{
                        margin: {{ l: 280, r: 30, t: 30, b: 50 }},
                        xaxis: {{ title: '<b>-log₁₀(p-value)</b>' }},
                        yaxis: {{ automargin: true }},
                        template: 'plotly_white',
                        paper_bgcolor: '#FFFFFF',
                        plot_bgcolor: '#FFFFFF',
                        height: 380
                    }});
                    Plotly.Plots.resize(goDiv);
                }} else {{
                    Plotly.react(goDiv, [], {{
                        title: "<b>Enrichr API query unavailable: network offline or server timeout</b>",
                        paper_bgcolor: "#FFFFFF", plot_bgcolor: "#F8FAFC", height: 260
                    }});
                }}
            }}

            if (keggDiv && window.Plotly) {{
                if (sortedKegg.length) {{
                    const yVal = sortedKegg.map(item => item.term.length > 55 ? item.term.slice(0,55) + '...' : item.term).reverse();
                    const xVal = sortedKegg.map(item => parseFloat(item.logP)).reverse();
                    const hoverText = sortedKegg.map(item => `<b>${{item.term}}</b><br>-log₁₀(p): ${{item.logP}}<br>p-value: ${{item.pvalue}} | FDR: ${{item.adjPvalue}}<br>Overlap: ${{item.overlap}}<br>Associated Genes: ${{item.genes}}`).reverse();
                    Plotly.react(keggDiv, [{{
                        x: xVal,
                        y: yVal,
                        type: 'bar',
                        orientation: 'h',
                        text: xVal.map(v => v.toFixed(1)),
                        textposition: 'auto',
                        marker: {{
                            color: '#8B5CF6',
                            line: {{ color: '#6D28D9', width: 1 }}
                        }},
                        hoverinfo: 'text',
                        hovertext: hoverText
                    }}], {{
                        margin: {{ l: 280, r: 30, t: 30, b: 50 }},
                        xaxis: {{ title: '<b>-log₁₀(p-value)</b>' }},
                        yaxis: {{ automargin: true }},
                        template: 'plotly_white',
                        paper_bgcolor: '#FFFFFF',
                        plot_bgcolor: '#FFFFFF',
                        height: 380
                    }});
                    Plotly.Plots.resize(keggDiv);
                }} else {{
                    Plotly.react(keggDiv, [], {{
                        title: "<b>Enrichr API query unavailable: network offline or server timeout</b>",
                        paper_bgcolor: "#FFFFFF", plot_bgcolor: "#F8FAFC", height: 260
                    }});
                }}
            }}

            let noticeText = `✔ GO & KEGG Enrichment Generated for |ΔPSI| ≥ ${{psiCut}} & FDR ≤ ${{fdrCut}} • Active Splicing Targets: ${{activeGenes.length}}`;
            if (!sortedGo.length && !sortedKegg.length) {{
                noticeText = "⚠️ Enrichr API query unavailable: network offline or server timeout (no synthetic mock pathways displayed)";
            }}
            const goNot = document.getElementById('go-notice');
            const keggNot = document.getElementById('kegg-notice');
            if (goNot) goNot.textContent = noticeText;
            if (keggNot) keggNot.textContent = noticeText;
        }}

        function triggerEnrichmentWithAction() {{
            if (!generateBtn) return;
            generateBtn.innerHTML = '⚡ Generating GO & KEGG...';
            isEnrichmentGenerated = true;
            
            const goCard = document.getElementById('go-card-section');
            const keggCard = document.getElementById('kegg-card-section');
            if (goCard) {{
                goCard.classList.remove('card-updating');
                void goCard.offsetWidth;
                goCard.classList.add('card-updating');
            }}
            if (keggCard) {{
                keggCard.classList.remove('card-updating');
                void keggCard.offsetWidth;
                keggCard.classList.add('card-updating');
            }}

            setTimeout(() => {{
                generateEnrichment();
                generateBtn.innerHTML = '🚀 Generate GO / KEGG';
                if (goCard) {{
                    goCard.scrollIntoView({{ behavior: 'smooth', block: 'start' }});
                }}
            }}, 200);
        }}

        function downloadFilteredCSV() {{
            if (!currentFilteredGenes.length) return;

            const headers = ["geneSymbol", "gene_id", "splicing_status", "delta_psi", "as_fdr", "as_pvalue", "event_type", "coordinates", "inc_counts", "exc_counts", "log2FoldChange", "deg_fdr"];
            let csvContent = headers.join(",") + "\\n";

            currentFilteredGenes.forEach(g => {{
                const row = [
                    `"${{g.geneSymbol || ''}}"`,
                    `"${{g.gene_id || ''}}"`,
                    `"${{g.current_splicing_status || ''}}"`,
                    (g.delta_psi !== undefined ? Number(g.delta_psi) : 0).toFixed(4),
                    (g.as_fdr !== undefined ? Number(g.as_fdr) : 1.0).toExponential(3),
                    (g.as_pvalue !== undefined ? Number(g.as_pvalue) : 1.0).toExponential(3),
                    `"${{g.event_type || 'None'}}"`,
                    `"${{g.coordinates || 'N/A'}}"` ,
                    g.inc_counts || 0,
                    g.exc_counts || 0,
                    (g.log2FoldChange !== undefined ? Number(g.log2FoldChange) : 0).toFixed(4),
                    (g.deg_fdr !== undefined ? Number(g.deg_fdr) : 1.0).toExponential(3)
                ];
                csvContent += row.join(",") + "\\n";
            }});

            const psiCut = psiSlider.value;
            const fdrCut = fdrSlider.value;
            const blob = new Blob([csvContent], {{ type: 'text/csv;charset=utf-8;' }});
            const link = document.createElement("a");
            const url = URL.createObjectURL(blob);
            
            link.setAttribute("href", url);
            link.setAttribute("download", `gensplice_splicing_targets_dPSI${{psiCut}}_FDR${{fdrCut}}.csv`);
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }}

        function downloadGoCSV() {{
            if (!isEnrichmentGenerated || !currentGoData || !currentGoData.length) {{
                alert("Please click '🚀 Generate GO / KEGG' first to compute enrichment data before downloading.");
                return;
            }}
            let csvContent = "GO Term,Overlap,P-value,Adjusted P-value,Associated Genes\\n";
            currentGoData.forEach(r => {{
                csvContent += `"${{r.term}}","${{r.overlap}}","${{r.pvalue}}","${{r.adjPvalue}}","${{r.genes}}"\\n`;
            }});
            const blob = new Blob([csvContent], {{ type: 'text/csv;charset=utf-8;' }});
            const link = document.createElement("a");
            link.href = URL.createObjectURL(blob);
            link.setAttribute("download", `go_enrichment_splicing_targets_dPSI${{psiSlider.value}}_FDR${{fdrSlider.value}}.csv`);
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }}

        function downloadKeggCSV() {{
            if (!isEnrichmentGenerated || !currentKeggData || !currentKeggData.length) {{
                alert("Please click '🚀 Generate GO / KEGG' first to compute enrichment data before downloading.");
                return;
            }}
            let csvContent = "KEGG Pathway,Overlap,P-value,Adjusted P-value,Associated Genes\\n";
            currentKeggData.forEach(r => {{
                csvContent += `"${{r.term}}","${{r.overlap}}","${{r.pvalue}}","${{r.adjPvalue}}","${{r.genes}}"\\n`;
            }});
            const blob = new Blob([csvContent], {{ type: 'text/csv;charset=utf-8;' }});
            const link = document.createElement("a");
            link.href = URL.createObjectURL(blob);
            link.setAttribute("download", `kegg_enrichment_splicing_targets_dPSI${{psiSlider.value}}_FDR${{fdrSlider.value}}.csv`);
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }}

        psiSlider.addEventListener('input', updateThresholds);
        fdrSlider.addEventListener('input', updateThresholds);
        downloadBtn.addEventListener('click', downloadFilteredCSV);
        if (downloadGoBtn) downloadGoBtn.addEventListener('click', downloadGoCSV);
        if (downloadKeggBtn) downloadKeggBtn.addEventListener('click', downloadKeggCSV);
        if (ncbiSelect) {{
            ncbiSelect.addEventListener('change', (e) => handleMasterGeneSelection(e.target.value));
        }}
        if (sashimiGeneSelect) {{
            sashimiGeneSelect.addEventListener('change', (e) => handleMasterGeneSelection(e.target.value));
        }}
        if (sashimiEventSelect) {{
            sashimiEventSelect.addEventListener('change', (e) => {{
                const curGene = sashimiGeneSelect ? sashimiGeneSelect.value : (ncbiSelect ? ncbiSelect.value : '');
                const evIdx = parseInt(e.target.value, 10);
                renderSashimiPlot(curGene, evIdx >= 0 ? evIdx : 0);
                renderPrimerDetails(curGene, evIdx);
            }});
        }}
        const isoformTbody = document.getElementById('isoform-table-body');
        if (isoformTbody) {{
            isoformTbody.addEventListener('click', (e) => {{
                const tr = e.target.closest('tr.isoform-data-row');
                if (tr && tr.dataset.gene) {{
                    const g = tr.dataset.gene;
                    const coords = tr.dataset.coords;
                    let targetIdx = 0;
                    const events = (sashimiEventsMap && sashimiEventsMap[g]) || [];
                    if (coords && events.length) {{
                        const found = events.findIndex(ev => ev.coordinates === coords);
                        if (found >= 0) targetIdx = found;
                    }}
                    handleMasterGeneSelection(g, targetIdx);
                }}
            }});
        }}
        if (downloadIsoformBtn) {{
            downloadIsoformBtn.addEventListener('click', downloadIsoformCSV);
        }}
        if (toggleIsoformBtn) {{
            toggleIsoformBtn.addEventListener('click', toggleIsoformView);
        }}
        if (generateBtn) {{
            generateBtn.addEventListener('click', triggerEnrichmentWithAction);
        }}

        updateThresholds();
    </script>
</body>
</html>
"""

    with open(output_html_path, "w", encoding="utf-8") as f:
        f.write(full_html)

    try:
        os.chmod(output_html_path, 0o644)
    except Exception:
        pass

    print(f"  ✔ Standalone Light Mode HTML Report updated: {output_html_path}")
    return output_html_path
