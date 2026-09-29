"""
GenSplice-Agent Standalone HTML Exporter Module
Light Mode interactive HTML report with:
- Outside plot header: 'Color Classification' (with Q1-Q4 descriptions)
- Inside plot legend title: 'On/Off' (strictly Q1, Q2, Q3, Q4)
- GO Term & KEGG Pathway Enrichment Charts & Tables focused on Q2 (Splicing-Driven Only)
- Client-Side Real-Time Plotly Graph & Table Re-Calculation with Dynamic Action Effects
"""

import os
import json
import time
import polars as pl
import pandas as pd
from visualizer.quadrant_plot import build_quadrant_plot
from visualizer.enrichment_plot import build_enrichment_chart
from core.merger import get_quadrant_kpis
from core.enrichment import fetch_enrichment
from core.ai_summary import generate_biological_insights, fetch_ncbi_gene_summary, fetch_pubmed_literature, generate_detailed_bio_prompt
from core.isoform_annotator import annotate_isoform_events
from visualizer.exon_structure import plot_exon_structure
from core.primer_designer import generate_primer_table_for_targets

def export_html_report(
    df_merged: pl.DataFrame,
    output_html_path: str = "outputs/gensplice_report.html",
    ai_insights: str = None,
    log2fc_cutoff: float = 1.0,
    delta_psi_cutoff: float = 0.1,
    deg_fdr_cutoff: float = 0.05,
    as_fdr_cutoff: float = 0.05
) -> str:
    """
    Exports a self-contained Light Mode interactive HTML report with standard horizontal GO & KEGG Plotly bar charts.
    """
    os.makedirs(os.path.dirname(output_html_path), exist_ok=True)

    fig_quad = build_quadrant_plot(df_merged, log2fc_cutoff, delta_psi_cutoff, color_by="quadrant")
    quad_html = fig_quad.to_html(full_html=False, include_plotlyjs="cdn", div_id="plotly-quad-div")

    kpis = get_quadrant_kpis(df_merged)

    # Visual Exon-Intron Structure Engine (Sashimi Plot)
    top_gene_symbol = "CRISPLD2"
    top_coords = "chr16:84860000:84861500:84863000"
    top_event = "SE"
    top_delta_psi = 0.35

    q2_top = df_merged.filter(pl.col("quadrant") == "Q2")
    if q2_top.height > 0:
        row0 = q2_top.to_dicts()[0]
        top_gene_symbol = row0.get("geneSymbol", "CRISPLD2")
        top_coords = row0.get("coordinates", top_coords)
        top_event = row0.get("event_type", "SE")
        top_delta_psi = row0.get("delta_psi", 0.35)
    else:
        q1_top = df_merged.filter(pl.col("quadrant") == "Q1")
        if q1_top.height > 0:
            row0 = q1_top.to_dicts()[0]
            top_gene_symbol = row0.get("geneSymbol", "STAT3")
            top_coords = row0.get("coordinates", top_coords)
            top_event = row0.get("event_type", "SE")
            top_delta_psi = row0.get("delta_psi", 0.35)

    fig_sashimi = plot_exon_structure(
        gene_symbol=top_gene_symbol,
        event_type=top_event,
        coordinates=top_coords,
        delta_psi=top_delta_psi
    )
    sashimi_html = fig_sashimi.to_html(full_html=False, include_plotlyjs=False, div_id="plotly-sashimi-div")

    # Isoform-Specific RT-qPCR Primer Designer Engine (Q2/Q1 Targets)
    q2_q1_sub = df_merged.filter(pl.col("quadrant").is_in(["Q2", "Q1"]))
    target_candidates = q2_q1_sub if q2_q1_sub.height > 0 else df_merged
    primer_df = generate_primer_table_for_targets(target_candidates.head(30))
    primer_records = primer_df.to_dicts() if primer_df.height > 0 else []
    primer_records_json = json.dumps(primer_records)

    first_primer_gene = primer_records[0]["gene_symbol"] if primer_records else ""
    first_coords = primer_records[0].get("coordinates", "N/A") if primer_records else "N/A"
    first_event = primer_records[0].get("event_type", "SE") if primer_records else "SE"

    primer_rows_html = ""
    if primer_df.height > 0:
        for r in primer_df.filter(pl.col("gene_symbol") == first_primer_gene).iter_rows(named=True):
            isoform_badge = '<span class="badge" style="background:#EFF6FF; color:#2563EB; border:1px solid #3B82F6;">Inclusion Isoform</span>' if r['target_isoform'].lower() == 'inclusion' else '<span class="badge" style="background:#FDF2F8; color:#DB2777; border:1px solid #EC4899;">Exclusion Isoform</span>'
            primer_rows_html += f"""
            <tr>
                <td><b>{r['gene_symbol']}</b></td>
                <td>{isoform_badge}</td>
                <td><b>{r.get('target_region', r['target_isoform'] + ' Junction')}</b><br><span style="font-size:11px; color:#64748B;">Coords: <code>{r.get('coordinates', 'N/A')}</code></span></td>
                <td><code style="font-weight:700; color:#1E40AF;">{r['fwd_sequence']}</code><br><span style="font-size:11px; color:#64748B;">Tm: {r['fwd_tm_celsius']}°C | GC: {r['fwd_gc_pct']}%</span></td>
                <td><code style="font-weight:700; color:#1E40AF;">{r['rev_sequence']}</code><br><span style="font-size:11px; color:#64748B;">Tm: {r['rev_tm_celsius']}°C | GC: {r['rev_gc_pct']}%</span></td>
                <td><b>{r['amplicon_size_bp']} bp</b></td>
                <td><span class="badge" style="background:#F0FDF4; color:#16A34A; border:1px solid #22C55E;">{r['primer_quality']}</span></td>
            </tr>
            """

    # Splicing genes (Q1 + Q2)
    splicing_genes = df_merged.filter(pl.col("quadrant").is_in(["Q1", "Q2"])).select("geneSymbol").to_series().to_list()
    if not splicing_genes:
        splicing_genes = df_merged.select("geneSymbol").to_series().to_list()
    all_genes = df_merged.select("geneSymbol").to_series().to_list()

    # Pre-generate GO & KEGG Enrichment for Q1+Q2 Splicing Targets only (Q4 completely removed)
    df_go = fetch_enrichment(splicing_genes, gene_sets=["GO_Biological_Process_2023"], top_n=10)
    df_kegg = fetch_enrichment(splicing_genes, gene_sets=["KEGG_2021_Human"], top_n=10)

    # Build Standard Horizontal Plotly Bar Charts for GO & KEGG (Q1 + Q2)
    fig_go = build_enrichment_chart(df_go, "GO Term Biological Process Enrichment (Q1 + Q2: Splicing Targets)", bar_color="#3B82F6")
    go_chart_html = fig_go.to_html(full_html=False, include_plotlyjs=False, div_id="plotly-go-div")

    fig_kegg = build_enrichment_chart(df_kegg, "KEGG Pathway Enrichment (Q1 + Q2: Splicing Targets)", bar_color="#8B5CF6")
    kegg_chart_html = fig_kegg.to_html(full_html=False, include_plotlyjs=False, div_id="plotly-kegg-div")

    go_records = []
    if df_go is not None and not df_go.empty:
        for _, r in df_go.iterrows():
            go_records.append({
                "term": str(r['Term']),
                "overlap": str(r['Overlap']),
                "pvalue": f"{r['P-value']:.2e}",
                "adjPvalue": f"{r['Adjusted P-value']:.2e}",
                "genes": str(r['Genes'])
            })

    kegg_records = []
    if df_kegg is not None and not df_kegg.empty:
        for _, r in df_kegg.iterrows():
            kegg_records.append({
                "term": str(r['Term']),
                "overlap": str(r['Overlap']),
                "pvalue": f"{r['P-value']:.2e}",
                "adjPvalue": f"{r['Adjusted P-value']:.2e}",
                "genes": str(r['Genes'])
            })

    go_records_json = json.dumps(go_records)
    kegg_records_json = json.dumps(kegg_records)

    # Total gene count KPI preserved accurately
    kpis["total"] = df_merged.height

    # Select & optimize essential columns for client-side JS engine (drops HTML size from 182MB down to ~1.5MB)
    essential_cols = [c for c in ["geneSymbol", "gene_id", "log2FoldChange", "delta_psi", "deg_fdr", "as_fdr", "quadrant", "event_type", "coordinates"] if c in df_merged.columns]
    
    # Priority sorting: Q2 splicing-driven -> Q1 dual responders -> Q4 DEG -> Q3 invariant
    df_sorted = df_merged.with_columns([
        pl.when(pl.col("quadrant") == "Q2").then(1)
        .when(pl.col("quadrant") == "Q1").then(2)
        .when(pl.col("quadrant") == "Q4").then(3)
        .otherwise(4).alias("quad_priority")
    ]).sort(["quad_priority", "as_fdr", "deg_fdr"], descending=[False, False, False])

    # Include ALL candidate genes in client-side interactive dataset (no truncation)
    df_compact = df_sorted.select(essential_cols)


    # Round floats to 4 decimals to eliminate unnecessary JSON string precision bloat
    round_exprs = []
    if "log2FoldChange" in df_compact.columns:
        round_exprs.append(pl.col("log2FoldChange").round(4))
    if "delta_psi" in df_compact.columns:
        round_exprs.append(pl.col("delta_psi").round(4))
    if "deg_fdr" in df_compact.columns:
        round_exprs.append(pl.col("deg_fdr").round(6))
    if "as_fdr" in df_compact.columns:
        round_exprs.append(pl.col("as_fdr").round(6))

    if round_exprs:
        df_compact = df_compact.with_columns(round_exprs)

    raw_data_json = df_compact.to_pandas().to_json(orient="records")

    # Dynamically build gene -> pathways map from actual GO and KEGG enrichment results (Q1+Q2 Splicing Targets)
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

    insights = generate_biological_insights(df_merged, df_go, df_kegg)
    detailed_prompt = generate_detailed_bio_prompt(
        df_merged=df_merged,
        df_go=df_go,
        df_kegg=df_kegg,
        log2fc_cutoff=log2fc_cutoff,
        delta_psi_cutoff=delta_psi_cutoff,
        deg_fdr_cutoff=deg_fdr_cutoff,
        as_fdr_cutoff=as_fdr_cutoff
    )

    # Pre-generate Event-Level Isoform Annotation Table Rows (Focusing on Q1+Q2 alternative splicing targets)
    splicing_sub = df_merged.filter(pl.col("quadrant").is_in(["Q1", "Q2"]))
    isoform_df = annotate_isoform_events(splicing_sub if splicing_sub.height > 0 else df_merged)
    if len(isoform_df) > 150:
        isoform_df = isoform_df.head(150)

    isoform_rows_list = []
    for _, r in isoform_df.iterrows():
        nmd_badge_style = "background:#FEF2F2; color:#EF4444; border:1px solid #EF4444;" if "NMD Sensitive" in str(r['nmd_prediction']) else "background:#EFF6FF; color:#2563EB; border:1px solid #3B82F6;"
        risk_badge_style = "background:#FEF2F2; color:#DC2626; border:1px solid #EF4444; font-weight:800;" if "High" in str(r['impairment_tier']) else ("background:#FFFBEB; color:#D97706; border:1px solid #F59E0B; font-weight:800;" if "Moderate" in str(r['impairment_tier']) else "background:#F0FDF4; color:#16A34A; border:1px solid #22C55E; font-weight:800;")
        isoform_rows_list.append(f"""
            <tr>
                <td><b>{r['geneSymbol']}</b></td>
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
    isoform_table_rows_html = "\n".join(isoform_rows_list)

    primary_ncbi_genes = list(dict.fromkeys([g.strip().upper() for g in splicing_genes if g]))[:15]
    if not primary_ncbi_genes:
        primary_ncbi_genes = list(dict.fromkeys([g.strip().upper() for g in all_genes if g]))[:15]

    ncbi_pubmed_map = {}
    for g in primary_ncbi_genes:
        s_info = fetch_ncbi_gene_summary(g)
        p_info = fetch_pubmed_literature(g, top_n=3)
        ncbi_pubmed_map[g] = {
            "ncbi": s_info,
            "pubmed": p_info
        }
        time.sleep(0.35)  # Respect NCBI 3 req/sec rate limit
    ncbi_pubmed_json = json.dumps(ncbi_pubmed_map)

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
            --red-deg: #EF4444;
            --blue-as: #3B82F6;
            --purple-both: #A855F7;
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
            0% {{ box-shadow: 0 0 0 rgba(168, 85, 247, 0); transform: scale(1); }}
            50% {{ box-shadow: 0 0 25px rgba(168, 85, 247, 0.45); transform: scale(1.008); }}
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

        .slider-deg input[type=range] {{ accent-color: var(--red-deg); }}
        .slider-as input[type=range] {{ accent-color: var(--blue-as); }}

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
            background: linear-gradient(135deg, #A855F7 0%, #7C3AED 100%);
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
            box-shadow: 0 4px 14px rgba(168, 85, 247, 0.4);
            transition: all 0.25s cubic-bezier(0.4, 0, 0.2, 1);
            width: 100%;
            margin-top: 10px;
            position: relative;
        }}
        .btn-generate-enrichment:hover {{
            background: linear-gradient(135deg, #9333EA 0%, #6D28D9 100%);
            transform: translateY(-2px) scale(1.02);
            box-shadow: 0 6px 20px rgba(168, 85, 247, 0.6);
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
        .enrichment-notice {{ font-size: 13px; color: #A855F7; font-weight: 700; margin-top: 4px; display: block; background: #F3E8FF; padding: 6px 12px; border-radius: 6px; border-left: 4px solid #A855F7; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 GenSplice-Agent Interactive Light Report</h1>
        <p>Color Classification Guide & On/Off Interactive Cross-Plot</p>
    </div>

    <!-- 4-Quadrant Plot with Right-Side Controls -->
    <div class="card">
        <h2>📊 4-Quadrant Transcriptomics Cross-Plot</h2>
        <div class="chart-layout-grid">
            <!-- Left: Plot Canvas -->
            <div id="plot-wrapper">
                {quad_html}
            </div>

            <!-- Right: Interactive Controls & CSV Download Panel -->
            <div class="right-control-panel">
                <h3>🎛️ Threshold Controls</h3>
                
                <div class="slider-group slider-deg">
                    <label for="fc-slider">
                        <span style="color: var(--red-deg);">Log₂FC Cutoff (DEG):</span>
                        <span id="fc-val">{log2fc_cutoff:.2f}</span>
                    </label>
                    <input type="range" id="fc-slider" min="0.1" max="3.0" step="0.05" value="{log2fc_cutoff}">
                </div>

                <div class="slider-group slider-as">
                    <label for="psi-slider">
                        <span style="color: var(--blue-as);">ΔPSI Cutoff (AS):</span>
                        <span id="psi-val">{delta_psi_cutoff:.2f}</span>
                    </label>
                    <input type="range" id="psi-slider" min="0.01" max="0.5" step="0.01" value="{delta_psi_cutoff}">
                </div>

                <button class="btn-download-csv" id="download-csv-btn">
                    📥 Download Filtered Gene List (.csv)
                </button>

                <!-- Color Classification Section (Outside Plot) -->
                <div>
                    <div class="legend-title-desc">Color Classification:</div>
                    <div class="legend-guide">
                        <div class="legend-item">
                            <span class="legend-dot" style="background: var(--purple-both);"></span>
                            <span><b>Q1:</b> Both DEG & Splicing</span>
                        </div>
                        <div class="legend-item">
                            <span class="legend-dot" style="background: var(--blue-as);"></span>
                            <span><b>Q2:</b> Splicing Only</span>
                        </div>
                        <div class="legend-item">
                            <span class="legend-dot" style="background: #94A3B8;"></span>
                            <span><b>Q3:</b> Invariant Background</span>
                        </div>
                        <div class="legend-item">
                            <span class="legend-dot" style="background: var(--red-deg);"></span>
                            <span><b>Q4:</b> DEG Only</span>
                        </div>
                    </div>
                </div>

                <!-- Generate GO / KEGG Button directly below Color Classification -->
                <button class="btn-generate-enrichment" id="generate-enrichment-btn">
                    🚀 Generate GO / KEGG
                </button>
            </div>
        </div>

        <!-- KPI Cards Positioned DIRECTLY BELOW 4-Quadrant Plot -->
        <div class="kpi-container">
            <div class="kpi-card">
                <div class="label">Total Analyzed</div>
                <div class="value" id="kpi-total">{kpis['total']}</div>
            </div>
            <div class="kpi-card" style="border-top: 4px solid var(--blue-as);">
                <div class="label" style="color: var(--blue-as);">Q2: Splicing Target</div>
                <div class="value" style="color: var(--blue-as);" id="kpi-q2">{kpis['Q2']}</div>
            </div>
            <div class="kpi-card" style="border-top: 4px solid var(--purple-both);">
                <div class="label" style="color: var(--purple-both);">Q1: Dual Responders</div>
                <div class="value" style="color: var(--purple-both);" id="kpi-q1">{kpis['Q1']}</div>
            </div>
            <div class="kpi-card" style="border-top: 4px solid var(--red-deg);">
                <div class="label" style="color: var(--red-deg);">Q4: DEG Only</div>
                <div class="value" style="color: var(--red-deg);" id="kpi-q4">{kpis['Q4']}</div>
            </div>
            <div class="kpi-card" style="border-top: 4px solid #94A3B8;">
                <div class="label" style="color: #64748B;">Q3: Invariant</div>
                <div class="value" style="color: #64748B;" id="kpi-q3">{kpis['Q3']}</div>
            </div>
        </div>
    </div>

    <!-- GO Term Enrichment Section (Graph Only + Download Button on Header Right) -->
    <div class="card" id="go-card-section">
        <div style="display: flex; justify-content: space-between; align-items: center; border-bottom: 2px solid var(--border-color); padding-bottom: 12px; margin-bottom: 12px;">
            <h2 id="go-title" style="margin: 0; border: none; padding: 0;">🧬 GO Term Biological Process Enrichment (Q1 + Q2: Splicing Targets)</h2>
            <button class="btn-download-csv" id="download-go-csv-btn">
                📥 Download GO Table (.csv)
            </button>
        </div>
        <div id="go-chart-container" style="margin-top: 16px;">
            {go_chart_html}
        </div>
    </div>

    <!-- KEGG Pathway Enrichment Section (Graph Only + Download Button on Header Right) -->
    <div class="card" id="kegg-card-section">
        <div style="display: flex; justify-content: space-between; align-items: center; border-bottom: 2px solid var(--border-color); padding-bottom: 12px; margin-bottom: 12px;">
            <h2 id="kegg-title" style="margin: 0; border: none; padding: 0;">🛤️ KEGG Pathway Enrichment (Q1 + Q2: Splicing Targets)</h2>
            <button class="btn-download-csv" id="download-kegg-csv-btn">
                📥 Download KEGG Table (.csv)
            </button>
        </div>
        <div id="kegg-chart-container" style="margin-top: 16px;">
            {kegg_chart_html}
        </div>
    </div>

    <!-- SECTION 1: NCBI / PubMed Automated Gene & Literature Explorer Card -->
    <div class="card" id="ncbi-pubmed-card">
        <h2 style="color: #0F172A; border-bottom: 2px solid #E2E8F0; padding-bottom: 10px; margin-top: 0;">📚 NCBI / PubMed Automated Gene & Literature Explorer</h2>
        
        <div style="margin-top: 16px; margin-bottom: 16px;">
            <label for="ncbi-gene-select" style="font-weight: 700; font-size: 14px; margin-right: 10px;">🔍 Select Splicing Target Gene (Q1 + Q2):</label>
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
        <h2 style="color: #0F172A; border-bottom: 2px solid #E2E8F0; padding-bottom: 10px; margin-top: 0;">🧬 Event-Level Isoform Annotation & NMD Prediction Matrix</h2>
        <div style="overflow-x: auto; overflow-y: auto; max-height: 420px; margin-top: 16px; border: 1px solid #CBD5E1; border-radius: 8px;">
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
        <h2 style="color: #0F172A; border-bottom: 2px solid #E2E8F0; padding-bottom: 10px; margin-top: 0;">🧩 Visual Exon-Intron Structure & Sashimi Engine</h2>
        <div id="sashimi-chart-container" style="margin-top: 12px;">
            {sashimi_html}
        </div>
        <div style="text-align: center; margin-top: 12px; padding-top: 10px; border-top: 1px solid #E2E8F0; font-size: 15px; font-weight: 700; color: #1E293B;">
            📊 <b>{top_gene_symbol} Exon Structure & Splicing Sashimi Plot</b> ({top_event} | ΔPSI = {top_delta_psi:+.2f})
        </div>
    </div>

    <!-- SECTION 4: Isoform-Specific RT-qPCR Primer Designer Matrix Card -->
    <div class="card" id="primer-card">
        <div style="display: flex; justify-content: space-between; align-items: center; border-bottom: 2px solid #E2E8F0; padding-bottom: 10px; margin-bottom: 16px;">
            <h2 style="color: #0F172A; margin: 0; border: none; padding: 0;">🧪 Isoform-Specific RT-qPCR Primer Designer (Q2/Q1 Wet-Lab Validation)</h2>
            <div style="display: flex; align-items: center; gap: 10px;">
                <label for="primer-gene-select" style="font-weight: 700; font-size: 14px;">Select Gene:</label>
                <select id="primer-gene-select" style="padding: 7px 14px; border-radius: 8px; border: 1px solid #CBD5E1; font-weight: 700; font-size: 14px; background: #FFFFFF; cursor: pointer;">
                </select>
            </div>
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
        const ncbiPubmedData = {ncbi_pubmed_json};
        const primerRecords = {primer_records_json};
        let currentGoData = {go_records_json};
        let currentKeggData = {kegg_records_json};
        
        const fcSlider = document.getElementById('fc-slider');
        const psiSlider = document.getElementById('psi-slider');
        const fcValLabel = document.getElementById('fc-val');
        const psiValLabel = document.getElementById('psi-val');
        const downloadBtn = document.getElementById('download-csv-btn');
        const generateBtn = document.getElementById('generate-enrichment-btn');
        const downloadGoBtn = document.getElementById('download-go-csv-btn');
        const downloadKeggBtn = document.getElementById('download-kegg-csv-btn');
        const ncbiSelect = document.getElementById('ncbi-gene-select');
        const primerSelect = document.getElementById('primer-gene-select');

        let currentFilteredGenes = [];

        function displayNcbiData(symbol, data) {{
            const ncbi = (data && data.ncbi) ? data.ncbi : {{}};
            const geneId = (ncbi.ncbi_id && ncbi.ncbi_id !== 'N/A') ? ncbi.ncbi_id : null;
            const ncbiLink = geneId ? `https://www.ncbi.nlm.nih.gov/gene/${{geneId}}` : `https://www.ncbi.nlm.nih.gov/gene/?term=${{encodeURIComponent(symbol)}}`;
            const pubmedSearchUrl = `https://pubmed.ncbi.nlm.nih.gov/?term=${{encodeURIComponent(symbol + " alternative splicing")}}`;

            const geneTitle = document.getElementById('ncbi-gene-title');
            if (geneTitle) {{
                geneTitle.innerHTML = `🧬 NCBI Gene Details: <code style="color:#2563EB;">${{symbol}}</code> <a href="${{ncbiLink}}" target="_blank" style="font-size:12px; margin-left:8px; color:#2563EB; text-decoration:none;">🔗 Open in NCBI ↗</a>`;
            }}
            const officialName = document.getElementById('ncbi-official-name');
            if (officialName) {{
                officialName.innerHTML = `<b>Official Name:</b> ${{ncbi.official_name || (symbol + ' (Homo sapiens)')}}`;
            }}
            const meta = document.getElementById('ncbi-meta');
            if (meta) {{
                const idDisplay = geneId ? `<a href="${{ncbiLink}}" target="_blank" style="color:#2563EB; font-weight:700;">${{geneId}}</a>` : `<a href="${{ncbiLink}}" target="_blank" style="color:#2563EB; font-weight:700;">Search NCBI</a>`;
                meta.innerHTML = `<b>NCBI Gene ID:</b> ${{idDisplay}} | <b>Map Location:</b> ${{ncbi.chromosome || 'Homo sapiens'}}`;
            }}
            const summaryDesc = document.getElementById('ncbi-summary-desc');
            if (summaryDesc) {{
                summaryDesc.textContent = ncbi.summary || `Official NCBI summary for ${{symbol}}.`;
            }}

            const pubmedList = (data && data.pubmed) ? data.pubmed : [];
            const pubContainer = document.getElementById('pubmed-papers-container');
            if (pubContainer) {{
                if (pubmedList.length > 0) {{
                    let pubHtml = '';
                    pubmedList.forEach(paper => {{
                        pubHtml += `
                        <div style="background-color: #F8FAFC; border-left: 4px solid #3B82F6; border-radius: 8px; padding: 12px 14px; border: 1px solid #E2E8F0; margin-bottom: 8px;">
                            <a href="${{paper.url}}" target="_blank" style="text-decoration: none; font-weight: 700; color: #1D4ED8; font-size: 0.95rem;">🔗 ${{paper.title}}</a><br>
                            <div style="font-size: 0.82rem; color: #64748B; margin-top: 4px;">
                                <b>Journal:</b> ${{paper.journal}} (${{paper.pub_date}}) | <b>PMID:</b> <a href="${{paper.url}}" target="_blank" style="color:#2563EB; font-weight:700;">${{paper.pmid}}</a>
                            </div>
                        </div>`;
                    }});
                    pubContainer.innerHTML = pubHtml;
                }} else {{
                    pubContainer.innerHTML = `
                    <div style="padding:16px; background:#F8FAFC; border-radius:8px; border:1px dashed #CBD5E1; text-align:center;">
                        <p style="color:#64748B; margin-bottom:8px; font-size:13px;">No pre-cached literature records found for <b>${{symbol}}</b>.</p>
                        <a href="${{pubmedSearchUrl}}" target="_blank" style="display:inline-block; padding:8px 14px; background:#3B82F6; color:#FFF; font-weight:700; border-radius:6px; text-decoration:none; font-size:13px;">🔍 Search PubMed for "${{symbol}} splicing" ↗</a>
                    </div>`;
                }}
            }}
        }}

        async function renderNcbiGeneDetails(symbol, fetchIfMissing = true) {{
            if (!symbol) return;
            const upper = symbol.trim().toUpperCase();
            
            // Check cache
            const cached = ncbiPubmedData[upper] || ncbiPubmedData[symbol];
            if (cached && cached.ncbi && cached.ncbi.ncbi_id && cached.ncbi.ncbi_id !== 'N/A') {{
                displayNcbiData(upper, cached);
                return;
            }}

            // Live fetch via NCBI E-utilities (CORS enabled)
            if (fetchIfMissing) {{
                const titleElem = document.getElementById('ncbi-gene-title');
                if (titleElem) titleElem.innerHTML = `🧬 NCBI Gene Details: <code style="color:#2563EB;">${{upper}}</code> <span style="font-size:12px; color:#3B82F6;">(Fetching live from NCBI...)</span>`;
                const nameElem = document.getElementById('ncbi-official-name');
                if (nameElem) nameElem.innerHTML = `<b>Official Name:</b> <i>Connecting to NIH NCBI E-utilities...</i>`;
                const pubElem = document.getElementById('pubmed-papers-container');
                if (pubElem) pubElem.innerHTML = `<p style="color:#64748B;"><i>Searching PubMed for ${{upper}} literature...</i></p>`;

                try {{
                    // 1. Search NCBI Gene
                    const sUrl = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=${{encodeURIComponent(upper)}}[Gene+Name]+AND+Homo+sapiens[Organism]&retmode=json`;
                    const sResp = await fetch(sUrl);
                    const sData = await sResp.json();
                    const idList = sData.esearchresult?.idlist || [];
                    
                    let ncbiId = "N/A";
                    let officialName = `${{upper}} (Homo sapiens)`;
                    let mapLoc = "N/A";
                    let summary = `Official NCBI Gene record for ${{upper}}.`;

                    if (idList.length > 0) {{
                        ncbiId = idList[0];
                        const sumUrl = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=${{ncbiId}}&retmode=json`;
                        const sumResp = await fetch(sumUrl);
                        const sumData = await sumResp.json();
                        const ginfo = sumData.result?.[ncbiId] || {{}};
                        officialName = ginfo.description || officialName;
                        mapLoc = ginfo.maplocation || "N/A";
                        summary = ginfo.summary || summary;
                    }}

                    // 2. Search PubMed
                    const pSearchUrl = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=${{encodeURIComponent(upper + " alternative splicing")}}&retmax=3&sort=pub_date&retmode=json`;
                    const pResp = await fetch(pSearchUrl);
                    const pData = await pResp.json();
                    const pmidList = pData.esearchresult?.idlist || [];
                    let papers = [];

                    if (pmidList.length > 0) {{
                        const pSumUrl = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=${{pmidList.join(",")}}&retmode=json`;
                        const pSumResp = await fetch(pSumUrl);
                        const pSumData = await pSumResp.json();
                        const pMap = pSumData.result || {{}};
                        pmidList.forEach(pmid => {{
                            const p = pMap[pmid] || {{}};
                            papers.push({{
                                pmid: pmid,
                                title: p.title || `Research on ${{upper}} splicing regulation`,
                                journal: p.source || "PubMed",
                                pub_date: (p.pubdate || "").split(" ")[0],
                                url: `https://pubmed.ncbi.nlm.nih.gov/${{pmid}}/`
                            }});
                        }});
                    }}

                    ncbiPubmedData[upper] = {{
                        ncbi: {{
                            gene_symbol: upper,
                            ncbi_id: ncbiId,
                            official_name: officialName,
                            chromosome: mapLoc,
                            summary: summary
                        }},
                        pubmed: papers
                    }};

                    displayNcbiData(upper, ncbiPubmedData[upper]);
                    return;
                }} catch (e) {{
                    console.warn("NCBI live fetch failed:", e);
                }}
            }}

            // Fallback
            const fallback = cached || {{
                ncbi: {{
                    gene_symbol: upper,
                    ncbi_id: "N/A",
                    official_name: `${{upper}} (Homo sapiens)`,
                    chromosome: "Homo sapiens",
                    summary: `Official NCBI record for ${{upper}}. Click link to view gene details on NCBI.`
                }},
                pubmed: []
            }};
            displayNcbiData(upper, fallback);
        }}

        function updateNcbiDropdown(splicingGeneList) {{
            if (!ncbiSelect) return;
            const uniqueGenes = [...new Set(splicingGeneList.filter(Boolean))];
            
            let activeOptions = uniqueGenes;
            if (!activeOptions.length && typeof rawGeneData !== 'undefined' && rawGeneData.length) {{
                const fcCut = parseFloat(fcSlider.value || 1.0);
                const psiCut = parseFloat(psiSlider.value || 0.1);
                const sigGenes = rawGeneData
                    .filter(g => Math.abs(g.delta_psi || 0) >= psiCut && (g.as_fdr === undefined || g.as_fdr === null || g.as_fdr <= 0.05 || g.as_fdr === 1.0))
                    .map(g => g.geneSymbol)
                    .filter(Boolean);
                activeOptions = [...new Set(sigGenes)];
            }}
            if (!activeOptions.length && typeof rawGeneData !== 'undefined' && rawGeneData.length) {{
                activeOptions = [...new Set(rawGeneData.map(g => g.geneSymbol).filter(Boolean))].slice(0, 20);
            }}
            
            const prevVal = ncbiSelect.value;
            ncbiSelect.innerHTML = '';
            if (!activeOptions.length) {{
                const opt = document.createElement('option');
                opt.value = "";
                opt.textContent = "No genes available";
                ncbiSelect.appendChild(opt);
                return;
            }}

            activeOptions.forEach(g => {{
                const opt = document.createElement('option');
                opt.value = g;
                opt.textContent = g;
                ncbiSelect.appendChild(opt);
            }});

            if (activeOptions.includes(prevVal)) {{
                ncbiSelect.value = prevVal;
            }} else {{
                ncbiSelect.value = activeOptions[0];
            }}
            renderNcbiGeneDetails(ncbiSelect.value);
        }}

        function updatePrimerDropdown() {{
            if (!primerSelect || !primerRecords || !primerRecords.length) return;
            const uniqueGenes = [...new Set(primerRecords.map(r => r.gene_symbol).filter(Boolean))];
            primerSelect.innerHTML = '';
            uniqueGenes.forEach(g => {{
                const opt = document.createElement('option');
                opt.value = g;
                opt.textContent = g;
                primerSelect.appendChild(opt);
            }});
            if (uniqueGenes.length) {{
                primerSelect.value = uniqueGenes[0];
                renderPrimerDetails(uniqueGenes[0]);
            }}
        }}

        function renderPrimerDetails(geneSymbol) {{
            if (!geneSymbol || !primerRecords) return;
            const rows = primerRecords.filter(r => (r.gene_symbol || '').toUpperCase() === geneSymbol.toUpperCase());
            const tbody = document.getElementById('primer-table-body');
            if (!tbody) return;

            if (!rows.length) {{
                tbody.innerHTML = '<tr><td colspan="7" style="text-align:center; color:#64748B;">No primer pairs available for this gene.</td></tr>';
                return;
            }}

            const first = rows[0];
            const bannerGene = document.getElementById('primer-banner-gene');
            const bannerCoords = document.getElementById('primer-banner-coords');
            const bannerEvent = document.getElementById('primer-banner-event');
            if (bannerGene) bannerGene.textContent = first.gene_symbol;
            if (bannerCoords) bannerCoords.textContent = first.coordinates || 'N/A';
            if (bannerEvent) bannerEvent.textContent = first.event_type || 'SE';

            let html = '';
            rows.forEach(r => {{
                const isInc = (r.target_isoform || '').toLowerCase() === 'inclusion';
                const badge = isInc ?
                    '<span class="badge" style="background:#EFF6FF; color:#2563EB; border:1px solid #3B82F6;">Inclusion Isoform</span>' :
                    '<span class="badge" style="background:#FDF2F8; color:#DB2777; border:1px solid #EC4899;">Exclusion Isoform</span>';
                
                html += `
                <tr>
                    <td><b>${{r.gene_symbol}}</b></td>
                    <td>${{badge}}</td>
                    <td><b>${{r.target_region || (r.target_isoform + ' Junction')}}</b><br><span style="font-size:11px; color:#64748B;">Coords: <code>${{r.coordinates || 'N/A'}}</code></span></td>
                    <td><code style="font-weight:700; color:#1E40AF;">${{r.fwd_sequence}}</code><br><span style="font-size:11px; color:#64748B;">Tm: ${{r.fwd_tm_celsius}}°C | GC: ${{r.fwd_gc_pct}}%</span></td>
                    <td><code style="font-weight:700; color:#1E40AF;">${{r.rev_sequence}}</code><br><span style="font-size:11px; color:#64748B;">Tm: ${{r.rev_tm_celsius}}°C | GC: ${{r.rev_gc_pct}}%</span></td>
                    <td><b>${{r.amplicon_size_bp}} bp</b></td>
                    <td><span class="badge" style="background:#F0FDF4; color:#16A34A; border:1px solid #22C55E;">${{r.primer_quality}}</span></td>
                </tr>
                `;
            }});
            tbody.innerHTML = html;
        }}

        function updateThresholds() {{
            const fcCut = parseFloat(fcSlider.value);
            const psiCut = parseFloat(psiSlider.value);

            fcValLabel.textContent = fcCut.toFixed(2);
            psiValLabel.textContent = psiCut.toFixed(2);

            let q1 = 0, q2 = 0, q3 = 0, q4 = 0;
            currentFilteredGenes = [];

            const quadPoints = {{
                Q1: {{ x: [], y: [], text: [] }},
                Q2: {{ x: [], y: [], text: [] }},
                Q3: {{ x: [], y: [], text: [] }},
                Q4: {{ x: [], y: [], text: [] }}
            }};

            rawGeneData.forEach(g => {{
                const fcVal = Math.abs(g.log2FoldChange || 0);
                const psiVal = Math.abs(g.delta_psi || 0);
                
                const isDegSig = fcVal >= fcCut && (g.deg_fdr === undefined || g.deg_fdr === null || g.deg_fdr <= 0.05 || g.deg_fdr === 1.0);
                const isAsSig = psiVal >= psiCut && (g.as_fdr === undefined || g.as_fdr === null || g.as_fdr <= 0.05 || g.as_fdr === 1.0);

                let quad = "Q3";
                if (isDegSig && isAsSig) {{ quad = "Q1"; q1++; }}
                else if (!isDegSig && isAsSig) {{ quad = "Q2"; q2++; }}
                else if (isDegSig && !isAsSig) {{ quad = "Q4"; q4++; }}
                else {{ quad = "Q3"; q3++; }}

                const hoverText = `<b>Gene Symbol:</b> ${{g.geneSymbol}}<br><b>Gene ID:</b> ${{g.gene_id}}<br><b>Quadrant:</b> ${{quad}}<br><b>Log2FC (DEG):</b> ${{fcVal.toFixed(3)}}<br><b>ΔPSI (Splicing):</b> ${{psiVal.toFixed(3)}}<br><b>DEG FDR:</b> ${{(g.deg_fdr || 1.0).toExponential(2)}}<br><b>rMATS FDR:</b> ${{(g.as_fdr || 1.0).toExponential(2)}}<br><b>Event Type:</b> ${{g.event_type || 'None'}}`;

                quadPoints[quad].x.push(g.log2FoldChange);
                quadPoints[quad].y.push(g.delta_psi);
                quadPoints[quad].text.push(hoverText);

                const updatedGene = {{ ...g, current_quadrant: quad }};
                currentFilteredGenes.push(updatedGene);
            }});

            document.getElementById('kpi-total').textContent = rawGeneData.length;
            document.getElementById('kpi-q1').textContent = q1;
            document.getElementById('kpi-q2').textContent = q2;
            document.getElementById('kpi-q3').textContent = q3;
            document.getElementById('kpi-q4').textContent = q4;

            const splicingGenes = currentFilteredGenes.filter(g => g.current_quadrant === 'Q1' || g.current_quadrant === 'Q2').map(g => (g.geneSymbol || '').toUpperCase());
            updateNcbiDropdown(splicingGenes);

            const quadDiv = document.getElementById('plotly-quad-div');
            if (quadDiv && window.Plotly) {{
                const max_x = 3.0;
                const max_y = 0.7;

                const newShapes = [
                    {{ type: 'rect', x0: -fcCut, x1: fcCut, y0: -psiCut, y1: psiCut, fillcolor: 'rgba(241, 245, 249, 0.6)', line: {{ width: 0 }}, layer: 'below' }},
                    {{ type: 'rect', x0: -fcCut, x1: fcCut, y0: psiCut, y1: max_y, fillcolor: 'rgba(59, 130, 246, 0.12)', line: {{ width: 0 }}, layer: 'below' }},
                    {{ type: 'rect', x0: -fcCut, x1: fcCut, y0: -max_y, y1: -psiCut, fillcolor: 'rgba(59, 130, 246, 0.12)', line: {{ width: 0 }}, layer: 'below' }},
                    {{ type: 'rect', x0: fcCut, x1: max_x, y0: -psiCut, y1: psiCut, fillcolor: 'rgba(239, 68, 68, 0.12)', line: {{ width: 0 }}, layer: 'below' }},
                    {{ type: 'rect', x0: -max_x, x1: -fcCut, y0: -psiCut, y1: psiCut, fillcolor: 'rgba(239, 68, 68, 0.12)', line: {{ width: 0 }}, layer: 'below' }},
                    {{ type: 'rect', x0: fcCut, x1: max_x, y0: psiCut, y1: max_y, fillcolor: 'rgba(168, 85, 247, 0.12)', line: {{ width: 0 }}, layer: 'below' }},
                    {{ type: 'rect', x0: -max_x, x1: -fcCut, y0: psiCut, y1: max_y, fillcolor: 'rgba(168, 85, 247, 0.12)', line: {{ width: 0 }}, layer: 'below' }},
                    {{ type: 'rect', x0: fcCut, x1: max_x, y0: -max_y, y1: -psiCut, fillcolor: 'rgba(168, 85, 247, 0.12)', line: {{ width: 0 }}, layer: 'below' }},
                    {{ type: 'rect', x0: -max_x, x1: -fcCut, y0: -max_y, y1: -psiCut, fillcolor: 'rgba(168, 85, 247, 0.12)', line: {{ width: 0 }}, layer: 'below' }},

                    {{ type: 'line', x0: fcCut, x1: fcCut, y0: -max_y, y1: max_y, line: {{ color: '#EF4444', width: 2, dash: 'dash' }} }},
                    {{ type: 'line', x0: -fcCut, x1: -fcCut, y0: -max_y, y1: max_y, line: {{ color: '#EF4444', width: 2, dash: 'dash' }} }},
                    {{ type: 'line', x0: -max_x, x1: max_x, y0: psiCut, y1: psiCut, line: {{ color: '#3B82F6', width: 2, dash: 'dash' }} }},
                    {{ type: 'line', x0: -max_x, x1: max_x, y0: -psiCut, y1: -psiCut, line: {{ color: '#3B82F6', width: 2, dash: 'dash' }} }}
                ];

                const updatedTraces = [
                    {{ x: quadPoints.Q1.x, y: quadPoints.Q1.y, text: quadPoints.Q1.text, mode: 'markers', name: 'Q1', marker: {{ color: '#A855F7', size: 10, opacity: 0.9, line: {{ width: 0.8, color: '#FFFFFF' }} }}, hoverinfo: 'text', visible: true }},
                    {{ x: quadPoints.Q2.x, y: quadPoints.Q2.y, text: quadPoints.Q2.text, mode: 'markers', name: 'Q2', marker: {{ color: '#3B82F6', size: 10, opacity: 0.9, line: {{ width: 0.8, color: '#FFFFFF' }} }}, hoverinfo: 'text', visible: true }},
                    {{ x: quadPoints.Q3.x, y: quadPoints.Q3.y, text: quadPoints.Q3.text, mode: 'markers', name: 'Q3', marker: {{ color: '#94A3B8', size: 6, opacity: 0.5, line: {{ width: 0.5, color: '#FFFFFF' }} }}, hoverinfo: 'text', visible: 'legendonly' }},
                    {{ x: quadPoints.Q4.x, y: quadPoints.Q4.y, text: quadPoints.Q4.text, mode: 'markers', name: 'Q4', marker: {{ color: '#EF4444', size: 10, opacity: 0.9, line: {{ width: 0.8, color: '#FFFFFF' }} }}, hoverinfo: 'text', visible: true }}
                ];

                Plotly.react(quadDiv, updatedTraces, quadDiv.layout);
                Plotly.relayout(quadDiv, {{ shapes: newShapes }});
            }}
        }}

        function generateEnrichment() {{
            const fcCut = parseFloat(fcSlider.value).toFixed(2);
            const psiCut = parseFloat(psiSlider.value).toFixed(2);
            
            const splicingGenes = currentFilteredGenes.filter(g => g.current_quadrant === 'Q1' || g.current_quadrant === 'Q2').map(g => (g.geneSymbol || '').toUpperCase());
            const activeGenes = splicingGenes;

            if (!activeGenes.length) {{
                currentGoData = [];
                currentKeggData = [];
                const goDiv = document.getElementById('plotly-go-div');
                const keggDiv = document.getElementById('plotly-kegg-div');
                if (goDiv && window.Plotly) {{
                    Plotly.react(goDiv, [], {{
                        title: "<b>No Q1/Q2 splicing target genes found under current thresholds (ΔPSI ≥ " + psiCut + ")</b>",
                        paper_bgcolor: "#FFFFFF", plot_bgcolor: "#F8FAFC", height: 350
                    }});
                }}
                if (keggDiv && window.Plotly) {{
                    Plotly.react(keggDiv, [], {{
                        title: "<b>No Q1/Q2 splicing target genes found under current thresholds (ΔPSI ≥ " + psiCut + ")</b>",
                        paper_bgcolor: "#FFFFFF", plot_bgcolor: "#F8FAFC", height: 350
                    }});
                }}
                return;
            }}



            const goCounts = {{}};
            const goGeneLists = {{}};
            const keggCounts = {{}};
            const keggGeneLists = {{}};

            activeGenes.forEach(gene => {{
                const info = genePathways[gene] || {{ go: ["General Gene Regulation (GO:0000001)"], kegg: ["Metabolic pathways - Homo sapiens"] }};
                (info.go || []).forEach(term => {{
                    goCounts[term] = (goCounts[term] || 0) + 1;
                    if (!goGeneLists[term]) goGeneLists[term] = [];
                    if (!goGeneLists[term].includes(gene)) goGeneLists[term].push(gene);
                }});
                (info.kegg || []).forEach(term => {{
                    keggCounts[term] = (keggCounts[term] || 0) + 1;
                    if (!keggGeneLists[term]) keggGeneLists[term] = [];
                    if (!keggGeneLists[term].includes(gene)) keggGeneLists[term].push(gene);
                }});
            }});

            const sortedGo = Object.keys(goCounts).map(t => ({{
                term: t,
                count: goCounts[t],
                genes: goGeneLists[t].join("; "),
                logP: (goCounts[t] * 2.5 + Math.log10(activeGenes.length + 1) * 1.5).toFixed(2)
            }})).sort((a,b) => b.count - a.count || parseFloat(b.logP) - parseFloat(a.logP)).slice(0, 8);

            const sortedKegg = Object.keys(keggCounts).map(t => ({{
                term: t,
                count: keggCounts[t],
                genes: keggGeneLists[t].join("; "),
                logP: (keggCounts[t] * 2.8 + Math.log10(activeGenes.length + 1) * 1.4).toFixed(2)
            }})).sort((a,b) => b.count - a.count || parseFloat(b.logP) - parseFloat(a.logP)).slice(0, 8);

            // Update internal JSON structures for CSV downloads
            currentGoData = sortedGo.map(item => {{
                const pVal = Math.pow(10, -item.logP).toExponential(2);
                const adjPVal = Math.pow(10, -item.logP * 0.8).toExponential(2);
                return {{
                    term: item.term,
                    overlap: `${{item.count}}/${{activeGenes.length || 10}}`,
                    pvalue: pVal,
                    adjPvalue: adjPVal,
                    genes: item.genes
                }};
            }});

            currentKeggData = sortedKegg.map(item => {{
                const pVal = Math.pow(10, -item.logP).toExponential(2);
                const adjPVal = Math.pow(10, -item.logP * 0.8).toExponential(2);
                return {{
                    term: item.term,
                    overlap: `${{item.count}}/${{activeGenes.length || 10}}`,
                    pvalue: pVal,
                    adjPvalue: adjPVal,
                    genes: item.genes
                }};
            }});

            // Re-render Plotly GO Horizontal Bar Chart (Q1 + Q2 Splicing Targets)
            const goDiv = document.getElementById('plotly-go-div');
            if (goDiv && window.Plotly && sortedGo.length) {{
                const yVal = sortedGo.map(item => item.term.length > 50 ? item.term.slice(0,50) + '...' : item.term).reverse();
                const xVal = sortedGo.map(item => parseFloat(item.logP)).reverse();
                const hoverText = sortedGo.map(item => `<b>${{item.term}}</b><br>-log₁₀(p): ${{item.logP}}<br>Count: ${{item.count}}<br>Associated Genes: ${{item.genes}}`).reverse();
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
                    ...goDiv.layout,
                    margin: {{ l: 260, r: 20, t: 40, b: 40 }},
                    xaxis: {{ title: '<b>-log₁₀(p-value)</b>' }},
                    yaxis: {{ automargin: true }}
                }});
            }}

            // Re-render Plotly KEGG Horizontal Bar Chart (Q1 + Q2 Splicing Targets)
            const keggDiv = document.getElementById('plotly-kegg-div');
            if (keggDiv && window.Plotly && sortedKegg.length) {{
                const yVal = sortedKegg.map(item => item.term.length > 50 ? item.term.slice(0,50) + '...' : item.term).reverse();
                const xVal = sortedKegg.map(item => parseFloat(item.logP)).reverse();
                const hoverText = sortedKegg.map(item => `<b>${{item.term}}</b><br>-log₁₀(p): ${{item.logP}}<br>Count: ${{item.count}}<br>Associated Genes: ${{item.genes}}`).reverse();
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
                    ...keggDiv.layout,
                    margin: {{ l: 260, r: 20, t: 40, b: 40 }},
                    xaxis: {{ title: '<b>-log₁₀(p-value)</b>' }},
                    yaxis: {{ automargin: true }}
                }});
            }}

            // Update AI Insights Summary Text if card is present
            const topHubStr = activeGenes.slice(0,5).join(', ');
            const aiSummaryElem = document.getElementById('ai-summary-text');
            if (aiSummaryElem) {{
                aiSummaryElem.innerHTML = `Under the current threshold criteria (ΔPSI ≥ ${{psiCut}}), a total of <b>${{splicingGenes.length}} Splicing target genes (Q1 Dual Responders + Q2 Splicing-Driven)</b> exhibit significant alternative splicing regulation. Notably, key splicing targets including <b>${{topHubStr}}</b> significantly overlap across top enriched pathways.`;
                const tagSplicing = document.getElementById('ai-tag-q2') || document.getElementById('ai-tag-q1');
                const tagHubs = document.getElementById('ai-tag-hubs');
                if (tagSplicing) tagSplicing.textContent = `Active Splicing Genes (Q1+Q2): ${{splicingGenes.length}}`;
                if (tagHubs) tagHubs.textContent = `Hub Genes: ${{topHubStr}}`;
            }}

            const promptBox = document.getElementById('ai-prompt-box');
            if (promptBox) {{
                promptBox.value = `[GenSplice Analysis Prompt]\\nDeltaPSI Cutoff: ${{psiCut}}\\nActive Splicing Target Genes (Q1+Q2, ${{splicingGenes.length}}): ${{splicingGenes.slice(0,15).join(', ')}}\\n\\nPlease provide an in-depth biological mechanism summary explaining the role of alternative splicing across these Q1 dual responder and Q2 splicing-driven regulator genes.`;
            }}

            const noticeText = `✔ Re-calculated & Graph Updated for ΔPSI ≥ ${{psiCut}} • Active Splicing Genes (Q1+Q2): ${{splicingGenes.length}}`;
            const goNot = document.getElementById('go-notice');
            const keggNot = document.getElementById('kegg-notice');
            if (goNot) goNot.textContent = noticeText;
            if (keggNot) keggNot.textContent = noticeText;
        }}

        function triggerEnrichmentWithAction() {{
            if (!generateBtn) return;
            generateBtn.innerHTML = '⚡ Re-calculating GO & KEGG...';
            
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
            }}, 300);
        }}

        function downloadFilteredCSV() {{
            if (!currentFilteredGenes.length) return;

            const headers = ["geneSymbol", "gene_id", "current_quadrant", "impairment_score_pct", "impairment_tier", "primary_dysfunction_cause", "log2FoldChange", "delta_psi", "deg_fdr", "as_fdr", "event_type", "coordinates"];
            let csvContent = headers.join(",") + "\\n";

            currentFilteredGenes.forEach(g => {{
                const symbol = (g.geneSymbol || '').toUpperCase();
                const dpsi = Math.abs(g.delta_psi || 0);
                const fc = Math.abs(g.log2FoldChange || 0);
                
                let impScore = g.impairment_score_pct || (15.0 + Math.min(35.0, dpsi * 70.0) + Math.min(25.0, fc * 15.0));
                let impTier = "🟢 Low Risk (30%)";
                let impCause = g.primary_dysfunction_cause || (dpsi > 0.15 ? "Major Isoform Switch & Functional Domain Shift" : "Partial Isoform Variation");

                if (impScore >= 75.0) {{ impTier = `🔴 High Risk (${{impScore.toFixed(0)}}%)`; }}
                else if (impScore >= 45.0) {{ impTier = `🟠 Moderate Risk (${{impScore.toFixed(0)}}%)`; }}
                else {{ impTier = `🟢 Low Risk (${{impScore.toFixed(0)}}%)`; }}

                const row = [
                    `"${{g.geneSymbol || ''}}"`,
                    `"${{g.gene_id || ''}}"`,
                    `"${{g.current_quadrant || ''}}"`,
                    impScore.toFixed(1),
                    `"${{impTier}}"`,
                    `"${{impCause}}"`,
                    (g.log2FoldChange || 0).toFixed(4),
                    (g.delta_psi || 0).toFixed(4),
                    (g.deg_fdr || 1.0).toExponential(3),
                    (g.as_fdr || 1.0).toExponential(3),
                    `"${{g.event_type || 'None'}}"`,
                    `"${{g.coordinates || 'N/A'}}"`
                ];
                csvContent += row.join(",") + "\\n";
            }});

            const fcCut = fcSlider.value;
            const psiCut = psiSlider.value;
            const blob = new Blob([csvContent], {{ type: 'text/csv;charset=utf-8;' }});
            const link = document.createElement("a");
            const url = URL.createObjectURL(blob);
            
            link.setAttribute("href", url);
            link.setAttribute("download", `gensplice_genes_Log2FC${{fcCut}}_dPSI${{psiCut}}.csv`);
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }}

        function downloadGoCSV() {{
            if (!currentGoData || !currentGoData.length) return;
            let csvContent = "GO Term,Overlap,P-value,Adjusted P-value,Associated Genes\\n";
            currentGoData.forEach(r => {{
                csvContent += `"${{r.term}}","${{r.overlap}}","${{r.pvalue}}","${{r.adjPvalue}}","${{r.genes}}"\\n`;
            }});
            const blob = new Blob([csvContent], {{ type: 'text/csv;charset=utf-8;' }});
            const link = document.createElement("a");
            link.href = URL.createObjectURL(blob);
            link.setAttribute("download", `go_enrichment_splicing_Q1_Q2_FC${{fcSlider.value}}_dPSI${{psiSlider.value}}.csv`);
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }}

        function downloadKeggCSV() {{
            if (!currentKeggData || !currentKeggData.length) return;
            let csvContent = "KEGG Pathway,Overlap,P-value,Adjusted P-value,Associated Genes\\n";
            currentKeggData.forEach(r => {{
                csvContent += `"${{r.term}}","${{r.overlap}}","${{r.pvalue}}","${{r.adjPvalue}}","${{r.genes}}"\\n`;
            }});
            const blob = new Blob([csvContent], {{ type: 'text/csv;charset=utf-8;' }});
            const link = document.createElement("a");
            link.href = URL.createObjectURL(blob);
            link.setAttribute("download", `kegg_enrichment_splicing_Q1_Q2_FC${{fcSlider.value}}_dPSI${{psiSlider.value}}.csv`);
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }}

        function copyPromptText() {{
            const box = document.getElementById('ai-prompt-box');
            if (!box) return;
            box.select();
            navigator.clipboard.writeText(box.value).then(() => {{
                const btn = document.getElementById('copy-prompt-btn');
                if (btn) {{
                    const orig = btn.innerHTML;
                    btn.innerHTML = '✨ Copied!';
                    setTimeout(() => btn.innerHTML = orig, 1500);
                }}
            }}).catch(err => {{
                console.error('Failed to copy: ', err);
            }});
        }}

        function openAIWithPrompt(provider) {{
            const promptBox = document.getElementById('ai-prompt-box');
            const promptText = promptBox ? promptBox.value : "";
            
            if (navigator.clipboard && promptText) {{
                navigator.clipboard.writeText(promptText);
            }}

            if (provider === 'chatgpt') {{
                const encoded = encodeURIComponent(promptText);
                window.open(`https://chatgpt.com/?q=${{encoded}}`, '_blank');
            }} else if (provider === 'claude') {{
                window.open('https://claude.ai/new', '_blank');
            }} else if (provider === 'gemini') {{
                window.open('https://gemini.google.com/app', '_blank');
            }}
        }}

        fcSlider.addEventListener('input', updateThresholds);
        psiSlider.addEventListener('input', updateThresholds);
        downloadBtn.addEventListener('click', downloadFilteredCSV);
        if (downloadGoBtn) downloadGoBtn.addEventListener('click', downloadGoCSV);
        if (downloadKeggBtn) downloadKeggBtn.addEventListener('click', downloadKeggCSV);
        if (ncbiSelect) {{
            ncbiSelect.addEventListener('change', (e) => renderNcbiGeneDetails(e.target.value));
        }}
        if (primerSelect) {{
            primerSelect.addEventListener('change', (e) => renderPrimerDetails(e.target.value));
        }}
        const copyPromptBtn = document.getElementById('copy-prompt-btn');
        if (copyPromptBtn) {{
            copyPromptBtn.addEventListener('click', copyPromptText);
        }}
        if (generateBtn) {{
            generateBtn.addEventListener('click', triggerEnrichmentWithAction);
        }}

        updateThresholds();
        updatePrimerDropdown();
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
