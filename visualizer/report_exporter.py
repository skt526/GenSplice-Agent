"""
GenSplice-Agent Standalone HTML Exporter Module
Generates self-contained, interactive HTML reports openable in Chrome / browsers.
"""

import os
import polars as pl
import pandas as pd
from visualizer.quadrant_plot import build_quadrant_plot
from visualizer.dual_volcano import build_dual_volcano_plot
from core.merger import get_quadrant_kpis

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
    Exports a self-contained interactive HTML report containing:
    1. Summary KPI Stat Cards
    2. Interactive 4-Quadrant Cross-Plot
    3. Dual Volcano Parallel Plot
    4. Top Q2 Splicing-Driven Target Genes Table
    5. Gemini AI Mechanism Evaluation (if available)
    """
    os.makedirs(os.path.dirname(output_html_path), exist_ok=True)

    # 1. Build Figures
    fig_quad = build_quadrant_plot(df_merged, log2fc_cutoff, delta_psi_cutoff, color_by="quadrant")
    fig_volc = build_dual_volcano_plot(df_merged, log2fc_cutoff, delta_psi_cutoff, deg_fdr_cutoff, as_fdr_cutoff)

    # Extract Plotly HTML snippets (include plotly js inline for offline viewing)
    quad_html = fig_quad.to_html(full_html=False, include_plotlyjs="cdn")
    volc_html = fig_volc.to_html(full_html=False, include_plotlyjs=False)

    # 2. Extract KPI Stats
    kpis = get_quadrant_kpis(df_merged)

    # 3. Format Top Q2 Genes Table (up to 20 top Q2 target genes)
    q2_genes_html = ""
    if df_merged.height > 0:
        q2_df = df_merged.filter(pl.col("quadrant") == "Q2").sort(pl.col("delta_psi").abs(), descending=True).head(20)
        if q2_df.height > 0:
            q2_pdf = q2_df.to_pandas()
            rows_html = ""
            for _, r in q2_pdf.iterrows():
                rows_html += f"""
                <tr>
                    <td><strong>{r['geneSymbol']}</strong></td>
                    <td><code>{r['gene_id']}</code></td>
                    <td><span class="badge badge-purple">{r['event_type']}</span></td>
                    <td><strong style="color: #8E44AD;">{r['delta_psi']:.3f}</strong></td>
                    <td>{r['log2FoldChange']:.3f}</td>
                    <td>{r['as_fdr']:.2e}</td>
                    <td><small>{r['coordinates']}</small></td>
                </tr>
                """
            q2_genes_html = f"""
            <table class="data-table">
                <thead>
                    <tr>
                        <th>Gene Symbol</th>
                        <th>Gene ID</th>
                        <th>Event Type</th>
                        <th>ΔPSI (Splicing)</th>
                        <th>Log₂FC (Expression)</th>
                        <th>rMATS FDR</th>
                        <th>Exon Coordinates</th>
                    </tr>
                </thead>
                <tbody>
                    {rows_html}
                </tbody>
            </table>
            """
        else:
            q2_genes_html = "<p class='no-data'>No Q2 target genes detected with current cutoffs.</p>"

    # 4. Format AI Section
    ai_section_html = ""
    if ai_insights:
        import markdown
        ai_html_content = markdown.markdown(ai_insights)
        ai_section_html = f"""
        <div class="card ai-card">
            <h2>🧬 Google Gemini AI Biological Mechanism Evaluation</h2>
            <div class="ai-content">
                {ai_html_content}
            </div>
        </div>
        """

    # 5. Full HTML Template
    full_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GenSplice-Agent Transcriptomics Report</title>
    <style>
        :root {{
            --primary: #8E44AD;
            --q1-color: #E63946;
            --q2-color: #8E44AD;
            --q3-color: #95A5A6;
            --q4-color: #2980B9;
            --bg-body: #F8FAFC;
            --bg-card: #FFFFFF;
            --text-main: #1E293B;
            --text-muted: #64748B;
            --border-color: #E2E8F0;
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
            background: linear-gradient(135deg, #4A0E4E 0%, #8E44AD 100%);
            color: white;
            padding: 32px;
            border-radius: 16px;
            margin-bottom: 24px;
            box-shadow: 0 10px 25px -5px rgba(142, 68, 173, 0.3);
        }}
        .header h1 {{
            margin: 0 0 8px 0;
            font-size: 28px;
            font-weight: 700;
        }}
        .header p {{
            margin: 0;
            opacity: 0.9;
            font-size: 15px;
        }}
        .kpi-container {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 16px;
            margin-bottom: 24px;
        }}
        .kpi-card {{
            background: var(--bg-card);
            border-radius: 12px;
            padding: 20px;
            border: 1px solid var(--border-color);
            box-shadow: 0 1px 3px rgba(0,0,0,0.05);
            text-align: center;
        }}
        .kpi-card .value {{
            font-size: 32px;
            font-weight: 800;
            margin-top: 4px;
        }}
        .kpi-card .label {{
            font-size: 13px;
            color: var(--text-muted);
            text-transform: uppercase;
            letter-spacing: 0.5px;
            font-weight: 600;
        }}
        .card {{
            background: var(--bg-card);
            border-radius: 16px;
            padding: 24px;
            margin-bottom: 24px;
            border: 1px solid var(--border-color);
            box-shadow: 0 4px 6px -1px rgba(0,0,0,0.05);
        }}
        .card h2 {{
            margin-top: 0;
            font-size: 20px;
            color: var(--text-main);
            border-bottom: 2px solid var(--border-color);
            padding-bottom: 12px;
        }}
        .data-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 16px;
            font-size: 14px;
        }}
        .data-table th, .data-table td {{
            padding: 12px 16px;
            text-align: left;
            border-bottom: 1px solid var(--border-color);
        }}
        .data-table th {{
            background-color: #F1F5F9;
            font-weight: 600;
            color: #475569;
        }}
        .badge {{
            display: inline-block;
            padding: 3px 8px;
            border-radius: 6px;
            font-size: 12px;
            font-weight: 600;
        }}
        .badge-purple {{
            background-color: rgba(142, 68, 173, 0.15);
            color: #8E44AD;
        }}
        .ai-card {{
            border-left: 6px solid #8E44AD;
            background: linear-gradient(180deg, #FFFFFF 0%, #FAF5FF 100%);
        }}
        .ai-content {{
            line-height: 1.7;
            font-size: 15px;
        }}
        .footer {{
            text-align: center;
            color: var(--text-muted);
            font-size: 13px;
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid var(--border-color);
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 GenSplice-Agent Transcriptomics Analysis Report</h1>
        <p>Integrative Quantity (DEG) vs Quality (Alternative Splicing 5 Event Types) Report</p>
    </div>

    <!-- Summary KPI Cards -->
    <div class="kpi-container">
        <div class="kpi-card">
            <div class="label">Total Analyzed Genes</div>
            <div class="value">{kpis['total']}</div>
        </div>
        <div class="kpi-card" style="border-top: 4px solid #8E44AD;">
            <div class="label" style="color: #8E44AD;">Q2: Splicing-Driven Target</div>
            <div class="value" style="color: #8E44AD;">{kpis['Q2']}</div>
        </div>
        <div class="kpi-card" style="border-top: 4px solid #E63946;">
            <div class="label" style="color: #E63946;">Q1: Dual Responders</div>
            <div class="value" style="color: #E63946;">{kpis['Q1']}</div>
        </div>
        <div class="kpi-card" style="border-top: 4px solid #2980B9;">
            <div class="label" style="color: #2980B9;">Q4: Abundance-Driven</div>
            <div class="value" style="color: #2980B9;">{kpis['Q4']}</div>
        </div>
        <div class="kpi-card" style="border-top: 4px solid #95A5A6;">
            <div class="label" style="color: #95A5A6;">Q3: Invariant Background</div>
            <div class="value" style="color: #95A5A6;">{kpis['Q3']}</div>
        </div>
    </div>

    <!-- 4-Quadrant Plot -->
    <div class="card">
        <h2>📊 4-Quadrant Transcriptomics Cross-Plot</h2>
        {quad_html}
    </div>

    <!-- Dual Volcano Plot -->
    <div class="card">
        <h2>🌋 Dual Volcano Parallel View</h2>
        {volc_html}
    </div>

    <!-- Q2 Primary Target Genes Table -->
    <div class="card">
        <h2>🎯 Primary Target Genes (Q2: Splicing-Driven / Masked)</h2>
        <p style="color: var(--text-muted); font-size: 14px;">
            Genes with statistically significant splicing variation (|ΔPSI| ≥ {delta_psi_cutoff}) but unchanged overall gene expression (|Log₂FC| &lt; {log2fc_cutoff}).
        </p>
        {q2_genes_html}
    </div>

    <!-- Gemini AI Section -->
    {ai_section_html}

    <div class="footer">
        Generated automatically by GenSplice-Agent Pipeline & Dashboard • Openable directly in Google Chrome / Browsers
    </div>
</body>
</html>
"""

    with open(output_html_path, "w", encoding="utf-8") as f:
        f.write(full_html)

    print(f"  ✔ Standalone HTML Report exported to: {output_html_path}")
    return output_html_path
