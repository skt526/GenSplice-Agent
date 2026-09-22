"""
GenSplice-Agent Standalone HTML Exporter Module
Light Mode interactive HTML report with:
- Dual Volcano Parallel View removed per user request.
- Legend Header explaining color meanings.
- One-line concise descriptions for Q1, Q2, Q3, Q4.
"""

import os
import polars as pl
import pandas as pd
from visualizer.quadrant_plot import build_quadrant_plot
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
    Exports a self-contained Light Mode interactive HTML report:
    - 4-Quadrant plot at top.
    - Dual Volcano view removed per user request.
    - Legend title explains color system with concise one-line Q1-Q4 descriptions.
    """
    os.makedirs(os.path.dirname(output_html_path), exist_ok=True)

    fig_quad = build_quadrant_plot(df_merged, log2fc_cutoff, delta_psi_cutoff, color_by="quadrant")
    quad_html = fig_quad.to_html(full_html=False, include_plotlyjs="cdn", div_id="plotly-quad-div")

    kpis = get_quadrant_kpis(df_merged)

    # Format Top Q2 Genes Table
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
                    <td><span class="badge badge-blue">{r['event_type']}</span></td>
                    <td><strong style="color: #3B82F6;">{r['delta_psi']:.3f}</strong></td>
                    <td>{r['log2FoldChange']:.3f}</td>
                    <td>{r['as_fdr']:.2e}</td>
                    <td><small>{r['coordinates']}</small></td>
                </tr>
                """
            q2_genes_html = f"""
            <table class="data-table" id="q2-table">
                <thead>
                    <tr>
                        <th>Gene Symbol</th>
                        <th>Gene ID</th>
                        <th>Event Type</th>
                        <th>ΔPSI (Splicing: Blue)</th>
                        <th>Log₂FC (DEG: Red)</th>
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
            q2_genes_html = "<p class='no-data' id='q2-table' style='color:#64748B;'>No Q2 target genes detected with current cutoffs.</p>"

    # Format AI Section
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

    raw_data_json = df_merged.to_pandas().to_json(orient="records")

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
        }}
        .card h2 {{ margin-top: 0; font-size: 20px; color: var(--text-main); border-bottom: 2px solid var(--border-color); padding-bottom: 12px; }}

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

        .legend-title-desc {{
            font-size: 13.5px;
            font-weight: 700;
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

        .data-table {{ width: 100%; border-collapse: collapse; margin-top: 16px; font-size: 14px; color: var(--text-main); }}
        .data-table th, .data-table td {{ padding: 12px 16px; text-align: left; border-bottom: 1px solid var(--border-color); }}
        .data-table th {{ background-color: #F1F5F9; font-weight: 600; color: #475569; }}
        .badge {{ display: inline-block; padding: 3px 8px; border-radius: 6px; font-size: 12px; font-weight: 600; }}
        .badge-blue {{ background-color: rgba(59, 130, 246, 0.15); color: #3B82F6; border: 1px solid #3B82F6; }}
        .ai-card {{ border-left: 6px solid var(--purple-both); background: #FFFFFF; }}
        .ai-content {{ line-height: 1.7; font-size: 15px; color: #334155; }}
        .footer {{ text-align: center; color: var(--text-muted); font-size: 13px; margin-top: 40px; padding-top: 20px; border-top: 1px solid var(--border-color); }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 GenSplice-Agent Interactive Light Report</h1>
        <p>Quadrant Classification & Color Meaning System</p>
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

                <!-- Color System Title & One-line descriptions -->
                <div>
                    <div class="legend-title-desc">Quadrant Classification & Color System:</div>
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

    <!-- Q2 Primary Target Genes Table -->
    <div class="card">
        <h2>🎯 Target Genes (Q2: Alternative Splicing Only - Blue)</h2>
        <div id="q2-table-container">
            {q2_genes_html}
        </div>
    </div>

    <!-- Gemini AI Section -->
    {ai_section_html}

    <div class="footer">
        Generated automatically by GenSplice-Agent Pipeline • Light Mode Theme
    </div>

    <!-- Real-Time Threshold & Dynamic CSV Script -->
    <script>
        const rawGeneData = {raw_data_json};
        
        const fcSlider = document.getElementById('fc-slider');
        const psiSlider = document.getElementById('psi-slider');
        const fcValLabel = document.getElementById('fc-val');
        const psiValLabel = document.getElementById('psi-val');
        const downloadBtn = document.getElementById('download-csv-btn');

        let currentFilteredGenes = [];

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
                const isDegSig = Math.abs(g.log2FoldChange) >= fcCut && (g.deg_fdr || 1.0) <= 0.05;
                const isAsSig = Math.abs(g.delta_psi) >= psiCut && (g.as_fdr || 1.0) <= 0.05;

                let quad = "Q3";
                if (isDegSig && isAsSig) {{ quad = "Q1"; q1++; }}
                else if (!isDegSig && isAsSig) {{ quad = "Q2"; q2++; }}
                else if (isDegSig && !isAsSig) {{ quad = "Q4"; q4++; }}
                else {{ quad = "Q3"; q3++; }}

                const hoverText = `<b>Gene Symbol:</b> ${{g.geneSymbol}}<br><b>Gene ID:</b> ${{g.gene_id}}<br><b>Quadrant:</b> ${{quad}}<br><b>Log2FC (DEG):</b> ${{g.log2FoldChange.toFixed(3)}}<br><b>ΔPSI (Splicing):</b> ${{g.delta_psi.toFixed(3)}}<br><b>DEG FDR:</b> ${{(g.deg_fdr || 1.0).toExponential(2)}}<br><b>rMATS FDR:</b> ${{(g.as_fdr || 1.0).toExponential(2)}}<br><b>Event Type:</b> ${{g.event_type || 'None'}}`;

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
                    {{ x: quadPoints.Q1.x, y: quadPoints.Q1.y, text: quadPoints.Q1.text, mode: 'markers', name: 'Q1: Both DEG & Splicing', marker: {{ color: '#A855F7', size: 10, opacity: 0.9, line: {{ width: 0.8, color: '#FFFFFF' }} }}, hoverinfo: 'text' }},
                    {{ x: quadPoints.Q2.x, y: quadPoints.Q2.y, text: quadPoints.Q2.text, mode: 'markers', name: 'Q2: Splicing Only', marker: {{ color: '#3B82F6', size: 10, opacity: 0.9, line: {{ width: 0.8, color: '#FFFFFF' }} }}, hoverinfo: 'text' }},
                    {{ x: quadPoints.Q3.x, y: quadPoints.Q3.y, text: quadPoints.Q3.text, mode: 'markers', name: 'Q3: Invariant Background', marker: {{ color: '#94A3B8', size: 6, opacity: 0.5, line: {{ width: 0.5, color: '#FFFFFF' }} }}, hoverinfo: 'text' }},
                    {{ x: quadPoints.Q4.x, y: quadPoints.Q4.y, text: quadPoints.Q4.text, mode: 'markers', name: 'Q4: DEG Only', marker: {{ color: '#EF4444', size: 10, opacity: 0.9, line: {{ width: 0.8, color: '#FFFFFF' }} }}, hoverinfo: 'text' }}
                ];

                Plotly.react(quadDiv, updatedTraces, quadDiv.layout);
                Plotly.relayout(quadDiv, {{ shapes: newShapes }});
            }}
        }}

        function downloadFilteredCSV() {{
            if (!currentFilteredGenes.length) return;

            const headers = ["geneSymbol", "gene_id", "current_quadrant", "log2FoldChange", "delta_psi", "deg_fdr", "as_fdr", "event_type", "coordinates"];
            let csvContent = headers.join(",") + "\\n";

            currentFilteredGenes.forEach(g => {{
                const row = [
                    `"${{g.geneSymbol || ''}}"`,
                    `"${{g.gene_id || ''}}"`,
                    `"${{g.current_quadrant || ''}}"`,
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

        fcSlider.addEventListener('input', updateThresholds);
        psiSlider.addEventListener('input', updateThresholds);
        downloadBtn.addEventListener('click', downloadFilteredCSV);

        updateThresholds();
    </script>
</body>
</html>
"""

    with open(output_html_path, "w", encoding="utf-8") as f:
        f.write(full_html)

    print(f"  ✔ Standalone Light Mode HTML Report updated with color title & one-line Q1-Q4 descriptions: {output_html_path}")
    return output_html_path
