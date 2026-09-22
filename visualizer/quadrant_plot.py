"""
GenSplice-Agent Quadrant Cross-Plot Visualizer (Plotly)
Light Mode Theme with Shaded Translucent Threshold Regions
- Canvas corner description boxes removed.
- Legend title explains color system.
- One-line descriptions for Q1, Q2, Q3, Q4 in legend.
"""

import polars as pl
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from config import (
    QUADRANT_COLORS, EVENT_COLORS,
    DEG_THRESHOLD_COLOR, AS_THRESHOLD_COLOR, REGION_FILL_COLORS
)

# Concise one-line labels for Q1 - Q4
QUADRANT_SHORT_DESCS = {
    "Q1": "Q1: Both DEG & Splicing",
    "Q2": "Q2: Splicing Only (Target)",
    "Q3": "Q3: Invariant Background",
    "Q4": "Q4: DEG Only"
}

def build_quadrant_plot(
    df_merged: pl.DataFrame,
    log2fc_cutoff: float = 1.0,
    delta_psi_cutoff: float = 0.1,
    color_by: str = "quadrant"  # 'quadrant' or 'event_type'
) -> go.Figure:
    """
    Creates a 2D scatter plot with translucent shaded regions and Q1 - Q4 color descriptions.
    """
    pdf = df_merged.to_pandas() if df_merged.height > 0 else pd.DataFrame(columns=[
        "geneSymbol", "gene_id", "quadrant", "log2FoldChange", "delta_psi",
        "deg_fdr", "as_fdr", "event_type", "coordinates"
    ])

    # Hover text generator
    if len(pdf) > 0:
        hover_texts = []
        for _, row in pdf.iterrows():
            ht = (
                f"<b>Gene Symbol:</b> {row['geneSymbol']}<br>"
                f"<b>Gene ID:</b> {row['gene_id']}<br>"
                f"<b>Quadrant:</b> {row['quadrant']}<br>"
                f"<b>Log2FC (DEG):</b> {row['log2FoldChange']:.3f}<br>"
                f"<b>ΔPSI (Splicing):</b> {row['delta_psi']:.3f}<br>"
                f"<b>DEG FDR:</b> {row['deg_fdr']:.2e}<br>"
                f"<b>rMATS FDR:</b> {row['as_fdr']:.2e}<br>"
                f"<b>Event Type:</b> {row['event_type']}<br>"
                f"<b>Coordinates:</b> {row['coordinates']}"
            )
            hover_texts.append(ht)
        pdf["hover_text"] = hover_texts
    else:
        pdf["hover_text"] = []

    fig = go.Figure()

    # Determine plot bounds
    max_x = max(abs(pdf["log2FoldChange"].max() if len(pdf) > 0 else 2), 2.5) + 0.5
    max_y = max(abs(pdf["delta_psi"].max() if len(pdf) > 0 else 0.5), 0.6) + 0.1

    # 1. Translucent Background Shaded Quadrant Regions
    shapes = [
        # Q3 Center Soft Gray
        dict(
            type="rect", x0=-log2fc_cutoff, x1=log2fc_cutoff, y0=-delta_psi_cutoff, y1=delta_psi_cutoff,
            fillcolor=REGION_FILL_COLORS["Q3"], line=dict(width=0), layer="below"
        ),
        # Q2 Splicing Only Translucent Blue
        dict(
            type="rect", x0=-log2fc_cutoff, x1=log2fc_cutoff, y0=delta_psi_cutoff, y1=max_y,
            fillcolor=REGION_FILL_COLORS["Q2"], line=dict(width=0), layer="below"
        ),
        dict(
            type="rect", x0=-log2fc_cutoff, x1=log2fc_cutoff, y0=-max_y, y1=-delta_psi_cutoff,
            fillcolor=REGION_FILL_COLORS["Q2"], line=dict(width=0), layer="below"
        ),
        # Q4 DEG Only Translucent Red
        dict(
            type="rect", x0=log2fc_cutoff, x1=max_x, y0=-delta_psi_cutoff, y1=delta_psi_cutoff,
            fillcolor=REGION_FILL_COLORS["Q4"], line=dict(width=0), layer="below"
        ),
        dict(
            type="rect", x0=-max_x, x1=-log2fc_cutoff, y0=-delta_psi_cutoff, y1=delta_psi_cutoff,
            fillcolor=REGION_FILL_COLORS["Q4"], line=dict(width=0), layer="below"
        ),
        # Q1 Dual Responders Translucent Purple
        dict(
            type="rect", x0=log2fc_cutoff, x1=max_x, y0=delta_psi_cutoff, y1=max_y,
            fillcolor=REGION_FILL_COLORS["Q1"], line=dict(width=0), layer="below"
        ),
        dict(
            type="rect", x0=-max_x, x1=-log2fc_cutoff, y0=delta_psi_cutoff, y1=max_y,
            fillcolor=REGION_FILL_COLORS["Q1"], line=dict(width=0), layer="below"
        ),
        dict(
            type="rect", x0=log2fc_cutoff, x1=max_x, y0=-max_y, y1=-delta_psi_cutoff,
            fillcolor=REGION_FILL_COLORS["Q1"], line=dict(width=0), layer="below"
        ),
        dict(
            type="rect", x0=-max_x, x1=-log2fc_cutoff, y0=-max_y, y1=-delta_psi_cutoff,
            fillcolor=REGION_FILL_COLORS["Q1"], line=dict(width=0), layer="below"
        ),

        # 2. Dashed Threshold Lines
        dict(
            type="line", x0=log2fc_cutoff, x1=log2fc_cutoff, y0=-max_y, y1=max_y,
            line=dict(color=DEG_THRESHOLD_COLOR, width=2, dash="dash")
        ),
        dict(
            type="line", x0=-log2fc_cutoff, x1=-log2fc_cutoff, y0=-max_y, y1=max_y,
            line=dict(color=DEG_THRESHOLD_COLOR, width=2, dash="dash")
        ),
        dict(
            type="line", x0=-max_x, x1=max_x, y0=delta_psi_cutoff, y1=delta_psi_cutoff,
            line=dict(color=AS_THRESHOLD_COLOR, width=2, dash="dash")
        ),
        dict(
            type="line", x0=-max_x, x1=max_x, y0=-delta_psi_cutoff, y1=-delta_psi_cutoff,
            line=dict(color=AS_THRESHOLD_COLOR, width=2, dash="dash")
        )
    ]

    # 3. Add Traces for all 4 quadrants (Q1, Q2, Q3, Q4)
    if color_by == "quadrant":
        all_quadrants = ["Q1", "Q2", "Q3", "Q4"]
        for quad in all_quadrants:
            sub = pdf[pdf["quadrant"] == quad] if len(pdf) > 0 and "quadrant" in pdf.columns else pd.DataFrame()
            fig.add_trace(
                go.Scatter(
                    x=sub["log2FoldChange"] if len(sub) > 0 else [],
                    y=sub["delta_psi"] if len(sub) > 0 else [],
                    mode="markers",
                    name=QUADRANT_SHORT_DESCS.get(quad, quad),
                    marker=dict(
                        color=QUADRANT_COLORS.get(quad, "#94A3B8"),
                        size=10 if quad in ["Q1", "Q2", "Q4"] else 6,
                        opacity=0.9 if quad in ["Q1", "Q2", "Q4"] else 0.5,
                        line=dict(width=0.8, color="#FFFFFF")
                    ),
                    text=sub["hover_text"] if len(sub) > 0 else [],
                    hoverinfo="text"
                )
            )
    else:
        event_types = ["SE", "RI", "MXE", "A5SS", "A3SS"]
        for et in event_types:
            sub = pdf[pdf["event_type"] == et] if len(pdf) > 0 and "event_type" in pdf.columns else pd.DataFrame()
            fig.add_trace(
                go.Scatter(
                    x=sub["log2FoldChange"] if len(sub) > 0 else [],
                    y=sub["delta_psi"] if len(sub) > 0 else [],
                    mode="markers",
                    name=et,
                    marker=dict(
                        color=EVENT_COLORS.get(et, "#94A3B8"),
                        size=8,
                        opacity=0.85,
                        line=dict(width=0.5, color="#FFFFFF")
                    ),
                    text=sub["hover_text"] if len(sub) > 0 else [],
                    hoverinfo="text"
                )
            )

    fig.update_layout(
        title=dict(
            text="<b>GenSplice 4-Quadrant Transcriptomics Cross-Plot</b>",
            x=0.02,
            font=dict(size=20, family="sans-serif", color="#0F172A")
        ),
        xaxis_title=dict(text="<b>Log₂ Fold Change (DEG: Red Threshold)</b>", font=dict(size=14, color="#DC2626")),
        yaxis_title=dict(text="<b>ΔPSI (Alternative Splicing: Blue Threshold)</b>", font=dict(size=14, color="#2563EB")),
        xaxis=dict(
            range=[-max_x, max_x],
            zeroline=True,
            zerolinecolor="#CBD5E1",
            gridcolor="#F1F5F9",
            tickfont=dict(color="#0F172A")
        ),
        yaxis=dict(
            range=[-max_y, max_y],
            zeroline=True,
            zerolinecolor="#CBD5E1",
            gridcolor="#F1F5F9",
            tickfont=dict(color="#0F172A")
        ),
        shapes=shapes,
        annotations=[],
        legend=dict(
            title=dict(text="<b>Quadrant Classification & Color System</b>", font=dict(size=13, color="#0F172A")),
            orientation="v",
            yanchor="top",
            y=1.0,
            xanchor="left",
            x=1.02,
            bgcolor="#FFFFFF",
            bordercolor="#E2E8F0",
            borderwidth=1,
            font=dict(color="#0F172A", size=12)
        ),
        template="plotly_white",
        paper_bgcolor="#FFFFFF",
        plot_bgcolor="#FFFFFF",
        margin=dict(l=60, r=260, t=80, b=60),
        height=680
    )

    return fig
