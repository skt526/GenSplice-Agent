"""
GenSplice-Agent Quadrant Cross-Plot Visualizer (Plotly)
Light Mode Theme with Shaded Translucent Threshold Regions & Right Legend
"""

import polars as pl
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from config import (
    QUADRANT_COLORS, QUADRANT_LABELS, EVENT_COLORS,
    DEG_THRESHOLD_COLOR, AS_THRESHOLD_COLOR, REGION_FILL_COLORS
)

def build_quadrant_plot(
    df_merged: pl.DataFrame,
    log2fc_cutoff: float = 1.0,
    delta_psi_cutoff: float = 0.1,
    color_by: str = "quadrant"  # 'quadrant' or 'event_type'
) -> go.Figure:
    """
    Creates a Light Mode 2D scatter plot with translucent shaded quadrant regions:
    - Q1 (Both DEG & Splicing): Purple Translucent Region (rgba(168, 85, 247, 0.12))
    - Q2 (Splicing Only): Blue Translucent Region (rgba(59, 130, 246, 0.12))
    - Q4 (DEG Only): Red Translucent Region (rgba(239, 68, 68, 0.12))
    - Q3 (Invariant): Soft Gray Translucent Region
    - Threshold Lines: Red vertical (DEG), Blue horizontal (Splicing)
    - Legend: Right side
    """
    if df_merged.height == 0:
        fig = go.Figure()
        fig.update_layout(
            title="No data available for Quadrant Plot",
            template="plotly_white",
            paper_bgcolor="#FFFFFF",
            plot_bgcolor="#FFFFFF"
        )
        return fig

    pdf = df_merged.to_pandas()

    # Hover text
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

    fig = go.Figure()

    # Determine plot bounds
    max_x = max(abs(pdf["log2FoldChange"].max() if len(pdf) > 0 else 2), 2.5) + 0.5
    max_y = max(abs(pdf["delta_psi"].max() if len(pdf) > 0 else 0.5), 0.6) + 0.1

    # -------------------------------------------------------------
    # 1. Translucent Quadrant Shaded Regions (Background Layer)
    # -------------------------------------------------------------
    shapes = [
        # Q3 (Center Invariant: Soft Gray)
        dict(
            type="rect", x0=-log2fc_cutoff, x1=log2fc_cutoff, y0=-delta_psi_cutoff, y1=delta_psi_cutoff,
            fillcolor=REGION_FILL_COLORS["Q3"], line=dict(width=0), layer="below"
        ),
        # Q2 (Top & Bottom Center Splicing Only: Translucent Blue)
        dict(
            type="rect", x0=-log2fc_cutoff, x1=log2fc_cutoff, y0=delta_psi_cutoff, y1=max_y,
            fillcolor=REGION_FILL_COLORS["Q2"], line=dict(width=0), layer="below"
        ),
        dict(
            type="rect", x0=-log2fc_cutoff, x1=log2fc_cutoff, y0=-max_y, y1=-delta_psi_cutoff,
            fillcolor=REGION_FILL_COLORS["Q2"], line=dict(width=0), layer="below"
        ),
        # Q4 (Left & Right Center DEG Only: Translucent Red)
        dict(
            type="rect", x0=log2fc_cutoff, x1=max_x, y0=-delta_psi_cutoff, y1=delta_psi_cutoff,
            fillcolor=REGION_FILL_COLORS["Q4"], line=dict(width=0), layer="below"
        ),
        dict(
            type="rect", x0=-max_x, x1=-log2fc_cutoff, y0=-delta_psi_cutoff, y1=delta_psi_cutoff,
            fillcolor=REGION_FILL_COLORS["Q4"], line=dict(width=0), layer="below"
        ),
        # Q1 (4 Corner Dual Responders: Translucent Purple)
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

        # -------------------------------------------------------------
        # 2. Dashed Threshold Lines
        # -------------------------------------------------------------
        # Vertical DEG Threshold Lines (RED)
        dict(
            type="line", x0=log2fc_cutoff, x1=log2fc_cutoff, y0=-max_y, y1=max_y,
            line=dict(color=DEG_THRESHOLD_COLOR, width=2, dash="dash")
        ),
        dict(
            type="line", x0=-log2fc_cutoff, x1=-log2fc_cutoff, y0=-max_y, y1=max_y,
            line=dict(color=DEG_THRESHOLD_COLOR, width=2, dash="dash")
        ),
        # Horizontal Splicing Threshold Lines (BLUE)
        dict(
            type="line", x0=-max_x, x1=max_x, y0=delta_psi_cutoff, y1=delta_psi_cutoff,
            line=dict(color=AS_THRESHOLD_COLOR, width=2, dash="dash")
        ),
        dict(
            type="line", x0=-max_x, x1=max_x, y0=-delta_psi_cutoff, y1=-delta_psi_cutoff,
            line=dict(color=AS_THRESHOLD_COLOR, width=2, dash="dash")
        )
    ]

    # Add Scatter Traces
    if color_by == "quadrant":
        for quad in ["Q1", "Q2", "Q4", "Q3"]:
            sub = pdf[pdf["quadrant"] == quad]
            if len(sub) == 0:
                continue
            fig.add_trace(
                go.Scatter(
                    x=sub["log2FoldChange"],
                    y=sub["delta_psi"],
                    mode="markers",
                    name=QUADRANT_LABELS.get(quad, quad),
                    marker=dict(
                        color=QUADRANT_COLORS.get(quad, "#94A3B8"),
                        size=10 if quad in ["Q1", "Q2", "Q4"] else 6,
                        opacity=0.9 if quad in ["Q1", "Q2", "Q4"] else 0.5,
                        line=dict(width=0.8, color="#FFFFFF")
                    ),
                    text=sub["hover_text"],
                    hoverinfo="text"
                )
            )
    else:
        event_types = pdf["event_type"].unique()
        for et in event_types:
            sub = pdf[pdf["event_type"] == et]
            fig.add_trace(
                go.Scatter(
                    x=sub["log2FoldChange"],
                    y=sub["delta_psi"],
                    mode="markers",
                    name=et,
                    marker=dict(
                        color=EVENT_COLORS.get(et, "#94A3B8"),
                        size=8,
                        opacity=0.85,
                        line=dict(width=0.5, color="#FFFFFF")
                    ),
                    text=sub["hover_text"],
                    hoverinfo="text"
                )
            )

    # Quadrant Badges
    annotations = [
        # Q2 (Top-Left): Splicing-Driven Target (BLUE)
        dict(
            x=-max_x * 0.72, y=max_y * 0.85,
            text="<b>Q2: Splicing-Driven</b><br>(Alternative Splicing Only)",
            showarrow=False,
            font=dict(size=11, color="#2563EB"),
            align="center",
            bordercolor="#3B82F6", borderwidth=1, borderpad=4, bgcolor="rgba(255, 255, 255, 0.9)"
        ),
        # Q1 (Top-Right): Dual Responders (PURPLE)
        dict(
            x=max_x * 0.72, y=max_y * 0.85,
            text="<b>Q1: Dual Responders</b><br>(Both DEG & Splicing)",
            showarrow=False,
            font=dict(size=11, color="#7C3AED"),
            align="center",
            bordercolor="#A855F7", borderwidth=1, borderpad=4, bgcolor="rgba(255, 255, 255, 0.9)"
        ),
        # Q3 (Bottom-Left): Invariant (GRAY)
        dict(
            x=-max_x * 0.72, y=-max_y * 0.85,
            text="<b>Q3: Background Invariant</b>",
            showarrow=False,
            font=dict(size=10, color="#64748B"),
            align="center",
            bordercolor="#CBD5E1", borderwidth=1, borderpad=4, bgcolor="rgba(255, 255, 255, 0.9)"
        ),
        # Q4 (Bottom-Right): Expression-Driven (RED)
        dict(
            x=max_x * 0.72, y=-max_y * 0.85,
            text="<b>Q4: Expression-Driven</b><br>(DEG Only)",
            showarrow=False,
            font=dict(size=10, color="#DC2626"),
            align="center",
            bordercolor="#EF4444", borderwidth=1, borderpad=4, bgcolor="rgba(255, 255, 255, 0.9)"
        ),
    ]

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
        annotations=annotations,
        # Legend on the RIGHT side
        legend=dict(
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
        margin=dict(l=60, r=220, t=80, b=60), # Right margin for legend & sliders
        height=680
    )

    return fig
