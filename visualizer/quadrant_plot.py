"""
GenSplice-Agent Quadrant Cross-Plot Visualizer (Plotly)
"""

import polars as pl
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from config import QUADRANT_COLORS, QUADRANT_LABELS, EVENT_COLORS

def build_quadrant_plot(
    df_merged: pl.DataFrame,
    log2fc_cutoff: float = 1.0,
    delta_psi_cutoff: float = 0.1,
    color_by: str = "quadrant"  # 'quadrant' or 'event_type'
) -> go.Figure:
    """
    Creates an interactive 2D scatter plot of Log2FC (X-axis) vs Delta_PSI (Y-axis)
    divided into 4 functional quadrants with dashed thresholds and interactive tooltips.
    """
    if df_merged.height == 0:
        fig = go.Figure()
        fig.update_layout(title="No data available for Quadrant Plot")
        return fig

    pdf = df_merged.to_pandas()

    # Create hover text
    hover_texts = []
    for _, row in pdf.iterrows():
        ht = (
            f"<b>Gene Symbol:</b> {row['geneSymbol']}<br>"
            f"<b>Gene ID:</b> {row['gene_id']}<br>"
            f"<b>Quadrant:</b> {row['quadrant']}<br>"
            f"<b>Log2FC (Expression):</b> {row['log2FoldChange']:.3f}<br>"
            f"<b>ΔPSI (Splicing):</b> {row['delta_psi']:.3f}<br>"
            f"<b>DEG FDR:</b> {row['deg_fdr']:.2e}<br>"
            f"<b>AS FDR:</b> {row['as_fdr']:.2e}<br>"
            f"<b>Event Type:</b> {row['event_type']}<br>"
            f"<b>Coordinates:</b> {row['coordinates']}"
        )
        hover_texts.append(ht)

    pdf["hover_text"] = hover_texts

    fig = go.Figure()

    if color_by == "quadrant":
        for quad in ["Q1", "Q2", "Q3", "Q4"]:
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
                        color=QUADRANT_COLORS.get(quad, "#999999"),
                        size=9 if quad in ["Q1", "Q2"] else 6,
                        opacity=0.85 if quad in ["Q1", "Q2"] else 0.5,
                        line=dict(width=0.5, color="white")
                    ),
                    text=sub["hover_text"],
                    hoverinfo="text"
                )
            )
    else: # color by event_type
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
                        color=EVENT_COLORS.get(et, "#7F8C8D"),
                        size=7,
                        opacity=0.8,
                        line=dict(width=0.5, color="white")
                    ),
                    text=sub["hover_text"],
                    hoverinfo="text"
                )
            )

    # Determine plot bounds
    max_x = max(abs(pdf["log2FoldChange"].max() if len(pdf) > 0 else 2), 2.5) + 0.5
    max_y = max(abs(pdf["delta_psi"].max() if len(pdf) > 0 else 0.5), 0.6) + 0.1

    # Add Dashed Threshold Lines (Draggable shapes)
    shapes = [
        # Right Log2FC Threshold Line
        dict(
            type="line", x0=log2fc_cutoff, x1=log2fc_cutoff, y0=-max_y, y1=max_y,
            line=dict(color="#E74C3C", width=2, dash="dash"),
            name="log2fc_pos"
        ),
        # Left Log2FC Threshold Line
        dict(
            type="line", x0=-log2fc_cutoff, x1=-log2fc_cutoff, y0=-max_y, y1=max_y,
            line=dict(color="#E74C3C", width=2, dash="dash"),
            name="log2fc_neg"
        ),
        # Upper Delta PSI Threshold Line
        dict(
            type="line", x0=-max_x, x1=max_x, y0=delta_psi_cutoff, y1=delta_psi_cutoff,
            line=dict(color="#9B59B6", width=2, dash="dash"),
            name="psi_pos"
        ),
        # Lower Delta PSI Threshold Line
        dict(
            type="line", x0=-max_x, x1=max_x, y0=-delta_psi_cutoff, y1=-delta_psi_cutoff,
            line=dict(color="#9B59B6", width=2, dash="dash"),
            name="psi_neg"
        )
    ]

    # Add Quadrant Corner Labels / Badges
    annotations = [
        # Q2 (Top-Left): Splicing-Driven Target
        dict(
            x=-max_x * 0.75, y=max_y * 0.85,
            text="<b>Q2: Splicing-Driven</b><br>(Primary Target: No DEG, Splicing Changed)",
            showarrow=False,
            font=dict(size=11, color="#8E44AD"),
            align="center",
            bordercolor="#8E44AD", borderwidth=1, borderpad=4, bgcolor="rgba(142, 68, 173, 0.1)"
        ),
        # Q1 (Top-Right): Dual Responders
        dict(
            x=max_x * 0.75, y=max_y * 0.85,
            text="<b>Q1: Dual Responders</b><br>(Both DEG & Splicing Changed)",
            showarrow=False,
            font=dict(size=11, color="#E63946"),
            align="center",
            bordercolor="#E63946", borderwidth=1, borderpad=4, bgcolor="rgba(230, 57, 70, 0.1)"
        ),
        # Q3 (Bottom-Left): Invariant
        dict(
            x=-max_x * 0.75, y=-max_y * 0.85,
            text="<b>Q3: Background Invariant</b>",
            showarrow=False,
            font=dict(size=10, color="#7F8C8D"),
            align="center",
            bordercolor="#7F8C8D", borderwidth=1, borderpad=4, bgcolor="rgba(127, 140, 141, 0.1)"
        ),
        # Q4 (Bottom-Right): Abundance-Driven
        dict(
            x=max_x * 0.75, y=-max_y * 0.85,
            text="<b>Q4: Expression-Driven</b><br>(DEG Only, No Splicing Change)",
            showarrow=False,
            font=dict(size=10, color="#2980B9"),
            align="center",
            bordercolor="#2980B9", borderwidth=1, borderpad=4, bgcolor="rgba(41, 128, 185, 0.1)"
        ),
    ]

    fig.update_layout(
        title=dict(
            text="<b>GenSplice 4-Quadrant Cross-Plot (Real-Time Interactive Thresholds)</b>",
            x=0.5,
            font=dict(size=18, family="sans-serif")
        ),
        xaxis_title=dict(text="<b>Log₂ Fold Change (Expression Quantity)</b>", font=dict(size=14)),
        yaxis_title=dict(text="<b>ΔPSI (Alternative Splicing Quality)</b>", font=dict(size=14)),
        xaxis=dict(range=[-max_x, max_x], zeroline=True, zerolinecolor="#CBD5E1", gridcolor="#F1F5F9"),
        yaxis=dict(range=[-max_y, max_y], zeroline=True, zerolinecolor="#CBD5E1", gridcolor="#F1F5F9"),
        shapes=shapes,
        annotations=annotations,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="center",
            x=0.5,
            bgcolor="rgba(255,255,255,0.8)",
            bordercolor="#E2E8F0",
            borderwidth=1
        ),
        template="plotly_white",
        margin=dict(l=60, r=60, t=100, b=60),
        height=680
    )

    return fig
