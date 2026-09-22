"""
GenSplice-Agent Quadrant Cross-Plot Visualizer (Plotly)
Black & White Dark Theme with Red/Blue/Purple Marker & Threshold System
"""

import polars as pl
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from config import (
    QUADRANT_COLORS, QUADRANT_LABELS, EVENT_COLORS,
    DEG_THRESHOLD_COLOR, AS_THRESHOLD_COLOR
)

def build_quadrant_plot(
    df_merged: pl.DataFrame,
    log2fc_cutoff: float = 1.0,
    delta_psi_cutoff: float = 0.1,
    color_by: str = "quadrant"  # 'quadrant' or 'event_type'
) -> go.Figure:
    """
    Creates a Black & White dark mode 2D scatter plot of Log2FC (X-axis) vs Delta_PSI (Y-axis)
    - DEG Only: Red (#EF4444)
    - Splicing Only: Blue (#3B82F6)
    - Both DEG & Splicing: Purple (#A855F7)
    - Threshold lines: Red for Log2FC, Blue for Delta_PSI
    - Legend & Controls: Right side
    """
    if df_merged.height == 0:
        fig = go.Figure()
        fig.update_layout(
            title="No data available for Quadrant Plot",
            template="plotly_dark",
            paper_bgcolor="#0B0F19",
            plot_bgcolor="#111827"
        )
        return fig

    pdf = df_merged.to_pandas()

    # Create rich dark-mode hover text
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

    if color_by == "quadrant":
        for quad in ["Q1", "Q2", "Q4", "Q3"]: # Put targets on top
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
                        color=QUADRANT_COLORS.get(quad, "#4B5563"),
                        size=10 if quad in ["Q1", "Q2", "Q4"] else 6,
                        opacity=0.9 if quad in ["Q1", "Q2", "Q4"] else 0.4,
                        line=dict(width=0.8, color="#FFFFFF")
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
                        color=EVENT_COLORS.get(et, "#4B5563"),
                        size=8,
                        opacity=0.85,
                        line=dict(width=0.5, color="#FFFFFF")
                    ),
                    text=sub["hover_text"],
                    hoverinfo="text"
                )
            )

    # Determine plot bounds
    max_x = max(abs(pdf["log2FoldChange"].max() if len(pdf) > 0 else 2), 2.5) + 0.5
    max_y = max(abs(pdf["delta_psi"].max() if len(pdf) > 0 else 0.5), 0.6) + 0.1

    # Threshold Lines matching axes (Red for Log2FC DEG, Blue for Delta PSI Splicing)
    shapes = [
        # Vertical Log2FC DEG Thresholds (RED)
        dict(
            type="line", x0=log2fc_cutoff, x1=log2fc_cutoff, y0=-max_y, y1=max_y,
            line=dict(color=DEG_THRESHOLD_COLOR, width=2, dash="dash"),
            name="log2fc_pos"
        ),
        dict(
            type="line", x0=-log2fc_cutoff, x1=-log2fc_cutoff, y0=-max_y, y1=max_y,
            line=dict(color=DEG_THRESHOLD_COLOR, width=2, dash="dash"),
            name="log2fc_neg"
        ),
        # Horizontal Delta PSI Splicing Thresholds (BLUE)
        dict(
            type="line", x0=-max_x, x1=max_x, y0=delta_psi_cutoff, y1=delta_psi_cutoff,
            line=dict(color=AS_THRESHOLD_COLOR, width=2, dash="dash"),
            name="psi_pos"
        ),
        dict(
            type="line", x0=-max_x, x1=max_x, y0=-delta_psi_cutoff, y1=-delta_psi_cutoff,
            line=dict(color=AS_THRESHOLD_COLOR, width=2, dash="dash"),
            name="psi_neg"
        )
    ]

    # Add Quadrant Badges
    annotations = [
        # Q2 (Top-Left): Splicing-Driven Target (BLUE)
        dict(
            x=-max_x * 0.72, y=max_y * 0.85,
            text="<b>Q2: Splicing-Driven</b><br>(Alternative Splicing Only)",
            showarrow=False,
            font=dict(size=11, color="#3B82F6"),
            align="center",
            bordercolor="#3B82F6", borderwidth=1, borderpad=4, bgcolor="rgba(59, 130, 246, 0.15)"
        ),
        # Q1 (Top-Right): Dual Responders (PURPLE)
        dict(
            x=max_x * 0.72, y=max_y * 0.85,
            text="<b>Q1: Dual Responders</b><br>(Both DEG & Splicing)",
            showarrow=False,
            font=dict(size=11, color="#A855F7"),
            align="center",
            bordercolor="#A855F7", borderwidth=1, borderpad=4, bgcolor="rgba(168, 85, 247, 0.15)"
        ),
        # Q3 (Bottom-Left): Invariant (GRAY)
        dict(
            x=-max_x * 0.72, y=-max_y * 0.85,
            text="<b>Q3: Background Invariant</b>",
            showarrow=False,
            font=dict(size=10, color="#9CA3AF"),
            align="center",
            bordercolor="#4B5563", borderwidth=1, borderpad=4, bgcolor="rgba(75, 85, 99, 0.2)"
        ),
        # Q4 (Bottom-Right): Expression-Driven (RED)
        dict(
            x=max_x * 0.72, y=-max_y * 0.85,
            text="<b>Q4: Expression-Driven</b><br>(DEG Only)",
            showarrow=False,
            font=dict(size=10, color="#EF4444"),
            align="center",
            bordercolor="#EF4444", borderwidth=1, borderpad=4, bgcolor="rgba(239, 68, 68, 0.15)"
        ),
    ]

    fig.update_layout(
        title=dict(
            text="<b>GenSplice 4-Quadrant Transcriptomics Cross-Plot</b>",
            x=0.02,
            font=dict(size=20, family="sans-serif", color="#F9FAFB")
        ),
        xaxis_title=dict(text="<b>Log₂ Fold Change (DEG: Red Threshold)</b>", font=dict(size=14, color="#EF4444")),
        yaxis_title=dict(text="<b>ΔPSI (Alternative Splicing: Blue Threshold)</b>", font=dict(size=14, color="#3B82F6")),
        xaxis=dict(
            range=[-max_x, max_x],
            zeroline=True,
            zerolinecolor="#374151",
            gridcolor="#1F2937",
            tickfont=dict(color="#F9FAFB")
        ),
        yaxis=dict(
            range=[-max_y, max_y],
            zeroline=True,
            zerolinecolor="#374151",
            gridcolor="#1F2937",
            tickfont=dict(color="#F9FAFB")
        ),
        shapes=shapes,
        annotations=annotations,
        # Legend positioned on the RIGHT side
        legend=dict(
            orientation="v",
            yanchor="top",
            y=1.0,
            xanchor="left",
            x=1.02,
            bgcolor="#111827",
            bordercolor="#374151",
            borderwidth=1,
            font=dict(color="#F9FAFB", size=12)
        ),
        template="plotly_dark",
        paper_bgcolor="#0B0F19",
        plot_bgcolor="#111827",
        margin=dict(l=60, r=220, t=80, b=60), # Right margin expanded for legend & sliders
        height=680
    )

    return fig
