"""
GenSplice-Agent Dual Volcano Parallel Visualizer (Plotly)
"""

import numpy as np
import polars as pl
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from config import QUADRANT_COLORS

def build_dual_volcano_plot(
    df_merged: pl.DataFrame,
    log2fc_cutoff: float = 1.0,
    delta_psi_cutoff: float = 0.1,
    deg_fdr_cutoff: float = 0.05,
    as_fdr_cutoff: float = 0.05
) -> go.Figure:
    """
    Builds side-by-side parallel volcano plots:
    - Subplot 1 (Left): DEG Volcano Plot (Log2FC vs -log10(DEG_FDR))
    - Subplot 2 (Right): Splicing Volcano Plot (Delta_PSI vs -log10(AS_FDR))
    """
    if df_merged.height == 0:
        fig = go.Figure()
        fig.update_layout(title="No data available for Dual Volcano Plot")
        return fig

    pdf = df_merged.to_pandas()

    # Calculate -log10(FDR) with floating point protection
    pdf["neg_log10_deg_fdr"] = -np.log10(np.clip(pdf["deg_fdr"], 1e-15, 1.0))
    pdf["neg_log10_as_fdr"] = -np.log10(np.clip(pdf["as_fdr"], 1e-15, 1.0))

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=(
            "<b>Differential Expression (DEG) Volcano</b>",
            "<b>Alternative Splicing (rMATS) Volcano</b>"
        ),
        horizontal_spacing=0.12
    )

    # --------------------------
    # Subplot 1: DEG Volcano
    # --------------------------
    deg_sig = pdf["is_deg_sig"]
    
    # Non-sig DEG points
    fig.add_trace(
        go.Scatter(
            x=pdf[~deg_sig]["log2FoldChange"],
            y=pdf[~deg_sig]["neg_log10_deg_fdr"],
            mode="markers",
            name="Not Significant",
            marker=dict(color="#BDC3C7", size=6, opacity=0.5),
            text=pdf[~deg_sig]["geneSymbol"],
            hovertemplate="<b>%{text}</b><br>Log2FC: %{x:.3f}<br>-log10(FDR): %{y:.2f}<extra></extra>"
        ),
        row=1, col=1
    )
    # Significant DEG points
    fig.add_trace(
        go.Scatter(
            x=pdf[deg_sig]["log2FoldChange"],
            y=pdf[deg_sig]["neg_log10_deg_fdr"],
            mode="markers",
            name="Significant DEG",
            marker=dict(color="#E63946", size=8, opacity=0.85, line=dict(width=0.5, color="white")),
            text=pdf[deg_sig]["geneSymbol"],
            hovertemplate="<b>%{text}</b><br>Log2FC: %{x:.3f}<br>-log10(FDR): %{y:.2f}<extra></extra>"
        ),
        row=1, col=1
    )

    # --------------------------
    # Subplot 2: AS Volcano
    # --------------------------
    as_sig = pdf["is_as_sig"]
    
    # Non-sig AS points
    fig.add_trace(
        go.Scatter(
            x=pdf[~as_sig]["delta_psi"],
            y=pdf[~as_sig]["neg_log10_as_fdr"],
            mode="markers",
            name="Not Significant",
            showlegend=False,
            marker=dict(color="#BDC3C7", size=6, opacity=0.5),
            text=pdf[~as_sig]["geneSymbol"],
            hovertemplate="<b>%{text}</b><br>ΔPSI: %{x:.3f}<br>-log10(FDR): %{y:.2f}<extra></extra>"
        ),
        row=1, col=2
    )
    # Significant AS points
    fig.add_trace(
        go.Scatter(
            x=pdf[as_sig]["delta_psi"],
            y=pdf[as_sig]["neg_log10_as_fdr"],
            mode="markers",
            name="Significant Splicing (AS)",
            marker=dict(color="#8E44AD", size=8, opacity=0.85, line=dict(width=0.5, color="white")),
            text=pdf[as_sig]["geneSymbol"],
            hovertemplate="<b>%{text}</b><br>ΔPSI: %{x:.3f}<br>-log10(FDR): %{y:.2f}<extra></extra>"
        ),
        row=1, col=2
    )

    # Add Cutoff lines for Subplot 1
    deg_y_cutoff = -np.log10(deg_fdr_cutoff)
    fig.add_vline(x=log2fc_cutoff, line_dash="dash", line_color="#E74C3C", row=1, col=1)
    fig.add_vline(x=-log2fc_cutoff, line_dash="dash", line_color="#E74C3C", row=1, col=1)
    fig.add_hline(y=deg_y_cutoff, line_dash="dash", line_color="#E74C3C", row=1, col=1)

    # Add Cutoff lines for Subplot 2
    as_y_cutoff = -np.log10(as_fdr_cutoff)
    fig.add_vline(x=delta_psi_cutoff, line_dash="dash", line_color="#9B59B6", row=1, col=2)
    fig.add_vline(x=-delta_psi_cutoff, line_dash="dash", line_color="#9B59B6", row=1, col=2)
    fig.add_hline(y=as_y_cutoff, line_dash="dash", line_color="#9B59B6", row=1, col=2)

    fig.update_xaxes(title_text="<b>Log₂ Fold Change</b>", row=1, col=1, zerolinecolor="#CBD5E1")
    fig.update_yaxes(title_text="<b>-log₁₀(DEG FDR)</b>", row=1, col=1, zerolinecolor="#CBD5E1")

    fig.update_xaxes(title_text="<b>ΔPSI (Percent Spliced In Difference)</b>", row=1, col=2, zerolinecolor="#CBD5E1")
    fig.update_yaxes(title_text="<b>-log₁₀(rMATS FDR)</b>", row=1, col=2, zerolinecolor="#CBD5E1")

    fig.update_layout(
        template="plotly_white",
        height=550,
        legend=dict(orientation="h", yanchor="bottom", y=1.05, xanchor="center", x=0.5),
        margin=dict(l=60, r=60, t=80, b=60)
    )

    return fig
