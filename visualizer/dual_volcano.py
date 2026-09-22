"""
GenSplice-Agent Dual Volcano Parallel Visualizer (Plotly)
Black & White Dark Mode Theme
"""

import numpy as np
import polars as pl
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from config import DEG_THRESHOLD_COLOR, AS_THRESHOLD_COLOR

def build_dual_volcano_plot(
    df_merged: pl.DataFrame,
    log2fc_cutoff: float = 1.0,
    delta_psi_cutoff: float = 0.1,
    deg_fdr_cutoff: float = 0.05,
    as_fdr_cutoff: float = 0.05
) -> go.Figure:
    """
    Builds side-by-side parallel volcano plots in Black & White dark theme:
    - Subplot 1 (Left): DEG Volcano Plot (Red points & threshold)
    - Subplot 2 (Right): Splicing Volcano Plot (Blue points & threshold)
    """
    if df_merged.height == 0:
        fig = go.Figure()
        fig.update_layout(
            title="No data available for Dual Volcano Plot",
            template="plotly_dark",
            paper_bgcolor="#0B0F19",
            plot_bgcolor="#111827"
        )
        return fig

    pdf = df_merged.to_pandas()

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
    # Subplot 1: DEG Volcano (Red)
    # --------------------------
    deg_sig = pdf["is_deg_sig"]
    
    # Non-sig DEG points
    fig.add_trace(
        go.Scatter(
            x=pdf[~deg_sig]["log2FoldChange"],
            y=pdf[~deg_sig]["neg_log10_deg_fdr"],
            mode="markers",
            name="Not Significant DEG",
            marker=dict(color="#4B5563", size=6, opacity=0.4),
            text=pdf[~deg_sig]["geneSymbol"],
            hovertemplate="<b>%{text}</b><br>Log2FC: %{x:.3f}<br>-log10(FDR): %{y:.2f}<extra></extra>"
        ),
        row=1, col=1
    )
    # Significant DEG points (RED)
    fig.add_trace(
        go.Scatter(
            x=pdf[deg_sig]["log2FoldChange"],
            y=pdf[deg_sig]["neg_log10_deg_fdr"],
            mode="markers",
            name="Significant DEG (Red)",
            marker=dict(color="#EF4444", size=8, opacity=0.9, line=dict(width=0.5, color="#FFFFFF")),
            text=pdf[deg_sig]["geneSymbol"],
            hovertemplate="<b>%{text}</b><br>Log2FC: %{x:.3f}<br>-log10(FDR): %{y:.2f}<extra></extra>"
        ),
        row=1, col=1
    )

    # --------------------------
    # Subplot 2: AS Volcano (Blue)
    # --------------------------
    as_sig = pdf["is_as_sig"]
    
    # Non-sig AS points
    fig.add_trace(
        go.Scatter(
            x=pdf[~as_sig]["delta_psi"],
            y=pdf[~as_sig]["neg_log10_as_fdr"],
            mode="markers",
            name="Not Significant Splicing",
            showlegend=False,
            marker=dict(color="#4B5563", size=6, opacity=0.4),
            text=pdf[~as_sig]["geneSymbol"],
            hovertemplate="<b>%{text}</b><br>ΔPSI: %{x:.3f}<br>-log10(FDR): %{y:.2f}<extra></extra>"
        ),
        row=1, col=2
    )
    # Significant AS points (BLUE)
    fig.add_trace(
        go.Scatter(
            x=pdf[as_sig]["delta_psi"],
            y=pdf[as_sig]["neg_log10_as_fdr"],
            mode="markers",
            name="Significant Splicing (Blue)",
            marker=dict(color="#3B82F6", size=8, opacity=0.9, line=dict(width=0.5, color="#FFFFFF")),
            text=pdf[as_sig]["geneSymbol"],
            hovertemplate="<b>%{text}</b><br>ΔPSI: %{x:.3f}<br>-log10(FDR): %{y:.2f}<extra></extra>"
        ),
        row=1, col=2
    )

    # Add Cutoff lines for Subplot 1 (RED)
    deg_y_cutoff = -np.log10(deg_fdr_cutoff)
    fig.add_vline(x=log2fc_cutoff, line_dash="dash", line_color=DEG_THRESHOLD_COLOR, row=1, col=1)
    fig.add_vline(x=-log2fc_cutoff, line_dash="dash", line_color=DEG_THRESHOLD_COLOR, row=1, col=1)
    fig.add_hline(y=deg_y_cutoff, line_dash="dash", line_color=DEG_THRESHOLD_COLOR, row=1, col=1)

    # Add Cutoff lines for Subplot 2 (BLUE)
    as_y_cutoff = -np.log10(as_fdr_cutoff)
    fig.add_vline(x=delta_psi_cutoff, line_dash="dash", line_color=AS_THRESHOLD_COLOR, row=1, col=2)
    fig.add_vline(x=-delta_psi_cutoff, line_dash="dash", line_color=AS_THRESHOLD_COLOR, row=1, col=2)
    fig.add_hline(y=as_y_cutoff, line_dash="dash", line_color=AS_THRESHOLD_COLOR, row=1, col=2)

    fig.update_xaxes(title_text="<b>Log₂ Fold Change</b>", row=1, col=1, zerolinecolor="#374151", gridcolor="#1F2937")
    fig.update_yaxes(title_text="<b>-log₁₀(DEG FDR)</b>", row=1, col=1, zerolinecolor="#374151", gridcolor="#1F2937")

    fig.update_xaxes(title_text="<b>ΔPSI (Difference)</b>", row=1, col=2, zerolinecolor="#374151", gridcolor="#1F2937")
    fig.update_yaxes(title_text="<b>-log₁₀(rMATS FDR)</b>", row=1, col=2, zerolinecolor="#374151", gridcolor="#1F2937")

    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor="#0B0F19",
        plot_bgcolor="#111827",
        height=550,
        legend=dict(
            orientation="v",
            yanchor="top",
            y=1.0,
            xanchor="left",
            x=1.02,
            bgcolor="#111827",
            bordercolor="#374151",
            borderwidth=1,
            font=dict(color="#F9FAFB")
        ),
        margin=dict(l=60, r=200, t=80, b=60)
    )

    return fig
