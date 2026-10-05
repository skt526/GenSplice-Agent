"""
GenSplice-Agent Alternative Splicing Volcano Plot Visualizer (Plotly)
Light Mode Theme with Shaded Translucent Significance Regions
- X-axis: ΔPSI (Inclusion Difference: Treatment − Control)
- Y-axis: -log10(rMATS FDR) (Statistical Significance)
- In-canvas legend: Splicing Direction (Inclusion Favored, Exclusion Favored, Non-Significant)
"""

import math
import numpy as np
import polars as pl
import pandas as pd
import plotly.graph_objects as go
from config import VOLCANO_COLORS, VOLCANO_SHADING

def build_as_volcano_plot(
    df_splicing: pl.DataFrame,
    delta_psi_cutoff: float = 0.10,
    fdr_cutoff: float = 0.05,
    max_fdr_cap: float = 1e-20
) -> go.Figure:
    """
    Creates an interactive 2D Volcano Plot for Alternative Splicing events:
    - X-axis: ΔPSI (-1.0 to +1.0)
    - Y-axis: -log10(rMATS FDR)
    - Translucent shaded significance zones for Inclusion Favored and Exclusion Favored
    - Dynamic cutoff lines and real-time hover information
    """
    if df_splicing is None or df_splicing.height == 0:
        pdf = pd.DataFrame(columns=[
            "geneSymbol", "gene_id", "delta_psi", "as_fdr", "as_pvalue",
            "event_type", "coordinates", "inc_counts", "exc_counts"
        ])
    else:
        # Filter for rows that have alternative splicing event data
        if "event_type" in df_splicing.columns:
            valid_df = df_splicing.filter(
                pl.col("event_type").is_not_null() & 
                (pl.col("event_type") != "None")
            )
            pdf = valid_df.to_pandas() if valid_df.height > 0 else df_splicing.to_pandas()
        else:
            pdf = df_splicing.to_pandas()

    fig = go.Figure()

    if len(pdf) == 0:
        fig.update_layout(
            title=dict(text="<b>No Alternative Splicing Events Available</b>", x=0.02),
            template="plotly_white",
            paper_bgcolor="#FFFFFF",
            plot_bgcolor="#FFFFFF",
            height=680
        )
        return fig

    # Ensure required columns exist
    if "delta_psi" not in pdf.columns:
        pdf["delta_psi"] = 0.0
    if "as_fdr" not in pdf.columns:
        pdf["as_fdr"] = 1.0
    if "as_pvalue" not in pdf.columns:
        pdf["as_pvalue"] = pdf["as_fdr"]
    if "geneSymbol" not in pdf.columns:
        pdf["geneSymbol"] = pdf["gene_id"] if "gene_id" in pdf.columns else "Unknown"
    if "gene_id" not in pdf.columns:
        pdf["gene_id"] = pdf["geneSymbol"]
    if "event_type" not in pdf.columns:
        pdf["event_type"] = "SE"
    if "coordinates" not in pdf.columns:
        pdf["coordinates"] = "N/A"
    if "inc_counts" not in pdf.columns:
        pdf["inc_counts"] = 0
    if "exc_counts" not in pdf.columns:
        pdf["exc_counts"] = 0

    # Fill NaN values safely
    pdf["delta_psi"] = pd.to_numeric(pdf["delta_psi"], errors="coerce").fillna(0.0)
    pdf["as_fdr"] = pd.to_numeric(pdf["as_fdr"], errors="coerce").fillna(1.0)
    pdf["as_pvalue"] = pd.to_numeric(pdf["as_pvalue"], errors="coerce").fillna(1.0)

    # Calculate -log10(FDR)
    def calc_log_p(fdr_val):
        v = float(fdr_val)
        if not np.isfinite(v) or v >= 1.0:
            return 0.0
        v = max(v, max_fdr_cap)
        return -math.log10(v)

    pdf["log_fdr"] = pdf["as_fdr"].apply(calc_log_p)

    # Classify Splicing Status
    status_list = []
    for _, r in pdf.iterrows():
        dpsi = float(r["delta_psi"])
        fdr = float(r["as_fdr"])
        is_sig = (abs(dpsi) >= delta_psi_cutoff) and (fdr <= fdr_cutoff)
        if is_sig and dpsi > 0:
            status_list.append("Inclusion")
        elif is_sig and dpsi < 0:
            status_list.append("Exclusion")
        else:
            status_list.append("Non-Significant")
    pdf["splicing_status"] = status_list

    # Hover text formatting
    hover_texts = []
    for _, r in pdf.iterrows():
        status_label = (
            "Significant Inclusion Favored (ΔPSI ≥ +" + f"{delta_psi_cutoff:.2f})"
            if r["splicing_status"] == "Inclusion"
            else (
                "Significant Exclusion Favored (ΔPSI ≤ -" + f"{delta_psi_cutoff:.2f})"
                if r["splicing_status"] == "Exclusion"
                else "Non-Significant Background"
            )
        )
        sign = "+" if r["delta_psi"] >= 0 else ""
        hover_texts.append(
            f"<b>{r['geneSymbol']}</b> ({r['gene_id']})<br>"
            f"Status: <b>{status_label}</b><br>"
            f"Event: {r['event_type']} | Coords: {r['coordinates']}<br>"
            f"ΔPSI: {sign}{r['delta_psi']:.3f}<br>"
            f"rMATS FDR: {r['as_fdr']:.2e} (-log₁₀ FDR: {r['log_fdr']:.2f})<br>"
            f"rMATS p-value: {r['as_pvalue']:.2e}<br>"
            f"Read Counts: Inc={int(r['inc_counts'])} | Exc={int(r['exc_counts'])}"
        )
    pdf["hover_text"] = hover_texts

    # Calculate plot axis ranges
    psi_series = pdf["delta_psi"].dropna()
    log_fdr_series = pdf["log_fdr"].dropna()

    max_data_x = float(psi_series.abs().max()) if len(psi_series) > 0 else 0.8
    max_data_y = float(log_fdr_series.max()) if len(log_fdr_series) > 0 else 3.5

    max_x = max(round(max_data_x + 0.1, 2), delta_psi_cutoff + 0.1, 1.05)
    fdr_threshold_y = -math.log10(max(fdr_cutoff, 1e-20))
    max_y = max(round(max_data_y + 1.0, 1), round(fdr_threshold_y + 1.5, 1), 5.0)

    # 1. Translucent Background Shaded Significance Regions
    shapes = [
        # Center & Bottom: Non-significant background
        dict(
            type="rect", x0=-max_x, x1=max_x, y0=0, y1=fdr_threshold_y,
            fillcolor=VOLCANO_SHADING["Background"], line=dict(width=0), layer="below"
        ),
        dict(
            type="rect", x0=-delta_psi_cutoff, x1=delta_psi_cutoff, y0=fdr_threshold_y, y1=max_y,
            fillcolor=VOLCANO_SHADING["Background"], line=dict(width=0), layer="below"
        ),
        # Upper Right: Significant Inclusion Favored (Translucent Blue)
        dict(
            type="rect", x0=delta_psi_cutoff, x1=max_x, y0=fdr_threshold_y, y1=max_y,
            fillcolor=VOLCANO_SHADING["Inclusion"], line=dict(width=0), layer="below"
        ),
        # Upper Left: Significant Exclusion Favored (Translucent Red)
        dict(
            type="rect", x0=-max_x, x1=-delta_psi_cutoff, y0=fdr_threshold_y, y1=max_y,
            fillcolor=VOLCANO_SHADING["Exclusion"], line=dict(width=0), layer="below"
        ),
        # 2. Dashed Threshold Lines
        # Horizontal line: FDR Cutoff
        dict(
            type="line", x0=-max_x, x1=max_x, y0=fdr_threshold_y, y1=fdr_threshold_y,
            line=dict(color="#6366F1", width=2, dash="dash")
        ),
        # Vertical lines: ΔPSI Cutoffs
        dict(
            type="line", x0=delta_psi_cutoff, x1=delta_psi_cutoff, y0=0, y1=max_y,
            line=dict(color="#2563EB", width=2, dash="dash")
        ),
        dict(
            type="line", x0=-delta_psi_cutoff, x1=-delta_psi_cutoff, y0=0, y1=max_y,
            line=dict(color="#E11D48", width=2, dash="dash")
        )
    ]

    # 3. Add Scatter Traces for Splicing Categories
    categories = [
        ("Inclusion", "Inclusion Favored", VOLCANO_COLORS["Inclusion"], 10, 0.9, True),
        ("Exclusion", "Exclusion Favored", VOLCANO_COLORS["Exclusion"], 10, 0.9, True),
        ("Non-Significant", "Non-Significant Background", VOLCANO_COLORS["Non-Significant"], 6, 0.5, True)
    ]

    for cat_id, cat_name, cat_color, pt_size, pt_opacity, pt_vis in categories:
        sub = pdf[pdf["splicing_status"] == cat_id] if "splicing_status" in pdf.columns else pd.DataFrame()
        # Downsample non-significant dots if >2500 to keep HTML lightweight
        if cat_id == "Non-Significant" and len(sub) > 2500:
            sub = sub.sample(n=2500, random_state=42)

        fig.add_trace(
            go.Scatter(
                x=sub["delta_psi"] if len(sub) > 0 else [],
                y=sub["log_fdr"] if len(sub) > 0 else [],
                mode="markers",
                name=cat_name,
                visible=pt_vis,
                marker=dict(
                    color=cat_color,
                    size=pt_size,
                    opacity=pt_opacity,
                    line=dict(width=0.8, color="#FFFFFF")
                ),
                text=sub["hover_text"] if len(sub) > 0 else [],
                hoverinfo="text"
            )
        )

    # 4. Regional In-plot Label Annotations
    annotations = [
        dict(
            x=-max_x * 0.65, y=max_y * 0.95,
            xref="x", yref="y",
            text="<b>Exclusion Favored</b> (ΔPSI < 0)",
            showarrow=False,
            font=dict(size=12, color="#BE123C", family="sans-serif"),
            bgcolor="rgba(255, 255, 255, 0.8)",
            bordercolor="#FECDD3",
            borderwidth=1,
            borderpad=4
        ),
        dict(
            x=max_x * 0.65, y=max_y * 0.95,
            xref="x", yref="y",
            text="<b>Inclusion Favored</b> (ΔPSI > 0)",
            showarrow=False,
            font=dict(size=12, color="#1D4ED8", family="sans-serif"),
            bgcolor="rgba(255, 255, 255, 0.8)",
            bordercolor="#BFDBFE",
            borderwidth=1,
            borderpad=4
        ),
        dict(
            x=max_x * 0.96, y=fdr_threshold_y + 0.15,
            xref="x", yref="y",
            text=f"FDR = {fdr_cutoff}",
            showarrow=False,
            font=dict(size=11, color="#4F46E5", family="sans-serif"),
            xanchor="right"
        )
    ]

    fig.update_layout(
        title=dict(
            text="<b>Alternative Splicing Volcano Plot (ΔPSI vs −log₁₀ FDR)</b>",
            x=0.02,
            font=dict(size=20, family="sans-serif", color="#0F172A")
        ),
        xaxis_title=dict(
            text="<b>ΔPSI (Exon Inclusion Difference: Treatment − Control)</b>",
            font=dict(size=14, color="#0F172A")
        ),
        yaxis_title=dict(
            text="<b>−log₁₀(rMATS FDR) (Statistical Significance)</b>",
            font=dict(size=14, color="#4F46E5")
        ),
        xaxis=dict(
            range=[-max_x, max_x],
            zeroline=True,
            zerolinecolor="#94A3B8",
            zerolinewidth=1.5,
            gridcolor="#F1F5F9",
            tickfont=dict(color="#0F172A")
        ),
        yaxis=dict(
            range=[0, max_y],
            zeroline=True,
            zerolinecolor="#CBD5E1",
            gridcolor="#F1F5F9",
            tickfont=dict(color="#0F172A")
        ),
        shapes=shapes,
        annotations=annotations,
        legend=dict(
            title=dict(text="<b>Splicing Status</b>", font=dict(size=13, color="#0F172A")),
            orientation="v",
            yanchor="top",
            y=1.0,
            xanchor="left",
            x=1.02,
            bgcolor="#FFFFFF",
            bordercolor="#E2E8F0",
            borderwidth=1,
            font=dict(color="#0F172A", size=13)
        ),
        template="plotly_white",
        paper_bgcolor="#FFFFFF",
        plot_bgcolor="#FFFFFF",
        margin=dict(l=60, r=220, t=80, b=60),
        height=680
    )

    return fig
