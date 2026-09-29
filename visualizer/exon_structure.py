"""
GenSplice-Agent Exon-Intron Isoform Structure Engine (Visualizer)
Renders interactive Sashimi-style exon-intron structure and junction splice arcs using Plotly.
"""

import plotly.graph_objects as go

def plot_exon_structure(
    gene_symbol: str = "STAT3",
    event_type: str = "SE",
    coordinates: str = "chr17:42300000:42301500:42303000",
    inc_counts: int = 145,
    exc_counts: int = 22,
    delta_psi: float = 0.35
) -> go.Figure:
    """
    Renders an interactive 2D Sashimi-style genomic structure plot for alternative splicing events.
    Displays Inclusion isoform (top) vs Exclusion isoform (bottom) with exon blocks and junction arcs.
    """
    fig = go.Figure()

    # Parse coordinates or fallback to mock coordinates if unparseable
    try:
        parts = coordinates.replace(":", "-").split("-")
        chrom = parts[0]
        ex1_start, ex1_end = int(parts[1]), int(parts[2])
        ex2_start, ex2_end = int(parts[3]), int(parts[4])
        ex3_start, ex3_end = int(parts[5]), int(parts[6])
    except Exception:
        chrom = "chr17"
        ex1_start, ex1_end = 1000, 1200
        ex2_start, ex2_end = 1600, 1750  # Skipped Exon
        ex3_start, ex3_end = 2200, 2400

    # Genomic coordinate span
    x_min = ex1_start - 200
    x_max = ex3_end + 200

    # Colors
    color_constitutive = "#2b5c8f"  # Blue
    color_skipped = "#e74c3c"       # Red
    color_inc_arc = "#27ae60"        # Green
    color_exc_arc = "#8e44ad"        # Purple

    # -------------------------------------------------------------
    # 1. Draw Intron Backbone Lines
    # -------------------------------------------------------------
    # Inclusion Track (Y = 2)
    fig.add_trace(go.Scatter(
        x=[ex1_start, ex3_end], y=[2, 2],
        mode="lines", line=dict(color="#7f8c8d", width=2, dash="dash"),
        name="Intron (Inclusion Track)", hoverinfo="skip"
    ))

    # Exclusion Track (Y = -2)
    fig.add_trace(go.Scatter(
        x=[ex1_start, ex3_end], y=[-2, -2],
        mode="lines", line=dict(color="#7f8c8d", width=2, dash="dash"),
        name="Intron (Exclusion Track)", hoverinfo="skip"
    ))

    # -------------------------------------------------------------
    # 2. Draw Exon Boxes
    # -------------------------------------------------------------
    def add_exon_box(x0, x1, y_center, height, color, name):
        fig.add_shape(
            type="rect",
            x0=x0, x1=x1,
            y0=y_center - height/2, y1=y_center + height/2,
            fillcolor=color, line=dict(color="black", width=1),
            layer="above"
        )
        # Invisible trace for hover label
        fig.add_trace(go.Scatter(
            x=[(x0 + x1)/2], y=[y_center],
            mode="markers", marker=dict(size=1, color="rgba(0,0,0,0)"),
            name=name, hovertemplate=f"<b>{name}</b><br>Coords: {chrom}:{x0}-{x1}<extra></extra>"
        ))

    # Exon 1 & Exon 3 (Constitutive)
    add_exon_box(ex1_start, ex1_end, 2, 0.8, color_constitutive, "Upstream Exon 1")
    add_exon_box(ex3_start, ex3_end, 2, 0.8, color_constitutive, "Downstream Exon 3")

    add_exon_box(ex1_start, ex1_end, -2, 0.8, color_constitutive, "Upstream Exon 1")
    add_exon_box(ex3_start, ex3_end, -2, 0.8, color_constitutive, "Downstream Exon 3")

    # Exon 2 (Target Skipped / Retained Exon)
    if event_type in ["SE", "MXE", "A5SS", "A3SS"]:
        add_exon_box(ex2_start, ex2_end, 2, 0.8, color_skipped, f"Alternative Exon 2 ({event_type})")

    # -------------------------------------------------------------
    # 3. Draw Junction Arcs & Read Counts
    # -------------------------------------------------------------
    # Inclusion Junction Arc (Exon 1 -> Exon 2 -> Exon 3)
    if event_type == "SE":
        inc_mid_x = (ex1_end + ex2_start) / 2
        fig.add_trace(go.Scatter(
            x=[ex1_end, inc_mid_x, ex2_start], y=[2.4, 3.2, 2.4],
            mode="lines+text", line=dict(color=color_inc_arc, width=3),
            text=["", f"Inclusion Reads: {inc_counts}", ""],
            textposition="top center", name="Inclusion Junction J1"
        ))
        
        inc_mid_x2 = (ex2_end + ex3_start) / 2
        fig.add_trace(go.Scatter(
            x=[ex2_end, inc_mid_x2, ex3_start], y=[2.4, 3.2, 2.4],
            mode="lines", line=dict(color=color_inc_arc, width=3),
            name="Inclusion Junction J2"
        ))

    # Exclusion Junction Arc (Exon 1 -> Exon 3 Skipping Exon 2)
    exc_mid_x = (ex1_end + ex3_start) / 2
    fig.add_trace(go.Scatter(
        x=[ex1_end, exc_mid_x, ex3_start], y=[-2.4, -3.4, -2.4],
        mode="lines+text", line=dict(color=color_exc_arc, width=3),
        text=["", f"Exclusion Reads: {exc_counts}", ""],
        textposition="bottom center", name="Exclusion Junction J3 (Skipping)"
    ))

    # -------------------------------------------------------------
    # Layout Formatting (Title positioned at the bottom below the plot)
    # -------------------------------------------------------------
    fig.update_layout(
        title=dict(
            text=f"<b>{gene_symbol} Exon Structure & Splicing Sashimi Plot</b> ({event_type} | ΔPSI = {delta_psi:+.2f})",
            x=0.5, y=0.01, xanchor="center", yanchor="bottom",
            font=dict(size=14, color="#1e293b")
        ),
        xaxis=dict(
            title=f"Genomic Coordinate ({chrom})",
            range=[x_min, x_max], showgrid=True, gridcolor="#f1f5f9"
        ),
        yaxis=dict(
            showticklabels=True,
            tickvals=[2, -2],
            ticktext=["Inclusion Isoform", "Exclusion Isoform"],
            range=[-4.5, 4.5], showgrid=False
        ),
        plot_bgcolor="white",
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(l=40, r=40, t=30, b=75),
        height=390
    )

    return fig

if __name__ == "__main__":
    fig = plot_exon_structure()
    print("Exon structure visualizer unit test OK!")
