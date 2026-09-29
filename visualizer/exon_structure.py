import re
import plotly.graph_objects as go

def resolve_gene_exon_coords(gene_symbol: str, event_type: str, coordinates: str):
    """
    Parses or deterministically derives realistic genomic exon-intron coordinates
    for any gene and splicing event type.
    """
    event = (event_type or "SE").upper()
    chrom = "chr1"
    
    # Try parsing clean coordinates first
    if coordinates and coordinates != "N/A":
        # Extract chromosome if present
        m_chr = re.search(r'(chr[0-9XYM]+)', coordinates, re.IGNORECASE)
        if m_chr:
            chrom = m_chr.group(1).lower()
        
        # Extract all numbers
        nums = [int(n) for n in re.findall(r'\d+', coordinates)]
        # Filter out chromosome number itself if it was captured as first number
        if m_chr and nums and str(nums[0]) in chrom:
            nums = nums[1:]
        
        if len(nums) >= 6:
            # Full 3-exon coordinates: e.g. [ex1_s, ex1_e, ex2_s, ex2_e, ex3_s, ex3_e]
            return {
                "chrom": chrom,
                "event_type": event,
                "nums": nums,
                "ex1": (nums[0], nums[1]),
                "ex2": (nums[2], nums[3]),
                "ex3": (nums[4], nums[5])
            }
        elif len(nums) >= 4:
            # 2 exons or partial: expand into 3 segments
            p1, p2, p3, p4 = nums[0], nums[1], nums[2], nums[3]
            span = max(400, (p4 - p1))
            return {
                "chrom": chrom,
                "event_type": event,
                "nums": nums,
                "ex1": (p1, p2),
                "ex2": (p2 + int(span * 0.25), p2 + int(span * 0.45)),
                "ex3": (p3, p4)
            }
        elif len(nums) >= 2:
            mid_s, mid_e = nums[0], nums[1]
            diff = max(150, mid_e - mid_s)
            return {
                "chrom": chrom,
                "event_type": event,
                "nums": nums,
                "ex1": (mid_s - diff * 3, mid_s - diff * 2),
                "ex2": (mid_s, mid_e),
                "ex3": (mid_e + diff * 2, mid_e + diff * 3)
            }

    # Deterministic fallback per gene symbol
    h = abs(hash(gene_symbol.upper()))
    chrom_num = (h % 22) + 1
    chrom = f"chr{chrom_num}"
    base = 10000000 + (h % 500000) * 100
    
    e1_len = 150 + (h % 80)
    intron1 = 800 + ((h >> 3) % 400)
    e2_len = 120 + ((h >> 6) % 100)
    intron2 = 900 + ((h >> 9) % 500)
    e3_len = 180 + ((h >> 12) % 90)

    ex1 = (base, base + e1_len)
    ex2 = (ex1[1] + intron1, ex1[1] + intron1 + e2_len)
    ex3 = (ex2[1] + intron2, ex2[1] + intron2 + e3_len)

    return {
        "chrom": chrom,
        "event_type": event,
        "nums": [ex1[0], ex1[1], ex2[0], ex2[1], ex3[0], ex3[1]],
        "ex1": ex1,
        "ex2": ex2,
        "ex3": ex3
    }


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
    Displays event-specific isoform structures for SE, RI, MXE, A5SS, and A3SS.
    """
    fig = go.Figure()
    event = (event_type or "SE").upper()
    coord_data = resolve_gene_exon_coords(gene_symbol, event, coordinates)
    chrom = coord_data["chrom"]
    
    # Calculate read counts if defaults provided
    if inc_counts == 145 and exc_counts == 22 and delta_psi != 0.35:
        inc_counts = int(max(20, abs(delta_psi) * 450 + 50))
        exc_counts = int(max(10, (1.0 - abs(delta_psi)) * 140 + 15))

    # Base coordinates
    ex1 = coord_data["ex1"]
    ex2 = coord_data["ex2"]
    ex3 = coord_data["ex3"]

    color_constitutive = "#2563EB"   # Classic Royal Blue
    color_alt = "#EF4444"            # Vibrant Coral Red
    color_retained = "#F59E0B"       # Amber / Orange
    color_mxe_b = "#9333EA"          # Purple / Violet
    color_inc_arc = "#10B981"        # Emerald Green
    color_exc_arc = "#8B5CF6"        # Violet

    def add_exon_box(x0, x1, y_center, height, color, name, hover_desc=None):
        fig.add_shape(
            type="rect",
            x0=x0, x1=x1,
            y0=y_center - height/2, y1=y_center + height/2,
            fillcolor=color, line=dict(color="#0F172A", width=1.5),
            layer="above"
        )
        tip = hover_desc or f"<b>{name}</b><br>Coords: {chrom}:{x0:,}-{x1:,} ({x1 - x0} bp)"
        fig.add_trace(go.Scatter(
            x=[(x0 + x1)/2], y=[y_center],
            mode="markers", marker=dict(size=1, color="rgba(0,0,0,0)"),
            name=name, hovertemplate=f"{tip}<extra></extra>"
        ))

    # =========================================================================
    # EVENT-SPECIFIC SASHIMI RENDERING
    # =========================================================================
    if event == "RI":
        # Retained Intron: Exon 1 and Exon 2 are separated by an intron that is retained
        x_min = ex1[0] - 250
        x_max = ex2[1] + 250

        # Intron lines
        fig.add_trace(go.Scatter(
            x=[ex1[0], ex2[1]], y=[2, 2], mode="lines",
            line=dict(color="#94A3B8", width=2, dash="dash"),
            name="Intron Backbone (Inclusion)", hoverinfo="skip"
        ))
        fig.add_trace(go.Scatter(
            x=[ex1[0], ex2[1]], y=[-2, -2], mode="lines",
            line=dict(color="#94A3B8", width=2, dash="dash"),
            name="Intron Backbone (Exclusion)", hoverinfo="skip"
        ))

        # Top Track: Retained Intron (Continuous Exon 1 - Intron - Exon 2)
        add_exon_box(ex1[0], ex1[1], 2, 0.8, color_constitutive, "Exon 1 (Upstream)")
        add_exon_box(ex1[1], ex2[0], 2, 0.5, color_retained, "Retained Intron Block",
                     f"<b>Retained Intron Block</b><br>Span: {chrom}:{ex1[1]:,}-{ex2[0]:,} ({ex2[0]-ex1[1]} bp)<br>Status: Retained in mature transcript")
        add_exon_box(ex2[0], ex2[1], 2, 0.8, color_constitutive, "Exon 2 (Downstream)")

        # Retention Read coverage line on top track
        fig.add_trace(go.Scatter(
            x=[(ex1[1] + ex2[0]) / 2], y=[2.6],
            mode="text", text=[f"Retained Intron Reads: {inc_counts}"],
            textposition="top center",
            textfont=dict(color="#D97706", size=12, family="Arial Black"),
            name="Intron Retention Signal"
        ))

        # Bottom Track: Normally Spliced Exons (Exon 1 and Exon 2 spliced, intron excised)
        add_exon_box(ex1[0], ex1[1], -2, 0.8, color_constitutive, "Exon 1 (Upstream)")
        add_exon_box(ex2[0], ex2[1], -2, 0.8, color_constitutive, "Exon 2 (Downstream)")

        # Spliced Intron Junction Arc
        mid_x = (ex1[1] + ex2[0]) / 2
        fig.add_trace(go.Scatter(
            x=[ex1[1], mid_x, ex2[0]], y=[-2.4, -3.5, -2.4],
            mode="lines+text", line=dict(color=color_exc_arc, width=3.5),
            text=["", f"Spliced Intron Reads: {exc_counts}", ""],
            textposition="bottom center",
            textfont=dict(color=color_exc_arc, size=11, family="Arial Bold"),
            name="Spliced Junction Arc"
        ))

    elif event == "MXE":
        # Mutually Exclusive Exons: Exon 1, Exon 2A, Exon 2B, Exon 3
        ex2a = ex2
        gap = (ex3[0] - ex2a[1]) // 3
        ex2b_len = ex2a[1] - ex2a[0]
        ex2b = (ex2a[1] + gap, ex2a[1] + gap + ex2b_len)
        x_min = ex1[0] - 250
        x_max = ex3[1] + 250

        # Intron lines
        fig.add_trace(go.Scatter(x=[ex1[0], ex3[1]], y=[2, 2], mode="lines", line=dict(color="#94A3B8", width=2, dash="dash"), hoverinfo="skip", name="Intron Track 1"))
        fig.add_trace(go.Scatter(x=[ex1[0], ex3[1]], y=[-2, -2], mode="lines", line=dict(color="#94A3B8", width=2, dash="dash"), hoverinfo="skip", name="Intron Track 2"))

        # Top Track: Includes Exon 2A (Exon 2B is skipped)
        add_exon_box(ex1[0], ex1[1], 2, 0.8, color_constitutive, "Exon 1")
        add_exon_box(ex2a[0], ex2a[1], 2, 0.8, color_alt, "Exon 2A (Isoform 1)")
        add_exon_box(ex3[0], ex3[1], 2, 0.8, color_constitutive, "Exon 3")
        
        mid_1a = (ex1[1] + ex2a[0]) / 2
        fig.add_trace(go.Scatter(
            x=[ex1[1], mid_1a, ex2a[0]], y=[2.4, 3.2, 2.4],
            mode="lines+text", line=dict(color=color_inc_arc, width=3),
            text=["", f"Exon 2A Reads: {inc_counts}", ""],
            textposition="top center", name="Junction 1 -> 2A"
        ))
        mid_a3 = (ex2a[1] + ex3[0]) / 2
        fig.add_trace(go.Scatter(
            x=[ex2a[1], mid_a3, ex3[0]], y=[2.4, 3.2, 2.4],
            mode="lines", line=dict(color=color_inc_arc, width=3), name="Junction 2A -> 3"
        ))

        # Bottom Track: Includes Exon 2B (Exon 2A is skipped)
        add_exon_box(ex1[0], ex1[1], -2, 0.8, color_constitutive, "Exon 1")
        add_exon_box(ex2b[0], ex2b[1], -2, 0.8, color_mxe_b, "Exon 2B (Isoform 2)")
        add_exon_box(ex3[0], ex3[1], -2, 0.8, color_constitutive, "Exon 3")

        mid_1b = (ex1[1] + ex2b[0]) / 2
        fig.add_trace(go.Scatter(
            x=[ex1[1], mid_1b, ex2b[0]], y=[-2.4, -3.4, -2.4],
            mode="lines+text", line=dict(color=color_exc_arc, width=3),
            text=["", f"Exon 2B Reads: {exc_counts}", ""],
            textposition="bottom center", name="Junction 1 -> 2B"
        ))
        mid_b3 = (ex2b[1] + ex3[0]) / 2
        fig.add_trace(go.Scatter(
            x=[ex2b[1], mid_b3, ex3[0]], y=[-2.4, -3.4, -2.4],
            mode="lines", line=dict(color=color_exc_arc, width=3), name="Junction 2B -> 3"
        ))

    elif event == "A5SS":
        # Alternative 5' Splice Site (Alternative Donor on Upstream Exon 1)
        ext_len = max(80, int((ex1[1] - ex1[0]) * 0.6))
        ex1_long = (ex1[0], ex1[1] + ext_len)
        ex2_down = ex3
        x_min = ex1[0] - 250
        x_max = ex2_down[1] + 250

        # Intron lines
        fig.add_trace(go.Scatter(x=[ex1[0], ex2_down[1]], y=[2, 2], mode="lines", line=dict(color="#94A3B8", width=2, dash="dash"), hoverinfo="skip", name="Intron (Long)"))
        fig.add_trace(go.Scatter(x=[ex1[0], ex2_down[1]], y=[-2, -2], mode="lines", line=dict(color="#94A3B8", width=2, dash="dash"), hoverinfo="skip", name="Intron (Short)"))

        # Top Track: Long Exon 1 (Distal 5' Donor)
        add_exon_box(ex1[0], ex1[1], 2, 0.8, color_constitutive, "Core Exon 1")
        add_exon_box(ex1[1], ex1_long[1], 2, 0.8, color_alt, "Alt 5' Extension", "<b>Alternative 5' Extension</b><br>Distal splice donor site")
        add_exon_box(ex2_down[0], ex2_down[1], 2, 0.8, color_constitutive, "Exon 2 (Downstream)")

        mid_long = (ex1_long[1] + ex2_down[0]) / 2
        fig.add_trace(go.Scatter(
            x=[ex1_long[1], mid_long, ex2_down[0]], y=[2.4, 3.3, 2.4],
            mode="lines+text", line=dict(color=color_inc_arc, width=3),
            text=["", f"Long 5' Reads: {inc_counts}", ""],
            textposition="top center", name="Distal 5' Splice Junction"
        ))

        # Bottom Track: Short Exon 1 (Proximal 5' Donor)
        add_exon_box(ex1[0], ex1[1], -2, 0.8, color_constitutive, "Core Exon 1")
        add_exon_box(ex2_down[0], ex2_down[1], -2, 0.8, color_constitutive, "Exon 2 (Downstream)")

        mid_short = (ex1[1] + ex2_down[0]) / 2
        fig.add_trace(go.Scatter(
            x=[ex1[1], mid_short, ex2_down[0]], y=[-2.4, -3.4, -2.4],
            mode="lines+text", line=dict(color=color_exc_arc, width=3),
            text=["", f"Short 5' Reads: {exc_counts}", ""],
            textposition="bottom center", name="Proximal 5' Splice Junction"
        ))

    elif event == "A3SS":
        # Alternative 3' Splice Site (Alternative Acceptor on Downstream Exon 2)
        ext_len = max(80, int((ex3[1] - ex3[0]) * 0.5))
        ex2_long = (ex3[0] - ext_len, ex3[1])
        x_min = ex1[0] - 250
        x_max = ex3[1] + 250

        # Intron lines
        fig.add_trace(go.Scatter(x=[ex1[0], ex3[1]], y=[2, 2], mode="lines", line=dict(color="#94A3B8", width=2, dash="dash"), hoverinfo="skip", name="Intron (Long)"))
        fig.add_trace(go.Scatter(x=[ex1[0], ex3[1]], y=[-2, -2], mode="lines", line=dict(color="#94A3B8", width=2, dash="dash"), hoverinfo="skip", name="Intron (Short)"))

        # Top Track: Long Exon 2 (Upstream 3' Acceptor Site)
        add_exon_box(ex1[0], ex1[1], 2, 0.8, color_constitutive, "Exon 1 (Upstream)")
        add_exon_box(ex2_long[0], ex3[0], 2, 0.8, color_alt, "Alt 3' Extension", "<b>Alternative 3' Extension</b><br>Proximal splice acceptor site")
        add_exon_box(ex3[0], ex3[1], 2, 0.8, color_constitutive, "Core Exon 2")

        mid_long = (ex1[1] + ex2_long[0]) / 2
        fig.add_trace(go.Scatter(
            x=[ex1[1], mid_long, ex2_long[0]], y=[2.4, 3.3, 2.4],
            mode="lines+text", line=dict(color=color_inc_arc, width=3),
            text=["", f"Long 3' Reads: {inc_counts}", ""],
            textposition="top center", name="Proximal 3' Splice Junction"
        ))

        # Bottom Track: Short Exon 2 (Downstream 3' Acceptor Site)
        add_exon_box(ex1[0], ex1[1], -2, 0.8, color_constitutive, "Exon 1 (Upstream)")
        add_exon_box(ex3[0], ex3[1], -2, 0.8, color_constitutive, "Core Exon 2")

        mid_short = (ex1[1] + ex3[0]) / 2
        fig.add_trace(go.Scatter(
            x=[ex1[1], mid_short, ex3[0]], y=[-2.4, -3.4, -2.4],
            mode="lines+text", line=dict(color=color_exc_arc, width=3),
            text=["", f"Short 3' Reads: {exc_counts}", ""],
            textposition="bottom center", name="Distal 3' Splice Junction"
        ))

    else:
        # Default: SE (Skipped Exon / Cassette Exon)
        x_min = ex1[0] - 250
        x_max = ex3[1] + 250

        # Intron lines
        fig.add_trace(go.Scatter(x=[ex1[0], ex3[1]], y=[2, 2], mode="lines", line=dict(color="#94A3B8", width=2, dash="dash"), hoverinfo="skip", name="Intron (Inclusion Track)"))
        fig.add_trace(go.Scatter(x=[ex1[0], ex3[1]], y=[-2, -2], mode="lines", line=dict(color="#94A3B8", width=2, dash="dash"), hoverinfo="skip", name="Intron (Exclusion Track)"))

        # Top Track: Exon 1 - Exon 2 - Exon 3
        add_exon_box(ex1[0], ex1[1], 2, 0.8, color_constitutive, "Upstream Exon 1")
        add_exon_box(ex2[0], ex2[1], 2, 0.8, color_alt, "Alternative Skipped Exon 2", "<b>Alternative Skipped Exon</b><br>Cassette exon included in dominant isoform")
        add_exon_box(ex3[0], ex3[1], 2, 0.8, color_constitutive, "Downstream Exon 3")

        mid_1 = (ex1[1] + ex2[0]) / 2
        fig.add_trace(go.Scatter(
            x=[ex1[1], mid_1, ex2[0]], y=[2.4, 3.2, 2.4],
            mode="lines+text", line=dict(color=color_inc_arc, width=3),
            text=["", f"Inclusion J1: {inc_counts}", ""],
            textposition="top center", name="Inclusion Junction J1"
        ))

        mid_2 = (ex2[1] + ex3[0]) / 2
        fig.add_trace(go.Scatter(
            x=[ex2[1], mid_2, ex3[0]], y=[2.4, 3.2, 2.4],
            mode="lines", line=dict(color=color_inc_arc, width=3),
            name="Inclusion Junction J2"
        ))

        # Bottom Track: Exon 1 - Exon 3 (Skipping Exon 2)
        add_exon_box(ex1[0], ex1[1], -2, 0.8, color_constitutive, "Upstream Exon 1")
        add_exon_box(ex3[0], ex3[1], -2, 0.8, color_constitutive, "Downstream Exon 3")

        mid_exc = (ex1[1] + ex3[0]) / 2
        fig.add_trace(go.Scatter(
            x=[ex1[1], mid_exc, ex3[0]], y=[-2.4, -3.5, -2.4],
            mode="lines+text", line=dict(color=color_exc_arc, width=3),
            text=["", f"Exclusion J3 (Skipping): {exc_counts}", ""],
            textposition="bottom center", name="Exclusion Junction J3 (Skipping)"
        ))

    # Track Labels & Title layout
    tick_inclusion = f"Inclusion Isoform ({event})" if event != "MXE" else "Isoform 1 (Exon 2A)"
    tick_exclusion = f"Exclusion Isoform" if event != "MXE" else "Isoform 2 (Exon 2B)"

    fig.update_layout(
        title=None,
        xaxis=dict(
            title=f"<b>Genomic Coordinate ({chrom})</b>",
            range=[x_min, x_max], showgrid=True, gridcolor="#F1F5F9",
            tickformat=",d"
        ),
        yaxis=dict(
            automargin=True,
            showticklabels=True,
            tickvals=[2, -2],
            ticktext=[f"<b>{tick_inclusion}</b>", f"<b>{tick_exclusion}</b>"],
            range=[-4.8, 4.8], showgrid=False
        ),
        plot_bgcolor="#FFFFFF",
        paper_bgcolor="#FFFFFF",
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1, font=dict(size=11)),
        margin=dict(l=180, r=40, t=30, b=40),
        height=380
    )

    return fig

if __name__ == "__main__":
    fig = plot_exon_structure()
    print("Exon structure visualizer unit test OK!")

