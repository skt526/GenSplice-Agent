"""
GenSplice-Agent Plotly Enrichment Bar & Dot Plot Visualizers
Renders light-mode dynamic horizontal bar charts and bubble dot plots for GO Terms & KEGG Pathways.
"""

import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

def build_enrichment_chart(df_enr: pd.DataFrame, title: str, bar_color: str = "#3B82F6") -> go.Figure:
    """
    Builds a horizontal Plotly bar chart for enrichment terms sorted by -log10(p-value).
    """
    if df_enr is None or df_enr.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No significant enriched terms found for the selected quadrant.",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="#64748B")
        )
        fig.update_layout(
            paper_bgcolor="#FFFFFF",
            plot_bgcolor="#F8FAFC",
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            height=300
        )
        return fig
        
    df_sorted = df_enr.sort_values(by="P-value", ascending=True).head(10).iloc[::-1]
    
    df_sorted["Display_Term"] = df_sorted["Term"].apply(
        lambda t: t[:50] + "..." if len(str(t)) > 53 else str(t)
    )
    
    fig = px.bar(
        df_sorted,
        x="log_p",
        y="Display_Term",
        orientation="h",
        labels={"log_p": "-log₁₀(p-value)", "Display_Term": "Enriched Term / Pathway"},
        hover_data={"Term": True, "P-value": ":.4f", "Adjusted P-value": ":.4f", "Overlap": True, "Genes": True},
        title=f"<b>{title}</b>"
    )
    
    fig.update_traces(
        marker_color=bar_color,
        marker_line_color="#1E40AF",
        marker_line_width=1,
        opacity=0.85
    )
    
    fig.update_layout(
        paper_bgcolor="#FFFFFF",
        plot_bgcolor="#F8FAFC",
        font=dict(family="Inter, Roboto, sans-serif", color="#0F172A", size=12),
        margin=dict(l=20, r=20, t=50, b=40),
        xaxis=dict(
            title="<b>-log₁₀ (p-value)</b>",
            gridcolor="#E2E8F0",
            zerolinecolor="#CBD5E1"
        ),
        yaxis=dict(
            title="",
            gridcolor="#E2E8F0"
        ),
        height=380
    )
    
    return fig

def build_enrichment_dot_plot(
    df_enr: pd.DataFrame, 
    title: str,
    color_scale: str = "Purples_r",
    border_color: str = "#6B21A8"
) -> go.Figure:
    """
    Builds a dynamic Dot / Bubble plot for single-quadrant enrichment results.
    X-axis: Gene Ratio (Overlap Ratio)
    Y-axis: Enriched Term
    Size: Gene Count Overlap
    Color: Adjusted P-value
    """
    if df_enr is None or df_enr.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No significant enriched terms found for this quadrant.",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="#64748B")
        )
        fig.update_layout(
            paper_bgcolor="#FFFFFF",
            plot_bgcolor="#F8FAFC",
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            height=300
        )
        return fig

    df_sorted = df_enr.sort_values(by="P-value", ascending=True).head(10).iloc[::-1].copy()
    
    def parse_overlap_num(val):
        try:
            return int(str(val).split("/")[0])
        except Exception:
            return 5

    def parse_gene_ratio(val):
        try:
            parts = str(val).split("/")
            return round(float(parts[0]) / float(parts[1]), 3)
        except Exception:
            return 0.1

    df_sorted["Gene_Count"] = df_sorted["Overlap"].apply(parse_overlap_num)
    df_sorted["Gene_Ratio"] = df_sorted["Overlap"].apply(parse_gene_ratio)
    df_sorted["Display_Term"] = df_sorted["Term"].apply(
        lambda t: t[:50] + "..." if len(str(t)) > 53 else str(t)
    )

    fig = px.scatter(
        df_sorted,
        x="Gene_Ratio",
        y="Display_Term",
        size="Gene_Count",
        color="Adjusted P-value",
        color_continuous_scale=color_scale,
        size_max=22,
        labels={"Gene_Ratio": "Gene Ratio (Overlap Ratio)", "Display_Term": "Enriched Term / Pathway"},
        hover_data={"Term": True, "P-value": ":.4f", "Adjusted P-value": ":.4f", "Overlap": True, "Genes": True},
        title=f"<b>{title}</b>"
    )

    fig.update_traces(
        marker=dict(line=dict(width=1.5, color=border_color))
    )

    fig.update_layout(
        paper_bgcolor="#FFFFFF",
        plot_bgcolor="#F8FAFC",
        font=dict(family="Inter, Roboto, sans-serif", color="#0F172A", size=12),
        margin=dict(l=20, r=20, t=50, b=40),
        xaxis=dict(
            title="<b>Gene Ratio (Overlap Ratio)</b>",
            gridcolor="#E2E8F0",
            zerolinecolor="#CBD5E1"
        ),
        yaxis=dict(
            title="",
            gridcolor="#E2E8F0"
        ),
        height=380
    )

    return fig

def build_combined_quadrant_dot_plot(
    df_q1: pd.DataFrame, 
    df_q2: pd.DataFrame, 
    df_q4: pd.DataFrame, 
    title: str = "Comparative GO Term Dot + Bubble Plot across Quadrants"
) -> go.Figure:
    """
    Builds a multi-quadrant comparative Dot / Bubble Plot:
    X-axis: Quadrants (Q1, Q2, Q4)
    Y-axis: Enriched GO Terms
    Bubble Size: Gene Count Overlap
    Bubble Color: Adjusted P-value
    """
    frames = []
    for df, q_label in [
        (df_q1, "Q1"),
        (df_q2, "Q2"),
        (df_q4, "Q4")
    ]:
        if df is not None and not df.empty:
            sub = df.head(8).copy()
            sub["Quadrant"] = q_label
            frames.append(sub)
            
    if not frames:
        fig = go.Figure()
        fig.add_annotation(
            text="No enriched GO terms found across Q1, Q2, Q4.",
            xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="#64748B")
        )
        fig.update_layout(paper_bgcolor="#FFFFFF", plot_bgcolor="#F8FAFC", height=350)
        return fig
        
    combined = pd.concat(frames, ignore_index=True)
    
    def parse_overlap_num(val):
        try:
            return int(str(val).split("/")[0])
        except Exception:
            return 5

    combined["Gene_Count"] = combined["Overlap"].apply(parse_overlap_num)
    combined["Display_Term"] = combined["Term"].apply(
        lambda t: t[:45] + "..." if len(str(t)) > 48 else str(t)
    )
    
    # Sort terms by minimum p-value so most significant are at top
    term_order = combined.groupby("Display_Term")["P-value"].min().sort_values(ascending=False).index.tolist()
    
    fig = px.scatter(
        combined,
        x="Quadrant",
        y="Display_Term",
        size="Gene_Count",
        color="Adjusted P-value",
        color_continuous_scale="Purples_r",
        size_max=24,
        category_orders={"Display_Term": term_order, "Quadrant": ["Q1", "Q2", "Q4"]},
        labels={"Adjusted P-value": "Adjusted P-value", "Display_Term": "GO Biological Process Term", "Quadrant": "Quadrant Category", "Gene_Count": "Gene Count"},
        hover_data={"Term": True, "P-value": ":.4f", "Adjusted P-value": ":.4f", "Overlap": True, "Genes": True},
        title=f"<b>{title}</b>"
    )

    fig.update_traces(
        marker=dict(line=dict(width=1.5, color="#4C1D95"))
    )

    fig.update_layout(
        paper_bgcolor="#FFFFFF",
        plot_bgcolor="#F8FAFC",
        font=dict(family="Inter, Roboto, sans-serif", color="#0F172A", size=12),
        margin=dict(l=20, r=20, t=50, b=50),
        xaxis=dict(
            title="<b>Quadrant Category (X-axis: Q1, Q2, Q4)</b>",
            gridcolor="#E2E8F0",
            zerolinecolor="#CBD5E1"
        ),
        yaxis=dict(
            title="",
            gridcolor="#E2E8F0"
        ),
        height=520
    )

    return fig
