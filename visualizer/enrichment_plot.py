"""
GenSplice-Agent Plotly Enrichment Bar Visualizer
Renders light-mode dynamic horizontal bar charts for GO Terms & KEGG Pathways.
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
            text="No significant enriched terms found for significant splicing targets.",
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
