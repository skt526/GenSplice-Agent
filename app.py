"""
GenSplice-Agent Streamlit Dashboard
Light Mode Theme with Dynamic Threshold-Reactive GO Term & KEGG Pathway Analysis
- Outside plot header: 'Color Classification' (with Q1-Q4 descriptions)
- Inside plot legend title: 'On/Off' (strictly Q1, Q2, Q3, Q4)
- GO Term & KEGG Pathway Sections focused on Q1 (Both DEG & Splicing)
- Download Table CSV Buttons on the right side of GO & KEGG headers
- Main Result Page displays ONLY the interactive graphs
- Generate GO/KEGG Button re-calculates GO & KEGG with gradient button animations & toast visual feedback
- Prominent Chart Type Radio Selectors (Bar & Dot/Bubble Plot)
- Always-Visible KEGG Pathway Explorer Dropdown
"""

import os
import streamlit as st
import polars as pl
import pandas as pd

from config import (
    DEFAULT_LOG2FC_CUTOFF, DEFAULT_DELTA_PSI_CUTOFF,
    DEFAULT_DEG_FDR_CUTOFF, DEFAULT_AS_FDR_CUTOFF,
    THEME_BG, THEME_CARD_BG, THEME_TEXT, THEME_MUTED, THEME_BORDER,
    QUADRANT_COLORS
)
from core.deg_loader import load_deg_data
from core.rmats_loader import load_rmats_data, select_primary_splicing_events
from core.merger import merge_deg_and_rmats, get_quadrant_kpis
from core.enrichment import fetch_enrichment
from visualizer.quadrant_plot import build_quadrant_plot
from visualizer.enrichment_plot import build_enrichment_chart, build_enrichment_dot_plot
import urllib.parse
from visualizer.report_exporter import export_html_report
from core.ai_summary import generate_biological_insights, fetch_ncbi_gene_summary, fetch_pubmed_literature, generate_detailed_bio_prompt
from core.isoform_annotator import annotate_isoform_events

# Streamlit Page Config
st.set_page_config(
    page_title="GenSplice-Agent Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for Premium Light Mode Theme & Dynamic Button Action Effects
st.markdown(f"""
<style>
    .stApp {{
        background-color: {THEME_BG};
        color: {THEME_TEXT};
    }}
    .main-header {{
        font-size: 2.4rem;
        font-weight: 900;
        letter-spacing: -0.5px;
        color: #0F172A;
        margin-bottom: 0.2rem;
    }}
    .sub-header {{
        color: #64748B;
        font-size: 1rem;
        margin-bottom: 1.5rem;
    }}
    .kpi-box {{
        background-color: {THEME_CARD_BG};
        border-radius: 12px;
        padding: 16px;
        border: 1px solid {THEME_BORDER};
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
        text-align: center;
        margin-top: 16px;
    }}
    .kpi-val {{
        font-size: 2.2rem;
        font-weight: 800;
        margin-top: 4px;
    }}
    .kpi-lbl {{
        font-size: 0.8rem;
        color: #64748B;
        font-weight: 700;
        text-transform: uppercase;
        letter-spacing: 0.5px;
    }}
    .control-panel-right {{
        background-color: {THEME_CARD_BG};
        border-radius: 14px;
        padding: 20px;
        border: 2px solid #3B82F6;
        box-shadow: 0 4px 12px rgba(59, 130, 246, 0.1);
    }}
    .stButton>button {{
        border-radius: 8px;
        font-weight: 700;
        background-color: #FFFFFF;
        color: #0F172A;
        border: 1px solid #CBD5E1;
    }}
    .stButton>button:hover {{
        background-color: #F1F5F9;
        border-color: #3B82F6;
        color: #3B82F6;
    }}
    .btn-generate-enrichment>button {{
        background: linear-gradient(135deg, #A855F7 0%, #7C3AED 100%) !important;
        color: #FFFFFF !important;
        border: none !important;
        font-weight: 800 !important;
        font-size: 1.05rem !important;
        border-radius: 10px !important;
        padding: 12px 18px !important;
        box-shadow: 0 4px 14px rgba(168, 85, 247, 0.4) !important;
        transition: all 0.25s cubic-bezier(0.4, 0, 0.2, 1) !important;
    }}
    .btn-generate-enrichment>button:hover {{
        background: linear-gradient(135deg, #9333EA 0%, #6D28D9 100%) !important;
        transform: translateY(-2px) scale(1.02) !important;
        box-shadow: 0 6px 20px rgba(168, 85, 247, 0.6) !important;
    }}
    .btn-generate-enrichment>button:active {{
        transform: translateY(1px) scale(0.97) !important;
    }}
    .pathway-box {{
        background-color: #F8FAFC;
        border-left: 5px solid #A855F7;
        padding: 16px;
        border-radius: 8px;
        margin-top: 12px;
        margin-bottom: 16px;
    }}
</style>
""", unsafe_allow_html=True)

# Sidebar Configuration
st.sidebar.image("https://img.icons8.com/color/96/dna-helix.png", width=64)
st.sidebar.title("GenSplice-Agent")
st.sidebar.caption("Light Theme Edition")
st.sidebar.markdown("---")

# 1. Data Input Paths
st.sidebar.header("📁 Data Inputs")
deg_path = st.sidebar.text_input("DEG File Path (CSV/TSV)", value="outputs/03_deg/deg_result.csv")
rmats_dir = st.sidebar.text_input("rMATS Directory Path", value="outputs/04_rmats")

# Fallback to test data if default outputs don't exist yet
if not os.path.exists(deg_path) and os.path.exists("test_data/deg_result.csv"):
    deg_path = "test_data/deg_result.csv"
if not os.path.exists(rmats_dir) and os.path.exists("test_data"):
    rmats_dir = "test_data"

st.sidebar.markdown("---")

# 2. Sidebar Color Options
st.sidebar.header("🎨 Plot Color Settings")
color_option = st.sidebar.radio("Color Palette", ["By Quadrant (Red/Blue/Purple)", "By Splicing Event Type"])
color_by = "quadrant" if color_option.startswith("By Quadrant") else "event_type"

st.sidebar.markdown("---")

# Data Loading Engine
@st.cache_data(ttl=60)
def load_and_process_raw(deg_f, rmats_d):
    if not os.path.exists(deg_f) or not os.path.exists(rmats_d):
        return pl.DataFrame(), pl.DataFrame()
    df_deg = load_deg_data(deg_f)
    df_rmats = load_rmats_data(rmats_d)
    primary_rmats = select_primary_splicing_events(df_rmats)
    return df_deg, primary_rmats

df_deg_raw, df_rmats_raw = load_and_process_raw(deg_path, rmats_dir)

# Main Dashboard Layout
st.markdown('<div class="main-header">🧬 GenSplice-Agent Platform</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-header">Light Theme • Q1 (Both DEG & Splicing) Focus • Red (DEG) | Blue (Splicing) | Purple (Both)</div>', unsafe_allow_html=True)

if df_deg_raw.height == 0 and df_rmats_raw.height == 0:
    st.warning("⚠️ No valid DEG or rMATS data found. Please run `./GenSplice` or `./test` to generate output data.")
    st.stop()

# Main Plot Area and Right-Side Controls Layout
plot_col, control_col = st.columns([3, 1])

# Threshold Sliders on the Right Side
with control_col:
    st.markdown('<div class="control-panel-right">', unsafe_allow_html=True)
    st.markdown("### 🎛️ Threshold Controls (Right)")
    
    log2fc_cutoff = st.slider("Log₂FC Cutoff (DEG: Red)", 0.1, 3.0, DEFAULT_LOG2FC_CUTOFF, 0.05, key="log2fc_slider")
    delta_psi_cutoff = st.slider("ΔPSI Cutoff (Splicing: Blue)", 0.01, 0.5, DEFAULT_DELTA_PSI_CUTOFF, 0.01, key="psi_slider")
    deg_fdr_cutoff = st.select_slider("DEG FDR Cutoff", options=[0.001, 0.01, 0.05, 0.1], value=DEFAULT_DEG_FDR_CUTOFF, key="deg_fdr_slider")
    as_fdr_cutoff = st.select_slider("rMATS FDR Cutoff", options=[0.001, 0.01, 0.05, 0.1], value=DEFAULT_AS_FDR_CUTOFF, key="as_fdr_slider")
    
    st.markdown("---")
    
    # Dynamic Merger calculation based on active sliders
    df_merged = merge_deg_and_rmats(
        df_deg=df_deg_raw,
        df_rmats=df_rmats_raw,
        log2fc_cutoff=log2fc_cutoff,
        delta_psi_cutoff=delta_psi_cutoff,
        deg_fdr_cutoff=deg_fdr_cutoff,
        as_fdr_cutoff=as_fdr_cutoff
    )

    # Dynamic CSV Export Button placed directly below Threshold Sliders (Includes Functional Impairment Risk Scores)
    df_merged_annotated = annotate_isoform_events(df_merged)
    csv_bytes = df_merged_annotated.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="📥 Download Filtered Gene List (.csv)",
        data=csv_bytes,
        file_name=f"gensplice_genes_FC{log2fc_cutoff}_dPSI{delta_psi_cutoff}.csv",
        mime="text/csv",
        help="Click to download the current filtered gene list matching these exact threshold cutoffs"
    )

    st.markdown("---")
    # Outside plot header strictly: Color Classification
    st.markdown("### 📌 Color Classification")
    st.markdown("""
    - 🟣 **Q1:** Both DEG & Splicing
    - 🔵 **Q2:** Splicing Only
    - ⚪ **Q3:** Invariant Background
    - 🔴 **Q4:** DEG Only
    """)
    
    st.markdown("---")
    # Generate GO/KEGG Button directly below Color Classification
    st.markdown('<div class="btn-generate-enrichment">', unsafe_allow_html=True)
    gen_enrich_btn = st.button(
        "🚀 Generate GO / KEGG",
        use_container_width=True,
        help="Click to compute/refresh GO Term and KEGG Pathway analysis based on the updated Q1 gene list matching current threshold cutoffs"
    )
    st.markdown('</div>', unsafe_allow_html=True)
    
    st.markdown('</div>', unsafe_allow_html=True)

# Manage Session State for GO & KEGG Enrichment Calculation
q1_current_genes = df_merged.filter(pl.col("quadrant") == "Q1").select("geneSymbol").to_series().to_list()
if not q1_current_genes:
    q1_current_genes = df_merged.select("geneSymbol").to_series().to_list()

if "go_results_df" not in st.session_state or gen_enrich_btn:
    if gen_enrich_btn:
        st.toast("⚡ Recalculating GO Term & KEGG Pathway Enrichments...", icon="🧬")
    with st.spinner("⚡ Re-calculating GO & KEGG Pathways for active Q1 genes..."):
        st.session_state["go_results_df"] = fetch_enrichment(q1_current_genes, gene_sets=["GO_Biological_Process_2023"], top_n=10)
        st.session_state["kegg_results_df"] = fetch_enrichment(q1_current_genes, gene_sets=["KEGG_2021_Human"], top_n=10)
        st.session_state["enrichment_q1_count"] = len(q1_current_genes)
        st.session_state["enrichment_cutoffs_str"] = f"Log₂FC ≥ {log2fc_cutoff:.2f}, ΔPSI ≥ {delta_psi_cutoff:.2f}"
    if gen_enrich_btn:
        st.toast(f"✨ GO & KEGG Enrichment Updated! ({len(q1_current_genes)} Q1 genes)", icon="🚀")

# Export HTML Report Button
export_col1, export_col2 = st.columns([3, 1])
with export_col2:
    if st.button("🌐 Export Interactive Chrome HTML Report"):
        html_out_path = export_html_report(
            df_merged=df_merged,
            output_html_path="outputs/gensplice_report.html",
            log2fc_cutoff=log2fc_cutoff,
            delta_psi_cutoff=delta_psi_cutoff,
            deg_fdr_cutoff=deg_fdr_cutoff,
            as_fdr_cutoff=as_fdr_cutoff
        )
        st.success("✔ Light Mode HTML Report generated!")
        with open(html_out_path, "r", encoding="utf-8") as f:
            html_bytes = f.read().encode("utf-8")
        st.download_button(
            label="💾 Download HTML Report",
            data=html_bytes,
            file_name="gensplice_report.html",
            mime="text/html"
        )

# Main Navigation Tabs (Updated Section Titles highlighting Q1)
tab1, tab2, tab3, tab4 = st.tabs([
    "📊 4-Quadrant Cross-Plot",
    "🧬 GO Term Analysis (Q1: Both DEG & Splicing)",
    "🛤️ KEGG Pathway Analysis (Q1: Both DEG & Splicing)",
    "📋 Data Explorer & Target Selector"
])

# Tab 1: Quadrant Plot & Gene Count KPI Cards BELOW Plot
with tab1:
    st.subheader("Interactive 4-Quadrant Cross-Plot")
    fig_quad = build_quadrant_plot(df_merged, log2fc_cutoff, delta_psi_cutoff, color_by=color_by)
    st.plotly_chart(
        fig_quad,
        use_container_width=True,
        config={
            "editable": True,
            "edits": {"shapePosition": True},
            "displayModeBar": True
        }
    )

    # KPI Summary Cards Positioned DIRECTLY BELOW 4-Quadrant Plot
    st.markdown("---")
    st.markdown("### 📊 Gene Count Summary Metrics (Updated in Real Time)")
    kpis = get_quadrant_kpis(df_merged)
    col1, col2, col3, col4, col5 = st.columns(5)

    with col1:
        st.markdown(f'<div class="kpi-box"><div class="kpi-lbl">Total Analyzed</div><div class="kpi-val" style="color:#0F172A;">{kpis["total"]}</div></div>', unsafe_allow_html=True)
    with col2:
        st.markdown(f'<div class="kpi-box" style="border-top:4px solid #3B82F6;"><div class="kpi-lbl" style="color:#3B82F6;">Q2: Splicing Target</div><div class="kpi-val" style="color:#3B82F6;">{kpis["Q2"]}</div></div>', unsafe_allow_html=True)
    with col3:
        st.markdown(f'<div class="kpi-box" style="border-top:4px solid #A855F7;"><div class="kpi-lbl" style="color:#A855F7;">Q1: Dual Responders</div><div class="kpi-val" style="color:#A855F7;">{kpis["Q1"]}</div></div>', unsafe_allow_html=True)
    with col4:
        st.markdown(f'<div class="kpi-box" style="border-top:4px solid #EF4444;"><div class="kpi-lbl" style="color:#EF4444;">Q4: DEG Only</div><div class="kpi-val" style="color:#EF4444;">{kpis["Q4"]}</div></div>', unsafe_allow_html=True)
    with col5:
        st.markdown(f'<div class="kpi-box" style="border-top:4px solid #94A3B8;"><div class="kpi-lbl" style="color:#64748B;">Q3: Invariant</div><div class="kpi-val" style="color:#64748B;">{kpis["Q3"]}</div></div>', unsafe_allow_html=True)

# Helper function to render GO Term Section (Default Q1, ONLY graph on page + Download button on header right)
def render_go_section(df_merged_data):
    st.subheader("🧬 GO Term Analysis (Q1: Both DEG & Splicing)")
    st.markdown("Perform Gene Ontology enrichment focusing primarily on **Q1 (Dual Responders: Both DEG & Splicing)** genes.")
    
    st.info(f"💡 **Current GO Analysis Scope:** {st.session_state.get('enrichment_q1_count', 0)} Q1 genes ({st.session_state.get('enrichment_cutoffs_str', '')}). Adjust threshold sliders on the right panel and click **🚀 Generate GO / KEGG** to re-calculate.")
    
    df_go = st.session_state.get("go_results_df", pd.DataFrame())

    # Prominently Visible Controls Bar with Download CSV Button on the right
    ctrl_col1, ctrl_col2, ctrl_col3 = st.columns([2, 1, 1])
    with ctrl_col1:
        q_filter_go = st.multiselect(
            "Select Quadrants for GO Analysis",
            options=["Q1", "Q2", "Q3", "Q4"],
            default=["Q1"],
            key="q_filter_go_select"
        )
    with ctrl_col2:
        chart_style = st.radio(
            "📊 Chart Style",
            ["Bar Chart", "Dot / Bubble Plot"],
            horizontal=True,
            key="go_chart_type"
        )
    with ctrl_col3:
        if df_go is not None and not df_go.empty:
            go_csv_bytes = df_go[["Term", "Overlap", "P-value", "Adjusted P-value", "Genes"]].to_csv(index=False).encode('utf-8')
            st.download_button(
                label="📥 Download GO Table (.csv)",
                data=go_csv_bytes,
                file_name=f"go_terms_Q1_FC{log2fc_cutoff:.2f}_dPSI{delta_psi_cutoff:.2f}.csv",
                mime="text/csv",
                key="dl_go_csv_btn",
                use_container_width=True
            )
    
    if not q_filter_go:
        st.info("Please select at least one quadrant to analyze GO terms.")
        return
        
    if chart_style == "Dot / Bubble Plot":
        fig_go = build_enrichment_dot_plot(df_go, f"GO Biological Process Dot Plot ({', '.join(q_filter_go)})")
    else:
        fig_go = build_enrichment_chart(df_go, f"Top GO Biological Processes ({', '.join(q_filter_go)})", bar_color="#A855F7")
        
    # ONLY display the graph on the result page
    st.plotly_chart(fig_go, use_container_width=True)

# Helper function to render KEGG Section (Default Q1, ONLY graph on page + Download button on header right)
def render_kegg_section(df_merged_data):
    st.subheader("🛤️ KEGG Pathway Analysis (Q1: Both DEG & Splicing)")
    st.markdown("Perform KEGG pathway enrichment focusing on **Q1 (Dual Responders: Both DEG & Splicing)** target genes.")
    
    st.info(f"💡 **Current KEGG Analysis Scope:** {st.session_state.get('enrichment_q1_count', 0)} Q1 genes ({st.session_state.get('enrichment_cutoffs_str', '')}). Adjust threshold sliders on the right panel and click **🚀 Generate GO / KEGG** to re-calculate.")
    
    df_kegg = st.session_state.get("kegg_results_df", pd.DataFrame())

    # Prominently Visible Controls Bar with Download CSV Button on the right
    ctrl_col1, ctrl_col2, ctrl_col3 = st.columns([2, 1, 1])
    with ctrl_col1:
        q_filter_kegg = st.multiselect(
            "Select Quadrants for KEGG Analysis",
            options=["Q1", "Q2", "Q3", "Q4"],
            default=["Q1"],
            key="q_filter_kegg_select"
        )
    with ctrl_col2:
        chart_style = st.radio(
            "📊 Chart Style",
            ["Bar Chart", "Dot / Bubble Plot"],
            horizontal=True,
            key="kegg_chart_type"
        )
    with ctrl_col3:
        if df_kegg is not None and not df_kegg.empty:
            kegg_csv_bytes = df_kegg[["Term", "Overlap", "P-value", "Adjusted P-value", "Genes"]].to_csv(index=False).encode('utf-8')
            st.download_button(
                label="📥 Download KEGG Table (.csv)",
                data=kegg_csv_bytes,
                file_name=f"kegg_pathways_Q1_FC{log2fc_cutoff:.2f}_dPSI{delta_psi_cutoff:.2f}.csv",
                mime="text/csv",
                key="dl_kegg_csv_btn",
                use_container_width=True
            )
    
    if not q_filter_kegg:
        st.info("Please select at least one quadrant to analyze KEGG pathways.")
        return
        
    if chart_style == "Dot / Bubble Plot":
        fig_kegg = build_enrichment_dot_plot(df_kegg, f"KEGG Pathway Dot Plot ({', '.join(q_filter_kegg)})")
    else:
        fig_kegg = build_enrichment_chart(df_kegg, f"Top KEGG Pathways ({', '.join(q_filter_kegg)})", bar_color="#A855F7")
        
    # ONLY display the graph on the result page
    st.plotly_chart(fig_kegg, use_container_width=True)

    # Prominently Visible KEGG Pathway Explorer Box on Main Canvas
    st.markdown("---")
    st.markdown("### 🎯 KEGG Pathway Interactive Gene Explorer")
    st.markdown("Select an enriched pathway below to view all corresponding target genes and transcriptomics metrics:")
    
    kegg_terms = df_kegg["Term"].tolist() if (df_kegg is not None and not df_kegg.empty) else ["No enriched pathways found"]
    selected_pathway = st.selectbox(
        "🔍 Choose KEGG Pathway to Inspect Target Genes",
        options=kegg_terms,
        key="kegg_pathway_dropdown_main"
    )
    
    if selected_pathway and df_kegg is not None and not df_kegg.empty:
        matched_rows = df_kegg[df_kegg["Term"] == selected_pathway]
        if not matched_rows.empty:
            matched_row = matched_rows.iloc[0]
            genes_str = str(matched_row.get("Genes", ""))
            g_list = [g.strip() for g in genes_str.split(";") if g.strip() and g.strip() != "N/A"]
            
            st.markdown(f'<div class="pathway-box"><b>Selected Pathway:</b> {selected_pathway}<br><b>P-value:</b> {matched_row["P-value"]:.2e} | <b>Adjusted P-value:</b> {matched_row["Adjusted P-value"]:.2e} | <b>Overlap Count:</b> {matched_row["Overlap"]}</div>', unsafe_allow_html=True)
            
            if g_list:
                pathway_genes_df = df_merged_data.filter(pl.col("geneSymbol").is_in(g_list))
                if pathway_genes_df.height > 0:
                    st.markdown(f"**Target Genes belonging to `{selected_pathway}`:**")
                    st.dataframe(
                        pathway_genes_df.to_pandas()[[
                            "geneSymbol", "gene_id", "quadrant", "event_type",
                            "delta_psi", "log2FoldChange", "as_fdr", "deg_fdr", "coordinates"
                        ]],
                        use_container_width=True,
                        hide_index=True
                    )

# Tab 2: GO Term Analysis (Q1 Focused)
with tab2:
    render_go_section(df_merged)

# Tab 3: KEGG Pathway Analysis (Q1 Focused + Interactive Explorer)
with tab3:
    render_kegg_section(df_merged)

# Tab 4: Data Explorer
with tab4:
    st.subheader("Transcriptomics Data Explorer")
    quad_filter = st.multiselect("Filter by Quadrant", options=["Q1", "Q2", "Q3", "Q4"], default=["Q1", "Q2"], key="explorer_quad_filter")
    
    filtered_df = df_merged.filter(pl.col("quadrant").is_in(quad_filter)) if quad_filter else df_merged
    pdf_view = filtered_df.to_pandas()
    
    st.dataframe(
        pdf_view[[
            "geneSymbol", "gene_id", "quadrant", "event_type", "delta_psi",
            "log2FoldChange", "as_fdr", "deg_fdr", "coordinates"
        ]],
        use_container_width=True,
        hide_index=True
    )

    st.markdown("---")
    st.subheader("🧬 Event-Level Isoform Annotation & NMD Prediction Matrix")
    st.caption("Includes Transcript ID, Event Type (SE/RI/MXE/A5SS/A3SS), Exon Coordinates, CDS Frame, PTC Position, NMD Prediction, Protein Domain Overlap & Subcellular Localization Consequences.")

    isoform_annot_df = annotate_isoform_events(filtered_df)
    
    st.dataframe(
        isoform_annot_df[[
            "geneSymbol", "impairment_tier", "primary_dysfunction_cause", "event_type", "transcript_id", "coordinates",
            "delta_psi", "log2FoldChange", "cds_frame", "ptc_position",
            "nmd_prediction", "protein_domain", "localization_consequence"
        ]],
        use_container_width=True,
        hide_index=True,
        height=360
    )

    st.markdown("---")
    st.subheader("⚠️ Predicted Functional Impairment & Loss-of-Function (LoF) Risk Ranking")
    st.caption("Quantitative Loss-of-Function (LoF) Risk Score (0-100%) based on Frame-Shift, NMD Degradation, Catalytic Domain Loss & ΔPSI Magnitude.")

    high_risk_df = isoform_annot_df[isoform_annot_df["impairment_score_pct"] >= 75.0]
    if not high_risk_df.empty:
        st.error(f"🚨 High Functional Impairment Risk Detected in {len(high_risk_df)} target genes! (NMD Degradation / Critical Domain Loss)")
        for _, r in high_risk_df.head(4).iterrows():
            st.markdown(f"• **{r['geneSymbol']}** (`{r['event_type']}`): **{r['impairment_tier']}** — {r['primary_dysfunction_cause']} (Transcript: `{r['transcript_id']}`)")

    isoform_csv_bytes = isoform_annot_df.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="📥 Download Event-Level Isoform & Impairment Table (.csv)",
        data=isoform_csv_bytes,
        file_name="gensplice_event_level_isoform_annotations.csv",
        mime="text/csv",
        key="dl_isoform_csv_btn"
    )

# ==============================================================================
# SECTION 1: 📚 NCBI / PubMed Automated Gene & Literature Explorer (100% Free NIH API)
# ==============================================================================
st.markdown("---")
st.subheader("📚 NCBI / PubMed Automated Gene & Literature Explorer")
st.caption("100% Free NIH E-utilities API Integration • Official NCBI Descriptions & Clickable PubMed Papers")

target_gene_options = q1_current_genes if q1_current_genes else df_merged.select("geneSymbol").to_series().to_list()
selected_ncbi_gene = st.selectbox(
    "🔍 Select Q1 Target Gene to Inspect NCBI Details & PubMed Papers:",
    options=target_gene_options,
    key="ncbi_gene_selector_main"
)

if selected_ncbi_gene:
    ncbi_info = fetch_ncbi_gene_summary(selected_ncbi_gene)
    pubmed_papers = fetch_pubmed_literature(selected_ncbi_gene, top_n=3)

    g_col1, g_col2 = st.columns([1, 1])
    with g_col1:
        st.markdown(f"""
        <div style="background-color: #FFFFFF; border: 1px solid #CBD5E1; border-radius: 12px; padding: 20px; height: 100%; box-shadow: 0 2px 6px rgba(0,0,0,0.04);">
            <h4 style="margin-top: 0; color: #0F172A; border-bottom: 2px solid #E2E8F0; padding-bottom: 8px;">🧬 NCBI Official Gene Summary: <code style="color:#A855F7;">{ncbi_info['gene_symbol']}</code></h4>
            <p style="margin-bottom: 6px;"><b>Official Name:</b> {ncbi_info['official_name']}</p>
            <p style="margin-bottom: 6px;"><b>NCBI Gene ID:</b> <code>{ncbi_info['ncbi_id']}</code> | <b>Map Location:</b> {ncbi_info['chromosome']}</p>
            <p style="font-size: 0.92rem; color: #475569; line-height: 1.55; margin-top: 10px;">{ncbi_info['summary']}</p>
        </div>
        """, unsafe_allow_html=True)

    with g_col2:
        st.markdown(f"#### 📖 Recent PubMed Literature for `{selected_ncbi_gene}`")
        for paper in pubmed_papers:
            st.markdown(f"""
            <div style="background-color: #F8FAFC; border-left: 4px solid #3B82F6; border-radius: 8px; padding: 12px 16px; margin-bottom: 10px; border: 1px solid #E2E8F0; border-left: 4px solid #3B82F6;">
                <a href="{paper['url']}" target="_blank" style="text-decoration: none; font-weight: 700; color: #1D4ED8; font-size: 0.95rem;">🔗 {paper['title']}</a><br>
                <div style="font-size: 0.82rem; color: #64748B; margin-top: 4px;">
                    <b>Journal:</b> {paper['journal']} ({paper['pub_date']}) | <b>PMID:</b> <code style="color:#2563EB;">{paper['pmid']}</code>
                </div>
            </div>
            """, unsafe_allow_html=True)
