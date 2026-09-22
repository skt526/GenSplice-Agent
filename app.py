"""
GenSplice-Agent Streamlit Dashboard & AI Agent
Integrative Quantity (DEG) & Quality (Alternative Splicing 5 Event Types) Platform
"""

import os
import streamlit as st
import polars as pl
import pandas as pd

from config import (
    DEFAULT_LOG2FC_CUTOFF, DEFAULT_DELTA_PSI_CUTOFF,
    DEFAULT_DEG_FDR_CUTOFF, DEFAULT_AS_FDR_CUTOFF
)
from core.deg_loader import load_deg_data
from core.rmats_loader import load_rmats_data, select_primary_splicing_events
from core.merger import merge_deg_and_rmats, get_quadrant_kpis
from visualizer.quadrant_plot import build_quadrant_plot
from visualizer.dual_volcano import build_dual_volcano_plot
from visualizer.report_exporter import export_html_report
from ai.gemini_evaluator import evaluate_gene_with_gemini

# Streamlit Page Config
st.set_page_config(
    page_title="GenSplice-Agent Transcriptomics Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for Premium Design
st.markdown("""
<style>
    .main-header {
        font-size: 2.2rem;
        font-weight: 800;
        background: linear-gradient(135deg, #4A0E4E 0%, #8E44AD 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        margin-bottom: 0.5rem;
    }
    .kpi-box {
        background-color: #FFFFFF;
        border-radius: 12px;
        padding: 16px;
        border: 1px solid #E2E8F0;
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
        text-align: center;
    }
    .kpi-val {
        font-size: 2rem;
        font-weight: 800;
        margin-top: 4px;
    }
    .kpi-lbl {
        font-size: 0.85rem;
        color: #64748B;
        font-weight: 600;
        text-transform: uppercase;
    }
    .stButton>button {
        border-radius: 8px;
        font-weight: 600;
    }
    .threshold-container {
        background-color: #F8FAFC;
        border-radius: 12px;
        padding: 16px 24px;
        border: 2px solid #8E44AD;
        margin-bottom: 16px;
    }
</style>
""", unsafe_allow_html=True)

# Sidebar Configuration
st.sidebar.image("https://img.icons8.com/color/96/dna-helix.png", width=64)
st.sidebar.title("GenSplice-Agent")
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

# 2. Sidebar Quick Controls
st.sidebar.header("⚙️ Global Settings")
color_option = st.sidebar.radio("Plot Color Palette", ["By Quadrant", "By Splicing Event Type"])
color_by = "quadrant" if color_option == "By Quadrant" else "event_type"

st.sidebar.markdown("---")

# 3. Gemini AI Key
st.sidebar.header("🔑 Google Gemini AI Key")
api_key = st.sidebar.text_input("Gemini API Key (BYOK)", type="password", help="Enter your Gemini API key for live mechanism analysis")

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
st.markdown("Integrative Quantitative (DEG) & Qualitative (Alternative Splicing) Transcriptomics Dashboard")

if df_deg_raw.height == 0 and df_rmats_raw.height == 0:
    st.warning("⚠️ No valid DEG or rMATS data found. Please run `./GenSplice` or `./test` to generate output data.")
    st.stop()

# Real-Time Interactive Threshold Adjustment Panel directly at top
st.markdown('<div class="threshold-container">', unsafe_allow_html=True)
st.subheader("🎛️ Real-Time Interactive Threshold Adjustment")
tcol1, tcol2, tcol3, tcol4 = st.columns(4)

with tcol1:
    log2fc_cutoff = st.slider("Log₂FC Cutoff (|Log₂FC|)", 0.1, 3.0, DEFAULT_LOG2FC_CUTOFF, 0.05, key="log2fc_slider")
with tcol2:
    delta_psi_cutoff = st.slider("ΔPSI Cutoff (|ΔPSI|)", 0.01, 0.5, DEFAULT_DELTA_PSI_CUTOFF, 0.01, key="psi_slider")
with tcol3:
    deg_fdr_cutoff = st.select_slider("DEG FDR Cutoff", options=[0.001, 0.01, 0.05, 0.1], value=DEFAULT_DEG_FDR_CUTOFF, key="deg_fdr_slider")
with tcol4:
    as_fdr_cutoff = st.select_slider("rMATS FDR Cutoff", options=[0.001, 0.01, 0.05, 0.1], value=DEFAULT_AS_FDR_CUTOFF, key="as_fdr_slider")
st.markdown('</div>', unsafe_allow_html=True)

# Merge dataset dynamically based on live slider values
df_merged = merge_deg_and_rmats(
    df_deg=df_deg_raw,
    df_rmats=df_rmats_raw,
    log2fc_cutoff=log2fc_cutoff,
    delta_psi_cutoff=delta_psi_cutoff,
    deg_fdr_cutoff=deg_fdr_cutoff,
    as_fdr_cutoff=as_fdr_cutoff
)

# Summary KPI Cards
kpis = get_quadrant_kpis(df_merged)
col1, col2, col3, col4, col5 = st.columns(5)

with col1:
    st.markdown(f'<div class="kpi-box"><div class="kpi-lbl">Total Analyzed</div><div class="kpi-val" style="color:#1E293B;">{kpis["total"]}</div></div>', unsafe_allow_html=True)
with col2:
    st.markdown(f'<div class="kpi-box" style="border-top:4px solid #8E44AD;"><div class="kpi-lbl" style="color:#8E44AD;">Q2: Splicing Target</div><div class="kpi-val" style="color:#8E44AD;">{kpis["Q2"]}</div></div>', unsafe_allow_html=True)
with col3:
    st.markdown(f'<div class="kpi-box" style="border-top:4px solid #E63946;"><div class="kpi-lbl" style="color:#E63946;">Q1: Dual Responders</div><div class="kpi-val" style="color:#E63946;">{kpis["Q1"]}</div></div>', unsafe_allow_html=True)
with col4:
    st.markdown(f'<div class="kpi-box" style="border-top:4px solid #2980B9;"><div class="kpi-lbl" style="color:#2980B9;">Q4: DEG Only</div><div class="kpi-val" style="color:#2980B9;">{kpis["Q4"]}</div></div>', unsafe_allow_html=True)
with col5:
    st.markdown(f'<div class="kpi-box" style="border-top:4px solid #95A5A6;"><div class="kpi-lbl" style="color:#95A5A6;">Q3: Invariant</div><div class="kpi-val" style="color:#95A5A6;">{kpis["Q3"]}</div></div>', unsafe_allow_html=True)

st.markdown("<br>", unsafe_allow_html=True)

# Export HTML Report Button in Top Action Bar
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
        st.success("✔ Report generated! Openable in Chrome:")
        with open(html_out_path, "r", encoding="utf-8") as f:
            html_bytes = f.read().encode("utf-8")
        st.download_button(
            label="💾 Download Chrome HTML Report",
            data=html_bytes,
            file_name="gensplice_report.html",
            mime="text/html"
        )

# Main Navigation Tabs
tab1, tab2, tab3, tab4 = st.tabs([
    "📊 4-Quadrant Cross-Plot",
    "🌋 Dual Volcano View",
    "📋 Data Explorer & Target Selector",
    "🤖 Gemini AI Biological Analyst"
])

# Tab 1: Quadrant Plot
with tab1:
    st.subheader("Interactive 4-Quadrant Cross-Plot (Real-Time Drag & Update)")
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

# Tab 2: Dual Volcano
with tab2:
    st.subheader("Parallel Dual Volcano View (DEG vs rMATS Splicing)")
    fig_volc = build_dual_volcano_plot(df_merged, log2fc_cutoff, delta_psi_cutoff, deg_fdr_cutoff, as_fdr_cutoff)
    st.plotly_chart(fig_volc, use_container_width=True)

# Tab 3: Data Explorer
with tab3:
    st.subheader("Transcriptomics Data Explorer")
    quad_filter = st.multiselect("Filter by Quadrant", options=["Q1", "Q2", "Q3", "Q4"], default=["Q2", "Q1"])
    
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

# Tab 4: Gemini AI Agent
with tab4:
    st.subheader("🤖 Google Gemini AI Biological Evaluator")
    st.markdown("Select a target gene (especially **Q2 Splicing-Driven** genes) to analyze functional domain loss, NMD, and qRT-PCR validation primers.")
    
    target_genes = df_merged.filter(pl.col("quadrant").is_in(["Q2", "Q1"])).select("geneSymbol").to_series().to_list()
    if not target_genes:
        target_genes = df_merged.select("geneSymbol").to_series().to_list()
        
    selected_gene_symbol = st.selectbox("Select Target Gene for AI Evaluation", options=target_genes)
    
    if selected_gene_symbol:
        sub_df = df_merged.filter(pl.col("geneSymbol") == selected_gene_symbol)
        if sub_df.height > 0:
            gene_dict = sub_df.to_pandas().iloc[0].to_dict()
            
            st.info(f"**Target Gene:** `{gene_dict.get('geneSymbol')}` | **Quadrant:** `{gene_dict.get('quadrant')}` | **ΔPSI:** `{gene_dict.get('delta_psi'):.3f}` | **Log₂FC:** `{gene_dict.get('log2FoldChange'):.3f}`")
            
            if st.button("🚀 Run Gemini AI Analysis"):
                with st.spinner(f"Evaluating biological mechanism for {selected_gene_symbol} via Gemini..."):
                    ai_result = evaluate_gene_with_gemini(gene_dict, api_key=api_key)
                    st.markdown(ai_result)
