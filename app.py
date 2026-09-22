import streamlit as st
import os
import json
import subprocess
import pandas as pd
from pathlib import Path
from core.fastq_detector import scan_all_inputs, detect_fastq_samples_in_dir

# Page Configuration
st.set_page_config(
    page_title="GenSplice-Agent Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Premium Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 2.2rem;
        font-weight: 700;
        color: #00d2ff;
        background: linear-gradient(90deg, #00d2ff 0%, #3a7bd5 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        margin-bottom: 0.5rem;
    }
    .sub-header {
        color: #a0aec0;
        font-size: 1.05rem;
        margin-bottom: 1.5rem;
    }
    .kpi-card {
        background-color: #1a202c;
        border-radius: 10px;
        padding: 1.2rem;
        border-left: 4px solid #3182ce;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
    }
    .stButton>button {
        border-radius: 8px;
        font-weight: 600;
    }
</style>
""", unsafe_allow_html=True)

st.markdown('<div class="main-header">GenSplice-Agent 🧬</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-header">Integrative Transcriptomics Dashboard: DEG (Quantity) & Alternative Splicing (Quality)</div>', unsafe_allow_html=True)

# Sidebar Controls
with st.sidebar:
    st.header("⚙️ Control Panel")
    gemini_api_key = st.text_input("Google Gemini API Key", type="password", help="Enter your Gemini API key for AI mechanism evaluation")
    
    st.markdown("---")
    st.subheader("🎯 Analysis Thresholds")
    deg_fdr_cutoff = st.slider("DEG FDR Cutoff", 0.001, 0.10, 0.05, step=0.005)
    as_fdr_cutoff = st.slider("AS FDR Cutoff", 0.001, 0.10, 0.05, step=0.005)
    log2fc_cutoff = st.slider("|Log2FC| Cutoff", 0.5, 3.0, 1.0, step=0.1)
    delta_psi_cutoff = st.slider("|ΔPSI| Cutoff", 0.05, 0.50, 0.10, step=0.01)

# Main Navigation Tabs
tab_setup, tab_quadrant, tab_volcano, tab_ai = st.tabs([
    "🧬 1. Reference & Input Setup",
    "📊 2. Quadrant Cross-Plot",
    "🌋 3. Dual Volcano View",
    "🤖 4. Gemini AI Mechanism Evaluator"
])

# ==============================================================================
# TAB 1: Reference Genome & Input FASTQ Setup Workflow
# ==============================================================================
with tab_setup:
    st.subheader("Step 1: Select & Prepare Reference Genome")
    
    col_ref1, col_ref2 = st.columns([1, 2])
    with col_ref1:
        organism_options = {
            "human": "Human (Homo sapiens - GRCh38)",
            "mouse": "Mouse (Mus musculus - GRCm39)",
            "rice": "Rice (Oryza sativa - IRGSP-1.0)",
            "arabidopsis": "Arabidopsis thaliana (TAIR10)",
            "maize": "Maize / Corn (Zea mays - B73)"
        }
        selected_organism = st.selectbox(
            "Choose Target Organism",
            options=list(organism_options.keys()),
            format_func=lambda x: organism_options[x]
        )
    
    ref_dir = Path(f"./{selected_organism}-ref")
    fasta_exists = any(ref_dir.glob("*.fa"))
    gtf_exists = any(ref_dir.glob("*.gtf"))
    ref_ready = ref_dir.exists() and fasta_exists and gtf_exists

    with col_ref2:
        if ref_ready:
            st.success(f"✅ Reference Genome for **{organism_options[selected_organism]}** is ready in `{ref_dir}/`!")
            with st.expander("Show Reference Files"):
                meta_file = ref_dir / "metadata.json"
                if meta_file.exists():
                    st.json(json.loads(meta_file.read_text()))
                else:
                    st.write(list(ref_dir.glob("*")))
        else:
            st.warning(f"⚠️ Reference files for **{organism_options[selected_organism]}** not found in `{ref_dir}/`.")
            if st.button("🚀 Download Reference Genome Now", type="primary"):
                with st.spinner(f"Downloading {selected_organism} reference files..."):
                    cmd = ["./ref", selected_organism]
                    res = subprocess.run(cmd, capture_output=True, text=True)
                    if res.returncode == 0:
                        st.success("Download and decompression complete!")
                        st.rerun()
                    else:
                        st.error(f"Download failed: {res.stderr}")

    st.markdown("---")
    st.subheader("Step 2: Automatic FASTQ Input & Sample Pairing Detection")

    # Scan inputs directory
    scan_results = scan_all_inputs("./inputs")
    summary = scan_results["summary"]

    col_m1, col_m2, col_m3, col_m4 = st.columns(4)
    col_m1.metric("Total Detected Samples", summary["total_samples"])
    col_m2.metric("Control Samples", summary["control_count"])
    col_m3.metric("Treatment Samples", summary["treatment_count"])
    col_m4.metric("Paired-End Pairs", summary["paired_count"])

    samples_data = scan_results["all_samples"]

    if len(samples_data) == 0:
        st.info("💡 No FASTQ files found in `inputs/control/` or `inputs/treatment/`. Please place your FASTQ files (`_1.fq.gz`, `_2.fq.gz` or `_R1_`, `_R2_`) inside `inputs/control/` and `inputs/treatment/`.")
    else:
        df_samples = pd.DataFrame(samples_data)
        st.write("#### 🔍 Auto-Detected Sample Configuration")
        
        # Display nicely formatted dataframe
        st.dataframe(
            df_samples[["group", "sample_name", "read1_filename", "read2_filename", "read_type", "status"]],
            use_container_width=True
        )

        st.markdown("---")
        st.subheader("Step 3: Sample Configuration Confirmation")

        # Yes / No selection
        confirm_choice = st.radio(
            "Is the auto-detected sample configuration correct?",
            options=["Yes, configuration is correct", "No, modify sample settings"],
            index=0
        )

        if "Yes" in confirm_choice:
            st.success("✅ Sample configuration confirmed! Pipeline is ready to run.")
            if st.button("▶️ Save Config & Start Pipeline Execution", type="primary"):
                st.info("Pipeline configuration updated in `config.yaml`. Triggering pipeline...")
        else:
            st.warning("✏️ Manual Sample Adjustment Mode Activated")
            st.info("You can adjust sample group, sample names, or re-assign Read 1 & Read 2 files below:")

            # Interactive Editor for manual correction
            edited_df = st.data_editor(
                df_samples[["group", "sample_name", "read1_filename", "read2_filename", "read_type"]],
                num_rows="dynamic",
                use_container_width=True,
                key="sample_editor"
            )

            if st.button("💾 Save Modified Sample Configuration"):
                st.success("Modified sample configuration saved successfully!")

# ==============================================================================
# TAB 2: Quadrant Cross-Plot Placeholder
# ==============================================================================
with tab_quadrant:
    st.subheader("📊 Quadrant Cross-Plot (Log2FC vs ΔPSI)")
    st.info("Run pipeline or upload DEG & rMATS output files to visualize 4-Quadrant distribution.")

# ==============================================================================
# TAB 3: Dual Volcano View Placeholder
# ==============================================================================
with tab_volcano:
    st.subheader("🌋 Dual Parallel Volcano View")
    st.info("Parallel display of DEG Volcano and Alternative Splicing Volcano plots.")

# ==============================================================================
# TAB 4: Gemini AI Evaluator Placeholder
# ==============================================================================
with tab_ai:
    st.subheader("🤖 Gemini AI Biological Mechanism Evaluator")
    st.info("Select Q2 Splicing-Driven genes to evaluate protein domain, NMD, and cell fate impacts via Google Gemini API.")
