"""
GenSplice-Agent Gemini AI Biological Evaluator Engine
Uses google-genai SDK for BYOK mechanism interpretation.
"""

import os

try:
    from google import genai
    GEMINI_SDK_AVAILABLE = True
except ImportError:
    GEMINI_SDK_AVAILABLE = False

SYSTEM_INSTRUCTION = """
You are a senior bioinformatics research scientist specializing in transcriptomics, RNA alternative splicing, and molecular genetics.
Your job is to analyze quantitative gene expression data (DEG, Log2FC) alongside qualitative splicing data (rMATS 5 event types, ΔPSI).
Focus especially on 'Splicing-Driven' genes (Q2: total expression remains unchanged, but splicing isoform structure changes drastically).
Deliver rigorous biological insights, domain loss/NMD hypotheses, and experimental qRT-PCR primer validation strategies.
Write your analysis clearly in GitHub-flavored Markdown.
"""

PROMPT_TEMPLATE = """
### Target Gene Analysis Request

**Gene Information:**
- Gene Symbol: {gene_symbol}
- Gene ID: {gene_id}
- Quadrant Classification: {quadrant} ({quadrant_desc})
- Expression Change (Log₂FC): {log2fc:.3f} (DEG FDR: {deg_fdr:.2e})
- Splicing Difference (ΔPSI): {delta_psi:.3f} (rMATS FDR: {as_fdr:.2e})
- Splicing Event Type: {event_type}
- Exon Coordinates: {coordinates}

**Analytical Directives:**
1. **Quantity vs Quality Evaluation**: Explain why evaluating total expression (Log₂FC = {log2fc:.3f}) alone fails to capture the physiological impact of this gene, and how the splicing change (ΔPSI = {delta_psi:.3f}, {event_type}) alters functional isoform stoichiometry.
2. **Molecular Mechanism Hypotheses**: Propose specific molecular mechanism hypotheses (e.g., protein domain loss, nonsense-mediated mRNA decay (NMD), nuclear localization signal disruption, or isoform switching).
3. **Experimental Validation Strategy**: Outline a step-by-step experimental validation protocol using isoform-specific qRT-PCR, detailing specific primer design strategies targeting exon-exon junctions (included vs. skipped isoforms).
"""

def evaluate_gene_with_gemini(
    gene_info: dict,
    api_key: str = None,
    model_name: str = "gemini-2.5-flash"
) -> str:
    """
    Calls Google Gemini AI API to evaluate a target gene's splicing mechanism.
    """
    effective_api_key = api_key or os.environ.get("GEMINI_API_KEY")

    if not effective_api_key:
        return (
            "⚠️ **Gemini API Key Required**\n\n"
            "Please enter your Google Gemini API Key in the sidebar or set the `GEMINI_API_KEY` environment variable "
            "to activate automated biological mechanism evaluation."
        )

    if not GEMINI_SDK_AVAILABLE:
        return (
            "⚠️ **google-genai SDK Not Installed**\n\n"
            "Please install the SDK via `pip install google-genai` to enable live AI evaluations."
        )

    quadrant_desc = {
        "Q1": "Dual Responders (Expression & Splicing Both Changed)",
        "Q2": "Splicing-Driven / Masked (Expression Unchanged, Splicing Changed - Primary Target)",
        "Q3": "Invariant Background",
        "Q4": "Abundance-Driven (Expression Changed Only)"
    }.get(gene_info.get("quadrant", "Q2"), "Target Gene")

    prompt = PROMPT_TEMPLATE.format(
        gene_symbol=gene_info.get("geneSymbol", "N/A"),
        gene_id=gene_info.get("gene_id", "N/A"),
        quadrant=gene_info.get("quadrant", "Q2"),
        quadrant_desc=quadrant_desc,
        log2fc=float(gene_info.get("log2FoldChange", 0.0)),
        deg_fdr=float(gene_info.get("deg_fdr", 1.0)),
        delta_psi=float(gene_info.get("delta_psi", 0.0)),
        as_fdr=float(gene_info.get("as_fdr", 1.0)),
        event_type=gene_info.get("event_type", "SE"),
        coordinates=gene_info.get("coordinates", "N/A")
    )

    try:
        client = genai.Client(api_key=effective_api_key)
        response = client.models.generate_content(
            model=model_name,
            contents=[SYSTEM_INSTRUCTION, prompt]
        )
        return response.text
    except Exception as e:
        return f"❌ **Gemini API Error:** {str(e)}"
