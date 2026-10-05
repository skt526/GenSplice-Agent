"""
GenSplice-Agent Configuration & Global Design System
Light Mode Theme with Shaded Translucent Threshold Regions
"""

# Default Statistical Cutoffs
DEFAULT_LOG2FC_CUTOFF = 1.0
DEFAULT_DELTA_PSI_CUTOFF = 0.1
DEFAULT_DEG_FDR_CUTOFF = 0.05
DEFAULT_AS_FDR_CUTOFF = 0.05
DEFAULT_MIN_JUNCTION_READS = 10  # Minimum total junction reads across replicates to filter low-coverage false positives

# Background Noise Exclusion Cutoffs (Excludes unperturbed genes/events with zero biological change)
DEFAULT_FILTER_NOISE = True
NOISE_DELTA_PSI_CUTOFF = 0.05  # |dPSI| <= 0.05 represents unperturbed / subtle technical noise
NOISE_LOG2FC_CUTOFF = 0.01     # |log2FC| <= 0.01 represents unperturbed / zero expression shift
NOISE_FDR_CUTOFF = 0.90        # FDR >= 0.90 indicates non-significant / zero-confidence background noise

# Premium Light Mode Theme Colors
THEME_BG = "#F8FAFC"
THEME_CARD_BG = "#FFFFFF"
THEME_TEXT = "#0F172A"
THEME_MUTED = "#64748B"
THEME_BORDER = "#E2E8F0"

# Threshold Line Colors
AS_THRESHOLD_COLOR = "#3B82F6"   # Blue horizontal lines for Splicing

# 5 rMATS Splicing Event Type Colors
EVENT_COLORS = {
    "SE": "#EF4444",    # Exon Skipping (Red)
    "RI": "#10B981",    # Retained Intron (Green)
    "MXE": "#F59E0B",   # Mutually Exclusive Exons (Amber)
    "A5SS": "#A855F7",  # Alt 5' Splice Site (Purple)
    "A3SS": "#06B6D4",  # Alt 3' Splice Site (Cyan)
    "None": "#94A3B8"   # No Splicing Event (Gray)
}

# Alternative Splicing Volcano Plot Theme Colors
VOLCANO_COLORS = {
    "Inclusion": "#2563EB",        # Royal Blue (Inclusion Favored, ΔPSI > 0)
    "Exclusion": "#E11D48",        # Vibrant Crimson / Rose Red (Exclusion Favored, ΔPSI < 0)
    "Non-Significant": "#94A3B8"   # Slate Gray (Invariant / Background)
}

VOLCANO_SHADING = {
    "Inclusion": "rgba(37, 99, 235, 0.08)",    # Soft Translucent Blue
    "Exclusion": "rgba(225, 29, 72, 0.08)",     # Soft Translucent Rose/Red
    "Background": "rgba(241, 245, 249, 0.50)"  # Neutral Translucent Slate
}
