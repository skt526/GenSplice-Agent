"""
GenSplice-Agent Configuration & Global Design System
"""

import os

# Default Statistical Cutoffs
DEFAULT_LOG2FC_CUTOFF = 1.0
DEFAULT_DELTA_PSI_CUTOFF = 0.1
DEFAULT_DEG_FDR_CUTOFF = 0.05
DEFAULT_AS_FDR_CUTOFF = 0.05

# Quadrant Palette (Harmonious Modern Theme)
QUADRANT_COLORS = {
    "Q1": "#E63946",  # Vibrant Coral Red (Dual Responders)
    "Q2": "#8E44AD",  # Amethyst Purple (Splicing-Driven / Primary Target)
    "Q3": "#95A5A6",  # Muted Cool Gray (Background Invariant)
    "Q4": "#2980B9",  # Deep Ocean Blue (Expression-Driven)
}

QUADRANT_LABELS = {
    "Q1": "Q1: Dual Responders (DEG & AS)",
    "Q2": "Q2: Splicing-Driven (Masked / No DEG Change)",
    "Q3": "Q3: Invariant Background",
    "Q4": "Q4: Abundance-Driven (DEG Only)",
}

# 5 rMATS Splicing Event Type Colors
EVENT_COLORS = {
    "SE": "#E74C3C",    # Exon Skipping (Red)
    "RI": "#2ECC71",    # Retained Intron (Green)
    "MXE": "#F39C12",   # Mutually Exclusive Exons (Orange)
    "A5SS": "#9B59B6",  # Alt 5' Splice Site (Purple)
    "A3SS": "#1ABC9C",  # Alt 3' Splice Site (Teal)
    "None": "#BDC3C7"   # No Splicing Event (Gray)
}

EVENT_NAMES = {
    "SE": "Skipping Exon (SE)",
    "RI": "Retained Intron (RI)",
    "MXE": "Mutually Exclusive Exons (MXE)",
    "A5SS": "Alternative 5' Splice Site (A5SS)",
    "A3SS": "Alternative 3' Splice Site (A3SS)",
    "None": "No Splicing Event"
}
