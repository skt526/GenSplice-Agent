"""
GenSplice-Agent Configuration & Global Design System
Black & White Dark Mode Theme with High-Contrast Biological Color System
"""

import os

# Default Statistical Cutoffs
DEFAULT_LOG2FC_CUTOFF = 1.0
DEFAULT_DELTA_PSI_CUTOFF = 0.1
DEFAULT_DEG_FDR_CUTOFF = 0.05
DEFAULT_AS_FDR_CUTOFF = 0.05

# Black & White Dark Mode Theme Colors
THEME_BG = "#0B0F19"
THEME_CARD_BG = "#111827"
THEME_TEXT = "#F9FAFB"
THEME_MUTED = "#9CA3AF"
THEME_BORDER = "#1F2937"

# Quadrant Palette (DEG = Red, AS = Blue, Both = Purple, Invariant = Gray)
QUADRANT_COLORS = {
    "Q1": "#A855F7",  # Purple (Both DEG & Alternative Splicing)
    "Q2": "#3B82F6",  # Blue (Alternative Splicing Only - Target)
    "Q3": "#4B5563",  # Dark Muted Gray (Invariant Background)
    "Q4": "#EF4444",  # Red (DEG Expression Only)
}

QUADRANT_LABELS = {
    "Q1": "Q1: Dual Responders (Both DEG & Splicing)",
    "Q2": "Q2: Splicing-Driven (Alternative Splicing Only)",
    "Q3": "Q3: Invariant Background",
    "Q4": "Q4: Expression-Driven (DEG Only)",
}

# Threshold Line Colors (Matching Axes)
DEG_THRESHOLD_COLOR = "#EF4444"  # Red for Log2FC DEG threshold
AS_THRESHOLD_COLOR = "#3B82F6"   # Blue for Delta PSI Splicing threshold

# 5 rMATS Splicing Event Type Colors (Harmonized with Dark Theme)
EVENT_COLORS = {
    "SE": "#EF4444",    # Exon Skipping (Red)
    "RI": "#10B981",    # Retained Intron (Green)
    "MXE": "#F59E0B",   # Mutually Exclusive Exons (Amber/Orange)
    "A5SS": "#A855F7",  # Alt 5' Splice Site (Purple)
    "A3SS": "#06B6D4",  # Alt 3' Splice Site (Cyan)
    "None": "#4B5563"   # No Splicing Event (Gray)
}

EVENT_NAMES = {
    "SE": "Skipping Exon (SE)",
    "RI": "Retained Intron (RI)",
    "MXE": "Mutually Exclusive Exons (MXE)",
    "A5SS": "Alternative 5' Splice Site (A5SS)",
    "A3SS": "Alternative 3' Splice Site (A3SS)",
    "None": "No Splicing Event"
}
