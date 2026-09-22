"""
GenSplice-Agent Configuration & Global Design System
Light Mode Theme with Shaded Translucent Threshold Regions
"""

import os

# Default Statistical Cutoffs
DEFAULT_LOG2FC_CUTOFF = 1.0
DEFAULT_DELTA_PSI_CUTOFF = 0.1
DEFAULT_DEG_FDR_CUTOFF = 0.05
DEFAULT_AS_FDR_CUTOFF = 0.05

# Premium Light Mode Theme Colors
THEME_BG = "#F8FAFC"
THEME_CARD_BG = "#FFFFFF"
THEME_TEXT = "#0F172A"
THEME_MUTED = "#64748B"
THEME_BORDER = "#E2E8F0"

# Quadrant Palette (DEG = Red, AS = Blue, Both = Purple, Invariant = Slate Gray)
QUADRANT_COLORS = {
    "Q1": "#A855F7",  # Purple (Both DEG & Alternative Splicing)
    "Q2": "#3B82F6",  # Blue (Alternative Splicing Only - Primary Target)
    "Q3": "#94A3B8",  # Slate Gray (Invariant Background)
    "Q4": "#EF4444",  # Red (DEG Expression Only)
}

# Translucent Region Shading Colors (Layered below scatter markers)
REGION_FILL_COLORS = {
    "Q1": "rgba(168, 85, 247, 0.12)", # Translucent Purple for Both
    "Q2": "rgba(59, 130, 246, 0.12)",  # Translucent Blue for Splicing Only
    "Q3": "rgba(241, 245, 249, 0.60)", # Soft Gray for Invariant
    "Q4": "rgba(239, 68, 68, 0.12)"   # Translucent Red for DEG Only
}

QUADRANT_LABELS = {
    "Q1": "Q1: Dual Responders (Both DEG & Splicing)",
    "Q2": "Q2: Splicing-Driven (Alternative Splicing Only)",
    "Q3": "Q3: Invariant Background",
    "Q4": "Q4: Expression-Driven (DEG Only)",
}

# Threshold Line Colors
DEG_THRESHOLD_COLOR = "#EF4444"  # Red vertical lines for DEG
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
