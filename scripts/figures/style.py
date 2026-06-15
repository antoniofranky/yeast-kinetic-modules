"""
Shared matplotlib style for all thesis figures.

Usage (from a script in a subdirectory):
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[N] / "results" / "scripts"))
    import thesis_style
    thesis_style.apply()
"""

import matplotlib as mpl

# ---------------------------------------------------------------------------
# Font sizes
# ---------------------------------------------------------------------------
FS_LABEL  = 20   # axis labels
FS_TICK   = 16   # tick labels
FS_LEGEND = 15   # legend entries
FS_LTITLE = 16   # legend title
FS_PANEL  = 22   # A / B panel labels
FS_ANNOT  = 13   # in-plot annotations

# ---------------------------------------------------------------------------
# Colors — Wong (2011) / Okabe-Ito colorblind-safe palette
# ---------------------------------------------------------------------------
# Use these for categorical / accent colors across all figures.
# Order chosen so the first few are maximally distinct.
CB_TEAL      = "#009E73"   # teal-green   → primary accent, PGLS lines
CB_BLUE      = "#0072B2"   # deep blue    → OLS lines, second series
CB_ORANGE    = "#E69F00"   # amber-orange → third series / bar fill
CB_VERMILLON = "#D55E00"   # vermillion   → outliers / warning
CB_SKYBLUE   = "#56B4E9"   # sky blue     → neutral bar fill
CB_YELLOW    = "#F0E442"   # yellow       → fourth series (use sparingly)
CB_PINK      = "#CC79A7"   # pink-mauve   → fifth series
CB_BLACK     = "#000000"   # black        → data points, axes

# Convenient aliases used in scripts
ACCENT        = CB_TEAL
COLOR_A       = CB_BLUE
COLOR_B       = CB_ORANGE
COLOR_C       = CB_VERMILLON
COLOR_D       = CB_SKYBLUE
COLOR_OLS     = CB_BLUE
COLOR_PGLS    = CB_TEAL
COLOR_OUTLIER = CB_VERMILLON
COLOR_BAR     = CB_SKYBLUE

# Full ordered palette list (drop-in for sns.color_palette)
PALETTE = [CB_TEAL, CB_BLUE, CB_ORANGE, CB_VERMILLON,
           CB_SKYBLUE, CB_YELLOW, CB_PINK, CB_BLACK]

# ---------------------------------------------------------------------------
# Apply to matplotlib rcParams
# ---------------------------------------------------------------------------
def apply():
    """Call once at the top of a plotting script."""
    mpl.rcParams.update({
        "font.family":            "sans-serif",
        "font.size":              FS_TICK,
        "axes.titlesize":         FS_PANEL,
        "axes.labelsize":         FS_LABEL,
        "xtick.labelsize":        FS_TICK,
        "ytick.labelsize":        FS_TICK,
        "legend.fontsize":        FS_LEGEND,
        "legend.title_fontsize":  FS_LTITLE,
        "legend.framealpha":      0.85,
        "legend.edgecolor":       "0.7",
        "axes.grid":              True,
        "grid.alpha":             0.3,
        "grid.color":             "#d0d0d0",
        "axes.spines.top":        False,
        "axes.spines.right":      False,
        "figure.dpi":             150,
        "savefig.dpi":            300,
        "savefig.bbox":           "tight",
    })
