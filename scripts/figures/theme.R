# =============================================================================
# thesis_theme.R — Shared ggplot2 style for all thesis figures
#
# Colors: Wong (2011) / Okabe-Ito colorblind-safe palette
#
# Usage:
#   source(file.path(script_dir, "thesis_theme.R"))
# =============================================================================

# ---------------------------------------------------------------------------
# Font sizes
# ---------------------------------------------------------------------------
THESIS_BASE_SIZE  <- 16   # ggplot2 base_size (scales axes, legend, etc.)
THESIS_TITLE_SIZE <- 20   # panel labels A / B (plot.title)
THESIS_ANNOT_SIZE <- 4.5  # geom_text / annotate size (ggplot mm ≈ pt/2.845)

# ---------------------------------------------------------------------------
# Colors — Wong (2011) / Okabe-Ito colorblind-safe palette
# Matches thesis_style.py exactly.
# ---------------------------------------------------------------------------
CB_TEAL      <- "#009E73"   # teal-green   → primary accent / PGLS lines
CB_BLUE      <- "#0072B2"   # deep blue    → OLS lines / second series
CB_ORANGE    <- "#E69F00"   # amber-orange → third series / bar fill
CB_VERMILLON <- "#D55E00"   # vermillion   → outliers / warning
CB_SKYBLUE   <- "#56B4E9"   # sky blue     → neutral bar fill
CB_YELLOW    <- "#F0E442"   # yellow       → fourth series (use sparingly)
CB_PINK      <- "#CC79A7"   # pink-mauve   → fifth series
CB_BLACK     <- "#000000"   # black        → data points, axes

# Aliases
ACCENT_COLOR  <- CB_TEAL
COLOR_OLS     <- CB_BLUE
COLOR_PGLS    <- CB_TEAL
COLOR_OUTLIER <- CB_VERMILLON
COLOR_BAR     <- CB_SKYBLUE

# Full ordered palette vector
THESIS_PALETTE <- c(CB_TEAL, CB_BLUE, CB_ORANGE, CB_VERMILLON,
                    CB_SKYBLUE, CB_YELLOW, CB_PINK, CB_BLACK)

# ---------------------------------------------------------------------------
# Base theme (for scatter / line / bar plots)
# ---------------------------------------------------------------------------
theme_thesis <- function(base_size = THESIS_BASE_SIZE, ...) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
  ggplot2::theme(
    axis.line         = ggplot2::element_line(linewidth = 0.8, color = "black"),
    axis.ticks        = ggplot2::element_line(linewidth = 0.6, color = "black"),
    axis.ticks.length = ggplot2::unit(0.12, "cm"),
    axis.text         = ggplot2::element_text(color = "black", face = "plain"),
    axis.title        = ggplot2::element_text(face = "plain"),
    axis.title.x      = ggplot2::element_text(face = "plain", margin = ggplot2::margin(t = 10)),
    panel.grid.major  = ggplot2::element_line(color = "gray92", linewidth = 0.3),
    panel.grid.minor  = ggplot2::element_blank(),
    legend.background = ggplot2::element_rect(fill = "white", color = "gray80",
                                              linewidth = 0.3),
    legend.title      = ggplot2::element_text(face = "bold"),
    legend.key.size   = ggplot2::unit(0.45, "cm"),
    plot.title        = ggplot2::element_text(face = "bold",
                                              size = THESIS_TITLE_SIZE),
    ...
  )
}

# Variant without y-axis text (for multi-panel rows sharing a y label)
theme_thesis_noy <- function(base_size = THESIS_BASE_SIZE, ...) {
  theme_thesis(base_size = base_size, ...) +
  ggplot2::theme(
    axis.text.y  = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank()
  )
}

# Compact theme for dense plots (heatmaps) — readable but space-efficient
theme_thesis_heatmap <- function(base_size = 12, ...) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "sans") +
  ggplot2::theme(
    panel.grid      = ggplot2::element_blank(),
    panel.border    = ggplot2::element_rect(color = "black", fill = NA,
                                            linewidth = 0.5),
    legend.title    = ggplot2::element_text(face = "bold",
                                            size = base_size * 0.9),
    legend.text     = ggplot2::element_text(size = base_size * 0.85),
    legend.key.height = ggplot2::unit(1.2, "cm"),
    legend.key.width  = ggplot2::unit(0.3, "cm"),
    axis.ticks      = ggplot2::element_line(linewidth = 0.3),
    plot.margin     = ggplot2::margin(5, 5, 70, 5),
    ...
  )
}
