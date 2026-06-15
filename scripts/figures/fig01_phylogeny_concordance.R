library(ape)
library(phytools)
library(tidyverse)
library(ggtree)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(scales)

cat("\n", rep("=", 80), "\n", sep = "")
cat("CREATING PHYLOGENETIC FIGURE: NETWORK STRUCTURE & ROBUSTNESS\n")
cat(rep("=", 80), "\n\n", sep = "")

library(here)
REPO_ROOT   <- here::here()
DATA_DIR    <- file.path(REPO_ROOT, "data")
RESULTS_DIR <- file.path(REPO_ROOT, "results")
FIGURES_DIR <- file.path(REPO_ROOT, "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)
source(file.path(REPO_ROOT, "scripts/figures/theme.R"))

# --- GLOBAL SETTINGS FOR CONSISTENCY ---
font_size <- 34  # 20" figure displayed at ~5.9" (29.5% scale) → 34*0.295 ≈ 10pt
bar_width <- 0.7
box_line_width <- 0.7
# ---------------------------------------

# Load data
df_matched <- read_csv(file.path(DATA_DIR, "annotations/cocoa_results_merged.csv"), show_col_types = FALSE)
df_kinetic <- read_csv(file.path(RESULTS_DIR, "kinetic/species_kinetic_summary.csv.gz"), show_col_types = FALSE)
df_reactions <- read_csv(file.path(RESULTS_DIR, "kinetic/giant_module_reaction_summary.csv.gz"), show_col_types = FALSE)
df_split <- read_csv(file.path(RESULTS_DIR, "kinetic/split_model_sizes.csv.gz"), show_col_types = FALSE)

# Join kinetic data
df_matched <- df_matched %>%
  left_join(df_kinetic %>% select(species_id, n_acr, largest_module_size),
            by = c("old_species_id" = "species_id")) %>%
  left_join(df_reactions %>% select(species_id, n_giant_reactions),
            by = c("old_species_id" = "species_id")) %>%
  left_join(df_split %>% select(species_id, n_total_split_reactions),
            by = c("old_species_id" = "species_id"))

df_matched <- df_matched %>%
  mutate(giant_reaction_proportion = n_giant_reactions / n_total_split_reactions)

################################################################################
## 1. PREPARE TREE
################################################################################

cat("Preparing phylogenetic tree...\n")

tree_file <- file.path(DATA_DIR, "phylogeny/clade_tree_manual.nwk")
clade_tree <- read.tree(tree_file)
clade_tree$tip.label <- gsub("_", " ", clade_tree$tip.label)
clade_tree$tip.label <- gsub("Dipodascaceae Trichomonascaceae", "Dipodascaceae/\nTrichomonascaceae", clade_tree$tip.label)

n_tips <- length(clade_tree$tip.label)

p_tree <- ggtree(clade_tree, ladderize = FALSE, linewidth = 1) +
  theme_tree2() +
  coord_cartesian(xlim = c(-0.5, 3.5), ylim = c(0.5, n_tips + 0.5), clip = "off") +
  # Suppress the "Scale for y is already present" message by setting the scale explicitly once
  scale_y_continuous(expand = expansion(mult = c(0, 0)), breaks = 1:n_tips, labels = NULL) +
  theme(
    plot.margin = margin(15, -10, 10, 10),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

tree_data <- p_tree$data
tree_tips <- tree_data %>% filter(isTip) %>% arrange(y) %>% pull(label)

################################################################################
## 2. PREPARE DATA
################################################################################

df_matched <- df_matched %>%
  mutate(`Major clade` = gsub("Dipodascaceae/Trichomonascaceae", "Dipodascaceae/\nTrichomonascaceae", `Major clade`))

df_plot <- df_matched %>%
  filter(`Major clade` %in% tree_tips) %>%
  mutate(clade_factor = factor(`Major clade`, levels = tree_tips))

################################################################################
## 3. CREATE PLOTS
################################################################################

cat("Creating plots...\n")

species_counts <- df_plot %>%
  group_by(clade_factor) %>%
  summarise(n_species = n(), .groups = "drop")

# --- THEMES ---
theme_with_labels <- theme_thesis(font_size) +
  theme(
    axis.text.y  = element_text(size = font_size, hjust = 1, margin = margin(r = 5)),
    axis.text.x  = element_text(size = font_size - 10, angle = 45, hjust = 1),
    axis.title.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(15, 5, 10, -15),
    panel.grid.major.x = element_line(color = "gray92", linewidth = 0.4),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank()
  )

theme_no_labels <- theme_thesis_noy(font_size) +
  theme(
    axis.text.x  = element_text(size = font_size - 10, angle = 45, hjust = 1),
    legend.position = "none",
    plot.margin = margin(15, 5, 10, 5),
    panel.grid.major.x = element_line(color = "gray92", linewidth = 0.4),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank()
  )

# --- PANEL B: SPECIES COUNTS ---
p_species <- ggplot(species_counts, aes(x = clade_factor, y = n_species)) +
  geom_col(fill = COLOR_BAR, width = bar_width) +
  geom_text(aes(label = n_species), hjust = -0.2, size = 5, color = "black", fontface = "bold") +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 110), expand = expansion(mult = c(0, 0.05))) +
  coord_flip(xlim = c(0.5, n_tips + 0.5), clip = "off") +
  labs(x = NULL, y = "Species count") +
  theme_with_labels

# --- PANEL C: REACTIONS ---
p_reactions <- ggplot(df_plot, aes(x = clade_factor, y = n_reactions)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, 
               linewidth = box_line_width,          # NEW: replaces size
               median.linewidth = box_line_width,   # NEW: replaces fatten
               width = bar_width, fill = "white") +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), n.breaks = 4) +
  coord_flip(xlim = c(0.5, n_tips + 0.5), clip = "off") +
  labs(x = NULL, y = "Reactions") +
  theme_no_labels

# --- PANEL: METABOLITES ---
p_metabolites <- ggplot(df_plot, aes(x = clade_factor, y = n_metabolites)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, 
               linewidth = box_line_width, 
               median.linewidth = box_line_width,
               width = bar_width, fill = "white") +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), n.breaks = 4) +
  coord_flip(xlim = c(0.5, n_tips + 0.5), clip = "off") +
  labs(x = NULL, y = "Metabolites") +
  theme_no_labels

# --- PANEL: COMPLEXES ---
p_complexes <- ggplot(df_plot, aes(x = clade_factor, y = n_complexes)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, 
               linewidth = box_line_width, 
               median.linewidth = box_line_width,
               width = bar_width, fill = "white") +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), n.breaks = 4) +
  coord_flip(xlim = c(0.5, n_tips + 0.5), clip = "off") +
  labs(x = NULL, y = "Complexes") +
  theme_no_labels

# --- PANEL: BALANCED ---
p_balanced <- ggplot(df_plot, aes(x = clade_factor, y = n_balanced_complexes)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, 
               linewidth = box_line_width, 
               median.linewidth = box_line_width,
               width = bar_width, fill = "white") +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), n.breaks = 4) +
  coord_flip(xlim = c(0.5, n_tips + 0.5), clip = "off") +
  labs(x = NULL, y = "Balanced\ncomplexes") +
  theme_with_labels

# --- PANEL: CONCORDANT ---
p_concordant <- ggplot(df_plot, aes(x = clade_factor, y = n_concordant_total)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, 
               linewidth = box_line_width, 
               median.linewidth = box_line_width,
               width = bar_width, fill = "white") +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), n.breaks = 4) +
  coord_flip(xlim = c(0.5, n_tips + 0.5), clip = "off") +
  labs(x = NULL, y = "Concordant\npairs") +
  theme_no_labels

# --- PANEL: GIANT MODULE REACTION PROPORTION ---
p_module <- ggplot(df_plot, aes(x = clade_factor, y = giant_reaction_proportion)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2,
               linewidth = box_line_width,
               median.linewidth = box_line_width,
               width = bar_width, fill = "white",
               na.rm = TRUE) +  # EXPLICITLY IGNORE MISSING DATA
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), n.breaks = 4,
                     labels = scales::percent_format(accuracy = 1)) +
  coord_flip(xlim = c(0.5, n_tips + 0.5), clip = "off") +
  labs(x = NULL, y = "Reactions in\ngiant module (%)") +
  theme_no_labels

# --- PANEL: ACR ---
p_acr <- ggplot(df_plot, aes(x = clade_factor, y = n_acr)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, 
               linewidth = box_line_width, 
               median.linewidth = box_line_width,
               width = bar_width, fill = "white",
               na.rm = TRUE) + # EXPLICITLY IGNORE MISSING DATA
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)), n.breaks = 4) +
  coord_flip(xlim = c(0.5, n_tips + 0.5), clip = "off") +
  labs(x = NULL, y = "ACR\nmetabolites") +
  theme_no_labels

################################################################################
## 4. COMBINE AND SAVE
################################################################################

cat("Combining panels...\n")

p_tree_a <- p_tree + ggtitle("A") + theme(plot.title = element_text(face = "bold", size = 28, hjust = 0))
p_tree_b <- p_tree + ggtitle("B") + theme(plot.title = element_text(face = "bold", size = 28, hjust = 0))

# Layout: Tree (0.6) | Data (1) | Data (1) | Data (1) | Data (1)
layout_spec <- c(0.6, 1, 1, 1, 1)

row_a <- p_tree_a + p_species + p_reactions + p_metabolites + p_complexes +
  plot_layout(widths = layout_spec)

row_b <- p_tree_b + p_balanced + p_concordant + p_module + p_acr +
  plot_layout(widths = layout_spec)

combined <- row_a / row_b + plot_layout(heights = c(1, 1))

cat("Saving figure...\n")
dir.create("../../docs/thesis/figures", showWarnings = FALSE, recursive = TRUE)

fig_height <- max(14, length(tree_tips) * 1.1) * 2

# Lower panel only: tree once on the left, then all data panels in a single row
theme_poster_labels <- theme_with_labels +
  theme(
    axis.text.y   = element_text(size = font_size * 1.25, hjust = 1, margin = margin(r = 5)),
    axis.text.x   = element_text(size = font_size * 1.15),
    axis.title.x  = element_text(face = "plain", size = font_size * 1.35),
    axis.ticks.length = unit(0.25, "cm")
  )
p_module_labeled <- p_module + theme(
  axis.text.x   = element_text(size = font_size * 1.15),
  axis.title.x  = element_text(face = "plain", size = font_size * 1.35),
  axis.ticks.length = unit(0.25, "cm")
)
p_acr_poster <- p_acr + theme(
  axis.text.x   = element_text(size = font_size * 1.15),
  axis.title.x  = element_text(face = "plain", size = font_size * 1.35),
  axis.ticks.length = unit(0.25, "cm")
)
p_balanced_poster <- p_balanced + theme_poster_labels
p_concordant_poster <- p_concordant + theme(
  axis.text.x   = element_text(size = font_size * 1.15),
  axis.title.x  = element_text(face = "plain", size = font_size * 1.35),
  axis.ticks.length = unit(0.25, "cm")
)
lower_panel <- p_tree + p_balanced + p_concordant + p_module + p_acr +
  plot_layout(widths = c(0.6, 1.2, 1.2, 1.2, 1.2))

panel_w <- 22
panel_h <- fig_height / 2

# Upper panel (row A: tree + species + reactions + metabolites + complexes)
ggsave(file.path(FIGURES_DIR, "phylo_network_robustness_upper_panel.svg"), row_a,
       width = panel_w, height = panel_h, units = "in", device = svg)
ggsave(file.path(FIGURES_DIR, "phylo_network_robustness_upper_panel.pdf"), row_a,
       width = panel_w, height = panel_h, units = "in", device = cairo_pdf)
ggsave(file.path(FIGURES_DIR, "phylo_network_robustness_upper_panel.png"), row_a,
       width = panel_w, height = panel_h, units = "in", dpi = 300, bg = "white")

# Lower panel (row B: tree + balanced + concordant + module + acr)
ggsave(file.path(FIGURES_DIR, "phylo_network_robustness_lower_panel.svg"), lower_panel,
       width = panel_w, height = panel_h, units = "in", device = svg)
ggsave(file.path(FIGURES_DIR, "phylo_network_robustness_lower_panel.pdf"), lower_panel,
       width = panel_w, height = panel_h, units = "in", device = cairo_pdf)
ggsave(file.path(FIGURES_DIR, "phylo_network_robustness_lower_panel.png"), lower_panel,
       width = panel_w, height = panel_h, units = "in", dpi = 300, bg = "white")

ggsave(file.path(FIGURES_DIR, "phylo_network_robustness.pdf"), combined,
       width = 20, height = 28, units = "in", device = cairo_pdf)

ggsave(file.path(FIGURES_DIR, "phylo_network_robustness.svg"), combined,
       width = 20, height = 28, units = "in", device = svg)

ggsave(file.path(FIGURES_DIR, "phylo_network_robustness.png"), combined,
       width = 20, height = 28, units = "in", dpi = 300, bg = "white")

cat("Done.\n")

