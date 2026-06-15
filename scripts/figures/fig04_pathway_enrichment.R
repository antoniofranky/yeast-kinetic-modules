library(tidyverse)
library(seriation)

cat("\n", rep("=", 80), "\n", sep = "")
cat("PATHWAY ENRICHMENT HEATMAP - THESIS VERSION (30 pathways)\n")
cat(rep("=", 80), "\n\n", sep = "")

library(here)
REPO_ROOT   <- here::here()
DATA_DIR    <- file.path(REPO_ROOT, "data")
RESULTS_DIR <- file.path(REPO_ROOT, "results")
FIGURES_DIR <- file.path(REPO_ROOT, "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)
source(file.path(REPO_ROOT, "scripts/figures/theme.R"))

# =============================================================================
# 1. DEFINE THESIS PATHWAY SELECTION (30 pathways aligned with Langary et al.)
# =============================================================================

cat("Defining thesis pathway selection (30 pathways)...\n\n")

# 10 ENRICHED (in giant kinetic module → concentration robustness)
enriched_pathways <- c(
  "Autophagy - other",
  "Autophagy - yeast",
  "Phagosome",
  "Glycerophospholipid metabolism",
  "Inositol phosphate metabolism",
  "Glycerolipid metabolism",
  "Steroid biosynthesis",
  "C5-Branched dibasic acid metabolism",
  "MAPK signaling pathway - yeast",
  "Oxidative phosphorylation"
)

# 10 HIGH VARIANCE (evolutionary adaptation)
# Clade-specific (5) + Species-specific (5)
high_variance_pathways <- c(
  # Clade-specific
  "Thiamine metabolism",
  "Pentose phosphate pathway",
  "Peroxisome",
  "Fatty acid biosynthesis",
  "One carbon pool by folate",
  # Species-specific
  "Meiosis - yeast",
  "Sphingolipid metabolism",
  "Butanoate metabolism",
  "Nicotinate and nicotinamide metabolism",
  "Pantothenate and CoA biosynthesis"
)

# 10 DEPLETED (outside giant module → less robust)
depleted_pathways <- c(
  "Fatty acid degradation",
  "Lysine biosynthesis",
  "Arginine biosynthesis",
  "Phenylalanine, tyrosine and tryptophan biosynthesis",
  "Tryptophan metabolism",
  "Citrate cycle (TCA cycle)",
  "Glyoxylate and dicarboxylate metabolism",
  "Pyruvate metabolism",
  "Purine metabolism",
  "Pyrimidine metabolism"
)

# Combine all 30 pathways
thesis_pathways <- c(enriched_pathways, high_variance_pathways, depleted_pathways)

cat("Selected pathway categories:\n")
cat("  - Enriched (in giant module):     ", length(enriched_pathways), "\n")
cat("  - High variance (evolution):      ", length(high_variance_pathways), "\n")
cat("  - Depleted (outside giant module):", length(depleted_pathways), "\n")
cat("  - TOTAL:                          ", length(thesis_pathways), "\n\n")

# =============================================================================
# 2. LOAD DATA
# =============================================================================

cat("Loading data...\n")

# Metadata for clade information
df_meta <- read_csv(file.path(DATA_DIR, "annotations/cocoa_results_merged.csv"), show_col_types = FALSE) %>%
  select(old_species_id, `Major clade`)

# Pre-computed enrichment (same formula: observed_rate - background_rate,
# background_rate = n_giant_reactions / n_total_split_reactions per species)
df_pathway_enrichment <- read_csv(file.path(RESULTS_DIR, "pathway/pathway_enrichment_detailed.csv"),
                                  show_col_types = FALSE)

cat("  Loaded", nrow(df_pathway_enrichment), "species×pathway enrichment values\n")
cat("  Across", n_distinct(df_pathway_enrichment$species_id), "species\n\n")

# =============================================================================
# 3. ENRICHMENT ALREADY COMPUTED
# =============================================================================

cat("Using pre-computed pathway enrichment...\n")

cat("  Calculated enrichment for", n_distinct(df_pathway_enrichment$pathway),
    "pathways\n")

# Filter to thesis pathways only
df_thesis <- df_pathway_enrichment %>%
  filter(pathway %in% thesis_pathways)

cat("  Filtered to", n_distinct(df_thesis$pathway), "thesis pathways\n")

# Check which pathways are missing
missing_pathways <- setdiff(thesis_pathways, unique(df_thesis$pathway))
if (length(missing_pathways) > 0) {
  cat("\n  WARNING: These pathways not found in data:\n")
  for (p in missing_pathways) {
    cat("    -", p, "\n")
  }
  cat("\n")
}

# =============================================================================
# 4. ORDER SPECIES BY CLADE (PHYLOGENETIC)
# =============================================================================

cat("Ordering species by phylogenetic clades...\n")

# Define clade order
clade_order <- c(
  "Saccharomycetaceae",
  "Saccharomycodaceae",
  "Phaffomycetaceae",
  "CUG-Ser2",
  "CUG-Ser1",
  "CUG-Ala",
  "Pichiaceae",
  "Sporopachydermia",  # NOTE: metadata has "Sporopachydermia" not "Sporopachydermia clade"
  "Alloascoideaceae",
  "Dipodascaceae/Trichomonascaceae",
  "Trigonopsidaceae",
  "Lipomycetaceae",
  "Outgroup"
)

# Create matrix for ordering
temp_mat <- df_thesis %>%
  select(pathway, species_id, enrichment) %>%
  pivot_wider(names_from = species_id, values_from = enrichment, values_fill = 0) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

enrichment_matrix_for_clustering <- t(temp_mat)

# Function to order species within clade
order_species_by_similarity <- function(species_ids, enrichment_matrix) {
  if (length(species_ids) <= 2) {
    return(species_ids)
  }

  species_matrix <- enrichment_matrix[species_ids, , drop = FALSE]
  dist_matrix <- dist(species_matrix, method = "manhattan")
  hc <- hclust(dist_matrix, method = "average")
  species_order <- seriate(dist_matrix, method = "OLO", control = list(hclust = hc))
  species_order_idx <- get_order(species_order)

  species_ids[species_order_idx]
}

# Order species within each clade
species_in_matrix <- colnames(temp_mat)
species_with_clade <- df_meta %>%
  filter(old_species_id %in% species_in_matrix) %>%
  select(old_species_id, `Major clade`) %>%
  mutate(`Major clade` = factor(`Major clade`, levels = clade_order))

ordered_species_list <- list()
for (clade in clade_order) {
  clade_species <- species_with_clade %>%
    filter(`Major clade` == clade) %>%
    pull(old_species_id)

  if (length(clade_species) > 0) {
    ordered_clade_species <- order_species_by_similarity(clade_species, enrichment_matrix_for_clustering)
    ordered_species_list[[clade]] <- ordered_clade_species
  }
}

species_order_final <- unlist(ordered_species_list, use.names = FALSE)

cat("  Ordered", length(species_order_final), "species across",
    length(ordered_species_list), "clades\n\n")

# Calculate clade boundaries
species_ordered <- df_meta %>%
  filter(old_species_id %in% species_order_final) %>%
  mutate(old_species_id = factor(old_species_id, levels = species_order_final)) %>%
  arrange(old_species_id) %>%
  mutate(old_species_id = as.character(old_species_id))

clade_boundaries <- species_ordered %>%
  mutate(species_idx = row_number()) %>%
  group_by(`Major clade`) %>%
  summarise(
    start = min(species_idx),
    end = max(species_idx),
    midpoint = (min(species_idx) + max(species_idx)) / 2,
    n_species = n(),
    .groups = "drop"
  ) %>%
  arrange(start)

# =============================================================================
# 5. ORDER PATHWAYS BY HIERARCHICAL CLUSTERING (like enrichment heatmap)
# =============================================================================

cat("Ordering pathways by hierarchical clustering...\n")

# First, create full matrix for clustering
temp_mat_full <- df_thesis %>%
  select(pathway, species_id, enrichment) %>%
  pivot_wider(names_from = pathway, values_from = enrichment, values_fill = 0) %>%
  column_to_rownames("species_id") %>%
  as.matrix()

# Only keep pathways in our thesis selection
pathways_in_data <- intersect(thesis_pathways, colnames(temp_mat_full))
temp_mat_full <- temp_mat_full[, pathways_in_data, drop = FALSE]

# Transpose for pathway clustering (pathways as rows)
mat_for_clustering <- t(temp_mat_full)

# Handle NA values
for (i in 1:nrow(mat_for_clustering)) {
  na_idx <- is.na(mat_for_clustering[i, ])
  if (any(na_idx)) {
    mat_for_clustering[i, na_idx] <- 0
  }
}

# Hierarchical clustering with optimal leaf ordering
pathway_dist <- dist(mat_for_clustering, method = "manhattan")
pathway_hc <- hclust(pathway_dist, method = "average")

cat("  Applying optimal leaf ordering...\n")
pathway_order_seriate <- seriate(pathway_dist, method = "OLO", control = list(hclust = pathway_hc))
pathway_order_idx <- get_order(pathway_order_seriate)

# Get ordered pathway names
pathway_order <- rownames(mat_for_clustering)[pathway_order_idx]

# Reverse order so high enrichment pathways are at top
pathway_order <- rev(pathway_order)

cat("  Final pathway count:", length(pathway_order), "\n\n")

# =============================================================================
# 6. CREATE MATRIX AND PREPARE DATA
# =============================================================================

cat("Creating pathway × species matrix...\n")

# Matrix: pathways as rows, species as columns
mat <- df_thesis %>%
  select(pathway, species_id, enrichment) %>%
  pivot_wider(names_from = species_id, values_from = enrichment, values_fill = NA) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

# Reorder
mat <- mat[pathway_order, species_order_final, drop = FALSE]

# Clean pathway labels: remove organism suffixes redundant in a yeast thesis
rownames(mat) <- rownames(mat) %>%
  str_remove(" - yeast$") %>%
  str_replace(" - other$", " (other)") %>%
  str_replace("Phenylalanine, tyrosine and tryptophan biosynthesis",
              "Phe, Tyr & Trp biosynthesis") %>%
  str_replace("Glyoxylate and dicarboxylate metabolism",
              "Glyoxylate & dicarboxylate metabolism")

cat("  Matrix:", nrow(mat), "pathways ×", ncol(mat), "species\n")

# Convert to long format for ggplot
mat_long <- mat %>%
  as.data.frame() %>%
  rownames_to_column("pathway") %>%
  pivot_longer(-pathway, names_to = "species_id", values_to = "enrichment") %>%
  mutate(
    pathway = factor(pathway, levels = rownames(mat)),
    species_id = factor(species_id, levels = colnames(mat)),
    pathway_idx = as.numeric(pathway),
    species_idx = as.numeric(species_id)
  )

# Determine symmetric color scale range around 0
max_abs_enrichment <- max(abs(range(mat, na.rm = TRUE)))
enrichment_limits <- c(-max_abs_enrichment, max_abs_enrichment)

cat("  Enrichment range: [", round(enrichment_limits[1], 3), ", ",
    round(enrichment_limits[2], 3), "]\n\n", sep = "")

# =============================================================================
# 7. CREATE HEATMAP WITH GGPLOT2
# =============================================================================

cat("Creating heatmap...\n")

dir.create("../../docs/thesis/figures", showWarnings = FALSE, recursive = TRUE)

# Create heatmap
p_heatmap <- ggplot(mat_long, aes(x = species_idx, y = pathway_idx)) +

  # Main heatmap tiles
  geom_tile(aes(fill = enrichment), color = NA) +

  # Diverging color scale: Wong/Okabe-Ito blue–white–vermillion
  # Both endpoints are colorblind-safe (deuteranopia, protanopia, tritanopia)
  scale_fill_gradient2(
    name = "Enrichment",
    low = CB_BLUE,       # #0072B2 — Wong blue (depleted)
    mid = "#FFFFFF",
    high = CB_VERMILLON, # #D55E00 — Wong vermillion (enriched)
    midpoint = 0,
    limits = enrichment_limits,
    breaks = seq(enrichment_limits[1], enrichment_limits[2], length.out = 5),
    labels = sprintf("%.2f", seq(enrichment_limits[1], enrichment_limits[2], length.out = 5)),
    na.value = "#696969"
  ) +

  # Thin black clade separators (drawn after tiles, minimal overdraw at this width)
  geom_vline(
    data = clade_boundaries %>% filter(row_number() != 1),
    aes(xintercept = start - 0.5),
    color = "black",
    linewidth = 0.1
  ) +

  # Axes - keep original breaks, add labels manually
  scale_x_continuous(
    expand = c(0, 0),
    breaks = clade_boundaries$midpoint,
    labels = NULL,  # Remove default labels
    position = "bottom"
  ) +

  # Add clade labels manually with offset for Alloascoideaceae
  geom_text(
    data = {
      clade_labels <- clade_boundaries
      clade_labels$x_pos <- clade_labels$midpoint
      # Add line break for Dipodascaceae/Trichomonascaceae
      clade_labels$`Major clade` <- gsub("Dipodascaceae/Trichomonascaceae",
                                          "Dipodascaceae/\nTrichomonascaceae",
                                          clade_labels$`Major clade`)
      # Offset Sporopachydermia label left
      sporo_idx <- which(grepl("Sporopachydermia", clade_labels$`Major clade`))
      if (length(sporo_idx) > 0) {
        clade_labels$x_pos[sporo_idx] <- clade_labels$x_pos[sporo_idx] - 3
      }
      # Offset Alloascoideaceae label right
      allo_idx <- which(grepl("Alloascoideaceae", clade_labels$`Major clade`))
      if (length(allo_idx) > 0) {
        clade_labels$x_pos[allo_idx] <- clade_labels$x_pos[allo_idx] + 6
      }
      clade_labels
    },
    aes(x = x_pos, y = 0, label = `Major clade`),
    angle = 45,
    hjust = 1,
    vjust = 1,
    size = 3.0,
    lineheight = 0.5,
    inherit.aes = FALSE
  ) +

  scale_y_continuous(
    expand = c(0, 0),
    breaks = seq(1, nrow(mat)),
    labels = str_wrap(rownames(mat), width = 38),
    position = "left"
  ) +

  coord_cartesian(ylim = c(0.5, nrow(mat) + 0.5), xlim = c(0.5, ncol(mat) + 0.5), clip = "off") +

  labs(x = NULL, y = NULL) +

  # Theme
  theme_thesis_heatmap() +
  theme(
    axis.text.x        = element_blank(),
    axis.text.y        = element_text(angle = 0, hjust = 1, size = 7.5, lineheight = 0.65),
    legend.position    = "top",
    legend.title       = element_text(margin = margin(r = 12, unit = "mm")),
    legend.text        = element_text(margin = margin(t = 2, l = 1, r = 2, unit = "mm")),
    legend.key.height  = unit(4, "mm"),
    legend.key.width   = unit(20, "mm"),
    plot.margin        = margin(t = 2, r = 2, b = 25, l = 2, unit = "mm"),
  ) +
  guides(fill = guide_colorbar(
    title.position = "left",
    title.vjust    = 0.85,
    barwidth       = unit(40, "mm"),
    barheight      = unit(4, "mm"),
    ticks          = TRUE,
    frame.colour   = "grey40",
    ticks.colour   = "grey40"
  ))

cat("  Heatmap created\n")

# Save PDF
ggsave(
  file.path(FIGURES_DIR, "pathway_enrichment_thesis.pdf"),
  p_heatmap,
  width = 8,
  height = 5.5,
  device = cairo_pdf
)
cat("  ✓ Saved PDF (8\" × 5.5\" - optimized for A4 landscape)\n")

# Save PNG
ggsave(
  file.path(FIGURES_DIR, "pathway_enrichment_thesis.png"),
  p_heatmap,
  width = 8,
  height = 5.5,
  dpi = 300,
  bg = "white"
)
cat("  ✓ Saved PNG (8\" × 5.5\" - optimized for A4 landscape)\n")

# SVG for web/external use
ggsave(
  file.path(FIGURES_DIR, "pathway_enrichment_thesis.svg"),
  p_heatmap,
  width = 8,
  height = 5.5,
  device = svg
)
cat("  ✓ Saved SVG (7.5\" × 8\" - vector for web/external use)\n")

# =============================================================================
# 8. SUMMARY STATISTICS
# =============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("SUMMARY STATISTICS\n")
cat(rep("=", 80), "\n\n", sep = "")

# Calculate enrichment statistics by category
pathway_categories <- tibble(
  pathway = c(enriched_pathways, high_variance_pathways, depleted_pathways),
  category = c(
    rep("Enriched", length(enriched_pathways)),
    rep("High Variance", length(high_variance_pathways)),
    rep("Depleted", length(depleted_pathways))
  )
)

category_stats <- df_thesis %>%
  left_join(pathway_categories, by = "pathway") %>%
  filter(!is.na(category)) %>%
  group_by(category) %>%
  summarise(
    n_pathways = n_distinct(pathway),
    mean_enrichment = mean(enrichment, na.rm = TRUE),
    median_enrichment = median(enrichment, na.rm = TRUE),
    sd_enrichment = sd(enrichment, na.rm = TRUE),
    min_enrichment = min(enrichment, na.rm = TRUE),
    max_enrichment = max(enrichment, na.rm = TRUE),
    .groups = "drop"
  )

cat("Pathway category statistics:\n\n")
print(category_stats)

# Save pathway selection
pathway_selection <- pathway_categories %>%
  mutate(order = 1:n())

write_csv(
  pathway_selection,
  file.path(RESULTS_DIR, "kinetic/thesis_pathway_selection_final.csv")
)

cat("\n\nSaved:\n")
cat("  - ../../docs/thesis/figures/pathway_enrichment_thesis.pdf (7.5\" × 6\" - A4 optimized)\n")
cat("  - ../../docs/thesis/figures/pathway_enrichment_thesis.png (7.5\" × 6\" - A4 optimized)\n")
cat("  - results/kinetic_phylogenetic/thesis_pathway_selection_final.csv\n")

cat("\nDone!\n")

