library(tidyverse)
library(seriation)

cat("\n", rep("=", 80), "\n", sep = "")
cat("PATHWAY ENRICHMENT HEATMAP - SUPPLEMENTARY (all pathways)\n")
cat(rep("=", 80), "\n\n", sep = "")

library(here)
REPO_ROOT   <- here::here()
DATA_DIR    <- file.path(REPO_ROOT, "data")
RESULTS_DIR <- file.path(REPO_ROOT, "results")
FIGURES_DIR <- file.path(REPO_ROOT, "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)
source(file.path(REPO_ROOT, "scripts/figures/theme.R"))

# =============================================================================
# 1. LOAD DATA
# =============================================================================

cat("Loading data...\n")

df_meta <- read_csv(file.path(DATA_DIR, "annotations/cocoa_results_merged.csv"), show_col_types = FALSE) %>%
  select(old_species_id, `Major clade`)

df_pathway_enrichment <- read_csv(file.path(RESULTS_DIR, "pathway/pathway_enrichment_detailed.csv"),
                                  show_col_types = FALSE)

cat("  Loaded", nrow(df_pathway_enrichment), "species×pathway enrichment values\n")
cat("  Across", n_distinct(df_pathway_enrichment$species_id), "species\n")
cat("  Across", n_distinct(df_pathway_enrichment$pathway), "pathways\n\n")

# Use all pathways (no curated selection)
df_thesis <- df_pathway_enrichment

# =============================================================================
# 2. ORDER SPECIES BY CLADE (PHYLOGENETIC)
# =============================================================================

cat("Ordering species by phylogenetic clades...\n")

clade_order <- c(
  "Saccharomycetaceae",
  "Saccharomycodaceae",
  "Phaffomycetaceae",
  "CUG-Ser2",
  "CUG-Ser1",
  "CUG-Ala",
  "Pichiaceae",
  "Sporopachydermia",
  "Alloascoideaceae",
  "Dipodascaceae/Trichomonascaceae",
  "Trigonopsidaceae",
  "Lipomycetaceae",
  "Outgroup"
)

temp_mat <- df_thesis %>%
  select(pathway, species_id, enrichment) %>%
  pivot_wider(names_from = species_id, values_from = enrichment, values_fill = 0) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

enrichment_matrix_for_clustering <- t(temp_mat)

order_species_by_similarity <- function(species_ids, enrichment_matrix) {
  if (length(species_ids) <= 2) return(species_ids)
  sub_mat <- enrichment_matrix[species_ids, , drop = FALSE]
  sub_mat[is.na(sub_mat)] <- 0
  d <- dist(sub_mat, method = "euclidean")
  if (any(is.na(d)) || any(is.infinite(d))) return(species_ids)
  hc <- hclust(d, method = "average")
  species_ids[hc$order]
}

species_with_clade <- df_meta %>%
  filter(old_species_id %in% colnames(temp_mat)) %>%
  select(old_species_id, `Major clade`) %>%
  mutate(`Major clade` = factor(`Major clade`, levels = clade_order))

ordered_species_list <- list()
for (clade in clade_order) {
  clade_species <- species_with_clade %>%
    filter(`Major clade` == clade) %>%
    pull(old_species_id)
  if (length(clade_species) > 0) {
    ordered_species_list[[clade]] <- order_species_by_similarity(clade_species, enrichment_matrix_for_clustering)
  }
}

species_order_final <- unlist(ordered_species_list)
cat("  Ordered", length(species_order_final), "species across",
    length(ordered_species_list), "clades\n\n")

clade_boundaries <- species_with_clade %>%
  mutate(old_species_id = as.character(old_species_id)) %>%
  slice(match(species_order_final, old_species_id)) %>%
  mutate(species_idx = row_number()) %>%
  group_by(`Major clade`) %>%
  summarise(
    start    = min(species_idx),
    end      = max(species_idx),
    midpoint = (min(species_idx) + max(species_idx)) / 2,
    n_species = n(),
    .groups  = "drop"
  ) %>%
  arrange(start)

# =============================================================================
# 3. ORDER PATHWAYS BY HIERARCHICAL CLUSTERING
# =============================================================================

cat("Ordering pathways by hierarchical clustering...\n")

temp_mat_full <- df_thesis %>%
  select(pathway, species_id, enrichment) %>%
  pivot_wider(names_from = pathway, values_from = enrichment, values_fill = 0) %>%
  column_to_rownames("species_id") %>%
  as.matrix()

mat_for_clustering <- t(temp_mat_full)
for (i in seq_len(nrow(mat_for_clustering))) {
  na_idx <- is.na(mat_for_clustering[i, ])
  if (any(na_idx)) mat_for_clustering[i, na_idx] <- 0
}

pathway_dist   <- dist(mat_for_clustering, method = "manhattan")
pathway_hc     <- hclust(pathway_dist, method = "average")
pathway_order_seriate <- seriate(pathway_dist, method = "OLO", control = list(hclust = pathway_hc))
pathway_order_idx     <- get_order(pathway_order_seriate)
pathway_order  <- rev(rownames(mat_for_clustering)[pathway_order_idx])

cat("  Final pathway count:", length(pathway_order), "\n\n")

# =============================================================================
# 4. BUILD MATRIX
# =============================================================================

cat("Creating pathway × species matrix...\n")

mat <- df_thesis %>%
  select(pathway, species_id, enrichment) %>%
  pivot_wider(names_from = species_id, values_from = enrichment, values_fill = NA) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

mat <- mat[pathway_order, species_order_final, drop = FALSE]

rownames(mat) <- rownames(mat) %>%
  str_remove(" - yeast$") %>%
  str_replace(" - other$", " (other)")

cat("  Matrix:", nrow(mat), "pathways ×", ncol(mat), "species\n")

mat_long <- mat %>%
  as.data.frame() %>%
  rownames_to_column("pathway") %>%
  pivot_longer(-pathway, names_to = "species_id", values_to = "enrichment") %>%
  mutate(
    pathway    = factor(pathway, levels = rownames(mat)),
    species_id = factor(species_id, levels = colnames(mat)),
    pathway_idx  = as.numeric(pathway),
    species_idx  = as.numeric(species_id)
  )

max_abs_enrichment <- max(abs(range(mat, na.rm = TRUE)))
enrichment_limits  <- c(-max_abs_enrichment, max_abs_enrichment)
cat("  Enrichment range: [", round(enrichment_limits[1], 3), ", ",
    round(enrichment_limits[2], 3), "]\n\n", sep = "")

# =============================================================================
# 5. CREATE HEATMAP
# =============================================================================

cat("Creating heatmap...\n")

dir.create(file.path(FIGURES_DIR, "supplementary", showWarnings = FALSE, recursive = TRUE)

p_heatmap <- ggplot(mat_long, aes(x = species_idx, y = pathway_idx)) +

  geom_tile(aes(fill = enrichment), color = NA) +

  scale_fill_gradient2(
    name     = "Enrichment",
    low      = CB_BLUE,
    mid      = "#FFFFFF",
    high     = CB_VERMILLON,
    midpoint = 0,
    limits   = enrichment_limits,
    breaks   = seq(enrichment_limits[1], enrichment_limits[2], length.out = 5),
    labels   = sprintf("%.2f", seq(enrichment_limits[1], enrichment_limits[2], length.out = 5)),
    na.value = "#696969"
  ) +

  geom_vline(
    data = clade_boundaries %>% filter(row_number() != 1),
    aes(xintercept = start - 0.5),
    color = "black", linewidth = 0.1
  ) +

  scale_x_continuous(
    expand = c(0, 0),
    breaks = clade_boundaries$midpoint,
    labels = NULL,
    position = "bottom"
  ) +

  geom_text(
    data = {
      clade_labels <- clade_boundaries
      clade_labels$x_pos <- clade_labels$midpoint
      clade_labels$`Major clade` <- gsub("Dipodascaceae/Trichomonascaceae",
                                         "Dipodascaceae/\nTrichomonascaceae",
                                         clade_labels$`Major clade`)
      sporo_idx <- which(grepl("Sporopachydermia", clade_labels$`Major clade`))
      if (length(sporo_idx) > 0)
        clade_labels$x_pos[sporo_idx] <- clade_labels$x_pos[sporo_idx] - 3
      allo_idx <- which(grepl("Alloascoideaceae", clade_labels$`Major clade`))
      if (length(allo_idx) > 0)
        clade_labels$x_pos[allo_idx] <- clade_labels$x_pos[allo_idx] + 6
      clade_labels
    },
    aes(x = x_pos, y = 0, label = `Major clade`),
    angle = 45, hjust = 1, vjust = 1,
    size = 3.0, lineheight = 0.5, inherit.aes = FALSE
  ) +

  scale_y_continuous(
    expand = c(0, 0),
    breaks = seq(1, nrow(mat)),
    labels = str_wrap(rownames(mat), width = 28),
    position = "left"
  ) +

  coord_cartesian(ylim = c(0.5, nrow(mat) + 0.5),
                  xlim = c(0.5, ncol(mat) + 0.5), clip = "off") +

  labs(x = NULL, y = NULL) +

  theme_thesis_heatmap() +
  theme(
    axis.text.x       = element_blank(),
    axis.text.y       = element_text(angle = 0, hjust = 1, size = 7, lineheight = 0.8),
    legend.position   = "top",
    legend.title      = element_text(margin = margin(r = 12, unit = "mm")),
    legend.text       = element_text(margin = margin(t = 2, l = 1, r = 2, unit = "mm")),
    legend.key.height = unit(4, "mm"),
    legend.key.width  = unit(20, "mm"),
    plot.margin       = margin(t = 2, r = 2, b = 25, l = 2, unit = "mm"),
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

ggsave(
  file.path(FIGURES_DIR, "supplementary/pathway_enrichment_supplementary.pdf"),
  p_heatmap,
  width  = 6,
  height = 14,
  device = cairo_pdf
)
cat("  ✓ Saved PDF\n")

ggsave(
  file.path(FIGURES_DIR, "supplementary/pathway_enrichment_supplementary.png"),
  p_heatmap,
  width  = 6,
  height = 14,
  dpi    = 300,
  bg     = "white"
)
cat("  ✓ Saved PNG\n")

