library(tidyverse)
library(ggplot2)
library(ape)  # For reading tree

cat("\n", rep("=", 80), "\n", sep = "")
cat("CORE COMPLEX PROPERTY IDENTITY HEATMAP\n")
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

# Load property matrix
prop_matrix <- read_csv(file.path(RESULTS_DIR, "concordance/core_complex_property_matrix.csv"),
                        show_col_types = FALSE)
cat("  Property matrix: ", nrow(prop_matrix), " complexes × ",
    ncol(prop_matrix) - 1, " species\n", sep = "")

# Load conservation scores
conservation <- read_csv(file.path(RESULTS_DIR, "concordance/core_complex_consistency.csv"),
                         show_col_types = FALSE)
cat("  Conservation scores: ", nrow(conservation), " complexes\n", sep = "")

# Load species metadata
metadata <- read_csv(file.path(DATA_DIR, "annotations/cocoa_results_merged.csv"),
                     show_col_types = FALSE)
cat("  Species metadata: ", nrow(metadata), " species\n\n", sep = "")

# =============================================================================
# 2. ORDER SPECIES BY CLADE (PHYLOGENETIC)
# =============================================================================

cat("Ordering species by phylogenetic clades...\n")

# Define clade order from tree (rooted, so Outgroup is last)
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

# First, identify core complexes (present in all species)
n_species_total <- ncol(prop_matrix) - 1
core_complex_mask <- apply(prop_matrix %>% select(-Complex) %>% as.matrix(), 1,
                           function(row) all(row != 5))
core_complexes_for_clustering <- prop_matrix$Complex[core_complex_mask]

cat("  Using ", length(core_complexes_for_clustering),
    " core complexes for species clustering\n", sep = "")

# Extract core complex matrix for clustering (species as rows, complexes as columns)
core_matrix_for_clustering <- prop_matrix %>%
  filter(Complex %in% core_complexes_for_clustering) %>%
  select(-Complex) %>%
  as.matrix() %>%
  t()  # Transpose so species are rows

# Function to order species within a clade by hierarchical clustering
order_species_by_similarity <- function(species_ids, core_matrix) {
  if (length(species_ids) <= 2) {
    return(species_ids)  # No need to reorder 1-2 species
  }

  # Get subset of matrix for these species
  species_matrix <- core_matrix[species_ids, , drop = FALSE]

  # Compute distance matrix (Manhattan distance works well for categorical-ish data)
  dist_matrix <- dist(species_matrix, method = "manhattan")

  # Hierarchical clustering
  hc <- hclust(dist_matrix, method = "average")

  # Return species in dendrogram order
  species_ids[hc$order]
}

# Order species within each clade by concordance pattern similarity
species_in_matrix <- colnames(prop_matrix)[-1]
species_with_clade <- metadata %>%
  filter(old_species_id %in% species_in_matrix) %>%
  select(old_species_id, `Major clade`) %>%
  mutate(`Major clade` = factor(`Major clade`, levels = clade_order))

# Process each clade
ordered_species_list <- list()
for (clade in clade_order) {
  clade_species <- species_with_clade %>%
    filter(`Major clade` == clade) %>%
    pull(old_species_id)

  if (length(clade_species) > 0) {
    ordered_clade_species <- order_species_by_similarity(clade_species, core_matrix_for_clustering)
    ordered_species_list[[clade]] <- ordered_clade_species
  }
}

# Combine into final order
species_order_final <- unlist(ordered_species_list, use.names = FALSE)

# Create species_ordered dataframe in the new order
species_ordered <- metadata %>%
  filter(old_species_id %in% species_order_final) %>%
  mutate(old_species_id = factor(old_species_id, levels = species_order_final)) %>%
  arrange(old_species_id) %>%
  mutate(old_species_id = as.character(old_species_id))

cat("  Ordered ", nrow(species_ordered), " species across ",
    length(unique(species_ordered$`Major clade`)), " clades\n", sep = "")
cat("  (species within each clade ordered by concordance pattern similarity)\n\n")

# Calculate clade boundaries for vertical lines (sorted by position in plot)
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
  arrange(start)  # Sort by position to get correct phylogenetic order

cat("Clade sizes (in phylogenetic order):\n")
for (i in 1:nrow(clade_boundaries)) {
  cat(sprintf("  %-35s: %3d species\n",
              as.character(clade_boundaries$`Major clade`[i]),
              clade_boundaries$n_species[i]))
}
cat("\n")

# =============================================================================
# 3. ORDER COMPLEXES BY CONCORDANCE PROPORTION
# =============================================================================

cat("Ordering complexes by proportion of concordance types...\n")

# Calculate presence (number of species where complex is present = not code 5)
# and proportion of concordance types (values 1-4, anything but None/0)
mat_values <- prop_matrix %>%
  select(-Complex) %>%
  as.matrix()

complex_stats <- tibble(
  Complex = prop_matrix$Complex,
  Presence = apply(mat_values, 1, function(row) sum(row != 5)),
  # Proportion of concordance types (1=Concordant, 2=Trivially_concordant, 3=Balanced, 4=Trivially_balanced)
  # among present values (excluding 5=Not present)
  Concordance_Prop = apply(mat_values, 1, function(row) {
    present <- row[row != 5]
    if (length(present) == 0) return(0)
    sum(present %in% c(1, 2, 3, 4)) / length(present)
  })
)

# Calculate dominant type from the matrix (including None as a possible dominant type)
# If None is the most prevalent value, return "None"
# Otherwise return the most prevalent concordance type
# New codes: 0=None, 1=Concordant, 2=Trivially concordant, 3=Balanced, 4=Trivially balanced, 5=Not present
calc_dominant_type <- function(row_values) {
  # Exclude "Not present" (5)
  present_values <- row_values[row_values != 5]
  if (length(present_values) == 0) return("None")

  # Count all types including None (0)
  counts <- table(present_values)
  dominant_code <- as.numeric(names(which.max(counts)))

  if (dominant_code == 0) return("None")
  else if (dominant_code == 1) return("Concordant")
  else if (dominant_code == 2) return("Trivially concordant")
  else if (dominant_code == 3) return("Balanced")
  else if (dominant_code == 4) return("Trivially balanced")
  else return("None")
}

matrix_dominant_type <- apply(mat_values, 1, calc_dominant_type)

complex_properties_df <- tibble(
  Complex = prop_matrix$Complex,
  Dominant_Type = matrix_dominant_type
)

# Join presence and property data
# Order: Trivially balanced > Balanced > Trivially concordant > Concordant > None
property_order <- c("Trivially balanced", "Balanced", "Trivially concordant", "Concordant", "None")

complex_ordered <- complex_stats %>%
  left_join(complex_properties_df, by = "Complex") %>%
  mutate(
    Dominant = factor(Dominant_Type, levels = property_order),
    Score = Concordance_Prop  # Use concordance proportion as the score
  ) %>%
  # Sort by: 1) Presence (descending), 2) Dominant type (incl. None), 3) Concordance proportion (descending)
  arrange(desc(Presence), Dominant, desc(Concordance_Prop)) %>%
  select(Complex, Presence, Dominant, Score, Concordance_Prop)

cat("  Presence range: ",
    min(complex_ordered$Presence), " to ",
    max(complex_ordered$Presence), " species\n", sep = "")

# Count complexes by presence category
cat("\n  Complexes by presence:\n")
cat(sprintf("    Core (all %d species): %d complexes\n",
            ncol(prop_matrix) - 1,
            sum(complex_ordered$Presence == (ncol(prop_matrix) - 1))))
cat(sprintf("    High (≥300 species): %d complexes\n",
            sum(complex_ordered$Presence >= 300)))
cat(sprintf("    Medium (100-299 species): %d complexes\n",
            sum(complex_ordered$Presence >= 100 & complex_ordered$Presence < 300)))
cat(sprintf("    Low (<100 species): %d complexes\n",
            sum(complex_ordered$Presence < 100)))

# Count complexes by property type
cat("\n  Complexes by dominant property:\n")
property_counts <- complex_ordered %>%
  count(Dominant)

for (i in 1:nrow(property_counts)) {
  cat(sprintf("    %-12s: %4d complexes\n",
              as.character(property_counts$Dominant[i]),
              property_counts$n[i]))
}
cat("\n")

# =============================================================================
# 4. FILTER TO CORE COMPLEXES AND PREPARE MATRIX FOR PLOTTING
# =============================================================================

cat("Filtering to core complexes (present in all species)...\n")

# Filter to core complexes (present in all species)
n_species_total <- ncol(prop_matrix) - 1
core_complexes <- complex_ordered %>%
  filter(Presence == n_species_total) %>%
  pull(Complex)

cat(sprintf("  Found %d core complexes (present in all %d species)\n",
            length(core_complexes), n_species_total))

# Reorder matrix columns (species) and rows (complexes) - CORE ONLY
species_cols <- as.character(species_ordered$old_species_id)
complex_rows <- complex_ordered %>%
  filter(Complex %in% core_complexes) %>%
  pull(Complex)

# Extract matrix values in correct order (core complexes only)
mat_data <- prop_matrix %>%
  filter(Complex %in% core_complexes) %>%
  select(Complex, all_of(species_cols)) %>%
  mutate(Complex = factor(Complex, levels = complex_rows)) %>%
  arrange(Complex)

# Convert to long format for ggplot
mat_long <- mat_data %>%
  pivot_longer(-Complex, names_to = "Species", values_to = "Property") %>%
  mutate(
    Species = factor(Species, levels = species_cols),
    Complex_idx = as.numeric(Complex),
    Species_idx = as.numeric(Species)
  )

# Convert property codes to factor with labels
# Codes: 0=None, 1=Concordant, 2=Trivially concordant, 3=Balanced, 4=Trivially balanced, 5=Not present
mat_long <- mat_long %>%
  mutate(Property_label = factor(
    Property,
    levels = c(0, 1, 2, 3, 4, 5),
    labels = c("None", "Concordant", "Trivially concordant", "Balanced", "Trivially balanced", "Not present")
  ))

cat("  Matrix dimensions: ", length(unique(mat_long$Complex)), " × ",
    length(unique(mat_long$Species)), "\n", sep = "")
cat("  Total cells: ", nrow(mat_long), "\n\n", sep = "")

# =============================================================================
# 5. CREATE HEATMAP
# =============================================================================

cat("Creating heatmap...\n")

# Colorblind-friendly palette (Wong / Okabe-Ito — matches thesis_theme.R)
property_colors <- c(
  "Trivially balanced"   = CB_YELLOW,    # #F0E442 — yellow
  "Balanced"             = CB_ORANGE,    # #E69F00 — amber-orange
  "Trivially concordant" = CB_SKYBLUE,   # #56B4E9 — sky blue
  "Concordant"           = CB_BLUE,      # #0072B2 — deep blue
  "None"                 = "#999999",    # gray
  "Not present"          = "#FFFFFF"     # white
)

# Create heatmap
p_heatmap <- ggplot(mat_long, aes(x = Species_idx, y = Complex_idx)) +
  geom_tile(aes(fill = Property_label), color = NA) +

  # Color scale
  scale_fill_manual(
    values = property_colors,
    name = NULL,
    na.value = "white"
  ) +

  # Add vertical lines for clade boundaries (thin black)
  geom_vline(
    data = clade_boundaries %>% filter(row_number() != 1),
    aes(xintercept = start - 0.5),
    color = "black",
    linewidth = 0.15
  ) +

  # Axes - matching enrichment heatmap style
  scale_x_continuous(
    expand = c(0, 0),
    breaks = clade_boundaries$midpoint,
    labels = NULL,  # Remove default labels, add manually below
    position = "bottom"
  ) +

  scale_y_continuous(
    expand = c(0, 0),
    breaks = NULL
  ) +

  # Labels
  labs(
    x = NULL,
    y = "Core Complexes"
  ) +

  # Theme
  theme_thesis_heatmap() +
  theme(
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank(),
    legend.position    = "top",
    legend.title       = element_blank(),
    legend.key.width   = unit(5, "mm"),
    legend.key.height  = unit(5, "mm"),
    legend.spacing.x   = unit(2, "mm"),
    legend.text        = element_text(margin = margin(l = 2, r = 6, unit = "mm")),
    plot.margin        = margin(t = 5, r = 5, b = 28, l = 5, unit = "mm"),
  ) +
  guides(fill = guide_legend(
    override.aes = list(color = "grey40", linewidth = 0.3)
  )) +

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
    aes(x = x_pos, y = -16, label = `Major clade`),
    angle = 45,
    hjust = 1,
    vjust = 1,
    size = 3.5,
    lineheight = 0.5,
    inherit.aes = FALSE
  ) +

  coord_cartesian(ylim = c(0.5, length(unique(mat_long$Complex)) + 0.5),
                  xlim = c(0.5, length(unique(mat_long$Species)) + 0.5),
                  clip = "off")

cat("  Heatmap created\n\n")

# =============================================================================
# 6. SAVE FIGURE
# =============================================================================

cat("Saving figures...\n")

dir.create("../../docs/thesis/figures", showWarnings = FALSE, recursive = TRUE)

# A4-optimized dimensions (similar to pathway enrichment heatmap)
fig_width <- 7.5   # inches (suitable for A4 width)
fig_height <- 7    # inches (compact format for A4)

# High-resolution PNG
ggsave(
  file.path(FIGURES_DIR, "core_property_heatmap1.png"),
  p_heatmap,
  width = fig_width,
  height = fig_height,
  dpi = 300,
  bg = "white"
)
cat(sprintf("  Saved: core_property_heatmap1.png (%.1f\" × %.1f\", 300 dpi - A4 optimized)\n",
            fig_width, fig_height))

# Vector PDF for publication
ggsave(
  file.path(FIGURES_DIR, "core_property_heatmap1.pdf"),
  p_heatmap,
  width = fig_width,
  height = fig_height,
  device = cairo_pdf
)
cat("  Saved: core_property_heatmap.pdf (vector - A4 optimized)\n\n")

# SVG for web/external use
ggsave(
  file.path(FIGURES_DIR, "core_property_heatmap1.svg"),
  p_heatmap,
  width = fig_width,
  height = fig_height,
  device = svg
)
cat("  Saved: core_property_heatmap1.svg (vector - web/external use)\n\n")

# =============================================================================
# 7. SUMMARY STATISTICS
# =============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("SUMMARY STATISTICS\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Matrix dimensions:\n")
cat("  Complexes: ", nrow(mat_data), "\n", sep = "")
cat("  Species: ", length(species_cols), "\n", sep = "")
cat("  Total cells: ", nrow(mat_long), "\n\n", sep = "")

cat("Property distribution (core complexes):\n")
prop_counts <- mat_long %>%
  filter(!is.na(Property_label)) %>%
  count(Property_label) %>%
  mutate(Percentage = 100 * n / sum(n))

for (i in 1:nrow(prop_counts)) {
  cat(sprintf("  %-22s: %7d (%5.2f%%)\n",
              as.character(prop_counts$Property_label[i]),
              prop_counts$n[i],
              prop_counts$Percentage[i]))
}
cat("\n")

cat("Clade structure:\n")
cat("  Total clades: ", nrow(clade_boundaries), "\n", sep = "")
cat("  Species per clade:\n")
cat(sprintf("    Min: %d\n", min(clade_boundaries$n_species)))
cat(sprintf("    Max: %d\n", max(clade_boundaries$n_species)))
cat(sprintf("    Mean: %.1f\n\n", mean(clade_boundaries$n_species)))

cat("Conservation scores (for core complexes):\n")
core_scores <- complex_ordered %>% filter(Complex %in% core_complexes)
cat(sprintf("  Highly conserved (≥0.99): %d\n",
            sum(core_scores$Score >= 0.99, na.rm = TRUE)))
cat(sprintf("  Moderately conserved (0.80-0.99): %d\n",
            sum(core_scores$Score >= 0.80 & core_scores$Score < 0.99, na.rm = TRUE)))
cat(sprintf("  Variable (<0.80): %d\n\n",
            sum(core_scores$Score < 0.80 & core_scores$Score > 0, na.rm = TRUE)))

# =============================================================================
# 8. SESSION INFO
# =============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("COMPLETE\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Output files:\n")
cat("  - results/figures/main/core_property_heatmap1.png (core complexes only)\n")
cat("  - results/figures/main/core_property_heatmap1.pdf (core complexes only)\n\n")

cat("Property codes (ConcordanceType enum):\n")
cat("  0 = None (singleton)\n")
cat("  1 = Concordant\n")
cat("  2 = Trivially_concordant\n")
cat("  3 = Balanced\n")
cat("  4 = Trivially_balanced\n")
cat("  5 = Not present\n\n")

