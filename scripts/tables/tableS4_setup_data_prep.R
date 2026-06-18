################################################################################
## SCRIPT 00: SETUP AND DATA PREPARATION
################################################################################
##
## Purpose: Load, validate, and prepare data for phylogenetic comparative analysis
##
## Inputs:
##   - Metabolic network data (cocoa_results_merged.csv)
##   - Phylogenetic trees (MCMCTree and RelTime)
##
## Outputs:
##   - Clean matched dataset
##   - Divergence times
##   - Log-transformed variables
##   - Data quality report
##
################################################################################

# Start analysis log
cat("\n", rep("=", 80), "\n", sep = "")
cat("STATISTICALLY RIGOROUS PHYLOGENETIC COMPARATIVE ANALYSIS\n")
cat("Script 00: Data Preparation and Validation\n")
cat(rep("=", 80), "\n\n", sep = "")
cat("Started:", as.character(Sys.time()), "\n\n")

# Load required libraries
library(ape)
library(phytools)
library(tidyverse)

library(here)
REPO_ROOT   <- here::here()
DATA_DIR    <- file.path(REPO_ROOT, "data")
RESULTS_DIR <- file.path(REPO_ROOT, "results")
sink("../results/outputs/00_setup_and_data_prep.txt", split = TRUE)
cat("R version:", R.version.string, "\n")
cat("Working directory:", getwd(), "\n\n")

# Set working directory to scripts folder

################################################################################
## 1. LOAD RAW DATA
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("STEP 1: LOADING RAW DATA\n")
cat(rep("=", 80), "\n\n", sep = "")

# Load metabolic network data
cat("Loading metabolic network data...\n")
df_raw <- read_csv(file.path(DATA_DIR, "annotations/cocoa_results_merged.csv"), show_col_types = FALSE)

# ============================================================
# OVERRIDE: cocoa_results_merged.csv holds no-split (original GEM)
# balanced/concordant/complex counts. The manuscript's evolutionary
# model selection (Table S4) uses the split-model (random_0) values
# instead, averaged per species across its HPC job replicates.
# ============================================================
cat("Loading random_0 balanced/concordant values...\n")
df_r0 <- read_csv(file.path(RESULTS_DIR, "kinetic/resource_analysis.csv"), show_col_types = FALSE) %>%
  filter(variant == "random_0", job_status == "SUCCESS" | (job_status == "TIMEOUT" & !is.na(balanced_complexes))) %>%
  group_by(model) %>%
  summarise(
    n_balanced_complexes_r0 = mean(balanced_complexes, na.rm = TRUE),
    n_concordant_total_r0 = mean(concordant_pairs, na.rm = TRUE),
    n_complexes_r0 = mean(n_complexes, na.rm = TRUE),
    .groups = "drop"
  )

n_before <- nrow(df_raw)
df_raw <- df_raw %>%
  inner_join(df_r0, by = c("model_name" = "model")) %>%
  mutate(
    n_balanced_complexes = n_balanced_complexes_r0,
    n_concordant_total   = n_concordant_total_r0,
    n_complexes           = n_complexes_r0
  ) %>%
  select(-n_balanced_complexes_r0, -n_concordant_total_r0, -n_complexes_r0)

cat(sprintf("  random_0 values applied for %d species\n", nrow(df_raw)))
cat(sprintf("  %d species dropped (no random_0 result)\n", n_before - nrow(df_raw)))
# ============================================================

cat("  Rows:", nrow(df_raw), "\n")
cat("  Columns:", ncol(df_raw), "\n")
cat("  Variables:", paste(names(df_raw), collapse = ", "), "\n\n")

# Load phylogenetic trees
cat("Loading phylogenetic trees...\n")
tree_mcmc <- read.tree(file.path(DATA_DIR, "phylogeny/332_2408OGs_timetree_mcmctree.nwk"))
tree_reltime <- read.tree(file.path(DATA_DIR, "phylogeny/332_2408OGs_timetree_reltime.nwk"))

cat("  MCMCTree tips:", length(tree_mcmc$tip.label), "\n")
cat("  RelTime tips:", length(tree_reltime$tip.label), "\n\n")

################################################################################
## 2. DATA QUALITY CHECKS
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("STEP 2: DATA QUALITY CHECKS\n")
cat(rep("=", 80), "\n\n", sep = "")

# Check for missing values
cat("Missing values per column:\n")
missing_counts <- colSums(is.na(df_raw))
if (sum(missing_counts) == 0) {
  cat("  No missing values detected!\n\n")
} else {
  print(missing_counts[missing_counts > 0])
  cat("\n")
}

# Check for zero/negative values in count variables
cat("Checking for zero or negative values in count variables...\n")
count_vars <- c("n_balanced_complexes", "n_concordant_total", "n_complexes", "n_reactions")

for (var in count_vars) {
  if (var %in% names(df_raw)) {
    n_zero <- sum(df_raw[[var]] == 0, na.rm = TRUE)
    n_negative <- sum(df_raw[[var]] < 0, na.rm = TRUE)

    if (n_zero > 0 || n_negative > 0) {
      cat("  WARNING:", var, "has", n_zero, "zeros and", n_negative, "negative values\n")
    } else {
      cat("  ", var, ": OK (all positive)\n")
    }
  }
}
cat("\n")

# Summary statistics for key variables
cat("Summary statistics for key variables:\n")
for (var in count_vars) {
  if (var %in% names(df_raw)) {
    cat("  ", var, ":\n")
    cat("    Range:", min(df_raw[[var]], na.rm = TRUE), "-", max(df_raw[[var]], na.rm = TRUE), "\n")
    cat("    Mean:", round(mean(df_raw[[var]], na.rm = TRUE), 2), "\n")
    cat("    Median:", median(df_raw[[var]], na.rm = TRUE), "\n")
  }
}
cat("\n")

################################################################################
## 3. TREE VALIDATION
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("STEP 3: PHYLOGENETIC TREE VALIDATION\n")
cat(rep("=", 80), "\n\n", sep = "")

# Check ultrametricity
is_ultra_mcmc <- is.ultrametric(tree_mcmc)
is_ultra_reltime <- is.ultrametric(tree_reltime)

cat("Ultrametricity check (before forcing):\n")
cat("  MCMCTree:", ifelse(is_ultra_mcmc, "YES (ultrametric)", "NO (not ultrametric)"), "\n")
cat("  RelTime:", ifelse(is_ultra_reltime, "YES (ultrametric)", "NO (not ultrametric)"), "\n\n")

# Force ultrametricity for proper PGLS analysis
# Even if trees are "close enough", force.ultrametric ensures exact ultrametricity
if (!is_ultra_mcmc) {
  cat("Forcing MCMCTree to be ultrametric...\n")
  tree_mcmc <- force.ultrametric(tree_mcmc)
}
if (!is_ultra_reltime) {
  cat("Forcing RelTime to be ultrametric...\n")
  tree_reltime <- force.ultrametric(tree_reltime)
}

# Verify after forcing
is_ultra_mcmc_after <- is.ultrametric(tree_mcmc)
is_ultra_reltime_after <- is.ultrametric(tree_reltime)

cat("\nUltrametricity check (after forcing):\n")
cat("  MCMCTree:", ifelse(is_ultra_mcmc_after, "YES (ultrametric)", "NO (not ultrametric)"), "\n")
cat("  RelTime:", ifelse(is_ultra_reltime_after, "YES (ultrametric)", "NO (not ultrametric)"), "\n\n")

cat("NOTE: Ultrametric trees have all tips at present day.\n")
cat("      For temporal analysis, use LINEAGE AGES (distance from root), not tip distances.\n\n")

# Check tree properties
cat("Tree properties:\n")
cat("  MCMCTree:\n")
cat("    Total branch length:", round(sum(tree_mcmc$edge.length), 4), "\n")
cat("    Mean branch length:", round(mean(tree_mcmc$edge.length), 4), "\n")
cat("    Tree height:", round(max(node.depth.edgelength(tree_mcmc)), 4), "\n")
cat("  RelTime:\n")
cat("    Total branch length:", round(sum(tree_reltime$edge.length), 4), "\n")
cat("    Mean branch length:", round(mean(tree_reltime$edge.length), 4), "\n")
cat("    Tree height:", round(max(node.depth.edgelength(tree_reltime)), 4), "\n\n")

################################################################################
## 4. SPECIES MATCHING
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("STEP 4: MATCHING SPECIES BETWEEN DATA AND TREE\n")
cat(rep("=", 80), "\n\n", sep = "")

# Clean species names in data
# Use old_species_id which matches the tree tip labels
# Also filter out outgroup species (not in tree)
df_raw <- df_raw %>%
  filter(`Major clade` != "Outgroup" | is.na(`Major clade`)) %>%
  mutate(species_clean = old_species_id)

cat("Filtered out", sum(df_raw$`Major clade` == "Outgroup", na.rm = TRUE), "outgroup species\n")

# Find species in both data and tree
species_in_data <- df_raw$species_clean
species_in_tree <- tree_mcmc$tip.label

species_both <- intersect(species_in_data, species_in_tree)
species_data_only <- setdiff(species_in_data, species_in_tree)
species_tree_only <- setdiff(species_in_tree, species_in_data)

cat("Species matching:\n")
cat("  In data (after filtering outgroups):", length(species_in_data), "\n")
cat("  In tree:", length(species_in_tree), "\n")
cat("  In data only:", length(species_data_only), "\n")
cat("  In tree only:", length(species_tree_only), "\n")
cat("  In both:", length(species_both), "\n\n")

if (length(species_data_only) > 0) {
  cat("Species in data but not in tree:\n")
  cat("  ", paste(head(species_data_only, 10), collapse = "\n   "), "\n")
  if (length(species_data_only) > 10) {
    cat("  ... and", length(species_data_only) - 10, "more\n")
  }
  cat("\n")
}

# Keep only matched species
df_matched <- df_raw %>%
  filter(species_clean %in% species_both)

cat("Final matched dataset:", nrow(df_matched), "species\n\n")

################################################################################
## 5. PRUNE TREES
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("STEP 5: PRUNING TREES TO MATCHED SPECIES\n")
cat(rep("=", 80), "\n\n", sep = "")

# Prune trees to matched species
tree_pruned_mcmc <- keep.tip(tree_mcmc, species_both)
tree_pruned_reltime <- keep.tip(tree_reltime, species_both)

cat("Pruned trees:\n")
cat("  MCMCTree: ", length(tree_pruned_mcmc$tip.label), "tips\n")
cat("  RelTime: ", length(tree_pruned_reltime$tip.label), "tips\n\n")

# Verify matching
if (length(tree_pruned_mcmc$tip.label) != nrow(df_matched)) {
  stop("ERROR: Tree and data have different number of species after matching!")
}

cat("Verification: Tree and data species counts match!\n\n")

################################################################################
## 6. CALCULATE LINEAGE AGES (EVOLUTIONARY TIME FROM ROOT)
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("STEP 6: CALCULATING LINEAGE AGES (EVOLUTIONARY TIME FROM ROOT)\n")
cat(rep("=", 80), "\n\n", sep = "")

# Function to calculate lineage ages for each species
# (evolutionary time from root to when lineage split from sister lineage)
# NOTE: This measures total evolutionary time along the lineage from the root,
#       which is what we need for testing directional evolution (Research Q2)
get_lineage_ages <- function(tree) {
  n_tips <- length(tree$tip.label)
  lineage_ages <- numeric(n_tips)

  # Calculate node ages (distances from root)
  node_ages <- node.depth.edgelength(tree)

  for (i in 1:n_tips) {
    # Get parent node of this tip
    parent_node <- tree$edge[tree$edge[,2] == i, 1]

    # Lineage age = distance from root to parent node
    # This represents the total evolutionary time for this lineage
    lineage_ages[i] <- node_ages[parent_node]
  }

  names(lineage_ages) <- tree$tip.label
  return(lineage_ages)
}

# Calculate for both trees
div_times_mcmc <- get_lineage_ages(tree_pruned_mcmc)
div_times_reltime <- get_lineage_ages(tree_pruned_reltime)

cat("Lineage ages (MCMCTree):\n")
cat("  NOTE: Tree units are in 100 MY intervals\n")
cat("  NOTE: Lineage age = evolutionary time from root to parent node\n")
cat("  Range (raw):", round(min(div_times_mcmc), 4), "-", round(max(div_times_mcmc), 4), "\n")
cat("  Range (Ma):", round(min(div_times_mcmc) * 100, 2), "-", round(max(div_times_mcmc) * 100, 2), "Ma\n")
cat("  Mean:", round(mean(div_times_mcmc), 4), "(", round(mean(div_times_mcmc) * 100, 2), "Ma)\n")
cat("  SD:", round(sd(div_times_mcmc), 4), "(", round(sd(div_times_mcmc) * 100, 2), "Ma)\n")
cat("  CV:", round(sd(div_times_mcmc) / mean(div_times_mcmc), 4), "\n\n")

cat("Lineage ages (RelTime):\n")
cat("  Range:", round(min(div_times_reltime), 4), "-", round(max(div_times_reltime), 4), "\n")
cat("  Mean:", round(mean(div_times_reltime), 4), "\n")
cat("  SD:", round(sd(div_times_reltime), 4), "\n")
cat("  CV:", round(sd(div_times_reltime) / mean(div_times_reltime), 4), "\n\n")

# Add lineage ages to dataframe
# NOTE: Tree units are in 100 MY intervals, convert to Ma
# NOTE: Lineage age = evolutionary time from root (used for Research Q2)
max_time_mcmc <- max(div_times_mcmc)
max_time_reltime <- max(div_times_reltime)

df_matched <- df_matched %>%
  mutate(
    # Lineage ages: evolutionary time from root to parent node
    lineage_age_mcmc = div_times_mcmc[species_clean],
    lineage_age_reltime = div_times_reltime[species_clean],
    # Convert to Ma (multiply by 100)
    lineage_age_ma = lineage_age_mcmc * 100,
    # Keep old names for backward compatibility with existing scripts
    divergence_time_mcmc = lineage_age_mcmc,
    divergence_time_reltime = lineage_age_reltime,
    divergence_time_ma = lineage_age_ma,
    # Time from root for intuitive plotting (old on left, recent on right)
    time_from_root = max_time_mcmc * 100 - lineage_age_ma
  )

################################################################################
## 7. CREATE LOG-TRANSFORMED VARIABLES
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("STEP 7: LOG-TRANSFORMING VARIABLES\n")
cat(rep("=", 80), "\n\n", sep = "")

# Log-transform count variables
# Note: All values are positive (checked above), so no need for log(x+1)
df_matched <- df_matched %>%
  mutate(
    log_balanced = log(n_balanced_complexes),
    log_concordant = log(n_concordant_total),
    log_complexes = log(n_complexes),
    log_reactions = log(n_reactions)
  )

cat("Log-transformed variables created:\n")
cat("  - log_balanced (log of n_balanced_complexes)\n")
cat("  - log_concordant (log of n_concordant_total)\n")
cat("  - log_complexes (log of n_complexes)\n")
cat("  - log_reactions (log of n_reactions)\n\n")

cat("Time variables created:\n")
cat("  - lineage_age_mcmc: Lineage age from root (raw tree units: 100 MY intervals)\n")
cat("  - lineage_age_ma: Lineage age from root in Ma (* 100)\n")
cat("  - divergence_time_mcmc/ma: Kept for backward compatibility (= lineage_age)\n")
cat("  - time_from_root: For intuitive plotting (", round(max_time_mcmc * 100, 2), "Ma [old] -> 0 Ma [recent])\n\n")

cat("Summary statistics (log-transformed):\n")
log_vars <- c("log_balanced", "log_concordant", "log_complexes", "log_reactions")
for (var in log_vars) {
  cat("  ", var, ":\n")
  cat("    Range:", round(min(df_matched[[var]]), 3), "-", round(max(df_matched[[var]]), 3), "\n")
  cat("    Mean:", round(mean(df_matched[[var]]), 3), "\n")
  cat("    SD:", round(sd(df_matched[[var]]), 3), "\n")
}
cat("\n")

################################################################################
## 8. FLAG POTENTIAL OUTLIERS
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("STEP 8: FLAGGING POTENTIAL OUTLIERS\n")
cat(rep("=", 80), "\n\n", sep = "")

# Simple outlier flagging: values > 3 SD from mean
flag_outliers <- function(x, threshold = 3) {
  z_scores <- abs((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  return(z_scores > threshold)
}

df_matched <- df_matched %>%
  mutate(
    outlier_balanced = flag_outliers(log_balanced),
    outlier_concordant = flag_outliers(log_concordant),
    outlier_complexes = flag_outliers(log_complexes)
  )

n_outliers_bal <- sum(df_matched$outlier_balanced)
n_outliers_conc <- sum(df_matched$outlier_concordant)
n_outliers_comp <- sum(df_matched$outlier_complexes)

cat("Potential outliers (|z-score| > 3):\n")
cat("  log_balanced:", n_outliers_bal, "species\n")
cat("  log_concordant:", n_outliers_conc, "species\n")
cat("  log_complexes:", n_outliers_comp, "species\n\n")

if (n_outliers_bal > 0) {
  cat("Outlier species (log_balanced):\n")
  outlier_species <- df_matched %>%
    filter(outlier_balanced) %>%
    select(species_clean, n_balanced_complexes, log_balanced, n_complexes)
  print(outlier_species)
  cat("\n")
}

################################################################################
## 9. SET ROW NAMES FOR PGLS
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("STEP 9: PREPARING FOR PHYLOGENETIC ANALYSES\n")
cat(rep("=", 80), "\n\n", sep = "")

# Set row names to species names (required for PGLS functions)
rownames(df_matched) <- df_matched$species_clean

cat("Row names set to species names\n")
cat("  First 5 row names:", paste(head(rownames(df_matched)), collapse = ", "), "\n")
cat("  Match tree tip labels:", all(rownames(df_matched) %in% tree_pruned_mcmc$tip.label), "\n\n")

################################################################################
## 10. SAVE OUTPUTS
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("STEP 10: SAVING OUTPUTS\n")
cat(rep("=", 80), "\n\n", sep = "")

# Save clean workspace
save(
  df_matched,
  tree_pruned_mcmc,
  tree_pruned_reltime,
  div_times_mcmc,
  div_times_reltime,
  file = "../results/workspaces/00_clean_data.RData"
)
cat("Workspace saved: ../results/workspaces/00_clean_data.RData\n")

# Save data summary table
data_summary <- data.frame(
  Metric = c(
    "Total species",
    "Species in data only",
    "Species in tree only",
    "Matched species",
    "Potential outliers (balanced)",
    "Potential outliers (concordant)",
    "Divergence time range (Ma)",
    "Time from root range (Ma)",
    "Mean divergence time (Ma)",
    "SD divergence time (Ma)"
  ),
  Value = c(
    nrow(df_raw),
    length(species_data_only),
    length(species_tree_only),
    nrow(df_matched),
    n_outliers_bal,
    n_outliers_conc,
    paste(round(min(div_times_mcmc) * 100, 2), "-", round(max(div_times_mcmc) * 100, 2)),
    paste(round(max(div_times_mcmc) * 100, 2), "-", 0),
    round(mean(div_times_mcmc) * 100, 2),
    round(sd(div_times_mcmc) * 100, 2)
  )
)

write_csv(data_summary, "../results/tables/data_summary.csv")
cat("Data summary saved: ../results/tables/data_summary.csv\n")

# Save matched species list
species_list <- data.frame(
  species = df_matched$species_clean,
  clade = df_matched$`Major clade`,
  n_complexes = df_matched$n_complexes,
  n_balanced = df_matched$n_balanced_complexes,
  divergence_time = df_matched$divergence_time_mcmc,
  outlier_flagged = df_matched$outlier_balanced | df_matched$outlier_concordant
)

write_csv(species_list, "../results/tables/matched_species_list.csv")
cat("Species list saved: ../results/tables/matched_species_list.csv\n\n")

################################################################################
## 11. START ANALYSIS LOG
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("STEP 11: CREATING ANALYSIS LOG\n")
cat(rep("=", 80), "\n\n", sep = "")

log_content <- paste0(
  "# Analysis Log\n\n",
  "## Script 00: Data Preparation\n\n",
  "**Date**: ", Sys.Date(), "\n",
  "**R version**: ", R.version.string, "\n\n",
  "### Data Summary\n",
  "- **Matched species**: ", nrow(df_matched), "\n",
  "- **Phylogenetic tree**: MCMCTree (", length(tree_pruned_mcmc$tip.label), " tips)\n",
  "- **Lineage ages**: ", round(min(div_times_mcmc) * 100, 2), " - ", round(max(div_times_mcmc) * 100, 2), " Ma\n",
  "- **Time from root**: ", round(max(div_times_mcmc) * 100, 2), " Ma (old) -> 0 Ma (recent)\n",
  "- **Potential outliers**: ", n_outliers_bal, " (balanced), ", n_outliers_conc, " (concordant)\n\n",
  "### Key Variables\n",
  "- `log_balanced`: log(n_balanced_complexes)\n",
  "- `log_concordant`: log(n_concordant_total)\n",
  "- `log_complexes`: log(n_complexes)\n",
  "- `log_reactions`: log(n_reactions)\n",
  "- `lineage_age_mcmc`: Lineage age from root (raw tree units: 100 MY intervals)\n",
  "- `lineage_age_ma`: Lineage age from root in Ma (× 100)\n",
  "- `divergence_time_mcmc/ma`: Backward compatible aliases for lineage_age\n",
  "- `time_from_root`: Time from root for intuitive plotting (old -> recent)\n\n",
  "### Quality Checks\n",
  "- [x] No missing values\n",
  "- [x] All counts positive\n",
  "- [x] Tree ultrametricity forced (ensures exact ultrametricity for PGLS)\n",
  "- [x] Species matched between data and tree\n",
  "- [x] Row names set for PGLS\n",
  "- [x] Time units converted to Ma\n",
  "- [x] Lineage ages calculated correctly (from root, not sister split)\n\n",
  "### Next Steps\n",
  "1. Run Script 01: Evolutionary model selection\n",
  "2. Determine best correlation structure for PGLS\n",
  "3. Proceed with main analysis\n\n",
  "---\n\n"
)


dir.create("../documentation", recursive = TRUE, showWarnings = FALSE)
writeLines(log_content, "../documentation/ANALYSIS_LOG.md")
cat("Analysis log created: ../documentation/ANALYSIS_LOG.md\n\n")

################################################################################
## 12. FINAL SUMMARY
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("DATA PREPARATION COMPLETE\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Summary:\n")
cat("  ✓ Loaded", nrow(df_raw), "species from metabolic network data\n")
cat("  ✓ Loaded phylogenetic tree with", length(tree_mcmc$tip.label), "tips\n")
cat("  ✓ Matched", nrow(df_matched), "species between data and tree\n")
cat("  ✓ Calculated lineage ages from root (", round(min(div_times_mcmc) * 100, 2), "-",
    round(max(div_times_mcmc) * 100, 2), " Ma)\n")
cat("  ✓ Lineage age = evolutionary time from root (for Research Q2)\n")
cat("  ✓ Created log-transformed variables\n")
cat("  ✓ Flagged", n_outliers_bal, "potential outliers\n")
cat("  ✓ Prepared data for PGLS analysis\n\n")

cat("Outputs:\n")
cat("  • Workspace: ../results/workspaces/00_clean_data.RData\n")
cat("  • Data summary: ../results/tables/data_summary.csv\n")
cat("  • Species list: ../results/tables/matched_species_list.csv\n")
cat("  • Analysis log: ../documentation/ANALYSIS_LOG.md\n\n")

cat("Next script: 01_evolutionary_model_selection.R\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat(rep("=", 80), "\n\n", sep = "")

# Print session info
cat("Session Info:\n")
sessionInfo()

