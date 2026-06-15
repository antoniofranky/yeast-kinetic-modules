################################################################################
## SCRIPT 01: EVOLUTIONARY MODEL SELECTION
################################################################################
##
## Purpose: Determine the best evolutionary model for each trait
##          This is CRITICAL - it determines the correlation structure for PGLS
##
## Models tested:
##   - Brownian Motion (BM): Neutral drift
##   - Ornstein-Uhlenbeck (OU): Stabilizing selection
##   - Early Burst (EB): Rapid early evolution
##   - Pagel's Lambda: Variable phylogenetic signal
##
## Note: Drift/mean_trend model excluded because it requires non-ultrametric trees.
##       Directional evolution is tested in Script 02 via lineage_age_mcmc covariate.
##
## Traits tested:
##   Count traits (log-transformed):
##     - log_balanced, log_concordant, log_reactions
##     - log_giant_reactions, log_acr
##   Proportion traits (logit-transformed):
##     - logit_prop_balanced, logit_prop_concordant
##     - logit_prop_giant, logit_prop_acr
##
## Selection criteria:
##   - AICc (small-sample corrected)
##   - Akaike weights
##   - Parameter plausibility
##   - If ΔAICc < 2, prefer simpler model (BM)
##
################################################################################

library(ape)
library(geiger)
library(tidyverse)

library(here)
REPO_ROOT   <- here::here()
DATA_DIR    <- file.path(REPO_ROOT, "data")
RESULTS_DIR <- file.path(REPO_ROOT, "results")

sink("../results/outputs/01_evolutionary_model_selection.txt", split = TRUE)

cat("\n", rep("=", 80), "\n", sep = "")
cat("Script 01: Evolutionary Model Selection\n")
cat(rep("=", 80), "\n\n", sep = "")
cat("Started:", as.character(Sys.time()), "\n\n")

################################################################################
## 1. LOAD DATA
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("STEP 1: LOADING DATA\n")
cat(rep("=", 80), "\n\n", sep = "")

load("../results/workspaces/00_clean_data.RData")

cat("Loaded workspace with:\n")
cat("  Species:", nrow(df_matched), "\n")
cat("  Tree tips:", length(tree_pruned_mcmc$tip.label), "\n\n")

# Load kinetic data
base_dir <- REPO_ROOT

kinetic <- read_csv(
  file.path(base_dir, "results/analysis/results/kinetic_phylogenetic/species_kinetic_summary.csv.gz"),
  show_col_types = FALSE
) %>%
  dplyr::select(species_id, n_acr)

giant_module <- read_csv(
  file.path(base_dir, "results/analysis/results/kinetic_phylogenetic/giant_module_reaction_summary.csv.gz"),
  show_col_types = FALSE
) %>%
  dplyr::select(species_id, n_giant_reactions)

split_sizes <- read_csv(
  file.path(base_dir, "results/analysis/results/kinetic_phylogenetic/split_model_sizes.csv.gz"),
  show_col_types = FALSE
) %>%
  dplyr::select(species_id, n_total_split_reactions, n_total_split_metabolites)

# Merge kinetic data with df_matched
kinetic_merged <- kinetic %>%
  inner_join(giant_module, by = "species_id") %>%
  inner_join(split_sizes, by = "species_id")

df_matched <- df_matched %>%
  left_join(kinetic_merged, by = c("species_clean" = "species_id"))

n_kinetic <- sum(!is.na(df_matched$n_giant_reactions))
cat("Merged kinetic data:\n")
cat("  Species with kinetic data:", n_kinetic, "of", nrow(df_matched), "\n\n")

# Create log-transformed kinetic traits
df_matched <- df_matched %>%
  mutate(
    log_giant_reactions = log(n_giant_reactions),
    log_acr = log(n_acr)
  )

# Create proportions
df_matched <- df_matched %>%
  mutate(
    prop_balanced = n_balanced_complexes / n_complexes,
    prop_concordant = n_concordant_total / choose(n_complexes, 2),
    prop_giant = n_giant_reactions / n_total_split_reactions,
    prop_acr = n_acr / n_total_split_metabolites
  )

# Logit-transform proportions (with offset to avoid Inf)
logit <- function(p, eps = 1e-6) {
  p <- pmin(pmax(p, eps), 1 - eps)
  log(p / (1 - p))
}

df_matched <- df_matched %>%
  mutate(
    logit_prop_balanced = logit(prop_balanced),
    logit_prop_concordant = logit(prop_concordant),
    logit_prop_giant = logit(prop_giant),
    logit_prop_acr = logit(prop_acr)
  )

cat("Computed proportions:\n")
cat("  prop_balanced:   mean =", round(mean(df_matched$prop_balanced, na.rm = TRUE), 4), "\n")
cat("  prop_concordant: mean =", round(mean(df_matched$prop_concordant, na.rm = TRUE), 6), "\n")
cat("  prop_giant:      mean =", round(mean(df_matched$prop_giant, na.rm = TRUE), 4), "\n")
cat("  prop_acr:        mean =", round(mean(df_matched$prop_acr, na.rm = TRUE), 4), "\n\n")

# Size-corrected residuals: regress out log(n_total_split_reactions) from log(n_giant_reactions)
# Uses OLS (not PGLS) for size correction, consistent with phylogenetic signal analysis
df_giant_resid <- df_matched %>%
  filter(!is.na(n_giant_reactions), !is.na(n_total_split_reactions)) %>%
  mutate(log_total_split = log(n_total_split_reactions))

lm_giant <- lm(log_giant_reactions ~ log_total_split, data = df_giant_resid)
cat("Size correction: log(n_giant_reactions) ~ log(n_total_split_reactions)\n")
cat("  Slope:", round(coef(lm_giant)[2], 4), "\n")
cat("  R²:   ", round(summary(lm_giant)$r.squared, 4), "\n\n")

df_giant_resid$resid_log_giant <- residuals(lm_giant)
df_matched <- df_matched %>%
  left_join(df_giant_resid %>% select(species_clean, resid_log_giant),
            by = "species_clean")

################################################################################
## 2. FIT EVOLUTIONARY MODELS
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("STEP 2: FITTING EVOLUTIONARY MODELS\n")
cat(rep("=", 80), "\n\n", sep = "")

# Traits to test: count traits + proportion traits
traits <- c(
  # Count traits (log-transformed)
  "log_balanced", "log_concordant", "log_reactions",
  "log_giant_reactions", "log_acr",
  # Proportion traits (logit-transformed)
  "logit_prop_balanced", "logit_prop_concordant",
  "logit_prop_giant", "logit_prop_acr",
  # Size-corrected residuals
  "resid_log_giant"
)

# Storage for results
all_results <- list()

for (trait in traits) {
  cat("\n", rep("-", 80), "\n", sep = "")
  cat("TRAIT:", trait, "\n")
  cat(rep("-", 80), "\n\n", sep = "")

  # Extract trait values
  trait_values <- df_matched[[trait]]
  names(trait_values) <- df_matched$species_clean

  # Remove NA values
  trait_values <- trait_values[!is.na(trait_values)]

  # Remove infinite values (can occur with logit of extreme proportions)
  trait_values <- trait_values[is.finite(trait_values)]

  cat("Fitting models for", length(trait_values), "species...\n\n")

  if (length(trait_values) < 10) {
    cat("  SKIPPING: Too few valid values\n")
    next
  }

  # Fit Brownian Motion
  cat("  Fitting BM...")
  fit_bm <- tryCatch(
    fitContinuous(tree_pruned_mcmc, trait_values, model = "BM"),
    error = function(e) { cat(" FAILED\n"); return(NULL) }
  )
  if (!is.null(fit_bm)) cat(" AICc =", round(fit_bm$opt$aicc, 2), "\n")

  # Fit Ornstein-Uhlenbeck
  cat("  Fitting OU...")
  fit_ou <- tryCatch(
    fitContinuous(tree_pruned_mcmc, trait_values, model = "OU"),
    error = function(e) { cat(" FAILED\n"); return(NULL) }
  )
  if (!is.null(fit_ou)) cat(" AICc =", round(fit_ou$opt$aicc, 2), "\n")

  # Fit Early Burst
  cat("  Fitting EB...")
  fit_eb <- tryCatch(
    fitContinuous(tree_pruned_mcmc, trait_values, model = "EB"),
    error = function(e) { cat(" FAILED\n"); return(NULL) }
  )
  if (!is.null(fit_eb)) cat(" AICc =", round(fit_eb$opt$aicc, 2), "\n")

  # Fit Pagel's Lambda
  cat("  Fitting Lambda...")
  fit_lambda <- tryCatch(
    fitContinuous(tree_pruned_mcmc, trait_values, model = "lambda"),
    error = function(e) { cat(" FAILED\n"); return(NULL) }
  )
  if (!is.null(fit_lambda)) cat(" AICc =", round(fit_lambda$opt$aicc, 2), "\n")

  # Compile results
  models_list <- list(BM = fit_bm, OU = fit_ou, EB = fit_eb,
                      Lambda = fit_lambda)
  models_list <- models_list[!sapply(models_list, is.null)]

  if (length(models_list) == 0) {
    cat("\n  WARNING: No models converged for this trait\n")
    next
  }

  # Extract AICc values
  aicc_values <- sapply(models_list, function(x) x$opt$aicc)

  # Calculate delta AICc and weights
  min_aicc <- min(aicc_values)
  delta_aicc <- aicc_values - min_aicc
  aicc_weights <- exp(-0.5 * delta_aicc) / sum(exp(-0.5 * delta_aicc))

  # Create results table
  results_df <- data.frame(
    Model = names(models_list),
    AICc = aicc_values,
    Delta_AICc = delta_aicc,
    Weight = aicc_weights,
    stringsAsFactors = FALSE
  )

  # Extract parameters
  results_df$Parameters <- sapply(names(models_list), function(m) {
    if (m == "BM") {
      return(paste("σ²=", round(models_list[[m]]$opt$sigsq, 4)))
    } else if (m == "OU") {
      return(paste("α=", round(models_list[[m]]$opt$alpha, 4),
                   ", σ²=", round(models_list[[m]]$opt$sigsq, 4)))
    } else if (m == "EB") {
      return(paste("r=", round(models_list[[m]]$opt$a, 4),
                   ", σ²=", round(models_list[[m]]$opt$sigsq, 4)))
    } else if (m == "Lambda") {
      return(paste("λ=", round(models_list[[m]]$opt$lambda, 4),
                   ", σ²=", round(models_list[[m]]$opt$sigsq, 4)))
    } else {
      return("NULL")
    }
  })

  # Sort by AICc
  results_df <- results_df[order(results_df$AICc), ]

  cat("\nModel Comparison:\n")
  print(results_df, row.names = FALSE)
  cat("\n")

  # Select best model with caution
  best_model <- results_df$Model[1]

  # Apply conservative selection rule: if ΔAICc < 2, prefer BM
  if ("BM" %in% results_df$Model && results_df$Delta_AICc[results_df$Model == "BM"] < 2) {
    selected_model <- "BM"
    cat("SELECTED: BM (ΔAICc < 2, prefer simpler model)\n")
  } else {
    selected_model <- best_model
    cat("SELECTED:", selected_model, "(clear best fit)\n")
  }

  # Store results
  all_results[[trait]] <- list(
    results_table = results_df,
    selected_model = selected_model,
    fits = models_list,
    n_species = length(trait_values)
  )
}

################################################################################
## 3. SUMMARY AND MODEL RECOMMENDATIONS
################################################################################

cat("\n\n", rep("=", 80), "\n", sep = "")
cat("STEP 3: MODEL SELECTION SUMMARY\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Selected models for each trait:\n")
for (trait in names(all_results)) {
  cat("  ", trait, ": ", all_results[[trait]]$selected_model,
      " (n = ", all_results[[trait]]$n_species, ")\n", sep = "")
}
cat("\n")

# Recommendation based on count traits (focal traits for PGLS)
count_traits <- c("log_balanced", "log_concordant", "log_reactions",
                  "log_giant_reactions", "log_acr")
count_traits <- count_traits[count_traits %in% names(all_results)]
selected_models <- sapply(all_results[count_traits], function(x) x$selected_model)

cat("Recommendation for PGLS correlation structure:\n")
cat("  Based on count traits:", paste(count_traits, collapse = ", "), "\n")

model_counts <- table(selected_models)
majority_model <- names(which.max(model_counts))

if (all(selected_models == "BM")) {
  cat("  → Use corBrownian() for all analyses\n")
  cat("  → Brownian motion is appropriate for all traits\n\n")
  recommended_cor <- "corBrownian"
} else if (all(selected_models == "Lambda")) {
  cat("  → Use corPagel() for all analyses\n")
  cat("  → Variable phylogenetic signal detected\n\n")
  recommended_cor <- "corPagel"
} else if (all(selected_models == "OU")) {
  cat("  → Use corMartins() for all analyses\n")
  cat("  → Stabilizing selection detected\n\n")
  recommended_cor <- "corMartins"
} else {
  cat("  → MIXED models selected\n")
  cat("  → Majority model:", majority_model, "(", model_counts[majority_model],
      "of", length(selected_models), "traits)\n")
  if (majority_model == "Lambda") {
    cat("  → Recommend corPagel() for all analyses\n\n")
    recommended_cor <- "corPagel"
  } else if (majority_model == "OU") {
    cat("  → Recommend corMartins() for all analyses\n\n")
    recommended_cor <- "corMartins"
  } else {
    cat("  → Recommend corPagel() as flexible default\n\n")
    recommended_cor <- "corPagel"
  }
}

################################################################################
## 4. SAVE OUTPUTS
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("STEP 4: SAVING OUTPUTS\n")
cat(rep("=", 80), "\n\n", sep = "")

# Save workspace
save(all_results, recommended_cor, selected_models, df_matched,
     file = "../results/workspaces/01_model_selection.RData")
cat("Workspace saved: ../results/workspaces/01_model_selection.RData\n")

# Save summary table
summary_table <- do.call(rbind, lapply(names(all_results), function(trait) {
  df <- all_results[[trait]]$results_table
  df$Trait <- trait
  df$Selected <- df$Model == all_results[[trait]]$selected_model
  df$N_species <- all_results[[trait]]$n_species
  return(df)
}))

write_csv(summary_table, "../results/tables/model_selection_summary.csv")
cat("Summary table saved: ../results/tables/model_selection_summary.csv\n")

# Create simple recommendation file
rec_df <- data.frame(
  Trait = names(all_results),
  Selected_Model = sapply(all_results, function(x) x$selected_model),
  N_species = sapply(all_results, function(x) x$n_species),
  Recommended_Correlation = recommended_cor
)
write_csv(rec_df, "../results/tables/model_recommendations.csv")
cat("Recommendations saved: ../results/tables/model_recommendations.csv\n\n")

# Update analysis log
log_update <- paste0(
  "\n## Script 01: Evolutionary Model Selection\n\n",
  "**Date**: ", Sys.Date(), "\n\n",
  "### Models Tested\n",
  "- Brownian Motion (BM)\n",
  "- Ornstein-Uhlenbeck (OU)\n",
  "- Early Burst (EB)\n",
  "- Pagel's Lambda\n\n",
  "Note: Drift/mean_trend excluded (requires non-ultrametric trees).\n\n",
  "### Count Traits Tested\n",
  paste0("- ", count_traits, ": ", selected_models, collapse = "\n"), "\n\n",
  "### Recommendation\n",
  "Use `", recommended_cor, "` as correlation structure for PGLS analyses.\n\n",
  "### Selection Rule\n",
  "- If ΔAICc < 2 between models → prefer simpler model (BM)\n",
  "- This protects against overfitting\n\n",
  "---\n\n"
)

cat(log_update, file = "../documentation/ANALYSIS_LOG.md", append = TRUE)
cat("Analysis log updated\n\n")

################################################################################
## 5. FINAL SUMMARY
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("MODEL SELECTION COMPLETE\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Summary:\n")
cat("  ✓ Tested 4 evolutionary models for", length(all_results), "traits\n")
cat("  ✓ Count traits:", length(count_traits), "\n")
cat("  ✓ Proportion traits:", length(all_results) - length(count_traits), "\n")
cat("  ✓ Selected best model for each trait\n")
cat("  ✓ Recommended correlation structure:", recommended_cor, "\n\n")

cat("Next steps:\n")
cat("  1. Run Script 02: PGLS main analysis\n")
cat("  2. Use", recommended_cor, "in all PGLS models\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat(rep("=", 80), "\n\n", sep = "")

sessionInfo()

sink()

