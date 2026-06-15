# ==============================================================================
# PGLS Analysis: Kinetic Properties, Model Size, and Phenotypes
# ==============================================================================
# Tests correlations between:
#   - Kinetic properties: giant module size (reactions), ACR metabolites
#   - Model size: total reactions
#   - Phenotypes: substrate usage, biomass yield
# Performed on full dataset and excluding species with defective kinetic models (<1%)
# ==============================================================================

library(tidyverse)
library(readxl)
library(ape)
library(caper)
library(ggrepel)

library(here)
REPO_ROOT   <- here::here()
DATA_DIR    <- file.path(REPO_ROOT, "data")
RESULTS_DIR <- file.path(REPO_ROOT, "results")
FIGURES_DIR <- file.path(REPO_ROOT, "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)
source(file.path(REPO_ROOT, "scripts/figures/theme.R"))

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading data...\n")

# Kinetic properties
kinetic <- read_csv(file.path(RESULTS_DIR, "kinetic/species_kinetic_summary.csv.gz")) %>%
    dplyr::select(species_id, n_acr, clade)

# Giant module data
giant_module <- read_csv(file.path(RESULTS_DIR, "kinetic/giant_module_reaction_summary.csv.gz")) %>%
    dplyr::select(species_id, n_giant_reactions) %>%
    left_join(
        read_csv(file.path(RESULTS_DIR, "kinetic/split_model_sizes.csv.gz")) %>%
            dplyr::select(species_id, n_total_split_reactions),
        by = "species_id"
    )

# Phenotype data
panel_a <- read_excel(file.path(DATA_DIR, "phenotypes/Panel_A.xlsx")) %>%
    rename(species_raw = `...1`, biomass_yield = `biomass yield`) %>%
    mutate(species_id = str_replace_all(species_raw, "'", "")) %>%
    dplyr::select(species_id, biomass_yield)

# Substrate usage (Biolog)
biolog <- read_tsv(file.path(DATA_DIR, "phenotypes/Biolog_Substrate.tsv"))
substrate_usage <- biolog %>%
    dplyr::select(5:ncol(.)) %>%
    pivot_longer(cols = everything(), names_to = "species_id", values_to = "growth") %>%
    filter(growth == "1") %>%
    group_by(species_id) %>%
    summarise(n_substrates = n())

# Phylogeny
tree <- read.tree(file.path(DATA_DIR, "phylogeny/332_2408OGs_timetree_mcmctree.nwk"))

# ==============================================================================
# 2. MERGE AND PREPARE DATA
# ==============================================================================

df <- kinetic %>%
    inner_join(giant_module, by = "species_id") %>%
    inner_join(panel_a, by = "species_id") %>%
    inner_join(substrate_usage, by = "species_id") %>%
    distinct(species_id, .keep_all = TRUE) %>%
    mutate(
        # Proportional measures
        giant_module_pct = (n_giant_reactions / n_total_split_reactions) * 100,
        giant_module_prop = n_giant_reactions / n_total_split_reactions,
        acr_prop = n_acr / n_total_split_reactions,
        # Natural-log-transformed network properties (consistent with rest
        # of analysis pipeline, which uses base-R log() = ln)
        log_total_split_reactions = log(n_total_split_reactions),
        log_giant_reactions       = log(n_giant_reactions),
        log_acr                   = log(pmax(n_acr, 1)),
        # Reactions whose substrate complex is NOT in the giant kinetic module
        # (complement of n_giant_reactions within n_total_split_reactions)
        n_non_giant_reactions     = n_total_split_reactions - n_giant_reactions,
        log_non_giant_reactions   = log(pmax(n_non_giant_reactions, 1))
    ) %>%
    as.data.frame()

# Match with phylogeny
row.names(df) <- df$species_id
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, df$species_id))
df <- df[tree_pruned$tip.label, ]

cat("\nTotal species:", nrow(df), "\n")

# ==============================================================================
# 3. IDENTIFY KINETIC OUTLIERS
# ==============================================================================

cat("\n==================================================\n")
cat("IDENTIFYING KINETIC OUTLIERS\n")
cat("==================================================\n\n")

cat("Giant module size (% reactions) - Summary:\n")
print(summary(df$giant_module_pct))

cat("\n\nSpecies with <1% reactions in giant module:\n")
outliers <- df %>%
    filter(giant_module_pct < 1) %>%
    arrange(giant_module_pct) %>%
    dplyr::select(species_id, clade, giant_module_pct, n_giant_reactions, n_total_split_reactions, n_acr)
print(outliers)

# Create datasets
df_full <- df
df_clean <- df %>% filter(giant_module_pct >= 1)

tree_full <- tree_pruned
tree_clean <- drop.tip(tree, setdiff(tree$tip.label, df_clean$species_id))
df_clean <- df_clean[tree_clean$tip.label, ]

cat("\n\nDataset sizes:")
cat("\n  Full dataset:  ", nrow(df_full), "species")
cat("\n  Clean dataset: ", nrow(df_clean), "species (excluding", nrow(df_full) - nrow(df_clean), "outliers)\n\n")

# ==============================================================================
# 4. PGLS ANALYSIS FUNCTION
# ==============================================================================

run_pgls <- function(formula_str, data, tree, dataset_name) {
    comp_data <- comparative.data(
        phy = tree, data = data, names.col = species_id,
        vcv = TRUE, na.omit = TRUE, warn.dropped = TRUE
    )

    model <- pgls(as.formula(formula_str), data = comp_data, lambda = "ML")

    # Extract results
    coef_table <- summary(model)$coefficients
    slope_val <- coef_table[2, 1]
    r2_val    <- summary(model)$r.squared
    # Signed correlation (Zoran's comment): r = sign(slope) * sqrt(R^2)
    signed_r  <- sign(slope_val) * sqrt(r2_val)
    results <- list(
        formula = formula_str,
        dataset = dataset_name,
        n = nrow(data),
        p_value = coef_table[2, 4],
        slope = slope_val,
        r_squared = r2_val,
        pearson_r = signed_r,
        adj_r_squared = summary(model)$adj.r.squared,
        lambda = as.numeric(summary(model)$param["lambda"])
    )

    return(results)
}

# ==============================================================================
# 5. RUN ALL ANALYSES
# ==============================================================================

cat("\n==================================================\n")
cat("RUNNING PGLS ANALYSES\n")
cat("==================================================\n\n")

# Define all model formulas
# Network properties are natural-log-transformed (consistent with the rest
# of the pipeline); biomass yield
# and substrate usage breadth remain on the linear scale (already bounded).
models <- list(
    # === SUBSTRATE USAGE CORRELATIONS ===
    # Kinetic properties (absolute, log)
    list(formula = "log_giant_reactions ~ n_substrates",
         label = "ln(Giant module reactions) ~ Substrate usage"),
    list(formula = "log_acr ~ n_substrates",
         label = "ln(ACR metabolites) ~ Substrate usage"),

    # Kinetic properties (proportional)
    list(formula = "giant_module_prop ~ n_substrates",
         label = "Giant module proportion ~ Substrate usage"),
    list(formula = "acr_prop ~ n_substrates",
         label = "ACR proportion ~ Substrate usage"),

    # Model size vs substrates (log)
    list(formula = "log_total_split_reactions ~ n_substrates",
         label = "ln(Total reactions) ~ Substrate usage"),

    # === BIOMASS YIELD CORRELATIONS ===
    # Kinetic properties (absolute, log)
    list(formula = "biomass_yield ~ log_giant_reactions",
         label = "Biomass yield ~ ln(Giant module reactions)"),
    list(formula = "biomass_yield ~ log_acr",
         label = "Biomass yield ~ ln(ACR metabolites)"),

    # Kinetic properties (proportional)
    list(formula = "biomass_yield ~ giant_module_prop",
         label = "Biomass yield ~ Giant module proportion"),
    list(formula = "biomass_yield ~ acr_prop",
         label = "Biomass yield ~ ACR proportion"),

    # Model size vs biomass yield (log)
    list(formula = "biomass_yield ~ log_total_split_reactions",
         label = "Biomass yield ~ ln(Total reactions)"),

    # === REACTIONS OUTSIDE THE GIANT MODULE (complement of giant) ===
    list(formula = "biomass_yield ~ log_non_giant_reactions",
         label = "Biomass yield ~ ln(Non-giant module reactions)"),
    list(formula = "log_non_giant_reactions ~ n_substrates",
         label = "ln(Non-giant module reactions) ~ Substrate usage")
)

# Run analyses on both datasets
all_results <- tibble()

for (model in models) {
    cat("\n--------------------------------------------------\n")
    cat(model$label, "\n")
    cat("--------------------------------------------------\n")

    # Full dataset
    result_full <- run_pgls(model$formula, df_full, tree_full, "Full")
    cat("Full dataset (n=", result_full$n, "): p = ",
        format(result_full$p_value, digits = 3, scientific = TRUE),
        ", R² = ", format(result_full$r_squared, digits = 3), "\n", sep = "")

    # Clean dataset
    result_clean <- run_pgls(model$formula, df_clean, tree_clean, "Clean")
    cat("Clean dataset (n=", result_clean$n, "): p = ",
        format(result_clean$p_value, digits = 3, scientific = TRUE),
        ", R² = ", format(result_clean$r_squared, digits = 3), "\n", sep = "")

    # Store results
    all_results <- bind_rows(
        all_results,
        tibble(
            comparison = model$label,
            dataset = "Full",
            n = result_full$n,
            p_value = result_full$p_value,
            slope = result_full$slope,
            r_squared = result_full$r_squared,
            pearson_r = result_full$pearson_r,
            lambda = result_full$lambda,
            adj_r_squared = result_full$adj_r_squared
        ),
        tibble(
            comparison = model$label,
            dataset = "Clean",
            n = result_clean$n,
            p_value = result_clean$p_value,
            slope = result_clean$slope,
            r_squared = result_clean$r_squared,
            pearson_r = result_clean$pearson_r,
            lambda = result_clean$lambda,
            adj_r_squared = result_clean$adj_r_squared
        )
    )
}

# ==============================================================================
# 6. SUMMARY TABLE
# ==============================================================================

cat("\n\n==================================================\n")
cat("SUMMARY TABLE\n")
cat("==================================================\n\n")

summary_table <- all_results %>%
    mutate(
        sig = case_when(
            p_value < 0.001 ~ "***",
            p_value < 0.01 ~ "**",
            p_value < 0.05 ~ "*",
            p_value < 0.1 ~ ".",
            TRUE ~ ""
        ),
        p_formatted = format(p_value, digits = 3, scientific = TRUE),
        r2_formatted = format(r_squared, digits = 3),
        r_formatted  = format(pearson_r, digits = 2),
        lambda_formatted = format(lambda, digits = 2)
    ) %>%
    dplyr::select(comparison, dataset, n, p_formatted, sig,
                  r_formatted, r2_formatted, lambda_formatted, slope)

print(summary_table, n = Inf)

# Save results
all_results_save <- all_results %>% dplyr::select(-adj_r_squared)  # Remove any problematic columns
write_csv(all_results_save, file.path(FIGURES_DIR, "pgls_results_kinetic_phenotype.csv"))
write_csv(summary_table, file.path(FIGURES_DIR, "pgls_results_summary.csv"))

cat("\n\nResults saved to:", FIGURES_DIR, "\n")

# ==============================================================================
# 7. KEY FINDINGS
# ==============================================================================

cat("\n\n==================================================\n")
cat("KEY FINDINGS\n")
cat("==================================================\n\n")

significant_results <- all_results %>%
    filter(p_value < 0.05) %>%
    arrange(p_value)

if (nrow(significant_results) > 0) {
    cat("Significant correlations (p < 0.05):\n\n")
    for (i in 1:nrow(significant_results)) {
        row <- significant_results[i, ]
        cat(sprintf("  %s [%s dataset, n=%d]\n    p = %.2g, R² = %.3f, slope = %.2g\n\n",
                    row$comparison, row$dataset, row$n,
                    row$p_value, row$r_squared, row$slope))
    }
} else {
    cat("No significant correlations found (all p >= 0.05).\n")
}

marginal_results <- all_results %>%
    filter(p_value >= 0.05 & p_value < 0.1) %>%
    arrange(p_value)

if (nrow(marginal_results) > 0) {
    cat("\nMarginal trends (0.05 ≤ p < 0.1):\n\n")
    for (i in 1:nrow(marginal_results)) {
        row <- marginal_results[i, ]
        cat(sprintf("  %s [%s dataset, n=%d]\n    p = %.2g, R² = %.3f, slope = %.2g\n\n",
                    row$comparison, row$dataset, row$n,
                    row$p_value, row$r_squared, row$slope))
    }
}

# ==============================================================================
# 7b. JOINT PGLS: PARTIAL EFFECTS OF SIZE AND GIANT MODULE
# ==============================================================================
# Reviewer asked whether giant module size predicts biomass yield
# independently of total reaction count. Tested with multiple-predictor PGLS
# (consistent with the univariate phylogenetic correction used elsewhere).

cat("\n==================================================\n")
cat("JOINT PGLS: SIZE + GIANT MODULE -> BIOMASS YIELD\n")
cat("==================================================\n\n")

comp_data_clean <- comparative.data(
    phy = tree_clean, data = df_clean, names.col = species_id,
    vcv = TRUE, na.omit = TRUE, warn.dropped = FALSE
)

# Collinearity diagnostic (VIF for 2 predictors = 1 / (1 - r^2))
r_pred  <- cor(df_clean$log_total_split_reactions, df_clean$log_giant_reactions)
vif_pred <- 1 / (1 - r_pred^2)
cat(sprintf("Predictor correlation (ln total reactions vs ln giant module): r = %.3f\n", r_pred))
cat(sprintf("VIF: %.2f", vif_pred))
if (vif_pred > 5) cat("  [>5: high collinearity, partial coefficients unstable]")
cat("\n\n")

# Nested models on the same dataset (log-transformed predictors)
m_size  <- pgls(biomass_yield ~ log_total_split_reactions,                         data = comp_data_clean, lambda = "ML")
m_giant <- pgls(biomass_yield ~ log_giant_reactions,                               data = comp_data_clean, lambda = "ML")
m_joint <- pgls(biomass_yield ~ log_total_split_reactions + log_giant_reactions,   data = comp_data_clean, lambda = "ML")

print_model <- function(model, label) {
    coefs <- summary(model)$coefficients
    cat(label, "\n")
    print(round(coefs, 6))
    cat(sprintf("  R^2 = %.3f, adj.R^2 = %.3f, AIC = %.2f, lambda = %.3f\n\n",
                summary(model)$r.squared, summary(model)$adj.r.squared,
                AIC(model), as.numeric(summary(model)$param["lambda"])))
}

print_model(m_size,  "biomass_yield ~ log(n_total_split_reactions)")
print_model(m_giant, "biomass_yield ~ log(n_giant_reactions)")
print_model(m_joint, "biomass_yield ~ log(n_total_split_reactions) + log(n_giant_reactions)")

# AIC comparison (lower is better; delta < 2 = indistinguishable)
aic_table <- data.frame(
    model      = c("size only", "giant only", "size + giant"),
    n_params   = c(2, 2, 3),
    AIC        = c(AIC(m_size), AIC(m_giant), AIC(m_joint)),
    R2         = c(summary(m_size)$r.squared,
                   summary(m_giant)$r.squared,
                   summary(m_joint)$r.squared)
)
aic_table$delta_AIC <- aic_table$AIC - min(aic_table$AIC)
cat("AIC comparison:\n")
print(aic_table, row.names = FALSE)
cat("\n")

# Likelihood ratio test: does adding giant module improve over size-only?
ll_size  <- as.numeric(logLik(m_size))
ll_joint <- as.numeric(logLik(m_joint))
lr_stat  <- 2 * (ll_joint - ll_size)
lr_p     <- pchisq(lr_stat, df = 1, lower.tail = FALSE)
cat(sprintf("LRT (size vs size+giant): chi^2 = %.3f, df = 1, p = %.4f\n\n", lr_stat, lr_p))

# Save joint model output
joint_coefs <- as.data.frame(summary(m_joint)$coefficients)
joint_coefs$predictor <- rownames(joint_coefs)
joint_coefs <- joint_coefs[, c("predictor", setdiff(colnames(joint_coefs), "predictor"))]
colnames(joint_coefs) <- c("predictor", "estimate", "std_error", "t_value", "p_value")
joint_coefs$vif <- c(NA, vif_pred, vif_pred)
write_csv(joint_coefs, file.path(FIGURES_DIR, "pgls_joint_biomass.csv"))
write_csv(aic_table,   file.path(FIGURES_DIR, "pgls_joint_aic_comparison.csv"))
cat("Saved: pgls_joint_biomass.csv, pgls_joint_aic_comparison.csv\n")

# ==============================================================================
# 7c. JOINT PGLS: SIZE + ACR -> BIOMASS YIELD
# ==============================================================================

cat("\n==================================================\n")
cat("JOINT PGLS: SIZE + ACR -> BIOMASS YIELD\n")
cat("==================================================\n\n")

r_pred_acr  <- cor(df_clean$log_total_split_reactions, df_clean$log_acr)
vif_pred_acr <- 1 / (1 - r_pred_acr^2)
cat(sprintf("Predictor correlation (ln total reactions vs ln ACR): r = %.3f\n", r_pred_acr))
cat(sprintf("VIF: %.2f", vif_pred_acr))
if (vif_pred_acr > 5) cat("  [>5: high collinearity]")
cat("\n\n")

m_acr       <- pgls(biomass_yield ~ log_acr,                                       data = comp_data_clean, lambda = "ML")
m_joint_acr <- pgls(biomass_yield ~ log_total_split_reactions + log_acr,           data = comp_data_clean, lambda = "ML")

print_model(m_acr,       "biomass_yield ~ log(n_acr)")
print_model(m_joint_acr, "biomass_yield ~ log(n_total_split_reactions) + log(n_acr)")

aic_table_acr <- data.frame(
    model    = c("size only", "ACR only", "size + ACR"),
    n_params = c(2, 2, 3),
    AIC      = c(AIC(m_size), AIC(m_acr), AIC(m_joint_acr)),
    R2       = c(summary(m_size)$r.squared,
                 summary(m_acr)$r.squared,
                 summary(m_joint_acr)$r.squared)
)
aic_table_acr$delta_AIC <- aic_table_acr$AIC - min(aic_table_acr$AIC)
cat("AIC comparison:\n")
print(aic_table_acr, row.names = FALSE)
cat("\n")

ll_joint_acr <- as.numeric(logLik(m_joint_acr))
lr_stat_acr <- 2 * (ll_joint_acr - ll_size)
lr_p_acr    <- pchisq(lr_stat_acr, df = 1, lower.tail = FALSE)
cat(sprintf("LRT (size vs size+ACR): chi^2 = %.3f, df = 1, p = %.4f\n\n", lr_stat_acr, lr_p_acr))

joint_coefs_acr <- as.data.frame(summary(m_joint_acr)$coefficients)
joint_coefs_acr$predictor <- rownames(joint_coefs_acr)
joint_coefs_acr <- joint_coefs_acr[, c("predictor", setdiff(colnames(joint_coefs_acr), "predictor"))]
colnames(joint_coefs_acr) <- c("predictor", "estimate", "std_error", "t_value", "p_value")
joint_coefs_acr$vif <- c(NA, vif_pred_acr, vif_pred_acr)
write_csv(joint_coefs_acr, file.path(FIGURES_DIR, "pgls_joint_biomass_acr.csv"))

# ==============================================================================
# 7d. JOINT PGLS: SIZE + GIANT MODULE -> SUBSTRATE USAGE
# ==============================================================================

cat("\n==================================================\n")
cat("JOINT PGLS: SIZE + GIANT MODULE -> SUBSTRATE USAGE\n")
cat("==================================================\n\n")
cat("Univariate models all p > 0.3; joint model checks for suppression.\n\n")

m_sub_size  <- pgls(n_substrates ~ log_total_split_reactions,                         data = comp_data_clean, lambda = "ML")
m_sub_giant <- pgls(n_substrates ~ log_giant_reactions,                               data = comp_data_clean, lambda = "ML")
m_sub_joint <- pgls(n_substrates ~ log_total_split_reactions + log_giant_reactions,   data = comp_data_clean, lambda = "ML")

print_model(m_sub_size,  "n_substrates ~ log(n_total_split_reactions)")
print_model(m_sub_giant, "n_substrates ~ log(n_giant_reactions)")
print_model(m_sub_joint, "n_substrates ~ log(n_total_split_reactions) + log(n_giant_reactions)")

aic_table_sub <- data.frame(
    model    = c("size only", "giant only", "size + giant"),
    n_params = c(2, 2, 3),
    AIC      = c(AIC(m_sub_size), AIC(m_sub_giant), AIC(m_sub_joint)),
    R2       = c(summary(m_sub_size)$r.squared,
                 summary(m_sub_giant)$r.squared,
                 summary(m_sub_joint)$r.squared)
)
aic_table_sub$delta_AIC <- aic_table_sub$AIC - min(aic_table_sub$AIC)
cat("AIC comparison:\n")
print(aic_table_sub, row.names = FALSE)
cat("\n")

ll_sub_size  <- as.numeric(logLik(m_sub_size))
ll_sub_joint <- as.numeric(logLik(m_sub_joint))
lr_stat_sub  <- 2 * (ll_sub_joint - ll_sub_size)
lr_p_sub     <- pchisq(lr_stat_sub, df = 1, lower.tail = FALSE)
cat(sprintf("LRT (size vs size+giant for substrate usage): chi^2 = %.3f, df = 1, p = %.4f\n\n",
            lr_stat_sub, lr_p_sub))

joint_coefs_sub <- as.data.frame(summary(m_sub_joint)$coefficients)
joint_coefs_sub$predictor <- rownames(joint_coefs_sub)
joint_coefs_sub <- joint_coefs_sub[, c("predictor", setdiff(colnames(joint_coefs_sub), "predictor"))]
colnames(joint_coefs_sub) <- c("predictor", "estimate", "std_error", "t_value", "p_value")
write_csv(joint_coefs_sub, file.path(FIGURES_DIR, "pgls_joint_substrate.csv"))

cat("Saved: pgls_joint_biomass_acr.csv, pgls_joint_substrate.csv\n")

# ==============================================================================
# 7e. JOINT PGLS: GIANT + NON-GIANT REACTIONS -> BIOMASS YIELD
# ==============================================================================
# Decompose total reactions into "in giant kinetic module" and "outside it".
# If biomass yield correlates with reactions outside the giant module too,
# the apparent giant-module effect just reflects general network-size effects.

cat("\n==================================================\n")
cat("JOINT PGLS: GIANT + NON-GIANT REACTIONS -> BIOMASS YIELD\n")
cat("==================================================\n\n")

r_pred_ng <- cor(df_clean$log_giant_reactions, df_clean$log_non_giant_reactions)
vif_pred_ng <- 1 / (1 - r_pred_ng^2)
cat(sprintf("Predictor correlation (ln giant vs ln non-giant reactions): r = %.3f\n", r_pred_ng))
cat(sprintf("VIF: %.2f", vif_pred_ng))
if (vif_pred_ng > 5) cat("  [>5: high collinearity]")
cat("\n\n")

m_nongiant   <- pgls(biomass_yield ~ log_non_giant_reactions,
                      data = comp_data_clean, lambda = "ML")
m_joint_ng   <- pgls(biomass_yield ~ log_giant_reactions + log_non_giant_reactions,
                      data = comp_data_clean, lambda = "ML")

print_model(m_nongiant, "biomass_yield ~ ln(non_giant_reactions)")
print_model(m_joint_ng, "biomass_yield ~ ln(giant_reactions) + ln(non_giant_reactions)")

aic_table_ng <- data.frame(
    model    = c("giant only", "non-giant only", "giant + non-giant"),
    n_params = c(2, 2, 3),
    AIC      = c(AIC(m_giant), AIC(m_nongiant), AIC(m_joint_ng)),
    R2       = c(summary(m_giant)$r.squared,
                 summary(m_nongiant)$r.squared,
                 summary(m_joint_ng)$r.squared)
)
aic_table_ng$delta_AIC <- aic_table_ng$AIC - min(aic_table_ng$AIC)
cat("AIC comparison:\n")
print(aic_table_ng, row.names = FALSE)
cat("\n")

ll_giant_only   <- as.numeric(logLik(m_giant))
ll_nongiant     <- as.numeric(logLik(m_nongiant))
ll_joint_ng     <- as.numeric(logLik(m_joint_ng))
lr_g_to_joint   <- 2 * (ll_joint_ng - ll_giant_only)
lr_ng_to_joint  <- 2 * (ll_joint_ng - ll_nongiant)
lr_p_g_to_j     <- pchisq(lr_g_to_joint,  df = 1, lower.tail = FALSE)
lr_p_ng_to_j    <- pchisq(lr_ng_to_joint, df = 1, lower.tail = FALSE)
cat(sprintf("LRT (giant-only -> giant+non-giant): chi^2 = %.3f, p = %.4f\n",
            lr_g_to_joint, lr_p_g_to_j))
cat(sprintf("LRT (non-giant-only -> giant+non-giant): chi^2 = %.3f, p = %.4f\n",
            lr_ng_to_joint, lr_p_ng_to_j))

joint_coefs_ng <- as.data.frame(summary(m_joint_ng)$coefficients)
joint_coefs_ng$predictor <- rownames(joint_coefs_ng)
joint_coefs_ng <- joint_coefs_ng[, c("predictor", setdiff(colnames(joint_coefs_ng), "predictor"))]
colnames(joint_coefs_ng) <- c("predictor", "estimate", "std_error", "t_value", "p_value")
joint_coefs_ng$vif <- c(NA, vif_pred_ng, vif_pred_ng)
write_csv(joint_coefs_ng, file.path(FIGURES_DIR, "pgls_joint_biomass_giant_nongiant.csv"))
write_csv(aic_table_ng,   file.path(FIGURES_DIR, "pgls_joint_giant_nongiant_aic.csv"))
cat("Saved: pgls_joint_biomass_giant_nongiant.csv, pgls_joint_giant_nongiant_aic.csv\n")

cat("\n==================================================\n")
cat("GENERATING FIGURE FOR THESIS\n")
cat("==================================================\n\n")

# Create combined figure for thesis
library(patchwork)

# Style settings (sourced from thesis_theme.R)
theme_style <- theme_thesis()

# Panel A: Biomass yield vs ln(Total reactions)
comp_data_a <- comparative.data(phy = tree_clean, data = df_clean,
                                names.col = species_id, vcv = TRUE,
                                na.omit = TRUE, warn.dropped = FALSE)
model_a <- pgls(biomass_yield ~ log_total_split_reactions, data = comp_data_a, lambda = "ML")

annot_a <- sprintf("PGLS: r = %.2f, p = %.1e",
                   sign(coef(model_a)[2]) * sqrt(summary(model_a)$r.squared),
                   summary(model_a)$coefficients[2, 4])

p_a <- ggplot(df_clean, aes(x = log_total_split_reactions, y = biomass_yield)) +
    geom_point(alpha = 0.5, size = 2.5, color = CB_BLACK) +
    geom_abline(intercept = coef(model_a)[1], slope = coef(model_a)[2],
                color = COLOR_PGLS, linewidth = 1.2) +
    annotate("text", x = Inf, y = -Inf,
             label = annot_a,
             hjust = 1.05, vjust = -0.4, size = THESIS_ANNOT_SIZE,
             fontface = "italic") +
    labs(x = "ln(Total reactions)",
         y = "Biomass yield", title = "A") +
    theme_style

# Panel B: Biomass yield vs ln(Giant module reactions)
model_b <- pgls(biomass_yield ~ log_giant_reactions, data = comp_data_a, lambda = "ML")

annot_b <- sprintf("PGLS: r = %.2f, p = %.1e",
                   sign(coef(model_b)[2]) * sqrt(summary(model_b)$r.squared),
                   summary(model_b)$coefficients[2, 4])

p_b <- ggplot(df_clean, aes(x = log_giant_reactions, y = biomass_yield)) +
    geom_point(alpha = 0.5, size = 2.5, color = CB_BLACK) +
    geom_abline(intercept = coef(model_b)[1], slope = coef(model_b)[2],
                color = COLOR_PGLS, linewidth = 1.2) +
    annotate("text", x = Inf, y = -Inf,
             label = annot_b,
             hjust = 1.05, vjust = -0.4, size = THESIS_ANNOT_SIZE,
             fontface = "italic") +
    labs(x = "ln(Giant module reactions)",
         y = NULL, title = "B") +
    theme_style +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Combine panels
p_combined <- p_a + p_b

ggsave(file.path(FIGURES_DIR, "kinetic_phenotype_combined.pdf"),
       plot = p_combined, width = 10, height = 4.5, device = cairo_pdf)
ggsave(file.path(FIGURES_DIR, "kinetic_phenotype_combined.svg"),
       plot = p_combined, width = 10, height = 4.5, device = svg)
ggsave(file.path(FIGURES_DIR, "kinetic_phenotype_combined.png"),
       plot = p_combined, width = 10, height = 4.5, dpi = 300, bg = "white")

cat("Figures saved: kinetic_phenotype_combined.{pdf,svg,png}\n")

cat("\n==================================================\n")
cat("Analysis complete!\n")
cat("==================================================\n")

