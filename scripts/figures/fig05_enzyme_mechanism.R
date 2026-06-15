#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  require(ggplot2)
  require(dplyr)
  require(readr)
  require(tidyr)
})

# Get script directory
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else getwd()
library(here)
REPO_ROOT   <- here::here()
RESULTS_DIR <- file.path(REPO_ROOT, "results")
FIGURES_DIR <- file.path(REPO_ROOT, "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)
source(file.path(REPO_ROOT, "scripts/figures/theme.R"))

# Input file paths
concordance_csv <- file.path(RESULTS_DIR, 'concordance/concordance_stats_all_variants.csv')
metrics_csv <- file.path(RESULTS_DIR, 'concordance/enzyme_mechanism_metrics.csv')

concordance_csv <- normalizePath(concordance_csv, mustWork = FALSE)
metrics_csv <- normalizePath(metrics_csv, mustWork = FALSE)

if(!file.exists(concordance_csv)) {
  stop(paste('Concordance CSV not found:', concordance_csv))
}

if(!file.exists(metrics_csv)) {
  stop(paste('Metrics CSV not found:', metrics_csv))
}

# Read data
concordance_df <- read_csv(concordance_csv, show_col_types = FALSE)
metrics_df <- read_csv(metrics_csv, show_col_types = FALSE)

# Filter for the variants we want
variants_of_interest <- c("random_0", "random_25", "random_50", "random_75", "random_100")
concordance_df <- concordance_df %>% filter(variant %in% variants_of_interest)
metrics_df <- metrics_df %>% filter(variant %in% variants_of_interest)

# Calculate percentages for metrics INCLUDING ACR as proportion
metrics_df <- metrics_df %>%
  mutate(
    pct_balanced = (balanced / n_complexes) * 100,
    pct_giant = (largest_module_size / n_complexes) * 100,
    pct_acr = (acr_metabolites / n_metabolites) * 100
  )

# Merge the datasets
combined_df <- concordance_df %>%
  inner_join(
    metrics_df %>% select(seed, model, variant, n_metabolites, acr_metabolites, largest_module_size, balanced, n_complexes, pct_balanced, pct_giant, pct_acr),
    by = c("seed", "model", "variant"),
    suffix = c("", "_met")
  )

# Hierarchical approach: First average across seeds within each model
model_means <- combined_df %>%
  group_by(variant, model) %>%
  summarize(
    mean_concordance = mean(pct_concordance, na.rm = TRUE),
    mean_balanced = mean(pct_balanced, na.rm = TRUE),
    mean_giant = mean(pct_giant, na.rm = TRUE),
    mean_acr = mean(pct_acr, na.rm = TRUE),
    .groups = 'drop'
  )


# Then calculate summary statistics across models (n=13)
# Using base R aggregate to avoid dplyr issues
summary_by_variant_list <- list(
  variant = sort(unique(model_means$variant)),
  n_models = numeric(length(unique(model_means$variant))),
  mean_concordance = numeric(length(unique(model_means$variant))),
  sd_concordance = numeric(length(unique(model_means$variant))),
  se_concordance = numeric(length(unique(model_means$variant))),
  mean_balanced = numeric(length(unique(model_means$variant))),
  sd_balanced = numeric(length(unique(model_means$variant))),
  se_balanced = numeric(length(unique(model_means$variant))),
  mean_giant = numeric(length(unique(model_means$variant))),
  sd_giant = numeric(length(unique(model_means$variant))),
  se_giant = numeric(length(unique(model_means$variant))),
  mean_acr = numeric(length(unique(model_means$variant))),
  sd_acr = numeric(length(unique(model_means$variant))),
  se_acr = numeric(length(unique(model_means$variant)))
)

for (i in seq_along(summary_by_variant_list$variant)) {
  v <- summary_by_variant_list$variant[i]
  vdata <- model_means[model_means$variant == v, ]

  summary_by_variant_list$n_models[i] <- nrow(vdata)
  summary_by_variant_list$mean_concordance[i] <- mean(vdata$mean_concordance)
  summary_by_variant_list$sd_concordance[i] <- sd(vdata$mean_concordance)
  summary_by_variant_list$se_concordance[i] <- sd(vdata$mean_concordance) / sqrt(nrow(vdata))
  summary_by_variant_list$mean_balanced[i] <- mean(vdata$mean_balanced)
  summary_by_variant_list$sd_balanced[i] <- sd(vdata$mean_balanced)
  summary_by_variant_list$se_balanced[i] <- sd(vdata$mean_balanced) / sqrt(nrow(vdata))
  summary_by_variant_list$mean_giant[i] <- mean(vdata$mean_giant)
  summary_by_variant_list$sd_giant[i] <- sd(vdata$mean_giant)
  summary_by_variant_list$se_giant[i] <- sd(vdata$mean_giant) / sqrt(nrow(vdata))
  summary_by_variant_list$mean_acr[i] <- mean(vdata$mean_acr)
  summary_by_variant_list$sd_acr[i] <- sd(vdata$mean_acr)
  summary_by_variant_list$se_acr[i] <- sd(vdata$mean_acr) / sqrt(nrow(vdata))
}

summary_by_variant <- as.data.frame(summary_by_variant_list)

# Convert variant to factor with proper ordering
summary_by_variant$variant <- factor(summary_by_variant$variant,
                                      levels = c("random_0", "random_25", "random_50", "random_75", "random_100"))

# Prepare data for plotting - ALL 4 METRICS AS PROPORTIONS
plot_data <- summary_by_variant %>%
  select(variant, mean_concordance, se_concordance, mean_balanced, se_balanced, mean_giant, se_giant, mean_acr, se_acr) %>%
  pivot_longer(
    cols = c(mean_concordance, mean_balanced, mean_giant, mean_acr),
    names_to = "metric",
    values_to = "mean",
    names_prefix = "mean_"
  ) %>%
  mutate(
    se = case_when(
      metric == "concordance" ~ summary_by_variant$se_concordance[match(variant, summary_by_variant$variant)],
      metric == "balanced" ~ summary_by_variant$se_balanced[match(variant, summary_by_variant$variant)],
      metric == "giant" ~ summary_by_variant$se_giant[match(variant, summary_by_variant$variant)],
      metric == "acr" ~ summary_by_variant$se_acr[match(variant, summary_by_variant$variant)]
    )
  )

# Set metric labels
plot_data$metric <- factor(plot_data$metric,
                            levels = c("concordance", "balanced", "giant", "acr"),
                            labels = c("Complexes in concordance",
                                      "Balanced",
                                      "Complexes in giant module",
                                      "ACR metabolites"))

# Colorblind-friendly palette (Wong/Okabe-Ito — from thesis_theme.R)
colors <- c("Complexes in concordance" = CB_BLUE,
            "Balanced"                 = CB_ORANGE,
            "Complexes in giant module"= CB_TEAL,
            "ACR metabolites"          = CB_PINK)

# Create BAR CHART - all metrics as proportions
p <- ggplot(plot_data, aes(x = variant, y = mean, fill = metric)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8, color = "black", linewidth = 0.25) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.85), width = 0.3, linewidth = 0.4) +
  scale_fill_manual(values = colors, name = NULL) +
  scale_y_continuous(
    name = "Proportion (%)",
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    labels = c("random_0" = "0%", "random_25" = "25%", "random_50" = "50%",
               "random_75" = "75%", "random_100" = "100%")
  ) +
  labs(x = "Random binding") +
  theme_thesis(base_size = 20) +
  theme(
    legend.position   = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    axis.text.x       = element_text(face = "plain"),
    panel.grid.major.x = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE))

# Save plot
output_file <- file.path(FIGURES_DIR, "enzyme_mechanism_concordance.png")

ggsave(output_file,                              plot = p, width = 10, height = 6, dpi = 300, bg = "white")
ggsave(sub("\\.png$", ".pdf", output_file),      plot = p, width = 10, height = 6, device = cairo_pdf)
ggsave(sub("\\.png$", ".svg", output_file),      plot = p, width = 10, height = 6, device = svg)
cat('Wrote', output_file, '\n')

# Print summary statistics
cat('\nSummary statistics by variant:\n')
for (v in levels(summary_by_variant$variant)) {
  row <- summary_by_variant %>% filter(variant == v)
  cat(sprintf('\n%s:\n', v))
  cat(sprintf('  Complexes in concordance: %.2f%% ± %.2f%%\n',
              row$mean_concordance, row$se_concordance))
  cat(sprintf('  Balanced complexes: %.2f%% ± %.2f%%\n',
              row$mean_balanced, row$se_balanced))
  cat(sprintf('  Complexes in giant module: %.2f%% ± %.2f%%\n',
              row$mean_giant, row$se_giant))
  cat(sprintf('  ACR metabolites: %.2f%% ± %.2f%%\n',
              row$mean_acr, row$se_acr))
}

# Save combined data
combined_output <- file.path(RESULTS_DIR, "concordance/concordance_comparison_all_variants.csv")
write_csv(summary_by_variant, combined_output)
cat('\nSaved summary data to:', combined_output, '\n')
