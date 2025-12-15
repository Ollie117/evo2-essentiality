#!/usr/bin/env Rscript
# Distribution plots for mutation size analysis
#
# Generates violin/boxplots showing delta log-likelihood distributions
# by essentiality class across mutation size levels.

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)

# ============================================================================
# CONFIGURATION
# ============================================================================

scores_file <- "../../../data/processed/essentiality_scores_mutation_size.csv"
tnseq_file <- "../../../data/processed/tnseq_clean.csv"
output_dir <- "../../../results/mutation_size"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# LOAD AND MERGE DATA
# ============================================================================

cat("Loading data...\n")
scores <- read_csv(scores_file, col_types = cols())
tnseq <- read_csv(tnseq_file, col_types = cols())

merged <- scores %>%
  inner_join(tnseq, by = "gene_id")

cat(sprintf("  Merged dataset: %d genes\n", nrow(merged)))
cat(sprintf("  Essential: %d\n", sum(merged$is_essential == 1)))
cat(sprintf("  Non-essential: %d\n", sum(merged$is_essential == 0)))

# ============================================================================
# RESHAPE DATA FOR PLOTTING
# ============================================================================

# Pivot to long format for faceted plotting
delta_ll_cols <- paste0("delta_ll_size", 1:5)
patterns <- c("TAA", "TAATAA", "TAATAATAA", "TAATAATAATAG", "TAATAATAATAGTGA")

long_data <- merged %>%
  select(gene_id, final_call_binary, all_of(delta_ll_cols)) %>%
  pivot_longer(
    cols = all_of(delta_ll_cols),
    names_to = "size_level",
    values_to = "delta_ll"
  ) %>%
  mutate(
    size = as.numeric(gsub("delta_ll_size", "", size_level)),
    pattern = factor(patterns[size], levels = patterns),
    label = factor(
      paste0("Size ", size, " (", patterns[size], ")"),
      levels = paste0("Size ", 1:5, " (", patterns, ")")
    )
  )

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\nDelta log-likelihood summary by size and essentiality:\n")
summary_stats <- long_data %>%
  group_by(size, pattern, final_call_binary) %>%
  summarise(
    n = n(),
    mean = mean(delta_ll, na.rm = TRUE),
    median = median(delta_ll, na.rm = TRUE),
    sd = sd(delta_ll, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_stats)

write_csv(summary_stats, file.path(output_dir, "delta_ll_summary_by_size_essentiality.csv"))

# ============================================================================
# PLOT 1: FACETED VIOLIN PLOT BY SIZE
# ============================================================================

cat("\nGenerating faceted violin plot...\n")

p_violin_facet <- ggplot(long_data, aes(x = final_call_binary, y = delta_ll, fill = final_call_binary)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.3, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Essential" = "#e34a33", "Non-Essential" = "#3182bd")) +
  facet_wrap(~label, nrow = 1) +
  labs(
    title = "Delta Log-Likelihood Distributions Across Mutation Sizes",
    x = "Essentiality Classification",
    y = "Delta Log-Likelihood (Δℓ)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(size = 8),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "delta_ll_violin_faceted_mutation_size.png"), p_violin_facet, 
       width = 14, height = 5, dpi = 300)
cat(sprintf("  Saved: %s/delta_ll_violin_faceted_mutation_size.png\n", output_dir))

# ============================================================================
# PLOT 2: BOXPLOT COMPARING SIZES (ESSENTIAL GENES ONLY)
# ============================================================================

cat("Generating boxplot for essential genes across sizes...\n")

essential_data <- long_data %>% filter(final_call_binary == "Essential")

p_box_essential <- ggplot(essential_data, aes(x = pattern, y = delta_ll, fill = pattern)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_brewer(palette = "OrRd") +
  labs(
    title = "Delta Log-Likelihood for Essential Genes Across Mutation Sizes",
    subtitle = "More negative values indicate stronger deleterious signal",
    x = "Stop Codon Pattern",
    y = "Delta Log-Likelihood (Δℓ)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "delta_ll_boxplot_essential_mutation_size.png"), p_box_essential, 
       width = 8, height = 6, dpi = 300)
cat(sprintf("  Saved: %s/delta_ll_boxplot_essential_mutation_size.png\n", output_dir))

# ============================================================================
# PLOT 3: DENSITY PLOT OVERLAY BY ESSENTIALITY (SIZE 1 vs SIZE 5)
# ============================================================================

cat("Generating density comparison (minimal vs baseline)...\n")

comparison_data <- long_data %>%
  filter(size %in% c(1, 5)) %>%
  mutate(size_label = ifelse(size == 1, "Minimal (TAA)", "Baseline (TAATAATAATAGTGA)"))

p_density_compare <- ggplot(comparison_data, aes(x = delta_ll, fill = final_call_binary)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Essential" = "#e34a33", "Non-Essential" = "#3182bd")) +
  facet_wrap(~size_label, ncol = 1) +
  labs(
    title = "Density of Delta Log-Likelihood: Minimal vs Baseline Perturbation",
    x = "Delta Log-Likelihood (Δℓ)",
    y = "Density",
    fill = "Classification"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "delta_ll_density_minimal_vs_baseline.png"), p_density_compare, 
       width = 9, height = 7, dpi = 300)
cat(sprintf("  Saved: %s/delta_ll_density_minimal_vs_baseline.png\n", output_dir))

# ============================================================================
# PLOT 4: SEPARATION BY SIZE (MEAN DELTA_LL)
# ============================================================================

cat("Generating separation plot...\n")

separation_data <- long_data %>%
  group_by(size, pattern, final_call_binary) %>%
  summarise(mean_delta_ll = mean(delta_ll, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = final_call_binary, values_from = mean_delta_ll) %>%
  mutate(separation = `Non-Essential` - Essential)

p_separation <- ggplot(separation_data, aes(x = pattern, y = separation)) +
  geom_col(fill = "#756bb1", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", separation)), vjust = -0.5, size = 3.5) +
  labs(
    title = "Mean Delta Log-Likelihood Separation by Mutation Size",
    subtitle = "Difference between non-essential and essential gene means",
    x = "Stop Codon Pattern",
    y = "Separation (Non-Essential - Essential)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "delta_ll_separation_mutation_size.png"), p_separation, 
       width = 8, height = 6, dpi = 300)
cat(sprintf("  Saved: %s/delta_ll_separation_mutation_size.png\n", output_dir))

cat("\n=== Distribution Plots Complete ===\n")