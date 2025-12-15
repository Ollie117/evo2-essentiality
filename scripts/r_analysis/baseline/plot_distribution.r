#!/usr/bin/env Rscript
# Distribution plot for baseline validation
#
# Generates violin plot with boxplot overlay showing delta log-likelihood
# distributions for essential vs non-essential genes.

library(ggplot2)
library(readr)
library(dplyr)

# ============================================================================
# CONFIGURATION
# ============================================================================

input_file <- "../../../data/processed/merged_data.csv"
output_dir <- "../../../results/validation"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading merged data...\n")
data <- read_csv(input_file, col_types = cols())
cat(sprintf("  Loaded %d genes\n", nrow(data)))
cat(sprintf("  Essential: %d\n", sum(data$is_essential == 1)))
cat(sprintf("  Non-essential: %d\n", sum(data$is_essential == 0)))

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\nDelta log-likelihood summary by essentiality:\n")
summary_stats <- data %>%
  group_by(final_call_binary) %>%
  summarise(
    n = n(),
    mean = mean(delta_ll, na.rm = TRUE),
    median = median(delta_ll, na.rm = TRUE),
    sd = sd(delta_ll, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_stats)

# Save summary stats
write_csv(summary_stats, file.path(output_dir, "delta_ll_summary_by_essentiality.csv"))

# ============================================================================
# VIOLIN PLOT WITH BOXPLOT
# ============================================================================

cat("\nGenerating violin plot...\n")

p_violin <- ggplot(data, aes(x = final_call_binary, y = delta_ll, fill = final_call_binary)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Essential" = "#e34a33", "Non-Essential" = "#3182bd")) +
  labs(
    title = "Distribution of Delta Log-Likelihood by Gene Essentiality",
    subtitle = "More negative values indicate greater predicted deleteriousness",
    x = "Essentiality Classification (TnSeq)",
    y = "Delta Log-Likelihood (Δℓ)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "delta_ll_violin_by_essentiality.png"), p_violin, 
       width = 7, height = 6, dpi = 300)
cat(sprintf("  Saved: %s/delta_ll_violin_by_essentiality.png\n", output_dir))

# ============================================================================
# BOXPLOT ONLY (ALTERNATIVE)
# ============================================================================

cat("Generating boxplot...\n")

p_box <- ggplot(data, aes(x = final_call_binary, y = delta_ll, fill = final_call_binary)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Essential" = "#e34a33", "Non-Essential" = "#3182bd")) +
  labs(
    title = "Delta Log-Likelihood by Gene Essentiality",
    x = "Essentiality Classification (TnSeq)",
    y = "Delta Log-Likelihood (Δℓ)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "delta_ll_boxplot_by_essentiality.png"), p_box, 
       width = 7, height = 6, dpi = 300)
cat(sprintf("  Saved: %s/delta_ll_boxplot_by_essentiality.png\n", output_dir))

# ============================================================================
# DENSITY PLOT (ALTERNATIVE)
# ============================================================================

cat("Generating density plot...\n")

p_density <- ggplot(data, aes(x = delta_ll, fill = final_call_binary)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Essential" = "#e34a33", "Non-Essential" = "#3182bd")) +
  labs(
    title = "Density of Delta Log-Likelihood by Gene Essentiality",
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

ggsave(file.path(output_dir, "delta_ll_density_by_essentiality.png"), p_density, 
       width = 8, height = 6, dpi = 300)
cat(sprintf("  Saved: %s/delta_ll_density_by_essentiality.png\n", output_dir))

cat("\n=== Distribution Plots Complete ===\n")