#!/usr/bin/env Rscript
# Publication-quality plots for mutation size analysis
#
# Style: Pastel colours, low-saturation, warm-cool contrast

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(patchwork)

# ============================================================================
# CONFIGURATION
# ============================================================================

scores_file <- "data/processed/essentiality_scores_mutation_size.csv"
tnseq_file <- "data/processed/tnseq_clean.csv"
data_dir <- "data/processed/mutation_size"
output_dir <- "results/mutation_size"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Mutation size levels
patterns <- c("TAA", "TAATAA", "TAATAATAA", "TAATAATAATAG", "TAATAATAATAGTGA")
pattern_lengths <- c("3 bp", "6 bp", "9 bp", "12 bp", "15 bp")

# ============================================================================
# COLOUR PALETTE - Pastel, low-saturation
# ============================================================================

colours <- list(
  essential = "#e08060",
  nonessential = "#6a9bc3",
  grey = "#8c8c8c"
)

# ROC curves - pastel, distinguishable
roc_colours <- c(
  "#e08060",  # pastel orange
  "#f0c070",  # pastel yellow/tan
  "#7ec8a0",  # pastel green
  "#6a9bc3",  # pastel blue
  "#c47a9a"   # pastel pink/mauve
)

# Bar chart colours - pastel
bar_colours <- list(
  auroc = "#9a7bb0",
  separation = "#7ab77a"
)

# ============================================================================
# CUSTOM THEME
# ============================================================================

theme_paper <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      strip.text = element_text(size = base_size, face = "bold"),
      panel.background = element_rect(fill = "white", color = "black", size = 0.8),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.5),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading data...\n")

# AUC results
auc_data <- read_csv(file.path(data_dir, "auc_comparison_mutation_size.csv"), col_types = cols())
auc_data$pattern_label <- factor(pattern_lengths, levels = pattern_lengths)

# ROC curve data
roc_data_list <- list()
for (i in 1:5) {
  roc_data_list[[i]] <- read_csv(file.path(data_dir, sprintf("roc_data_size%d.csv", i)), col_types = cols()) %>%
    mutate(size_level = i)
}
roc_data_all <- bind_rows(roc_data_list)
roc_data_all$pattern_label <- factor(
  pattern_lengths[roc_data_all$size_level],
  levels = pattern_lengths
)

# Scores and TnSeq for distribution plots
scores <- read_csv(scores_file, col_types = cols())
tnseq <- read_csv(tnseq_file, col_types = cols())
merged <- scores %>% inner_join(tnseq, by = "gene_id")

# Reshape to long format
delta_cols <- paste0("delta_ll_size", 1:5)
long_data <- merged %>%
  select(gene_id, final_call_binary, all_of(delta_cols)) %>%
  pivot_longer(
    cols = all_of(delta_cols),
    names_to = "size_level",
    values_to = "delta_ll"
  ) %>%
  mutate(
    level = as.numeric(gsub("delta_ll_size", "", size_level)),
    pattern = patterns[level],
    label = factor(pattern_lengths[level], levels = pattern_lengths)
  )

cat(sprintf("  Loaded %d genes\n", nrow(merged)))

# ============================================================================
# PLOT A: ROC CURVES - Pastel colours
# ============================================================================

cat("\nGenerating plots...\n")

p_roc <- ggplot(roc_data_all, aes(x = 1 - specificity, y = sensitivity, 
                                   color = pattern_label)) +
  geom_hline(yintercept = 0.5, color = "grey80", size = 0.3) +
  geom_vline(xintercept = 0.5, color = "grey80", size = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
              color = colours$grey, size = 0.5) +
  geom_line(size = 1, alpha = 0.85) +
  scale_color_manual(values = roc_colours, name = "Pattern length") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(
    x = "False positive rate (1 - Specificity)",
    y = "True positive rate (Sensitivity)"
  ) +
  theme_paper() +
  theme(
    legend.position = c(0.75, 0.25),
    legend.background = element_rect(fill = "white", color = "grey50", size = 0.3)
  ) +
  coord_equal()

ggsave(file.path(output_dir, "roc_curves_mutation_size.png"), p_roc,
       width = 6, height = 6, dpi = 300)
cat("  Saved roc_curves_mutation_size.png\n")

# ============================================================================
# PLOT B: BOXPLOT - Unfilled, pastel colours
# ============================================================================

set.seed(42)
jitter_data <- long_data %>%
  group_by(level, final_call_binary) %>%
  sample_frac(0.15) %>%
  ungroup()

p_boxplot <- ggplot(long_data, aes(x = label, y = delta_ll)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = colours$grey, size = 0.5) +
  geom_jitter(data = jitter_data,
              aes(color = final_call_binary),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              size = 1, alpha = 0.5) +
  geom_boxplot(aes(color = final_call_binary),
               fill = NA,
               outlier.shape = NA, 
               width = 0.5,
               position = position_dodge(width = 0.75),
               size = 0.6) +
  scale_color_manual(values = c("Essential" = colours$essential, 
                                 "Non-Essential" = colours$nonessential),
                     name = "Classification") +
  labs(
    x = "Stop codon pattern length",
    y = "Delta log-likelihood (Δℓ)"
  ) +
  theme_paper() +
  theme(legend.position = "right")

ggsave(file.path(output_dir, "delta_ll_boxplot_mutation_size.png"), p_boxplot,
       width = 10, height = 5, dpi = 300)
cat("  Saved delta_ll_boxplot_mutation_size.png\n")

# ============================================================================
# PLOT C: AUROC BAR PLOT - Pastel
# ============================================================================

p_auroc_bar <- ggplot(auc_data, aes(x = pattern_label, y = auc)) +
  geom_col(fill = bar_colours$auroc, width = 0.6, alpha = 0.8) +
  geom_errorbar(aes(ymin = auc_lower, ymax = auc_upper), 
                width = 0.15, size = 0.5, color = "black") +
  geom_text(aes(label = sprintf("%.3f", auc)), 
            vjust = -0.5, size = 2.8) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0.02)) +
  labs(
    x = "Pattern length",
    y = "AUROC"
  ) +
  theme_paper() +
  theme(axis.text.x = element_text(size = 9))

ggsave(file.path(output_dir, "auroc_barplot_mutation_size.png"), p_auroc_bar,
       width = 5, height = 4, dpi = 300)
cat("  Saved auroc_barplot_mutation_size.png\n")

# ============================================================================
# PLOT D: SEPARATION BAR PLOT - Pastel
# ============================================================================

separation_data <- long_data %>%
  group_by(level, label, final_call_binary) %>%
  summarise(
    mean_delta_ll = mean(delta_ll, na.rm = TRUE),
    sd_delta_ll = sd(delta_ll, na.rm = TRUE),
    n = n(),
    se = sd_delta_ll / sqrt(n),
    .groups = "drop"
  ) %>%
  select(level, label, final_call_binary, mean_delta_ll, se) %>%
  pivot_wider(names_from = final_call_binary, 
              values_from = c(mean_delta_ll, se)) %>%
  mutate(
    separation = `mean_delta_ll_Non-Essential` - `mean_delta_ll_Essential`,
    se_separation = sqrt(`se_Non-Essential`^2 + `se_Essential`^2),
    ci_lower = separation - 1.96 * se_separation,
    ci_upper = separation + 1.96 * se_separation
  )

separation_data$pattern_label <- factor(pattern_lengths, levels = pattern_lengths)

p_separation <- ggplot(separation_data, aes(x = pattern_label, y = separation)) +
  geom_col(fill = bar_colours$separation, width = 0.6, alpha = 0.8) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                width = 0.15, size = 0.5, color = "black") +
  geom_text(aes(label = sprintf("%.3f", separation)), 
            vjust = -0.5, size = 2.8) +
  labs(
    x = "Pattern length",
    y = "Mean Δℓ separation"
  ) +
  theme_paper() +
  theme(axis.text.x = element_text(size = 9))

ggsave(file.path(output_dir, "delta_ll_separation_mutation_size.png"), p_separation,
       width = 5, height = 4, dpi = 300)
cat("  Saved delta_ll_separation_mutation_size.png\n")

# ============================================================================
# PLOT E: DENSITY COMPARISON (MINIMAL vs BASELINE)
# ============================================================================

comparison_data <- long_data %>%
  filter(level %in% c(1, 5)) %>%
  mutate(size_label = factor(
    ifelse(level == 1, "Minimal (TAA, 3 bp)", "Baseline (15 bp)"),
    levels = c("Minimal (TAA, 3 bp)", "Baseline (15 bp)")
  ))

p_density <- ggplot(comparison_data, aes(x = delta_ll, fill = final_call_binary)) +
  geom_density(alpha = 0.5, color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", color = colours$grey, size = 0.5) +
  scale_fill_manual(values = c("Essential" = colours$essential, 
                                "Non-Essential" = colours$nonessential),
                    name = "Classification") +
  facet_wrap(~size_label, ncol = 1) +
  labs(
    x = "Delta log-likelihood (Δℓ)",
    y = "Density"
  ) +
  theme_paper() +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = NA, color = "black", size = 0.5)
  )

ggsave(file.path(output_dir, "delta_ll_density_comparison.png"), p_density,
       width = 7, height = 6, dpi = 300)
cat("  Saved delta_ll_density_comparison.png\n")

# ============================================================================
# COMBINED FIGURE: ROC top, boxplot middle, bar charts bottom
# ============================================================================

cat("\nGenerating combined figure...\n")

p_A <- p_roc + labs(tag = "A")
p_B <- p_boxplot + labs(tag = "B")
p_C <- p_auroc_bar + labs(tag = "C")
p_D <- p_separation + labs(tag = "D")

combined <- p_A / p_B / (p_C | p_D) +
  plot_layout(heights = c(1.5, 1, 0.8))

ggsave(file.path(output_dir, "figure_mutation_size_combined.png"), combined,
       width = 12, height = 14, dpi = 300)
cat("  Saved figure_mutation_size_combined.png\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n=== Summary ===\n")
cat("AUROC by pattern length:\n")
print(auc_data %>% select(pattern_label, auc, auc_lower, auc_upper))

cat("\n=== All plots saved to", output_dir, "===\n")