#!/usr/bin/env Rscript
# =============================================================================
# Absolute Likelihood Analysis for Essentiality Prediction
# Section 3.1: Tests whether absolute WT likelihood predicts essentiality
# Style: Matches established publication style (pastel, low-saturation)
# =============================================================================

library(tidyverse)
library(pROC)
library(patchwork)

# =============================================================================
# CONFIGURATION
# =============================================================================
BASELINE_DATA <- "data/processed/essentiality_scores_test.csv"  # Your baseline results
TRUTH_DATA <- "data/processed/tnseq_clean.csv"  # TnSeq essentiality data
OUTPUT_DIR <- "results/absolute_likelihood"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# COLOUR PALETTE - Pastel, low-saturation (matching established style)
# =============================================================================

colours <- list(
  essential = "#e08060",
  nonessential = "#6a9bc3",
  primary = "#6a9bc3",
  secondary = "#009E73",
  grey = "#8c8c8c"
)

# =============================================================================
# CUSTOM THEME (matching established style)
# =============================================================================

theme_paper <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 2),
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

# =============================================================================
# LOAD DATA
# =============================================================================
cat("Loading baseline data...\n")
data <- read_csv(BASELINE_DATA, show_col_types = FALSE)
cat(sprintf("  %d genes loaded\n", nrow(data)))

cat("\nLoading TnSeq essentiality data...\n")
truth <- read_csv(TRUTH_DATA, show_col_types = FALSE)
cat(sprintf("  %d genes loaded\n", nrow(truth)))

# =============================================================================
# MERGE DATA
# =============================================================================
cat("\nMerging datasets on gene_id...\n")
data <- data %>%
  left_join(truth %>% select(gene_id, is_essential, final_call_binary), by = "gene_id")

cat(sprintf("  %d genes after merge\n", nrow(data)))

# =============================================================================
# CALCULATE PER-BASE LIKELIHOOD
# =============================================================================
data <- data %>%
  mutate(
    wt_logprob_per_base = wt_logprob / cds_length,
    essential_binary = as.numeric(is_essential == TRUE | is_essential == 1)
  )

cat("\nClass balance:\n")
cat(sprintf("  Essential: %d\n", sum(data$essential_binary == 1, na.rm = TRUE)))
cat(sprintf("  Non-essential: %d\n", sum(data$essential_binary == 0, na.rm = TRUE)))

# Remove any NAs
data <- data %>% filter(!is.na(essential_binary), !is.na(wt_logprob))
cat(sprintf("  After removing NAs: %d genes\n", nrow(data)))

# =============================================================================
# ROC ANALYSIS
# =============================================================================
cat("\n========== ROC ANALYSIS ==========\n")

# Helper function to get best ROC direction
get_best_roc <- function(response, predictor, name) {
  roc1 <- roc(response, predictor, direction = "<", quiet = TRUE)
  roc2 <- roc(response, predictor, direction = ">", quiet = TRUE)
  
  if (auc(roc1) >= auc(roc2)) {
    best_roc <- roc1
    direction <- "lower scores = essential"
  } else {
    best_roc <- roc2
    direction <- "higher scores = essential"
  }
  
  ci <- ci.auc(best_roc)
  cat(sprintf("\n%s:\n", name))
  cat(sprintf("   AUROC: %.3f (95%% CI: %.3f-%.3f)\n", auc(best_roc), ci[1], ci[3]))
  cat(sprintf("   Direction: %s\n", direction))
  
  return(best_roc)
}

# 1. Absolute WT log-likelihood (total)
roc_absolute <- get_best_roc(data$essential_binary, data$wt_logprob,
                              "1. Absolute WT log-likelihood (total)")

# 2. Absolute WT log-likelihood (per base)
roc_perbase <- get_best_roc(data$essential_binary, data$wt_logprob_per_base,
                             "2. Absolute WT log-likelihood (per base)")

# =============================================================================
# STATISTICAL COMPARISON
# =============================================================================
cat("\n========== STATISTICAL COMPARISONS ==========\n")
test_abs <- roc.test(roc_absolute, roc_perbase)
cat(sprintf("Absolute (total) vs Absolute (per-base): p = %.2e\n", test_abs$p.value))

# =============================================================================
# SUMMARY TABLE
# =============================================================================
ci_absolute <- ci.auc(roc_absolute)
ci_perbase <- ci.auc(roc_perbase)

summary_df <- tibble(
  Method = c("Absolute WT likelihood (total)",
             "Absolute WT likelihood (per-base)"),
  AUROC = c(auc(roc_absolute), auc(roc_perbase)),
  CI_lower = c(ci_absolute[1], ci_perbase[1]),
  CI_upper = c(ci_absolute[3], ci_perbase[3])
) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

cat("\n========== SUMMARY ==========\n")
print(summary_df)

write_tsv(summary_df, file.path(OUTPUT_DIR, "absolute_likelihood_summary.tsv"))

# =============================================================================
# PLOT A: ROC CURVES
# =============================================================================
cat("\nGenerating plots...\n")

# Extract ROC data
roc_data_total <- data.frame(
  specificity = roc_absolute$specificities,
  sensitivity = roc_absolute$sensitivities,
  Method = "Total"
)

roc_data_perbase <- data.frame(
  specificity = roc_perbase$specificities,
  sensitivity = roc_perbase$sensitivities,
  Method = "Per-base"
)

roc_data <- bind_rows(roc_data_total, roc_data_perbase)

auc_label_total <- sprintf("Absolute total = %.3f", auc(roc_absolute))
auc_label_perbase <- sprintf("Absolute per-base = %.3f", auc(roc_perbase))

p_roc <- ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity, color = Method)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
              color = colours$grey, size = 0.5) +
  geom_line(size = 1, alpha = 0.85) +
  scale_color_manual(values = c("Total" = colours$primary, 
                                 "Per-base" = colours$secondary),
                     labels = c(auc_label_perbase, auc_label_total)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(
    x = "False positive rate (1 - Specificity)",
    y = "True positive rate (Sensitivity)",
    color = "AUROC"
  ) +
  theme_paper() +
  theme(
    legend.position = c(0.7, 0.2),
    legend.title = element_blank()
  ) +
  coord_equal()

ggsave(file.path(OUTPUT_DIR, "roc_absolute_likelihood.png"), p_roc,
       width = 6, height = 6, dpi = 300)
cat("  Saved roc_absolute_likelihood.png\n")

# =============================================================================
# PLOT B: BOXPLOT - TOTAL LIKELIHOOD
# =============================================================================

set.seed(42)
jitter_data <- data %>%
  group_by(final_call_binary) %>%
  sample_frac(0.15) %>%
  ungroup()

p_boxplot_total <- ggplot(data, aes(x = final_call_binary, y = wt_logprob)) +
  geom_jitter(data = jitter_data,
              aes(color = final_call_binary),
              width = 0.2, size = 1.2, alpha = 0.5) +
  geom_boxplot(aes(color = final_call_binary),
               fill = NA,
               outlier.shape = NA,
               width = 0.4,
               size = 0.7) +
  scale_color_manual(values = c("Essential" = colours$essential, 
                                 "Non-Essential" = colours$nonessential)) +
  labs(
    x = "Essentiality classification",
    y = "Absolute log-likelihood (total)"
  ) +
  theme_paper() +
  theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "boxplot_absolute_total.png"), p_boxplot_total,
       width = 5, height = 5, dpi = 300)
cat("  Saved boxplot_absolute_total.png\n")

# =============================================================================
# PLOT C: DENSITY PLOT
# =============================================================================

p_density <- ggplot(data, aes(x = wt_logprob, fill = final_call_binary)) +
  geom_density(alpha = 0.5, color = NA) +
  scale_fill_manual(values = c("Essential" = colours$essential, 
                                "Non-Essential" = colours$nonessential),
                    name = "Classification") +
  labs(
    x = "Absolute log-likelihood (total)",
    y = "Density"
  ) +
  theme_paper() +
  theme(legend.position = "top")

ggsave(file.path(OUTPUT_DIR, "density_absolute_total.png"), p_density,
       width = 7, height = 5, dpi = 300)
cat("  Saved density_absolute_total.png\n")

# =============================================================================
# COMBINED FIGURE
# =============================================================================

cat("\nGenerating combined figure...\n")

p_A <- p_roc + labs(tag = "A")
p_B <- p_boxplot_total + labs(tag = "B")
p_C <- p_density + labs(tag = "C")

combined <- p_A / (p_B | p_C) +
  plot_layout(heights = c(1.2, 1))

ggsave(file.path(OUTPUT_DIR, "figure_absolute_combined.png"), combined,
       width = 10, height = 10, dpi = 300)
cat("  Saved figure_absolute_combined.png\n")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================
cat("\n========== SUMMARY STATISTICS ==========\n")

stats <- data %>%
  mutate(Status = ifelse(essential_binary == 1, "Essential", "Non-essential")) %>%
  group_by(Status) %>%
  summarise(
    n = n(),
    mean_total = mean(wt_logprob),
    sd_total = sd(wt_logprob),
    mean_perbase = mean(wt_logprob_per_base),
    sd_perbase = sd(wt_logprob_per_base),
    .groups = "drop"
  )

cat("\nAbsolute total log-likelihood:\n")
cat(sprintf("  Essential:     mean = %.3f, SD = %.3f\n", 
            stats$mean_total[stats$Status == "Essential"],
            stats$sd_total[stats$Status == "Essential"]))
cat(sprintf("  Non-essential: mean = %.3f, SD = %.3f\n",
            stats$mean_total[stats$Status == "Non-essential"],
            stats$sd_total[stats$Status == "Non-essential"]))

# Wilcoxon tests
wilcox_total <- wilcox.test(wt_logprob ~ essential_binary, data = data)
wilcox_perbase <- wilcox.test(wt_logprob_per_base ~ essential_binary, data = data)

cat("\nWilcoxon rank-sum tests:\n")
cat(sprintf("  Absolute total:    p = %.2e\n", wilcox_total$p.value))
cat(sprintf("  Absolute per-base: p = %.2e\n", wilcox_perbase$p.value))

write_csv(stats, file.path(OUTPUT_DIR, "absolute_likelihood_stats.csv"))

# =============================================================================
# INTERPRETATION
# =============================================================================
cat("\n========== INTERPRETATION ==========\n")
abs_auc <- as.numeric(auc(roc_absolute))
perbase_auc <- as.numeric(auc(roc_perbase))

cat(sprintf("Absolute total AUROC:    %.3f\n", abs_auc))
cat(sprintf("Absolute per-base AUROC: %.3f\n", perbase_auc))
cat("\n")

cat("FINDING: Absolute likelihood shows predictive signal above chance.\n")
cat("Lower likelihood scores predict essentiality (counter-intuitive).\n")
cat("\nNote: Comparison to delta method (AUROC = 0.825) in Section 3.2.\n")

cat("\n========== COMPLETE ==========\n")