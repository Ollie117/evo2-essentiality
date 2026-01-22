#!/usr/bin/env Rscript
# Publication-quality plots for baseline validation
#
# Style: Pastel colours, low-saturation, warm-cool contrast

library(ggplot2)
library(readr)
library(dplyr)
library(pROC)
library(patchwork)

# ============================================================================
# CONFIGURATION
# ============================================================================

input_file <- "data/processed/merged_data.csv"
output_dir <- "results/validation"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# COLOUR PALETTE - Pastel, low-saturation
# ============================================================================

colours <- list(
  essential = "#e08060",
  nonessential = "#6a9bc3",
  primary = "#6a9bc3",
  grey = "#8c8c8c"
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
data <- read_csv(input_file, col_types = cols())
cat(sprintf("  Loaded %d genes\n", nrow(data)))
cat(sprintf("  Essential: %d\n", sum(data$is_essential == 1)))
cat(sprintf("  Non-essential: %d\n", sum(data$is_essential == 0)))

# ============================================================================
# ROC ANALYSIS
# ============================================================================

predictor <- -data$delta_ll
response <- data$is_essential

roc_obj <- roc(response, predictor, levels = c(0, 1), direction = "<")
auc_val <- auc(roc_obj)
ci_val <- ci.auc(roc_obj)

cat(sprintf("\nAUROC: %.3f (95%% CI: %.3f - %.3f)\n", auc_val, ci_val[1], ci_val[3]))

roc_data <- data.frame(
  specificity = roc_obj$specificities,
  sensitivity = roc_obj$sensitivities
)

# ============================================================================
# PLOT A: ROC CURVE
# ============================================================================

cat("\nGenerating plots...\n")

auc_label <- sprintf("AUROC = %.3f\n(95%% CI: %.3f–%.3f)", auc_val, ci_val[1], ci_val[3])

p_roc <- ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity)) +
  geom_hline(yintercept = 0.5, color = "grey80", size = 0.3) +
  geom_vline(xintercept = 0.5, color = "grey80", size = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
              color = colours$grey, size = 0.5) +
  geom_line(color = colours$primary, size = 1, alpha = 0.85) +
  annotate("text", x = 0.6, y = 0.2, label = auc_label, 
           size = 3.5, hjust = 0, fontface = "bold") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(
    x = "False positive rate (1 - Specificity)",
    y = "True positive rate (Sensitivity)"
  ) +
  theme_paper() +
  coord_equal()

ggsave(file.path(output_dir, "roc_curve_baseline.png"), p_roc,
       width = 6, height = 6, dpi = 300)
cat("  Saved roc_curve_baseline.png\n")

# ============================================================================
# PLOT B: BOXPLOT WITH JITTERED POINTS
# ============================================================================

set.seed(42)
jitter_data <- data %>%
  group_by(final_call_binary) %>%
  sample_frac(0.15) %>%
  ungroup()

p_boxplot <- ggplot(data, aes(x = final_call_binary, y = delta_ll)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = colours$grey, size = 0.5) +
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
    y = "Delta log-likelihood (Δℓ)"
  ) +
  theme_paper() +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "delta_ll_boxplot_baseline.png"), p_boxplot,
       width = 5, height = 5, dpi = 300)
cat("  Saved delta_ll_boxplot_baseline.png\n")

# ============================================================================
# PLOT C: DENSITY PLOT
# ============================================================================

p_density <- ggplot(data, aes(x = delta_ll, fill = final_call_binary)) +
  geom_density(alpha = 0.5, color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", color = colours$grey, size = 0.5) +
  scale_fill_manual(values = c("Essential" = colours$essential, 
                                "Non-Essential" = colours$nonessential),
                    name = "Classification") +
  labs(
    x = "Delta log-likelihood (Δℓ)",
    y = "Density"
  ) +
  theme_paper() +
  theme(legend.position = "top")

ggsave(file.path(output_dir, "delta_ll_density_baseline.png"), p_density,
       width = 7, height = 5, dpi = 300)
cat("  Saved delta_ll_density_baseline.png\n")

# ============================================================================
# COMBINED FIGURE
# ============================================================================

cat("\nGenerating combined figure...\n")

p_A <- p_roc + labs(tag = "A")
p_B <- p_boxplot + labs(tag = "B")
p_C <- p_density + labs(tag = "C")

combined <- p_A / (p_B | p_C) +
  plot_layout(heights = c(1.2, 1))

ggsave(file.path(output_dir, "figure_baseline_combined.png"), combined,
       width = 10, height = 10, dpi = 300)
cat("  Saved figure_baseline_combined.png\n")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

summary_stats <- data %>%
  group_by(final_call_binary) %>%
  summarise(
    n = n(),
    mean = mean(delta_ll, na.rm = TRUE),
    median = median(delta_ll, na.rm = TRUE),
    sd = sd(delta_ll, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(summary_stats, file.path(output_dir, "delta_ll_summary_baseline.csv"))

cat("\n=== Summary Statistics ===\n")
print(summary_stats)

cat("\n=== All plots saved to", output_dir, "===\n")