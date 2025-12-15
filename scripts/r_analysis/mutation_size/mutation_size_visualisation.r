#!/usr/bin/env Rscript
# Mutation Size Analysis - Visualization
#
# Generates publication-ready plots for mutation size AUROC comparison

library(ggplot2)
library(readr)
library(dplyr)

# ============================================================================
# CONFIGURATION
# ============================================================================

input_dir <- "mutation_size"
output_dir <- "../../results/mutation_size"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading results...\n")
auc_data <- read_csv(file.path(input_dir, "auc_comparison_mutation_size.csv"), col_types = cols())

# Load ROC curve data for all sizes
roc_data_all <- data.frame()
for (i in 1:5) {
  roc_data <- read_csv(file.path(input_dir, paste0("roc_data_size", i, ".csv")), col_types = cols())
  roc_data$Size <- i
  roc_data$Pattern <- auc_data$Pattern[i]
  roc_data_all <- rbind(roc_data_all, roc_data)
}

# Create ordered factor for plotting
roc_data_all$Label <- factor(
  paste0(roc_data_all$Pattern, " (AUC=", round(auc_data$AUC[roc_data_all$Size], 3), ")"),
  levels = paste0(auc_data$Pattern, " (AUC=", round(auc_data$AUC, 3), ")")
)

# ============================================================================
# SAVE PERFORMANCE METRICS
# ============================================================================

cat("Saving performance metrics...\n")

performance_metrics <- auc_data %>%
  select(Size, Pattern, Length_bp, Stop_codons, 
         AUC, AUC_CI_Lower, AUC_CI_Upper,
         Sensitivity, Specificity, PPV, NPV, Accuracy)

write_csv(performance_metrics, file.path(output_dir, "performance_metrics_mutation_size.csv"))
cat(sprintf("  Saved: %s/performance_metrics_mutation_size.csv\n", output_dir))

# ============================================================================
# PLOT 1: AUROC BAR CHART
# ============================================================================

cat("Generating AUROC bar chart...\n")

p_auc <- ggplot(auc_data, aes(x = factor(Size), y = AUC)) +
  geom_col(aes(fill = AUC), width = 0.7) +
  geom_errorbar(aes(ymin = AUC_CI_Lower, ymax = AUC_CI_Upper), width = 0.2) +
  geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -0.5, size = 3.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  scale_fill_gradient(low = "#fee8c8", high = "#e34a33") +
  scale_x_discrete(labels = auc_data$Pattern) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(
    title = "AUROC by Mutation Size",
    subtitle = "Stop codon pattern length: 3bp (TAA) to 15bp (TAATAATAATAGTGA)",
    x = "Stop Codon Pattern",
    y = "AUROC"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "auroc_barplot_mutation_size.png"), p_auc, 
       width = 8, height = 6, dpi = 300)

# ============================================================================
# PLOT 2: COMBINED ROC CURVES
# ============================================================================

cat("Generating combined ROC curves...\n")

p_roc <- ggplot(roc_data_all, aes(x = FPR, y = TPR, color = Label)) +
  geom_line(linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(
    values = c("#fee8c8", "#fdbb84", "#fc8d59", "#e34a33", "#b30000")
  ) +
  labs(
    title = "ROC Curves by Mutation Size",
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Pattern"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    aspect.ratio = 1,
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "roc_curves_mutation_size.png"), p_roc, 
       width = 9, height = 7, dpi = 300)

# ============================================================================
# PLOT 3: AUROC TREND LINE
# ============================================================================

cat("Generating AUROC trend plot...\n")

p_trend <- ggplot(auc_data, aes(x = Length_bp, y = AUC)) +
  geom_ribbon(aes(ymin = AUC_CI_Lower, ymax = AUC_CI_Upper), alpha = 0.2, fill = "#e34a33") +
  geom_line(color = "#e34a33", linewidth = 1) +
  geom_point(size = 3, color = "#e34a33") +
  geom_text(aes(label = Pattern), vjust = -1.5, size = 3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  scale_x_continuous(breaks = auc_data$Length_bp) +
  scale_y_continuous(limits = c(0.5, 0.9), breaks = seq(0.5, 0.9, 0.05)) +
  labs(
    title = "AUROC vs Mutation Size",
    subtitle = "Performance increases with pattern length, plateauing at ~12bp",
    x = "Stop Pattern Length (bp)",
    y = "AUROC"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "auroc_trend_mutation_size.png"), p_trend, 
       width = 8, height = 6, dpi = 300)

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n=== Mutation Size Results Summary ===\n\n")
cat("Performance Metrics:\n")
print(performance_metrics)

cat("\nPlots saved to:", output_dir, "\n")
cat("  - performance_metrics_mutation_size.csv\n")
cat("  - auroc_barplot_mutation_size.png\n")
cat("  - roc_curves_mutation_size.png\n")
cat("  - auroc_trend_mutation_size.png\n")