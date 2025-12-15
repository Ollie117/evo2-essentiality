#!/usr/bin/env Rscript
# Mutation Size Analysis - ROC analysis across 5 mutation size levels
#
# Computes ROC/AUC for each mutation size (TAA through TAATAATAATAGTGA)
# to assess whether Evo2 can detect minimal perturbations.

library(dplyr)
library(readr)
library(pROC)
library(ggplot2)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input files
scores_file <- "../../../data/processed/essentiality_mutation_sensitivity3.csv"
tnseq_file <- "../../../data/processed/tnseq_clean.csv"

# Output directory
output_dir <- "../../../data/processed/mutation_size"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading mutation size scores...\n")
scores <- read_csv(scores_file, col_types = cols())
cat(sprintf("  Loaded %d genes\n", nrow(scores)))

cat("Loading TnSeq classifications...\n")
tnseq <- read_csv(tnseq_file, col_types = cols())
cat(sprintf("  Loaded %d genes\n", nrow(tnseq)))

# ============================================================================
# MERGE DATA
# ============================================================================

cat("\nMerging datasets...\n")
merged <- scores %>%
  inner_join(tnseq, by = "gene_id")
cat(sprintf("  Merged dataset: %d genes\n", nrow(merged)))
cat(sprintf("  Essential: %d\n", sum(merged$is_essential == 1)))
cat(sprintf("  Non-essential: %d\n", sum(merged$is_essential == 0)))

# ============================================================================
# ROC ANALYSIS FOR EACH SIZE LEVEL
# ============================================================================

cat("\nPerforming ROC analysis for each mutation size...\n")

# Define size levels
size_levels <- data.frame(
  size = 1:5,
  pattern = c("TAA", "TAATAA", "TAATAATAA", "TAATAATAATAG", "TAATAATAATAGTGA"),
  length_bp = c(3, 6, 9, 12, 15),
  stop_codons = c(1, 2, 3, 4, 5)
)

# Storage for results
roc_results <- list()
auc_comparison <- data.frame()

for (i in 1:5) {
  delta_col <- paste0("delta_ll_size", i)
  
  # Create ROC object
  roc_obj <- roc(
    response = merged$is_essential,
    predictor = -merged[[delta_col]],
    levels = c(0, 1),
    direction = "<"
  )
  
  # Store ROC object
  roc_results[[i]] <- roc_obj
  
  # Get AUC and CI
  auc_val <- as.numeric(roc_obj$auc)
  auc_ci <- ci.auc(roc_obj)
  
  # Get performance at optimal threshold
  perf <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity", "ppv", "npv", "accuracy"))
  
  # Add to comparison table
  auc_comparison <- rbind(auc_comparison, data.frame(
    Size = i,
    Pattern = size_levels$pattern[i],
    Length_bp = size_levels$length_bp[i],
    Stop_codons = size_levels$stop_codons[i],
    AUC = auc_val,
    AUC_CI_Lower = auc_ci[1],
    AUC_CI_Upper = auc_ci[3],
    Sensitivity = perf$sensitivity,
    Specificity = perf$specificity,
    PPV = perf$ppv,
    NPV = perf$npv,
    Accuracy = perf$accuracy
  ))
  
  cat(sprintf("  Size %d (%s): AUC = %.4f [%.4f - %.4f]\n", 
              i, size_levels$pattern[i], auc_val, auc_ci[1], auc_ci[3]))
}

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\nSaving results...\n")

# Save AUC comparison table
write_csv(auc_comparison, file.path(output_dir, "auc_comparison_mutation_size.csv"))
cat(sprintf("  Saved: %s/auc_comparison_mutation_size.csv\n", output_dir))

# Save ROC curve data for each size
for (i in 1:5) {
  roc_data <- data.frame(
    FPR = 1 - roc_results[[i]]$specificities,
    TPR = roc_results[[i]]$sensitivities,
    Threshold = roc_results[[i]]$thresholds
  )
  write_csv(roc_data, file.path(output_dir, paste0("roc_data_size", i, ".csv")))
}
cat(sprintf("  Saved: ROC data for all 5 sizes\n"))

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n=== Mutation Size Analysis Complete ===\n\n")
cat("AUC Summary:\n")
print(auc_comparison %>% select(Size, Pattern, Length_bp, AUC, AUC_CI_Lower, AUC_CI_Upper))
cat(sprintf("\nResults saved to: %s/\n", output_dir))