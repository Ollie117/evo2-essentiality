#!/usr/bin/env Rscript
# ROC analysis for genomic context experiment
#
# Evaluates discriminative performance across 5 context window sizes.

library(pROC)
library(readr)
library(dplyr)
library(ggplot2)

# ============================================================================
# CONFIGURATION
# ============================================================================

scores_file <- "data/processed/essentiality_genomic_context.csv"
tnseq_file <- "data/processed/tnseq_clean.csv"
output_dir <- "results/genomic_context"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Context levels (bp per side)
context_levels <- c(512, 1024, 2048, 4096, 8192)

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
# ROC ANALYSIS FOR EACH CONTEXT LEVEL
# ============================================================================

cat("\nPerforming ROC analysis for each context level...\n")

auc_results <- data.frame(
  context_level = integer(),
  context_bp = integer(),
  auc = numeric(),
  auc_lower = numeric(),
  auc_upper = numeric(),
  stringsAsFactors = FALSE
)

roc_objects <- list()

for (i in 1:5) {
  context_bp <- context_levels[i]
  delta_col <- paste0("delta_ll_context", i)
  
  # Invert scores (more negative = more essential = higher rank)
  predictor <- -merged[[delta_col]]
  response <- merged$is_essential
  
  # Remove NAs
  valid_idx <- !is.na(predictor)
  predictor <- predictor[valid_idx]
  response <- response[valid_idx]
  
  # Compute ROC
  roc_obj <- roc(response, predictor, levels = c(0, 1), direction = "<")
  roc_objects[[i]] <- roc_obj
  
  # Get AUC with CI
  auc_val <- auc(roc_obj)
  ci_val <- ci.auc(roc_obj)
  
  auc_results <- rbind(auc_results, data.frame(
    context_level = i,
    context_bp = context_bp,
    auc = as.numeric(auc_val),
    auc_lower = ci_val[1],
    auc_upper = ci_val[3]
  ))
  
  cat(sprintf("  Context %d (%d bp): AUC = %.3f (95%% CI: %.3f - %.3f)\n",
              i, context_bp, auc_val, ci_val[1], ci_val[3]))
  
  # Save ROC curve data
  roc_data <- data.frame(
    specificity = roc_obj$specificities,
    sensitivity = roc_obj$sensitivities
  )
  write_csv(roc_data, file.path(output_dir, sprintf("roc_data_context%d.csv", i)))
}

# Save AUC comparison
write_csv(auc_results, file.path(output_dir, "auc_comparison_genomic_context.csv"))
cat(sprintf("\nSaved AUC results to %s/auc_comparison_genomic_context.csv\n", output_dir))

# ============================================================================
# PERFORMANCE METRICS AT OPTIMAL THRESHOLD
# ============================================================================

cat("\nComputing performance metrics at optimal threshold...\n")

metrics_results <- data.frame(
  context_level = integer(),
  context_bp = integer(),
  threshold = numeric(),
  sensitivity = numeric(),
  specificity = numeric(),
  ppv = numeric(),
  npv = numeric(),
  accuracy = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:5) {
  context_bp <- context_levels[i]
  roc_obj <- roc_objects[[i]]
  
  # Optimal threshold via Youden index
  coords_best <- coords(roc_obj, "best", best.method = "youden")
  
  # Calculate metrics
  threshold <- coords_best$threshold
  sens <- coords_best$sensitivity
  spec <- coords_best$specificity
  
  # Get confusion matrix values for PPV/NPV
  delta_col <- paste0("delta_ll_context", i)
  predictor <- -merged[[delta_col]]
  response <- merged$is_essential
  valid_idx <- !is.na(predictor)
  predictor <- predictor[valid_idx]
  response <- response[valid_idx]
  
  predicted <- ifelse(predictor >= threshold, 1, 0)
  
  tp <- sum(predicted == 1 & response == 1)
  fp <- sum(predicted == 1 & response == 0)
  tn <- sum(predicted == 0 & response == 0)
  fn <- sum(predicted == 0 & response == 1)
  
  ppv <- tp / (tp + fp)
  npv <- tn / (tn + fn)
  accuracy <- (tp + tn) / (tp + fp + tn + fn)
  
  metrics_results <- rbind(metrics_results, data.frame(
    context_level = i,
    context_bp = context_bp,
    threshold = threshold,
    sensitivity = sens,
    specificity = spec,
    ppv = ppv,
    npv = npv,
    accuracy = accuracy
  ))
}

write_csv(metrics_results, file.path(output_dir, "performance_metrics_genomic_context.csv"))
cat(sprintf("Saved performance metrics to %s/performance_metrics_genomic_context.csv\n", output_dir))

print(metrics_results)

cat("\n=== ROC Analysis Complete ===\n")