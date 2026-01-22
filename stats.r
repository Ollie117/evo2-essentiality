#!/usr/bin/env Rscript
# Extract all numbers for thesis
#
# Outputs all [X] placeholder values needed for Methods, Results, Discussion

library(readr)
library(dplyr)
library(tidyr)
library(pROC)

# ============================================================================
# OUTPUT FILE
# ============================================================================

output_file <- "results/thesis_numbers.txt"
sink(output_file)

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("THESIS NUMBERS - Generated", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ============================================================================
# 1. DATASET COUNTS (Methods 2.1)
# ============================================================================

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("1. DATASET COUNTS (Methods 2.1)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

tnseq <- read_csv("data/processed/tnseq_clean.csv", col_types = cols(), show_col_types = FALSE)

n_essential <- sum(tnseq$is_essential == 1)
n_nonessential <- sum(tnseq$is_essential == 0)
n_total <- nrow(tnseq)

cat(sprintf("Total genes: %d\n", n_total))
cat(sprintf("Essential genes: %d\n", n_essential))
cat(sprintf("Non-essential genes: %d\n", n_nonessential))
cat(sprintf("\nFor thesis: 'The final dataset comprised %d essential and %d non-essential genes.'\n", 
            n_essential, n_nonessential))

# ============================================================================
# 2. BASELINE VALIDATION (Results 3.1)
# ============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("2. BASELINE VALIDATION (Results 3.1)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

baseline <- read_csv("data/processed/merged_data.csv", col_types = cols(), show_col_types = FALSE)

# ROC analysis
predictor <- -baseline$delta_ll
response <- baseline$is_essential
roc_obj <- roc(response, predictor, levels = c(0, 1), direction = "<", quiet = TRUE)
auc_val <- auc(roc_obj)
ci_val <- ci.auc(roc_obj)

cat("ROC Performance:\n")
cat(sprintf("  AUROC: %.3f\n", auc_val))
cat(sprintf("  95%% CI: %.3f - %.3f\n", ci_val[1], ci_val[3]))

# Optimal threshold metrics
coords_best <- coords(roc_obj, "best", best.method = "youden")
threshold <- coords_best$threshold
sens <- coords_best$sensitivity
spec <- coords_best$specificity

# Confusion matrix
predicted <- ifelse(predictor >= threshold, 1, 0)
tp <- sum(predicted == 1 & response == 1)
fp <- sum(predicted == 1 & response == 0)
tn <- sum(predicted == 0 & response == 0)
fn <- sum(predicted == 0 & response == 1)

ppv <- tp / (tp + fp)
npv <- tn / (tn + fn)
accuracy <- (tp + tn) / (tp + fp + tn + fn)

cat("\nPerformance at optimal threshold:\n")
cat(sprintf("  Threshold: %.4f\n", threshold))
cat(sprintf("  Sensitivity: %.3f (%.1f%%)\n", sens, sens * 100))
cat(sprintf("  Specificity: %.3f (%.1f%%)\n", spec, spec * 100))
cat(sprintf("  PPV: %.3f (%.1f%%)\n", ppv, ppv * 100))
cat(sprintf("  NPV: %.3f (%.1f%%)\n", npv, npv * 100))
cat(sprintf("  Accuracy: %.3f (%.1f%%)\n", accuracy, accuracy * 100))

# Distribution statistics
essential_stats <- baseline %>% filter(is_essential == 1) %>% pull(delta_ll)
nonessential_stats <- baseline %>% filter(is_essential == 0) %>% pull(delta_ll)

cat("\nDelta log-likelihood distributions:\n")
cat(sprintf("  Essential: mean = %.4f, SD = %.4f, median = %.4f\n", 
            mean(essential_stats), sd(essential_stats), median(essential_stats)))
cat(sprintf("  Non-essential: mean = %.4f, SD = %.4f, median = %.4f\n", 
            mean(nonessential_stats), sd(nonessential_stats), median(nonessential_stats)))
cat(sprintf("  Separation (NE - E): %.4f\n", mean(nonessential_stats) - mean(essential_stats)))

# ============================================================================
# 3. MUTATION SIZE ANALYSIS (Results 3.2)
# ============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("3. MUTATION SIZE ANALYSIS (Results 3.2)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Load AUC comparison
if (file.exists("data/processed/mutation_size/auc_comparison_mutation_size.csv")) {
  auc_mutation <- read_csv("data/processed/mutation_size/auc_comparison_mutation_size.csv", 
                           col_types = cols(), show_col_types = FALSE)
  
  cat("AUROC by pattern length:\n")
  cat(sprintf("%-12s  %s  %s\n", "Pattern", "AUROC", "95% CI"))
  cat("-" |> rep(45) |> paste(collapse = ""), "\n")
  for (i in 1:nrow(auc_mutation)) {
    cat(sprintf("%-12s  %.3f  (%.3f - %.3f)\n", 
                auc_mutation$pattern_length[i],
                auc_mutation$auc[i],
                auc_mutation$auc_lower[i],
                auc_mutation$auc_upper[i]))
  }
} else {
  cat("AUC comparison file not found. Run mutation_size_roc_analysis.R first.\n")
}

# Load scores for separation analysis
if (file.exists("data/processed/essentiality_scores_mutation_size.csv")) {
  scores_mut <- read_csv("data/processed/essentiality_scores_mutation_size.csv", 
                         col_types = cols(), show_col_types = FALSE)
  scores_mut <- scores_mut %>% inner_join(tnseq, by = "gene_id")
  
  cat("\nMean separation by pattern length:\n")
  pattern_lengths <- c("3 bp", "6 bp", "9 bp", "12 bp", "15 bp")
  for (i in 1:5) {
    col_name <- paste0("delta_ll_size", i)
    essential_mean <- mean(scores_mut[[col_name]][scores_mut$is_essential == 1], na.rm = TRUE)
    nonessential_mean <- mean(scores_mut[[col_name]][scores_mut$is_essential == 0], na.rm = TRUE)
    separation <- nonessential_mean - essential_mean
    cat(sprintf("  %s: %.4f (NE mean: %.4f, E mean: %.4f)\n", 
                pattern_lengths[i], separation, nonessential_mean, essential_mean))
  }
}

# ============================================================================
# 4. GENOMIC CONTEXT ANALYSIS (Results 3.3)
# ============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("4. GENOMIC CONTEXT ANALYSIS (Results 3.3)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Load AUC comparison
if (file.exists("data/processed/genomic_context/auc_comparison_genomic_context.csv")) {
  auc_context <- read_csv("data/processed/genomic_context/auc_comparison_genomic_context.csv", 
                          col_types = cols(), show_col_types = FALSE)
  
  cat("AUROC by context length:\n")
  cat(sprintf("%-15s  %s  %s\n", "Context", "AUROC", "95% CI"))
  cat("-" |> rep(45) |> paste(collapse = ""), "\n")
  for (i in 1:nrow(auc_context)) {
    cat(sprintf("%-15s  %.3f  (%.3f - %.3f)\n", 
                paste0(auc_context$context_bp[i], " bp"),
                auc_context$auc[i],
                auc_context$auc_lower[i],
                auc_context$auc_upper[i]))
  }
} else {
  cat("AUC comparison file not found. Run genomic_context_roc_analysis.R first.\n")
}

# Load scores for separation analysis
# Check which filename exists
context_file <- if (file.exists("data/processed/essentiality_scores_genomic_context.csv")) {
  "data/processed/essentiality_scores_genomic_context.csv"
} else if (file.exists("data/processed/essentiality_genomic_context.csv")) {
  "data/processed/essentiality_genomic_context.csv"
} else {
  NULL
}

if (!is.null(context_file)) {
  scores_ctx <- read_csv(context_file, col_types = cols(), show_col_types = FALSE)
  scores_ctx <- scores_ctx %>% inner_join(tnseq, by = "gene_id")
  
  cat("\nMean separation by context length:\n")
  context_lengths <- c("512 bp", "1024 bp", "2048 bp", "4096 bp", "8192 bp")
  for (i in 1:5) {
    col_name <- paste0("delta_ll_context", i)
    if (col_name %in% names(scores_ctx)) {
      essential_mean <- mean(scores_ctx[[col_name]][scores_ctx$is_essential == 1], na.rm = TRUE)
      nonessential_mean <- mean(scores_ctx[[col_name]][scores_ctx$is_essential == 0], na.rm = TRUE)
      separation <- nonessential_mean - essential_mean
      cat(sprintf("  %s: %.4f (NE mean: %.4f, E mean: %.4f)\n", 
                  context_lengths[i], separation, nonessential_mean, essential_mean))
    }
  }
}

# ============================================================================
# 5. SUMMARY TABLE FOR COPY-PASTE
# ============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("5. COPY-PASTE READY TEXT\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("METHODS 2.1:\n")
cat(sprintf("'The final dataset comprised %d essential and %d non-essential genes.'\n\n", 
            n_essential, n_nonessential))

cat("RESULTS 3.1 (Baseline):\n")
cat(sprintf("'AUROC of %.3f (95%% CI: %.3f–%.3f)'\n", auc_val, ci_val[1], ci_val[3]))
cat(sprintf("'sensitivity of %.1f%%, specificity of %.1f%%, and overall accuracy of %.1f%%'\n", 
            sens * 100, spec * 100, accuracy * 100))
cat(sprintf("'Essential genes: mean = %.3f, SD = %.3f'\n", mean(essential_stats), sd(essential_stats)))
cat(sprintf("'Non-essential genes: mean = %.3f, SD = %.3f'\n", mean(nonessential_stats), sd(nonessential_stats)))

sink()

cat("\n\nResults saved to:", output_file, "\n")
cat("Open this file for all thesis numbers.\n")