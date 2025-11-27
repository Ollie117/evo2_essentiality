#!/usr/bin/env Rscript
# Perform ROC analysis and generate validation metrics.
#
# Calculates AUC, sensitivity, specificity, and other performance metrics
# for Evo2 predictions. Generates ROC curves and distribution plots.

library(dplyr)
library(readr)
library(pROC)

# ============================================================================
# CONFIGURATION
# ============================================================================

merged_data_file <- "data/processed/merged_data.csv"
output_dir <- "results/validation"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading merged dataset...\n")
merged_data <- read_csv(merged_data_file, col_types = cols())
cat(sprintf("  Loaded %d genes\n", nrow(merged_data)))

# ============================================================================
# ROC ANALYSIS
# ============================================================================

cat("\nPerforming ROC analysis...\n")

# Create ROC curve
# pROC expects: response = 0/1, predictor = continuous score
# We use negative delta_ll because negative values indicate essentiality
roc_obj <- roc(
  response = merged_data$is_essential,
  predictor = -merged_data$delta_ll,
  levels = c(0, 1),
  direction = "<"  # Higher negative delta_ll = essential
)

# Calculate AUC and confidence interval
auc_value <- as.numeric(roc_obj$auc)
auc_ci <- ci.auc(roc_obj)

cat(sprintf("\nROC Analysis Results:\n"))
cat(sprintf("  AUC: %.4f\n", auc_value))
cat(sprintf("  95%% CI: [%.4f - %.4f]\n", auc_ci[1], auc_ci[3]))

# Get optimal threshold (Youden index)
optimal_coords <- coords(roc_obj, "best", ret = "threshold")
optimal_threshold <- optimal_coords$threshold

cat(sprintf("  Optimal threshold (Youden): %.4f\n", optimal_threshold))

# ============================================================================
# PERFORMANCE METRICS
# ============================================================================

cat("\nPerformance Metrics at optimal threshold:\n")

# Get all metrics at optimal threshold
perf_at_threshold <- coords(
  roc_obj,
  optimal_threshold,
  ret = c("threshold", "sensitivity", "specificity", "ppv", "npv", "accuracy")
)

cat(sprintf("  Sensitivity: %.4f\n", perf_at_threshold$sensitivity))
cat(sprintf("  Specificity: %.4f\n", perf_at_threshold$specificity))
cat(sprintf("  PPV (Precision): %.4f\n", perf_at_threshold$ppv))
cat(sprintf("  NPV: %.4f\n", perf_at_threshold$npv))
cat(sprintf("  Accuracy: %.4f\n", perf_at_threshold$accuracy))

# Create results dataframe
performance_metrics <- data.frame(
  Metric = c("AUC", "Sensitivity", "Specificity", "PPV", "NPV", "Accuracy"),
  Value = c(
    auc_value,
    perf_at_threshold$sensitivity,
    perf_at_threshold$specificity,
    perf_at_threshold$ppv,
    perf_at_threshold$npv,
    perf_at_threshold$accuracy
  ),
  CI_Lower = c(auc_ci[1], rep(NA, 5)),
  CI_Upper = c(auc_ci[3], rep(NA, 5))
)

# Note: Plotting is handled in separate script (05_plot_roc_curve.R)

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\nSaving results...\n")

# Save performance metrics
write_csv(performance_metrics, file.path(output_dir, "performance_metrics.csv"))
cat(sprintf("  Saved: %s/performance_metrics.csv\n", output_dir))

# Save ROC curve data
roc_data <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities,
  Threshold = roc_obj$thresholds
)
write_csv(roc_data, file.path(output_dir, "roc_data.csv"))
cat(sprintf("  Saved: %s/roc_data.csv\n", output_dir))

# Save merged data with predictions
merged_with_pred <- merged_data %>%
  mutate(
    evo2_prediction = if_else(-delta_ll > optimal_threshold, 
                               "Essential", "Non-Essential"),
    prediction_correct = final_call_binary == evo2_prediction
  ) %>%
  select(gene_id, gene_name, cds_length, delta_ll, final_call, 
         final_call_binary, evo2_prediction, prediction_correct)

write_csv(merged_with_pred, file.path(output_dir, "merged_predictions.csv"))
cat(sprintf("  Saved: %s/merged_predictions.csv\n", output_dir))

cat("\n=== ROC Analysis Complete ===\n")
cat(sprintf("Results saved to: %s/\n", output_dir))