#!/usr/bin/env Rscript
# Perform ROC analysis for mutation sensitivity levels 1-4.
#
# Calculates AUC, sensitivity, specificity, and other performance metrics
# for each sensitivity level. Generates ROC curves and distribution plots.

library(dplyr)
library(readr)
library(pROC)

# ============================================================================
# CONFIGURATION
# ============================================================================

merged_data_file <- "C:/Users/ollfo/OneDrive/Desktop/Projects/evo2/evo2_essentiality/data/processed/merged_data_mutation_sensitivity.csv"
output_dir <- "C:/Users/ollfo/OneDrive/Desktop/Projects/evo2/evo2_essentiality/results/mutation_sensitivity"

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
# ROC ANALYSIS FOR ALL SENSITIVITY LEVELS
# ============================================================================

# Store results for all sensitivity levels
all_results <- list()

for (sens_level in 1:4) {
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat(sprintf("SENSITIVITY LEVEL %d\n", sens_level))
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  # Get the appropriate delta_ll column for this sensitivity level
  delta_ll_col <- paste0("delta_ll_sens", sens_level)
  
  cat(sprintf("\nPerforming ROC analysis for %s...\n", delta_ll_col))
  
  # Create ROC curve
  # pROC expects: response = 0/1, predictor = continuous score
  # We use negative delta_ll because negative values indicate essentiality
  roc_obj <- roc(
    response = merged_data$is_essential,
    predictor = -merged_data[[delta_ll_col]],
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
  
  # ========================================================================
  # PERFORMANCE METRICS
  # ========================================================================
  
  cat(sprintf("\nPerformance Metrics at optimal threshold:\n"))
  
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
  
  # Create results dataframe for this sensitivity level
  performance_metrics <- data.frame(
    Sensitivity_Level = sens_level,
    Metric = c("AUC", "Sensitivity", "Specificity", "PPV", "NPV", "Accuracy", "Optimal_Threshold"),
    Value = c(
      auc_value,
      perf_at_threshold$sensitivity,
      perf_at_threshold$specificity,
      perf_at_threshold$ppv,
      perf_at_threshold$npv,
      perf_at_threshold$accuracy,
      optimal_threshold
    ),
    CI_Lower = c(auc_ci[1], rep(NA, 6)),
    CI_Upper = c(auc_ci[3], rep(NA, 6))
  )
  
  # Store results
  all_results[[sens_level]] <- list(
    performance_metrics = performance_metrics,
    roc_obj = roc_obj,
    optimal_threshold = optimal_threshold,
    merged_data = merged_data,
    delta_ll_col = delta_ll_col
  )
  
}

# ========================================================================
# SAVE RESULTS FOR EACH SENSITIVITY LEVEL
# ========================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("SAVING RESULTS\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Combine all performance metrics into one dataframe
all_performance_metrics <- bind_rows(lapply(all_results, function(x) x$performance_metrics))

write_csv(all_performance_metrics, 
          file.path(output_dir, "performance_metrics_all_levels.csv"))
cat(sprintf("\n  Saved: %s/performance_metrics_all_levels.csv\n", output_dir))

# Save individual results for each sensitivity level
for (sens_level in 1:4) {
  
  results <- all_results[[sens_level]]
  roc_obj <- results$roc_obj
  optimal_threshold <- results$optimal_threshold
  merged_data <- results$merged_data
  delta_ll_col <- results$delta_ll_col
  perf_metrics <- results$performance_metrics
  
  # Save performance metrics
  perf_file <- file.path(output_dir, sprintf("performance_metrics_sens%d.csv", sens_level))
  write_csv(perf_metrics, perf_file)
  cat(sprintf("  Saved: %s\n", perf_file))
  
  # Save ROC curve data
  roc_data <- data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    Threshold = roc_obj$thresholds
  )
  roc_file <- file.path(output_dir, sprintf("roc_data_sens%d.csv", sens_level))
  write_csv(roc_data, roc_file)
  cat(sprintf("  Saved: %s\n", roc_file))
  
  # Save merged data with predictions for this sensitivity level
  merged_with_pred <- merged_data %>%
    mutate(
      evo2_prediction = if_else(-!!sym(delta_ll_col) > optimal_threshold, 
                                 "Essential", "Non-Essential"),
      prediction_correct = final_call_binary == evo2_prediction
    ) %>%
    select(gene_id, gene_name, cds_length, !!sym(delta_ll_col), 
           final_call, final_call_binary, evo2_prediction, prediction_correct)
  
  pred_file <- file.path(output_dir, sprintf("merged_predictions_sens%d.csv", sens_level))
  write_csv(merged_with_pred, pred_file)
  cat(sprintf("  Saved: %s\n", pred_file))
  
}

cat("\n=== ROC Analysis Complete ===\n")
cat(sprintf("Results saved to: %s/\n", output_dir))