#!/usr/bin/env Rscript
# Generate ROC plots for mutation sensitivity analysis
#
# Creates individual ROC curves for each sensitivity level and a 
# combined comparison plot showing all four levels

library(ggplot2)
library(readr)
library(dplyr)

# ============================================================================
# CONFIGURATION
# ============================================================================

data_dir <- "C:/Users/ollfo/OneDrive/Desktop/Projects/evo2/evo2_essentiality/results/mutation_sensitivity"
output_dir <- "C:/Users/ollfo/OneDrive/Desktop/Projects/evo2/evo2_essentiality/results/mutation_sensitivity"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# LOAD PERFORMANCE METRICS
# ============================================================================

cat("Loading performance metrics...\n")
perf_all <- read_csv(file.path(data_dir, "performance_metrics_all_levels.csv"), 
                     col_types = cols())

# Extract AUC values for each sensitivity level
auc_df <- perf_all %>%
  filter(Metric == "AUC") %>%
  select(Sensitivity_Level, Value)

auc_values <- setNames(auc_df$Value, auc_df$Sensitivity_Level)

cat(sprintf("  Sensitivity Level 1 AUC: %.4f\n", auc_values[1]))
cat(sprintf("  Sensitivity Level 2 AUC: %.4f\n", auc_values[2]))
cat(sprintf("  Sensitivity Level 3 AUC: %.4f\n", auc_values[3]))
cat(sprintf("  Sensitivity Level 4 AUC: %.4f\n", auc_values[4]))

# ============================================================================
# CREATE INDIVIDUAL ROC PLOTS
# ============================================================================

cat("\nCreating individual ROC plots...\n")

# Define colors for each sensitivity level
colors <- c("steelblue", "orange", "red", "darkred")
sens_labels <- c("Sensitivity 1 (Small)", "Sensitivity 2 (Medium)", 
                 "Sensitivity 3 (Large)", "Sensitivity 4 (Very Large)")

for (sens_level in 1:4) {
  
  cat(sprintf("  Plotting Sensitivity Level %d...\n", sens_level))
  
  # Load ROC data for this sensitivity level
  roc_file <- file.path(data_dir, sprintf("roc_data_sens%d.csv", sens_level))
  roc_data <- read_csv(roc_file, col_types = cols())
  
  # Get AUC for this level
  auc_val <- auc_values[sens_level]
  
  # Create ROC plot
  p_roc <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
    geom_line(size = 1.2, color = colors[sens_level]) +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color = "grey50", linetype = "dashed", size = 0.8) +
    labs(
      title = sprintf("ROC Curve - Evo2 Gene Essentiality Prediction\n%s\nAUC = %.4f", 
                      sens_labels[sens_level], auc_val),
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    xlim(0, 1) + ylim(0, 1) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      aspect.ratio = 1
    )
  
  # Save plot
  plot_file <- file.path(output_dir, sprintf("roc_curve_sens%d.png", sens_level))
  png(plot_file, width = 800, height = 800, res = 120)
  print(p_roc)
  dev.off()
  
  cat(sprintf("    Saved: %s\n", plot_file))
  
}

# ============================================================================
# CREATE COMPARISON PLOT
# ============================================================================

cat("\nCreating comparison plot...\n")

# Load all ROC data and add sensitivity level indicator
roc_data_all <- tibble()

for (sens_level in 1:4) {
  roc_file <- file.path(data_dir, sprintf("roc_data_sens%d.csv", sens_level))
  roc_data <- read_csv(roc_file, col_types = cols())
  roc_data$Sensitivity_Level <- sens_level
  roc_data$Label <- sens_labels[sens_level]
  roc_data_all <- bind_rows(roc_data_all, roc_data)
}

# Create comparison plot
p_comparison <- ggplot(roc_data_all, aes(x = FPR, y = TPR, color = Label, linetype = Label)) +
  geom_line(size = 1.2) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1, color = NA, linetype = NA), 
               color = "grey50", linetype = "dashed", size = 0.8) +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) +
  labs(
    title = "ROC Curves Comparison - Evo2 Performance Across Mutation Sensitivity Levels",
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Mutation Size",
    linetype = "Mutation Size"
  ) +
  xlim(0, 1) + ylim(0, 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "lower right",
    legend.background = element_rect(fill = "white", color = "grey50"),
    aspect.ratio = 1
  )

# Save comparison plot
comparison_file <- file.path(output_dir, "roc_curve_comparison.png")
png(comparison_file, width = 1000, height = 800, res = 120)
print(p_comparison)
dev.off()

cat(sprintf("  Saved: %s\n", comparison_file))

# ============================================================================
# CREATE AUC COMPARISON BAR PLOT
# ============================================================================

cat("\nCreating AUC comparison bar plot...\n")

auc_df <- data.frame(
  Sensitivity_Level = factor(1:4, labels = sens_labels),
  AUC = c(auc_values[1], auc_values[2], auc_values[3], auc_values[4])
)

p_auc <- ggplot(auc_df, aes(x = Sensitivity_Level, y = AUC, fill = Sensitivity_Level)) +
  geom_bar(stat = "identity", color = "black", size = 0.8) +
  scale_fill_manual(values = colors) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50", size = 0.8) +
  geom_text(aes(label = sprintf("%.4f", AUC)), vjust = -0.5, size = 4, fontface = "bold") +
  ylim(0, 1) +
  labs(
    title = "AUC Comparison Across Mutation Sensitivity Levels",
    x = "Mutation Sensitivity Level",
    y = "AUC"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Save AUC bar plot
auc_file <- file.path(output_dir, "auc_comparison.png")
png(auc_file, width = 900, height = 700, res = 120)
print(p_auc)
dev.off()

cat(sprintf("  Saved: %s\n", auc_file))

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n=== Plotting Complete ===\n")
cat(sprintf("\nGenerated plots:\n"))
cat(sprintf("  - roc_curve_sens1.png (AUC = %.4f)\n", auc_values[1]))
cat(sprintf("  - roc_curve_sens2.png (AUC = %.4f)\n", auc_values[2]))
cat(sprintf("  - roc_curve_sens3.png (AUC = %.4f)\n", auc_values[3]))
cat(sprintf("  - roc_curve_sens4.png (AUC = %.4f)\n", auc_values[4]))
cat(sprintf("  - roc_curve_comparison.png\n"))
cat(sprintf("  - auc_comparison.png\n"))
cat(sprintf("\nAll plots saved to: %s\n", output_dir))