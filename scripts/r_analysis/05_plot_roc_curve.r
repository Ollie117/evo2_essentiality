#!/usr/bin/env Rscript
# Regenerate ROC plot from saved ROC data

library(ggplot2)
library(readr)

# Load ROC data
roc_data <- read_csv("results/validation/roc_data.csv")

# Load performance metrics
perf <- read_csv("results/validation/performance_metrics.csv")
auc_value <- perf$Value[1]

# Create ROC plot
p_roc <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
  geom_line(size = 1.2, color = "steelblue") +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color = "grey50", linetype = "dashed", size = 0.8) +
  labs(
    title = sprintf("ROC Curve - Evo2 Gene Essentiality Prediction\nAUC = %.4f", auc_value),
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
png("results/validation/roc_curve.png", width = 800, height = 800, res = 120)
print(p_roc)
dev.off()

cat("ROC plot saved to: results/validation/roc_curve.png\n")