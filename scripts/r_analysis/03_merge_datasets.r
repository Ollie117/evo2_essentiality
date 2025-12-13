#!/usr/bin/env Rscript
# Merge Evo2 predictions with TnSeq experimental data.
#
# Combines cleaned Evo2 scores and TnSeq essentiality classifications
# on gene ID to create a unified dataset for validation.

library(dplyr)
library(readr)

# ============================================================================
# CONFIGURATION
# ============================================================================

evo2_file <- "data/processed/evo2_scores_clean.csv"
tnseq_file <- "data/processed/tnseq_clean.csv"
output_file <- "data/processed/merged_data.csv"

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading cleaned data files...\n")
evo2_scores <- read_csv(evo2_file, col_types = cols())
tnseq_data <- read_csv(tnseq_file, col_types = cols())

cat(sprintf("  Evo2 scores: %d genes\n", nrow(evo2_scores)))
cat(sprintf("  TnSeq data: %d genes\n", nrow(tnseq_data)))

# ============================================================================
# MERGE DATASETS
# ============================================================================

cat("\nMerging datasets on gene_id...\n")

merged_data <- evo2_scores %>%
  inner_join(tnseq_data, by = "gene_id")

cat(sprintf("  Merged dataset: %d genes\n", nrow(merged_data)))

# Check for genes that didn't match
evo2_unmatched <- setdiff(evo2_scores$gene_id, merged_data$gene_id)
tnseq_unmatched <- setdiff(tnseq_data$gene_id, merged_data$gene_id)

if (length(evo2_unmatched) > 0) {
  cat(sprintf("\n  Warning: %d Evo2 genes not found in TnSeq data\n", length(evo2_unmatched)))
}

if (length(tnseq_unmatched) > 0) {
  cat(sprintf("  Warning: %d TnSeq genes not found in Evo2 data\n", length(tnseq_unmatched)))
}

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\nMerged dataset summary:\n")
cat(sprintf("  Essential genes: %d\n", sum(merged_data$is_essential == 1)))
cat(sprintf("  Non-Essential genes: %d\n", sum(merged_data$is_essential == 0)))

cat("\nDelta_ll statistics by essentiality:\n")
summary_stats <- merged_data %>%
  group_by(final_call_binary) %>%
  summarise(
    n = n(),
    mean_delta_ll = mean(delta_ll, na.rm = TRUE),
    median_delta_ll = median(delta_ll, na.rm = TRUE),
    sd_delta_ll = sd(delta_ll, na.rm = TRUE),
    min_delta_ll = min(delta_ll),
    max_delta_ll = max(delta_ll),
    .groups = "drop"
  )

print(summary_stats)

# ============================================================================
# SAVE MERGED DATA
# ============================================================================

cat(sprintf("\nSaving merged dataset to: %s\n", output_file))
write_csv(merged_data, output_file)

cat("=== Merge Complete ===\n")