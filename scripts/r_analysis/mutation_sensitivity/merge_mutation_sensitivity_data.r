#!/usr/bin/env Rscript
# Merge Evo2 mutation sensitivity predictions with TnSeq experimental data.
#
# Combines cleaned Evo2 mutation sensitivity scores and TnSeq essentiality 
# classifications on gene ID to create unified datasets for each sensitivity level.

library(dplyr)
library(readr)

# ============================================================================
# CONFIGURATION
# ============================================================================

evo2_sens_file <- "C:/Users/ollfo/OneDrive/Desktop/Projects/evo2/evo2_essentiality/data/processed/evo2_mutation_sensitivity_clean.csv"
tnseq_file <- "C:/Users/ollfo/OneDrive/Desktop/Projects/evo2/evo2_essentiality/data/processed/tnseq_clean.csv"
output_dir <- "C:/Users/ollfo/OneDrive/Desktop/Projects/evo2/evo2_essentiality/data/processed"

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading cleaned data files...\n")
evo2_sens_scores <- read_csv(evo2_sens_file, col_types = cols())
tnseq_data <- read_csv(tnseq_file, col_types = cols())

cat(sprintf("  Evo2 mutation sensitivity scores: %d genes\n", nrow(evo2_sens_scores)))
cat(sprintf("  TnSeq experimental data: %d genes\n", nrow(tnseq_data)))

# ============================================================================
# MERGE DATASETS
# ============================================================================

cat("\nMerging datasets on gene_id...\n")

merged_data <- evo2_sens_scores %>%
  inner_join(tnseq_data, by = "gene_id")

cat(sprintf("  Merged dataset: %d genes\n", nrow(merged_data)))

# Check for genes that didn't match
evo2_unmatched <- setdiff(evo2_sens_scores$gene_id, merged_data$gene_id)
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

# ============================================================================
# DELTA_LL STATISTICS BY ESSENTIALITY AND SENSITIVITY LEVEL
# ============================================================================

cat("\n=== Delta_ll Statistics by Essentiality and Sensitivity Level ===\n\n")

for (sens_level in 1:4) {
  col_name <- paste0("delta_ll_sens", sens_level)
  
  cat(sprintf("Sensitivity Level %d:\n", sens_level))
  
  summary_stats <- merged_data %>%
    group_by(final_call_binary) %>%
    summarise(
      n = n(),
      mean_delta_ll = mean(.data[[col_name]], na.rm = TRUE),
      median_delta_ll = median(.data[[col_name]], na.rm = TRUE),
      sd_delta_ll = sd(.data[[col_name]], na.rm = TRUE),
      min_delta_ll = min(.data[[col_name]]),
      max_delta_ll = max(.data[[col_name]]),
      .groups = "drop"
    )
  
  print(summary_stats)
  cat("\n")
}

# ============================================================================
# SAVE MERGED DATA
# ============================================================================

output_file <- file.path(output_dir, "merged_data_mutation_sensitivity.csv")
cat(sprintf("Saving merged dataset to: %s\n", output_file))
write_csv(merged_data, output_file)

cat("=== Merge Complete ===\n")