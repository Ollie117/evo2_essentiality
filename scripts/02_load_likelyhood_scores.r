#!/usr/bin/env Rscript
# Load Evo2 gene essentiality prediction scores.
#
# Reads the CSV file containing delta_ll scores for all genes.

library(dplyr)
library(readr)

# ============================================================================
# CONFIGURATION
# ============================================================================

evo2_scores_file <- "data/processed/essentiality_scores_test.csv"
output_file <- "data/processed/evo2_scores_clean.csv"

# Create output directory if it doesn't exist
if (!dir.exists("data/processed")) {
  dir.create("data/processed", recursive = TRUE)
}

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading Evo2 scores...\n")
evo2_scores <- read_csv(evo2_scores_file, col_types = cols())
cat(sprintf("  Loaded %d genes\n", nrow(evo2_scores)))

# ============================================================================
# DATA CLEANING
# ============================================================================

cat("\nCleaning Evo2 data...\n")

# Clean gene IDs (remove whitespace)
evo2_clean <- evo2_scores %>%
  mutate(gene_id = trimws(gene_id))

# Check for missing values in key columns
missing_summary <- evo2_clean %>%
  summarise(
    missing_gene_id = sum(is.na(gene_id)),
    missing_delta_ll = sum(is.na(delta_ll)),
    missing_wt_logprob = sum(is.na(wt_logprob)),
    missing_mut_logprob = sum(is.na(mut_logprob))
  )

cat("\n  Missing values:\n")
print(missing_summary)

# Remove genes with missing delta_ll scores
evo2_clean <- evo2_clean %>%
  filter(!is.na(delta_ll))

cat(sprintf("\n  Total genes after removing missing values: %d\n", nrow(evo2_clean)))

# Summary statistics
cat("\nDelta_ll statistics:\n")
cat(sprintf("  Mean: %.6f\n", mean(evo2_clean$delta_ll, na.rm = TRUE)))
cat(sprintf("  Median: %.6f\n", median(evo2_clean$delta_ll, na.rm = TRUE)))
cat(sprintf("  Std Dev: %.6f\n", sd(evo2_clean$delta_ll, na.rm = TRUE)))
cat(sprintf("  Min: %.6f\n", min(evo2_clean$delta_ll, na.rm = TRUE)))
cat(sprintf("  Max: %.6f\n", max(evo2_clean$delta_ll, na.rm = TRUE)))

# ============================================================================
# SAVE CLEANED DATA
# ============================================================================

cat(sprintf("\nSaving cleaned Evo2 scores to: %s\n", output_file))
write_csv(evo2_clean, output_file)

cat("=== Evo2 Scores Loading Complete ===\n")