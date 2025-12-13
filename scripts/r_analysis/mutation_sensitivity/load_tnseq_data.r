#!/usr/bin/env Rscript
# Load Evo2 gene essentiality prediction scores for mutation sensitivity analysis.
#
# Reads the CSV file containing delta_ll scores for all four sensitivity levels
# (sens1-sens4 representing increasing mutation sizes).

library(dplyr)
library(readr)

# ============================================================================
# CONFIGURATION
# ============================================================================

evo2_scores_file <- "C:/Users/ollfo/OneDrive/Desktop/Projects/evo2/evo2_essentiality/data/processed/essentiality_mutation_sensitivity.csv"
tnseq_file <- "C:/Users/ollfo/OneDrive/Desktop/Projects/evo2/evo2_essentiality/data/processed/tnseq_clean.csv"
output_dir <- "C:/Users/ollfo/OneDrive/Desktop/Projects/evo2/evo2_essentiality/data/processed"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading Evo2 mutation sensitivity scores...\n")
evo2_sens_scores <- read_csv(evo2_scores_file, col_types = cols())
cat(sprintf("  Loaded %d genes\n", nrow(evo2_sens_scores)))

# ============================================================================
# DATA CLEANING
# ============================================================================

cat("\nCleaning Evo2 mutation sensitivity data...\n")

# Clean gene IDs (remove whitespace)
evo2_sens_clean <- evo2_sens_scores %>%
  mutate(gene_id = trimws(gene_id))

# Select relevant columns: gene_id and all delta_ll columns for sensitivity levels
evo2_sens_clean <- evo2_sens_clean %>%
  select(gene_id, locus_tag, gene_name, cds_length,
         delta_ll_sens1, delta_ll_sens2, delta_ll_sens3, delta_ll_sens4,
         starts_with("mutation_total_stops"))

cat(sprintf("\n  Selected columns:\n"))
cat(sprintf("    - gene_id, locus_tag, gene_name, cds_length\n"))
cat(sprintf("    - delta_ll_sens1 through delta_ll_sens4\n"))
cat(sprintf("    - mutation_total_stops (for reference)\n"))

# Check for missing values in key columns
missing_summary <- evo2_sens_clean %>%
  summarise(
    missing_gene_id = sum(is.na(gene_id)),
    missing_sens1 = sum(is.na(delta_ll_sens1)),
    missing_sens2 = sum(is.na(delta_ll_sens2)),
    missing_sens3 = sum(is.na(delta_ll_sens3)),
    missing_sens4 = sum(is.na(delta_ll_sens4))
  )

cat("\n  Missing values by sensitivity level:\n")
print(missing_summary)

# Remove genes with any missing delta_ll values across sensitivity levels
genes_before <- nrow(evo2_sens_clean)
evo2_sens_clean <- evo2_sens_clean %>%
  filter(!is.na(delta_ll_sens1) & 
         !is.na(delta_ll_sens2) & 
         !is.na(delta_ll_sens3) & 
         !is.na(delta_ll_sens4))

genes_removed <- genes_before - nrow(evo2_sens_clean)
cat(sprintf("\n  Genes removed due to missing values: %d\n", genes_removed))
cat(sprintf("  Total genes for analysis: %d\n", nrow(evo2_sens_clean)))

# ============================================================================
# SUMMARY STATISTICS BY SENSITIVITY LEVEL
# ============================================================================

cat("\n=== Delta_ll Statistics by Sensitivity Level ===\n\n")

for (sens_level in 1:4) {
  col_name <- paste0("delta_ll_sens", sens_level)
  data <- evo2_sens_clean[[col_name]]
  
  cat(sprintf("Sensitivity Level %d (%s):\n", sens_level, col_name))
  cat(sprintf("  Mean:     %.6f\n", mean(data, na.rm = TRUE)))
  cat(sprintf("  Median:   %.6f\n", median(data, na.rm = TRUE)))
  cat(sprintf("  Std Dev:  %.6f\n", sd(data, na.rm = TRUE)))
  cat(sprintf("  Min:      %.6f\n", min(data, na.rm = TRUE)))
  cat(sprintf("  Max:      %.6f\n", max(data, na.rm = TRUE)))
  cat("\n")
}

# ============================================================================
# SAVE CLEANED DATA
# ============================================================================

output_file <- file.path(output_dir, "evo2_mutation_sensitivity_clean.csv")
cat(sprintf("Saving cleaned Evo2 mutation sensitivity scores to: %s\n", output_file))
write_csv(evo2_sens_clean, output_file)

cat("=== Evo2 Mutation Sensitivity Scores Loading Complete ===\n")