#!/usr/bin/env Rscript
# Load and clean TnSeq essentiality data from Excel file.
#
# Reads TnSeq data, extracts gene IDs and essentiality classifications,
# and applies data cleaning (handling uncertain values, standardizing classifications).

library(readxl)
library(dplyr)

# ============================================================================
# CONFIGURATION
# ============================================================================

tnseq_excel_file <- "data/raw/mbo002173137st3.xlsx"
output_file <- "data/processed/tnseq_clean.csv"

# Create output directory if it doesn't exist
if (!dir.exists("data/processed")) {
  dir.create("data/processed", recursive = TRUE)
}

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading TnSeq data from Excel...\n")
tnseq_data <- read_excel(tnseq_excel_file, sheet = 1)
cat(sprintf("  Loaded %d genes\n", nrow(tnseq_data)))

# ============================================================================
# DATA CLEANING
# ============================================================================

cat("\nCleaning TnSeq data...\n")

# Extract gene ID and final call columns
# Column 1 = ORF ID, last column = Final Call
tnseq_clean <- tnseq_data %>%
  select(1, ncol(tnseq_data)) %>%
  setNames(c("gene_id", "final_call"))

# Clean up gene IDs and final calls (remove whitespace)
tnseq_clean <- tnseq_clean %>%
  mutate(
    gene_id = trimws(gene_id),
    final_call = trimws(as.character(final_call))
  )

# Display classification distribution
cat("  Essentiality classifications found:\n")
classification_table <- table(tnseq_clean$final_call, useNA = "ifany")
print(classification_table)

# Create binary classification
# ES = Essential, ESD = Essential Domain (treat as Essential)
# NE = Non-Essential, GA = Growth-Advantage, GD = Growth-Defect (treat as Non-Essential)
# Uncertain/other = NA (exclude)

tnseq_clean <- tnseq_clean %>%
  mutate(
    is_essential = case_when(
      final_call %in% c("ES", "ESD") ~ 1,  # Essential
      final_call %in% c("NE", "GA", "GD") ~ 0,  # Non-Essential
      TRUE ~ NA_integer_  # Uncertain - exclude
    ),
    final_call_binary = case_when(
      is_essential == 1 ~ "Essential",
      is_essential == 0 ~ "Non-Essential",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(is_essential))  # Remove genes with uncertain classifications

# Summary statistics
cat(sprintf("\n  After cleaning:\n"))
cat(sprintf("    Total genes: %d\n", nrow(tnseq_clean)))
cat(sprintf("    Essential genes (ES/ESD): %d\n", sum(tnseq_clean$is_essential == 1, na.rm = TRUE)))
cat(sprintf("    Non-Essential genes (NE/GA/GD): %d\n", sum(tnseq_clean$is_essential == 0, na.rm = TRUE)))
cat(sprintf("    Excluded (Uncertain): %d\n", nrow(tnseq_data) - nrow(tnseq_clean)))

# ============================================================================
# SAVE CLEANED DATA
# ============================================================================

cat(sprintf("\nSaving cleaned TnSeq data to: %s\n", output_file))
write.csv(tnseq_clean, output_file, row.names = FALSE)

cat("=== TnSeq Data Loading Complete ===\n")