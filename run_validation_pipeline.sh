#!/bin/bash
echo "Running Evo2 validation pipeline..."
Rscript 01_load_tnseq_data.r && \
Rscript 02_load_likelyhood_scores.r && \
Rscript 03_merge_datasets.r && \
Rscript 04_roc_analysis.r && \
echo "Pipeline complete!"
