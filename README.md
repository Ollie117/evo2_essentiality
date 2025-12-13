# Evo2 Gene Essentiality Prediction



## Overview



## Key Results

- **Main finding 1**: value/statistic
- **Main finding 2**: value/statistic
- **Main finding 3**: value/statistic

## Installation

### Requirements
- Python 3.8+
- R 4.0+
- NVIDIA GPU (optional, for faster Evo2 inference)

### Setup

```bash
# Python dependencies
pip install -r scripts/model_setup/requirements.txt

# R dependencies
install.packages(c("readxl", "dplyr", "tidyr", "pROC", "ggplot2", "readr"))
```
 
## Usage

### Quick Start

```bash
# Phase 1: Prepare sequences and score with Evo2
python scripts/model_setup/01_parse_gtf.py
python scripts/model_setup/02_extract_sequences.py
python scripts/model_setup/03.1_generate_peturbated_sequences_evo_paper_strategy.py
python scripts/model_setup/04_evo_scoring.py

# Phase 2: Validate against experimental data
bash run_validation_pipeline.sh
```

### Input Data

- Genome: `data/raw/GCF_000195955.2_ASM19595v2_genomic.fna`
- Annotations: `data/raw/genomic.gtf`
- Experimental: `data/raw/mbo002173137st3.xlsx` (TnSeq essentiality data)

### Output

Results saved to `results/validation/`:
- `roc_curve.png` - ROC curve plot
- `performance_metrics.csv` - AUC, sensitivity, specificity, etc.
- `merged_predictions.csv` - Gene-level predictions

## Project Structure

```
evo2-essentiality/
├── data/
│   ├── raw/              # Original genomic and experimental data
│   └── processed/        # Intermediate and final outputs
├── scripts/
│   ├── model_setup/      # Sequence preparation & Evo2 scoring
│   └── analysis/         # Validation analysis & ROC curves
├── results/validation/   # Final results and plots
├── slurm/               # HPC job submission scripts
└── README.md
```

## Methodology

### Sequence Preparation
- Extract CDS from GTF annotations (3,906 genes)
- Add 8 kb genomic context (4 kb on each side)
- Generate perturbed sequences with multi-stop-codon insertions

### Evo2 Scoring
- Model: Evo2 7B (pretrained DNA language model)
- Metric: Delta log-likelihood (Δℓ) = log P(mutant) - log P(wild-type)
- Interpretation: negative Δℓ = essential, positive Δℓ = non-essential

### Validation
- Ground truth: TnSeq essentiality classifications (Sassetti et al. 2013)
- Evaluation: ROC analysis with AUC metric
- Optimal threshold: Youden index

## Results

| Metric | Value |
|--------|-------|
| AUC | 0.8246 |
| Sensitivity | 81.0% |
| Specificity | 72.1% |
| NPV | 96.2% |
| Accuracy | 73.2% |

## File Descriptions

### Python Scripts (`scripts/model_setup/`)
- `01_parse_gtf.py` - Extract CDS from GTF
- `02_extract_sequences.py` - Add genomic context
- `03.1_generate_peturbated_sequences_evo_paper_strategy.py` - Generate perturbed sequences
- `04_evo_scoring.py` - Score with Evo2

### R Scripts (`scripts/analysis/`)
- `01_load_tnseq_data.r` - Load and clean TnSeq data
- `02_load_likelyhood_scores.r` - Load Evo2 scores
- `03_merge_datasets.r` - Merge predictions with experimental data
- `04_roc_analysis.r` - Calculate ROC metrics
- `05_plot_roc_curve.r` - Generate ROC plot

## Technical Details

- **Organism**: *Mycobacterium tuberculosis* H37Rv
- **Genome**: GCF_000195955.2 (4.4 Mb, 3,906 genes)
- **Model**: Evo2 7B parameters
- **Runtime**: ~6 hours for full dataset on NVIDIA L40S

## Future Work

- [ ] Stratified analysis by functional categories
- [ ] Comparison with other prediction methods
- [ ] Conditional essentiality analysis
- [ ] Integration with pathway databases

## References

- Sassetti, C.M. & Rubin, E.J. (2013). "Mycobacterial Persistence and Pathogenesis through the Lens of Bacterial Genetics." *Nature Reviews Microbiology*.
- Arc Institute & Stanford. (2024). "Evo2: Evolutionary-inspired Large Language Models for Genomics."

## License

This project is provided for academic and research purposes.