#!/bin/bash

Rscript load_tnseq_data.r
Rscript merge_mutation_sensitivity_data.r
Rscript roc_analysis_mutation_sensitivity.r