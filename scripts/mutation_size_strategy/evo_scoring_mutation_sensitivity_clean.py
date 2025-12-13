#!/usr/bin/env python3
"""
Score mutation sensitivity perturbations with Evo2 model on cluster.

Tests Evo2's sensitivity to mutation burden by scoring mutations with 
1, 2, 3, or 4 stop codons inserted at 12bp intervals.

Questions: 
- Does Evo2 show linear response to mutation burden?
- At what point does it become insensitive or saturate?

Outputs delta_ll scores for all four sensitivity levels.
"""

import pandas as pd
import numpy as np
import torch
from evo2 import Evo2
from tqdm import tqdm
import sys
import gc
import warnings
warnings.filterwarnings('ignore')

def score_sequence(sequence: str, model, device) -> float:
    """Score a single sequence using Evo2 model."""
    try:
        if len(sequence) == 0 or pd.isna(sequence):
            return np.nan
            
        tokens = model.tokenizer.tokenize(sequence)
        
        if not tokens or len(tokens) == 0:
            return np.nan
        
        input_ids = torch.tensor(tokens, dtype=torch.long).unsqueeze(0).to(device)
        
        with torch.no_grad():
            output = model(input_ids)
            if isinstance(output, tuple):
                logits = output[0]
            else:
                logits = output
        
        if isinstance(logits, tuple):
            logits = logits[0]
        
        logprobs = torch.log_softmax(logits, dim=-1)
        logprobs = logprobs[:, :-1]
        input_ids_shifted = input_ids[:, 1:]
        
        seq_logprobs = torch.gather(
            logprobs, 2, input_ids_shifted.unsqueeze(-1)
        ).squeeze(-1)
        
        mean_logprob = seq_logprobs.mean().item()
        return mean_logprob
    
    except Exception as e:
        print(f"Warning - Error scoring sequence: {e}")
        return np.nan


def score_mutation_sensitivity_variants(perturbed_csv: str, output_file: str, device: str = 'cuda:0'):
    """
    Score mutations with varying stop codon densities (1-4 stops per 12bp interval).
    
    Args:
        perturbed_csv: Path to perturbed sequences CSV
        output_file: Output CSV file name
        device: Device to use (cuda:0, cpu, etc.)
    """
    
    print("=" * 80)
    print("EVO2 MUTATION SENSITIVITY PERTURBATION SCORING")
    print("=" * 80)
    
    # Load perturbed sequences
    print(f"\n[1/4] Loading mutation sensitivity perturbations...")
    print(f"      Input: {perturbed_csv}")
    genes_df = pd.read_csv(perturbed_csv)
    print(f"      Loaded {len(genes_df)} genes")
    
    # Load Evo2 model
    print(f"\n[2/4] Loading Evo2 model...")
    device = torch.device(device if torch.cuda.is_available() else 'cpu')
    print(f"      Device: {device}")
    
    model = Evo2('evo2_7b')
    print(f"      Model loaded successfully")
    
    # Process each gene
    print(f"\n[3/4] Scoring sequences...")
    results = []
    
    for idx, row in tqdm(genes_df.iterrows(), total=len(genes_df), desc="Progress"):
        gene_id = row['gene_id']
        wt_sequence = row['wt_sequence']
        mut_sens1 = row['mut_sens1_sequence']
        mut_sens2 = row['mut_sens2_sequence']
        mut_sens3 = row['mut_sens3_sequence']
        mut_sens4 = row['mut_sens4_sequence']
        
        try:
            # Score wild-type
            wt_score = score_sequence(wt_sequence, model, device)
            
            # Score mutations of different sensitivities
            sens1_score = score_sequence(mut_sens1, model, device) if pd.notna(mut_sens1) and mut_sens1 != 'N/A' else np.nan
            sens2_score = score_sequence(mut_sens2, model, device) if pd.notna(mut_sens2) and mut_sens2 != 'N/A' else np.nan
            sens3_score = score_sequence(mut_sens3, model, device) if pd.notna(mut_sens3) and mut_sens3 != 'N/A' else np.nan
            sens4_score = score_sequence(mut_sens4, model, device) if pd.notna(mut_sens4) and mut_sens4 != 'N/A' else np.nan
            
            # Calculate delta log-likelihood
            delta_ll_sens1 = sens1_score - wt_score if not np.isnan(sens1_score) else np.nan
            delta_ll_sens2 = sens2_score - wt_score if not np.isnan(sens2_score) else np.nan
            delta_ll_sens3 = sens3_score - wt_score if not np.isnan(sens3_score) else np.nan
            delta_ll_sens4 = sens4_score - wt_score if not np.isnan(sens4_score) else np.nan
            
            results.append({
                'gene_id': gene_id,
                'locus_tag': row['locus_tag'],
                'gene_name': row['gene_name'],
                'cds_length': row['cds_length'],
                'wt_logprob': wt_score,
                'mut_sens1_logprob': sens1_score,
                'mut_sens2_logprob': sens2_score,
                'mut_sens3_logprob': sens3_score,
                'mut_sens4_logprob': sens4_score,
                'delta_ll_sens1': delta_ll_sens1,
                'delta_ll_sens2': delta_ll_sens2,
                'delta_ll_sens3': delta_ll_sens3,
                'delta_ll_sens4': delta_ll_sens4,
                'mutation_total_stops_percent_cds_sens1': row['mutation_total_stops_percent_cds_sens1'],
                'mutation_total_stops_percent_cds_sens2': row['mutation_total_stops_percent_cds_sens2'],
                'mutation_total_stops_percent_cds_sens3': row['mutation_total_stops_percent_cds_sens3'],
                'mutation_total_stops_percent_cds_sens4': row['mutation_total_stops_percent_cds_sens4']
            })
        
        except Exception as e:
            print(f"      Warning - Error processing {gene_id}: {e}")
            results.append({
                'gene_id': gene_id,
                'locus_tag': row.get('locus_tag', ''),
                'gene_name': row.get('gene_name', ''),
                'cds_length': row.get('cds_length', np.nan),
                'wt_logprob': np.nan,
                'mut_sens1_logprob': np.nan,
                'mut_sens2_logprob': np.nan,
                'mut_sens3_logprob': np.nan,
                'mut_sens4_logprob': np.nan,
                'delta_ll_sens1': np.nan,
                'delta_ll_sens2': np.nan,
                'delta_ll_sens3': np.nan,
                'delta_ll_sens4': np.nan,
                'mutation_total_stops_percent_cds_sens1': row.get('mutation_total_stops_percent_cds_sens1', np.nan),
                'mutation_total_stops_percent_cds_sens2': row.get('mutation_total_stops_percent_cds_sens2', np.nan),
                'mutation_total_stops_percent_cds_sens3': row.get('mutation_total_stops_percent_cds_sens3', np.nan),
                'mutation_total_stops_percent_cds_sens4': row.get('mutation_total_stops_percent_cds_sens4', np.nan)
            })
        
        # Memory cleanup
        if idx % 100 == 0:
            gc.collect()
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
    
    # Create output dataframe
    output_df = pd.DataFrame(results)
    
    # Save results
    print(f"\n[4/4] Saving results...")
    output_df.to_csv(output_file, index=False)
    print(f"      Output: {output_file}")
    
    # Print summary statistics
    print(f"\n" + "=" * 80)
    print(f"SUMMARY STATISTICS")
    print(f"=" * 80)
    print(f"Total genes processed: {len(output_df)}")
    
    valid_sens1 = output_df['delta_ll_sens1'].notna().sum()
    valid_sens2 = output_df['delta_ll_sens2'].notna().sum()
    valid_sens3 = output_df['delta_ll_sens3'].notna().sum()
    valid_sens4 = output_df['delta_ll_sens4'].notna().sum()
    
    print(f"Genes with valid scores:")
    print(f"  Sensitivity 1: {valid_sens1}")
    print(f"  Sensitivity 2: {valid_sens2}")
    print(f"  Sensitivity 3: {valid_sens3}")
    print(f"  Sensitivity 4: {valid_sens4}")
    
    if valid_sens1 > 0:
        print(f"\n" + "-" * 80)
        print(f"DELTA LOG-LIKELIHOOD BY MUTATION BURDEN")
        print(f"-" * 80)
        
        sens1_vals = output_df['delta_ll_sens1'].dropna()
        sens2_vals = output_df['delta_ll_sens2'].dropna()
        sens3_vals = output_df['delta_ll_sens3'].dropna()
        sens4_vals = output_df['delta_ll_sens4'].dropna()
        
        print(f"\nSENSITIVITY 1 (1 stop codon per 12bp interval):")
        print(f"  Mean:   {sens1_vals.mean():.6f}")
        print(f"  Median: {sens1_vals.median():.6f}")
        print(f"  Std:    {sens1_vals.std():.6f}")
        
        print(f"\nSENSITIVITY 2 (2 stop codons per 12bp interval):")
        print(f"  Mean:   {sens2_vals.mean():.6f}")
        print(f"  Median: {sens2_vals.median():.6f}")
        print(f"  Std:    {sens2_vals.std():.6f}")
        
        print(f"\nSENSITIVITY 3 (3 stop codons per 12bp interval):")
        print(f"  Mean:   {sens3_vals.mean():.6f}")
        print(f"  Median: {sens3_vals.median():.6f}")
        print(f"  Std:    {sens3_vals.std():.6f}")
        
        print(f"\nSENSITIVITY 4 (4 stop codons per 12bp interval):")
        print(f"  Mean:   {sens4_vals.mean():.6f}")
        print(f"  Median: {sens4_vals.median():.6f}")
        print(f"  Std:    {sens4_vals.std():.6f}")
        
        # Statistical interpretation
        print(f"\n" + "-" * 80)
        print(f"INTERPRETATION")
        print(f"-" * 80)
        
        means = [sens1_vals.mean(), sens2_vals.mean(), sens3_vals.mean(), sens4_vals.mean()]
        
        # Check for monotonic relationship
        is_monotonic = (means[0] > means[1] > means[2] > means[3]) or (means[0] < means[1] < means[2] < means[3])
        
        if is_monotonic and means[0] > means[1]:
            print(f"✓ LINEAR RESPONSE: Delta_ll decreases with mutation burden")
            print(f"  Sens1: {means[0]:.6f}")
            print(f"  Sens2: {means[1]:.6f}")
            print(f"  Sens3: {means[2]:.6f}")
            print(f"  Sens4: {means[3]:.6f}")
            print(f"  → Evo2 shows sensitivity and linear response to mutation size")
        elif means[3] < means[2] and means[2] < means[1] and means[1] < means[0]:
            print(f"✓ MONOTONIC DECREASING: Consistent penalty increase with mutation burden")
            print(f"  → Evo2 recognizes larger mutations are worse")
        else:
            print(f"✗ NON-LINEAR or SATURATING RESPONSE")
            print(f"  Sens1: {means[0]:.6f}")
            print(f"  Sens2: {means[1]:.6f}")
            print(f"  Sens3: {means[2]:.6f}")
            print(f"  Sens4: {means[3]:.6f}")
            print(f"  → Evo2 may become insensitive or saturate at high mutation burden")
        
        # Spearman correlation with mutation burden
        from scipy.stats import spearmanr
        all_delta_ll = pd.concat([sens1_vals, sens2_vals, sens3_vals, sens4_vals])
        all_burdens = pd.concat([
            pd.Series([1] * len(sens1_vals)),
            pd.Series([2] * len(sens2_vals)),
            pd.Series([3] * len(sens3_vals)),
            pd.Series([4] * len(sens4_vals))
        ]).reset_index(drop=True)
        
        corr, p_value = spearmanr(all_burdens, all_delta_ll)
        print(f"\nSpearman correlation (stops per interval vs delta_ll):")
        print(f"  ρ = {corr:.4f}, p-value = {p_value:.4e}")
        if p_value < 0.05 and corr < 0:
            print(f"  Result: SIGNIFICANT negative correlation (more stops = worse fitness)")
        elif p_value < 0.05:
            print(f"  Result: SIGNIFICANT but unexpected direction")
        else:
            print(f"  Result: NOT significant")
    
    print(f"\n" + "=" * 80)
    print(f"SCORING COMPLETE")
    print(f"=" * 80 + "\n")
    
    return output_df


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python score_mutation_sensitivity_variants.py <input.csv> [output.csv] [device]")
        print("Example: python score_mutation_sensitivity_variants.py genes_mutation_sensitivity_perturbed.csv results.csv cuda:0")
        sys.exit(1)
    
    perturbed_csv = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'essentiality_scores_mutation_sensitivity.csv'
    device = sys.argv[3] if len(sys.argv) > 3 else 'cuda:0'
    
    scores_df = score_mutation_sensitivity_variants(perturbed_csv, output_file, device=device)