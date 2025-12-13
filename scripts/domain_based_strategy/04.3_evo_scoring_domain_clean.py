#!/usr/bin/env python3
"""
Score domain-based and non-domain-based perturbations with Evo2 model on cluster.

Outputs delta_ll scores for both domain and non-domain mutations to compare
whether Evo2 recognizes that domain disruptions are more harmful.
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


def score_domain_nondomain_variants(perturbed_csv: str, output_file: str, device: str = 'cuda:0'):
    """
    Score wild-type and domain/non-domain mutant sequences.
    
    Args:
        perturbed_csv: Path to perturbed sequences CSV
        output_file: Output CSV file name
        device: Device to use (cuda:0, cpu, etc.)
    """
    
    print("=" * 80)
    print("EVO2 DOMAIN VS NON-DOMAIN PERTURBATION SCORING")
    print("=" * 80)
    
    # Load perturbed sequences
    print(f"\n[1/4] Loading domain vs non-domain perturbations...")
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
        mut_domain = row['mut_domain_sequence']
        mut_nondomain = row['mut_nondomain_sequence']
        
        try:
            # Score wild-type
            wt_score = score_sequence(wt_sequence, model, device)
            
            # Score domain mutation
            domain_score = score_sequence(mut_domain, model, device) if pd.notna(mut_domain) and mut_domain != 'N/A' else np.nan
            
            # Score non-domain mutation
            nondomain_score = score_sequence(mut_nondomain, model, device) if pd.notna(mut_nondomain) and mut_nondomain != 'N/A' else np.nan
            
            # Calculate delta log-likelihood
            delta_ll_domain = domain_score - wt_score if not np.isnan(domain_score) else np.nan
            delta_ll_nondomain = nondomain_score - wt_score if not np.isnan(nondomain_score) else np.nan
            
            results.append({
                'gene_id': gene_id,
                'locus_tag': row['locus_tag'],
                'gene_name': row['gene_name'],
                'cds_length': row['cds_length'],
                'num_domains': row['num_domains'],
                'domain_coverage_bp': row['domain_coverage_bp'],
                'wt_logprob': wt_score,
                'mut_domain_logprob': domain_score,
                'mut_nondomain_logprob': nondomain_score,
                'delta_ll_domain': delta_ll_domain,
                'delta_ll_nondomain': delta_ll_nondomain,
                'delta_ll_difference': delta_ll_domain - delta_ll_nondomain if not (np.isnan(delta_ll_domain) or np.isnan(delta_ll_nondomain)) else np.nan
            })
        
        except Exception as e:
            print(f"      Warning - Error processing {gene_id}: {e}")
            results.append({
                'gene_id': gene_id,
                'locus_tag': row.get('locus_tag', ''),
                'gene_name': row.get('gene_name', ''),
                'cds_length': row.get('cds_length', np.nan),
                'num_domains': row.get('num_domains', np.nan),
                'domain_coverage_bp': row.get('domain_coverage_bp', np.nan),
                'wt_logprob': np.nan,
                'mut_domain_logprob': np.nan,
                'mut_nondomain_logprob': np.nan,
                'delta_ll_domain': np.nan,
                'delta_ll_nondomain': np.nan,
                'delta_ll_difference': np.nan
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
    
    valid_domain = output_df['delta_ll_domain'].notna().sum()
    valid_nondomain = output_df['delta_ll_nondomain'].notna().sum()
    valid_both = ((output_df['delta_ll_domain'].notna()) & (output_df['delta_ll_nondomain'].notna())).sum()
    
    print(f"Genes with valid domain mutations: {valid_domain}")
    print(f"Genes with valid non-domain mutations: {valid_nondomain}")
    print(f"Genes with valid BOTH mutations: {valid_both}")
    
    if valid_both > 0:
        print(f"\n" + "-" * 80)
        print(f"DELTA LOG-LIKELIHOOD BY MUTATION TYPE")
        print(f"-" * 80)
        
        domain_vals = output_df['delta_ll_domain'].dropna()
        nondomain_vals = output_df['delta_ll_nondomain'].dropna()
        difference_vals = output_df['delta_ll_difference'].dropna()
        
        print(f"\nDOMAIN MUTATIONS (disrupting conserved regions):")
        print(f"  Mean:   {domain_vals.mean():.6f}")
        print(f"  Median: {domain_vals.median():.6f}")
        print(f"  Std:    {domain_vals.std():.6f}")
        print(f"  Min/Max: {domain_vals.min():.6f} / {domain_vals.max():.6f}")
        
        print(f"\nNON-DOMAIN MUTATIONS (disrupting flexible regions):")
        print(f"  Mean:   {nondomain_vals.mean():.6f}")
        print(f"  Median: {nondomain_vals.median():.6f}")
        print(f"  Std:    {nondomain_vals.std():.6f}")
        print(f"  Min/Max: {nondomain_vals.min():.6f} / {nondomain_vals.max():.6f}")
        
        print(f"\nDIFFERENCE (Domain - Non-Domain):")
        print(f"  Mean:   {difference_vals.mean():.6f}")
        print(f"  Median: {difference_vals.median():.6f}")
        print(f"  Std:    {difference_vals.std():.6f}")
        
        # Statistical interpretation
        print(f"\n" + "-" * 80)
        print(f"INTERPRETATION")
        print(f"-" * 80)
        
        if domain_vals.mean() < nondomain_vals.mean():
            print(f"✓ Domain mutations are MORE DAMAGING than non-domain mutations")
            print(f"  → Evo2 appears to recognize functional importance of domains")
            print(f"  → Difference: {abs(domain_vals.mean() - nondomain_vals.mean()):.6f}")
        else:
            print(f"✗ Domain mutations are LESS/EQUALLY damaging than non-domain mutations")
            print(f"  → May indicate Evo2 is pattern-matching rather than understanding biology")
            print(f"  → Difference: {abs(domain_vals.mean() - nondomain_vals.mean()):.6f}")
        
        # T-test
        from scipy import stats
        t_stat, p_value = stats.ttest_rel(domain_vals, nondomain_vals)
        print(f"\nPaired t-test (Domain vs Non-Domain):")
        print(f"  t-statistic: {t_stat:.4f}")
        print(f"  p-value: {p_value:.4e}")
        if p_value < 0.05:
            print(f"  Result: SIGNIFICANT difference (p < 0.05)")
        else:
            print(f"  Result: NOT significant (p >= 0.05)")
    
    print(f"\n" + "=" * 80)
    print(f"SCORING COMPLETE")
    print(f"=" * 80 + "\n")
    
    return output_df


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python score_domain_nondomain_variants.py <input.csv> [output.csv] [device]")
        print("Example: python score_domain_nondomain_variants.py genes_domain_nondomain_perturbed.csv results.csv cuda:0")
        sys.exit(1)
    
    perturbed_csv = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'essentiality_scores_domain_nondomain.csv'
    device = sys.argv[3] if len(sys.argv) > 3 else 'cuda:0'
    
    scores_df = score_domain_nondomain_variants(perturbed_csv, output_file, device=device)