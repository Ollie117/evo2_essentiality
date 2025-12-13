#!/usr/bin/env python3
"""
Score perturbed and wild-type sequences using Evo2 model.

Loads sequences (wild-type and mutant) and calculates log-likelihood scores
for each using the Evo2 7B model. Computes delta log-likelihood (Δℓ) which
serves as a proxy for fitness impact and essentiality prediction.

Output includes:
- Wild-type log-likelihood score
- Mutant log-likelihood score
- Delta log-likelihood (Δℓ = mutant - wild-type)
  Negative Δℓ suggests mutation is deleterious (gene likely essential)
  Positive Δℓ suggests mutation is tolerated (gene likely non-essential)
"""

import pandas as pdf
import numpy as np
import torch
from evo2 import Evo2
from tqdm import tqdm
import sys

def score_sequence(sequence: str, model: Evo2, device: str = 'cuda:0') -> float:
    """
    Score a single sequence using Evo2 model.
    
    Args:
        sequence: DNA sequence to score
        model: Loaded Evo2 model
        device: Device to use (cuda:0, cpu, etc.)
    
    Returns:
        Mean log-likelihood score for the sequence
    """
    try:
        # Tokenize sequence
        tokens = model.tokenizer.tokenize(sequence)
        input_ids = torch.tensor(tokens, dtype=torch.long).unsqueeze(0).to(device)
        
        # Get logits
        with torch.no_grad():
            output = model(input_ids)
            # Handle tuple output - logits is typically first element
            if isinstance(output, tuple):
                logits = output[0]
            else:
                logits = output
        
        # Convert to log probabilities
        logprobs = torch.log_softmax(logits, dim=-1)
        
        # Get log probability of actual tokens
        # Remove last position (lookahead) and first token (BOS)
        logprobs = logprobs[:, :-1]
        input_ids = input_ids[:, 1:]
        
        # Gather log probabilities for actual tokens
        seq_logprobs = torch.gather(
            logprobs,
            2,
            input_ids.unsqueeze(-1)
        ).squeeze(-1)
        
        # Calculate mean log-likelihood
        mean_logprob = seq_logprobs.mean().item()
        return mean_logprob
    
    except Exception as e:
        print(f"Error scoring sequence: {e}")
        return np.nan

def score_perturbed_sequences(perturbed_csv: str, output_file: str = 'essentiality_scores.csv',
                              batch_size: int = 1, device: str = 'cuda:0'):
    """
    Score wild-type and mutant sequences, calculate delta log-likelihood.
    
    Args:
        perturbed_csv: Path to perturbed sequences CSV (from generate_perturbed_sequences.py)
        output_file: Output CSV file name
        batch_size: Batch size for processing (currently processes one at a time)
        device: Device to use (cuda:0, cpu, etc.)
    """
    
    # Load perturbed sequences
    print(f"Loading perturbed sequences from: {perturbed_csv}")
    genes_df = pd.read_csv(perturbed_csv)
    print(f"Loaded {len(genes_df)} genes")
    
    # Load Evo2 model
    print("\nLoading Evo2 model (7B variant)...")
    device = torch.device(device if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    model = Evo2('evo2_7b')
    print("Model loaded successfully\n")
    
    # Process each gene
    results = []
    
    for idx, row in tqdm(genes_df.iterrows(), total=len(genes_df), desc="Scoring sequences"):
        gene_id = row['gene_id']
        wt_sequence = row['wt_sequence']
        mut_sequence = row['mut_sequence']
        
        try:
            # Score wild-type
            wt_score = score_sequence(wt_sequence, model, device=str(device))
            
            # Score mutant
            mut_score = score_sequence(mut_sequence, model, device=str(device))
            
            # Calculate delta log-likelihood
            delta_ll = mut_score - wt_score
            
            results.append({
                'gene_id': gene_id,
                'locus_tag': row['locus_tag'],
                'gene_name': row['gene_name'],
                'cds_length': row['cds_length'],
                'wt_logprob': wt_score,
                'mut_logprob': mut_score,
                'delta_ll': delta_ll
            })
        
        except Exception as e:
            print(f"Error processing {gene_id}: {e}")
            results.append({
                'gene_id': gene_id,
                'locus_tag': row['locus_tag'],
                'gene_name': row['gene_name'],
                'cds_length': row['cds_length'],
                'wt_logprob': np.nan,
                'mut_logprob': np.nan,
                'delta_ll': np.nan
            })
    
    # Create output dataframe
    output_df = pd.DataFrame(results)
    
    # Save results
    output_df.to_csv(output_file, index=False)
    print(f"\nResults saved to: {output_file}")
    
    # Print summary statistics
    valid_scores = output_df['delta_ll'].notna().sum()
    print(f"\nSummary Statistics:")
    print(f"  Total genes processed: {len(output_df)}")
    print(f"  Genes with valid scores: {valid_scores}")
    
    if valid_scores > 0:
        print(f"  Mean Δℓ: {output_df['delta_ll'].mean():.6f}")
        print(f"  Median Δℓ: {output_df['delta_ll'].median():.6f}")
        print(f"  Std Dev Δℓ: {output_df['delta_ll'].std():.6f}")
        print(f"  Min Δℓ: {output_df['delta_ll'].min():.6f}")
        print(f"  Max Δℓ: {output_df['delta_ll'].max():.6f}")
        
        print(f"\nScore distribution:")
        print(f"  Negative Δℓ (deleterious): {(output_df['delta_ll'] < 0).sum()} genes")
        print(f"  Positive Δℓ (tolerated): {(output_df['delta_ll'] > 0).sum()} genes")
        print(f"  Zero Δℓ: {(output_df['delta_ll'] == 0).sum()} genes")
    
    # Show example predictions
    print(f"\nExample predictions (first 10 genes):")
    print(output_df.head(10)[['gene_id', 'gene_name', 'wt_logprob', 'mut_logprob', 'delta_ll']].to_string())
    
    return output_df

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python score_evo2.py <perturbed_sequences.csv> [output_file] [device]")
        print("Example: python score_evo2.py genes_evo_strategy_peturbed.csv essentiality_scores.csv cuda:0")
        sys.exit(1)
    
    perturbed_csv = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'essentiality_scores.csv'
    device = sys.argv[3] if len(sys.argv) > 3 else 'cuda:0'
    
    scores_df = score_perturbed_sequences(perturbed_csv, output_file, device=device)