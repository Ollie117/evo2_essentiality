#!/usr/bin/env python3
"""
Generate mutations with varying stop codon density to test model sensitivity.

Creates perturbations where stop codons are inserted at 12bp intervals,
but with different densities (1, 2, 3, or 4 stop codons per interval).

Tests: At what mutation burden does Evo2 become insensitive or saturate?
Does it show linear response to mutation burden?

Uses the Evo2 paper methodology (12bp intervals, TAATAATAATAGTGA pattern).
"""

import pandas as pd
import numpy as np
import sys

def insert_stops_at_intervals(sequence: str, offset: int = 12, num_stops_per_position: int = 1,
                             stop_codon: str = 'TAATAATAATAGTGA') -> str:
    """
    Insert stop codons at regular intervals throughout sequence.
    
    Args:
        sequence: DNA sequence
        offset: Interval between insertions (default: 12bp)
        num_stops_per_position: Number of stop codons to insert at each position
        stop_codon: Stop codon pattern
    
    Returns:
        Sequence with stops inserted at intervals
    """
    if len(sequence) < offset:
        return stop_codon * num_stops_per_position + sequence
    
    stop_pattern = stop_codon * num_stops_per_position
    mutated = []
    
    for i in range(0, len(sequence), offset):
        mutated.append(stop_pattern)
        mutated.append(sequence[i:i+offset])
    
    return ''.join(mutated)

def generate_mutation_sensitivity_perturbations(genes_csv: str, output_file: str = 'genes_mutation_sensitivity_perturbed.csv',
                                               stop_codon: str = 'TAATAATAATAGTGA'):
    """
    Generate mutations with varying stop codon density to test model sensitivity.
    
    For each gene, create four perturbations with stops inserted every 12bp:
    - Sensitivity 1: 1 stop codon at each 12bp interval (mild)
    - Sensitivity 2: 2 stop codons at each 12bp interval (moderate)
    - Sensitivity 3: 3 stop codons at each 12bp interval (strong, like paper)
    - Sensitivity 4: 4 stop codons at each 12bp interval (extreme)
    
    Tests: At what mutation burden does Evo2 become insensitive or saturate?
    
    Args:
        original_csv: Path to original genes_sequences_8k.csv
        perturbed_csv: Path to perturbed sequences CSV (with wt_sequence)
        scores_csv: Path to essentiality scores CSV
        output_prefix: Prefix for output files
    """
    
    print("=" * 80)
    print("MUTATION SENSITIVITY PERTURBATION GENERATION")
    print("=" * 80)
    
    # Load gene sequences
    print(f"\n[1/2] Loading gene sequences from: {genes_csv}")
    genes_df = pd.read_csv(genes_csv)
    print(f"      Loaded {len(genes_df)} genes")
    
    # Generate perturbations
    print(f"\n[2/2] Generating mutation sensitivity perturbations...")
    results = []
    
    for idx, row in genes_df.iterrows():
        if idx % 500 == 0:
            print(f"      Processing gene {idx}/{len(genes_df)}")
        
        gene_id = row['gene_id']
        wt_sequence = row['sequence']
        
        try:
            cds_start = int(row['cds_start_in_seq'])
            cds_end = int(row['cds_end_in_seq'])
            cds_sequence = wt_sequence[cds_start:cds_end]
            cds_length = cds_end - cds_start
            
            # Extract just the CDS portion
            cds_seq = wt_sequence[cds_start:cds_end]
            
            # Create mutations with different numbers of stop codons at 12bp intervals
            mut_sens1 = insert_stops_at_intervals(cds_seq, offset=12, num_stops_per_position=1, stop_codon=stop_codon)
            mut_sens2 = insert_stops_at_intervals(cds_seq, offset=12, num_stops_per_position=2, stop_codon=stop_codon)
            mut_sens3 = insert_stops_at_intervals(cds_seq, offset=12, num_stops_per_position=3, stop_codon=stop_codon)
            mut_sens4 = insert_stops_at_intervals(cds_seq, offset=12, num_stops_per_position=4, stop_codon=stop_codon)
            
            # Reconstruct full sequences with genomic context
            context_left = wt_sequence[:cds_start]
            context_right = wt_sequence[cds_end:]
            
            mut_sens1_full = context_left + mut_sens1 + context_right
            mut_sens2_full = context_left + mut_sens2 + context_right
            mut_sens3_full = context_left + mut_sens3 + context_right
            mut_sens4_full = context_left + mut_sens4 + context_right
            
            # Calculate mutation burdens
            stops_per_interval_1 = 1
            stops_per_interval_2 = 2
            stops_per_interval_3 = 3
            stops_per_interval_4 = 4
            
            num_intervals = (cds_length // 12) + 1
            
            total_stops_bp_1 = num_intervals * len(stop_codon) * stops_per_interval_1
            total_stops_bp_2 = num_intervals * len(stop_codon) * stops_per_interval_2
            total_stops_bp_3 = num_intervals * len(stop_codon) * stops_per_interval_3
            total_stops_bp_4 = num_intervals * len(stop_codon) * stops_per_interval_4
            
            results.append({
                'gene_id': gene_id,
                'locus_tag': row['locus_tag'],
                'gene_name': row['gene_name'],
                'cds_length': row['cds_length'],
                'wt_sequence': wt_sequence,
                'mut_sens1_sequence': mut_sens1_full,
                'mut_sens2_sequence': mut_sens2_full,
                'mut_sens3_sequence': mut_sens3_full,
                'mut_sens4_sequence': mut_sens4_full,
                'mutation_stops_per_interval': [1, 2, 3, 4],
                'mutation_total_stops_bp_sens1': total_stops_bp_1,
                'mutation_total_stops_bp_sens2': total_stops_bp_2,
                'mutation_total_stops_bp_sens3': total_stops_bp_3,
                'mutation_total_stops_bp_sens4': total_stops_bp_4,
                'mutation_total_stops_percent_cds_sens1': (total_stops_bp_1 / row['cds_length']) * 100,
                'mutation_total_stops_percent_cds_sens2': (total_stops_bp_2 / row['cds_length']) * 100,
                'mutation_total_stops_percent_cds_sens3': (total_stops_bp_3 / row['cds_length']) * 100,
                'mutation_total_stops_percent_cds_sens4': (total_stops_bp_4 / row['cds_length']) * 100,
                'perturbation_method': 'mutation_sensitivity',
                'perturbation_interval': 12,
                'stop_codon_pattern': stop_codon
            })
        
        except Exception as e:
            print(f"      Warning - Error processing {gene_id}: {e}")
    
    # Create output dataframe
    output_df = pd.DataFrame(results)
    
    # Save results
    output_df.to_csv(output_file, index=False)
    print(f"\n      Results saved to: {output_file}")
    
    # Print summary statistics
    print(f"\n" + "=" * 80)
    print(f"SUMMARY STATISTICS")
    print(f"=" * 80)
    print(f"Total genes processed: {len(output_df)}")
    print(f"Genes successfully perturbed: {len(output_df)}")
    
    print(f"\nMutation burden by sensitivity level (average across genes):")
    print(f"  Sensitivity 1 (1 stop/12bp): {output_df['mutation_total_stops_percent_cds_sens1'].mean():.2f}% of CDS")
    print(f"  Sensitivity 2 (2 stops/12bp): {output_df['mutation_total_stops_percent_cds_sens2'].mean():.2f}% of CDS")
    print(f"  Sensitivity 3 (3 stops/12bp): {output_df['mutation_total_stops_percent_cds_sens3'].mean():.2f}% of CDS")
    print(f"  Sensitivity 4 (4 stops/12bp): {output_df['mutation_total_stops_percent_cds_sens4'].mean():.2f}% of CDS")
    
    print(f"\nExample genes:")
    for idx, row in output_df.head(3).iterrows():
        print(f"\n  {row['gene_id']} ({row['gene_name']}):")
        print(f"    CDS length: {row['cds_length']} bp")
        print(f"    Sensitivity 1: {row['mutation_total_stops_percent_cds_sens1']:.2f}% of CDS")
        print(f"    Sensitivity 2: {row['mutation_total_stops_percent_cds_sens2']:.2f}% of CDS")
        print(f"    Sensitivity 3: {row['mutation_total_stops_percent_cds_sens3']:.2f}% of CDS")
        print(f"    Sensitivity 4: {row['mutation_total_stops_percent_cds_sens4']:.2f}% of CDS")
    
    print(f"\n" + "=" * 80)
    print(f"GENERATION COMPLETE")
    print(f"=" * 80 + "\n")
    
    return output_df

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python generate_mutation_sensitivity_perturbations.py <genes_sequences_8k.csv> [output_file] [stop_codon]")
        print("Example: python generate_mutation_sensitivity_perturbations.py genes_sequences_8k.csv genes_mutation_sensitivity_perturbed.csv TAATAATAATAGTGA")
        sys.exit(1)
    
    genes_csv = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'genes_mutation_sensitivity_perturbed.csv'
    stop_codon = sys.argv[3] if len(sys.argv) > 3 else 'TAATAATAATAGTGA'
    
    perturbed_df = generate_mutation_sensitivity_perturbations(genes_csv, output_file, stop_codon)