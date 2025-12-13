#!/usr/bin/env python3
"""
Generate perturbed sequences using Evo2 paper methodology.

Insert multiple stop codons (TAATAATAATAGTGA) at 12-nucleotide intervals
throughout the coding sequence to create mutations.

Reference: Evo2 paper section 4.3.10 Gene essentiality
"""

import pandas as pd
import sys

def insert_multiple_stop_codons(sequence: str, offset: int = 12, 
                                stop_codon_pattern: str = 'TAATAATAATAGTGA') -> str:
    """
    Insert multiple stop codons at regular intervals throughout the sequence.
    
    Args:
        sequence: DNA sequence (CDS only, not including context)
        offset: Nucleotide offset between insertions (default: 12)
        stop_codon_pattern: Pattern to insert (default: TAATAATAATAGTGA)
    
    Returns:
        Mutated sequence with stop codons inserted at intervals
    """
    if len(sequence) < offset:
        # Sequence too short to insert multiple stop codons, just insert one at start
        return stop_codon_pattern + sequence
    
    # Build mutated sequence by inserting stop codon pattern at regular intervals
    mutated = []
    for i in range(0, len(sequence), offset):
        mutated.append(stop_codon_pattern)
        mutated.append(sequence[i:i+offset])
    
    mutated_seq = ''.join(mutated)
    return mutated_seq

def generate_paper_perturbed_sequences(genes_csv: str, output_file: str = 'genes_paper_perturbed.csv',
                                      offset: int = 12, 
                                      stop_codon_pattern: str = 'TAATAATAATAGTGA'):
    """
    Generate perturbed sequences using the Evo2 paper's methodology.
    
    Args:
        genes_csv: Path to genes_sequences.csv from extract_sequences.py
        output_file: Output CSV file name
        offset: Nucleotide offset between stop codon insertions (default: 12)
        stop_codon_pattern: Stop codon pattern to insert (default: TAATAATAATAGTGA)
    """
    
    # Load gene sequences
    print(f"Loading sequences from: {genes_csv}")
    genes_df = pd.read_csv(genes_csv)
    print(f"Loaded {len(genes_df)} genes")
    
    print(f"\nPerturbation strategy: Evo2 paper methodology")
    print(f"  Stop codon pattern: {stop_codon_pattern}")
    print(f"  Insertion offset: {offset} bp")
    
    # Process each gene
    results = []
    errors = []
    
    for idx, row in genes_df.iterrows():
        if idx % 500 == 0:
            print(f"Processing gene {idx}/{len(genes_df)}")
        
        gene_id = row['gene_id']
        wt_sequence = row['sequence']
        
        # Extract just the CDS portion (not the genomic context)
        try:
            cds_start = int(row['cds_start_in_seq'])
            cds_end = int(row['cds_end_in_seq'])
            cds_sequence = wt_sequence[cds_start:cds_end]
            
            # Generate mutant sequence with multiple stop codons
            mut_sequence = insert_multiple_stop_codons(cds_sequence, offset=offset, 
                                                      stop_codon_pattern=stop_codon_pattern)
            
            # Create full sequence: context + mutant CDS + context
            context_left = wt_sequence[:cds_start]
            context_right = wt_sequence[cds_end:]
            full_mut_sequence = context_left + mut_sequence + context_right
            
            results.append({
                'gene_id': gene_id,
                'locus_tag': row['locus_tag'],
                'gene_name': row['gene_name'],
                'cds_start': row['cds_start'],
                'cds_end': row['cds_end'],
                'strand': row['strand'],
                'cds_length': row['cds_length'],
                'context_bp': row['context_bp'],
                'wt_sequence': wt_sequence,
                'mut_sequence': full_mut_sequence,
                'perturbation_method': 'paper_method',
                'stop_codon_pattern': stop_codon_pattern,
                'insertion_offset': offset
            })
        
        except Exception as e:
            print(f"Error processing {gene_id}: {e}")
            errors.append({'gene_id': gene_id, 'error': str(e)})
    
    # Create output dataframe
    output_df = pd.DataFrame(results)
    
    # Save results
    output_df.to_csv(output_file, index=False)
    print(f"\nResults saved to: {output_file}")
    
    # Print summary statistics
    print(f"\nSummary:")
    print(f"  Total genes processed: {len(genes_df)}")
    print(f"  Successfully perturbed: {len(results)}")
    print(f"  Errors: {len(errors)}")
    
    if errors:
        print(f"\nFirst 5 errors:")
        for err in errors[:5]:
            print(f"  {err['gene_id']}: {err['error']}")
    
    # Show example perturbations
    print(f"\nExample perturbations (first 3 genes):")
    for idx, row in output_df.head(3).iterrows():
        wt_display = row['wt_sequence'][:60] + '...' if len(row['wt_sequence']) > 60 else row['wt_sequence']
        mut_display = row['mut_sequence'][:60] + '...' if len(row['mut_sequence']) > 60 else row['mut_sequence']
        print(f"\n  {row['gene_id']} ({row['gene_name']}):")
        print(f"    WT:  {wt_display}")
        print(f"    MUT: {mut_display}")
        print(f"    WT length: {len(row['wt_sequence'])}, MUT length: {len(row['mut_sequence'])}")
    
    return output_df

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python generate_paper_perturbed_sequences.py <genes_sequences.csv> [output_file] [offset] [stop_codon_pattern]")
        print("Example: python generate_paper_perturbed_sequences.py genes_sequences_8k.csv genes_paper_perturbed.csv 12 TAATAATAATAGTGA")
        sys.exit(1)
    
    genes_csv = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'genes_paper_perturbed.csv'
    offset = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    stop_codon_pattern = sys.argv[4] if len(sys.argv) > 4 else 'TAATAATAATAGTGA'
    
    perturbed_df = generate_paper_perturbed_sequences(genes_csv, output_file, offset, stop_codon_pattern)