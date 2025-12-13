#!/usr/bin/env python3
"""
Generate perturbed sequences by inserting stop codons at codon position 10.

Takes wild-type gene sequences and creates mutant versions with stop codon insertions.
Output includes both wild-type and mutant sequences for downstream analysis.
"""

import pandas as pd
import sys

def insert_stop_codon(sequence: str, codon_position: int = 10, stop_codon: str = 'TAG') -> str:
    """
    Insert a stop codon at a specific codon position.
    
    Args:
        sequence: DNA sequence (must be multiple of 3 for valid codons)
        codon_position: Which codon to replace (1-indexed)
        stop_codon: Stop codon to insert (TAG, TAA, or TGA)
    
    Returns:
        Mutated sequence with stop codon insertion
        
    Raises:
        ValueError: If sequence length is not divisible by 3 or codon position is invalid
    """
    if len(sequence) % 3 != 0:
        raise ValueError(f"Sequence length {len(sequence)} is not divisible by 3")
    
    if codon_position < 1:
        raise ValueError("Codon position must be >= 1")
    
    if codon_position > len(sequence) // 3:
        raise ValueError(f"Codon position {codon_position} exceeds sequence length ({len(sequence) // 3} codons)")
    
    # Convert to 0-indexed position in the sequence
    start_pos = (codon_position - 1) * 3
    end_pos = start_pos + 3
    
    # Replace the codon with stop codon
    mutated = sequence[:start_pos] + stop_codon + sequence[end_pos:]
    
    return mutated

def generate_perturbed_sequences(genes_csv: str, output_file: str = 'genes_perturbed.csv',
                                codon_position: int = 10, stop_codon: str = 'TAG'):
    """
    Generate perturbed sequences with stop codon insertions.
    
    Args:
        genes_csv: Path to genes_sequences.csv from extract_sequences.py
        output_file: Output CSV file name
        codon_position: Position to insert stop codon (1-indexed)
        stop_codon: Stop codon to use (TAG, TAA, or TGA)
    """
    
    # Load gene sequences
    print(f"Loading sequences from: {genes_csv}")
    genes_df = pd.read_csv(genes_csv)
    print(f"Loaded {len(genes_df)} genes")
    
    # Process each gene
    results = []
    errors = []
    
    for idx, row in genes_df.iterrows():
        if idx % 500 == 0:
            print(f"Processing gene {idx}/{len(genes_df)}")
        
        gene_id = row['gene_id']
        wt_sequence = row['sequence']
        
        try:
            # Generate mutant sequence
            mut_sequence = insert_stop_codon(wt_sequence, codon_position=codon_position, 
                                            stop_codon=stop_codon)
            
            results.append({
                'gene_id': gene_id,
                'locus_tag': row['locus_tag'],
                'gene_name': row['gene_name'],
                'cds_start': row['cds_start'],
                'cds_end': row['cds_end'],
                'strand': row['strand'],
                'cds_length': row['cds_length'],
                'wt_sequence': wt_sequence,
                'mut_sequence': mut_sequence,
                'stop_codon_position': codon_position
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
        wt_display = row['wt_sequence'][:30] + '...' if len(row['wt_sequence']) > 30 else row['wt_sequence']
        mut_display = row['mut_sequence'][:30] + '...' if len(row['mut_sequence']) > 30 else row['mut_sequence']
        print(f"\n  {row['gene_id']} ({row['gene_name']}):")
        print(f"    WT:  {wt_display}")
        print(f"    MUT: {mut_display}")
    
    return output_df

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python generate_perturbed_sequences.py <genes_sequences.csv> [output_file] [codon_position] [stop_codon]")
        print("Example: python generate_perturbed_sequences.py genes_sequences.csv genes_perturbed.csv 10 TAG")
        sys.exit(1)
    
    genes_csv = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'genes_perturbed.csv'
    codon_position = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    stop_codon = sys.argv[4] if len(sys.argv) > 4 else 'TAG'
    
    perturbed_df = generate_perturbed_sequences(genes_csv, output_file, codon_position, stop_codon)