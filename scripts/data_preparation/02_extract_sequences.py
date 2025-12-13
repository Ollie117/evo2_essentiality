#!/usr/bin/env python3
"""
Extract DNA sequences from FNA (FASTA) file based on CDS coordinates.
Uses the genes_cds.csv output from parse_gtf.py to extract sequences for each gene.
"""

import pandas as pd
import sys

def read_fasta(fasta_file):
    """
    Read FASTA file and return sequence as a single string.
    Handles multi-line FASTA format.
    
    Args:
        fasta_file: Path to FASTA/FNA file
    
    Returns:
        Tuple of (sequence_id, full_sequence)
    """
    print(f"Reading FASTA file: {fasta_file}")
    
    sequence_id = None
    sequence = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                # Header line
                sequence_id = line[1:].split()[0]  # Get first part of header
                print(f"Found sequence: {sequence_id}")
            else:
                # Sequence line
                sequence.append(line)
    
    full_sequence = ''.join(sequence).upper()
    print(f"Total sequence length: {len(full_sequence)} bp")
    
    return sequence_id, full_sequence

def reverse_complement(seq):
    """
    Generate reverse complement of DNA sequence.
    
    Args:
        seq: DNA sequence string
    
    Returns:
        Reverse complement sequence
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))

def extract_sequences(csv_file, fasta_file, output_file='genes_sequences.csv', context_bp: int = 0):
    """
    Extract DNA sequences for each gene from FASTA file with optional genomic context.
    
    Args:
        csv_file: Path to genes_cds.csv from parse_gtf.py
        fasta_file: Path to FNA/FASTA file
        output_file: Output CSV with sequences
        context_bp: Amount of genomic context to include on each side (default: 0)
                   Total sequence length will be: context_bp + cds_length + context_bp
    
    Returns:
        DataFrame with genes and their sequences
    """
    
    # Read gene coordinates
    print(f"Reading gene coordinates from: {csv_file}")
    genes_df = pd.read_csv(csv_file)
    print(f"Loaded {len(genes_df)} genes")
    
    # Read genome sequence
    seq_id, genome_seq = read_fasta(fasta_file)
    genome_length = len(genome_seq)
    print(f"Genome length: {genome_length} bp")
    
    if context_bp > 0:
        print(f"Including {context_bp} bp of genomic context on each side")
    
    # Extract sequences for each gene
    sequences = []
    errors = []
    
    for idx, row in genes_df.iterrows():
        if idx % 500 == 0:
            print(f"Processing gene {idx}/{len(genes_df)}")
        
        gene_id = row['gene_id']
        start = row['cds_start']
        end = row['cds_end']
        strand = row['strand']
        
        try:
            # Calculate context boundaries
            context_start = max(0, start - 1 - context_bp)  # Convert to 0-based and subtract context
            context_end = min(genome_length, end + context_bp)  # end is already exclusive in 0-based
            
            # Extract sequence with context (convert to 0-based indexing)
            # GTF uses 1-based indexing, Python uses 0-based
            seq = genome_seq[context_start:context_end]
            
            # Track the actual CDS position within the extracted sequence
            cds_start_in_seq = (start - 1) - context_start
            cds_end_in_seq = cds_start_in_seq + (end - start + 1)
            
            # Reverse complement if on negative strand
            if strand == '-':
                seq = reverse_complement(seq)
                # Adjust CDS positions for reverse complement
                seq_len = len(seq)
                cds_start_in_seq, cds_end_in_seq = seq_len - cds_end_in_seq, seq_len - cds_start_in_seq
            
            sequences.append({
                'gene_id': gene_id,
                'locus_tag': row['locus_tag'],
                'gene_name': row['gene_name'],
                'cds_start': start,
                'cds_end': end,
                'strand': strand,
                'cds_length': row['cds_length'],
                'cds_start_in_seq': cds_start_in_seq,
                'cds_end_in_seq': cds_end_in_seq,
                'sequence_length': len(seq),
                'context_bp': context_bp,
                'sequence': seq
            })
        
        except Exception as e:
            errors.append({
                'gene_id': gene_id,
                'error': str(e)
            })
    
    # Create output dataframe
    output_df = pd.DataFrame(sequences)
    
    print(f"\nSuccessfully extracted: {len(sequences)} sequences")
    if errors:
        print(f"Errors encountered: {len(errors)}")
        for err in errors[:5]:  # Show first 5 errors
            print(f"  {err['gene_id']}: {err['error']}")
    
    # Save to CSV
    output_df.to_csv(output_file, index=False)
    print(f"Saved to: {output_file}")
    
    # Print summary
    print(f"\nSummary:")
    print(f"  Total genes: {len(output_df)}")
    print(f"  Average CDS length: {output_df['cds_length'].mean():.0f} bp")
    print(f"  Average total sequence length (with {context_bp} bp context): {output_df['sequence_length'].mean():.0f} bp")
    print(f"  Min/Max CDS length: {output_df['cds_length'].min()} / {output_df['cds_length'].max()} bp")
    print(f"  Min/Max total sequence length: {output_df['sequence_length'].min()} / {output_df['sequence_length'].max()} bp")
    
    # Show first few sequences (truncated)
    print(f"\nFirst 3 genes:")
    for idx, row in output_df.head(3).iterrows():
        seq_display = row['sequence'][:50] + '...' if len(row['sequence']) > 50 else row['sequence']
        print(f"  {row['gene_id']}: {seq_display}")
    
    return output_df

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python extract_sequences.py <genes_cds.csv> <fasta_file> [output_file] [context_bp]")
        print("Example: python extract_sequences.py genes_cds.csv NC_000962.3_genomic.fna genes_sequences.csv 4096")
        print("         (4096 bp context = 8192 bp total with 4096 on each side)")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else 'genes_sequences.csv'
    context_bp = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    
    sequences_df = extract_sequences(csv_file, fasta_file, output_file, context_bp=context_bp)