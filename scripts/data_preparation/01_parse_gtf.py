#!/usr/bin/env python3
"""
Parse GTF file to extract CDS (coding sequence) features for H37Rv genome.
Outputs a clean CSV with gene information and coordinates.
"""

import pandas as pd
import sys

def parse_gtf(gtf_file, output_file='genes_cds.csv'):
    """
    Parse GTF file and extract CDS features.
    
    Args:
        gtf_file: Path to GTF file
        output_file: Output CSV file name
    
    Returns:
        DataFrame with CDS features
    """
    
    # Column names for GTF format
    gtf_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 
                   'strand', 'frame', 'attribute']
    
    # Read GTF file, skipping comment lines
    print(f"Reading GTF file: {gtf_file}")
    df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, 
                     names=gtf_columns, dtype={'start': int, 'end': int})
    
    # Filter for CDS features only
    cds_df = df[df['feature'] == 'CDS'].copy()
    print(f"Found {len(cds_df)} CDS features")
    
    # Extract key information from the attribute column
    # Format: gene_id "Rv0001"; transcript_id "unassigned_transcript_1"; ... locus_tag "Rv0001"; ...
    
    cds_df['gene_id'] = cds_df['attribute'].str.extract(r'gene_id "([^"]+)"')
    cds_df['locus_tag'] = cds_df['attribute'].str.extract(r'locus_tag "([^"]+)"')
    cds_df['gene_name'] = cds_df['attribute'].str.extract(r'gene "([^"]+)"')
    cds_df['product'] = cds_df['attribute'].str.extract(r'product "([^"]+)"')
    cds_df['protein_id'] = cds_df['attribute'].str.extract(r'protein_id "([^"]+)"')
    
    # Select relevant columns
    cds_clean = cds_df[['seqname', 'gene_id', 'locus_tag', 'gene_name', 
                        'start', 'end', 'strand', 'product', 'protein_id']].copy()
    
    # Rename for clarity
    cds_clean.columns = ['chromosome', 'gene_id', 'locus_tag', 'gene_name', 
                         'cds_start', 'cds_end', 'strand', 'product', 'protein_id']
    
    # Calculate sequence length
    cds_clean['cds_length'] = cds_clean['cds_end'] - cds_clean['cds_start'] + 1
    
    # Remove duplicates (keep first occurrence of each gene)
    cds_clean = cds_clean.drop_duplicates(subset=['gene_id'], keep='first')
    
    # Sort by genomic position
    cds_clean = cds_clean.sort_values('cds_start').reset_index(drop=True)
    
    print(f"After deduplication: {len(cds_clean)} unique genes")
    
    # Save to CSV
    cds_clean.to_csv(output_file, index=False)
    print(f"Saved to: {output_file}")
    
    # Print summary statistics
    print(f"\nSummary:")
    print(f"  Total genes: {len(cds_clean)}")
    print(f"  Forward strand (+): {len(cds_clean[cds_clean['strand'] == '+'])}")
    print(f"  Reverse strand (-): {len(cds_clean[cds_clean['strand'] == '-'])}")
    print(f"  Average CDS length: {cds_clean['cds_length'].mean():.0f} bp")
    print(f"  Min/Max CDS length: {cds_clean['cds_length'].min()} / {cds_clean['cds_length'].max()} bp")
    
    # Show first few genes
    print(f"\nFirst 5 genes:")
    print(cds_clean.head())
    
    return cds_clean

if __name__ == "__main__":
    # Check command line arguments
    if len(sys.argv) < 2:
        print("Usage: python parse_gtf.py <gtf_file> [output_file]")
        print("Example: python parse_gtf.py NC_000962.3_genomic.gtf genes_cds.csv")
        sys.exit(1)
    
    gtf_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'genes_cds.csv'
    
    cds_dataframe = parse_gtf(gtf_file, output_file)
