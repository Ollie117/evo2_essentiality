#!/usr/bin/env python3
"""
Generate domain-based and non-domain-based perturbations.

For genes with Domain [FT] annotations, create two separate perturbed versions:
1. Stop codons inserted only in domain regions
2. Stop codons inserted only in non-domain regions

Uses the Evo2 paper methodology (12bp offset, TAATAATAATAGTGA pattern).
"""

import pandas as pd
import re
import sys

def parse_domain_regions(domain_ft_str: str) -> list:
    """
    Parse Domain [FT] string to extract domain coordinates.
    
    Format: "DOMAIN 16..428; /note="..."; DOMAIN 898..1182; /note="...";"
    
    Returns:
        List of tuples: [(start, end), (start, end), ...]
    """
    if pd.isna(domain_ft_str) or domain_ft_str == '':
        return []
    
    domains = []
    # Find all DOMAIN entries with coordinates
    pattern = r'DOMAIN\s+(\d+)\.\.(\d+);'
    matches = re.findall(pattern, str(domain_ft_str))
    
    for match in matches:
        start = int(match[0]) - 1  # Convert to 0-indexed
        end = int(match[1])
        domains.append((start, end))
    
    return domains

def get_domain_mask(cds_length: int, domain_regions: list) -> list:
    """
    Create a boolean mask indicating which positions are in domains.
    
    Args:
        cds_length: Length of CDS
        domain_regions: List of (start, end) tuples for domains
    
    Returns:
        Boolean list where True = in domain, False = not in domain
    """
    mask = [False] * cds_length
    
    for start, end in domain_regions:
        # Clamp to CDS bounds
        start = max(0, min(start, cds_length - 1))
        end = max(0, min(end, cds_length))
        mask[start:end] = [True] * (end - start)
    
    return mask

def insert_stops_in_region(sequence: str, region_mask: list, offset: int = 12,
                          stop_codon: str = 'TAATAATAATAGTGA') -> str:
    """
    Insert stop codons only at positions matching the region mask.
    
    Args:
        sequence: DNA sequence
        region_mask: Boolean list indicating target positions
        offset: Insertion interval (12bp)
        stop_codon: Pattern to insert
    
    Returns:
        Mutated sequence with stops inserted in specified regions
    """
    if len(region_mask) != len(sequence):
        raise ValueError("Region mask length must match sequence length")
    
    # Find all positions where we should insert (every 12bp within the region)
    insert_positions = []
    for i in range(0, len(sequence), offset):
        # Check if this position is in the target region
        if i < len(region_mask) and region_mask[i]:
            insert_positions.append(i)
    
    # If no positions to insert, return original
    if not insert_positions:
        return sequence
    
    # Build mutated sequence by inserting stops at selected positions
    mutated = []
    last_pos = 0
    
    for pos in insert_positions:
        mutated.append(stop_codon)
        mutated.append(sequence[last_pos:pos])
        last_pos = pos
    
    mutated.append(sequence[last_pos:])
    return ''.join(mutated)

def generate_domain_nondomain_perturbations(genes_csv: str, uniprot_tsv: str,
                                           output_file: str = 'genes_domain_nondomain_perturbed.csv',
                                           offset: int = 12,
                                           stop_codon: str = 'TAATAATAATAGTGA'):
    """
    Generate domain-based and non-domain-based perturbations.
    
    Args:
        genes_csv: Path to genes_sequences.csv
        uniprot_tsv: Path to uniprot_tb_domains.tsv
        output_file: Output CSV file name
        offset: Nucleotide offset for stop insertions (default: 12)
        stop_codon: Stop codon pattern (default: TAATAATAATAGTGA)
    """
    
    print("=" * 80)
    print("DOMAIN-BASED AND NON-DOMAIN-BASED PERTURBATION GENERATION")
    print("=" * 80)
    
    # Load gene sequences
    print(f"\n[1/3] Loading gene sequences from: {genes_csv}")
    genes_df = pd.read_csv(genes_csv)
    print(f"      Loaded {len(genes_df)} genes")
    
    # Load UniProt domain annotations
    print(f"\n[2/3] Loading UniProt annotations from: {uniprot_tsv}")
    uniprot_df = pd.read_csv(uniprot_tsv, sep='\t')
    
    # Create lookup: gene_id -> Domain [FT]
    # Map Gene Names (contains Rv IDs) to Domain [FT]
    domain_lookup = {}
    for idx, row in uniprot_df.iterrows():
        gene_names = str(row['Gene Names']) if pd.notna(row['Gene Names']) else ''
        domain_ft = row['Domain [FT]']
        
        # Gene Names might be space or semicolon separated, try to find Rv IDs
        # Format examples: "pknA MRA_0017" or "msl3 pks3 pks4 Rv1180/Rv1181"
        if 'Rv' in gene_names:
            # Extract all Rv IDs (format: Rv followed by digits, optionally with 'c' for complement)
            import re
            rv_ids = re.findall(r'Rv\d+c?', gene_names)
            for rv_id in rv_ids:
                domain_lookup[rv_id] = domain_ft
    
    print(f"      Loaded {len(domain_lookup)} UniProt entries mapped to Rv IDs")
    
    # Process each gene
    print(f"\n[3/3] Generating perturbations...")
    results = []
    skipped = []
    
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
            
            # Try to find domain annotations for this gene
            domain_ft = domain_lookup.get(gene_id, None)
            
            if pd.isna(domain_ft) or domain_ft == '':
                # No domain annotation - skip this gene
                skipped.append({
                    'gene_id': gene_id,
                    'reason': 'no_domain_annotation'
                })
                continue
            
            # Parse domain regions
            domain_regions = parse_domain_regions(domain_ft)
            
            if not domain_regions:
                # Domain field exists but couldn't parse - skip
                skipped.append({
                    'gene_id': gene_id,
                    'reason': 'failed_to_parse_domains'
                })
                continue
            
            # Create domain mask
            domain_mask = get_domain_mask(cds_length, domain_regions)
            has_domain = any(domain_mask)
            
            if not has_domain:
                # Domains parsed but all outside CDS bounds
                skipped.append({
                    'gene_id': gene_id,
                    'reason': 'domains_outside_cds'
                })
                continue
            
            # Create non-domain mask (inverse)
            nondomain_mask = [not x for x in domain_mask]
            
            # Generate domain mutation
            cds_mut_domain = insert_stops_in_region(cds_sequence, domain_mask, offset, stop_codon)
            
            # Generate non-domain mutation
            cds_mut_nondomain = insert_stops_in_region(cds_sequence, nondomain_mask, offset, stop_codon)
            
            # Create full sequences with genomic context
            context_left = wt_sequence[:cds_start]
            context_right = wt_sequence[cds_end:]
            
            mut_domain_full = context_left + cds_mut_domain + context_right
            mut_nondomain_full = context_left + cds_mut_nondomain + context_right
            
            results.append({
                'gene_id': gene_id,
                'locus_tag': row['locus_tag'],
                'gene_name': row['gene_name'],
                'cds_length': cds_length,
                'num_domains': len(domain_regions),
                'domain_coverage_bp': sum(end - start for start, end in domain_regions),
                'wt_sequence': wt_sequence,
                'mut_domain_sequence': mut_domain_full,
                'mut_nondomain_sequence': mut_nondomain_full,
                'perturbation_method': 'domain_vs_nondomain',
                'stop_codon_pattern': stop_codon,
                'insertion_offset': offset,
                'domain_ft': domain_ft
            })
        
        except Exception as e:
            skipped.append({
                'gene_id': gene_id,
                'reason': f'error: {str(e)}'
            })
    
    # Create output dataframe
    output_df = pd.DataFrame(results)
    
    # Save results
    output_df.to_csv(output_file, index=False)
    print(f"\n      Results saved to: {output_file}")
    
    # Print summary statistics
    print(f"\n" + "=" * 80)
    print(f"SUMMARY STATISTICS")
    print(f"=" * 80)
    print(f"Total genes in input: {len(genes_df)}")
    print(f"Genes with domain annotations: {len(output_df)}")
    print(f"Genes skipped: {len(skipped)}")
    print(f"Success rate: {(len(output_df) / len(genes_df)) * 100:.1f}%")
    
    if skipped:
        print(f"\nSkip reasons:")
        skip_reasons = {}
        for skip in skipped:
            reason = skip['reason']
            skip_reasons[reason] = skip_reasons.get(reason, 0) + 1
        for reason, count in sorted(skip_reasons.items(), key=lambda x: -x[1]):
            print(f"  {reason}: {count}")
    
    if len(output_df) > 0:
        print(f"\nDomain coverage statistics:")
        print(f"  Mean domains per gene: {output_df['num_domains'].mean():.2f}")
        print(f"  Mean domain coverage: {(output_df['domain_coverage_bp'].mean() / output_df['cds_length'].mean() * 100):.1f}% of CDS")
        
        print(f"\nExample genes:")
        for idx, row in output_df.head(3).iterrows():
            print(f"\n  {row['gene_id']} ({row['gene_name']}):")
            print(f"    CDS length: {row['cds_length']} bp")
            print(f"    Domains: {row['num_domains']} (coverage: {row['domain_coverage_bp']} bp)")
            print(f"    WT length: {len(row['wt_sequence'])}")
            print(f"    Domain mutation length: {len(row['mut_domain_sequence'])}")
            print(f"    Non-domain mutation length: {len(row['mut_nondomain_sequence'])}")
    else:
        print(f"\nNo genes with domain annotations found. Check gene ID mapping.")
    
    print(f"\n" + "=" * 80)
    print(f"GENERATION COMPLETE")
    print(f"=" * 80 + "\n")
    
    return output_df

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python generate_domain_nondomain_perturbations.py <genes_sequences.csv> <uniprot_tb_domains.tsv> [output_file] [offset] [stop_codon]")
        print("Example: python generate_domain_nondomain_perturbations.py genes_sequences_8k.csv uniprot_tb_domains.tsv genes_domain_nondomain_perturbed.csv 12 TAATAATAATAGTGA")
        sys.exit(1)
    
    genes_csv = sys.argv[1]
    uniprot_tsv = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else 'genes_domain_nondomain_perturbed.csv'
    offset = int(sys.argv[4]) if len(sys.argv) > 4 else 12
    stop_codon = sys.argv[5] if len(sys.argv) > 5 else 'TAATAATAATAGTGA'
    
    perturbed_df = generate_domain_nondomain_perturbations(genes_csv, uniprot_tsv, output_file, offset, stop_codon)