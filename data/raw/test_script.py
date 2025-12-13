import pandas as pd

df = pd.read_csv('uniprot_tb_domains.tsv', sep='\t')

# Check for "Linker" in the Region column
has_linker = df['Region'].notna() & df['Region'].str.contains('Linker', case=False, na=False)
has_domain = df['Domain [FT]'].notna() & (df['Domain [FT]'] != '')

linker_count = has_linker.sum()
domain_count = has_domain.sum()
both_count = (has_linker & has_domain).sum()
total_genes = len(df)

print(f"Total genes: {total_genes}")
print(f"Genes with Domain [FT]: {domain_count} ({(domain_count/total_genes)*100:.1f}%)")
print(f"Genes with Linker in Region: {linker_count} ({(linker_count/total_genes)*100:.1f}%)")
print(f"Genes with BOTH Domain [FT] and Linker: {both_count} ({(both_count/total_genes)*100:.1f}%)")