import pandas as pd

df = pd.read_csv('genes_domain_based_pertubed.csv')

# Count stop codons in each mutation
df['stops_in_domain'] = df['mut_domain_sequence'].str.count('TAATAATAATAGTGA')
df['stops_in_nondomain'] = df['mut_nondomain_sequence'].str.count('TAATAATAATAGTGA')

print(f"Mean stops in domain mutations: {df['stops_in_domain'].mean():.2f}")
print(f"Mean stops in non-domain mutations: {df['stops_in_nondomain'].mean():.2f}")
print(f"Ratio (non-domain / domain): {(df['stops_in_nondomain'].mean() / df['stops_in_domain'].mean()):.2f}x")