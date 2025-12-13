import pandas as pd

# Load data
df = pd.read_excel('mbo002173137st3.xlsx', sheet_name=0, header=1)

# Get ESD genes
esd_genes = df[df['Final Call'] == 'ESD'].copy()

# Check the TA site distribution
esd_genes['Essential Sites'] = esd_genes['Number of Sites Belonging to Essential State']
esd_genes['Non-Essential Sites'] = esd_genes['Number of Sites Belonging to Non-Essential State']

print("ESD genes with TA site distribution:")
print(esd_genes[['ORF ID', 'Name', 'Essential Sites', 'Non-Essential Sites']].to_string())