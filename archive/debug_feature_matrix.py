#!/usr/bin/env python
"""Debug: Check what gene IDs are actually in the feature matrix."""
import pandas as pd

df = pd.read_csv('results/node_features_matrix.csv', index_col=0, nrows=100)
print("Sample of gene IDs in feature matrix (first 100):")
print(df.index.tolist())
print(f"\nTotal unique IDs: {len(df.index)}")

# Check if any match our mock genes
mock_genes = ['TP53', 'PTEN', 'EGFR', 'FOXO3', 'DICER1', 'DROSHA']
print(f"\nSearching for mock genes in feature matrix:")
for gene in mock_genes:
    matches = [idx for idx in df.index if gene.upper() in str(idx).upper()]
    print(f"  {gene}: {len(matches)} matches - {matches[:3]}")
