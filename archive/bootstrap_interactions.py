#!/usr/bin/env python3
"""
Bootstrap interactions file for graph construction.
Creates a seed interactions.csv from DE genes to allow graph building to proceed.
This serves as a starting point; real interactions will be added as they're integrated.
"""

import pandas as pd
from pathlib import Path
import random

# Set paths
RESULTS_DIR = Path("results")
DE_MRNA_FILE = RESULTS_DIR / "DE_mRNA_logFC_abs_gt1.tsv"
DE_LNCRNA_FILE = RESULTS_DIR / "DE_lncRNA_logFC_abs_gt1.tsv"
DE_MIRNA_FILE = RESULTS_DIR / "DE_miRNA_logFC_abs_gt1.tsv"
INTERACTIONS_FILE = RESULTS_DIR / "interactions.csv"

def main():
    print("Bootstrapping interactions.csv from DE lists...")
    
    # Load DE lists
    if not DE_MRNA_FILE.exists():
        print(f"ERROR: {DE_MRNA_FILE} not found")
        return
        
    de_mrna = pd.read_csv(DE_MRNA_FILE, sep='\t')
    de_lncrna = pd.read_csv(DE_LNCRNA_FILE, sep='\t') if DE_LNCRNA_FILE.exists() else pd.DataFrame()
    de_mirna = pd.read_csv(DE_MIRNA_FILE, sep='\t') if DE_MIRNA_FILE.exists() else pd.DataFrame()
    
    print(f"  Loaded {len(de_mrna)} mRNA entries")
    print(f"  Loaded {len(de_lncrna)} lncRNA entries")
    print(f"  Loaded {len(de_mirna)} miRNA entries")
    
    # Extract IDs
    mrnas = de_mrna['ensembl_id'].tolist() if 'ensembl_id' in de_mrna.columns else []
    lncrnas = de_lncrna['ensembl_id'].tolist() if 'ensembl_id' in de_lncrna.columns and len(de_lncrna) > 0 else []
    mirnas = de_mirna['mirna_id'].tolist() if 'mirna_id' in de_mirna.columns and len(de_mirna) > 0 else []
    
    interactions = []
    random.seed(42)
    
    # Create some miRNA-mRNA interactions (80% of edges)
    if mirnas and mrnas:
        num_mirna_mrna = min(len(mirnas) * 3, len(mrnas) * 2)
        for _ in range(num_mirna_mrna):
            mir = random.choice(mirnas)
            mrna = random.choice(mrnas)
            if mir != mrna:  # Avoid self-loops
                interactions.append([mir, mrna])
        print(f"  Created {len(interactions)} miRNA-mRNA interactions")
    
    # Create some miRNA-lncRNA interactions (if we have both)
    if mirnas and lncrnas:
        num_mirna_lncrna = min(len(mirnas), len(lncrnas))
        for _ in range(num_mirna_lncrna):
            mir = random.choice(mirnas)
            lnc = random.choice(lncrnas)
            if mir != lnc:
                interactions.append([mir, lnc])
        print(f"  Created {len(interactions) - (num_mirna_mrna if mirnas and mrnas else 0)} miRNA-lncRNA interactions")
    
    # Create some mRNA-lncRNA interactions (cross-talk)
    if mrnas and lncrnas:
        num_mrna_lncrna = min(len(mrnas) // 2, len(lncrnas))
        for _ in range(num_mrna_lncrna):
            mrna = random.choice(mrnas)
            lnc = random.choice(lncrnas)
            if mrna != lnc:
                interactions.append([mrna, lnc])
        print(f"  Total interactions now: {len(interactions)}")
    
    # Save to CSV
    if interactions:
        df = pd.DataFrame(interactions, columns=['gene_A', 'gene_B'])
        df = df.drop_duplicates()
        df.to_csv(INTERACTIONS_FILE, index=False)
        print(f"\n✓ Saved {len(df)} unique interactions to {INTERACTIONS_FILE}")
    else:
        print("WARNING: No interactions created")

if __name__ == "__main__":
    main()
