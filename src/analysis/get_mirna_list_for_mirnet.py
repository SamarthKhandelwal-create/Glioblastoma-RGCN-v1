import pandas as pd
import os

# Paths
RESULTS_DIR = r"c:\Users\Samarth\OneDrive\Documents\AP Research\results"
DE_MIRNA_PATH = os.path.join(RESULTS_DIR, "DE_miRNA_logFC_abs_gt1.tsv")
OUTPUT_TXT_PATH = os.path.join(RESULTS_DIR, "mirna_list_for_mirnet.txt")

def main():
    if not os.path.exists(DE_MIRNA_PATH):
        print(f"Error: Could not find {DE_MIRNA_PATH}")
        return

    print(f"Reading miRNAs from {DE_MIRNA_PATH}...")
    df = pd.read_csv(DE_MIRNA_PATH, sep='\t')
    
    # Extract IDs
    # miRNet usually accepts 'hsa-mir-21' or 'hsa-miR-21-5p'
    # Our IDs are likely 'hsa-mir-XXX' or 'hsa-let-XXX'
    # We will output them as-is (clean string)
    
    mirnas = df['mirna_id'].dropna().astype(str).unique().tolist()
    
    # Clean/Sort
    mirnas = sorted([m.strip() for m in mirnas])
    
    print(f"Found {len(mirnas)} unique miRNAs.")
    
    with open(OUTPUT_TXT_PATH, 'w') as f:
        for m in mirnas:
            f.write(m + "\n")
            
    print(f"Saved list to {OUTPUT_TXT_PATH}")
    print("You can upload this file directly to https://www.mirnet.ca/")

if __name__ == "__main__":
    main()
