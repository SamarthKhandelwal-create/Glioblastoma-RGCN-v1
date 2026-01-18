import pandas as pd
import os

# Paths
RESULTS_DIR = r"c:\Users\Samarth\OneDrive\Documents\AP Research\results"
MIRNET_PATH = os.path.join(RESULTS_DIR, "mirnet_interactions.csv")
INTERACTIONS_PATH = os.path.join(RESULTS_DIR, "interactions.csv")

DE_MRNA_PATH = os.path.join(RESULTS_DIR, "DE_mRNA_logFC_abs_gt1.tsv")
DE_LNCRNA_PATH = os.path.join(RESULTS_DIR, "DE_lncRNA_logFC_abs_gt1.tsv")
DE_MIRNA_PATH = os.path.join(RESULTS_DIR, "DE_miRNA_logFC_abs_gt1.tsv")

def main():
    print("--- Integrating miRNet ---")
    
    if not os.path.exists(MIRNET_PATH):
        print(f"Error: {MIRNET_PATH} not found.")
        return
        
    df = pd.read_csv(MIRNET_PATH)
    print(f"Loaded miRNet data: {len(df)} rows")
    # Expected cols: 'Source', 'Target', 'Type'?
    # Usually Source=miRNA, Target=Gene
    
    # 2. Load Mappings
    de_mrna = pd.read_csv(DE_MRNA_PATH, sep='\t')
    de_lncrna = pd.read_csv(DE_LNCRNA_PATH, sep='\t')
    de_mirna = pd.read_csv(DE_MIRNA_PATH, sep='\t')
    
    mrna_map = dict(zip(de_mrna['gene_symbol'], de_mrna['ensembl_id']))
    lncrna_map = dict(zip(de_lncrna['gene_symbol'], de_lncrna['ensembl_id']))
    
    # Normalize our miRNAs for matching
    # Our DE: 'hsa-mir-21'
    # miRNet: 'hsa-mir-21' (usually good match if we used the generated list)
    my_mirnas = set(de_mirna['mirna_id'])
    
    new_edges = []
    
    col_source = 'Source' 
    col_target = 'Target'
    # Fallback if names differ
    if col_source not in df.columns:
        col_source = df.columns[0]
        col_target = df.columns[1]
        print(f"Guessing columns: Source={col_source}, Target={col_target}")
        
    for _, row in df.iterrows():
        src = str(row[col_source]).strip()
        dst = str(row[col_target]).strip()
        
        # Check source (miRNA)
        # Try raw match
        mir_id = None
        if src in my_mirnas:
            mir_id = src
        # Try lower match
        elif src.lower() in {m.lower() for m in my_mirnas}:
            # Recover original case? simpler to just map
            pass # TODO: Robust map if needed
            
        # If src matches a miRNA
        if mir_id:
            # Check Target (mRNA or lncRNA)
            gene_id = mrna_map.get(dst) or lncrna_map.get(dst)
            if gene_id:
                new_edges.append([mir_id, gene_id])
                
    print(f"Found {len(new_edges)} valid edges from miRNet.")
    
    if new_edges:
        new_df = pd.DataFrame(new_edges, columns=['gene_A', 'gene_B'])
        if os.path.exists(INTERACTIONS_PATH):
            old = pd.read_csv(INTERACTIONS_PATH)
            combined = pd.concat([old, new_df], ignore_index=True).drop_duplicates()
        else:
            combined = new_df.drop_duplicates()
            
        combined.to_csv(INTERACTIONS_PATH, index=False)
        print(f"Saved merged: {len(combined)} total rows.")

if __name__ == "__main__":
    main()
