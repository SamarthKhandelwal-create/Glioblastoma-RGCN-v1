import pandas as pd
import os

RESULTS_DIR = "results"
PREDS_PATH = os.path.join(RESULTS_DIR, "novel_predictions.csv")
DE_MRNA = os.path.join(RESULTS_DIR, "DE_mRNA_logFC_abs_gt1.tsv")
DE_LNCRNA = os.path.join(RESULTS_DIR, "DE_lncRNA_logFC_abs_gt1.tsv")

def main():
    # 1. Load Predictions
    df = pd.read_csv(PREDS_PATH)
    
    # 2. Build Mapping (Ensembl -> Symbol)
    symbol_map = {}
    
    if os.path.exists(DE_MRNA):
        m_df = pd.read_csv(DE_MRNA, sep='\t')
        for _, row in m_df.iterrows():
            symbol_map[row['ensembl_id']] = row['gene_symbol']
            
    if os.path.exists(DE_LNCRNA):
        l_df = pd.read_csv(DE_LNCRNA, sep='\t')
        for _, row in l_df.iterrows():
            symbol_map[row['ensembl_id']] = row['gene_symbol']
            
    # 3. Map Targets
    def get_symbol(ens_id):
        # Handle version stripping if needed
        # Our map keys likely have versions if loaded from that file.
        # But prediction script uses keys from node_mapping.
        # Let's try exact match first.
        if ens_id in symbol_map:
            return symbol_map[ens_id]
        
        # Try stripping version (ENSG...12 -> ENSG...)
        # But our map keys from DE file usually HAVE versions?
        # Let's check keys in map.
        # If map has versions, and ens_id has versions, it should match.
        return ens_id # Fallback
        
    df['Gene Symbol'] = df['Target'].apply(get_symbol)
    
    # Filter for unique genes (best score for each gene)
    unique_genes = df.sort_values('Score', ascending=False).drop_duplicates(subset='Target')
    
    # 4. Display
    print("--- Top 10 High-Confidence Gene Targets ---")
    print(unique_genes[['Gene Symbol', 'miRNA', 'Score']].head(10).to_string(index=False))

if __name__ == "__main__":
    main()
