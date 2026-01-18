import pandas as pd
import os

# Paths
RESULTS_DIR = r"c:\Users\Samarth\OneDrive\Documents\AP Research\results"
MIRTARBASE_PATH = os.path.join(RESULTS_DIR, "mirtarbase.csv") 
# Note: User might save as .xlsx, code should handle both if likely
INTERACTIONS_PATH = os.path.join(RESULTS_DIR, "interactions.csv")

DE_MRNA_PATH = os.path.join(RESULTS_DIR, "DE_mRNA_logFC_abs_gt1.tsv")
DE_MIRNA_PATH = os.path.join(RESULTS_DIR, "DE_miRNA_logFC_abs_gt1.tsv")

def normalize_id(gene_id):
    s = str(gene_id).strip()
    if s.startswith('ENSG'):
        return s.split('.')[0]
    return s.lower()

def main():
    print("--- Integrating miRTarBase ---")
    
    # 1. Load Data
    if not os.path.exists(MIRTARBASE_PATH):
        # Try Excel
        xlsx_path = MIRTARBASE_PATH.replace('.csv', '.xlsx')
        if os.path.exists(xlsx_path):
            print(f"Loading {xlsx_path}...")
            df = pd.read_excel(xlsx_path)
        else:
            print(f"Error: {MIRTARBASE_PATH} (or .xlsx) not found. Please download it first.")
            return
    else:
        print(f"Loading {MIRTARBASE_PATH}...")
        try:
            df = pd.read_csv(MIRTARBASE_PATH, encoding='utf-8', on_bad_lines='skip')
        except:
            df = pd.read_csv(MIRTARBASE_PATH, encoding='latin1', on_bad_lines='skip')

    print(f"Total miRTarBase rows: {len(df)}")
    
    # Check Columns (Standard miRTarBase Release 9.0 headers)
    # miRNA, Target Gene, Target Gene ID, Experiments, Support Type, References, PubMed ID
    # We need: miRNA, Target Gene (Symbol usually)
    
    # 2. Load Mapping Dictionaries
    print("Loading DE lists...")
    de_mrna = pd.read_csv(DE_MRNA_PATH, sep='\t')
    de_mirna = pd.read_csv(DE_MIRNA_PATH, sep='\t')
    
    # Map Symbol -> ID (Version Stripped for safety)
    mrna_map = dict(zip(de_mrna['gene_symbol'], de_mrna['ensembl_id']))
    
    # miRNA Set (Lowercased/Normalized for matching)
    mirna_set = set(de_mirna['mirna_id'].apply(normalize_id))
    
    # 3. Filter and Extract
    new_edges = []
    
    # Identifiers in miRTarBase
    # 'miRNA' (e.g. hsa-miR-20a-5p)
    # 'Target Gene' (e.g. HSPA4)
    
    c_mir = 'miRNA'
    c_gene = 'Target Gene'
    
    # Verify cols
    if c_mir not in df.columns or c_gene not in df.columns:
        # Try to guess
        print("Standard columns not found. Available:", df.columns.tolist())
        # Add fallback logic if needed
        return

    print("Filtering interactions...")
    for _, row in df.iterrows():
        mir_raw = str(row[c_mir]).strip()
        gene_raw = str(row[c_gene]).strip()
        
        # Check miRNA (Normalize)
        # miRTarBase: hsa-miR-21-5p
        # DE: hsa-mir-21 (stem) or hsa-let-7a
        # Fuzzy match logic
        mir_norm = mir_raw.lower()
        
        # Check exact or stem match
        # Try stripping '-5p'/'-3p'?
        # A simple check: is 'mir_norm' in our set? or is 'mir_norm' contained in?
        # Actually our DE is likely stem 'hsa-mir-21'. miRTarBase is mature 'hsa-miR-21-5p'.
        # 'hsa-mir-21' stem produces 'hsa-miR-21-5p'.
        # So if our DE ID is a substring of the miRTarBase ID (ignoring case), it's a match-ish.
        # But 'hsa-mir-21' (stem) is NOT a substring of 'hsa-miR-21-5p' strictly if dash differs?
        # Let's try simple 'hsa-mir-21' vs 'hsa-mir-21-5p'.
        
        mir_id = None
        
        # List based check (slow but safe for 300 items)
        for my_mir in mirna_set:
            # my_mir: hsa-mir-21
            # mir_norm: hsa-mir-21-5p
            if my_mir in mir_norm or mir_norm in my_mir:
                mir_id = my_mir # Map to OUR format
                break
                
        if not mir_id:
            continue
            
        # Check mRNA
        gene_id = mrna_map.get(gene_raw)
        if not gene_id:
            continue
            
        new_edges.append([mir_id, gene_id])
        
    print(f"Found {len(new_edges)} valid interactions involving our DE list.")
    
    if not new_edges:
        return

    # 4. Merge
    new_df = pd.DataFrame(new_edges, columns=['gene_A', 'gene_B'])
    
    if os.path.exists(INTERACTIONS_PATH):
        old_df = pd.read_csv(INTERACTIONS_PATH)
        combined = pd.concat([old_df, new_df], ignore_index=True)
    else:
        combined = new_df
        
    before = len(combined)
    combined.drop_duplicates(inplace=True)
    after = len(combined)
    
    combined.to_csv(INTERACTIONS_PATH, index=False)
    print(f"Merged. Total interactions: {after} (Added {after - len(old_df) if os.path.exists(INTERACTIONS_PATH) else after})")

if __name__ == "__main__":
    main()
