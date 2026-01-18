import pandas as pd
import os

# Paths
RESULTS_DIR = r"c:\Users\Samarth\OneDrive\Documents\AP Research\results"
M2G_NAME = "mir2gene.csv"
M2L_NAME = "mir2lnc.csv"

# Assumes script is run from root or src/preprocessing. 
# Safe approach: Check both root/data/raw and ../../data/raw
if os.path.exists(os.path.join("data", "raw", M2G_NAME)):
    DATA_DIR = os.path.join("data", "raw")
else:
    # If run from src/preprocessing
    DATA_DIR = os.path.join("..", "..", "data", "raw")

MIR2GENE_PATH = os.path.join(DATA_DIR, M2G_NAME)
MIR2LNC_PATH = os.path.join(DATA_DIR, M2L_NAME)
INTERACTIONS_PATH = os.path.join(RESULTS_DIR, "interactions.csv")

DE_MRNA_PATH = os.path.join(RESULTS_DIR, "DE_mRNA_logFC_abs_gt1.tsv")
DE_LNCRNA_PATH = os.path.join(RESULTS_DIR, "DE_lncRNA_logFC_abs_gt1.tsv")
DE_MIRNA_PATH = os.path.join(RESULTS_DIR, "DE_miRNA_logFC_abs_gt1.tsv")

def normalize_mir_id(s):
    """Normalize miRNA ID to help matching."""
    return s.strip().lower()

MIN_REF_COUNT = 1 # Set to 2 or higher to filter for higher confidence

def count_refs(lit_str):
    if pd.isna(lit_str) or str(lit_str) == 'nan':
        return 0
    return len(str(lit_str).split('|'))

def main():
    print("--- Integrating ENCORI (Manual Downloads) ---")
    print(f"Filtering for interactions with >= {MIN_REF_COUNT} literature references.")
    
    # 1. Load Mappings
    print("Loading DE lists...")
    de_mrna = pd.read_csv(DE_MRNA_PATH, sep='\t')
    de_lncrna = pd.read_csv(DE_LNCRNA_PATH, sep='\t')
    de_mirna = pd.read_csv(DE_MIRNA_PATH, sep='\t')
    
    # Maps: Symbol -> Ensembl ID
    mrna_map = dict(zip(de_mrna['gene_symbol'], de_mrna['ensembl_id']))
    lncrna_map = dict(zip(de_lncrna['gene_symbol'], de_lncrna['ensembl_id']))
    
    # Set of our miRNA IDs (normalized)
    my_mirnas = set(de_mirna['mirna_id'].apply(normalize_mir_id))
    # Keep original map for retrieval if needed, but we output valid IDs
    # Actually, we should output the ID format that matches our node list (which is the one in DE_miRNA)
    # So we need a way to go from "hsa-miR-16-5p" (ENCORI) -> "hsa-mir-16-1" (DE/Ensembl)
    # This is tricky because one mature can map to multiple precursors.
    # We will map if the *mature* sequence or ID is contained.
    # Heuristic: If our DE ID (e.g. hsa-mir-21) is in the ENCORI ID (hsa-miR-21-5p), it's a match.
    # OR if the ENCORI ID (hsa-let-7d) is equal to ours.
    
    # Create a lookup list
    my_mirna_list = de_mirna['mirna_id'].tolist()
    
    def resolve_mirna(encori_id):
        encori_norm = normalize_mir_id(encori_id)
        # Exact match check
        if encori_norm in my_mirnas:
            return encori_norm # Matches our format
        
        # Substring check
        # Check if any of our DE IDs is a substring of ENCORI ID?
        # e.g. DE: hsa-mir-21, Enc: hsa-miR-21-5p -> Match
        # But caution: hsa-mir-21 vs hsa-mir-210
        for my_m in my_mirna_list:
            norm_m = normalize_mir_id(my_m)
            # Boundary safe check?
            if norm_m == encori_norm:
                return my_m
            # If our stem is in their mature (and usually their mature is longer or equal)
            # Actually, user's list often has 'hsa-mir-xxx' (stem).
            # ENCORI has 'hsa-miR-xxx-5p'.
            # If 'hsa-mir-21' is in 'hsa-mir-21-5p', it's valid.
            if norm_m in encori_norm:
                 # Verify it's not a partial suffix match like mir-21 inside mir-210
                 # Simple check: encori starts with norm_m?
                 if encori_norm.startswith(norm_m):
                     return my_m
        return None

    new_edges = []
    
    # 2. Process mir2gene (miRNA -> mRNA)
    if os.path.exists(MIR2GENE_PATH):
        print(f"Processing {MIR2GENE_PATH}...")
        df_genes = pd.read_csv(MIR2GENE_PATH, on_bad_lines='skip')
        # Cols: ID (miRNA), Target (Gene Symbol)
        count = 0
        for _, row in df_genes.iterrows():
            mir_raw = str(row['ID'])
            gene_sym = str(row['Target'])
            
            # Resolve miRNA
            my_mir = resolve_mirna(mir_raw)
            if not my_mir:
                continue
                
            # Resolve mRNA
            my_gene = mrna_map.get(gene_sym)
            if not my_gene:
                continue
            
            # Filter by Refs
            refs = count_refs(row.get('Literature'))
            if refs < MIN_REF_COUNT:
                continue
                
            new_edges.append([my_mir, my_gene])
            count += 1
            
        print(f"  Added {count} miRNA-mRNA interactions.")
        
    else:
        print("  mir2gene.csv not found.")

    # 3. Process mir2lnc (miRNA -> lncRNA)
    if os.path.exists(MIR2LNC_PATH):
        print(f"Processing {MIR2LNC_PATH}...")
        df_lnc = pd.read_csv(MIR2LNC_PATH, on_bad_lines='skip')
        count = 0
        for _, row in df_lnc.iterrows():
            mir_raw = str(row['ID'])
            target_sym = str(row['Target']) 
            
            # Resolve miRNA
            my_mir = resolve_mirna(mir_raw)
            if not my_mir:
                continue
            
            # Resolve lncRNA
            my_lnc = lncrna_map.get(target_sym)
            if not my_lnc:
                continue
                
            # Note: Graph direction.
            # My schema: lncRNA (1) -> miRNA (2) (Sequestering)
            # ENCORI lists interaction. We map pairs.
            # Order in CSV: gene_A, gene_B. 
            # build_graph.py handles direction if Types are correct.
            # But usually we store [Source, Target].
            # For sequestering: lncRNA is source? Or miRNA is source?
            # Biological: lncRNA binds miRNA. 
            # build_graph.py logic:
            # if type_a == 1 (lnc) and type_b == 2 (mir): src=a, dst=b, rel=1
            # So if we output [lnc_id, mir_id], build_graph will pick it up.
            # If we output [mir_id, lnc_id], build_graph handles:
            # elif type_b == 1 and type_a == 2: src=b, dst=a, rel=1
            # So order doesn't matter for build_graph detection!
            
            new_edges.append([my_lnc, my_mir])
            count += 1
            
        print(f"  Added {count} lncRNA-miRNA interactions.")
    else:
        print("  mir2lnc.csv not found.")

    # 4. Save
    if new_edges:
        print(f"Merging {len(new_edges)} total new edges...")
        new_df = pd.DataFrame(new_edges, columns=['gene_A', 'gene_B'])
        
        if os.path.exists(INTERACTIONS_PATH):
            old_df = pd.read_csv(INTERACTIONS_PATH)
            combined = pd.concat([old_df, new_df], ignore_index=True)
        else:
            combined = new_df
            
        final_df = combined.drop_duplicates()
        final_df.to_csv(INTERACTIONS_PATH, index=False)
        print(f"Saved to {INTERACTIONS_PATH}. Total interactions: {len(final_df)}")
    else:
        print("No new interactions found matching our DE lists.")

if __name__ == "__main__":
    main()
