import torch
import pandas as pd
import numpy as np
import os
import json
from torch_geometric.data import Data

# Paths
RESULTS_DIR = r"c:\Users\Samarth\OneDrive\Documents\AP Research\results"
FEATURE_MATRIX_PATH = os.path.join(RESULTS_DIR, "node_features_matrix.csv")
INTERACTIONS_PATH = os.path.join(RESULTS_DIR, "interactions.csv")
OUTPUT_GRAPH_PATH = os.path.join(RESULTS_DIR, "hetero_graph_GBM.pt")
OUTPUT_MAPPING_PATH = os.path.join(RESULTS_DIR, "node_mapping.json")

# DE Files for Node Type Assignment
DE_MRNA_PATH = os.path.join(RESULTS_DIR, "DE_mRNA_logFC_abs_gt1.tsv")
DE_LNCRNA_PATH = os.path.join(RESULTS_DIR, "DE_lncRNA_logFC_abs_gt1.tsv")
DE_MIRNA_PATH = os.path.join(RESULTS_DIR, "DE_miRNA_logFC_abs_gt1.tsv")

def normalize_id(gene_id):
    """Strips Ensembl version numbers and lowercases miRNAs."""
    s = str(gene_id).strip()
    # Ensembl Version Strip
    if s.startswith('ENSG'):
        return s.split('.')[0]
    # miRNA Normalization
    return s.lower()

def main():
    print("--- Starting Phase 2: Heterogeneous Graph Construction (Robust ID Matching) ---")
    
    # 1. Load Feature Matrix
    print(f"Loading feature matrix from {FEATURE_MATRIX_PATH}...")
    if not os.path.exists(FEATURE_MATRIX_PATH):
        raise FileNotFoundError(f"{FEATURE_MATRIX_PATH} not found. Run Phase 1 first.")
        
    node_features_df = pd.read_csv(FEATURE_MATRIX_PATH, index_col=0)
    print(f"Loaded Feature Matrix: Shape {node_features_df.shape}")
    
    # 2. Load Gene Lists for Type Assignment & ID Mapping
    try:
        de_mrna = pd.read_csv(DE_MRNA_PATH, sep='\t') # ID col: ensembl_id
        de_lncrna = pd.read_csv(DE_LNCRNA_PATH, sep='\t') # ID col: ensembl_id
        de_mirna = pd.read_csv(DE_MIRNA_PATH, sep='\t') # ID col: mirna_id
    except FileNotFoundError as e:
        print(f"Critical Error: DE files not found. {e}")
        return

    # Create Sets and Version Maps (Base ID -> Versioned ID in Matrix)
    # The Matrix uses Versioned IDs (e.g. ENSG...14).
    # We need to map interactions (which might be Unversioned ENSG... or ENSG...12) to these Matrix IDs.
    
    # Helper to map CleanID -> MatrixID
    matrix_ids = set(node_features_df.index)
    id_map = {} # Clean/Normalized -> MatrixID
    
    for mid in matrix_ids:
        clean = normalize_id(mid)
        id_map[clean] = mid
        # Also map exact just in case
        id_map[mid] = mid
        
    print(f"Built ID Map with {len(id_map)} entries from {len(matrix_ids)} matrix nodes.")

    # Sets for typing (using Matrix IDs if possible, or Clean IDs for check)
    # We will trace back to Matrix IDs for the final node list.
    
    # Determine Types based on DE membership
    mrna_clean_set = set(de_mrna['ensembl_id'].apply(normalize_id))
    lncrna_clean_set = set(de_lncrna['ensembl_id'].apply(normalize_id))
    mirna_clean_set = set(de_mirna['mirna_id'].apply(normalize_id))
    
    # 2b. Inject Zero-Features for missing miRNAs (REMOVED per user request)
    # The user requested to remove node injections. 
    # NOTE: This means miRNAs without features in the CSV will NOT be in the graph.
    # Since all supported edges involve miRNAs, this may result in 0 edges unless 
    # miRNAs are already present in the feature matrix.
    
    # skipped injection logic
    
    # Re-build ID Map (just in case)
    matrix_ids = set(node_features_df.index)
    id_map = {} 
    for mid in matrix_ids:
        clean = normalize_id(mid)
        id_map[clean] = mid
        id_map[mid] = mid # exact match fallback
        
    # Re-build ID Map with injected nodes
    matrix_ids = set(node_features_df.index)
    id_map = {} 
    for mid in matrix_ids:
        clean = normalize_id(mid)
        id_map[clean] = mid
        id_map[mid] = mid # exact match fallback

    
    # 3. Align Matrix
    master_node_list = node_features_df.index.tolist()
    gene_to_index = {gene: idx for idx, gene in enumerate(master_node_list)}
    
    # Save Mapping
    with open(OUTPUT_MAPPING_PATH, 'w') as f:
        json.dump(gene_to_index, f, indent=4)
    
    # Tensor
    x = torch.tensor(node_features_df.values, dtype=torch.float)
    print(f"Feature Tensor x shape: {x.shape}")
    
    # 4. Assign Node Types
    # Map: 0=mRNA, 1=lncRNA, 2=miRNA
    gene_type_map = {}
    type_counts = {0: 0, 1: 0, 2: 0, -1: 0}
    
    for gene in master_node_list:
        clean = normalize_id(gene)
        g_type = -1
        
        if clean in mrna_clean_set:
            g_type = 0
        elif clean in lncrna_clean_set:
            g_type = 1
        elif clean in mirna_clean_set:
            g_type = 2
        else:
            # Fallback by pattern if not in DE list (shouldn't happen if filtered correctly)
            if 'ENSG' in gene:
                g_type = 0 # Assume mRNA/lncRNA default? Dangerous.
                # Check lncRNA again?
                pass
            elif 'mir' in gene.lower() or 'let' in gene.lower():
                g_type = 2
        
        gene_type_map[gene] = g_type
        type_counts[g_type] += 1
        
    print(f"Node Type Assignment: mRNA={type_counts[0]}, lncRNA={type_counts[1]}, miRNA={type_counts[2]}, Unknown={type_counts[-1]}")

    # 5. Process Interactions with Fuzzy Matching
    print(f"Loading interactions from {INTERACTIONS_PATH}...")
    interactions_df = pd.read_csv(INTERACTIONS_PATH)
    
    edge_list = []
    edge_types = []
    unique_edges = set()
    
    processed = 0
    accepted = 0
    
    for _, row in interactions_df.iterrows():
        raw_a = str(row['gene_A']).strip()
        raw_b = str(row['gene_B']).strip()
        processed += 1
        
        # Try to resolve to Matrix ID
        gene_a = id_map.get(normalize_id(raw_a))
        gene_b = id_map.get(normalize_id(raw_b))
        
        if not gene_a or not gene_b:
            # Debug sample failure
            if processed < 5:
                print(f"Skipping edge {raw_a}-{raw_b}: resolved ({gene_a}, {gene_b})")
            continue
            
        if gene_a not in gene_to_index or gene_b not in gene_to_index:
            continue
            
        idx_a = gene_to_index[gene_a]
        idx_b = gene_to_index[gene_b]
        
        type_a = gene_type_map[gene_a]
        type_b = gene_type_map[gene_b]
        
        src, dst, rel = None, None, None
        
        # 0: miRNA -> mRNA (Silencing)
        if type_a == 2 and type_b == 0:
            src, dst, rel = idx_a, idx_b, 0
        elif type_b == 2 and type_a == 0:
            src, dst, rel = idx_b, idx_a, 0
            
        # 1: lncRNA -> miRNA (Sequestering)
        elif type_a == 1 and type_b == 2:
            src, dst, rel = idx_a, idx_b, 1
        elif type_b == 1 and type_a == 2:
            src, dst, rel = idx_b, idx_a, 1
            
        if src is not None:
            if (src, dst, rel) not in unique_edges:
                unique_edges.add((src, dst, rel))
                edge_list.append([src, dst])
                edge_types.append(rel)
                accepted += 1
                
    print(f"Processed {processed} rows. Accepted {accepted} valid edges.")
    
    if accepted == 0:
        print("WARNING: No edges generated. Check ID mapping.")
        # Proceeding to save anyway to avoid crash, but object will be empty
    
    # 6. Save Graph
    if edge_list:
        edge_index = torch.tensor(edge_list, dtype=torch.long).t().contiguous()
        edge_attr = torch.tensor(edge_types, dtype=torch.long)
    else:
        edge_index = torch.empty((2, 0), dtype=torch.long)
        edge_attr = torch.empty((0,), dtype=torch.long)
        
    data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)
    torch.save(data, OUTPUT_GRAPH_PATH)
    print(f"Saved Heterogeneous Graph to {OUTPUT_GRAPH_PATH}")
    print(data)

if __name__ == "__main__":
    main()
