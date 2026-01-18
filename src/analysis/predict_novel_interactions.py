import torch
import torch.nn.functional as F
import pandas as pd
import json
import os
import sys

# Allow importing from src
sys.path.append(os.getcwd())
from src.model import RGCN_LinkPredictor

import torch_geometric.transforms as T

# --- Configuration ---
GRAPH_PATH = "results/hetero_graph_GBM.pt"
MODEL_PATH = "results/gnn_model.pth"
MAPPING_PATH = "results/node_mapping.json"
OUTPUT_PATH = "results/novel_predictions.csv"
TOP_K = 100

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# --- Model Definition ---
# Imported from src.model

def main():
    print("--- Starting Novel Interaction Discovery ---")
    
    # 1. Load Data
    data = torch.load(GRAPH_PATH, weights_only=False)
    data = data.to(device)
    print(f"Loaded Graph: {data.num_nodes} nodes, {data.edge_index.size(1)} edges")
    
    with open(MAPPING_PATH, 'r') as f:
        gene_to_idx = json.load(f)
    idx_to_gene = {v: k for k, v in gene_to_idx.items()}
    
    # 2. Load Model
    model = RGCN_LinkPredictor(
        num_nodes=data.num_nodes,
        num_features=data.num_features,
        num_relations=2,
        hidden_channels=32,
        out_channels=16
    ).to(device)
    
    model.load_state_dict(torch.load(MODEL_PATH))
    model.eval()
    print("Loaded Model.")
    
    # 3. Identify Node Types
    # We need to distinguish miRNAs from Genes. 
    # Heuristic: 'mir' in name (or check ID format)
    # Better: Use the mapping.
    
    mirna_indices = []
    gene_indices = []
    
    for idx in range(data.num_nodes):
        name = idx_to_gene[idx]
        if 'mir' in name.lower() or 'let' in name.lower():
            mirna_indices.append(idx)
        else:
            gene_indices.append(idx)
            
    print(f"Nodes: {len(mirna_indices)} miRNAs, {len(gene_indices)} Genes (mRNA/lncRNA)")
    
    # 4. Generate Candidate Edges (miRNA -> Gene)
    # We want to predict NEW regulation.
    # Total possible: len(mirna) * len(gene)
    # Filter out existing edges.
    
    existing_edges = set()
    ei = data.edge_index.cpu().numpy()
    for i in range(ei.shape[1]):
        u, v = ei[0, i], ei[1, i]
        existing_edges.add((u, v))
        
    candidates = []
    for m_idx in mirna_indices:
        for g_idx in gene_indices:
            # Check edge direction 0: miRNA -> mRNA.
            # In our graph, source=miRNA, target=Gene is type 0? 
            # Or is it lncRNA -> miRNA?
            # Let's just predict connectivity in general.
            # But strictly we want miRNA -> Gene.
            if (m_idx, g_idx) not in existing_edges:
                candidates.append([m_idx, g_idx])
                
    print(f"Generated {len(candidates)} candidate pairs (excluding existing edges).")
    
    # 5. Batch Prediction
    # Running all at once might OOM if too large? 
    # 314 * 1700 = ~500k pairs. Should fit in memory easily.
    
    candidate_tensor = torch.tensor(candidates, dtype=torch.long).t().to(device) # [2, NumCandidates]
    
    with torch.no_grad():
        # Encode whole graph first using EXISTING structure
        z = model.encode(data.x, data.edge_index, data.edge_attr)
        
        # Decode candidates
        scores = model.decode(z, candidate_tensor).sigmoid()
        
    # 6. Rank and Save
    scores_np = scores.cpu().numpy()
    
    results = []
    for i in range(len(candidates)):
        score = scores_np[i]
        m_idx = candidates[i][0]
        g_idx = candidates[i][1]
        results.append({
            "miRNA": idx_to_gene[m_idx],
            "Target": idx_to_gene[g_idx],
            "Score": score
        })
        
    df = pd.DataFrame(results)
    df = df.sort_values(by="Score", ascending=False)
    
    top_df = df.head(TOP_K)
    top_df.to_csv(OUTPUT_PATH, index=False)
    print(f"Saved top {TOP_K} novel predictions to {OUTPUT_PATH}")
    
    # 7. Novel miRNA Analysis
    # Which miRNAs appear most frequently in the top predictions?
    top_mirnas = df.head(1000)['miRNA'].value_counts().head(10)
    print("\nTop 'Novel' Regulators (miRNAs with most high-confidence new targets):")
    print(top_mirnas)

if __name__ == "__main__":
    main()
