"""
Create a minimal test graph for demonstration.
This creates a simple homogeneous graph with random interactions
for testing the CV pipeline without requiring full data integration.
"""

import torch
import pandas as pd
from pathlib import Path
from torch_geometric.data import Data

def create_test_graph():
    """Create a minimal test graph from existing node features."""
    results_dir = Path("results")
    
    # Load node features
    features_path = results_dir / "node_features_matrix.csv"
    features_df = pd.read_csv(features_path, index_col=0)
    
    num_nodes = len(features_df)
    x = torch.tensor(features_df.values, dtype=torch.float)
    
    print(f"✓ Loaded node features: {x.shape}")
    
    # Create random edges for testing (10% of possible edges)
    num_edges = int(num_nodes * 0.1)
    edge_indices = []
    for _ in range(num_edges):
        src = torch.randint(0, num_nodes, (1,)).item()
        dst = torch.randint(0, num_nodes, (1,)).item()
        if src != dst:
            edge_indices.append([src, dst])
    
    if edge_indices:
        edge_index = torch.tensor(edge_indices, dtype=torch.long).t().contiguous()
    else:
        edge_index = torch.empty((2, 0), dtype=torch.long)
    
    edge_attr = torch.zeros(edge_index.shape[1], dtype=torch.long)
    
    data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)
    
    output_path = results_dir / "hetero_graph_GBM.pt"
    torch.save(data, output_path)
    
    print(f"✓ Test graph created: {output_path}")
    print(f"  Nodes: {data.num_nodes}")
    print(f"  Edges: {data.num_edges}")
    print(f"  Features: {data.x.shape[1]}")
    
    return output_path

if __name__ == "__main__":
    create_test_graph()
