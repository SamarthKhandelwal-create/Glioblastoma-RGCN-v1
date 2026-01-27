#!/usr/bin/env python
"""Verify graph build and display statistics."""
import torch
from pathlib import Path
import pandas as pd
import sys


def main():
    print("=" * 60)
    print("GRAPH VERIFICATION")
    print("=" * 60)
    
    try:
        from torch_geometric.data.data import Data
        torch.serialization.add_safe_globals([Data])
        data = torch.load('results/hetero_graph_GBM.pt', map_location='cpu', weights_only=False)
    except Exception as e:
        print(f"ERROR: Could not load graph: {e}")
        return 1
    
    n_nodes = data.x.shape[0]
    n_edges = data.edge_index.shape[1]
    
    print(f"\nNodes: {n_nodes}")
    print(f"Edges: {n_edges}")
    
    if data.edge_attr is not None and len(data.edge_attr) > 0:
        edge_attr_list = data.edge_attr.tolist()
        type_counts = {}
        for rel_type in edge_attr_list:
            type_counts[rel_type] = type_counts.get(rel_type, 0) + 1
        
        print(f"\nEdges by type:")
        print(f"  Type 0 (miRNA->mRNA):   {type_counts.get(0, 0)}")
        print(f"  Type 1 (lncRNA->miRNA): {type_counts.get(1, 0)}")
    
    try:
        metadata = pd.read_csv('results/edge_metadata.csv')
        print(f"\nEdge provenance:")
        for source, count in metadata['provenance'].value_counts().items():
            print(f"  {source}: {count}")
    except FileNotFoundError:
        pass
    
    print("\n" + "=" * 60)
    return 0


if __name__ == '__main__':
    sys.exit(main())
