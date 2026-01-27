#!/usr/bin/env python
"""
Test script: Verify graph build with enhanced features and show statistics.
"""
import torch
from pathlib import Path
import pandas as pd
import sys

def main():
    print("=" * 70)
    print("GRAPH BUILD VERIFICATION")
    print("=" * 70)
    
    # Load graph
    try:
        from torch_geometric.data.data import Data
        torch.serialization.add_safe_globals([Data])
        data = torch.load('results/hetero_graph_GBM.pt', map_location='cpu', weights_only=False)
    except Exception as e:
        print(f"ERROR: Could not load graph: {e}")
        return 1
    
    n_nodes = data.x.shape[0]
    n_edges = data.edge_index.shape[1]
    
    print(f"\n📊 GRAPH STATISTICS:")
    print(f"   Total Nodes: {n_nodes}")
    print(f"   Total Edges: {n_edges}")
    
    # Breakdown by relation type
    if data.edge_attr is not None and len(data.edge_attr) > 0:
        edge_attr_list = data.edge_attr.tolist()
        type_counts = {}
        for rel_type in edge_attr_list:
            type_counts[rel_type] = type_counts.get(rel_type, 0) + 1
        
        print(f"\n🔗 EDGES BY RELATION TYPE:")
        print(f"   Type 0 (miRNA → mRNA silencing):  {type_counts.get(0, 0):>6} edges")
        print(f"   Type 1 (lncRNA → miRNA sponging): {type_counts.get(1, 0):>6} edges")
    
    # Load and display edge metadata
    try:
        metadata = pd.read_csv('results/edge_metadata.csv')
        
        print(f"\n📋 EDGE PROVENANCE SUMMARY:")
        prov_counts = metadata['provenance'].value_counts()
        for source, count in prov_counts.items():
            print(f"   {source:.<40} {count:>6} edges")
        
        print(f"\n📈 CONFIDENCE SCORE STATISTICS:")
        conf_stats = metadata['confidence']
        print(f"   Mean:   {conf_stats.mean():.4f}")
        print(f"   Median: {conf_stats.median():.4f}")
        print(f"   Min:    {conf_stats.min():.4f}")
        print(f"   Max:    {conf_stats.max():.4f}")
        print(f"   Std:    {conf_stats.std():.4f}")
        
    except FileNotFoundError:
        print("\n⚠️  Edge metadata file not found (optional)")
    except Exception as e:
        print(f"\n⚠️  Could not load edge metadata: {e}")
    
    print("\n" + "=" * 70)
    print("✅ GRAPH BUILD SUCCESSFUL")
    print("=" * 70)
    print("\nNext steps:")
    print("  1. Download real TargetScan predictions from http://www.targetscan.org/")
    print("  2. Integrate: python src/preprocessing/integrate_targetscan.py --input YOUR_FILE.csv")
    print("  3. Rebuild:   python src/graph/build_graph.py")
    print("  4. Train:     python src/training/train_model.py")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
