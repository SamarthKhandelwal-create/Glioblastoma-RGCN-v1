#!/usr/bin/env python
"""Quick script to rebuild graph and count edges."""
import subprocess
import sys
import os

os.environ['INJECT_MIRNA'] = 'true'

print("=" * 60)
print("REBUILDING GRAPH WITH ceRNA INFERENCE")
print("=" * 60)

# Run graph builder
result = subprocess.run([sys.executable, 'src/graph/build_graph.py'], cwd='.')
if result.returncode != 0:
    print(f"\nERROR: build_graph.py failed with code {result.returncode}")
    sys.exit(1)

print("\n" + "=" * 60)
print("COUNTING NODES AND EDGES")
print("=" * 60)

# Count nodes/edges
count_code = """
import torch
from pathlib import Path

p = Path('results/hetero_graph_GBM.pt')
if not p.exists():
    print('ERROR: Graph file not found')
    raise SystemExit(1)

from torch_geometric.data.data import Data
torch.serialization.add_safe_globals([Data])
data = torch.load(p, map_location='cpu', weights_only=False)

n_nodes = data.x.shape[0]
n_edges = data.edge_index.shape[1]

print(f'NUM_NODES: {n_nodes}')
print(f'NUM_EDGES: {n_edges}')

# Count by relation type
if data.edge_attr is not None and len(data.edge_attr) > 0:
    unique_types = {}
    for rel_type in data.edge_attr.tolist():
        unique_types[rel_type] = unique_types.get(rel_type, 0) + 1
    print(f'\\nEdges by relation type:')
    for rel_type in sorted(unique_types.keys()):
        print(f'  Type {rel_type}: {unique_types[rel_type]} edges')
"""

result = subprocess.run([sys.executable, '-c', count_code], cwd='.')
sys.exit(result.returncode)
