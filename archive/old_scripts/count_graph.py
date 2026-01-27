# Usage: python scripts/count_graph.py

import torch
from pathlib import Path
p = Path('results/hetero_graph_GBM.pt')
if not p.exists():
    print('ERROR: File not found:', p)
    raise SystemExit(1)
try:
    from torch_geometric.data.data import Data
    torch.serialization.add_safe_globals([Data])
    data = torch.load(p, map_location='cpu', weights_only=False)
except Exception as e:
    print('ERROR loading file:', e)
    raise

print('NUM_NODES', getattr(data, 'x').shape[0])
print('NUM_EDGES', getattr(data, 'edge_index').shape[1])
