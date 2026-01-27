import torch
from pathlib import Path
p = Path('results/hetero_graph_GBM.pt')
if not p.exists():
    print('ERROR: File not found:', p)
    raise SystemExit(1)
try:
    # Allowlist PyG Data class for safe deserialization (PyTorch 2.6+)
    from torch_geometric.data.data import Data
    torch.serialization.add_safe_globals([Data])
    data = torch.load(p, map_location='cpu', weights_only=False)
except Exception as e:
    print('ERROR loading file:', e)
    raise

print('Loaded object type:', type(data))
# node count
try:
    xshape = getattr(data, 'x').shape
    n_nodes = xshape[0]
except Exception as e:
    n_nodes = None
    print('Could not read x:', e)
# edge count
try:
    edge_shape = getattr(data, 'edge_index').shape
    n_edges = edge_shape[1]
except Exception as e:
    n_edges = None
    print('Could not read edge_index:', e)

print('NUM_NODES', n_nodes)
print('NUM_EDGES', n_edges)
print('edge_index shape:', getattr(data, 'edge_index', None))
print('edge_attr shape:', getattr(data, 'edge_attr', None))
# temporary placeholder
print('noop')
