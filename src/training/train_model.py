import torch
import torch.nn.functional as F
from torch_geometric.data import HeteroData
from torch_geometric.loader import LinkNeighborLoader
import torch_geometric.transforms as T
from sklearn.metrics import roc_auc_score
import os
import sys

# Allow importing from src
sys.path.append(os.getcwd())
from src.model import RGCN_LinkPredictor

# Device
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

# Paths
GRAPH_PATH = "results/hetero_graph_GBM.pt"

# --- 1. Load Graph ---
if not os.path.exists(GRAPH_PATH):
    raise FileNotFoundError(f"{GRAPH_PATH} not found.")

data = torch.load(GRAPH_PATH, weights_only=False)
print("Loaded Graph:", data)

# --- 2. Preprocessing ---
# Add reverse edges for message passing if not present
# ceRNA edges are directed: miRNA->mRNA, lncRNA->miRNA. 
# RGCN needs to propagate info both ways for representation learning.
# We will use ToUndirected or manually add rev edges. 
# Better: Transform to undirected for message passing, or add explicit reverse relations.
# Let's use T.ToUndirected() but keep edge types separate?
# Actually, HeteroConv handles relations.
# Let's inspect relations first.
# edge_attr is our relation type.
# We need to reshape into HeteroData object if it's currently a homogeneous Data object with edge_types.
# The build_graph.py saved a Data object with edge_index and edge_attr.
# We should convert this to a true HeteroData object for better PyG support.

def convert_to_hetero(data):
    # data.x is [NumNodes, 10]
    # data.edge_index is [2, NumEdges]
    # data.edge_attr is [NumEdges] (relation type 0 or 1)
    
    # We implicitly have different node types but they are all in one matrix X.
    # We can keep it homogeneous or split it.
    # RGCNConv expects homogeneous node features but can handle edge types.
    # If we use fast_rgcn, we pass edge_type tensor.
    
    return data # Keep as homogeneous Data object for RGCNConv

data = convert_to_hetero(data)
data = data.to(device)

# Split edges for Link Prediction
# We need positive edges (existing) and negative edges (random non-existing)
transform = T.RandomLinkSplit(
    num_val=0.1,
    num_test=0.1,
    is_undirected=False, # Directed graph
    add_negative_train_samples=True,
    edge_types=None, # Split all edge types
    rev_edge_types=None, 
)

train_data, val_data, test_data = transform(data)
print(f"Train Edges: {train_data.edge_index.size(1)}")
print(f"Val Edges: {val_data.edge_index.size(1)}")
print(f"Test Edges: {test_data.edge_index.size(1)}")

# --- 3. Model Definition ---
# Imported from src.model import RGCN_LinkPredictor

model = RGCN_LinkPredictor(
    num_nodes=data.num_nodes,
    num_features=data.num_features,
    num_relations=2, # 0: miRNA->mRNA, 1: lncRNA->miRNA
    hidden_channels=32,
    out_channels=16
).to(device)

optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4)
criterion = torch.nn.BCEWithLogitsLoss()

# --- 4. Training Loop ---
def train():
    model.train()
    optimizer.zero_grad()
    
    # Forward Pass (Encode all nodes using message passing on Train Edges)
    # Note: RGCN uses edge_type.
    z = model.encode(train_data.x, train_data.edge_index, train_data.edge_attr)
    
    # Predict on Train Positive Edges
    # train_data.edge_label_index contains the split edges (pos + neg)
    # train_data.edge_label contains 1 or 0
    
    out = model.decode(z, train_data.edge_label_index)
    loss = criterion(out, train_data.edge_label)
    
    loss.backward()
    optimizer.step()
    return loss.item()

@torch.no_grad()
def test(data):
    model.eval()
    z = model.encode(train_data.x, train_data.edge_index, train_data.edge_attr) # Encode using train edges
    out = model.decode(z, data.edge_label_index) # Decode on test pair
    
    return roc_auc_score(data.edge_label.cpu().numpy(), out.sigmoid().cpu().numpy())

print("Starting Training...")
for epoch in range(1, 101):
    loss = train()
    if epoch % 10 == 0:
        val_auc = test(val_data)
        print(f'Epoch: {epoch:03d}, Loss: {loss:.4f}, Val AUC: {val_auc:.4f}')

test_auc = test(test_data)
print(f'Final Test AUC: {test_auc:.4f}')

# Save Model
torch.save(model.state_dict(), "results/gnn_model.pth")
print("Model saved to results/gnn_model.pth")
