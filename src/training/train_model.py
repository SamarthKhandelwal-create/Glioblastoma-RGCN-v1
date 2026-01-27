import torch
import torch.nn.functional as F
import torch_geometric.transforms as T
from sklearn.metrics import roc_auc_score
import os
import sys
from pathlib import Path

sys.path.append(os.getcwd())
from src.model import RGCN_LinkPredictor

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

GRAPH_PATH = "results/hetero_graph_GBM.pt"

if not os.path.exists(GRAPH_PATH):
    raise FileNotFoundError(f"{GRAPH_PATH} not found.")

data = torch.load(GRAPH_PATH, weights_only=False)
print("Loaded Graph:", data)

data = data.to(device)

transform = T.RandomLinkSplit(
    num_val=0.1,
    num_test=0.1,
    is_undirected=False,
    add_negative_train_samples=True,
)

train_data, val_data, test_data = transform(data)
print(f"Train Edges: {train_data.edge_index.size(1)}")
print(f"Val Edges: {val_data.edge_index.size(1)}")
print(f"Test Edges: {test_data.edge_index.size(1)}")

num_relations = int(data.edge_attr.max().item()) + 1 if data.edge_attr is not None else 3
print(f"Number of relation types: {num_relations}")

model = RGCN_LinkPredictor(
    num_nodes=data.x.shape[0],
    in_channels=data.x.shape[1], 
    hidden_channels=64,
    out_channels=32,
    num_relations=num_relations
).to(device)

optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4)
criterion = torch.nn.BCEWithLogitsLoss()


def train():
    model.train()
    optimizer.zero_grad()
    
    z = model.encode(train_data.x, train_data.edge_index, train_data.edge_attr)
    out = model.decode(z, train_data.edge_label_index)
    loss = criterion(out, train_data.edge_label)
    
    loss.backward()
    optimizer.step()
    return loss.item()


@torch.no_grad()
def test(data):
    model.eval()
    z = model.encode(train_data.x, train_data.edge_index, train_data.edge_attr)
    out = model.decode(z, data.edge_label_index)
    return roc_auc_score(data.edge_label.cpu().numpy(), out.sigmoid().cpu().numpy())


print("Starting Training...")
for epoch in range(1, 101):
    loss = train()
    if epoch % 10 == 0:
        val_auc = test(val_data)
        print(f'Epoch: {epoch:03d}, Loss: {loss:.4f}, Val AUC: {val_auc:.4f}')

test_auc = test(test_data)
print(f'Final Test AUC: {test_auc:.4f}')

torch.save(model.state_dict(), "results/model.pt")
print("Model saved to results/model.pt")
