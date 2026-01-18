import torch
import torch.nn.functional as F
from torch_geometric.nn import RGCNConv

class RGCN_LinkPredictor(torch.nn.Module):
    def __init__(self, num_nodes, num_features, num_relations, hidden_channels, out_channels):
        super().__init__()
        # Encoder: RGCN
        self.conv1 = RGCNConv(num_features, hidden_channels, num_relations)
        self.conv2 = RGCNConv(hidden_channels, out_channels, num_relations)
        
    def encode(self, x, edge_index, edge_type):
        x = self.conv1(x, edge_index, edge_type).relu()
        # Optional: Add dropout here if needed consistently
        # x = F.dropout(x, p=0.5, training=self.training) 
        # Note: Previous train script had dropout only in training loop? 
        # No, it had it in encode() with training=self.training.
        # Let's add it back relative to the training flag.
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.conv2(x, edge_index, edge_type)
        return x

    def decode(self, z, edge_index):
        # Dot product
        src_idx = edge_index[0]
        dst_idx = edge_index[1]
        
        src_z = z[src_idx]
        dst_z = z[dst_idx]
        
        return (src_z * dst_z).sum(dim=1)
