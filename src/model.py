import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import RGCNConv, HeteroConv, SAGEConv, Linear


class RGCN_LinkPredictor(nn.Module):
    
    def __init__(self, num_nodes, in_channels=None, hidden_channels=64, out_channels=32, 
                 num_relations=3, num_features=None):
        super().__init__()
        self.num_nodes = num_nodes
        self.num_relations = num_relations
        
        feat_dim = in_channels if in_channels is not None else num_features
        
        self.conv1 = RGCNConv(feat_dim, hidden_channels, num_relations)
        self.conv2 = RGCNConv(hidden_channels, out_channels, num_relations)
        
        self.decoder = nn.Sequential(
            nn.Linear(out_channels * 2, hidden_channels),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(hidden_channels, 1)
        )
    
    def encode(self, x, edge_index, edge_type):
        x = self.conv1(x, edge_index, edge_type).relu()
        x = F.dropout(x, p=0.3, training=self.training)
        x = self.conv2(x, edge_index, edge_type)
        return x
    
    def decode(self, z, edge_index):
        src, dst = edge_index
        edge_features = torch.cat([z[src], z[dst]], dim=1)
        return self.decoder(edge_features).squeeze()
    
    def forward(self, x, edge_index, edge_type, pred_edges):
        z = self.encode(x, edge_index, edge_type)
        return self.decode(z, pred_edges)


class HeteroRGCN_LinkPredictor(nn.Module):

    def __init__(self, metadata, in_channels_dict, hidden_channels=64, out_channels=64, num_layers=2):
        super().__init__()
        self.metadata = metadata
        self.node_types = metadata[0]
        self.edge_types = metadata[1]

        self.input_proj = nn.ModuleDict()
        for ntype in self.node_types:
            in_ch = in_channels_dict.get(ntype, 1)
            self.input_proj[ntype] = Linear(in_ch, hidden_channels)

        self.convs = nn.ModuleList()
        for _ in range(num_layers):
            conv_dict = {}
            for src, rel, dst in self.edge_types:
                conv_dict[(src, rel, dst)] = SAGEConv((-1, -1), hidden_channels)
            self.convs.append(HeteroConv(conv_dict, aggr='mean'))

        self.final_proj = nn.ModuleDict({nt: Linear(hidden_channels, out_channels) for nt in self.node_types})

    def forward(self, x_dict, edge_index_dict):
        h = {}
        for ntype, x in x_dict.items():
            h[ntype] = self.input_proj[ntype](x).relu()

        for conv in self.convs:
            h = conv(h, edge_index_dict)
            for ntype in h:
                h[ntype] = F.relu(h[ntype])
                h[ntype] = F.dropout(h[ntype], p=0.5, training=self.training)

        out = {nt: self.final_proj[nt](h[nt]) for nt in h}
        return out

    def decode_relation(self, node_emb_dict, edge_index):
        src_idx = edge_index[0]
        dst_idx = edge_index[1]
        src_emb = node_emb_dict[0][src_idx]
        dst_emb = node_emb_dict[1][dst_idx]
        return (src_emb * dst_emb).sum(dim=1)

    def decode_all(self, node_emb_dict, edge_index_dict):
        scores = {}
        for key, eidx in edge_index_dict.items():
            src_type, _, dst_type = key
            src_emb = node_emb_dict[src_type]
            dst_emb = node_emb_dict[dst_type]
            scores[key] = (src_emb[eidx[0]] * dst_emb[eidx[1]]).sum(dim=1)
        return scores

