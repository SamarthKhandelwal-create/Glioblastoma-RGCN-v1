import torch
import torch.nn as nn
from torch_geometric.nn import RGCNConv, HeteroConv, SAGEConv, Linear


class RGCN_LinkPredictor(nn.Module):
    """RGCN-based link predictor for heterogeneous ceRNA graphs."""
    
    def __init__(self, num_nodes, in_channels, hidden_channels=64, out_channels=32, num_relations=3):
        super().__init__()
        self.num_nodes = num_nodes
        self.num_relations = num_relations
        
        # RGCN encoder layers
        self.conv1 = RGCNConv(in_channels, hidden_channels, num_relations)
        self.conv2 = RGCNConv(hidden_channels, out_channels, num_relations)
        
        # Link prediction decoder
        self.decoder = nn.Sequential(
            nn.Linear(out_channels * 2, hidden_channels),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(hidden_channels, 1)
        )
    
    def encode(self, x, edge_index, edge_type):
        """Encode node features using RGCN layers."""
        x = self.conv1(x, edge_index, edge_type).relu()
        x = torch.nn.functional.dropout(x, p=0.3, training=self.training)
        x = self.conv2(x, edge_index, edge_type)
        return x
    
    def decode(self, z, edge_index):
        """Decode edge probabilities from node embeddings."""
        src, dst = edge_index
        edge_features = torch.cat([z[src], z[dst]], dim=1)
        return self.decoder(edge_features).squeeze()
    
    def forward(self, x, edge_index, edge_type, pred_edges):
        """Full forward pass: encode + decode."""
        z = self.encode(x, edge_index, edge_type)
        return self.decode(z, pred_edges)


class HeteroRGCN_LinkPredictor(torch.nn.Module):
    """Heterogeneous link predictor using type projections and HeteroConv of SAGEConv per relation.

    - Projects each node type to a common hidden dimension.
    - Uses multiple HeteroConv layers composed of SAGEConv for each relation.
    - Provides helpers to compute dot-product scores for relation-specific edge index tensors.
    """

    def __init__(self, metadata, in_channels_dict, hidden_channels=64, out_channels=64, num_layers=2):
        super().__init__()
        # metadata: (node_types, edge_types)
        self.metadata = metadata
        self.node_types = metadata[0]
        self.edge_types = metadata[1]

        # Per-node-type input projection to common hidden dim
        self.input_proj = torch.nn.ModuleDict()
        for ntype in self.node_types:
            in_ch = in_channels_dict.get(ntype, None)
            if in_ch is None:
                # Assume scalar feature fallback
                in_ch = 1
            self.input_proj[ntype] = Linear(in_ch, hidden_channels)

        # Build hetero conv layers
        self.convs = torch.nn.ModuleList()
        for _ in range(num_layers):
            conv_dict = {}
            for src, rel, dst in self.edge_types:
                # Use SAGEConv per relation (src->dst)
                conv_dict[(src, rel, dst)] = SAGEConv((-1, -1), hidden_channels)
            self.convs.append(HeteroConv(conv_dict, aggr='mean'))

        self.final_proj = torch.nn.ModuleDict({nt: Linear(hidden_channels, out_channels) for nt in self.node_types})

    def forward(self, x_dict, edge_index_dict):
        # x_dict: {ntype: tensor}
        # edge_index_dict: {(src, rel, dst): edge_index}
        # Project inputs
        h = {}
        for ntype, x in x_dict.items():
            h[ntype] = self.input_proj[ntype](x).relu()

        for conv in self.convs:
            h = conv(h, edge_index_dict)
            # Apply activation and dropout per node type
            for ntype in h:
                h[ntype] = F.relu(h[ntype])
                h[ntype] = F.dropout(h[ntype], p=0.5, training=self.training)

        # Final projection
        out = {nt: self.final_proj[nt](h[nt]) for nt in h}
        return out

    def decode_relation(self, node_emb_dict, edge_index):
        # edge_index: [2, num_edges] with source and target indices for a single relation
        src_idx = edge_index[0]
        dst_idx = edge_index[1]
        # Assume node_emb_dict corresponds to the correct node type for this edge index
        # Caller must select appropriate embeddings per relation
        src_emb = node_emb_dict[0][src_idx]
        dst_emb = node_emb_dict[1][dst_idx]
        return (src_emb * dst_emb).sum(dim=1)

    # Utility decoder for relation-aware dict: returns dict of scores per relation
    def decode_all(self, node_emb_dict, edge_index_dict):
        scores = {}
        for key, eidx in edge_index_dict.items():
            src_type, _, dst_type = key
            src_emb = node_emb_dict[src_type]
            dst_emb = node_emb_dict[dst_type]
            scores[key] = (src_emb[eidx[0]] * dst_emb[eidx[1]]).sum(dim=1)
        return scores

    # TODO: Add margin-ranking or classification heads per relation if desired

