import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import SAGEConv, HeteroConv, Linear

class AdaptiveEdgeDenoiser(nn.Module):
    """
    Adaptive Perturbation Module:
    Learns which edges are noisy/weak and safely drops them to create robust contrastive views.
    """
    def __init__(self, hidden_channels):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(hidden_channels * 2, hidden_channels),
            nn.ReLU(),
            nn.Linear(hidden_channels, 1),
            nn.Sigmoid()
        )
    
    def forward(self, edge_index, src_features, dst_features, threshold=0.5):
        src_nodes = src_features[edge_index[0]]
        dst_nodes = dst_features[edge_index[1]]
        z = torch.cat([src_nodes, dst_nodes], dim=-1)
        edge_scores = self.mlp(z).squeeze()
        mask = edge_scores > threshold
        return edge_index[:, mask]

class AdaptiveMultimodalFusion(nn.Module):
    """
    Attention-Gated Cross-Modal Fusion:
    Instead of simply concatenating an RNA embedding and a Disease/Drug embedding,
    this module uses Multi-Head Attention to align their semantic spaces, and a Sigmoid
    gate to adaptively weight the most important cross-modal features.
    """
    def __init__(self, embed_dim, num_heads=4):
        super().__init__()
        self.attention = nn.MultiheadAttention(embed_dim, num_heads, batch_first=True)
        self.gate = nn.Sequential(
            nn.Linear(embed_dim * 2, embed_dim),
            nn.Sigmoid()
        )
        
    def forward(self, emb1, emb2):
        # Stack into shape: [batch_size, 2, embed_dim]
        stacked = torch.stack([emb1, emb2], dim=1)
        
        # Self-attention across the two modalities
        attn_out, _ = self.attention(stacked, stacked, stacked)
        
        # Gate weights calculated from their raw concatenation
        concat_emb = torch.cat([emb1, emb2], dim=-1)
        gate_weights = self.gate(concat_emb)
        
        # Multiply attention output by adaptive gate weights
        fused_emb = gate_weights * attn_out.mean(dim=1)
        return fused_emb

class ASMSG(torch.nn.Module):
    def __init__(self, hidden_channels=128, out_channels=64, num_layers=2):
        super(ASMSG, self).__init__()
        
        # 1. Feature Projection Layers
        self.lin_dict = torch.nn.ModuleDict({
            'rna': Linear(640, hidden_channels),
            'drug': Linear(384, hidden_channels),
            'disease': Linear(128, hidden_channels)
        })
        
        # 2. Heterogeneous Graph Convolutional Layers
        self.convs = torch.nn.ModuleList()
        for _ in range(num_layers):
            conv = HeteroConv({
                ('rna', 'associates_with', 'disease'): SAGEConv(hidden_channels, hidden_channels),
                ('rna', 'targets', 'drug'): SAGEConv(hidden_channels, hidden_channels),
                ('disease', 'rev_associates_with', 'rna'): SAGEConv(hidden_channels, hidden_channels),
                ('drug', 'rev_targets', 'rna'): SAGEConv(hidden_channels, hidden_channels),
            }, aggr='mean')
            self.convs.append(conv)
            
        # 3. Adaptive Components
        self.edge_denoiser = AdaptiveEdgeDenoiser(hidden_channels)
        self.fusion_module = AdaptiveMultimodalFusion(hidden_channels)
            
        # 4. Final Prediction Head (Link Prediction)
        # It now takes a single hidden_channels vector because the fusion module combined them
        self.predictor = nn.Sequential(
            nn.Linear(hidden_channels, hidden_channels // 2),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(hidden_channels // 2, 1)
        )

    def forward(self, x_dict, edge_index_dict, create_contrastive_view=False):
        h_dict = {}
        for node_type, x in x_dict.items():
            h = F.relu(self.lin_dict[node_type](x))
            if create_contrastive_view:
                h = F.dropout(h, p=0.2, training=self.training)
            h_dict[node_type] = h
            
        working_edge_index_dict = edge_index_dict.copy()
        if create_contrastive_view:
            for edge_type, edge_index in edge_index_dict.items():
                src_type, _, dst_type = edge_type
                working_edge_index_dict[edge_type] = self.edge_denoiser(
                    edge_index, h_dict[src_type], h_dict[dst_type]
                )

        for conv in self.convs:
            h_dict = conv(h_dict, working_edge_index_dict)
            h_dict = {key: F.relu(x) for key, x in h_dict.items()}
            
        return h_dict

    def predict_link(self, rna_embedding, target_embedding):
        """
        Uses Attention-Gated Multimodal Fusion to combine the RNA and Target embeddings,
        then predicts the probability of an association.
        """
        fused_embedding = self.fusion_module(rna_embedding, target_embedding)
        return torch.sigmoid(self.predictor(fused_embedding))

# --- Contrastive Loss Functions ---

def info_nce_loss(z1, z2, temperature=0.5):
    z1 = F.normalize(z1, dim=1)
    z2 = F.normalize(z2, dim=1)
    sim_matrix = torch.matmul(z1, z2.T) / temperature
    labels = torch.arange(z1.size(0)).to(z1.device)
    return F.cross_entropy(sim_matrix, labels)

def margin_triplet_loss(anchor, positive, negative, margin=1.0):
    """
    Forces the anchor node to be closer to its positive augmented view
    than to any randomly sampled negative view.
    """
    return F.triplet_margin_loss(anchor, positive, negative, margin=margin)

if __name__ == "__main__":
    print("Testing ASMSG Architecture instantiation...")
    model = ASMSG(hidden_channels=128)
    print("ASMSG Model with Adaptive Fusion Successfully Built!")
