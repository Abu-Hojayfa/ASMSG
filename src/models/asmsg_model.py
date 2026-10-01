import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import SAGEConv, HeteroConv, Linear

class ASMSG(torch.nn.Module):
    def __init__(self, hidden_channels=128, out_channels=64, num_layers=2):
        super(ASMSG, self).__init__()
        
        # 1. Feature Projection Layers
        # Our raw features have different dimensions (RNA: 640, Drug: 384, Disease: 128)
        # We need to project them all into a shared hidden dimension (e.g., 128)
        self.lin_dict = torch.nn.ModuleDict({
            'rna': Linear(640, hidden_channels),
            'drug': Linear(384, hidden_channels),
            'disease': Linear(128, hidden_channels)
        })
        
        # 2. Heterogeneous Graph Convolutional Layers
        # We use HeteroConv to apply message passing across the different edge types
        self.convs = torch.nn.ModuleList()
        for _ in range(num_layers):
            conv = HeteroConv({
                ('rna', 'associates_with', 'disease'): SAGEConv(hidden_channels, hidden_channels),
                ('rna', 'targets', 'drug'): SAGEConv(hidden_channels, hidden_channels),
                # We also want messages flowing backwards so the graph is undirected mathematically
                ('disease', 'rev_associates_with', 'rna'): SAGEConv(hidden_channels, hidden_channels),
                ('drug', 'rev_targets', 'rna'): SAGEConv(hidden_channels, hidden_channels),
            }, aggr='mean')
            self.convs.append(conv)
            
        # 3. Final Prediction Head (Link Prediction)
        # Given an RNA embedding and a Disease embedding, predict if an edge exists
        self.predictor = nn.Sequential(
            nn.Linear(hidden_channels * 2, hidden_channels),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(hidden_channels, 1) # Outputs a single score
        )

    def forward(self, x_dict, edge_index_dict):
        # Step 1: Project node features to the same dimensionality
        h_dict = {}
        for node_type, x in x_dict.items():
            h_dict[node_type] = F.relu(self.lin_dict[node_type](x))
            
        # Step 2: Message Passing through the Graph
        for conv in self.convs:
            h_dict = conv(h_dict, edge_index_dict)
            # Apply non-linearity
            h_dict = {key: F.relu(x) for key, x in h_dict.items()}
            
        return h_dict

    def predict_link(self, rna_embedding, target_embedding):
        """
        Combines the RNA and Target (Disease/Drug) embeddings to predict edge probability.
        """
        # Concatenate the two embeddings
        z = torch.cat([rna_embedding, target_embedding], dim=-1)
        # Pass through MLP predictor
        return torch.sigmoid(self.predictor(z))

# If you run this script directly, it will perform a quick architecture test
if __name__ == "__main__":
    print("Testing ASMSG Architecture instantiation...")
    model = ASMSG(hidden_channels=128)
    print("ASMSG Model Successfully Built!")
    print(model)
