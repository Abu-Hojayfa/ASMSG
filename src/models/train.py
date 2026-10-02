import torch
import torch.nn.functional as F
import torch_geometric.transforms as T
import copy
from sklearn.metrics import roc_auc_score
from torch_geometric.utils import negative_sampling
from asmsg_model import ASMSG, info_nce_loss, margin_triplet_loss

# 1. Setup
print("1. Loading the ASMSG Heterogeneous Graph...")
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
graph = torch.load('d:/fydp/features/asmsg_hetero_graph.pt', weights_only=False).to(device)
graph = T.ToUndirected()(graph)

print(f"2. Performing Inductive Data Splitting on {device}...")
transform = T.RandomLinkSplit(
    num_val=0.1,
    num_test=0.1,
    disjoint_train_ratio=0.3, 
    neg_sampling_ratio=1.0,
    add_negative_train_samples=False,
    edge_types=[('rna', 'associates_with', 'disease'), ('rna', 'targets', 'drug')],
    rev_edge_types=None
)
train_data, val_data, test_data = transform(graph)

# 3. Model & Optimizer setup
print("3. Building ASMSG Model with Adaptive Perturbation and Multimodal Fusion...")
model = ASMSG(hidden_channels=128, num_layers=2).to(device)
# Reduced learning rate for better convergence (Hyperparameter Tweak)
optimizer = torch.optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-5)

lambda_info_nce = 0.1 
lambda_triplet = 0.05

print("\n--- Starting Training Loop (Max Epochs: 200) ---")

def train():
    model.train()
    optimizer.zero_grad()
    
    x_dict = train_data.x_dict
    edge_index_dict = train_data.edge_index_dict

    # Standard vs Augmented Views
    z1_dict = model(x_dict, edge_index_dict, create_contrastive_view=False)
    z2_dict = model(x_dict, edge_index_dict, create_contrastive_view=True)
    
    contrastive_loss = 0
    triplet_loss = 0
    for node_type in ['rna', 'disease', 'drug']:
        anchor = z1_dict[node_type]
        positive = z2_dict[node_type]
        # Shuffle embeddings to create a random "negative" contrastive view
        indices = torch.randperm(anchor.size(0))
        negative = z1_dict[node_type][indices]
        
        contrastive_loss += info_nce_loss(anchor, positive)
        triplet_loss += margin_triplet_loss(anchor, positive, negative)
    
    link_loss = 0
    for edge_type in [('rna', 'associates_with', 'disease'), ('rna', 'targets', 'drug')]:
        edge_label_index = train_data[edge_type].edge_label_index
        num_src, num_dst = x_dict[edge_type[0]].size(0), x_dict[edge_type[2]].size(0)
        
        neg_edge_index = negative_sampling(
            edge_index=edge_index_dict[edge_type], 
            num_nodes=(num_src, num_dst),
            num_neg_samples=edge_label_index.size(1), 
            method='sparse'
        )
        
        train_edge_index = torch.cat([edge_label_index, neg_edge_index], dim=-1)
        train_labels = torch.cat([torch.ones(edge_label_index.size(1)), torch.zeros(neg_edge_index.size(1))]).to(device)
        
        rna_emb = z1_dict[edge_type[0]][train_edge_index[0]]
        target_emb = z1_dict[edge_type[2]][train_edge_index[1]]
        predictions = model.predict_link(rna_emb, target_emb).squeeze()
        
        link_loss += F.binary_cross_entropy(predictions, train_labels)

    total_loss = link_loss + (lambda_info_nce * contrastive_loss) + (lambda_triplet * triplet_loss)
    total_loss.backward()
    optimizer.step()
    
    return total_loss.item(), link_loss.item(), contrastive_loss.item(), triplet_loss.item()

@torch.no_grad()
def test(data):
    model.eval()
    x_dict = data.x_dict
    edge_index_dict = data.edge_index_dict
    z_dict = model(x_dict, edge_index_dict, create_contrastive_view=False)
    
    metrics = {}
    for edge_type in [('rna', 'associates_with', 'disease'), ('rna', 'targets', 'drug')]:
        edge_label_index = data[edge_type].edge_label_index
        edge_label = data[edge_type].edge_label
        
        rna_emb = z_dict[edge_type[0]][edge_label_index[0]]
        target_emb = z_dict[edge_type[2]][edge_label_index[1]]
        
        preds = model.predict_link(rna_emb, target_emb).squeeze().cpu().numpy()
        labels = edge_label.cpu().numpy()
        metrics[edge_type[1]] = roc_auc_score(labels, preds)
        
    return metrics

# Early Stopping Variables
best_val_auc = 0.0
patience = 25
patience_counter = 0
best_model_state = None

for epoch in range(1, 201):
    loss, link_l, cont_l, trip_l = train()
    
    if epoch % 5 == 0:
        val_metrics = test(val_data)
        avg_val_auc = (val_metrics['associates_with'] + val_metrics['targets']) / 2.0
        
        print(f"Epoch {epoch:03d} | Loss: {loss:.4f} (Link: {link_l:.4f}) | "
              f"Val Disease AUC: {val_metrics['associates_with']:.4f} | Val Drug AUC: {val_metrics['targets']:.4f}")
        
        # Early Stopping Logic
        min_delta = 0.001  # Requires at least a 0.1% improvement
        
        if avg_val_auc > (best_val_auc + min_delta):
            best_val_auc = avg_val_auc
            patience_counter = 0
            best_model_state = copy.deepcopy(model.state_dict())
        else:
            patience_counter += 5
            
        if patience_counter >= patience:
            print(f"\n[Early Stopping] Triggered at epoch {epoch}. No improvement for {patience} epochs.")
            break

print(f"\nTraining Complete! Best Average Validation AUC: {best_val_auc:.4f}")
# Load the best model weights
model.load_state_dict(best_model_state)
torch.save(model.state_dict(), 'd:/fydp/features/asmsg_best_model.pth')
print("Best model saved to 'features/asmsg_best_model.pth'.")
