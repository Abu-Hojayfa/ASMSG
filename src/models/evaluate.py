import os
import torch
import torch_geometric.transforms as T
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix, ConfusionMatrixDisplay, precision_recall_curve, average_precision_score, classification_report
from asmsg_model import ASMSG

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

print("1. Loading Graph and Creating Unseen Test Split...")
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
graph = torch.load('d:/fydp/features/asmsg_hetero_graph.pt', weights_only=False).to(device)
graph = T.ToUndirected()(graph)

# We use a fixed seed here so the test split is deterministic for evaluation
torch.manual_seed(42)
transform = T.RandomLinkSplit(
    num_val=0.1,
    num_test=0.1,
    disjoint_train_ratio=0.3,
    neg_sampling_ratio=1.0,
    add_negative_train_samples=False,
    edge_types=[('rna', 'associates_with', 'disease'), ('rna', 'targets', 'drug')],
    rev_edge_types=None
)
_, _, test_data = transform(graph)

print("2. Loading Trained ASMSG Model...")
model = ASMSG(hidden_channels=128, num_layers=2).to(device)
model.load_state_dict(torch.load('d:/fydp/features/asmsg_best_model.pth', weights_only=True))
model.eval()

# Helper function to evaluate and plot
def evaluate_and_plot(edge_type, title_prefix, save_dir):
    x_dict = test_data.x_dict
    edge_index_dict = test_data.edge_index_dict
    
    with torch.no_grad():
        z_dict = model(x_dict, edge_index_dict, create_contrastive_view=False)
        
    edge_label_index = test_data[edge_type].edge_label_index
    edge_label = test_data[edge_type].edge_label
    
    rna_emb = z_dict[edge_type[0]][edge_label_index[0]]
    target_emb = z_dict[edge_type[2]][edge_label_index[1]]
    
    preds = model.predict_link(rna_emb, target_emb).squeeze().detach().cpu().numpy()
    labels = edge_label.cpu().numpy()
    
    # Metrics
    auc = roc_auc_score(labels, preds)
    aupr = average_precision_score(labels, preds)
    
    binary_preds = (preds > 0.5).astype(int)
    print(f"\n--- {title_prefix} Test Set Results ---")
    print(f"Test AUC:  {auc:.4f}")
    print(f"Test AUPR: {aupr:.4f}")
    print(classification_report(labels, binary_preds, target_names=['No Link (0)', 'True Link (1)']))
    
    # 1. ROC Curve
    fpr, tpr, _ = roc_curve(labels, preds)
    plt.figure(figsize=(6, 5))
    plt.plot(fpr, tpr, color='blue', lw=2, label=f'AUC = {auc:.4f}')
    plt.plot([0, 1], [0, 1], color='gray', linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'{title_prefix} - ROC Curve')
    plt.legend(loc='lower right')
    plt.savefig(os.path.join(save_dir, f'{title_prefix.replace(" ", "_")}_ROC_Curve.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. PR Curve
    precision, recall, _ = precision_recall_curve(labels, preds)
    plt.figure(figsize=(6, 5))
    plt.plot(recall, precision, color='green', lw=2, label=f'AUPR = {aupr:.4f}')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(f'{title_prefix} - Precision-Recall Curve')
    plt.legend(loc='lower left')
    plt.savefig(os.path.join(save_dir, f'{title_prefix.replace(" ", "_")}_PR_Curve.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Confusion Matrix
    cm = confusion_matrix(labels, binary_preds)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=['No Link', 'True Link'])
    fig, ax = plt.subplots(figsize=(6, 5))
    disp.plot(ax=ax, cmap='Blues', colorbar=False)
    plt.title(f'{title_prefix} - Confusion Matrix')
    plt.savefig(os.path.join(save_dir, f'{title_prefix.replace(" ", "_")}_Confusion_Matrix.png'), dpi=300, bbox_inches='tight')
    plt.close()

vis_dir = 'd:/fydp/visualizer'
ensure_dir(vis_dir)

print("\n3. Evaluating on Unseen Test Set...")
evaluate_and_plot(('rna', 'associates_with', 'disease'), 'ncRNA_Disease', vis_dir)
evaluate_and_plot(('rna', 'targets', 'drug'), 'ncRNA_Drug', vis_dir)

print(f"\nEvaluation Complete! All plots saved to {vis_dir}")
