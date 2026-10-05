import os
import json
import torch
import pandas as pd
from torch_geometric.data import HeteroData

def build_hetero_graph():
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    clean_dir = os.path.join(base_dir, 'datasets', 'asmsg_clean')
    drugs_csv = os.path.join(base_dir, 'datasets', 'drugbank', 'drug_smiles.csv')
    out_graph_pt = os.path.join(clean_dir, 'asmsg_hetero_graph.pt')
    
    print("\n--- Constructing PyTorch Geometric HeteroData Graph ---")
    
    nodes_df = pd.read_csv(os.path.join(clean_dir, 'asmsg_nodes.csv'))
    diseases_df = pd.read_csv(os.path.join(clean_dir, 'asmsg_diseases.csv'))
    edges_df = pd.read_csv(os.path.join(clean_dir, 'asmsg_edges.csv'))
    drugs_df = pd.read_csv(drugs_csv) if os.path.exists(drugs_csv) else None
    
    ncrna_to_idx = {sym: i for i, sym in enumerate(nodes_df['RNA Symbol'])}
    disease_to_idx = {dis: i for i, dis in enumerate(diseases_df['Disease_ID'])}
    
    graph = HeteroData()
    
    # 1. Set Node Counts & Map Features
    graph['ncRNA'].num_nodes = len(nodes_df)
    graph['Disease'].num_nodes = len(diseases_df)
    
    # Load feature tensors if available
    rna_fm_path = os.path.join(clean_dir, 'rna_fm_embeddings.pt')
    kmer_path = os.path.join(clean_dir, 'ablation_3mer_features.pt')
    biobert_path = os.path.join(clean_dir, 'disease_biobert_embeddings.pt')
    chemberta_path = os.path.join(clean_dir, 'drug_chemberta_embeddings.pt')
    
    if os.path.exists(rna_fm_path):
        graph['ncRNA'].x = torch.load(rna_fm_path)
        print(f" - Attached RNA-FM embeddings {graph['ncRNA'].x.shape} to ncRNA nodes")
    elif os.path.exists(kmer_path):
        graph['ncRNA'].x = torch.load(kmer_path)
        print(f" - Attached 3-mer ablation features {graph['ncRNA'].x.shape} to ncRNA nodes")
    else:
        # Fallback to random initialization for graph structure compilation
        graph['ncRNA'].x = torch.randn(len(nodes_df), 640)
        print(f" - Initialized fallback random features {graph['ncRNA'].x.shape} for ncRNA nodes")
        
    if os.path.exists(biobert_path):
        graph['Disease'].x = torch.load(biobert_path)
        print(f" - Attached BioBERT embeddings {graph['Disease'].x.shape} to Disease nodes")
    else:
        graph['Disease'].x = torch.randn(len(diseases_df), 768)
        print(f" - Initialized fallback random features {graph['Disease'].x.shape} for Disease nodes")
        
    if drugs_df is not None:
        graph['Drug'].num_nodes = len(drugs_df)
        if os.path.exists(chemberta_path):
            graph['Drug'].x = torch.load(chemberta_path)
            print(f" - Attached ChemBERTa-2 embeddings {graph['Drug'].x.shape} to Drug nodes")
        else:
            graph['Drug'].x = torch.randn(len(drugs_df), 384)
            print(f" - Initialized fallback random features {graph['Drug'].x.shape} for Drug nodes")
            
    # 2. Build Edge Indices and Continuous Weights
    src_indices = [ncrna_to_idx[r] for r in edges_df['RNA Symbol']]
    dst_indices = [disease_to_idx[d] for d in edges_df['Disease_ID']]
    
    edge_index = torch.tensor([src_indices, dst_indices], dtype=torch.long)
    edge_weight = torch.tensor(edges_df['Edge_Weight'].values, dtype=torch.float32)
    conf_score = torch.tensor(edges_df['Confidence_Score'].values, dtype=torch.float32)
    
    # Forward edges: (ncRNA -> associated_with -> Disease)
    graph['ncRNA', 'associated_with', 'Disease'].edge_index = edge_index
    graph['ncRNA', 'associated_with', 'Disease'].edge_weight = edge_weight
    graph['ncRNA', 'associated_with', 'Disease'].confidence_score = conf_score
    
    # Reverse edges: (Disease -> rev_associated_with -> ncRNA)
    rev_edge_index = torch.stack([edge_index[1], edge_index[0]], dim=0)
    graph['Disease', 'rev_associated_with', 'ncRNA'].edge_index = rev_edge_index
    graph['Disease', 'rev_associated_with', 'ncRNA'].edge_weight = edge_weight
    graph['Disease', 'rev_associated_with', 'ncRNA'].confidence_score = conf_score
    
    print("\n" + "="*50)
    print("PYG HETERODATA GRAPH SUMMARY:")
    print(graph)
    print(f"Edge Index Shape: {edge_index.shape}")
    print(f"Edge Weights Range: [{edge_weight.min():.4f}, {edge_weight.max():.4f}]")
    print(f"Confidence Scores Range: [{conf_score.min():.4f}, {conf_score.max():.4f}]")
    print("="*50)
    
    torch.save(graph, out_graph_pt)
    print(f"Saved HeteroData graph object to {out_graph_pt}")
    return graph

if __name__ == "__main__":
    build_hetero_graph()
