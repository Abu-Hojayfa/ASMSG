import os
import json
import torch
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score

def evaluate_trivial_baselines(seed=42):
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    split_file = os.path.join(base_dir, 'datasets', 'asmsg_clean', f'asmsg_splits_seed{seed}.json')
    nodes_csv = os.path.join(base_dir, 'datasets', 'asmsg_clean', 'asmsg_nodes.csv')
    diseases_csv = os.path.join(base_dir, 'datasets', 'asmsg_clean', 'asmsg_diseases.csv')
    edges_csv = os.path.join(base_dir, 'datasets', 'asmsg_clean', 'asmsg_edges.csv')
    kmer_pt = os.path.join(base_dir, 'datasets', 'asmsg_clean', 'ablation_3mer_features.pt')
    
    if not os.path.exists(split_file):
        raise FileNotFoundError(f"Split file for seed {seed} not found at {split_file}")
        
    with open(split_file, 'r', encoding='utf-8') as f:
        split_data = json.load(f)
        
    nodes_df = pd.read_csv(nodes_csv)
    diseases_df = pd.read_csv(diseases_csv)
    edges_df = pd.read_csv(edges_csv)
    
    # Load 3-mer k-mer features if available
    kmer_feats = torch.load(kmer_pt).numpy() if os.path.exists(kmer_pt) else None
    
    node_to_idx = {sym: i for i, sym in enumerate(nodes_df['RNA Symbol'])}
    disease_to_idx = {dis: i for i, dis in enumerate(diseases_df['Disease_ID'])}
    
    # Build training degrees and biotype association probabilities
    ss_edge_indices = split_data['edges_ss']
    train_edges = edges_df.iloc[ss_edge_indices]
    
    ncrna_degrees = train_edges['RNA Symbol'].value_counts().to_dict()
    disease_degrees = train_edges['Disease_ID'].value_counts().to_dict()
    
    # Biotype link probability
    biotype_map = nodes_df.set_index('RNA Symbol')['RNA Type'].to_dict()
    train_biotypes = train_edges['RNA Symbol'].map(biotype_map).value_counts(normalize=True).to_dict()
    
    seq_lengths = nodes_df.set_index('RNA Symbol')['Seq_Length'].to_dict()
    
    results = {}
    quadrants = {
        'S-S': split_data['edges_ss'],
        'U-S': split_data['edges_us'],
        'S-U': split_data['edges_su'],
        'U-U': split_data['edges_uu']
    }
    
    np.random.seed(seed)
    
    for quad_name, pos_indices in quadrants.items():
        if len(pos_indices) == 0:
            continue
            
        pos_df = edges_df.iloc[pos_indices]
        
        # Sample equal number of negative edges
        neg_rnas = np.random.choice(nodes_df['RNA Symbol'].values, size=len(pos_df))
        neg_dis = np.random.choice(diseases_df['Disease_ID'].values, size=len(pos_df))
        
        y_true = np.array([1] * len(pos_df) + [0] * len(neg_rnas))
        
        eval_pairs = list(zip(pos_df['RNA Symbol'], pos_df['Disease_ID'])) + list(zip(neg_rnas, neg_dis))
        
        # 1. Degree-Product Baseline
        scores_deg = [ncrna_degrees.get(r, 0) * disease_degrees.get(d, 0) for r, d in eval_pairs]
        
        # 2. Biotype-Only Baseline
        scores_biotype = [train_biotypes.get(biotype_map.get(r, 'other'), 0.0) for r, d in eval_pairs]
        
        # 3. Sequence-Length-Only Baseline
        scores_len = [seq_lengths.get(r, 0) for r, d in eval_pairs]
        
        # 4. Uniform Random Baseline
        scores_rand = np.random.rand(len(eval_pairs))
        
        results[quad_name] = {
            'n_samples': len(y_true),
            'Degree-Product': {
                'AUROC': round(float(roc_auc_score(y_true, scores_deg)), 4),
                'AUPR': round(float(average_precision_score(y_true, scores_deg)), 4)
            },
            'Biotype-Only': {
                'AUROC': round(float(roc_auc_score(y_true, scores_biotype)), 4),
                'AUPR': round(float(average_precision_score(y_true, scores_biotype)), 4)
            },
            'Sequence-Length-Only': {
                'AUROC': round(float(roc_auc_score(y_true, scores_len)), 4),
                'AUPR': round(float(average_precision_score(y_true, scores_len)), 4)
            },
            'Uniform-Random': {
                'AUROC': round(float(roc_auc_score(y_true, scores_rand)), 4),
                'AUPR': round(float(average_precision_score(y_true, scores_rand)), 4)
            }
        }
        
    out_eval_file = os.path.join(base_dir, 'datasets', 'asmsg_clean', f'trivial_baselines_eval_seed{seed}.json')
    with open(out_eval_file, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2)
        
    print(f"\n--- Trivial Baselines Evaluation Results (Seed={seed}) ---")
    for q_name, metrics in results.items():
        print(f"Quadrant [{q_name}]:")
        for b_name, b_metrics in metrics.items():
            if isinstance(b_metrics, dict):
                print(f"  * {b_name:22s} -> AUROC: {b_metrics['AUROC']:.4f} | AUPR: {b_metrics['AUPR']:.4f}")
                
    return results

if __name__ == "__main__":
    evaluate_trivial_baselines(seed=42)
