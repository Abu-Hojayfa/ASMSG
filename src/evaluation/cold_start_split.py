import os
import json
import re
import numpy as np
import pandas as pd

def get_ncrna_cluster_key(symbol, rna_type):
    symbol_str = str(symbol).strip()
    rna_type_clean = str(rna_type).strip().lower()
    
    # 1. MicroRNA family clustering (miRBase seed family)
    if 'mir' in rna_type_clean or symbol_str.lower().startswith(('hsa-mir', 'mir-', 'hsa-let', 'let-')):
        # Extract base family (e.g., hsa-let-7a-5p -> let-7; hsa-mir-17-3p -> mir-17)
        match = re.search(r'(let-\d+|mir-\d+|mirlet-\d+)', symbol_str.lower())
        if match:
            return f"family_mir_{match.group(1)}"
        return f"family_mir_{symbol_str.split('-')[1]}" if '-' in symbol_str else f"family_mir_{symbol_str}"
        
    # 2. LncRNA gene family prefix (e.g., LINC00123 -> LINC; HOXA-AS1 -> HOXA)
    if 'lnc' in rna_type_clean:
        prefix = re.split(r'[-_0-9]', symbol_str)[0].upper()
        if len(prefix) >= 3:
            return f"family_lnc_{prefix}"
        return f"family_lnc_generic_{hash(symbol_str) % 50}"
        
    # 3. CircRNA / piRNA / snoRNA prefix family
    prefix = symbol_str.split('-')[0].upper() if '-' in symbol_str else symbol_str[:5].upper()
    return f"family_{rna_type_clean}_{prefix}"

def get_disease_cluster_key(disease_id, disease_name):
    dis_id_str = str(disease_id).strip()
    dis_name_clean = str(disease_name).strip().lower()
    
    # High-level disease subtrees
    if any(k in dis_name_clean for k in ['cancer', 'carcinoma', 'neoplasm', 'tumor', 'leukemia', 'lymphoma', 'melanoma', 'glioma', 'sarcoma', 'adenoma']):
        return "subtree_oncology"
    if any(k in dis_name_clean for k in ['cardiac', 'heart', 'vascular', 'hypertension', 'stroke', 'cardiopathy', 'atherosclerosis', 'myocardial']):
        return "subtree_cardiovascular"
    if any(k in dis_name_clean for k in ['brain', 'neuro', 'alzheimer', 'parkinson', 'sclerosis', 'seizure', 'epilepsy', 'dementia', 'huntington']):
        return "subtree_neurology"
    if any(k in dis_name_clean for k in ['liver', 'hepatic', 'hepatitis', 'cirrhosis', 'pancrea', 'biliary']):
        return "subtree_hepatology"
    if any(k in dis_name_clean for k in ['kidney', 'renal', 'nephro', 'glomerulo']):
        return "subtree_nephrology"
    if any(k in dis_name_clean for k in ['lung', 'pulmonary', 'respiratory', 'asthma', 'pneumonia', 'bronchitis', 'fibrosis']):
        return "subtree_pulmonology"
    if any(k in dis_name_clean for k in ['diabetes', 'metabolic', 'obesity', 'thyroid', 'gout']):
        return "subtree_metabolic"
    if any(k in dis_name_clean for k in ['infection', 'viral', 'bacterial', 'sepsis', 'hiv', 'flu']):
        return "subtree_infectious"
        
    # Default prefix cluster
    if dis_id_str.startswith('DOID:'):
        return f"doid_prefix_{dis_id_str[:7]}"
    if dis_id_str.startswith('MESH:'):
        return f"mesh_prefix_{dis_id_str[:7]}"
    return f"phenotype_cluster_{hash(dis_id_str) % 40}"

def create_leakage_controlled_splits(seed=42):
    np.random.seed(seed)
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    nodes_csv = os.path.join(base_dir, 'datasets', 'asmsg_clean', 'asmsg_nodes.csv')
    diseases_csv = os.path.join(base_dir, 'datasets', 'asmsg_clean', 'asmsg_diseases.csv')
    edges_csv = os.path.join(base_dir, 'datasets', 'asmsg_clean', 'asmsg_edges.csv')
    out_dir = os.path.join(base_dir, 'datasets', 'asmsg_clean')
    
    print(f"\n--- Generating Leakage-Controlled Cold-Start Splits (Seed={seed}) ---")
    nodes_df = pd.read_csv(nodes_csv)
    diseases_df = pd.read_csv(diseases_csv)
    edges_df = pd.read_csv(edges_csv)
    
    # 1. Cluster ncRNAs by sequence homology & seed family
    nodes_df['cluster_key'] = [get_ncrna_cluster_key(row['RNA Symbol'], row['RNA Type']) for _, row in nodes_df.iterrows()]
    unique_ncrna_clusters = list(nodes_df['cluster_key'].unique())
    np.random.shuffle(unique_ncrna_clusters)
    
    n_ncrna = len(unique_ncrna_clusters)
    train_ncrna_clusters = set(unique_ncrna_clusters[: int(0.8 * n_ncrna)])
    val_ncrna_clusters = set(unique_ncrna_clusters[int(0.8 * n_ncrna) : int(0.9 * n_ncrna)])
    test_ncrna_clusters = set(unique_ncrna_clusters[int(0.9 * n_ncrna) :])
    
    train_ncrna_set = set(nodes_df[nodes_df['cluster_key'].isin(train_ncrna_clusters)]['RNA Symbol'])
    val_ncrna_set = set(nodes_df[nodes_df['cluster_key'].isin(val_ncrna_clusters)]['RNA Symbol'])
    test_ncrna_set = set(nodes_df[nodes_df['cluster_key'].isin(test_ncrna_clusters)]['RNA Symbol'])
    
    # 2. Cluster Diseases by ontology subtrees
    diseases_df['cluster_key'] = [get_disease_cluster_key(row['Disease_ID'], row['Disease_Name']) for _, row in diseases_df.iterrows()]
    unique_dis_clusters = list(diseases_df['cluster_key'].unique())
    np.random.shuffle(unique_dis_clusters)
    
    n_dis = len(unique_dis_clusters)
    train_dis_clusters = set(unique_dis_clusters[: int(0.8 * n_dis)])
    val_dis_clusters = set(unique_dis_clusters[int(0.8 * n_dis) : int(0.9 * n_dis)])
    test_dis_clusters = set(unique_dis_clusters[int(0.9 * n_dis) :])
    
    train_disease_set = set(diseases_df[diseases_df['cluster_key'].isin(train_dis_clusters)]['Disease_ID'])
    val_disease_set = set(diseases_df[diseases_df['cluster_key'].isin(val_dis_clusters)]['Disease_ID'])
    test_disease_set = set(diseases_df[diseases_df['cluster_key'].isin(test_dis_clusters)]['Disease_ID'])
    
    # 3. Partition edges into 4 evaluation quadrants
    edges_ss, edges_us, edges_su, edges_uu = [], [], [], []
    
    for idx, row in edges_df.iterrows():
        rna = row['RNA Symbol']
        dis = row['Disease_ID']
        
        rna_in_train = rna in train_ncrna_set
        rna_in_test = rna in test_ncrna_set or rna in val_ncrna_set
        dis_in_train = dis in train_disease_set
        dis_in_test = dis in test_disease_set or dis in val_disease_set
        
        if rna_in_train and dis_in_train:
            edges_ss.append(idx)
        elif rna_in_test and dis_in_train:
            edges_us.append(idx)
        elif rna_in_train and dis_in_test:
            edges_su.append(idx)
        elif rna_in_test and dis_in_test:
            edges_uu.append(idx)
            
    split_info = {
        'seed': seed,
        'ncrna_counts': {
            'train': len(train_ncrna_set),
            'val': len(val_ncrna_set),
            'test': len(test_ncrna_set),
            'total': len(nodes_df)
        },
        'disease_counts': {
            'train': len(train_disease_set),
            'val': len(val_disease_set),
            'test': len(test_disease_set),
            'total': len(diseases_df)
        },
        'quadrant_edge_counts': {
            'S-S (Seen-Seen)': len(edges_ss),
            'U-S (Unseen ncRNA, Seen Disease)': len(edges_us),
            'S-U (Seen ncRNA, Unseen Disease)': len(edges_su),
            'U-U (Unseen-Unseen Dual Cold-Start)': len(edges_uu),
            'total_edges': len(edges_df)
        },
        'train_ncrnas': list(train_ncrna_set),
        'test_ncrnas': list(test_ncrna_set),
        'train_diseases': list(train_disease_set),
        'test_diseases': list(test_disease_set),
        'edges_ss': edges_ss,
        'edges_us': edges_us,
        'edges_su': edges_su,
        'edges_uu': edges_uu
    }
    
    out_split_file = os.path.join(out_dir, f'asmsg_splits_seed{seed}.json')
    with open(out_split_file, 'w', encoding='utf-8') as f:
        json.dump(split_info, f, indent=2)
        
    print(f" - Train ncRNAs: {len(train_ncrna_set)} | Test ncRNAs: {len(test_ncrna_set)}")
    print(f" - Train Diseases: {len(train_disease_set)} | Test Diseases: {len(test_disease_set)}")
    print(" - Quadrant Edges:")
    print(f"   * S-S (Seen-Seen): {len(edges_ss)} ({len(edges_ss)/len(edges_df)*100:.1f}%)")
    print(f"   * U-S (Unseen ncRNA): {len(edges_us)} ({len(edges_us)/len(edges_df)*100:.1f}%)")
    print(f"   * S-U (Unseen Disease): {len(edges_su)} ({len(edges_su)/len(edges_df)*100:.1f}%)")
    print(f"   * U-U (Dual Cold-Start): {len(edges_uu)} ({len(edges_uu)/len(edges_df)*100:.1f}%)")
    print(f"Saved split metadata to {out_split_file}")
    
    return split_info

if __name__ == "__main__":
    for s in [42, 43, 44, 45, 46]:
        create_leakage_controlled_splits(seed=s)
