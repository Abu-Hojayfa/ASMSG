import os
import pandas as pd
import numpy as np
import torch
from torch_geometric.data import HeteroData

# 1. Configuration & Paths
FEATURES_DIR = 'd:/fydp/features'
DATA_DIR = 'd:/fydp/datasets'

# File paths
HMDD_FILE = os.path.join(DATA_DIR, 'hmdd', 'alldata_v4.xlsx')
NCRNADRUG_FILE = os.path.join(DATA_DIR, 'ncrnadrug', 'DR_Curated.xlsx')
LNCRNA_FILE = os.path.join(DATA_DIR, 'lncrnadisease', 'website_alldata.tsv')

RNA_FEAT_FILE = os.path.join(FEATURES_DIR, 'rna_features_rnafm.npy')
RNA_NAMES_FILE = os.path.join(FEATURES_DIR, 'rna_names_order.csv')

DRUG_FEAT_FILE = os.path.join(FEATURES_DIR, 'drug_features_chemberta.npy')
DRUG_NAMES_FILE = os.path.join(FEATURES_DIR, 'drug_names_order.csv')

print("1. Loading Node Lists & Precomputed Features...")
# Load precomputed RNA and Drug lists (so their node IDs exactly match the embedding matrix rows)
rna_df = pd.read_csv(RNA_NAMES_FILE)
drug_df = pd.read_csv(DRUG_NAMES_FILE)

# Create dictionaries to quickly map string names to integer IDs
rna_map = {name.lower(): i for i, name in enumerate(rna_df['RNA_Name'])}
drug_map = {name.lower(): i for i, name in enumerate(drug_df['Drug_Name'])}

# Load embeddings
rna_features = torch.tensor(np.load(RNA_FEAT_FILE), dtype=torch.float)
drug_features = torch.tensor(np.load(DRUG_FEAT_FILE), dtype=torch.float)

print(f"Loaded {len(rna_map)} RNAs and {len(drug_map)} Drugs with features.")

# For diseases, we extract all unique disease names from HMDD and LncRNADisease
hmdd_df = pd.read_excel(HMDD_FILE)
hmdd_diseases = set(hmdd_df['disease'].dropna().astype(str).str.lower())

lnc_df = pd.read_csv(LNCRNA_FILE, sep='\t')
lnc_diseases = set(lnc_df['Disease Name'].dropna().astype(str).str.lower())

all_diseases = hmdd_diseases.union(lnc_diseases)
disease_list = list(all_diseases)
disease_map = {name: i for i, name in enumerate(disease_list)}
print(f"Extracted {len(disease_map)} unique Diseases.")

# 2. Extract Edges (Connections)
print("\n2. Building Graph Edges...")
rna_disease_src = []
rna_disease_dst = []

# Helper to find RNA ID with flexible matching (handles -3p, -5p suffixes)
def get_rna_id(rna_name, rna_map):
    if rna_name in rna_map:
        return rna_map[rna_name]
    # Try finding an RNA in the map that starts with this name (e.g. hsa-mir-21 -> hsa-mir-21-5p)
    for mr, idx in rna_map.items():
        if mr.startswith(rna_name) or rna_name.startswith(mr):
            return idx
    return None

# Process HMDD (RNA -> Disease)
for _, row in hmdd_df.iterrows():
    rna = str(row['miRNA']).lower()
    disease = str(row['disease']).lower()
    rna_id = get_rna_id(rna, rna_map)
    if rna_id is not None and disease in disease_map:
        rna_disease_src.append(rna_id)
        rna_disease_dst.append(disease_map[disease])

# Process LncRNADisease (RNA -> Disease)
for _, row in lnc_df.iterrows():
    rna = str(row['ncRNA Symbol']).lower()
    disease = str(row['Disease Name']).lower()
    rna_id = get_rna_id(rna, rna_map)
    if rna_id is not None and disease in disease_map:
        rna_disease_src.append(rna_id)
        rna_disease_dst.append(disease_map[disease])

# Process ncRNADrug (RNA -> Drug)
ncrnadrug_df = pd.read_excel(NCRNADRUG_FILE)
rna_drug_src = []
rna_drug_dst = []

for _, row in ncrnadrug_df.iterrows():
    rna = str(row['ncRNA_Name']).lower()
    drug = str(row['Drug_Name']).lower()
    rna_id = get_rna_id(rna, rna_map)
    if rna_id is not None and drug in drug_map:
        rna_drug_src.append(rna_id)
        rna_drug_dst.append(drug_map[drug])

print(f"Found {len(rna_disease_src)} RNA-Disease edges.")
print(f"Found {len(rna_drug_src)} RNA-Drug edges.")

# 3. Construct PyG HeteroData Object
print("\n3. Constructing PyTorch Geometric HeteroData...")
data = HeteroData()

# Add node features
data['rna'].x = rna_features
data['drug'].x = drug_features
# For diseases, if we don't have embeddings yet, we use a random feature or one-hot for now.
# Let's initialize with an empty feature matrix of size (num_diseases, 128) just so the GNN can train embeddings for them.
data['disease'].x = torch.randn((len(disease_map), 128)) 

# Add edge indices (Format: [2, num_edges])
edge_index_rna_disease = torch.tensor([rna_disease_src, rna_disease_dst], dtype=torch.long)
edge_index_rna_drug = torch.tensor([rna_drug_src, rna_drug_dst], dtype=torch.long)

data['rna', 'associates_with', 'disease'].edge_index = edge_index_rna_disease
data['rna', 'targets', 'drug'].edge_index = edge_index_rna_drug

# Save the graph
out_graph_path = os.path.join(FEATURES_DIR, 'asmsg_hetero_graph.pt')
torch.save(data, out_graph_path)

print(f"\n✅ Successfully built the Heterogeneous Graph!")
print(f"Saved to: {out_graph_path}")
print(data)
