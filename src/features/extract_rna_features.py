import os
import time
import numpy as np
import pandas as pd
import torch
import fm

FASTA_FILE = 'd:/fydp/datasets/miRBase/mature.fa'
HMDD_FILE = 'd:/fydp/datasets/hmdd/alldata_v4.xlsx'
NCRNADRUG_FILE = 'd:/fydp/datasets/ncrnadrug/DR_Curated.xlsx'
OUTPUT_DIR = 'd:/fydp/features'

print("Extracting unique RNAs needed for the graph...")
hmdd_df = pd.read_excel(HMDD_FILE)
ncrna_df = pd.read_excel(NCRNADRUG_FILE)

# Get unique miRNAs from the datasets (convert to lowercase for matching)
hmdd_rnas = set(hmdd_df['miRNA'].dropna().astype(str).str.lower())
ncrna_rnas = set(ncrna_df['ncRNA_Name'].dropna().astype(str).str.lower())
needed_rnas = hmdd_rnas.union(ncrna_rnas)
print(f"Total unique RNAs required for the graph: {len(needed_rnas)}")

print(f"Loading sequences from {FASTA_FILE}...")
sequences = []
names = []
with open(FASTA_FILE, 'r') as f:
    current_name = ""
    current_seq = ""
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_name:
                # Only keep if it is in our needed list
                # e.g., hsa-let-7a-5p
                if current_name.lower() in needed_rnas:
                    names.append(current_name)
                    sequences.append(current_seq)
            current_name = line.split(" ")[0][1:] 
            current_seq = ""
        else:
            current_seq += line
    if current_name and current_name.lower() in needed_rnas:
        names.append(current_name)
        sequences.append(current_seq)

print(f"Found {len(sequences)} matching miRNA sequences in miRBase.")

print("Loading RNA-FM model...")
model, alphabet = fm.pretrained.rna_fm_t12()
batch_converter = alphabet.get_batch_converter()
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = model.to(device)
model.eval()
print(f"RNA-FM loaded successfully on {device}!")

batch_size = 8
embeddings = []

print("Extracting RNA features... (this will be much faster now)")
start_time = time.time()

with torch.no_grad():
    for i in range(0, len(sequences), batch_size):
        batch_names = names[i:i+batch_size]
        batch_seqs = sequences[i:i+batch_size]
        
        data = list(zip(batch_names, batch_seqs))
        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        batch_tokens = batch_tokens.to(device)
        
        results = model(batch_tokens, repr_layers=[12])
        token_representations = results["representations"][12]
        
        for j, seq in enumerate(batch_seqs):
            seq_len = len(seq)
            seq_rep = token_representations[j, 1:seq_len+1].mean(dim=0).cpu().numpy()
            embeddings.append(seq_rep)
            
        if (i + batch_size) % 100 < batch_size:
            print(f"Processed {min(i + batch_size, len(sequences))}/{len(sequences)} RNAs...")

embeddings_matrix = np.array(embeddings)
print(f"\nFinished in {time.time() - start_time:.2f} seconds.")
print(f"Final RNA embeddings matrix shape: {embeddings_matrix.shape}") 

out_matrix_path = os.path.join(OUTPUT_DIR, 'rna_features_rnafm.npy')
out_names_path = os.path.join(OUTPUT_DIR, 'rna_names_order.csv')

np.save(out_matrix_path, embeddings_matrix)
with open(out_names_path, 'w') as f:
    f.write("RNA_Name\n")
    for name in names:
        f.write(f"{name}\n")

print("✅ Saved RNA feature matrix!")
