import pandas as pd
import numpy as np
import torch
from transformers import AutoTokenizer, AutoModel
import os
import time

# 1. Configuration
SMILES_FILE = 'd:/fydp/datasets/drugbank/drug_smiles.csv'
OUTPUT_DIR = 'd:/fydp/features'
# We use ChemBERTa-77M which naturally outputs a 384-dimensional embedding
MODEL_NAME = 'DeepChem/ChemBERTa-77M-MLM' 

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

print(f"Loading dataset from {SMILES_FILE}...")
df = pd.read_csv(SMILES_FILE)

# Filter out drugs that we couldn't find SMILES for
valid_drugs = df[df['Found'] == True].copy()
valid_drugs.reset_index(drop=True, inplace=True)
print(f"Found {len(valid_drugs)} drugs with valid SMILES strings.")

print(f"Loading {MODEL_NAME} from HuggingFace (this requires internet)...")
tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
model = AutoModel.from_pretrained(MODEL_NAME)

# Move model to GPU if available
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = model.to(device)
model.eval()
print(f"Model loaded successfully on {device}!")

embeddings = []

print("Extracting features... (this might take a minute)")
start_time = time.time()

# 2. Extract Embeddings
with torch.no_grad():
    for i, smiles in enumerate(valid_drugs['SMILES']):
        # Tokenize the SMILES string
        inputs = tokenizer(str(smiles), return_tensors='pt', padding=True, truncation=True, max_length=128)
        inputs = {k: v.to(device) for k, v in inputs.items()}
        
        # Pass through the model
        outputs = model(**inputs)
        
        # We take the embedding of the [CLS] token (the first token), which represents the whole sequence
        cls_embedding = outputs.last_hidden_state[0, 0, :].cpu().numpy()
        embeddings.append(cls_embedding)
        
        if (i + 1) % 100 == 0:
            print(f"Processed {i + 1}/{len(valid_drugs)} drugs...")

# Convert list of arrays to a single NumPy matrix
embeddings_matrix = np.array(embeddings)
print(f"\nFinished in {time.time() - start_time:.2f} seconds.")
print(f"Final embeddings matrix shape: {embeddings_matrix.shape}") 
# Expected shape: (Number of valid drugs, 384)

# 3. Save the results
out_matrix_path = os.path.join(OUTPUT_DIR, 'drug_features_chemberta.npy')
out_names_path = os.path.join(OUTPUT_DIR, 'drug_names_order.csv')

np.save(out_matrix_path, embeddings_matrix)
valid_drugs[['Drug_Name']].to_csv(out_names_path, index=False)

print(f"✅ Saved feature matrix to: {out_matrix_path}")
print(f"✅ Saved corresponding drug names to: {out_names_path}")
