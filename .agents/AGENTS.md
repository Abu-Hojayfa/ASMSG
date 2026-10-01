# ASMSG Agent Memory & Context

This file is maintained by agents to store critical state, paths, and decisions. **Do not deviate from these rules and paths.**

## 1. Project Constraints
- **Virtual Environment**: ALWAYS use the virtual environment located at `d:\fydp\venv`. When running Python scripts or pip installs, use `d:\fydp\venv\Scripts\python.exe` and `d:\fydp\venv\Scripts\pip.exe`.
- **Workspace**: All work is restricted to `d:\fydp` and its subfolders. Do NOT write files outside this directory.
- **Code Structure**: Store Python source code in `d:\fydp\src\`. Feature extractors are in `d:\fydp\src\features\`, graph building logic in `d:\fydp\src\graph\`, etc.

## 2. Current Project State
- **Phase 1: Data Collection (✅ Completed)**
  - All datasets successfully downloaded and stored in `d:\fydp\datasets\`.
  - HMDD v4.0, ncRNADrug, miRBase FASTA, and PubChem drug SMILES are ready.
- **Phase 2: Multimodal Feature Extraction (✅ Completed)**
  - 640-dim embeddings for RNA extracted using RNA-FM. Saved in `d:\fydp\features\rna_features_rnafm.npy`.
  - 384-dim embeddings for Drugs extracted using ChemBERTa-2. Saved in `d:\fydp\features\drug_features_chemberta.npy`.
- **Phase 3: Heterogeneous Graph Construction (✅ Completed)**
  - PyG `HeteroData` object successfully built and saved to `d:\fydp\features\asmsg_hetero_graph.pt`.
- **Phase 4: Core ASMSG Model Implementation (⏳ In Progress)**
  - Writing the Dual-View Graph Neural Network architecture in PyTorch.

## 3. Tech Stack & Dependencies (Installed in venv)
- `torch`, `torch-geometric` (PyG)
- `transformers` (for ChemBERTa-2)
- `pandas`, `numpy`, `scikit-learn`, `openpyxl`
- `rna-fm` (for RNA sequences)

## 4. Graph Architecture Notes
- **Nodes**: ncRNA, Disease, Drug.
- **Edges**: ncRNA-Disease (HMDD, LncRNADisease), ncRNA-Drug (ncRNADrug).
- **Features**: RNA-FM embeddings for RNAs, ChemBERTa-2 embeddings for Drugs, Semantic similarities for Diseases.
