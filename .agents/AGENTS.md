# ASMSG Agent Memory & Context

This file is maintained by agents to store critical state, paths, and decisions. **Do not deviate from these rules and paths.**

## 1. Project Constraints
- **Virtual Environment**: ALWAYS use the virtual environment located at `d:\fydp\venv`. When running Python scripts or pip installs, use `d:\fydp\venv\Scripts\python.exe` and `d:\fydp\venv\Scripts\pip.exe`.
- **Workspace**: All work is restricted to `d:\fydp` and its subfolders. Do NOT write files outside this directory.
- **Code Structure**: Store Python source code in `d:\fydp\src\`. Feature extractors are in `d:\fydp\src\features\`, graph building logic in `d:\fydp\src\graph\`, etc.

## 2. Current Project State
- **Phase 1: Data Collection (⚠️ Needs Update)**
  - Existing datasets: HMDD v4.0, ncRNADrug, miRBase FASTA, and PubChem drug SMILES.
  - **NEW GOAL**: Integrate **RNADisease v4.0** (to replace/supplement HMDD and LncRNADisease) and download bulk sequences from **RNAcentral** to achieve TRUE ncRNA coverage (miRNAs, lncRNAs, circRNAs, piRNAs, snoRNAs).
- **Phase 2: Multimodal Feature Extraction (⚠️ Needs Update)**
  - Current: 640-dim embeddings for 1,203 miRNAs using RNA-FM.
  - **NEW GOAL**: Expand to extract RNA-FM embeddings for lncRNAs and circRNAs (truncating >1024nt seqs).
  - Current: 384-dim embeddings for Drugs using ChemBERTa-2 (Completed).
- **Phase 3: Heterogeneous Graph Construction (⚠️ Needs Update)**
  - Current PyG `HeteroData` is miRNA-biased. Needs rebuilding to include lncRNA/circRNA nodes and their respective edges, plus disease name normalization.
- **Phase 4: Core ASMSG Model Implementation (⏳ In Progress)**
  - Writing the Dual-View Graph Neural Network architecture in PyTorch, adapting contrastive learning for scaled node counts.

## 3. Tech Stack & Dependencies (Installed in venv)
- `torch`, `torch-geometric` (PyG)
- `transformers` (for ChemBERTa-2)
- `pandas`, `numpy`, `scikit-learn`, `openpyxl`
- `rna-fm` (for RNA sequences)

## 4. Graph Architecture Notes
- **Nodes**: ncRNA (multi-type: miRNA, lncRNA, circRNA), Disease, Drug.
- **Edges**: ncRNA-Disease (from RNADisease v4.0, HMDD, LncRNADisease), ncRNA-Drug (from ncRNADrug).
- **Features**: RNA-FM embeddings for ALL ncRNAs (sourced from miRBase + RNAcentral), ChemBERTa-2 embeddings for Drugs, Semantic similarities for Diseases.
