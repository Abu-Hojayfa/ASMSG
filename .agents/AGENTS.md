# 🧬 ASMSG Agent Memory & Context

This file is maintained by agents to store critical state, paths, decisions, and system constraints. **Do not deviate from these rules and paths.**

---

## 1. Project Constraints & Execution Rules
- **Virtual Environment**: ALWAYS use the virtual environment located at `d:\fydp\venv`. When running Python scripts or pip installs, use `d:\fydp\venv\Scripts\python.exe` and `d:\fydp\venv\Scripts\pip.exe`.
- **Workspace Scope**: All work is restricted strictly to `d:\fydp` and its subfolders. Do NOT write files outside this directory.
- **Code Structure**: Store Python source code in `d:\fydp\src\`. Data processing scripts are in `d:\fydp\src\data\`, feature extractors in `d:\fydp\src\features\`, graph construction in `d:\fydp\src\graph\`, models in `d:\fydp\src\models\`, evaluation in `d:\fydp\src\evaluation\`.
- **Agent Mandatory Maintenance Rule**: Whenever a project phase is updated, completed, or revised, the acting agent MUST immediately update `AGENTS.md`, `Dataset_Details.md`, and `ASMSG_Action_Plan.md` to keep context ground truth 100% synchronized and eliminate chat hallucination.

---

## 2. Current Project State & Milestones

### Phase 1: Data Engine & Topology Overhaul (✅ COMPLETED & REVISED)
- **Source Dataset**: RNADisease v4.0 (Gold Standard) + RNAcentral + miRBase v22 + GENCODE v44 + Ensembl 110 + circBase + piRBase.
- **Node Resolution Engine (`match_ncrna_sequences.py`)**: Multi-key sequence indexing resolved **20,973 verified human ncRNA nodes** across 6 biotypes (miRNAs, lncRNAs, circRNAs, piRNAs, snoRNAs, tRNAs), retrieving true sequences (clamped to 1,022 nt).
- **Disease Ontology & Inheritance Engine (`build_asmsg_dataset.py`)**: Dataset-wide cross-row DO ID & MeSH inheritance increased DO ID row coverage to **211,120 association rows (89.95%)** and resolved raw disease strings into **2,749 canonical Disease Ontology nodes** (965 DO IDs, 345 MeSH IDs, 1,439 clean phenotypes).
- **Continuous Edge Weighting**: Computed log-scaled continuous **PMID literature evidence weights** \( w_{ij} \in (0, 1] \) across **132,021 unique edges** (derived from 234,698 literature association records; max evidence = 257 PMIDs/edge).
- **Output Artifacts**: Exported master clean dataset to `datasets/asmsg_clean/` (`asmsg_nodes.csv`, `asmsg_diseases.csv`, `asmsg_edges.csv`).

### Phase 2: Multimodal Embedding Extraction Engine (⏳ NEXT UP)
- **ncRNA Embeddings**: Run `src/features/run_rna_fm.py` to extract **640-dim RNA-FM embeddings** for 20,973 ncRNAs using chunked FP16 batching (`batch_size=16`).
- **Disease Embeddings**: Construct `src/features/extract_disease_biobert.py` to extract **768-dim BioBERT embeddings** for 2,749 DO ID terms.
- **Drug Embeddings**: Construct `src/features/extract_drug_chemberta.py` to extract **384-dim ChemBERTa-2 embeddings** for PubChem drug SMILES.

### Phase 3: PyG Heterogeneous Graph Builder (📋 PLANNED)
- Construct PyTorch Geometric `HeteroData` object containing multi-entity typed nodes (`ncRNA`, `Disease`, `Drug`) and edge indices with continuous PMID weights.

### Phase 4: Adaptive Dual-View GNN & Contrastive Learning (📋 PLANNED)
- Implement `src/models/asmsg_gnn.py` featuring learned confidence-aware edge-denoising gating MLP, dual-view GNN message passing, InfoNCE contrastive loss, and margin triplet loss.

### Phase 5: Inductive Cold-Start Evaluation & Benchmark Suite (📋 PLANNED)
- Partition graph into 3 evaluation regimes: Transductive 5-Fold CV, Inductive Cold-Start ncRNA, and Inductive Cold-Start Disease splits.
- Benchmark directly against SSCLMD, SSLGRDA, GSLRDA, MIFNDRA, and DMGAT baselines.

---

## 3. Tech Stack & Dependencies (Installed in venv)
- `torch`, `torch-geometric` (PyG)
- `transformers` (for ChemBERTa-2 and BioBERT)
- `pandas`, `numpy`, `scikit-learn`, `openpyxl`
- `rna-fm` (for sequence embeddings)

---

## 4. Benchmark Novelty Summary vs. 5 Baselines

| Baseline Paper | ncRNA Scope | Disease Scope | Feature Representation | Edge Topology | Evaluation Paradigm |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **SSCLMD** (2023) | miRNA (853) | 591 strings | GIP + 3-mer frequency | 5,424 binary | Transductive 5-fold CV |
| **SSLGRDA** (2024) | miRNA+lncRNA (1002) | 590 strings | Binary topology matrix | 9,122 binary | Transductive 5-fold CV |
| **GSLRDA** (2024) | miRNA+lncRNA (852) | 591 strings | GIP + Wang semantic | 5,424 binary | Transductive 5-fold CV |
| **MIFNDRA** (2023) | miRNA (788) | 441 strings | GIP + sequence length | 5,149 binary | Transductive 5-fold CV |
| **DMGAT** (2024) | miRNA (853) | 591 strings | Dual GAT on GIP | 5,424 binary | Transductive 5-fold CV |
| **ASMSG Clean (Ours)** | **All 6 Biotypes (20,973)** | **2,749 DO IDs (89.95% DOID row coverage)** | **RNA-FM (640d) + BioBERT (768d) + ChemBERTa-2 (384d)** | **132,021 edges w/ PMID weights \( w_{ij} \)** | **Disjoint Inductive Cold-Start** |
