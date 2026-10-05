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
- **Node Resolution Engine (`match_ncrna_sequences.py`)**: Multi-key sequence indexing resolved **20,973 verified human ncRNA nodes** across biotypes (74.1% miRNA, 14.6% lncRNA, 6.5% circRNA, 3.4% piRNA, 0.6% snoRNA, 0.03% tRNA), retrieving true sequences (clamped to 1,022 nt).
- **Disease Ontology & Inheritance Engine (`build_asmsg_dataset.py`)**: Dataset-wide cross-row DO ID & MeSH inheritance increased DO ID row coverage to **211,120 association rows (89.95%)** and resolved raw disease strings into **2,749 Disease Ontology nodes** (965 DO IDs, 345 MeSH IDs, 1,439 clean string phenotypes).
- **Continuous Edge Weighting**: Log-scaled **PMID evidence weights** \( w_{ij} \in (0, 1] \) across **132,021 unique edges** (max evidence = 257 PMIDs/edge; 71.3% single PMID).
- **Output Artifacts**: Exported master clean dataset to `datasets/asmsg_clean/` (`asmsg_nodes.csv`, `asmsg_diseases.csv`, `asmsg_edges.csv`).
- **Triage Additions Pending**: `CURATION_POLICY.md` (Task 1.4), adding RNADisease `score` column (Task 1.5), transparent per-biotype reporting (Task 1.6).

### Phase 2: Multimodal Embedding & Ablation Engine (⏳ NEXT UP)
- **ncRNA Embeddings**: Execute `src/features/run_rna_fm.py` to extract **640-dim RNA-FM embeddings** for 20,973 ncRNAs (`batch_size=16`, FP16).
- **Disease Embeddings**: Construct `src/features/extract_disease_biobert.py` for **768-dim BioBERT embeddings** and **768-dim SapBERT embeddings** (Task 2.2b ablation).
- **Drug Embeddings**: Construct `src/features/extract_drug_chemberta.py` for **384-dim ChemBERTa-2 embeddings**.
- **Ablation Feature Engine (Task 2.4)**: Extract 3-mer frequency vectors, one-hot biotype vectors, sequence length scalars, and random init baselines.

### Phase 3: PyG Heterogeneous Graph Builder (📋 PLANNED)
- Construct PyTorch Geometric `HeteroData` object with multi-entity typed nodes (`ncRNA`, `Disease`, `Drug`), edge indices, PMID evidence weights, and RNADisease confidence scores.

### Phase 4: Adaptive Dual-View GNN & Contrastive Learning (📋 PLANNED)
- Implement `src/models/asmsg_gnn.py` with learned confidence-aware edge-denoising gating MLP, dual-view GNN message passing, InfoNCE contrastive loss, and margin triplet loss.

### Phase 5: Leakage-Controlled Evaluation & Rigorous Benchmark Suite (📋 PLANNED)
- **Leakage Controls (Task 5.1/5.1b)**: Seed-family/CD-HIT 80% sequence identity ncRNA clustering & DO/MeSH subtree hierarchy-aware disease splits.
- **4 Evaluation Quadrants (Task 5.2)**: Transductive S-S (Seen-Seen), Inductive U-S (Unseen ncRNA), Inductive S-U (Unseen Disease), Inductive U-U (Unseen-Unseen).
- **Ranking & Classification Metrics (Task 5.2b)**: Hits@10, Hits@50, MRR, AUROC, AUPR.
- **5 Trivial Baselines (Task 5.3)**: Degree-only, Biotype-only, Sequence-length-only, k-NN-on-3-mer, Random.
- **Reimplemented Baselines (Task 5.4)**: Reimplement SSCLMD, SSLGRDA, GSLRDA, MIFNDRA, DMGAT architectures on identical ASMSG splits.
- **Statistical Protocol (Task 5.5)**: 5 random seeds, mean ± std, paired t-test / Wilcoxon signed-rank test.
- **Temporal Validation (Task 5.6)**: Split edges by PMID publication year before/after year T.

---

## 3. Tech Stack & Dependencies (Installed in venv)
- `torch`, `torch-geometric` (PyG)
- `transformers` (for ChemBERTa-2, BioBERT, SapBERT)
- `pandas`, `numpy`, `scikit-learn`, `openpyxl`
- `rna-fm` (for sequence embeddings)

---

## 4. Benchmark Methodological Summary & Limitations vs 5 Baselines

| Method | ncRNA Scope | Disease Scope | Feature Representation | Edge Topology | Evaluation Paradigm & Leakage Controls |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **SSCLMD** (2023) | miRNA (853) | 591 raw strings | GIP + 3-mer frequency | 5,424 binary | Transductive 5-fold CV (Leaky) |
| **SSLGRDA** (2024) | miRNA+lncRNA (1002) | 590 raw strings | Binary topology matrix | 9,122 binary | Transductive 5-fold CV (Leaky) |
| **GSLRDA** (2024) | miRNA+lncRNA (852) | 591 raw strings | GIP + Wang semantic | 5,424 binary | Transductive 5-fold CV (Leaky) |
| **MIFNDRA** (2023) | miRNA (788) | 441 raw strings | GIP + sequence length | 5,149 binary | Transductive 5-fold CV (Leaky) |
| **DMGAT** (2024) | miRNA (853) | 591 raw strings | Dual GAT on GIP | 5,424 binary | Transductive 5-fold CV (Leaky) |
| **ASMSG Clean (Ours)** | **Spanning 6 Biotypes (20,973 nodes; 74.1% miRNA)** | **2,749 DO/MeSH Nodes (89.95% DOID row coverage)** | **RNA-FM (640d) + BioBERT (768d) + ChemBERTa-2 (384d) + Ablations** | **132,021 edges w/ PMID weights \( w_{ij} \) & confidence scores** | **Identity-Clustered & Hierarchy-Aware Disjoint Inductive Cold-Start (S-S, U-S, S-U, U-U) + Temporal Validation** |
