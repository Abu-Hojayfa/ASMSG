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

### Phase 1: Data Engine & Topology Overhaul (⏳ IN PROGRESS — DATASET PERFECTION FOCUS)
- **Source Dataset**: RNADisease v4.0 (Gold Standard) + RNAcentral + miRBase v22 + GENCODE v44 + Ensembl 110 + circBase + piRBase.
- **Strict Disease Ontology Resolution (Option B)**: Mapped 1,761 DO IDs and 1,862 MeSH IDs, strictly filtering out unmapped strings to achieve **1,304 100% ontology-curated Disease Ontology nodes** (965 DO IDs, 339 MeSH IDs, 0 unmapped strings).
- **Sequence Resolution & Head-Tail Clamping**: Resolved **20,559 verified human ncRNA sequence nodes** across 6 biotype classes (lncRNA, circRNA, piRNA, miRNA, snoRNA, other_ncRNA) using Head-Tail dual-window clamping (first 511 nt + last 511 nt for long transcripts).
- **Dual-Weighted Continuous Associations**: Formulated composite continuous edge weights \( w_{ij} = 0.5 w_{\text{PMID}} + 0.5 S_{\text{score}} \) (\( w_{ij} \in [0.2269, 0.9924] \), mean 0.3773) across **121,768 high-confidence edges**.
- **ncRNA Topological Sequence Similarity Edges**: Generated **71,187 inter-ncRNA sequence similarity edges** (`ncrna_sequence_similarity_edges.csv`, 3-mer cosine similarity $\ge 0.85$) to connect degree-1 leaf nodes and enrich heterogeneous message passing.
- **Single Master Dataset Export**: Exported single perfected master dataset to `datasets/asmsg_clean/` (`asmsg_nodes.csv`, `asmsg_diseases.csv`, `asmsg_edges.csv`, `ncrna_sequence_similarity_edges.csv`, `CURATION_POLICY.md`).

### Phase 2: Multimodal Embedding Extraction Suite (📋 PLANNED)
- Extract RNA-FM sequence embeddings, BioBERT & SapBERT disease embeddings, ChemBERTa-2 drug SMILES embeddings, and ablation feature vectors.

### Phase 3: PyG Heterogeneous Graph Builder (📋 PLANNED)
- Assemble PyTorch Geometric `HeteroData` graph object (`asmsg_hetero_graph.pt`) containing multi-entity typed nodes (`ncRNA`, `Disease`, `Drug`), edge indices, and dual continuous edge weights.

### Phase 4: Adaptive Dual-View GNN & Contrastive Learning (📋 PLANNED)
- Implement `src/models/asmsg_gnn.py` with learned confidence-aware edge-denoising gating MLP, dual-view GNN message passing, InfoNCE contrastive loss, and margin triplet loss.

### Phase 5: Leakage-Controlled Evaluation & Benchmark Suite (📋 PLANNED)
- Implement identity-clustered & hierarchy-aware disjoint cold-start splits (S-S, U-S, S-U, U-U), 5 trivial baselines, SOTA baseline reimplementations, 5-seed statistical protocol, and temporal validation.

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
