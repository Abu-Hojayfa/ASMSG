# 🧬 ASMSG: Adaptive Self-Supervised Multimodal Graph Learning for ncRNA–Disease Association Prediction

> **Project Status:** Phase 1 (Data Engine & Topology Overhaul) ✅ COMPLETED | Phase 2 (Multimodal Foundation Embeddings) ⏳ IN PROGRESS  
> **Target Target Outlets:** *Bioinformatics*, *Briefings in Bioinformatics*, *IEEE/ACM TCBB*

---

## 🔬 Abstract & Core Scientific Contributions

Non-coding RNAs (miRNAs, lncRNAs, circRNAs, piRNAs, snoRNAs, tRNAs) are essential regulators in cellular biology and complex human pathologies. Accurate computational prediction of ncRNA–disease associations is vital for clinical biomarker discovery and targeted therapeutics.

Existing State-of-the-Art (SOTA) computational models (e.g., **SSCLMD**, **SSLGRDA**, **GSLRDA**, **MIFNDRA**, **DMGAT**) suffer from three critical methodological flaws:
1. **The Transductive Data Leakage Trap:** Over-reliance on label-derived Gaussian Interaction Profile (GIP) kernel similarities causes complete predictive collapse under cold-start inductive settings (predicting for newly discovered ncRNAs with 0 prior edges).
2. **Matrix Binarization Information Loss:** Reducing multi-study experimental confirmations into binary 0/1 edges treats a link confirmed by 200 independent PubMed publications identically to a single noisy screen.
3. **Restricted Biotype Scope:** Existing frameworks predict associations for miRNAs or lncRNAs only (~800–1,000 nodes), ignoring the vast multi-biotype non-coding interactome.

**ASMSG** resolves these bottlenecks through four key innovations:
* **Multimodal Foundation Feature Encoders:** Contextual 640-dim embeddings from **RNA-FM** for 20,973 ncRNAs, 768-dim embeddings from **BioBERT** for 2,749 Disease Ontology (DO ID) terms, and 384-dim embeddings from **ChemBERTa-2** for PubChem drug SMILES.
* **Continuous Literature Evidence Weights:** Log-scaled continuous edge weights \( w_{ij} \in (0, 1] \) derived from 234,698 literature association records in RNADisease v4.0 (up to 257 PubMed citations per edge).
* **Adaptive Dual-View Edge-Denoising GNN:** Learned gating MLP that dynamically prunes low-confidence edges during cross-view self-supervised InfoNCE contrastive learning.
* **Strict Disjoint Inductive Cold-Start Benchmark:** Disjoint train/test splits evaluating true zero-shot generalization on unseen ncRNAs and Diseases.

---

## 📊 Dataset Specifications (`datasets/asmsg_clean/`)

| Dataset Artifact | Entity / Element | Count | Description |
| :--- | :--- | :--- | :--- |
| **`asmsg_nodes.csv`** | Verified ncRNAs | **20,973** | Sourced across miRBase v22, GENCODE v44, Ensembl 110, circBase, and piRBase. Includes full nucleotide sequences (clamped to 1,022 nt). |
| **`asmsg_diseases.csv`** | Disease Nodes | **2,749** | Standardized onto Disease Ontology (DO ID) and MeSH ontologies, resolving string fragmentation. |
| **`asmsg_edges.csv`** | Association Edges | **132,021** | Unique ncRNA–Disease edges with continuous PubMed evidence weights \( w_{ij} \in (0, 1] \) (max 257 PMIDs/edge). |

---

## 🏆 SOTA Baseline Comparison Matrix

| Feature / Dimension | **SSCLMD** (2023) | **SSLGRDA** (2024) | **GSLRDA** (2024) | **MIFNDRA** (2023) | **DMGAT** (2024) | **ASMSG Clean (Ours)** |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Biotype Scope** | miRNA only | miRNA + lncRNA | miRNA + lncRNA | miRNA only | miRNA only | **All 6 ncRNA Biotypes (20,973 nodes)** |
| **Disease Standard** | 591 raw strings | 590 raw strings | 591 raw strings | 441 raw strings | 591 raw strings | **2,749 DO ID / MeSH standardized terms** |
| **Unique Edge Count** | 5,424 binary | 9,122 binary | 5,424 binary | 5,149 binary | 5,424 binary | **132,021 edges w/ PMID weights \( w_{ij} \)** |
| **Sequence Feature** | 3-mer frequency | None | None | Sequence length | None | **RNA-FM Foundation Model (640d vectors)** |
| **Disease Feature** | MeSH tree | None | Disease semantic | MeSH semantic | None | **BioBERT Language Model (768d vectors)** |
| **Evaluation Mode** | Transductive 5-fold | Transductive 5-fold | Transductive 5-fold | Transductive 5-fold | Transductive 5-fold | **Disjoint Inductive Cold-Start** |

---

## 📁 Repository & Code Structure

```
fydp/
├── .agents/
│   └── AGENTS.md               # Master agent memory, execution rules, & ground-truth state
├── datasets/
│   ├── asmsg_clean/            # Master clean dataset (20,973 ncRNAs, 2,749 DO IDs, 132,021 edges)
│   │   ├── asmsg_nodes.csv
│   │   ├── asmsg_diseases.csv
│   │   └── asmsg_edges.csv
│   ├── rnadisease_v4/          # RNADisease v4.0 source data & sequence FASTAs
│   └── Dataset_Details.md      # Detailed dataset specifications & reviewer defense rationales
├── src/
│   ├── data/
│   │   └── build_asmsg_dataset.py  # DO ID standardization & PMID evidence weighting engine
│   ├── features/
│   │   ├── match_ncrna_sequences.py# Multi-key sequence indexing across genomic FASTAs
│   │   └── run_rna_fm.py           # RNA-FM 640-dim embedding extractor
│   ├── graph/                      # PyG HeteroData graph construction scripts
│   ├── models/                     # Dual-View GNN architecture & InfoNCE loss
│   └── evaluation/                 # Inductive cold-start split & benchmarking suite
├── ASMSG_Action_Plan.md        # Master technical roadmap & Oxford-tier action plan
└── README.md
```

---

## 🚀 Execution & Verification Protocol

### Environment Activation
All commands MUST be executed within the dedicated Python virtual environment:
```powershell
d:\fydp\venv\Scripts\python.exe --version
```

### Data Engine Pipeline (Phase 1)
```powershell
# 1. Execute Multi-Key Sequence Matching Engine
d:\fydp\venv\Scripts\python.exe d:\fydp\src\features\match_ncrna_sequences.py

# 2. Build Master Clean Dataset (DO ID Standardization & Evidence Weighting)
d:\fydp\venv\Scripts\python.exe d:\fydp\src\data\build_asmsg_dataset.py
```
