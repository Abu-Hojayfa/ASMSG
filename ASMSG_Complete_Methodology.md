# ASMSG Project: Complete Technical Methodology Report

**Project:** Adaptive Self-Supervised Multimodal Graph Learning for ncRNA–Disease Association Prediction  
**Institution:** United International University (UIU)  
**Target Outlets:** *Bioinformatics*, *Briefings in Bioinformatics*, *IEEE/ACM TCBB*

---

## Table of Contents

1. [Executive Summary & Novelty Framework](#1-executive-summary--novelty-framework)
2. [Phase 1: Multi-Biotype Data Engine & Topology Overhaul](#2-phase-1-multi-biotype-data-engine--topology-overhaul)
3. [Phase 2: Multimodal Embedding Extraction Engine](#3-phase-2-multimodal-embedding-extraction-engine)
4. [Phase 3: PyG Heterogeneous Graph Builder](#4-phase-3-pyg-heterogeneous-graph-builder)
5. [Phase 4: Adaptive Dual-View GNN & Self-Supervised Contrastive Learning](#5-phase-4-adaptive-dual-view-gnn--self-supervised-contrastive-learning)
6. [Phase 5: Inductive Cold-Start Evaluation & Benchmarking Suite](#6-phase-5-inductive-cold-start-evaluation--benchmarking-suite)
7. [Benchmark Comparison Matrix vs. SOTA Baselines](#7-benchmark-comparison-matrix-vs-sota-baselines)

---

## 1. Executive Summary & Novelty Framework

Existing State-of-the-Art (SOTA) computational models (e.g., **SSCLMD**, **SSLGRDA**, **GSLRDA**, **MIFNDRA**, **DMGAT**) suffer from three critical methodological flaws:
1. **The Transductive Bottleneck:** Over-reliance on label-derived Gaussian Interaction Profile (GIP) kernel similarities causes complete predictive collapse under cold-start inductive settings (predicting for newly discovered ncRNAs with 0 prior edges).
2. **Matrix Binarization Information Loss:** Reducing multi-study experimental confirmations into binary 0/1 edges treats a link confirmed by 200 independent PubMed publications identically to a single noisy screen.
3. **Restricted Biotype Scope:** Existing frameworks predict associations for miRNAs or lncRNAs only (~800–1,000 nodes), ignoring the vast multi-biotype non-coding interactome.

**ASMSG** resolves these bottlenecks through four key innovations:
* **Multimodal Foundation Feature Encoders:** Contextual 640-dim embeddings from **RNA-FM** for 20,973 ncRNAs, 768-dim embeddings from **BioBERT** for 2,749 Disease Ontology (DO ID) terms, and 384-dim embeddings from **ChemBERTa-2** for PubChem drug SMILES.
* **Continuous Literature Evidence Weights:** Log-scaled continuous edge weights \( w_{ij} \in (0, 1] \) derived from 234,698 literature association records in RNADisease v4.0 (up to 257 PubMed citations per edge).
* **Adaptive Dual-View Edge-Denoising GNN:** Learned gating MLP that dynamically prunes low-confidence edges during cross-view self-supervised InfoNCE contrastive learning.
* **Strict Disjoint Inductive Cold-Start Benchmark:** Disjoint train/test splits evaluating true zero-shot generalization on unseen ncRNAs and Diseases.

---

## 2. Phase 1: Multi-Biotype Data Engine & Topology Overhaul (✅ COMPLETED & PERFECTED)

### 2.1 RNADisease v4.0 Primary Backbone & Biotype Standardization
We ingested RNADisease v4.0 (343,273 raw entries) to build a unified ground truth across **6 standardized non-coding RNA biotype classes** (`lncRNA`: 6,211, `circRNA`: 5,908, `piRNA`: 3,888, `miRNA`: 3,861, `snoRNA`: 315, `other_ncRNA`: 376).

### 2.2 Multi-Key Sequence Resolution & Head-Tail Dual Window Clamping (`src/features/match_ncrna_sequences.py`)
To prevent the high node attrition of naive string matching, we built a multi-tier lookup engine across miRBase v22, GENCODE v44, Ensembl 110, circBase, and piRBase:
* **Outcome:** Successfully retrieved verified biological sequences for **20,559 human ncRNA nodes**.
* **Head-Tail Dual Window Clamping:** Sequences exceeding 1,022 nucleotides are clamped using a **Head-Tail dual-window strategy** (first 511 nt + last 511 nt). This preserves both 5'-cap seed interaction domains and 3'-UTR regulatory binding sites, overcoming simple prefix truncation loss.

### 2.3 Strict Disease Ontology (DO ID) Standardization — Option B (`src/data/build_asmsg_dataset.py`)
Dataset-wide cross-row synonym inheritance mapped 1,761 DO IDs and 1,862 MeSH IDs across database rows. Unmapped raw string phenotypes were strictly filtered out, achieving **1,304 100% ontology-curated Disease Ontology nodes** (965 DO IDs, 339 MeSH IDs, 0 unmapped strings).

### 2.4 Composite Dual Continuous Edge Weight Formulation
To eliminate edge weight flatlining where 71.3% of edges possess 1 PMID, ASMSG formulates a **composite continuous edge weight** combining PubMed literature popularity $w_{\text{PMID}}$ with experimental confidence score $S_{\text{score}} \in [0.3290, 1.0000]$:
$$w_{ij} = 0.5 \cdot \frac{\log(1 + E_{ij})}{\max \log(1 + E_{ij})} + 0.5 \cdot S_{\text{score}}(i, j)$$
yielding continuous edge weights ($w_{ij} \in [0.2269, 0.9924]$, mean $0.3773$, std $0.1605$) across **121,768 high-confidence edges**.

### 2.5 Inter-ncRNA Topological Sequence Similarity Topology
To resolve degree skew where 65.3% of ncRNAs connect to a single disease, ASMSG constructs **71,187 supplementary `(ncRNA, sequence_similar_to, ncRNA)` topological edges** (`ncrna_sequence_similarity_edges.csv`) based on 3-mer k-mer frequency cosine similarity ($\ge 0.85$). This transforms leaf ncRNAs into connected graph components capable of multi-hop GNN message passing.

---

## 3. Phase 2: Multimodal Embedding Extraction Suite (📋 PLANNED)

1. **ncRNA Sequence Embeddings (`src/features/run_rna_fm.py`):** Passes all 20,559 ncRNA sequences through **RNA-FM** in FP16 mini-batches (`batch_size=16`) to produce a 640-dimensional feature matrix.
2. **Disease Semantic Embeddings (`src/features/extract_disease_biobert.py`):** Passes DO ID definitions through **BioBERT** (`dmis-lab/biobert-base-cased-v1.2`) and **SapBERT** (`cambridgeltl/SapBERT-from-PubMedBERT-fulltext`) to produce 768-dimensional semantic matrices.
3. **Drug Chemical Embeddings (`src/features/extract_drug_chemberta.py`):** Passes PubChem SMILES strings through **ChemBERTa-2** (`DeepChem/ChemBERTa-77M-MTR`) to produce a 384-dimensional chemical matrix.
4. **Feature Ablation Vectors (`src/features/build_ablation_features.py`):** 3-mer frequency (64d), biotype one-hot (13d), sequence length scalar (1d), and random init baseline (640d).

---

## 4. Phase 3: PyG Heterogeneous Graph Builder (📋 PLANNED)

Constructs a PyTorch Geometric `HeteroData` graph object containing:
* **Nodes:** `ncRNA` (20,559 nodes, 640d), `Disease` (1,304 nodes, 768d), `Drug` (1,594 nodes, 384d).
* **Edges:** `(ncRNA, associated_with, Disease)` with dual weights \( w_{ij} \), `(ncRNA, sequence_similar_to, ncRNA)`, `(ncRNA, targets, Drug)`, `(Drug, treats, Disease)`.

---

## 5. Phase 4: Adaptive Dual-View GNN & Self-Supervised Contrastive Learning (📋 PLANNED)

### 5.1 Learned Edge-Denoising Gating MLP
Instead of random edge dropouts (GraphCL), ASMSG computes a learned gating probability \( p_{ij} \):
$$p_{ij} = \sigma \left( \mathbf{W}_g [\mathbf{h}_i \,||\, \mathbf{h}_j] + b_g \right) \cdot w_{ij}$$
Low-gating edges are masked in the secondary contrastive view.

### 5.2 Joint Multi-Task Objective
$$\mathcal{L}_{\text{total}} = \mathcal{L}_{\text{BCE}} + \lambda_1 \mathcal{L}_{\text{InfoNCE}} + \lambda_2 \mathcal{L}_{\text{Triplet}}$$

---

## 6. Phase 5: Inductive Cold-Start Evaluation & Benchmarking Suite (📋 PLANNED)

1. **Evidence-Weighted Negative Sampling:** Selects negative pairs primarily from unobserved pairs with zero PubMed co-mentions and low feature similarity.
2. **Leakage-Controlled Disjoint Evaluation Regimes:**
   * **Identity-Clustered ncRNA Cold-Start:** miRBase seed-family & 80% sequence identity clusters held out.
   * **Hierarchy-Aware Disease Cold-Start:** DO/MeSH subtrees held out.
   * **4 Quadrants:** S-S (Seen-Seen), U-S (Unseen ncRNA), S-U (Unseen Disease), U-U (Unseen-Unseen Dual Cold-Start).
3. **Temporal Validation:** Edges partitioned by PMID publication date ($T = 2020$).
4. **5-Seed Protocol & Trivial Baselines:** 5 random seeds with confidence intervals, evaluated against 5 trivial baselines and 5 SOTA baseline reimplementations.

---

## 7. Benchmark Comparison Matrix vs. SOTA Baselines

| Feature / Dimension | **SSCLMD** (2023) | **SSLGRDA** (2024) | **GSLRDA** (2024) | **MIFNDRA** (2023) | **DMGAT** (2024) | **ASMSG Clean (Ours)** |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Biotype Scope** | miRNA only | miRNA + lncRNA | miRNA + lncRNA | miRNA only | miRNA only | **Spanning 6 Biotypes (20,559 nodes)** |
| **Disease Standard** | 591 raw strings | 590 raw strings | 591 raw strings | 441 raw strings | 591 raw strings | **1,304 DO ID / MeSH Standardized Terms (100% Curated)** |
| **Unique Association Edges** | 5,424 binary | 9,122 binary | 5,424 binary | 5,149 binary | 5,424 binary | **121,768 Dual Continuous Weighted Edges** |
| **ncRNA Topology** | Bipartite only | Bipartite only | Bipartite only | Bipartite only | Bipartite only | **71,187 Inter-ncRNA Sequence Similarity Edges** |
| **Evaluation Paradigm** | Transductive 5-fold | Transductive 5-fold | Transductive 5-fold | Transductive 5-fold | Transductive 5-fold | **Disjoint Inductive Cold-Start (S-S, U-S, S-U, U-U) + Temporal Validation** |
| **Unique Edge Count** | 5,424 binary | 9,122 binary | 5,424 binary | 5,149 binary | 5,424 binary | **132,021 edges w/ PMID weights \( w_{ij} \)** |
| **Sequence Feature** | 3-mer frequency | None | None | Sequence length | None | **RNA-FM Foundation Model (640d vectors)** |
| **Disease Feature** | MeSH tree | None | Disease semantic | MeSH semantic | None | **BioBERT Language Model (768d vectors)** |
| **Evaluation Mode** | Transductive 5-fold | Transductive 5-fold | Transductive 5-fold | Transductive 5-fold | Transductive 5-fold | **Disjoint Inductive Cold-Start** |
