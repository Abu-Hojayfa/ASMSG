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

## 2. Phase 1: Multi-Biotype Data Engine & Topology Overhaul (✅ COMPLETED)

### 2.1 RNADisease v4.0 Primary Backbone
We ingested RNADisease v4.0 (343,273 raw entries) to build a unified ground truth across **6 non-coding RNA biotypes** (miRNAs, lncRNAs, circRNAs, piRNAs, snoRNAs, tRNAs).

### 2.2 Multi-Key Sequence Resolution (`src/features/match_ncrna_sequences.py`)
To prevent the high node attrition of naive string matching, we built a multi-tier lookup engine across miRBase v22, GENCODE v44, Ensembl 110, circBase, and piRBase:
* **Outcome:** Successfully retrieved verified biological sequences for **20,973 human ncRNA nodes** (up from 11,769 in legacy pipelines).
* **Token Guard:** Sequences are programmatically clamped to a maximum length of 1,022 nucleotides to respect RNA-FM's 1,024-token context window (`<cls>` and `<eos>` consume 2 slots).

### 2.3 Disease Ontology (DO ID) Standardization (`src/data/build_asmsg_dataset.py`)
Mapped 3,321 raw disease strings into **2,749 standardized Disease Ontology (DO ID) and MeSH terms**, eliminating string fragmentation.

### 2.4 Continuous PMID Evidence Weighting
Instead of binary 0/1 matrices, we calculate continuous log-scaled edge weights \( w_{ij} \in (0, 1] \) across **132,021 unique edges**:
$$w_{ij} = \frac{\log(1 + E_{ij})}{\max \log(1 + E_{ij})}$$
where \( E_{ij} \) is the count of independent literature PubMed citations supporting the edge (max 257 PMIDs/edge).

---

## 3. Phase 2: Multimodal Embedding Extraction Engine (⏳ NEXT UP)

1. **ncRNA Sequence Embeddings (`src/features/run_rna_fm.py`):** Passes all 20,973 ncRNA sequences through **RNA-FM** in FP16 mini-batches (`batch_size=16`) to produce a 640-dimensional feature matrix.
2. **Disease Semantic Embeddings (`src/features/extract_disease_biobert.py`):** Passes DO ID definitions through **BioBERT** (`dmis-lab/biobert-base-cased-v1.2`) to produce a 768-dimensional semantic matrix.
3. **Drug Chemical Embeddings (`src/features/extract_drug_chemberta.py`):** Passes PubChem SMILES strings through **ChemBERTa-2** (`DeepChem/ChemBERTa-77M-MTR`) to produce a 384-dimensional chemical matrix.

---

## 4. Phase 3: PyG Heterogeneous Graph Builder

Constructs a PyTorch Geometric `HeteroData` graph object containing:
* **Nodes:** `ncRNA` (20,973 nodes, 640d), `Disease` (2,749 nodes, 768d), `Drug` (384d).
* **Edges:** `(ncRNA, associated_with, Disease)` with continuous weights \( w_{ij} \), `(ncRNA, targets, Drug)`, `(Drug, treats, Disease)`.

---

## 5. Phase 4: Adaptive Dual-View GNN & Self-Supervised Contrastive Learning

### 5.1 Learned Edge-Denoising Gating MLP
Instead of random edge dropouts (GraphCL), ASMSG computes a learned gating probability \( p_{ij} \):
$$p_{ij} = \sigma \left( \mathbf{W}_g [\mathbf{h}_i \,||\, \mathbf{h}_j] + b_g \right) \cdot w_{ij}$$
Low-gating edges are masked in the secondary contrastive view.

### 5.2 Joint Multi-Task Objective
$$\mathcal{L}_{\text{total}} = \mathcal{L}_{\text{BCE}} + \lambda_1 \mathcal{L}_{\text{InfoNCE}} + \lambda_2 \mathcal{L}_{\text{Triplet}}$$

---

## 6. Phase 5: Inductive Cold-Start Evaluation & Benchmarking Suite

1. **Evidence-Weighted Negative Sampling:** Selects negative pairs primarily from unobserved pairs with zero PubMed co-mentions and low feature similarity.
2. **Disjoint Evaluation Regimes:**
   * Regime 1: Transductive 5-Fold CV.
   * Regime 2: Inductive Cold-Start ncRNA (20% held-out ncRNAs completely detached during training).
   * Regime 3: Inductive Cold-Start Disease (20% held-out Diseases completely detached during training).
3. **Benchmark Comparison:** Benchmarks against SSCLMD, SSLGRDA, GSLRDA, MIFNDRA, DMGAT, GCN, GATv2, and RGCN on ASMSG Clean and HMDD v3.2.

---

## 7. Benchmark Comparison Matrix vs. SOTA Baselines

| Feature / Dimension | **SSCLMD** (2023) | **SSLGRDA** (2024) | **GSLRDA** (2024) | **MIFNDRA** (2023) | **DMGAT** (2024) | **ASMSG Clean (Ours)** |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Biotype Scope** | miRNA only | miRNA + lncRNA | miRNA + lncRNA | miRNA only | miRNA only | **All 6 ncRNA Biotypes (20,973 nodes)** |
| **Disease Standard** | 591 raw strings | 590 raw strings | 591 raw strings | 441 raw strings | 591 raw strings | **2,749 DO ID / MeSH standardized terms** |
| **Unique Edge Count** | 5,424 binary | 9,122 binary | 5,424 binary | 5,149 binary | 5,424 binary | **132,021 edges w/ PMID weights \( w_{ij} \)** |
| **Sequence Feature** | 3-mer frequency | None | None | Sequence length | None | **RNA-FM Foundation Model (640d vectors)** |
| **Disease Feature** | MeSH tree | None | Disease semantic | MeSH semantic | None | **BioBERT Language Model (768d vectors)** |
| **Evaluation Mode** | Transductive 5-fold | Transductive 5-fold | Transductive 5-fold | Transductive 5-fold | Transductive 5-fold | **Disjoint Inductive Cold-Start** |
