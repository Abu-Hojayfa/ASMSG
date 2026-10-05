# 📊 ASMSG Datasets: Comprehensive Specification & Benchmark Analysis

This document provides the definitive specification of the **ASMSG Clean** master dataset, its empirical metrics, data processing pipeline, and comparative analysis against key baseline papers.

---

## 1. Comparative Novelty & Methodological Scope vs. 5 Baseline Papers

Existing benchmarks (SSCLMD, SSLGRDA, GSLRDA, MIFNDRA, DMGAT) derive from small, inherited miRNA-dominated datasets evaluated predominantly under standard transductive 5-fold cross-validation. ASMSG Clean expands the topological and biotype scope while introducing leakage-controlled inductive cold-start evaluation:

| Feature / Metric | **SSCLMD** (2023) | **SSLGRDA** (2024) | **GSLRDA** (2024) | **MIFNDRA** (2023) | **DMGAT** (2024) | **ASMSG Clean (Ours — Option B)** |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Primary Journal** | *IEEE/ACM TCBB* | *Brief. Bioinform.* | *Comput. Biol. Med.* | *Comput. Biol. Med.* | *IEEE JBHI* | **This Work** |
| **Biotype Scope** | miRNA only | miRNA + lncRNA | miRNA + lncRNA | miRNA only | miRNA only | **Spanning 6 Biotypes (74.1% miRNA, 14.6% lncRNA, 6.5% circRNA, 3.4% piRNA, 0.6% snoRNA, 0.8% other_ncRNA)** |
| **Sequence-Matched Nodes** | 853 | 1,002 | 852 | 788 | 853 | **20,559 sequence-matched nodes** |
| **Disease Scope** | 591 raw strings | 590 raw strings | 591 raw strings | 441 raw strings | 591 raw strings | **1,304 Standardized Disease Nodes (965 DO IDs, 339 MeSH IDs; 100% Ontology Curated, 0 Raw Strings)** |
| **Unique Edge Count** | 5,424 binary (0/1) | 9,122 binary (0/1) | 5,424 binary (0/1) | 5,149 binary (0/1) | 5,424 binary (0/1) | **121,768 unique edges (derived from 234,698 literature records)** |
| **Edge Weight Scheme** | Unweighted binary | Unweighted binary | Unweighted binary | Unweighted binary | Unweighted binary | **Composite Dual Weights \( w_{ij} = 0.5 w_{\text{PMID}} + 0.5 S_{\text{score}} \) ($w_{ij} \in [0.2269, 0.9924]$)** |
| **Sequence Feature Vectors** | k-mer (3-mer) | None | None | Sequence length | None | **RNA-FM Contextual Embeddings (640d, Head-Tail Clamping) + 3-mer / One-hot / Length Ablations** |
| **Disease Feature Vectors** | MeSH tree | None | Disease semantic | MeSH semantic | None | **BioBERT (768d) + SapBERT (768d) Semantic Embeddings** |
| **Evaluation Paradigm** | Transductive 5-fold | Transductive 5-fold | Transductive 5-fold | Transductive 5-fold | Transductive 5-fold | **Identity-Clustered & Hierarchy-Aware Disjoint Inductive Cold-Start (S-S, U-S, S-U, U-U) + Temporal Validation** |

---

## 2. Empirical Graph Statistics & Distributions (`datasets/asmsg_clean/`)

An empirical analysis of the exported Option B master dataset reveals key structural properties:

### **A. Biotype Distribution**
To prevent ultra-sparse biotypes (e.g., tRNA with only 42 edges) from creating unstable GNN message-passing channels, rare biotypes (`tRNA`, `snRNA`, `scRNA`, `rRNA`, `pseudo`) are standardized into a consolidated `other_ncRNA` class while preserving primary sequence strings:
- **lncRNA:** 6,211 nodes
- **circRNA:** 5,908 nodes
- **piRNA:** 3,888 nodes
- **miRNA:** 3,861 nodes
- **other_ncRNA:** 376 nodes
- **snoRNA:** 315 nodes

### **B. Composite Dual Edge Weight Scheme**
To resolve flatlining where 71.3% of edges possess 1 PMID, edge weights combine PubMed literature popularity AND experimental confidence scores $S_{\text{score}} \in [0.3290, 1.0000]$:
$$w_{ij} = 0.5 \cdot \frac{\log(1 + E_{ij})}{\max \log(1 + E_{ij})} + 0.5 \cdot S_{\text{score}}(i, j)$$
This produces continuous edge weights ($w_{ij} \in [0.2269, 0.9924]$, mean $0.3772$) across all 121,768 edges.

### **C. Inter-ncRNA Sequence Similarity Topology**
To resolve degree skew where 65.3% of ncRNAs connect to a single disease, ASMSG constructs supplementary `(ncRNA, sequence_similar_to, ncRNA)` topological edges (`ncrna_sequence_similarity_edges.csv`) based on 3-mer k-mer frequency cosine similarity ($\ge 0.85$). This connects leaf ncRNAs to sequence-similar neighbors.

---

## 3. ASMSG Master Dataset Specifications (`datasets/asmsg_clean/`)

### **`asmsg_nodes.csv` (ncRNA Master Table)**
* **Count:** **20,559 sequence-matched human ncRNA nodes**
* **Source Databases:** miRBase v22, GENCODE v44, Ensembl 110, circBase, piRBase.
* **Sequence Clamping:** Head-Tail dual-window clamping (first 511 nt + last 511 nt) for transcripts $>1,022$ nt.

### **`asmsg_diseases.csv` (Disease Master Table — 100% Curated)**
* **Count:** **1,304 Standardized Disease Ontology & MeSH nodes (0 Unmapped Strings)**
* **Resolution Engine:** Cross-row DO ID & MeSH synonym inheritance mapped 1,761 DO IDs and 1,862 MeSH IDs, raising ontology row coverage to **100%**.

### **`asmsg_edges.csv` (Continuous Evidence Edge Table)**
* **Count:** **121,768 unique association edges**
* **Columns:** `RNA Symbol`, `Disease_ID`, `Evidence_Count`, `PMID_List`, `Confidence_Score`, `Edge_Weight`.

### **`ncrna_sequence_similarity_edges.csv` (ncRNA Similarity Topology Table)**
* **Columns:** `RNA_Symbol_1`, `RNA_Symbol_2`, `Similarity_Score` ($\ge 0.85$).

---

## 4. Methodological Rationale & Reviewer Defense

### **A. Sequence Retrieval Rate (20,973 / 61,947)**
Multi-key sequence indexing matched 20,973 ncRNAs (**33.86% matching rate**, retaining **81.8% of total association records**). Unmatched candidate symbols are unannotated high-throughput transcript IDs or deprecated aliases that lack entry in gold-standard genomic reference databases. Retaining nodes without sequence data would force zero-padding or random noise initialization, corrupting feature space.

### **B. Disease Entity Normalization & 89.95% DOID Row Coverage**
Cross-row inheritance mapped raw disease string variants (*"breast cancer"*, *"breast carcinoma"*, *"breast malignant neoplasm"*) to canonical Disease Ontology IDs (`DOID:1612`), resolving string fragmentation and boosting DO ID coverage from 68% to 89.95% across all 234,698 literature rows.

### **C. Negative Sampling Strategy**
RNADisease v4 contains exclusively positive association records. To mitigate false-negative bias, ASMSG uses evidence-weighted negative sampling, sampling unobserved pairs primarily from entity pairs with zero PubMed literature co-mentions and low topological overlap.

---

## 5. Supplementary Datasets

### **`datasets/ncrnadrug/DR_Curated.xlsx` (ncRNA-Drug Associations)**
Curated ncRNA-drug resistance and target interactions used to introduce **Drug nodes** into the heterogeneous graph topology.

### **`datasets/drugbank/drug_smiles.csv` (PubChem Drug SMILES)**
Canonical SMILES chemical structure strings for drug nodes, processed via **ChemBERTa-2** (`DeepChem/ChemBERTa-77M-MTR`) to produce 384-dimensional chemical feature representations.

---

## 6. Benchmarking Strategy: Reimplementation Protocol

Legacy binary matrices (`rda.csv`) have been removed to prevent data contamination. Benchmarking is conducted by reimplementing baseline model architectures (SSCLMD, SSLGRDA, GSLRDA, MIFNDRA, DMGAT) and evaluating them under identical train/test split conditions on the ASMSG Clean dataset using 5 random seeds, reporting mean ± std and statistical significance tests.
