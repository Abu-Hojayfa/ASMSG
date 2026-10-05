# 📊 ASMSG Datasets: Comprehensive Specification & Benchmark Analysis

This document provides the definitive specification of the **ASMSG Clean** master dataset, its empirical metrics, data processing pipeline, and comparative analysis against key baseline papers.

---

## 1. Comparative Novelty & Methodological Scope vs. 5 Baseline Papers

Existing benchmarks (SSCLMD, SSLGRDA, GSLRDA, MIFNDRA, DMGAT) derive from small, inherited miRNA-dominated datasets evaluated predominantly under standard transductive 5-fold cross-validation. ASMSG Clean expands the topological and biotype scope while introducing leakage-controlled inductive cold-start evaluation:

| Feature / Metric | **SSCLMD** (2023) | **SSLGRDA** (2024) | **GSLRDA** (2024) | **MIFNDRA** (2023) | **DMGAT** (2024) | **ASMSG Clean (Ours)** |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Primary Journal** | *IEEE/ACM TCBB* | *Brief. Bioinform.* | *Comput. Biol. Med.* | *Comput. Biol. Med.* | *IEEE JBHI* | **This Work** |
| **Biotype Scope** | miRNA only | miRNA + lncRNA | miRNA + lncRNA | miRNA only | miRNA only | **Spanning 6 Biotypes (74.1% miRNA, 14.6% lncRNA, 6.5% circRNA, 3.4% piRNA, 0.6% snoRNA, 0.03% tRNA)** |
| **Sequence-Matched Nodes** | 853 | 1,002 | 852 | 788 | 853 | **20,973 sequence-matched nodes** |
| **Disease Scope** | 591 raw strings | 590 raw strings | 591 raw strings | 441 raw strings | 591 raw strings | **2,749 Disease Nodes (965 DO IDs, 345 MeSH IDs, 1,439 clean string phenotypes; 89.95% DOID row coverage)** |
| **Unique Edge Count** | 5,424 binary (0/1) | 9,122 binary (0/1) | 5,424 binary (0/1) | 5,149 binary (0/1) | 5,424 binary (0/1) | **132,021 unique edges (derived from 234,698 literature records)** |
| **Edge Weight Scheme** | Unweighted binary | Unweighted binary | Unweighted binary | Unweighted binary | Unweighted binary | **Continuous PMID evidence weights \( w_{ij} \in (0, 1] \) & confidence scores** |
| **Sequence Feature Vectors** | k-mer (3-mer) | None | None | Sequence length | None | **RNA-FM Contextual Embeddings (640d) + 3-mer / One-hot / Length Ablations** |
| **Disease Feature Vectors** | MeSH tree | None | Disease semantic | MeSH semantic | None | **BioBERT (768d) + SapBERT (768d) Semantic Embeddings** |
| **Evaluation Paradigm** | Transductive 5-fold | Transductive 5-fold | Transductive 5-fold | Transductive 5-fold | Transductive 5-fold | **Identity-Clustered & Hierarchy-Aware Disjoint Inductive Cold-Start (S-S, U-S, S-U, U-U) + Temporal Validation** |

---

## 2. Empirical Graph Statistics & Distributions

An empirical analysis of the exported `datasets/asmsg_clean/` dataset reveals key structural properties:

### **A. Biotype Distribution**
While the graph contains 6 biotype classes, the edge density is heavily skewed toward miRNAs:
- **miRNA:** 3,940 nodes | 97,825 edges (**74.1% of all edges**)
- **lncRNA:** 6,287 nodes | 19,272 edges (14.6%)
- **circRNA:** 6,145 nodes | 8,542 edges (6.5%)
- **piRNA:** 3,907 nodes | 4,461 edges (3.4%)
- **snoRNA:** 315 nodes | 829 edges (0.6%)
- **tRNA:** 18 nodes | 42 edges (0.03%)
- **snRNA & Others:** 361 nodes | 1,050 edges (0.8%)

### **B. Edge Weight Distribution**
- **71.3% of edges** (94,089 out of 132,021) have `Evidence_Count = 1`, taking minimum weight `0.1248`.
- **5.0% of edges** have `Evidence_Count ≥ 5`.
- **1.5% of edges** have `Evidence_Count ≥ 10`.
- Summary: `min=0.1248, median=0.1248, mean=0.1610, max=1.0000`.

### **C. Degree Distribution & Hubs**
- **65.2% of ncRNA nodes** (13,678 / 20,973) have **degree = 1** (median degree = 1.0).
- **31.6% of disease nodes** (868 / 2,749) have **degree = 1**.
- Top disease hub (*Neoplasms / Cancer*): **8,955 connected edges**.

### **D. Disease Node Normalization**
- **965 nodes** resolved to official Disease Ontology IDs (`DOID:XXXX`).
- **345 nodes** resolved to official MeSH IDs (`MESH:XXXX`).
- **1,439 nodes** (52.3%) are standardized lowercase phenotype names (`NAME:clean_string`) pending further ontology alignment.

---

## 3. ASMSG Clean Master Dataset (`datasets/asmsg_clean/`)

### **`asmsg_nodes.csv` (ncRNA Master Table)**
* **Count:** **20,973 sequence-matched human ncRNA nodes**
* **Source Databases:** miRBase v22, GENCODE v44, Ensembl 110, circBase, piRBase.
* **Columns:**
  1. `RNA Symbol`: Standardized gene symbol or accession ID.
  2. `RNA Type`: Biotype class (`miRNA`, `lncRNA`, `circRNA`, `piRNA`, `snoRNA`, `tRNA`, `snRNA`, `rRNA`).
  3. `Sequence`: Verified nucleotide sequence string (RNA format 'U', clamped to max 1,022 nt).
  4. `Seq_Length`: Sequence length in nucleotides (nt).

### **`asmsg_diseases.csv` (Disease Master Table)**
* **Count:** **2,749 Disease Ontology & Phenotype nodes**
* **Resolution Engine:** Cross-row DO ID & MeSH inheritance mapped 1,381 unique disease strings to DO IDs and 1,514 strings to MeSH IDs, raising DO ID row coverage to **211,120 association rows (89.95%)**.
* **Columns:**
  1. `Disease_ID`: Standardized ontology key (`DOID:XXXXX`, `MESH:XXXXX`, `NAME:clean_phenotype`).
  2. `Disease_Name`: Primary canonical clinical condition name.
  3. `DO_ID`: Official Disease Ontology identifier.
  4. `MeSH_ID`: Official MeSH CUI identifier.
  5. `Total_Edges`: Number of ncRNA association edges attached to this disease node.

### **`asmsg_edges.csv` (Continuous Evidence Edge Table)**
* **Count:** **132,021 unique association edges** (constructed from 234,698 literature association records in RNADisease v4.0).
* **PMID Evidence Weighting Formula:**
  Continuous edge weights \( w_{ij} \in (0, 1] \) are computed via logarithmic scaling:
  $$w_{ij} = \frac{\log(1 + E_{ij})}{\max \log(1 + E_{ij})}$$
  where \( E_{ij} \) is the count of independent literature PubMed citations supporting the edge.
* **Columns:**
  1. `RNA Symbol`: Source ncRNA node symbol.
  2. `Disease_ID`: Target Disease node ontology ID.
  3. `Evidence_Count`: Number of supporting PubMed literature citations \( E_{ij} \) (max single edge = 257 PMIDs).
  4. `PMID_List`: Semicolon-separated list of supporting PubMed IDs.
  5. `Edge_Weight`: Continuous log-scaled weight \( w_{ij} \).

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
