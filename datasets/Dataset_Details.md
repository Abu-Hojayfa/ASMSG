# 📊 ASMSG Datasets: Comprehensive Guide & Benchmark Analysis

This document provides the definitive specification of the **ASMSG Clean** master dataset, its empirical metrics, data processing pipeline, and comparative novelty against the 5 key baseline papers in the field.

---

## 1. Comparative Novelty Matrix vs. 5 Benchmark Papers

To establish state-of-the-art (SOTA) publication novelty for peer review (*Bioinformatics*, *Nucleic Acids Research*, *IEEE TPAMI*), ASMSG Clean was designed to resolve the fundamental structural flaws of existing benchmark datasets:

| Feature / Metric | **SSCLMD** (2023) | **SSLGRDA** (2024) | **GSLRDA** (2024) | **MIFNDRA** (2023) | **DMGAT** (2024) | **ASMSG Clean (Ours)** |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Primary Journal** | *IEEE/ACM TCBB* | *Brief. Bioinform.* | *Comput. Biol. Med.* | *Comput. Biol. Med.* | *IEEE JBHI* | **This Work** |
| **Biotype Scope** | miRNA only | miRNA + lncRNA | miRNA + lncRNA | miRNA only | miRNA only | **All 6 ncRNA Biotypes (miRNA, lncRNA, circRNA, piRNA, snoRNA, tRNA)** |
| **Verified ncRNA Nodes** | 853 | 1,002 | 852 | 788 | 853 | **20,973 verified sequence nodes** |
| **Disease Scope** | 591 raw strings | 590 raw strings | 591 raw strings | 441 raw strings | 591 raw strings | **2,749 DO ID / MeSH standardized terms** |
| **Unique Edge Count** | 5,424 binary (0/1) | 9,122 binary (0/1) | 5,424 binary (0/1) | 5,149 binary (0/1) | 5,424 binary (0/1) | **132,021 unique edges (from 234,698 literature records)** |
| **Edge Weight Scheme** | Unweighted binary | Unweighted binary | Unweighted binary | Unweighted binary | Unweighted binary | **Continuous PMID evidence weights \( w_{ij} \in (0, 1] \)** |
| **Sequence Encoding** | k-mer (3-mer) | None | None | Sequence length | None | **RNA-FM Foundation Model (640d contextual embeddings)** |
| **Disease Encoding** | MeSH tree | None | Disease semantic | MeSH semantic | None | **BioBERT Language Model (768d DO ID semantic vectors)** |
| **Evaluation Mode** | Transductive 5-fold CV | Transductive 5-fold CV | Transductive 5-fold CV | Transductive 5-fold CV | Transductive 5-fold CV | **Disjoint Inductive Cold-Start (Unseen ncRNAs & Diseases)** |

---

## 2. ASMSG Clean Master Dataset (`datasets/asmsg_clean/`)

### **`asmsg_nodes.csv` (ncRNA Master Table)**
* **Count:** **20,973 verified human ncRNA nodes**
* **Source Databases:** miRBase v22, GENCODE v44, Ensembl 110, circBase, piRBase.
* **Columns:**
  1. `RNA Symbol`: Standardized gene symbol or accession ID.
  2. `RNA Type`: Biotype class (`miRNA`, `lncRNA`, `circRNA`, `piRNA`, `snoRNA`, `tRNA`, `snRNA`, `rRNA`).
  3. `Sequence`: Verified nucleotide sequence string (converted to RNA format 'U' and clamped to max 1,022 nt).
  4. `Seq_Length`: Sequence length in nucleotides (nt).

### **`asmsg_diseases.csv` (Disease Master Table)**
* **Count:** **2,749 standardized Disease Ontology nodes**
* **Source:** Disease Ontology (DO ID) and Medical Subject Headings (MeSH).
* **Columns:**
  1. `Disease_ID`: Standardized ontology key (e.g., `DOID:1612`, `MESH:D001943`, `NAME:breast carcinoma`).
  2. `Disease_Name`: Primary clinical condition name.
  3. `DO_ID`: Official Disease Ontology identifier.
  4. `MeSH_ID`: Official MeSH CUI identifier.
  5. `Total_Edges`: Number of verified ncRNA association edges attached to this disease node.

### **`asmsg_edges.csv` (Continuous Evidence Edge Table)**
* **Count:** **132,021 unique association edges** (constructed from 234,698 literature association records in RNADisease v4.0).
* **PMID Evidence Weighting Formula:**
  To preserve quantitative biological confidence instead of binarizing edges into 0/1, continuous edge weights \( w_{ij} \in (0, 1] \) are computed via logarithmic scaling:
  $$w_{ij} = \frac{\log(1 + E_{ij})}{\max \log(1 + E_{ij})}$$
  where \( E_{ij} \) is the count of independent literature PubMed citations supporting the edge.
* **Columns:**
  1. `RNA Symbol`: Source ncRNA node symbol.
  2. `Disease_ID`: Target Disease node ontology ID.
  3. `Evidence_Count`: Number of supporting PubMed literature citations \( E_{ij} \) (max single edge = 257 PMIDs).
  4. `PMID_List`: Semicolon-separated list of supporting PubMed IDs.
  5. `Edge_Weight`: Continuous log-scaled weight \( w_{ij} \).

---

## 3. Methodological Rationale & Reviewer Defense

### **A. Why 20,973 Nodes Out of 61,947 Candidates?**
Out of 61,947 candidate RNA symbols in RNADisease v4, multi-key sequence indexing successfully retrieved verified biological sequences for 20,973 ncRNAs (**33.86% matching rate**, retaining **81.8% of total association records**).  
The remaining unmapped candidate symbols are unannotated high-throughput transcript IDs or deprecated aliases that lack entry in gold-standard genomic reference databases. Retaining nodes without sequence data would force zero-padding or random noise initialization, corrupting Transformer embedding space.

### **B. PMID Evidence-Weighted Negative Sampling**
RNADisease v4 contains exclusively validated positive associations. To prevent false-negative bias during model training, ASMSG employs **evidence-weighted negative sampling**, selecting unobserved (ncRNA, Disease) pairs primarily from entity combinations with zero PubMed literature co-mentions and low topological similarity.

---

## 4. Supplementary Multi-Entity Datasets

### **`datasets/ncrnadrug/DR_Curated.xlsx` (ncRNA-Drug Associations)**
Curated ncRNA-drug resistance and target interactions used to introduce **Drug nodes** into the heterogeneous graph topology.

### **`datasets/drugbank/drug_smiles.csv` (PubChem Drug SMILES)**
Canonical SMILES chemical structure strings for drug nodes, processed via **ChemBERTa-2** (`DeepChem/ChemBERTa-77M-MTR`) to produce 384-dimensional chemical feature representations.

---

## 5. Benchmarking Strategy: Algorithmic Architecture Baselines

Notice: Static legacy binary matrices (e.g., legacy `rda.csv`) have been **purged** from the repository to prevent data contamination. 

In modern computational biology literature, **benchmarking is conducted by re-implementing baseline model architectures** (SSCLMD, SSLGRDA, GSLRDA, MIFNDRA, DMGAT) and evaluating them under identical train/test split conditions on standard gold-standard datasets:
1. **Primary Evaluation Benchmark:** ASMSG Clean (20,973 ncRNAs × 2,749 DO IDs).
2. **Secondary Gold-Standard Benchmark:** HMDD v3.2 / v4.0 (35,547 associations across 1,206 miRNAs and 893 diseases).
