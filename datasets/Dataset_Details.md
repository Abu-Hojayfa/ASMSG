# 📊 ASMSG Datasets: Complete Explanatory Guide

This document provides a detailed explanation of each dataset we collected, what it contains, and exactly how it will be used in your ASMSG framework.

---

## 1. HMDD v4.0 (miRNA-Disease)
**File:** `hmdd/alldata_v4.xlsx`
**Size:** 53,530 rows, 5 columns

### What is this?
HMDD (Human microRNA Disease Database) is the gold standard for manually curated, experimentally supported miRNA-disease associations. 

### What do the columns mean?
* **`code`**: Internal database ID.
* **`PMID`**: The PubMed ID of the scientific paper that proved this association.
* **`miRNA`**: The name of the microRNA (e.g., *hsa-mir-21*).
* **`disease`**: The name of the disease (e.g., *Breast Neoplasms*).
* **`description`**: A text summary of how the miRNA affects the disease (e.g., "upregulated in tumor tissue").

### How do we use it in ASMSG?
This is the core of your **Heterogeneous Graph**. We will extract all the unique `miRNA` and `disease` names to create nodes, and every row in this file will become an **Edge** connecting a miRNA node to a disease node.

---

## 2. ncRNADrug (ncRNA-Drug Associations)
**Files:** `DR_Curated.xlsx` (Drug Resistance) & `DT_CMap.xlsx` (Drug Target)
**Sizes:** ~29,000 and ~19,000 rows

### What is this?
A massive database linking non-coding RNAs (miRNAs, lncRNAs, circRNAs) to specific drugs. It shows whether an ncRNA causes resistance to a drug, or if the drug targets the ncRNA.

### What do the columns mean?
* **`ncRNA_Name`**: The specific RNA (e.g., *hsa-mir-155*).
* **`Drug_Name`**: The name of the drug (e.g., *Cisplatin*).
* **`Effect` / `Phenotype`**: What actually happens (e.g., "Resistance", "Sensitivity", or "Target").
* **`ncRNA_Type`**: Classifies if it's a miRNA, lncRNA, etc.

### How do we use it in ASMSG?
This dataset adds the **Drug nodes** to our graph. We will map the `ncRNA_Name` to the ones we found in HMDD, and add edges linking them to the `Drug_Name`. This forms the critical "ncRNA-Drug" relationship in your multi-entity graph.

---

## 3. Drug Structures (PubChem)
**File:** `drugbank/drug_smiles.csv`
**Size:** 1,594 rows (1,382 successfully found)

### What is this?
For a machine learning model to understand a drug, it needs to see its chemical structure, not just its English name. We fetched this data directly from the PubChem API.

### What do the columns mean?
* **`Drug_Name`**: The name of the drug (e.g., *Panobinostat*).
* **`SMILES`**: A standard text representation of a 3D chemical structure (e.g., `CC1=C(C2=CC=CC=C2N1)CCNCC...`).
* **`Found`**: True if PubChem successfully returned the structure.

### How do we use it in ASMSG?
We will feed every single `SMILES` string into the **ChemBERTa-2 Foundation Model** (in Phase 2). ChemBERTa-2 will convert this text string into a 384-dimensional mathematical vector (an "embedding") that captures the drug's exact chemical properties.

---

## 4. LncRNADisease v2.0
**File:** `website_alldata.tsv`
**Size:** 25,440 rows

### What is this?
A database specifically focused on Long non-coding RNAs (lncRNAs) and diseases. 

### What do the columns mean?
* **`ncRNA Symbol`**: The lncRNA name.
* **`Disease Name`**: The associated disease.
* **`Dysfunction Pattern`**: How the RNA is broken (e.g., mutated, over-expressed).

### How do we use it in ASMSG?
This simply adds more RNA and Disease nodes to your graph, expanding it beyond just miRNAs (from HMDD) to include thousands of lncRNAs, making your model's predictions much more robust.

---

## 5. Baseline Datasets (SSLGRDA)
**Files:** `rda.csv` (852x591), `sr.csv` (852x852), `sd.csv` (590x590)

### What is this?
These are standard, pre-computed matrices from a previous paper (SSLGRDA) that serve as a baseline.
* **`rda.csv`**: A binary matrix of 852 RNAs and 591 Diseases. (1 means they are linked, 0 means they aren't).
* **`sr.csv`**: An RNA-to-RNA similarity matrix (Gaussian interaction profile similarity). 
* **`sd.csv`**: A Disease-to-Disease similarity matrix (Semantic similarity).

### Why do the column names look like decimals?
These are raw mathematical matrices saved as CSVs without headers. Pandas accidentally read the first row of similarities (like `0.5714`) as the column names. 

### How do we use it in ASMSG?
We use these matrices to **evaluate and compare** your new ASMSG model against older models. Because these matrices are standardized in the industry, proving your model works on this specific `rda.csv` matrix proves your model is state-of-the-art.

---

## 6. miRBase Sequences
**File:** `mature.fa`
**Size:** 38,589 total sequences

### What is this?
A FASTA file containing the actual biological nucleotide sequences (A, C, G, U) for every known microRNA.

### What does it look like?
```text
>hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p
UGAGGUAGUAGGUUGUAUAGUU
```

### How do we use it in ASMSG?
We will feed the sequences (like `UGAGGUAG...`) into the **RNA-FM Foundation Model** (in Phase 2). RNA-FM will convert the biological sequence into a 640-dimensional mathematical vector (embedding) that captures the RNA's evolutionary and structural properties.
