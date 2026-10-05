# ASMSG Datasets (v2)

This project uses a unified dataset pipeline based on RNADisease v4.0, moving away from legacy miRNA-only databases to support true multimodal multi-biotype ncRNA-disease predictions.

## Folder Structure

```
datasets/
├── asmsg_clean/           # Master Overhauled Dataset (Clean CSVs)
│   ├── asmsg_nodes.csv    # 20,973 verified ncRNA nodes with sequence strings
│   ├── asmsg_diseases.csv # 2,749 Disease Ontology (DO ID) standardized nodes
│   └── asmsg_edges.csv    # 132,021 unique edges with continuous PMID evidence weights
│
├── rnadisease_v4/         # Primary Edge Source (RNADisease v4.0)
│   └── alldata.xlsx       # Replaces HMDD and LncRNADisease (343k entries)
│
├── ncrna_sequences/       # RNA-FM Sequence Input
│   ├── gencode_lncrna.fa.gz     # GENCODE Human lncRNAs
│   ├── ensembl_ncrna.fa.gz      # Ensembl Human ncRNAs
│   ├── circbase_circrna.fa.gz   # circBase Human circRNAs
│   └── pirbase_pirna.fa.gz      # piRBase Human piRNAs
│
├── ncrnadrug/             # Drug Target Edges
│   └── DR_Curated.xlsx    # ncRNA-drug curated resistance
│
├── drugbank/              # Drug ChemBERTa Input
│   └── drug_smiles.csv    # PubChem CID and SMILES strings
│
└── mirbase_sequences/     # miRBase sequences for miRNAs
```

## Master Output Datasets (`asmsg_clean/`)

1. **`asmsg_nodes.csv`**: Contains **20,973 verified ncRNA sequence nodes** across miRNAs, lncRNAs, circRNAs, piRNAs, snoRNAs, and tRNAs.
2. **`asmsg_diseases.csv`**: Contains **2,749 standardized Disease Ontology (DO ID) nodes**.
3. **`asmsg_edges.csv`**: Contains **132,021 unique edges** with continuous PMID evidence weights \( w_{ij} \in (0, 1] \) derived from 234,698 literature records.
