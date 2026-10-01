# ASMSG Datasets

## Folder Structure

```
datasets/
├── hmdd/              # HMDD v4.0 — miRNA-disease associations
│   └── alldata_v4.xlsx    # Full miRNA-disease association data (2023.07)
│
├── ncrnadrug/         # ncRNADrug — ncRNA-drug resistance & targeting
│   ├── DR_Curated.xlsx    # Drug Resistance - Curated (low-throughput)
│   ├── DR_GEO.xlsx        # Drug Resistance - GEO (high-throughput)
│   ├── DR_NCI60.xlsx      # Drug Resistance - NCI-60 (high-throughput)
│   ├── DR_CCLE.xlsx       # Drug Resistance - CCLE (high-throughput)
│   ├── DT_Curated.xlsx    # Drug Target - Curated (low-throughput)
│   ├── DT_GEO.xlsx        # Drug Target - GEO (high-throughput)
│   └── DT_CMap.xlsx       # Drug Target - CMap (high-throughput)
│
├── miRBase/           # miRBase — miRNA sequences
│   ├── mature.fa          # All mature miRNA sequences (FASTA)
│   └── hsa.gff3           # Human miRNA genome coordinates
│
├── lncrnadisease/     # LncRNADisease v2.0 — lncRNA-disease (MANUAL)
│   └── (download manually - server unreliable)
│
├── drugbank/          # DrugBank — Drug structures (MANUAL)
│   └── (requires academic license - see below)
│
└── baseline/          # Your existing SSLGRDA baseline data
    ├── rda.csv            # ncRNA-disease association matrix
    ├── sr.csv             # ncRNA similarity matrix
    └── sd.csv             # Disease similarity matrix
```

## Sources

| Dataset | URL | Status |
|---------|-----|--------|
| HMDD v4.0 | https://www.cuilab.cn/hmdd | ✅ Downloaded |
| ncRNADrug | http://www.jianglab.cn/ncRNADrug/ | ✅ Downloaded |
| miRBase | https://mirbase.org/download/ | ✅ Downloaded |
| LncRNADisease v2.0 | https://www.rnanut.net/lncrnadisease/ | ⚠️ Server down - download manually |
| DrugBank | https://go.drugbank.com/ | ⚠️ Needs academic license signup |

## Manual Downloads Needed

### DrugBank (Required for drug SMILES)
1. Go to https://go.drugbank.com/releases/latest
2. Register for a free academic license
3. Download "All Drug Structures" (SDF or SMILES format)
4. Place in `datasets/drugbank/`

### LncRNADisease v2.0 (Required for lncRNA-disease associations)
1. Try: https://www.rnanut.net/lncrnadisease/ (server is sometimes down)
2. Alternative: https://www.cuilab.cn/lncrnadisease
3. Download the full association data
4. Place in `datasets/lncrnadisease/`
