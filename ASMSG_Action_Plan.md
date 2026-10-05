# 🧬 ASMSG: Oxford-Rigor Master Action & Implementation Plan
**Adaptive Self-Supervised Multimodal Graph Learning for ncRNA–Disease Association Prediction**

---

## 1. Executive Feasibility, Novelty & Publication Audit

Before executing each phase, every component is subjected to a rigorous 3-question peer-review evaluation:

| Phase / Module | **Is it Novel?** | **Is it Technically Possible?** | **Is it Worth It for Publication?** | **Primary Risk & Engineering Mitigation** |
| :--- | :--- | :--- | :--- | :--- |
| **Phase 1: Multi-Biotype Data Engine** | **HIGHLY NOVEL** (6 biotypes, 20,973 nodes, continuous PMID weights) | **YES (COMPLETED)** | **ESSENTIAL** | *Risk:* 66% candidate sequence gap. *Fix:* Documented genomic database coverage defense in `Dataset_Details.md`. |
| **Phase 2: Multimodal Foundation Embeddings** | **HIGHLY NOVEL** (RNA-FM 640d + BioBERT 768d + ChemBERTa-2 384d) | **YES** | **ESSENTIAL** | *Risk:* PyTorch CUDA OOM during RNA-FM batching. *Fix:* Mini-batching (`batch_size=16`) with `fp16` mixed precision and disk caching (`.pt`). |
| **Phase 3: PyG HeteroData Graph Builder** | **MODERATE** (Multi-entity topology: ncRNA, Disease, Drug) | **YES** | **HIGH** | *Risk:* Memory footprint during GNN message passing. *Fix:* Node features consume ~63MB RAM; fits easily on GPU. |
| **Phase 4: Adaptive Dual-View GNN & InfoNCE** | **EXTREMELY NOVEL** (Learned edge-denoising gating vs random GraphCL) | **YES** | **CORE THEORETICAL CONTRIBUTION** | *Risk:* Representation collapse in InfoNCE. *Fix:* Joint objective: BCE Link Prediction + InfoNCE + Margin Triplet Loss. |
| **Phase 5: Inductive Cold-Start Evaluation** | **HIGHLY NOVEL** (Disjoint ncRNA & Disease cold-start splits + PMID Negative Sampling) | **YES** | **ESSENTIAL** | *Risk:* Data leakage across PyG graph layers. *Fix:* Test nodes completely detached from training `edge_index` matrix. |

---

## 2. Technical Architecture & Mathematical Specification

```
                     ┌──────────────────────────────────────────────┐
                     │          MULTIMODAL FEATURE ENGINE           │
                     └──────────────────────┬───────────────────────┘
                                            │
         ┌──────────────────────────────────┼──────────────────────────────────┐
         │                                  │                                  │
         ▼                                  ▼                                  ▼
 ┌───────────────┐                  ┌───────────────┐                  ┌───────────────┐
 │    RNA-FM     │                  │    BioBERT    │                  │  ChemBERTa-2  │
 │ (640d Vector) │                  │ (768d Vector) │                  │ (384d Vector) │
 └───────┬───────┘                  └───────┬───────┘                  └───────┬───────┘
         │                                  │                                  │
         ▼                                  ▼                                  ▼
    [ ncRNA Node ]                   [ Disease Node ]                    [ Drug Node ]
         │                                  │                                  │
         └──────────────────────────────────┼──────────────────────────────────┘
                                            │
                                            ▼
                     ┌──────────────────────────────────────────────┐
                     │    HETEROGENEOUS GRAPH (PyG HeteroData)      │
                     │  - Continuous Literature PMID Edge Weights   │
                     └──────────────────────┬───────────────────────┘
                                            │
                                            ▼
                     ┌──────────────────────────────────────────────┐
                     │         DUAL-VIEW ADAPTIVE GNN ENCODER       │
                     │  - Original View & Denoised Adaptive View    │
                     └──────────────────────┬───────────────────────┘
                                            │
                                            ▼
                     ┌──────────────────────────────────────────────┐
                     │    SELF-SUPERVISED CONTRASTIVE LOSS (InfoNCE)│
                     │    + INDUCTIVE LINK PREDICTION HEAD           │
                     └──────────────────────────────────────────────┘
```

### A. Adaptive Edge-Denoising Gating MLP
Instead of random edge dropouts (GraphCL/SSCLMD), ASMSG computes a learned gating probability \( p_{ij} \) for each edge \((i, j)\):
$$p_{ij} = \sigma \left( \mathbf{W}_g [\mathbf{h}_i \,||\, \mathbf{h}_j] + b_g \right) \cdot w_{ij}$$
where \( w_{ij} \) is the PubMed evidence weight. Edges with low gating probability are masked in the secondary contrastive view, forcing the GNN to learn robust representations.

### B. Joint Multi-Task Loss Objective
The overall objective function combines link prediction supervised loss, self-supervised InfoNCE contrastive loss, and margin triplet loss:
$$\mathcal{L}_{\text{total}} = \mathcal{L}_{\text{BCE}} + \lambda_1 \mathcal{L}_{\text{InfoNCE}} + \lambda_2 \mathcal{L}_{\text{Triplet}}$$

---

## 3. Phase-by-Phase Technical Execution Plan

### Phase 1: Data Engine Repair & Mapping Overhaul (✅ COMPLETED)
- [x] **Task 1.1:** Overhaul `match_ncrna_sequences.py` with multi-key sequence indexing. Retrieved **20,973 verified ncRNA sequence nodes**.
- [x] **Task 1.2:** Implement DO ID standardization in `build_asmsg_dataset.py`. Standardized 3,321 raw strings into **2,749 Disease Ontology nodes**.
- [x] **Task 1.3:** Calculate continuous PMID evidence edge weights \( w_{ij} \in (0, 1] \) across **132,021 unique edges**. Exported to `datasets/asmsg_clean/`.

### Phase 2: Multimodal Embedding Extraction Engine (⏳ NEXT UP)
- [ ] **Task 2.1:** Execute `src/features/run_rna_fm.py` using `d:\fydp\venv\Scripts\python.exe` to generate 640-dim embeddings for 20,973 ncRNAs using chunked FP16 batching (`batch_size=16`). Save to `datasets/asmsg_clean/rna_fm_embeddings.pt`.
- [ ] **Task 2.2:** Build `src/features/extract_disease_biobert.py` using BioBERT (`dmis-lab/biobert-base-cased-v1.2`) to produce 768-dim disease embeddings for 2,749 DO IDs. Save to `datasets/asmsg_clean/disease_biobert_embeddings.pt`.
- [ ] **Task 2.3:** Build `src/features/extract_drug_chemberta.py` using ChemBERTa-2 (`DeepChem/ChemBERTa-77M-MTR`) to produce 384-dim drug embeddings. Save to `datasets/asmsg_clean/drug_chemberta_embeddings.pt`.

### Phase 3: PyG HeteroData Graph Builder
- [ ] **Task 3.1:** Create `src/graph/build_hetero_graph.py` to construct `torch_geometric.data.HeteroData` containing all node features, edge indices, and continuous edge weights. Save to `datasets/asmsg_clean/asmsg_hetero_graph.pt`.

### Phase 4: Core Model Architecture & Loss Functions
- [ ] **Task 4.1:** Implement `src/models/asmsg_gnn.py` with dual-view encoder, adaptive gating MLP, and heterogeneous message passing (`HeteroConv` with `GATv2Conv`).
- [ ] **Task 4.2:** Implement InfoNCE contrastive loss, margin triplet loss, and BCE link prediction head in `src/models/losses.py`.

### Phase 5: Inductive Cold-Start Benchmarking Suite
- [ ] **Task 5.1:** Construct `src/data/negative_sampling.py` for PMID-evidence-weighted negative sampling.
- [ ] **Task 5.2:** Create `src/evaluation/cold_start_split.py` for 3 disjoint evaluation regimes:
  1. Standard Transductive 5-Fold CV.
  2. Inductive Cold-Start ncRNA (20% held-out ncRNAs completely detached from training graph).
  3. Inductive Cold-Start Disease (20% held-out Diseases completely detached from training graph).
- [ ] **Task 5.3:** Create `src/evaluation/benchmark_baselines.py` to evaluate model baselines (SSCLMD, SSLGRDA, GSLRDA, MIFNDRA, DMGAT, GCN, GATv2, RGCN) on both ASMSG Clean and HMDD v3.2 datasets.

---

## 4. Multi-Agent Track Plan

```
                  ┌─────────────────────────────────────────┐
                  │          DATA & MAPPING AGENT           │
                  │  - Phase 1 Overhaul (✅ COMPLETED)       │
                  │  - Negative Sampling Strategy           │
                  └────────────────────┬────────────────────┘
                                       │
                                       ▼
                  ┌─────────────────────────────────────────┐
                  │        FOUNDATION FEATURE AGENT         │
                  │  - Runs RNA-FM FP16 extraction          │
                  │  - Runs BioBERT disease extraction      │
                  │  - Runs ChemBERTa-2 drug extraction     │
                  └────────────────────┬────────────────────┘
                                       │
                                       ▼
                  ┌─────────────────────────────────────────┐
                  │       GRAPH & MODEL DESIGN AGENT        │
                  │  - Constructs PyG HeteroData graph      │
                  │  - Implements Adaptive Edge-Denoising   │
                  │  - Implements Hetero GNN + InfoNCE      │
                  └────────────────────┬────────────────────┘
                                       │
                                       ▼
                  ┌─────────────────────────────────────────┐
                  │          BENCHMARKING AGENT             │
                  │  - Evaluates 3 split regimes            │
                  │  - Benchmarks vs SSCLMD, SSLGRDA, etc.  │
                  │  - Runs HMDD v3.2 external validation   │
                  └─────────────────────────────────────────┘
```

---

## 5. Rigorous Verification Protocol

Before declaring completion of any phase:
1. **Data Integrity:** Zero NaN/Inf values, node IDs must be strictly continuous integers `0..N-1`.
2. **Leakage Prevention:** Inductive cold-start test nodes must have **zero edges** in the training graph `edge_index` matrix.
3. **Execution Command:** All scripts MUST be run with `d:\fydp\venv\Scripts\python.exe` and verified via explicit CLI output logs.
