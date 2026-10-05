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

### Phase 1: Data Engine Repair & Mapping Overhaul (✅ COMPLETED & REVISED)
- [x] **Task 1.1:** Overhaul `match_ncrna_sequences.py` with multi-key sequence indexing. Retrieved **20,973 verified ncRNA sequence nodes**.
- [x] **Task 1.2:** Implement DO ID standardization in `build_asmsg_dataset.py`. Standardized 3,321 raw strings into **2,749 Disease Ontology nodes**.
- [x] **Task 1.3:** Calculate continuous PMID evidence edge weights \( w_{ij} \in (0, 1] \) across **132,021 unique edges**. Exported to `datasets/asmsg_clean/`.
- [ ] **Task 1.4 (NEW):** Create `datasets/asmsg_clean/CURATION_POLICY.md` documenting evidence types, Homo sapiens restriction, dedup rules, conflict resolution, and precise definition of sequence matching.
- [ ] **Task 1.5 (NEW):** Add `score` column from RNADisease v4 to `asmsg_edges.csv` as an independent confidence signal distinct from PMID count.
- [ ] **Task 1.6 (NEW):** Report per-biotype statistics transparently, explicitly acknowledging miRNA dominance (74.1%).

### Phase 2: Multimodal Embedding Extraction & Feature Ablation Suite (⏳ NEXT UP)
- [ ] **Task 2.1:** Execute `src/features/run_rna_fm.py` using `d:\fydp\venv\Scripts\python.exe` to generate 640-dim embeddings for 20,973 ncRNAs using chunked FP16 batching (`batch_size=16`). Save to `datasets/asmsg_clean/rna_fm_embeddings.pt`.
- [ ] **Task 2.2:** Build `src/features/extract_disease_biobert.py` using BioBERT (`dmis-lab/biobert-base-cased-v1.2`) to produce 768-dim disease embeddings for 2,749 DO IDs. Save to `datasets/asmsg_clean/disease_biobert_embeddings.pt`.
- [ ] **Task 2.2b (NEW):** Extract **SapBERT** (`cambridgeltl/SapBERT-from-PubMedBERT-fulltext`) 768-dim disease embeddings for ablation comparison. Save to `datasets/asmsg_clean/disease_sapbert_embeddings.pt`.
- [ ] **Task 2.3:** Build `src/features/extract_drug_chemberta.py` using ChemBERTa-2 (`DeepChem/ChemBERTa-77M-MTR`) to produce 384-dim drug embeddings. Save to `datasets/asmsg_clean/drug_chemberta_embeddings.pt`.
- [ ] **Task 2.4 (NEW):** Build ablation feature sets:
  - 3-mer frequency vectors for all ncRNAs.
  - One-hot biotype encoding vectors.
  - Sequence-length-only scalar features.
  - Random init baseline embeddings (learn from scratch).

### Phase 3: PyG HeteroData Graph Builder
- [ ] **Task 3.1:** Create `src/graph/build_hetero_graph.py` to construct `torch_geometric.data.HeteroData` containing all node features (RNA-FM, BioBERT, ChemBERTa-2, ablation sets), edge indices, continuous PMID weights, and confidence scores. Save to `datasets/asmsg_clean/asmsg_hetero_graph.pt`.

### Phase 4: Core Model Architecture & Loss Functions
- [ ] **Task 4.1:** Implement `src/models/asmsg_gnn.py` with dual-view encoder, adaptive gating MLP, and heterogeneous message passing (`HeteroConv` with `GATv2Conv`).
- [ ] **Task 4.2:** Implement InfoNCE contrastive loss, margin triplet loss, and BCE link prediction head in `src/models/losses.py`.

### Phase 5: Leakage-Controlled Inductive Benchmarking & Evaluation Suite
- [ ] **Task 5.1 (REVISED):** Implement **identity-clustered splits** for ncRNA cold-start in `src/evaluation/cold_start_split.py`:
  - Cluster miRNAs by seed-region family (miRBase family annotations).
  - Cluster lncRNAs/circRNAs by sequence identity (CD-HIT at 80% threshold).
  - Hold out entire sequence/family clusters, preventing sequence leakage.
- [ ] **Task 5.1b (NEW):** Implement **hierarchy-aware splits** for disease cold-start:
  - Hold out entire disease subtrees in the DO/MeSH hierarchy, ensuring no ancestor/descendant leaks into training graph.
- [ ] **Task 5.2 (REVISED):** Implement all 4 evaluation settings:
  1. **S-S (Seen-Seen):** Standard transductive link prediction.
  2. **U-S (Unseen ncRNA, Seen Disease):** Inductive ncRNA cold-start.
  3. **S-U (Seen ncRNA, Unseen Disease):** Inductive disease cold-start.
  4. **U-U (Unseen-Unseen):** Dual cold-start (both entities unseen).
- [ ] **Task 5.2b (NEW):** Report ranking metrics alongside classification metrics:
  - Hits@10, Hits@50, MRR, AUROC, AUPR per disease and per ncRNA.
- [ ] **Task 5.3 (NEW):** Implement 5 trivial baseline predictors:
  1. Degree-product predictor (ncRNA degree × Disease degree).
  2. RNA Type biotype-only predictor.
  3. Sequence-length-only predictor.
  4. k-NN on 3-mer frequency vectors.
  5. Uniform random baseline.
- [ ] **Task 5.4 (REVISED):** Reimplement SSCLMD, SSLGRDA, GSLRDA, MIFNDRA, DMGAT architectures:
  - Train and evaluate all 5 baselines under identical ASMSG Clean splits.
  - Document adaptations for cold-start evaluation.
- [ ] **Task 5.5 (NEW):** Implement 5-seed protocol with confidence intervals:
  - Run all experiments across 5 random seeds (report mean ± std).
  - Paired t-test / Wilcoxon signed-rank test for statistical significance.
- [ ] **Task 5.6 (NEW):** Implement temporal validation split:
  - Train on edges with PMIDs published before year T; test on edges published in year T or later.

---

## 4. Multi-Agent Track Plan

```
                  ┌─────────────────────────────────────────┐
                  │          DATA & CURATION AGENT          │
                  │  - Phase 1 Overhaul (✅ COMPLETED)       │
                  │  - Task 1.4 CURATION_POLICY.md          │
                  │  - Task 1.5 Add Score Column            │
                  └────────────────────┬────────────────────┘
                                       │
                                       ▼
                  ┌─────────────────────────────────────────┐
                  │   FOUNDATION FEATURE & ABLATION AGENT   │
                  │  - RNA-FM FP16 extraction               │
                  │  - BioBERT + SapBERT extraction         │
                  │  - ChemBERTa-2 drug extraction          │
                  │  - 3-mer / One-hot / Length ablations   │
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
                  │     LEAKAGE-CONTROLLED BENCHMARK AGENT  │
                  │  - Identity & Hierarchy-aware splits    │
                  │  - Evaluates 4 quadrants (S-S, U-S, etc)│
                  │  - 5 Trivial baselines + 5 SOTA re-impls│
                  │  - 5-seed statistical & temporal tests  │
                  └─────────────────────────────────────────┘
```

---

## 5. Rigorous Verification Protocol

Before declaring completion of any phase:
1. **Data Integrity:** Zero NaN/Inf values, node IDs must be strictly continuous integers `0..N-1`.
2. **Leakage Control:** Test nodes in U-S, S-U, and U-U splits must have **zero edges** in the training graph `edge_index` matrix and zero sequence/ontology cluster overlap.
3. **Execution Command:** All scripts MUST be run with `d:\fydp\venv\Scripts\python.exe` and verified via explicit CLI output logs.
4. **No Claim Without Test:** Never declare a model or split fixed without running verification execution and checking metric outputs.

