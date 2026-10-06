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

## 3. Phase 1 Rigor Overhaul & 8-Point Critique Resolution Plan

To resolve all 8 critique points with absolute scientific rigor, Phase 1 is extended to implement the following five dedicated execution modules:

### Module 1: Synonym Engine, Alias Resolution & Enriched Edge Metadata (Critiques 4 & 7)
- [ ] **Task 1.1:** Build `src/data/alias_resolution.py` using HGNC Helper tables, NCBI Gene Synonyms, and RNAcentral Accessions to map legacy/deprecated symbols (`NCRNA...`, `NONHSAT...`, `MIMAT...`) prior to sequence matching.
- [ ] **Task 1.2:** Update `build_asmsg_dataset.py` to retain `Publication_Year`, `Evidence_Type` (experimental vs predicted), `Stable_RNA_ID` (Ensembl/miRBase/HGNC), and `Dysregulation_Direction` (up/down) directly in `asmsg_edges.csv`.

### Module 2: Exact Disease Ontology Hierarchy & Funnel Verification (Critiques 5 & 6)
- [ ] **Task 2.1:** Replace loose regex string stripping (`clean_name()`) with exact Disease Ontology `doid.obo` DAG parsing and MeSH tree mappings to preserve fine-grained parent-child term distinctions.
- [ ] **Task 2.2:** Build `src/data/verify_dataset_funnel.py` to calculate all row drops, mapping rates, biotype breakdowns, and edge statistics dynamically from exported CSVs, eliminating manual textual discrepancies.

### Module 3: Leakage-Free Dynamic Topology Generator (Critique 1)
- [ ] **Task 3.1:** Remove static sequence similarity edges from `datasets/asmsg_clean/`.
- [ ] **Task 3.2:** Implement `src/graph/build_split_topology.py` to construct sequence similarity edges **dynamically strictly within training fold nodes ($\mathcal{V}_{\text{train}}$)** during PyG graph building. Enforce biotype-stratified & length-normalized TF-IDF cosine similarity.

### Module 4: Gold-Standard Experimental Evidence Benchmark (Critique 3)
- [ ] **Task 4.1:** Partition `asmsg_edges.csv` into `Experimental_Gold_Standard` vs `Predicted_or_TextMined` subsets.
- [ ] **Task 4.2:** Benchmark all baseline and GNN models on the `Experimental_Gold_Standard` subset to prevent circular training on legacy predictor outputs.

### Module 5: Leakage-Controlled Negative Sampling & Ranking Suite (Critique 2)
- [ ] **Task 5.1:** Implement **Degree-Matched Negative Sampling** and **Hard-Negative Sampling** in `src/evaluation/negative_sampling.py`.
- [ ] **Task 5.2:** Implement **Full Candidate Ranking Protocol** (Top-K Accuracy, MRR, NDCG@K, Hits@K) against all unobserved candidate diseases per RNA.

---

## 4. Present State vs. Gained State Comparison Matrix

| Critique # | Domain | Present Pipeline (Before Fix) | Gained State (After Action Plan) | Methodological Gain / Impact |
| :--- | :--- | :--- | :--- | :--- |
| **1. Topology Leakage** | Sequence Similarity Edges | Static `ncrna_sequence_similarity_edges.csv` loaded globally across train/test nodes; single 0.85 threshold. | Dynamic per-split topology built **strictly on $\mathcal{V}_{\text{train}}$**; biotype-stratified & length-normalized TF-IDF cosine similarity. | **Zero cold-start data leakage**; valid U-S and U-U evaluation. |
| **2. Negative Sampling** | Evaluation Protocol | Uniform random negative sampling from unobserved pairs; 1:1 balanced AUROC/AUPR. | **Degree-Matched & Hard Negative Sampling** + Full Candidate Ranking (MRR, NDCG@K, Hits@K against all unobserved diseases). | **Eliminates degree/popularity shortcut**; measures true biological prediction. |
| **3. Evidence Circularity** | Edge Weighting & Evidence Type | Single composite $w_{ij}$ containing mixture of experimental and computational prediction scores. | Explicit `Evidence_Type` tagging (`Experimental` vs `Predicted`); evaluation on **Experimental Gold-Standard subset**. | **Eliminates circular reasoning** (no training on legacy predictor outputs). |
| **4. Sequence Match Rate** | Alias & Synonym Resolution | Exact symbol string match against FASTA headers (8,421 symbols matched, 13.6%). | **HGNC / NCBI / RNAcentral Alias Resolution Engine** resolving deprecated symbols and transcript accessions. | Recovers long-tail ncRNAs; increases verified node coverage to estimated 12,000+ nodes. |
| **5. Count Consistency** | Documentation & Funnel | Static manual documentation subject to revision drift. | **Programmatic Funnel Verification Script (`verify_dataset_funnel.py`)** executing directly against CSVs. | **0 manual discrepancy**; 100% reproducible statistical pipeline. |
| **6. Disease Ontology** | Disease Normalization | Loose string regex stripping (`clean_name()`), risking false parent-child term collapses. | Exact **Disease Ontology `doid.obo` DAG & MeSH Tree mapping tables**. | **Preserves exact DO hierarchy**; enables hierarchy-aware cold-start splits. |
| **7. Edge Metadata** | Master Edge Table Columns | `RNA Symbol`, `Disease_ID`, `Evidence_Count`, `PMID_List`, `Confidence_Score`, `Edge_Weight`. | Adds `Publication_Year`, `Evidence_Type`, `Stable_RNA_ID`, and `Dysregulation_Direction`. | Natively enables temporal splits and experimental filtering without runtime API calls. |
| **8. Supplementary Modality** | Drug Interaction Data | Undocumented DrugBank snapshot. | **Explicit DrugBank provenance & licensing documentation**, or decoupled optional modality. | Compliance with open-science standards and DrugBank academic licensing terms. |

---

## 5. Phase-by-Phase Technical Execution Plan (Phases 2–5)

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

