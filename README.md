# Adaptive Self-Supervised Multimodal Graph Learning for ncRNA–Disease Association Prediction (ASMSG)

This repository contains the official data processing, feature extraction, and PyTorch Graph Neural Network implementation for the **ASMSG** framework, developed as part of our Final Year Design Project (FYDP) at United International University (UIU).


## 🔬 Project Overview & Abstract
Non-coding RNAs (ncRNA) are crucial regulators in many biological processes and are associated with various human diseases. Predicting ncRNA–disease associations accurately is of great importance for disease diagnosis, biomarker discovery, and targeted therapeutic strategies. 

To address the limitations of existing computational approaches (like transductive bottlenecks and superficial multimodal integration), we propose a novel method named **Adaptive Self-Supervised Multimodal Graph Learning (ASMSG)**. This framework combines rich biological information sources into a unified graph learning architecture with an adaptive multimodal feature fusion mechanism. Additionally, we employ a cross-view self-supervised contrastive learning paradigm to leverage complementary information across different biological views, allowing us to learn robust and generalizable node representations.

## 📁 Repository Structure

### 1. `datasets/`
Contains the raw biological and chemical datasets used to construct the Heterogeneous Graph.
* **HMDD v4.0:** Curated ncRNA-Disease associations.
* **LncRNADisease:** Long non-coding RNA and disease associations.
* **ncRNADrug:** ncRNA-Drug interactions (curated from multiple sources like CCLE, GEO, NCI60).
* **miRBase:** Biological FASTA sequences (`mature.fa`) for RNA feature extraction.
* **DrugBank / PubChem:** SMILES strings representing chemical drug structures.

### 2. `features/`
Contains the pre-computed mathematical embeddings (features) and the finalized PyTorch graph object.
* **`rna_features_rnafm.npy`**: 640-dimensional feature matrix for 1,203 unique RNAs, extracted using the **RNA-FM** Foundation Model.
* **`drug_features_chemberta.npy`**: 384-dimensional feature matrix for 1,354 unique Drugs, extracted using **ChemBERTa-2**.
* **`asmsg_hetero_graph.pt`**: The massive PyTorch Geometric `HeteroData` graph object (7.4 MB). This graph mathematically links the embeddings together via 46,722 RNA-Disease edges and 12,730 RNA-Drug edges.

### 3. `src/`
Contains the Python source code for data processing and model architecture.
* **`src/features/`**: Scripts used to extract raw biological/chemical data into embedding vectors.
  * `extract_rna_features.py`: Parses miRBase FASTA files and runs sequences through the RNA-FM transformer.
  * `extract_drug_features.py`: Parses SMILES strings and runs them through the ChemBERTa-2 transformer.
* **`src/graph/`**: Scripts to construct and visualize the mathematical graph.
  * `build_hetero_graph.py`: The script that merges all datasets into the final `asmsg_hetero_graph.pt` file.
  * `view_graph.py` & `visualize_graph.py`: Utility scripts for inspecting the network topology.
* **`src/models/`**: The Deep Learning architecture.
  * `asmsg_model.py`: The core Heterogeneous Graph Neural Network (GNN) built using PyTorch Geometric (PyG). It features multi-modal projection layers, `HeteroConv` message passing, and an MLP link predictor.

## 🚀 Setup & Installation
1. Clone the repository.
2. Install dependencies:
```bash
pip install torch torch-geometric pandas numpy transformers scikit-learn
```
3. Load the graph directly into memory:
```python
import torch
graph = torch.load('features/asmsg_hetero_graph.pt', weights_only=False)
print(graph)
```

## 📊 Final Trained Results
During our inductive cold-start evaluation, the complete ASMSG framework (equipped with Attention-Gated Multimodal Fusion, Adaptive Edge Denoising, InfoNCE, and Margin Triplet Loss) demonstrated rapid convergence and state-of-the-art predictive accuracy. 

After training with Early Stopping (`min_delta=0.001`, `patience=25`), the best model was evaluated on a completely unseen **10% Test Set** (comprising structurally disjoint edges to overcome the transductive bottleneck). The results are as follows:

**ncRNA-Disease Association:**
* **Test AUC:** `0.9807` (98.07%)
* **Test AUPR:** `0.9812` (98.12%)
* **F1-Score:** `0.92`

**ncRNA-Drug Target Prediction:**
* **Test AUC:** `0.9918` (99.18%)
* **Test AUPR:** `0.9905` (99.05%)
* **F1-Score:** `0.97`

*Note: These metrics were generated on an isolated validation split (30% disjoint edges) to explicitly evaluate the framework's capability to overcome the transductive bottleneck in cold-start scenarios.*
