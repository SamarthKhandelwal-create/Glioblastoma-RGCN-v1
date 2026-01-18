Glioblastoma ceRNA Network Analysis

This project constructs and analyzes a heterogeneous Competitive Endogenous RNA (ceRNA) network for Glioblastoma Multiforme (GBM) using TCGA data. It integrates data from ENCORI, miRTarBase, and miRNet, creates a graph representation, and trains a Relational Graph Neural Network (RGCN) to predict novel miRNA-lncRNA interactions.

 Project Structure

- `src/`: Source code.
    - `preprocessing/`: Data fetching and feature engineering.
    - `graph/`: Graph construction.
    - `training/`: Model definition and training.
    - `analysis/`: Evaluation.
- `data/`: Input data files.
- `results/`: output files.

 Installation

1. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
 Usage Pipeline

1. Feature Engineering:
   ```bash
   python src/preprocessing/build_node_features.py
   ```

2. Graph Construction:
   ```bash
   python src/graph/build_graph.py
   ```

3. Model Training:
   ```bash
   python src/training/train_model.py
   ```

4. Evaluation:
   ```bash
   python src/training/evaluate_model.py
   ```

5. Inference (Novel Discovery):
   ```bash
   python src/analysis/predict_novel_interactions.py
   ```
