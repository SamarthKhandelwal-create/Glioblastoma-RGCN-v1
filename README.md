# Glioblastoma ceRNA Network Analysis

RGCN-based link prediction for GBM ceRNA networks using TCGA data.

## Structure

```
src/
  preprocessing/   # Data fetching, feature engineering
  graph/           # Graph construction  
  training/        # Model training, cross-validation
  analysis/        # Biomarker discovery, visualization
results/           # Output files, trained models
```

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# 1. Build node features from TCGA data
python src/preprocessing/build_node_features.py

# 2. Construct heterogeneous graph
python src/graph/build_graph.py

# 3. Train model
python src/training/train_model.py

# 4. Run 10-fold cross-validation
python src/training/train_cross_validation.py --graph-path results/hetero_graph_GBM.pt --outdir results/cv_results --folds 10 --epochs 100

# 5. Analyze results
python src/analysis/find_biomarkers.py
python src/analysis/visualize_results.py
```

## Graph Structure

- **Nodes**: mRNA, lncRNA, miRNA (10-dimensional feature vectors)
- **Edges**: 
  - Type 0: miRNA → mRNA (silencing)
  - Type 1: lncRNA → miRNA (sponging)
  - Type 2: ceRNA co-regulation (inferred)

## Data Sources

- TCGA-GBM expression and clinical data (via Xena)
- TargetScan miRNA target predictions
- ENCORI validated interactions

