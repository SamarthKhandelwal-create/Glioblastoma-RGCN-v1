# Glioblastoma ceRNA Network Analysis

This project constructs and analyzes a heterogeneous Competitive Endogenous RNA (ceRNA) network for Glioblastoma Multiforme (GBM) using TCGA data. It integrates data from ENCORI, miRTarBase, and miRNet, creates a graph representation, and trains a Relational Graph Neural Network (RGCN) to predict novel miRNA-lncRNA interactions.

## Project Structure

```
├── src/                    # Source code
│   ├── preprocessing/      # Data fetching and feature engineering
│   ├── graph/              # Graph construction
│   ├── training/           # Model definition and training
│   └── analysis/           # Evaluation, biomarkers, survival analysis
├── docs/                   # Documentation
├── results/                # Output files and artifacts
├── archive/                # Old scripts and status files (not used)
├── run_pipeline.py         # Main orchestrator script
├── verify_graph.py         # Quick graph verification utility
└── requirements.txt        # Python dependencies
```

## Installation

```bash
pip install -r requirements.txt
```

## Quick Start

Run the full pipeline:
```bash
python run_pipeline.py
```

## Usage (Step-by-Step)

### 1. Feature Engineering
```bash
python src/preprocessing/build_node_features.py
```

### 2. Graph Construction
```bash
# Optional: enable miRNA injection
$env:INJECT_MIRNA='true'  # PowerShell
# export INJECT_MIRNA=true  # Bash

python src/graph/build_graph.py
```

### 3. Model Training
```bash
# Single training run
python src/training/train_model.py

# Cross-validation
python src/training/train_cross_validation.py --graph-path results/hetero_graph_GBM.pt --folds 10
```

### 4. Analysis

```bash
# Biomarker discovery
python src/analysis/find_biomarkers.py

# Survival analysis (Kaplan-Meier + Cox regression)
python src/analysis/survival_analysis.py --offline  # Use --offline for demo mode

# Novel interaction prediction
python src/analysis/predict_novel_interactions.py

# GO/KEGG enrichment
python src/analysis/run_enrichment.py
```

### 5. Verification
```bash
python verify_graph.py
```

## Key Outputs

| File | Description |
|------|-------------|
| `results/hetero_graph_GBM.pt` | PyTorch Geometric graph object |
| `results/node_features_matrix.csv` | Z-scored node feature matrix |
| `results/node_mapping.json` | Node index ↔ ID mapping |
| `results/edge_metadata.csv` | Edge provenance and confidence |
| `results/biomarkers/` | Biomarker discovery results |
| `results/survival/` | Survival analysis (KM curves, Cox results) |
| `results/cv_results/` | Cross-validation fold results |
