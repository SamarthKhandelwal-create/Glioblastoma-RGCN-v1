# Glioblastoma-RGCN-v1

Pipeline: expression data → ceRNA graph → RGCN link prediction

## Key Files
- `src/preprocessing/build_node_features.py` - Feature engineering from TCGA
- `src/graph/build_graph.py` - Graph construction
- `src/model.py` - RGCN model definition
- `src/training/train_model.py` - Single training run
- `src/training/train_cross_validation.py` - 10-fold CV

## Artifacts
- `results/node_features_matrix.csv` - Node features (z-scored)
- `results/interactions.csv` - Edge list
- `results/hetero_graph_GBM.pt` - PyG Data object
- `results/node_mapping.json` - Node ID mapping

## Conventions
- IDs normalized via `normalize_id()`: strip Ensembl versions, lowercase miRNAs
- Edge types: 0=miRNA→mRNA, 1=lncRNA→miRNA, 2=ceRNA
- Set `INJECT_MIRNA=true` to include miRNAs missing from feature matrix

## Commands
```bash
python src/preprocessing/build_node_features.py
python src/graph/build_graph.py
python src/training/train_cross_validation.py --graph-path results/hetero_graph_GBM.pt --outdir results/cv_results --folds 10 --epochs 100
```
