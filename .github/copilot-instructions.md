# Copilot / AI Agent Instructions — Glioblastoma-RGCN-v1

Brief: build a reproducible pipeline from expression data → heterogeneous ceRNA graph → RGCN link-predictor. Focus on data artifacts, ID normalization, and reproducible command sequences.

Key pipeline stages
- Preprocessing / features: [src/preprocessing/build_node_features.py](src/preprocessing/build_node_features.py)
- Graph construction: [src/graph/build_graph.py](src/graph/build_graph.py)
- Model & training: [src/model.py](src/model.py), [src/training/train_model.py](src/training/train_model.py)
- Analysis: [src/analysis/](src/analysis/) — biomarkers, survival, enrichment, visualization
- Orchestration: `run_pipeline.py` (end-to-end) and `src/training/train_cross_validation.py` (CV)

Important artifacts (results/)
- `results/node_features_matrix.csv`: indexed by normalized IDs; columns = numeric features (z-scored).
- `results/interactions.csv`: raw interaction list used to build edges.
- `results/hetero_graph_GBM.pt`: PyG `Data` object consumed by training scripts.
- `results/node_mapping.json`: node index ↔ original ID mapping.

Project-specific conventions (must follow)
- ID normalization: always use `normalize_id()` logic in [src/graph/build_graph.py](src/graph/build_graph.py) — strip Ensembl versions, lowercase miRNA names.
- Results directory: controlled by `RESULTS_DIR` env var (default `results`); prefer it when reading/writing.
- DE files naming: `DE_*_logFC_abs_gt1.tsv` are expected by graph builder; missing files may skip nodes.
- miRNA handling: miRNA injection is opt-in via `INJECT_MIRNA=true` and `results/miRNA_features.csv`; otherwise absent miRNAs will be excluded (zero edges possible).
- Feature scaling: features are z-scored before saving and models expect numeric `data.x` with shape `[N, F]`.

How training expects the graph
- `train_model.py` expects a PyG `Data` with `x` (node features), `edge_index` ([2, E]), and `edge_attr` (relation ids). Relation mapping: e.g., 0 = miRNA→mRNA, 1 = lncRNA→miRNA — verify in `src/graph/build_graph.py` and `results/edge_metadata.csv`.

Quick reproducible commands
```powershell
pip install -r requirements.txt
python src/preprocessing/build_node_features.py
$env:INJECT_MIRNA='true'  # optional
python src/graph/build_graph.py
python src/training/train_model.py
python src/training/train_cross_validation.py --graph-path results/hetero_graph_GBM.pt --outdir cv_results --folds 10 --epochs 100 --seed 42
```

Common checks & debugging hints
- If `hetero_graph_GBM.pt` has zero edges: inspect `results/interactions.csv` and `results/node_features_matrix.csv` (ID mismatch or missing miRNAs).
- Always normalize IDs when matching: inconsistent Ensembl versions are the main source of silent drops.
- Use `verify_graph.py` to get quick node/edge counts and provenance summary.
- Network fetches (Xena) are required by preprocessing; run from a machine with internet access or use cached files in `results/`.

Where to look in the code
- ID & matching logic: [src/graph/build_graph.py](src/graph/build_graph.py)
- Feature creation and DE naming: [src/preprocessing/build_node_features.py](src/preprocessing/build_node_features.py)
- Training & data expectations: [src/training/train_model.py](src/training/train_model.py)

If anything is unclear or you want examples (minimal Data shapes, a tiny runner, or extra checks), tell me which area to expand.
