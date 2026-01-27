# Graph Building with ceRNA Inference

## Summary of Changes to `src/graph/build_graph.py`

### 1. **miRNA Feature Merging**
- Looks for `results/miRNA_features.csv` and merges into node feature matrix if present
- Fallback: `INJECT_MIRNA=true` env var injects zero-valued rows for miRNAs found in interactions but missing from features
- **Purpose**: Ensures miRNAs are included in graph (previously excluded)

### 2. **ceRNA Co-targeting Edge Inference** (NEW)
- **Function**: `infer_cerna_edges(gene_to_index, gene_type_map, mirna_to_targets)`
- **Logic**:
  - Tracks which miRNAs regulate which mRNAs during interaction processing
  - Post-processing: finds lncRNA-mRNA pairs sharing miRNA regulators
  - Creates **relation type 2** edges (ceRNA sponging)
- **Biological**: Implements the paper's ceRNA topology where lncRNAs sequester miRNAs that would otherwise silence mRNAs
- **Output**: Reports inferred edge count separately

### 3. **Edge Tracking for Inference**
- During interaction loop: build `mirna_to_targets` dict mapping miRNA → set of mRNA targets
- Used by ceRNA inference to identify shared regulators

## Running the Rebuild

### Option A: With miRNA Injection (ensure all miRNAs get features)
```powershell
$env:INJECT_MIRNA='true'
python src/graph/build_graph.py
python scripts/count_graph.py
```

### Option B: With miRNA Feature CSV
1. Create `results/miRNA_features.csv` with miRNA expression/features
2. Run:
```powershell
python src/graph/build_graph.py
python scripts/count_graph.py
```

### Option C: Standard (uses existing feature matrix)
```powershell
python src/graph/build_graph.py
python scripts/count_graph.py
```

## Expected Output Changes

**Before ceRNA inference**:
- Edges: miRNA→mRNA (rel 0) + lncRNA→miRNA (rel 1)
- Example: ~1,256 edges (from baseline)

**After ceRNA inference**:
- Same edges PLUS lncRNA→mRNA co-targeting edges (rel 2)
- Expected increase: depends on miRNA target overlap

**Example console output**:
```
Processed 1200 rows. Accepted 1000 valid edges.
Inferring ceRNA edges from shared miRNA regulators...
Inferred 347 ceRNA co-targeting edges.
Saved Heterogeneous Graph to results/hetero_graph_GBM.pt
```

## New Relation Types

| Type | Edge         | Biological Meaning |
|------|---------------|--------------------|
| 0    | miRNA → mRNA  | Direct silencing |
| 1    | lncRNA → miRNA | Sponging (sequestering) |
| 2    | lncRNA → mRNA | ceRNA co-targeting (NEW) |

## Files Modified

- `src/graph/build_graph.py` — Added ceRNA inference + miRNA feature merging
- `scripts/rebuild_and_count.py` — Helper to rebuild and report stats
- `results/README_MIRNA_INJECTION.md` — Usage docs

## Next Steps

1. Run graph rebuild with option A or B above
2. Check `results/hetero_graph_GBM.pt` node/edge counts
3. Train model with expanded graph: `python src/training/train_model.py`
