# Enhanced Graph Building with TargetScan, ID Mapping, and Provenance

## Summary of Enhancements

### 1. **TargetScan Prediction Integration**
- **File**: `src/preprocessing/integrate_targetscan.py`
- **Functionality**:
  - Loads TargetScan miRNA target predictions from CSV
  - Filters by prediction score (default cutoff: 0.5)
  - Adds high-confidence predictions (score ≥ 0.7) to graph as supplementary edges
  - Marks edges with 'TargetScan' provenance source
- **Output**: `results/targetscan_predictions.csv`

### 2. **Enhanced ID Mapping (mygene + HGNC)**
- **File**: `src/utils/id_mapping.py`
- **Functionality**:
  - Multi-source ID resolution: HGNC symbols ↔ Ensembl IDs ↔ Entrez IDs
  - Fuzzy string matching for partial/alias matches
  - Fallback to common gene hardcodes if HGNC file unavailable
  - miRNA alias resolution
- **Class**: `IDMapper`
- **Usage**:
  ```python
  from src.utils.id_mapping import get_mapper
  mapper = get_mapper()
  ensembl = mapper.resolve_gene_id('TP53', target_type='ensembl')
  ```

### 3. **Edge Provenance & Confidence Annotation**
- **Tracking**: Each edge now has:
  - **Source**: 'curated' (from interactions.csv), 'TargetScan' (predictions), or 'curated+TargetScan' (both)
  - **Confidence**: 1.0 for curated, prediction score (0-1) for TargetScan
- **Output**: `results/edge_metadata.csv`
  - Columns: edge_idx, source_node, target_node, relation_type, provenance, confidence
- **Enables**: Model can weight/filter edges during training based on confidence

## Workflow

### Step 1: Create/Integrate TargetScan Predictions
```powershell
# Option A: Create mock TargetScan file (for testing)
python src/preprocessing/integrate_targetscan.py --create-mock

# Option B: Load real TargetScan data
python src/preprocessing/integrate_targetscan.py --input your_targetscan.csv --score-cutoff 0.5
```

### Step 2: Rebuild Graph (with all enhancements)
```powershell
$env:INJECT_MIRNA='true'
python src/graph/build_graph.py
```

### Step 3: Verify Edge Metadata
```powershell
python -c "
import pandas as pd
metadata = pd.read_csv('results/edge_metadata.csv')
print(metadata.head())
print('\nProvenance distribution:')
print(metadata['provenance'].value_counts())
print('\nConfidence statistics:')
print(metadata['confidence'].describe())
"
```

## Example Output

```
Loaded 156 TargetScan predictions from results/targetscan_predictions.csv
Processed 1200 rows. Accepted 1000 valid edges.
Added 89 TargetScan-only edges with high confidence (score >= 0.7).

Edge Provenance Summary:
  curated: 1000 edges
  curated+TargetScan: 45 edges
  TargetScan: 89 edges
```

## Using Edge Metadata in Training

The `edge_metadata.csv` can be loaded during model training to:
- **Weight edges**: High-confidence edges get higher loss weights
- **Filter edges**: Train only on curated interactions or only on predictions
- **Analyze**: Determine which sources contribute most to model performance

Example:
```python
metadata = pd.read_csv('results/edge_metadata.csv')

# Weight edges by confidence
edge_weights = torch.tensor(metadata['confidence'].values, dtype=torch.float)

# Or filter to only high-confidence edges
high_conf = metadata[metadata['confidence'] >= 0.8]
```

## Files Modified/Created

| File | Type | Purpose |
|------|------|---------|
| `src/preprocessing/integrate_targetscan.py` | NEW | Load TargetScan predictions |
| `src/utils/id_mapping.py` | NEW | HGNC/mygene-style ID mapping |
| `src/utils/__init__.py` | NEW | Utils package |
| `src/graph/build_graph.py` | MODIFIED | Integrate TargetScan + provenance tracking |
| `results/targetscan_predictions.csv` | OUTPUT | TargetScan predictions (CSV) |
| `results/edge_metadata.csv` | OUTPUT | Edge provenance/confidence (CSV) |

## Future Enhancements

- [ ] Integrate miRanda predictions (similar to TargetScan)
- [ ] Add ENCODE/miRNA-seq co-expression edges
- [ ] Download HGNC file automatically on first run
- [ ] Support edge weighting in model trainer
