# ceRNA Inference Testing Guide

## Quick Start

### 1. Basic Test (Existing Graph)
```powershell
python test_cerna_inference.py
```

### 2. Test with miRNA Injection Enabled
```powershell
$env:INJECT_MIRNA='true'
python test_cerna_inference.py --inject
```

### 3. Full Pipeline + Test
```powershell
$env:INJECT_MIRNA='true'
python test_cerna_inference.py --inject --run
```

---

## Test Coverage

### Test 1: File Existence
- ✓ Checks for required output files
- ✓ Feature matrix, interactions, graph, mappings

### Test 2: Node Features
- ✓ Feature matrix shape and content
- ✓ Numeric feature validation
- ✓ NaN value detection
- ✓ miRNA presence (when INJECT_MIRNA=true)

### Test 3: Interactions
- ✓ Interaction file structure
- ✓ miRNA-target pair categorization
- ✓ Interaction count validation

### Test 4: Graph Structure
- ✓ Node type presence (mRNA, lncRNA, miRNA)
- ✓ Edge type validation
- ✓ ceRNA edge detection
- ✓ Total edge count

### Test 5: Node Mapping
- ✓ Gene-to-index mapping
- ✓ Node type mapping
- ✓ Type distribution

### Test 6: ceRNA Inference Logic
- ✓ Direct miRNA-target edge count
- ✓ Inferred ceRNA edge count
- ✓ Logical consistency

### Test 7: Edge Provenance
- ✓ Edge metadata file validation
- ✓ Provenance source tracking (curated, TargetScan, ENCORI)
- ✓ Confidence score statistics

### Test 8: ID Normalization
- ✓ Ensembl IDs stripped of versions
- ✓ miRNA IDs lowercase
- ✓ Consistent ID formatting

### Test 9: INJECT_MIRNA Mode
- ✓ miRNA injection validation
- ✓ Zero-feature miRNA rows
- ✓ Count verification

---

## Command Examples

### Example 1: Quick Validation (No Injection)
```bash
python test_cerna_inference.py --verbose
```

### Example 2: Full Workflow with Injection
```powershell
# Set environment variable
$env:INJECT_MIRNA='true'

# Run full pipeline
python src/preprocessing/build_node_features.py
python src/graph/build_graph.py

# Test ceRNA inference
python test_cerna_inference.py --inject --verbose
```

### Example 3: Custom Results Directory
```powershell
$env:INJECT_MIRNA='true'
$env:RESULTS_DIR='custom_results'
python test_cerna_inference.py --results-dir custom_results --verbose
```

### Example 4: Build and Test
```powershell
$env:INJECT_MIRNA='true'
python test_cerna_inference.py --inject --run --verbose
```

---

## Understanding Test Output

### Symbols
- ✓ **PASSED** - Test succeeded
- ✗ **FAILED** - Test failed (exit code 1)
- ⚠ **WARNING** - Suboptimal but not critical
- → **DEBUG** - Informational message

### Example Output
```
→ [DEBUG] Feature matrix shape: (500, 100)
✓ [INFO] Non-empty feature matrix: PASSED
⚠ [WARN] Sufficient miRNAs injected: Expected miRNAs with INJECT_MIRNA=true, found 45
→ [DEBUG] Direct miRNA-target edges: 1256
→ [DEBUG] Inferred ceRNA edges: 150
```

### Exit Codes
- **0** - All tests passed
- **1** - One or more tests failed

---

## Troubleshooting

### "No ceRNA edges found"
**Cause:** ceRNA inference requires at least 2 genes sharing a miRNA regulator.

**Solution:**
```powershell
# Check interaction count
(Import-Csv results/interactions.csv).Count

# Enable miRNA injection for more coverage
$env:INJECT_MIRNA='true'
python src/graph/build_graph.py
python test_cerna_inference.py --inject
```

### "Versioned Ensembl IDs found"
**Cause:** ID normalization not applied correctly.

**Solution:** Check `normalize_id()` in `src/graph/build_graph.py`:
```python
def normalize_id(gene_id):
    s = str(gene_id).strip()
    if s.upper().startswith("ENSG"):
        return s.split(".")[0]  # Strip version
    return s.lower()
```

### "No miRNAs in matrix (INJECT_MIRNA=true)"
**Cause:** No miRNAs found in interactions or feature file.

**Solution:**
1. Check `results/interactions.csv` has miRNA rows
2. Verify miRNA names contain 'mir' or 'let'
3. Check `results/miRNA_features.csv` exists (optional)

---

## Integration with Training Pipeline

After validation, use the enhanced graph for training:

```powershell
$env:INJECT_MIRNA='true'
python src/preprocessing/build_node_features.py
python src/graph/build_graph.py

# Validate before training
python test_cerna_inference.py --inject

# Run cross-validation with enhanced graph
python src/training/train_cross_validation.py `
    --graph-path results/hetero_graph_GBM.pt `
    --outdir cv_results `
    --folds 10 `
    --epochs 100 `
    --seed 42
```

---

## Expected Results

### With INJECT_MIRNA=false (Default)
- Node count: ~400-500
- Edge count: ~1,200-1,400
- ceRNA edges: ~100-300
- miRNA count: Lower (only from feature matrix)

### With INJECT_MIRNA=true (Recommended)
- Node count: ~400-550 (more miRNAs)
- Edge count: ~1,400-1,600 (more inferred edges)
- ceRNA edges: ~200-400 (more co-regulation)
- miRNA count: ~100-150 (fuller coverage)

---

## Advanced: Running with Different Data Sources

```powershell
# With TargetScan predictions
$env:INJECT_MIRNA='true'
python src/preprocessing/build_node_features.py
python src/graph/integrate_targetscan.py  # If available
python src/graph/build_graph.py
python test_cerna_inference.py --inject --verbose

# With ENCORI data
python src/graph/fetch_encori.py --api
python src/graph/build_graph.py
python test_cerna_inference.py --verbose
```

---

## Metrics to Monitor

After running tests, check these key metrics:

1. **Edge Expansion Ratio**: `total_edges / direct_mirna_edges`
   - Expected: 1.2-1.5x (20-50% more edges from ceRNA)

2. **miRNA Coverage**: `mirna_count / total_nodes`
   - Expected: 20-30% with INJECT_MIRNA=true

3. **ceRNA Ratio**: `cerna_edges / total_edges`
   - Expected: 10-30% of edges are inferred

4. **ID Success Rate**: `(total_nodes - unknown_type) / total_nodes`
   - Expected: >95% properly typed nodes

---

## File Artifacts Generated

| File | Purpose | Format |
|------|---------|--------|
| `results/hetero_graph_GBM.pt` | PyG HeteroData object | PyTorch |
| `results/node_mapping.json` | Gene index + type mapping | JSON |
| `results/edge_metadata.csv` | Edge provenance & confidence | CSV |
| `results/node_features_matrix.csv` | Feature matrix (z-scored) | CSV |
| `results/interactions.csv` | Raw interaction pairs | CSV |

---

## Questions?

For more information:
- See `ENHANCED_GRAPH_BUILDING.md` for architecture details
- Check `src/graph/build_graph.py` for implementation
- Review `src/training/train_model.py` for model integration
