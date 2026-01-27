# Quick Reference Card: ceRNA Inference Testing

## Your Command
```powershell
$env:INJECT_MIRNA='true'
python test_cerna_inference.py --inject
```

---

## Run Variations

| Need | Command |
|------|---------|
| **Just test** | `python test_cerna_inference.py --inject` |
| **Build + test** | `python test_cerna_inference.py --inject --run` |
| **With details** | `python test_cerna_inference.py --inject --verbose` |
| **Full workflow** | `python test_cerna_inference.py --inject --run --verbose` |
| **PowerShell** | `.\run_exact_test.ps1 -Run -Verbose` |

---

## Environment Setup

### One-liner (Recommended)
```powershell
$env:INJECT_MIRNA='true'; python test_cerna_inference.py --inject --run
```

### Step-by-step
```powershell
# Step 1: Set environment
$env:INJECT_MIRNA='true'

# Step 2: Build features (if needed)
python src/preprocessing/build_node_features.py

# Step 3: Build graph (if needed)
python src/graph/build_graph.py

# Step 4: Test
python test_cerna_inference.py --inject
```

---

## Expected Results

✓ **SUCCESS** (Exit code: 0)
- 25+ tests passed
- 0 tests failed
- ceRNA edges: 100-400
- Total edges: 1,200-1,600

✗ **FAILURE** (Exit code: 1)
- Check output for specific failures
- Review TEST_CERNA_GUIDE.md troubleshooting

---

## Test Coverage

| Component | Tested | Method |
|-----------|--------|--------|
| Nodes | ✓ | Feature matrix validation |
| Direct edges | ✓ | Interaction count |
| ceRNA edges | ✓ | Inference logic check |
| ID normalization | ✓ | Ensembl/miRNA format |
| Types (mRNA/lncRNA/miRNA) | ✓ | DE file membership |
| Provenance | ✓ | Edge metadata |
| INJECT_MIRNA | ✓ | miRNA coverage |

---

## Files Generated

```
results/
├── hetero_graph_GBM.pt          ← Use this for training
├── node_mapping.json
├── edge_metadata.csv
├── node_features_matrix.csv
└── interactions.csv
```

---

## Integration with Training

```powershell
# After successful test:
python src/training/train_model.py \
    --graph-path results/hetero_graph_GBM.pt \
    --epochs 100 \
    --seed 42
```

---

## Troubleshooting Checklist

- [ ] Python installed: `python --version`
- [ ] Dependencies: `pip list | grep torch`
- [ ] Feature matrix exists: `ls results/node_features_matrix.csv`
- [ ] Interactions exist: `ls results/interactions.csv`
- [ ] Permissions: Can write to `results/` directory
- [ ] miRNA data: Check if INJECT_MIRNA needed

---

## Key Metrics

After test completes:

| Metric | Check |
|--------|-------|
| Nodes | Feature matrix rows |
| Features | Feature matrix columns |
| Direct edges | Should be 1,200+ |
| ceRNA edges | Should be 100-400 |
| miRNA coverage | 20-30% of nodes |

---

## Common Commands

### Build everything from scratch
```powershell
$env:INJECT_MIRNA='true'
python src/preprocessing/build_node_features.py
python src/graph/build_graph.py
python test_cerna_inference.py --inject --verbose
```

### Just validate existing graph
```powershell
python test_cerna_inference.py --inject
```

### Test with verbose output (debug)
```powershell
$env:INJECT_MIRNA='true'
python test_cerna_inference.py --inject --verbose
```

### Check graph statistics
```powershell
python verify_graph.py
```

---

## Documentation

| Doc | For |
|-----|-----|
| `README_CERNA_TESTING.md` | Getting started |
| `TEST_CERNA_GUIDE.md` | Detailed reference |
| `IMPLEMENTATION_SUMMARY.md` | Technical details |
| `test_cerna_inference.py` | Test code |

---

## Questions?

1. Check test output (most issues explained)
2. Review `TEST_CERNA_GUIDE.md` Troubleshooting section
3. Inspect results files directly

---

## Next Steps (After Tests Pass)

```powershell
# Option 1: Single training run
python src/training/train_model.py --graph-path results/hetero_graph_GBM.pt

# Option 2: Cross-validation
python src/training/train_cross_validation.py \
    --graph-path results/hetero_graph_GBM.pt \
    --folds 10 --epochs 100 --seed 42
```

---

**TL;DR: Run this now**
```powershell
$env:INJECT_MIRNA='true'
python test_cerna_inference.py --inject --run --verbose
```
