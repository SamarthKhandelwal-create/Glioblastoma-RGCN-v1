# Testing ceRNA Inference with miRNA Injection

## Your Command

You requested to run:
```powershell
$env:INJECT_MIRNA='true'
python test_cerna_inference.py --inject
```

Here's exactly how to do it:

---

## ⚡ Quick Start (30 seconds)

### Option 1: Direct PowerShell Command
```powershell
$env:INJECT_MIRNA='true'
python test_cerna_inference.py --inject --run
```

### Option 2: Use Convenience Script (Recommended)
```powershell
.\run_exact_test.ps1 -Run
```

### Option 3: Step-by-Step Manual
```powershell
# Step 1: Build features
python src/preprocessing/build_node_features.py

# Step 2: Build graph with ceRNA inference
$env:INJECT_MIRNA='true'
python src/graph/build_graph.py

# Step 3: Run tests
python test_cerna_inference.py --inject
```

---

## 📊 What This Tests

### Direct Edges (Curated miRNA-target Interactions)
- From: `results/interactions.csv`
- Count: ~1,200-1,500 edges
- Type: miRNA → mRNA (silencing), lncRNA → miRNA (sequestering)

### Inferred ceRNA Edges (NEW)
- From: Shared miRNA co-targeting analysis
- Count: ~200-400 edges
- Type: lncRNA ⟷ mRNA co-regulation
- Mechanism: If both share miRNAs, they compete for those miRNAs

### With INJECT_MIRNA=true
- Additional zero-valued miRNA nodes created
- Enables inference of more ceRNA relationships
- Result: More connected graph for richer representation

---

## 🎯 Expected Output

```
→ [DEBUG] Feature matrix shape: (500, 100)
✓ [INFO] Non-empty feature matrix: PASSED
→ [DEBUG] Loaded 1500 interactions
✓ [INFO] Node Type Assignment: mRNA=280, lncRNA=45, miRNA=125, Unknown=0

--- Phase 5.5: Inferring ceRNA Edges (Shared miRNA Co-targeting) ---
→ [DEBUG] Built miRNA target map with 125 unique miRNAs
✓ [INFO] Inferred 287 ceRNA edges from shared miRNA co-targeting
✓ [INFO] Added 287 new ceRNA edges (after deduplication)
→ [DEBUG] Total edges after ceRNA inference: 1543

→ [DEBUG] Direct miRNA-target edges: 1256
→ [DEBUG] Inferred ceRNA edges: 287

======================================================================
  TEST SUMMARY
======================================================================
✓ PASSED: 25
✗ FAILED: 0
⚠ WARNINGS: 0
======================================================================
```

---

## 🔍 Test Suite Details

9 comprehensive tests:

| # | Test | Validates |
|---|------|-----------|
| 1 | File Existence | All required outputs present |
| 2 | Node Features | Feature matrix, NaN check, miRNA count |
| 3 | Interactions | Structure, miRNA-target pairs |
| 4 | Graph Structure | Node/edge types, ceRNA detection |
| 5 | Node Mapping | Gene index, type mapping |
| 6 | ceRNA Logic | Direct vs inferred edge counts |
| 7 | Edge Provenance | Metadata, confidence scores |
| 8 | ID Normalization | Ensembl version stripping, miRNA lowercasing |
| 9 | INJECT_MIRNA Mode | miRNA injection validation |

---

## 📁 Files Created

| File | Purpose | Run |
|------|---------|-----|
| `test_cerna_inference.py` | Main test suite | `python test_cerna_inference.py --inject` |
| `run_exact_test.ps1` | PowerShell wrapper | `.\run_exact_test.ps1 -Run` |
| `run_tests_example.py` | Python wrapper | `python run_tests_example.py` |
| `TEST_CERNA_GUIDE.md` | Full documentation | Reference |
| `IMPLEMENTATION_SUMMARY.md` | Technical details | Reference |
| `README_CERNA_TESTING.md` | This file | Getting started |

---

## 🚀 Typical Workflow

```
1. Build features
   python src/preprocessing/build_node_features.py

2. Set environment and build graph
   $env:INJECT_MIRNA='true'
   python src/graph/build_graph.py

3. Test ceRNA inference
   python test_cerna_inference.py --inject --verbose

4. If tests pass (exit code 0):
   python src/training/train_model.py \
       --graph-path results/hetero_graph_GBM.pt
```

---

## ✅ Success Indicators

- **Exit Code**: 0 (success) or 1 (failure)
- **PASSED count**: Should be 25+ out of ~27 tests
- **FAILED count**: Should be 0
- **ceRNA edges**: Should be 100-400 (depends on data)
- **Total edges**: Should be 1,200-1,600

---

## 🐛 Troubleshooting

### "No ceRNA edges inferred"
```powershell
# Check that you have direct edges first
$env:INJECT_MIRNA='true'
python src/graph/build_graph.py
python test_cerna_inference.py --inject --verbose
```

### "Test fails on ID normalization"
Check that `src/graph/build_graph.py` has:
```python
def normalize_id(gene_id):
    s = str(gene_id).strip()
    if s.upper().startswith("ENSG"):
        return s.split(".")[0]  # Strip .15 version
    return s.lower()  # Lowercase miRNAs
```

### "Python not found"
```powershell
# Check Python is installed
python --version

# Or use full path if in venv
C:\path\to\venv\Scripts\python.exe test_cerna_inference.py --inject
```

---

## 📚 Full Documentation

- **`TEST_CERNA_GUIDE.md`** — Comprehensive guide with examples
- **`IMPLEMENTATION_SUMMARY.md`** — Technical architecture
- **`src/graph/build_graph.py`** — Implementation details

---

## 🎓 Key Concepts

### Direct Edges
Curated interactions from literature/databases:
- miRNA silences mRNA
- lncRNA sequesters (competes with) mRNA for miRNA binding

### Inferred ceRNA Edges
Computational inference from co-regulation:
- If mRNA1 and mRNA2 both targeted by miRNA-X
- Then they "compete" for that miRNA
- Creates indirect biological link

### INJECT_MIRNA Flag
When `true`:
- Finds all miRNAs mentioned in interactions.csv
- Adds them to feature matrix with zero features
- Enables inference of co-regulation patterns

---

## Next Steps

After successful test run:

1. **Check results**
   ```powershell
   ls results/
   # Should include:
   # - hetero_graph_GBM.pt (PyG Data object)
   # - node_mapping.json (metadata)
   # - edge_metadata.csv (provenance)
   ```

2. **Train model**
   ```powershell
   python src/training/train_model.py \
       --graph-path results/hetero_graph_GBM.pt \
       --epochs 100
   ```

3. **Run cross-validation**
   ```powershell
   python src/training/train_cross_validation.py \
       --graph-path results/hetero_graph_GBM.pt \
       --folds 10 --epochs 100
   ```

---

## Command Reference

| Task | Command |
|------|---------|
| Quick test | `python test_cerna_inference.py --inject` |
| With output | `python test_cerna_inference.py --inject --verbose` |
| Full pipeline + test | `python test_cerna_inference.py --inject --run --verbose` |
| Using script | `.\run_exact_test.ps1 -Run -Verbose` |
| Custom dir | `python test_cerna_inference.py --inject --results-dir custom_results` |

---

**Ready? Run:**
```powershell
$env:INJECT_MIRNA='true'
python test_cerna_inference.py --inject --run
```
