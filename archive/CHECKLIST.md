# ✅ Deployment Checklist

## Pre-Execution

- [ ] Python 3.8+ installed
  ```powershell
  python --version
  ```

- [ ] Dependencies installed
  ```powershell
  pip list | grep -E 'torch|pandas|numpy|torch_geometric'
  ```

- [ ] Project directory accessible
  ```powershell
  ls src/preprocessing
  ls src/graph
  ```

- [ ] Results directory writable
  ```powershell
  ls -d results/
  ```

- [ ] Sufficient disk space (~100 MB)

---

## Initial Setup (Do Once)

- [ ] Navigate to project root
  ```powershell
  cd path/to/Glioblastoma-RGCN-v1
  ```

- [ ] Read this file (you're doing it!)

- [ ] Read START_HERE.md
  ```powershell
  notepad START_HERE.md
  ```

---

## Execution

### Option 1: Direct Command (Recommended)
- [ ] Open PowerShell
- [ ] Set environment variable:
  ```powershell
  $env:INJECT_MIRNA='true'
  ```
- [ ] Run test suite:
  ```powershell
  python test_cerna_inference.py --inject --run --verbose
  ```
- [ ] Monitor output for progress

### Option 2: Using PowerShell Script
- [ ] Run:
  ```powershell
  .\run_exact_test.ps1 -Run -Verbose
  ```

### Option 3: Using Python Script
- [ ] Run:
  ```powershell
  python run_tests_example.py
  ```

---

## Monitoring Execution

- [ ] Wait 2-5 minutes for completion
- [ ] Watch for test progress indicators:
  ```
  ✓ = Test passed
  ✗ = Test failed
  ⚠ = Warning (non-critical)
  → = Debug information
  ```

- [ ] Look for these sections:
  - "Building node features..." (if running --run)
  - "Building graph with ceRNA inference..." (if running --run)
  - "=== Test N: ..." (for each of 9 tests)
  - "TEST SUMMARY" (final results)

---

## After Execution

### Check Exit Code
```powershell
# PowerShell shows exit code
$LASTEXITCODE

# Should be:
# 0 = Success ✓
# 1 = Failure ✗
```

### Success Case (Exit Code 0)
- [ ] All 9 tests passed or nearly all passed
- [ ] FAILED count = 0
- [ ] Graph ready for training
- [ ] Proceed to model training

### Failure Case (Exit Code 1)
- [ ] Note the failing tests
- [ ] Read TEST_CERNA_GUIDE.md troubleshooting section
- [ ] Check specific error messages in output
- [ ] Verify:
  - [ ] Feature matrix exists: `results/node_features_matrix.csv`
  - [ ] Interactions exist: `results/interactions.csv`
  - [ ] miRNA data is present (if INJECT_MIRNA=true)

---

## Verification

### Quick Stats Check
```powershell
# Check node count
(Import-Csv results/node_features_matrix.csv).Count

# Check interaction count
(Import-Csv results/interactions.csv).Count
```

### Expected Ranges
- [ ] Nodes: 400-550
- [ ] Direct edges: 1,200-1,500
- [ ] ceRNA edges: 100-400
- [ ] Total edges: 1,300-1,900

### Graph File
- [ ] Check size: `ls -la results/hetero_graph_GBM.pt` (~5-50 MB)
- [ ] Check timestamp: File should be recent
- [ ] Verify readable: `file results/hetero_graph_GBM.pt`

---

## Common Issues Checklist

### Issue: "Python not found"
- [ ] Check Python path: `where python`
- [ ] Or use full venv path if needed
- [ ] Solution: Install Python or activate venv

### Issue: "Module not found (torch, pandas, etc)"
- [ ] Install dependencies: `pip install -r requirements.txt`
- [ ] Check installation: `pip list`
- [ ] Solution: Run pip install again

### Issue: "File not found (interactions.csv, etc)"
- [ ] Check results/ directory exists
- [ ] Check file permissions
- [ ] Run with `--run` flag to build missing files
- [ ] Solution: `python test_cerna_inference.py --inject --run`

### Issue: "No ceRNA edges found"
- [ ] Check interaction data quality
- [ ] Enable INJECT_MIRNA for more coverage
- [ ] Verify miRNA names have 'mir' or 'let'
- [ ] Solution: `$env:INJECT_MIRNA='true'; python src/graph/build_graph.py`

### Issue: "Test fails on ID normalization"
- [ ] Check normalize_id() function in src/graph/build_graph.py
- [ ] Ensure Ensembl IDs don't have versions
- [ ] Ensure miRNA names are lowercase
- [ ] Solution: Review build_graph.py implementation

---

## Next Steps (After Success)

### Immediate (If Exit Code 0)
- [ ] Review graph statistics
- [ ] Check edge metadata: `head results/edge_metadata.csv`
- [ ] Verify mapping: `head results/node_mapping.json`

### Short Term (Model Training)
```powershell
# Single training run
python src/training/train_model.py \
    --graph-path results/hetero_graph_GBM.pt \
    --epochs 100

# Cross-validation
python src/training/train_cross_validation.py \
    --graph-path results/hetero_graph_GBM.pt \
    --folds 10 --epochs 100
```

### Medium Term (Analysis)
- [ ] Monitor training metrics
- [ ] Track edge contribution to predictions
- [ ] Compare with/without ceRNA edges

---

## Documentation Reference

If you need to:
- **Get started**: Read START_HERE.md
- **Find a command**: Check QUICK_REFERENCE.md
- **Understand tests**: See TEST_CERNA_GUIDE.md
- **Deep dive**: Read IMPLEMENTATION_SUMMARY.md
- **Navigate**: Use INDEX.md
- **Troubleshoot**: Check TEST_CERNA_GUIDE.md → Troubleshooting

---

## Rollback/Reset

If you need to start over:

```powershell
# Remove old graph
rm results/hetero_graph_GBM.pt
rm results/edge_metadata.csv

# Rebuild everything
$env:INJECT_MIRNA='true'
python src/preprocessing/build_node_features.py
python src/graph/build_graph.py
python test_cerna_inference.py --inject --verbose
```

---

## Support Resources

| Question | Resource |
|----------|----------|
| "How do I run this?" | START_HERE.md |
| "What command should I use?" | QUICK_REFERENCE.md |
| "Why did test fail?" | TEST_CERNA_GUIDE.md |
| "How does this work?" | IMPLEMENTATION_SUMMARY.md |
| "Where do I find...?" | INDEX.md |

---

## Quick Command Reference

```powershell
# Recommended: Full workflow with output
$env:INJECT_MIRNA='true'
python test_cerna_inference.py --inject --run --verbose

# Minimal: Just test (no build)
python test_cerna_inference.py --inject

# Using scripts
.\run_exact_test.ps1 -Run
python run_tests_example.py

# Build only (no test)
python src/graph/build_graph.py

# After success: Train model
python src/training/train_model.py --graph-path results/hetero_graph_GBM.pt
```

---

## Success Checklist (Final)

After running and getting exit code 0:

- [ ] Tests completed without errors
- [ ] Graph file created: `results/hetero_graph_GBM.pt`
- [ ] Node mapping saved: `results/node_mapping.json`
- [ ] Edge metadata exported: `results/edge_metadata.csv`
- [ ] Edge count appropriate: 1,300-1,900
- [ ] ceRNA edges present: 100-400
- [ ] All output files readable
- [ ] Ready to proceed with model training

---

## Go!

```powershell
# THE COMMAND:
$env:INJECT_MIRNA='true'
python test_cerna_inference.py --inject --run --verbose

# Expected: Success! ✓
# Time: 2-5 minutes
```

---

**Status: ✅ READY TO EXECUTE**

Start with the command above. Check back here if any issues.
