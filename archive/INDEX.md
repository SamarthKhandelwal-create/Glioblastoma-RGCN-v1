# ceRNA Testing Suite - Index & Navigation

## 📍 You Are Here

This is the complete testing suite for ceRNA (competing endogenous RNA) edge inference with miRNA injection support.

---

## 🎯 Your Task

```powershell
$env:INJECT_MIRNA='true'
python test_cerna_inference.py --inject --run --verbose
```

---

## 📚 Documentation Structure

### Entry Points

**New to this? Start here:**
1. **[START_HERE.md](START_HERE.md)** ← **READ THIS FIRST**
   - What you got
   - Your exact command
   - 30-second quick start
   - Success examples

### Reference Guides

**Quick lookup:**
- **[QUICK_REFERENCE.md](QUICK_REFERENCE.md)** — Commands, expected results, checklist
- **[README_CERNA_TESTING.md](README_CERNA_TESTING.md)** — Getting started, troubleshooting

**Detailed documentation:**
- **[TEST_CERNA_GUIDE.md](TEST_CERNA_GUIDE.md)** — Comprehensive guide with 9 test breakdown
- **[IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md)** — Technical architecture

---

## 🧪 Test Scripts

### Main Test Suite
```powershell
python test_cerna_inference.py --inject --run --verbose
```
- Runs 9 validation tests
- Supports INJECT_MIRNA mode
- Returns exit code (0=pass, 1=fail)

### Convenience Wrappers
```powershell
# PowerShell (Windows)
.\run_exact_test.ps1 -Run -Verbose

# Python (cross-platform)
python run_tests_example.py
```

---

## 📊 What Gets Tested

| Test | Validates | Location |
|------|-----------|----------|
| 1. File Existence | Output files present | results/ |
| 2. Node Features | Feature matrix quality | node_features_matrix.csv |
| 3. Interactions | Interaction structure | interactions.csv |
| 4. Graph Structure | Node/edge types | hetero_graph_GBM.pt |
| 5. Node Mapping | Gene index mapping | node_mapping.json |
| 6. ceRNA Logic | Direct + inferred edges | Graph analysis |
| 7. Edge Provenance | Metadata + confidence | edge_metadata.csv |
| 8. ID Normalization | Ensembl/miRNA format | Data validation |
| 9. INJECT_MIRNA | miRNA injection | Feature matrix |

---

## 🎬 Quick Start Scenarios

### "Just run it"
```powershell
$env:INJECT_MIRNA='true'
python test_cerna_inference.py --inject --run
```
**Time:** ~2-5 minutes

### "I want to see details"
```powershell
$env:INJECT_MIRNA='true'
python test_cerna_inference.py --inject --run --verbose
```
**Time:** ~2-5 minutes + detailed output

### "I have existing graph"
```powershell
python test_cerna_inference.py --inject
```
**Time:** ~10-30 seconds

### "I'm using PowerShell wrapper"
```powershell
.\run_exact_test.ps1 -Run -Verbose
```
**Time:** ~2-5 minutes

---

## 🔍 Key Files Generated

After running tests:

| File | Purpose | Size |
|------|---------|------|
| `hetero_graph_GBM.pt` | PyG HeteroData object | ~5-50 MB |
| `node_mapping.json` | Gene → index mapping | ~50-100 KB |
| `edge_metadata.csv` | Edge provenance + confidence | ~1-5 MB |
| `node_features_matrix.csv` | Feature matrix (z-scored) | ~5-20 MB |
| `interactions.csv` | Raw interactions | ~500 KB |

---

## 🎯 Expected Results

### With INJECT_MIRNA=true (Recommended)
- **Nodes**: 400-550
- **Direct edges**: 1,200-1,500
- **ceRNA edges**: 100-400 (NEW!)
- **Total edges**: 1,300-1,900
- **miRNA coverage**: 20-30%
- **Exit code**: 0 (if tests pass)

### With INJECT_MIRNA=false
- **Nodes**: 350-450
- **Direct edges**: 1,000-1,300
- **ceRNA edges**: 50-200
- **Total edges**: 1,050-1,500
- **miRNA coverage**: 10-20%

---

## 📋 Test Interpretation

### ✓ Success (Exit code: 0)
```
✓ PASSED: 25+
✗ FAILED: 0
⚠ WARNINGS: 0-2
```
→ Graph ready for training

### ✗ Failure (Exit code: 1)
```
✓ PASSED: 18-24
✗ FAILED: 1+
```
→ Check error messages, troubleshoot

### ⚠ Warnings (Non-fatal)
```
⚠ WARNINGS: 1-3
```
→ Suboptimal but working

---

## 🔧 Troubleshooting Fast Path

| Issue | Solution | Doc |
|-------|----------|-----|
| Exit code 1 | Read test failures | test_cerna_inference.py output |
| No ceRNA edges | Enable INJECT_MIRNA or check data | TEST_CERNA_GUIDE.md |
| ID errors | Check normalize_id() | src/graph/build_graph.py |
| File not found | Run `--run` to build | START_HERE.md |
| Python error | Install dependencies | TEST_CERNA_GUIDE.md |

---

## 🚀 Next Steps After Success

```powershell
# 1. Verify graph statistics
python verify_graph.py

# 2. Train model with enhanced graph
python src/training/train_model.py \
    --graph-path results/hetero_graph_GBM.pt \
    --epochs 100 \
    --seed 42

# 3. Or run cross-validation
python src/training/train_cross_validation.py \
    --graph-path results/hetero_graph_GBM.pt \
    --folds 10 \
    --epochs 100 \
    --seed 42
```

---

## 📍 Documentation Index

| File | Purpose | When to Read |
|------|---------|--------------|
| **START_HERE.md** | Overview & quick start | First time setup |
| **QUICK_REFERENCE.md** | Command cheat sheet | Fast lookup |
| **README_CERNA_TESTING.md** | Getting started guide | Learning how to run |
| **TEST_CERNA_GUIDE.md** | Comprehensive reference | Detailed information |
| **IMPLEMENTATION_SUMMARY.md** | Technical architecture | Understanding design |
| **INDEX.md** | This file | Navigation |

---

## 🎓 Concepts Explained

### Direct Edges
Curated interactions from literature:
- miRNA silences mRNA
- lncRNA sequesters miRNA

### Inferred ceRNA Edges
Computational inference:
- If mRNA1 and mRNA2 both regulated by miR-X
- Create edge between them
- They "compete" for miRNA

### INJECT_MIRNA
Inclusion strategy:
- Without: Only miRNAs with expression data
- With: All miRNAs mentioned in interactions
- Effect: More inferred edges

---

## ✅ Verification Checklist

Before running tests:
- [ ] Python 3.8+ installed
- [ ] Dependencies installed (`torch`, `pandas`, `numpy`, `torch_geometric`)
- [ ] `results/` directory writable
- [ ] `src/` directory with source code
- [ ] Enough disk space (~100 MB)

After running tests:
- [ ] Exit code is 0
- [ ] 25+ tests passed
- [ ] 0 tests failed
- [ ] ceRNA edges detected
- [ ] Graph saved to `results/hetero_graph_GBM.pt`

---

## 🆘 Getting Help

1. **Quick answers**: QUICK_REFERENCE.md
2. **How-to guide**: TEST_CERNA_GUIDE.md
3. **Troubleshooting**: TEST_CERNA_GUIDE.md → Troubleshooting section
4. **Technical details**: IMPLEMENTATION_SUMMARY.md
5. **Code inspection**: test_cerna_inference.py

---

## 🎬 Action Now

```powershell
# Run this:
$env:INJECT_MIRNA='true'
python test_cerna_inference.py --inject --run

# Then:
# - Check exit code
# - If 0: success ✓
# - If 1: review error output and troubleshoot
```

---

## 📊 Pipeline Overview

```
You → Script (this file)
↓
Set INJECT_MIRNA=true
↓
Build features (optional)
↓
Build graph + ceRNA inference
↓
Run 9 validation tests
↓
Exit code 0/1
↓
Ready for training (if 0)
```

---

## 🔗 Related Commands

```powershell
# Build only (no test)
python src/graph/build_graph.py

# Test only (no build)
python test_cerna_inference.py --inject

# Build + test (recommended)
python test_cerna_inference.py --inject --run

# Verbose output (for debugging)
python test_cerna_inference.py --inject --run --verbose

# Using wrapper scripts
.\run_exact_test.ps1 -Run
python run_tests_example.py
```

---

## 💾 Output Artifacts

After successful run, find in `results/`:

```
results/
├── hetero_graph_GBM.pt          ← Main graph (PyG HeteroData)
├── node_mapping.json             ← Gene index mapping
├── edge_metadata.csv             ← Edge provenance + confidence
├── node_features_matrix.csv      ← Feature matrix (z-scored)
├── interactions.csv              ← Raw interactions
├── DE_*.tsv                      ← Differential expression files
└── miRNA_features.csv (optional) ← miRNA-specific features
```

---

**Ready? Read [START_HERE.md](START_HERE.md) then run:**
```powershell
$env:INJECT_MIRNA='true'
python test_cerna_inference.py --inject --run --verbose
```
