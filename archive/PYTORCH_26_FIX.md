# PyTorch 2.6 Compatibility Fix

## Issue

When running `test_cerna_inference.py`, the test suite crashed on Test 4 (Graph Structure) with:

```
_pickle.UnpicklingError: Weights only load failed...
Unsupported global: GLOBAL torch_geometric.data.data.DataEdgeAttr was not an allowed global by default.
```

**Root cause:** PyTorch 2.6+ introduced stricter security defaults for `torch.load()`, requiring explicit allowlisting of custom classes.

---

## Solution Applied

### 1. Import Required Classes
Added to test_cerna_inference.py:
```python
from torch_geometric.data.data import DataEdgeAttr

# PyTorch 2.6+ security: allow torch_geometric objects
try:
    torch.serialization.add_safe_globals([DataEdgeAttr])
except (AttributeError, TypeError):
    # Fallback for older PyTorch versions
    pass
```

### 2. Updated torch.load() Calls
Changed all `torch.load()` calls to use `weights_only=False`:
```python
try:
    hetero_data = torch.load(graph_path, weights_only=False)
except TypeError:
    # Fallback for older PyTorch versions without weights_only parameter
    hetero_data = torch.load(graph_path)
```

---

## How to Use

### Option 1: Quick Re-run
```powershell
.\run_test_fixed.ps1 -Verbose
```

### Option 2: Direct Command
```powershell
$env:INJECT_MIRNA='true'
python test_cerna_inference.py --inject --verbose
```

---

## What Was Changed

| File | Change |
|------|--------|
| `test_cerna_inference.py` | Added PyTorch 2.6 compatibility fixes |
| `run_test_fixed.ps1` | NEW: Convenience wrapper |

---

## Expected Output

```
✓ [INFO]
=== Test 4: Graph Structure ===
→ [DEBUG] Loaded graph with Data(x=[2075, 10], edge_index=[2, 1344], edge_attr=[1344])
✓ [INFO] Has expected node types: PASSED
✓ [INFO] Graph has edges: PASSED

... (continues with remaining tests)
```

---

## Verification

After running the test:
- [ ] Exit code is 0 (success)
- [ ] All 9 tests pass
- [ ] Graph structure validated
- [ ] ceRNA edges detected

---

## Technical Details

**Why this happens:**
- PyTorch 2.6 changed `torch.load()` default to `weights_only=True`
- This prevents loading pickled Python objects (like torch_geometric's custom classes)
- Solution: Explicitly set `weights_only=False` for trusted files

**Why the try-except:**
- Older PyTorch versions don't have the `weights_only` parameter
- The fallback `torch.load(graph_path)` works with PyTorch <2.6

---

**Ready? Run:**
```powershell
.\run_test_fixed.ps1 -Verbose
```
