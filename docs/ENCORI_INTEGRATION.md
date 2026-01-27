# Using ENCORI for Real ceRNA Interactions

## Quick Start

### Option 1: Use Literature-Based ENCORI (Recommended - No API Call)
```powershell
# Create curated ENCORI interactions from literature
python src/preprocessing/fetch_encori.py --fallback

# Rebuild graph with ENCORI + TargetScan edges
$env:INJECT_MIRNA='true'
python src/graph/build_graph.py

# Verify
python verify_graph.py
```

### Option 2: Try to Fetch from ENCORI API (Experimental)
```powershell
# Attempt to fetch real-time data from ENCORI API
python src/preprocessing/fetch_encori.py --api

# Rebuild
$env:INJECT_MIRNA='true'
python src/graph/build_graph.py
```

## What ENCORI Provides

**ENCORI (starBase)** is a comprehensive ceRNA database with:
- ✅ Experimentally validated miRNA-target interactions
- ✅ CLIP-seq evidence (high confidence)
- ✅ Multi-tissue expression data
- ✅ Gene type annotations (mRNA, lncRNA, miRNA)

## Expected Results

After integrating ENCORI:
```
TargetScan predictions:      88 edges
ENCORI validated:           13+ edges (literature-based)
Total graph edges:       1,350-1,400+
```

## Data Sources Used

| Source | Type | Edges | Confidence |
|--------|------|-------|-----------|
| interactions.csv (curated) | Direct | 1,255 | 1.0 |
| ENCORI (validated) | Experimental | 13+ | 0.85-0.95 |
| TargetScan (predicted) | Computational | 88+ | 0.6-1.0 |

## Files Generated

- `results/encori_interactions.csv` — ENCORI miRNA-target pairs
- `results/hetero_graph_GBM.pt` — Updated graph with all sources
- `results/edge_metadata.csv` — Provenance tracking

## Literature Sources (Fallback)

The fallback dataset includes well-documented GBM interactions:
- **miR-21**: PTEN, PDCD4, CDKN1A, LATS2 (silencing targets)
- **miR-155**: TP53, SOCS1, FOXO3A
- **miR-10b**: CDKN1B, E2F5
- **miR-221**: Cell cycle regulation
- **miR-34a**: p53 pathway

These are curated from peer-reviewed GBM literature.
