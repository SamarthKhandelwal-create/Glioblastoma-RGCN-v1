#!/usr/bin/env python
"""
Complete pipeline for GBM link prediction model:
1. Build heterogeneous graph from interactions
2. Run 10-fold cross-validation
3. Extract top genes and perform enrichment
4. Generate comprehensive visualizations
"""

import os
import sys
import subprocess
from pathlib import Path

def run_command(cmd, description):
    """Run shell command and report status."""
    print(f"\n{'='*80}")
    print(f"STEP: {description}")
    print(f"{'='*80}")
    print(f"Command: {cmd}\n")
    
    result = subprocess.run(cmd, shell=True, cwd=os.getcwd())
    
    if result.returncode != 0:
        print(f"\n❌ FAILED: {description}")
        return False
    print(f"\n✓ SUCCESS: {description}")
    return True

def main():
    repo_root = Path(__file__).parent
    os.chdir(repo_root)
    
    print(f"Starting GBM Link Prediction Pipeline")
    print(f"Repository root: {repo_root}")
    
    # Check Python version
    if sys.version_info < (3, 8):
        print(f"❌ Python 3.8+ required, got {sys.version}")
        return 1
    
    # Step 1: Check for required files
    print("\n" + "="*80)
    print("CHECKING PREREQUISITES")
    print("="*80)
    
    required_de_files = [
        "results/DE_mRNA_logFC_abs_gt1.tsv",
        "results/DE_lncRNA_logFC_abs_gt1.tsv",
        "results/DE_miRNA_logFC_abs_gt1.tsv",
        "results/node_features_matrix.csv",
    ]
    
    missing = []
    for f in required_de_files:
        if not (repo_root / f).exists():
            missing.append(f)
            print(f"  ❌ Missing: {f}")
        else:
            print(f"  ✓ Found: {f}")
    
    if missing:
        print(f"\n❌ Missing {len(missing)} required files. Cannot proceed.")
        print("   Please run Phase 1 (differential expression and feature matrix generation) first.")
        return 1
    
    # Step 2: Build graph
    graph_path = repo_root / "results/hetero_graph_GBM.pt"
    if not graph_path.exists():
        print(f"\n  Graph not found. Building...")
        if not run_command(
            f"python src/graph/build_graph.py",
            "Build heterogeneous graph"
        ):
            return 1
    else:
        print(f"\n✓ Graph already exists at {graph_path}")
    
    # Step 3: Run cross-validation
    cv_results_dir = repo_root / "results/cv_results"
    cv_results_dir.mkdir(parents=True, exist_ok=True)
    
    if not run_command(
        f"python src/training/train_cross_validation.py --graph-path {graph_path} --outdir {cv_results_dir} --folds 10 --epochs 100 --seed 42",
        "Run 10-fold cross-validation"
    ):
        return 1
    
    # Step 4: Generate analysis and visualizations
    analysis_dir = repo_root / "results/analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    
    if not run_command(
        f"python src/analysis/analyze_cv_results.py --cv-results-dir {cv_results_dir} --outdir {analysis_dir}",
        "Generate visualizations and enrichment analysis"
    ):
        return 1
    
    # Summary
    print(f"\n\n{'='*80}")
    print("PIPELINE COMPLETE ✓")
    print(f"{'='*80}")
    print(f"\nKey outputs:")
    print(f"  📊 Cross-validation results:  {cv_results_dir}/fold_results.csv")
    print(f"  📈 ROC/PR curves:            {analysis_dir}/roc_pr_curves.png")
    print(f"  📋 Metrics plot:             {analysis_dir}/cv_metrics.png")
    print(f"  🔬 Enrichment analysis:      {analysis_dir}/enrichment_results.png")
    print(f"  📝 Summary report:           {analysis_dir}/analysis_report.txt")
    print(f"\nView results in: {analysis_dir}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
