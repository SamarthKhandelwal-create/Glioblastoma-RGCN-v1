#!/usr/bin/env python
"""Complete pipeline: graph construction, training, and analysis."""

import os
import sys
import subprocess
from pathlib import Path


def run_command(cmd, description):
    print(f"\n{'='*60}\n{description}\n{'='*60}")
    result = subprocess.run(cmd, shell=True, cwd=os.getcwd())
    if result.returncode != 0:
        print(f"FAILED: {description}")
        return False
    return True


def main():
    repo_root = Path(__file__).parent
    os.chdir(repo_root)
    
    required_files = [
        "results/DE_mRNA_logFC_abs_gt1.tsv",
        "results/DE_lncRNA_logFC_abs_gt1.tsv",
        "results/DE_miRNA_logFC_abs_gt1.tsv",
        "results/node_features_matrix.csv",
    ]
    
    missing = [f for f in required_files if not (repo_root / f).exists()]
    if missing:
        print(f"Missing required files: {missing}")
        print("Run preprocessing first.")
        return 1
    
    graph_path = repo_root / "results/hetero_graph_GBM.pt"
    if not graph_path.exists():
        if not run_command("python src/graph/build_graph.py", "Build graph"):
            return 1
    
    cv_dir = repo_root / "results/cv_results"
    cv_dir.mkdir(parents=True, exist_ok=True)
    
    if not run_command(
        f"python src/training/train_cross_validation.py --graph-path {graph_path} --outdir {cv_dir} --folds 10 --epochs 100",
        "Cross-validation"
    ):
        return 1
    
    print(f"\nPipeline complete. Results in: {cv_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
