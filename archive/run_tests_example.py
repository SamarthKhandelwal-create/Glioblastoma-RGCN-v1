#!/usr/bin/env python3
"""
Quick Reference: Running ceRNA Inference Tests
================================================

This script demonstrates the recommended workflow for testing
the ceRNA inference feature with miRNA injection.

Usage:
    python run_tests_example.py

Demonstrates:
    ✓ Setting INJECT_MIRNA=true
    ✓ Running full pipeline
    ✓ Running test suite
    ✓ Collecting results
"""

import subprocess
import os
import sys
from pathlib import Path


def run_command(cmd, description):
    """Run shell command with error handling."""
    print(f"\n{'='*70}")
    print(f"  {description}")
    print(f"{'='*70}")
    print(f"Command: {cmd}\n")
    
    result = subprocess.run(cmd, shell=True)
    return result.returncode == 0


def main():
    print("""
╔══════════════════════════════════════════════════════════════════════╗
║             ceRNA INFERENCE TEST WORKFLOW EXAMPLE                    ║
╚══════════════════════════════════════════════════════════════════════╝
    """)
    
    # Step 1: Set environment
    print("Step 1: Configuring environment...")
    os.environ['INJECT_MIRNA'] = 'true'
    print(f"  ✓ INJECT_MIRNA = {os.environ['INJECT_MIRNA']}")
    
    # Step 2: Run preprocessing
    success = run_command(
        'python src/preprocessing/build_node_features.py',
        'STEP 2: Build Node Features'
    )
    if not success:
        print("✗ Preprocessing failed")
        return 1
    
    # Step 3: Run graph building with ceRNA inference
    success = run_command(
        'python src/graph/build_graph.py',
        'STEP 3: Build Graph with ceRNA Inference'
    )
    if not success:
        print("✗ Graph building failed")
        return 1
    
    # Step 4: Run comprehensive tests
    success = run_command(
        'python test_cerna_inference.py --inject --verbose',
        'STEP 4: Run Test Suite'
    )
    
    # Step 5: Summary
    print(f"\n{'='*70}")
    print("  WORKFLOW COMPLETE")
    print(f"{'='*70}")
    
    if success:
        print("""
✓ Pipeline executed successfully!

Next steps:
  1. Review test results above
  2. Check edge statistics in results/hetero_graph_GBM.pt
  3. Use enhanced graph for training:
  
     python src/training/train_model.py \\
         --graph-path results/hetero_graph_GBM.pt \\
         --epochs 100
  
  4. Run cross-validation:
  
     python src/training/train_cross_validation.py \\
         --graph-path results/hetero_graph_GBM.pt \\
         --folds 10 --epochs 100
        """)
        return 0
    else:
        print("""
✗ One or more steps failed

Troubleshooting:
  1. Check results/ directory for partial outputs
  2. Review error messages above
  3. Ensure all required files exist
  4. Try running steps individually:
  
     $env:INJECT_MIRNA='true'
     python src/preprocessing/build_node_features.py
     python src/graph/build_graph.py
     python test_cerna_inference.py --inject
        """)
        return 1


if __name__ == '__main__':
    sys.exit(main())
