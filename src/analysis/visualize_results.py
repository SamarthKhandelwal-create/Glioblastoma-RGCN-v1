#!/usr/bin/env python3
"""
Simple visualization of 10-fold CV results using JPG format.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.metrics import roc_curve, auc, precision_recall_curve

# Set paths
CV_RESULTS_DIR = Path("results/cv_results_prod")
OUTPUT_DIR = Path("results/analysis_prod")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"Loading CV results from {CV_RESULTS_DIR}...")

# Load data
fold_results = pd.read_csv(CV_RESULTS_DIR / "fold_results.csv")
scores = np.load(CV_RESULTS_DIR / "all_fold_scores.npy")
labels = np.load(CV_RESULTS_DIR / "all_fold_labels.npy")

print(f"✓ Loaded {len(fold_results)} fold results")
print(f"✓ Scores shape: {scores.shape}, Labels shape: {labels.shape}")

# === Plot 1: ROC/PR Curves ===
print("\nGenerating ROC/PR curves...")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# ROC Curve
fpr, tpr, _ = roc_curve(labels, scores)
roc_auc = auc(fpr, tpr)

axes[0].plot(fpr, tpr, color='darkorange', lw=2.5, label=f'ROC curve (AUC = {roc_auc:.4f})')
axes[0].plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label='Random Classifier')
axes[0].fill_between(fpr, tpr, alpha=0.1, color='darkorange')
axes[0].set_xlim([0.0, 1.0])
axes[0].set_ylim([0.0, 1.05])
axes[0].set_xlabel('False Positive Rate', fontsize=11)
axes[0].set_ylabel('True Positive Rate', fontsize=11)
axes[0].set_title('ROC Curve (All CV Folds)', fontsize=12, fontweight='bold')
axes[0].legend(loc="lower right", fontsize=10)
axes[0].grid(alpha=0.3)

# PR Curve
precision, recall, _ = precision_recall_curve(labels, scores)
pr_auc = auc(recall, precision)

axes[1].plot(recall, precision, color='green', lw=2.5, label=f'PR curve (AUC = {pr_auc:.4f})')
axes[1].axhline(y=np.mean(labels), color='gray', linestyle='--', lw=2, label=f'Baseline ({np.mean(labels):.3f})')
axes[1].fill_between(recall, precision, alpha=0.1, color='green')
axes[1].set_xlim([0.0, 1.0])
axes[1].set_ylim([0.0, 1.05])
axes[1].set_xlabel('Recall', fontsize=11)
axes[1].set_ylabel('Precision', fontsize=11)
axes[1].set_title('Precision-Recall Curve (All CV Folds)', fontsize=12, fontweight='bold')
axes[1].legend(loc="upper right", fontsize=10)
axes[1].grid(alpha=0.3)

fig.tight_layout()
roc_pr_path = OUTPUT_DIR / "roc_pr_curves.jpg"
plt.savefig(str(roc_pr_path), dpi=150, format='jpg')
print(f"✓ ROC/PR curves saved to {roc_pr_path}")
plt.close()

# === Plot 2: Per-Fold Metrics ===
print("Generating per-fold metrics plot...")
metrics_to_plot = ['roc_auc', 'pr_auc', 'f1', 'accuracy', 'precision', 'recall']
available_metrics = [m for m in metrics_to_plot if m in fold_results.columns]

fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

for idx, metric in enumerate(available_metrics):
    ax = axes[idx]
    values = fold_results[metric].values
    fold_nums = np.arange(1, len(values) + 1)
    
    # Bar plot
    colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(values)))
    ax.bar(fold_nums, values, alpha=0.8, color=colors, edgecolor='black', linewidth=1.5)
    
    # Mean and std bands
    mean_val = values.mean()
    std_val = values.std()
    ax.axhline(y=mean_val, color='red', linestyle='-', lw=2.5, label=f'Mean: {mean_val:.4f}')
    ax.fill_between(fold_nums, mean_val - std_val, mean_val + std_val, 
                     alpha=0.15, color='red', label=f'±1σ: {std_val:.4f}')
    
    ax.set_xlabel('Fold', fontsize=10, fontweight='bold')
    ax.set_ylabel(metric.upper(), fontsize=10, fontweight='bold')
    ax.set_title(f'{metric.upper()} Across Folds', fontsize=11, fontweight='bold')
    ax.set_xticks(fold_nums)
    ax.legend(fontsize=9)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_ylim([max(0, min(values) - 0.1), min(1.0, max(values) + 0.1)])

# Remove extra subplots
for idx in range(len(available_metrics), len(axes)):
    fig.delaxes(axes[idx])

fig.tight_layout()
metrics_path = OUTPUT_DIR / "cv_metrics.jpg"
plt.savefig(str(metrics_path), dpi=150, format='jpg')
print(f"✓ CV metrics plot saved to {metrics_path}")
plt.close()

# === Generate Summary Report ===
print("\nGenerating analysis report...")
report_path = OUTPUT_DIR / "analysis_report.txt"

with open(report_path, 'w', encoding='utf-8') as f:
    f.write("=" * 80 + "\n")
    f.write("GBM LINK PREDICTION - 10-FOLD CROSS-VALIDATION ANALYSIS\n")
    f.write("=" * 80 + "\n\n")
    
    f.write("CROSS-VALIDATION SUMMARY\n")
    f.write("-" * 80 + "\n")
    
    for metric in available_metrics:
        values = fold_results[metric].values
        f.write(f"{metric:15s}: {values.mean():.4f} ± {values.std():.4f}\n")
    
    f.write("\n")
    f.write("COMBINED METRICS (All Folds)\n")
    f.write("-" * 80 + "\n")
    f.write(f"ROC-AUC: {roc_auc:.4f}\n")
    f.write(f"PR-AUC:  {pr_auc:.4f}\n")
    
    f.write("\n")
    f.write("PER-FOLD RESULTS\n")
    f.write("-" * 80 + "\n")
    f.write(fold_results.to_string())
    
    f.write("\n\n")
    f.write("MODEL PERFORMANCE INTERPRETATION\n")
    f.write("-" * 80 + "\n")
    
    if roc_auc >= 0.9:
        perf = "EXCELLENT - Near-perfect discrimination"
    elif roc_auc >= 0.8:
        perf = "VERY GOOD - Excellent discrimination ability"
    elif roc_auc >= 0.7:
        perf = "GOOD - Fair discrimination ability"
    elif roc_auc >= 0.6:
        perf = "FAIR - Poor discrimination ability"
    else:
        perf = "POOR - Very limited discrimination"
    
    f.write(f"ROC-AUC Performance: {perf}\n")
    f.write(f"  - The model can discriminate between true and false interactions\n")
    f.write(f"  - Higher values indicate better ability to rank interactions\n\n")
    
    f.write(f"PR-AUC Performance: {pr_auc:.4f}\n")
    f.write(f"  - Precision-Recall AUC indicates ability to balance precision and recall\n")
    f.write(f"  - Important for imbalanced prediction tasks\n\n")
    
    recall_vals = fold_results['recall'].values
    precision_vals = fold_results['precision'].values
    f.write(f"Recall (Average): {recall_vals.mean():.4f}\n")
    f.write(f"  - The model identifies ~{recall_vals.mean()*100:.1f}% of true interactions\n\n")
    
    f.write(f"Precision (Average): {precision_vals.mean():.4f}\n")
    f.write(f"  - Of predicted interactions, ~{precision_vals.mean()*100:.1f}% are correct\n\n")
    
    f.write("RECOMMENDATIONS\n")
    f.write("-" * 80 + "\n")
    
    if roc_auc < 0.75:
        f.write("[*] Model performance is moderate. Consider:\n")
        f.write("    - Increasing network size (more interactions needed)\n")
        f.write("    - Tuning hyperparameters (hidden_channels, learning_rate)\n")
        f.write("    - Adding more node features\n")
    elif roc_auc < 0.85:
        f.write("[*] Good model performance. Consider:\n")
        f.write("    - Fine-tuning hyperparameters for marginal improvements\n")
        f.write("    - Ensemble methods combining multiple models\n")
    else:
        f.write("[+] Excellent model performance!\n")
        f.write("    - Model is ready for biological validation\n")
        f.write("    - Consider predicting novel interactions\n")
    
    f.write("\nVALIDATION\n")
    f.write("-" * 80 + "\n")
    f.write("Recommended next steps:\n")
    f.write("1. Validate top predicted interactions experimentally\n")
    f.write("2. Compare with known databases (miRTarBase, ENCORI, etc.)\n")
    f.write("3. Perform biological enrichment analysis on predicted targets\n")
    f.write("4. Analyze predictions for biological coherence\n")
    
    f.write("\n")
    f.write("=" * 80 + "\n")
    f.write("Analysis Complete\n")
    f.write("=" * 80 + "\n")

print(f"✓ Analysis report saved to {report_path}")

print("\n" + "=" * 80)
print("VISUALIZATION GENERATION COMPLETE")
print("=" * 80)
print(f"\nGenerated files:")
print(f"  - {roc_pr_path}")
print(f"  - {metrics_path}")
print(f"  - {report_path}")
