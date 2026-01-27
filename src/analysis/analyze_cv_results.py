"""
Post-training analysis and visualization for GBM link prediction model.

Generates:
1. ROC/PR curves from CV fold results
2. Top genes extraction
3. GO and KEGG enrichment analysis
4. Performance summary visualizations
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import auc, precision_recall_curve, roc_curve

try:
    from gseapy import enrichr
    GSEAPY_AVAILABLE = True
except ImportError:
    GSEAPY_AVAILABLE = False
    print("Warning: gseapy not installed. GO/KEGG enrichment will be skipped.")


def load_cv_results(cv_results_dir):
    """Load fold results and scores from CV."""
    cv_results_dir = Path(cv_results_dir)
    
    # Load fold results CSV
    fold_results_df = pd.read_csv(cv_results_dir / "fold_results.csv")
    
    # Load scores and labels
    scores = np.load(cv_results_dir / "all_fold_scores.npy")
    labels = np.load(cv_results_dir / "all_fold_labels.npy")
    
    return fold_results_df, scores, labels


def plot_roc_pr_curves(scores, labels, outdir):
    """Generate ROC and PR curves from all CV fold scores."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # ROC Curve
    fpr, tpr, _ = roc_curve(labels, scores)
    roc_auc = auc(fpr, tpr)
    
    axes[0].plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {roc_auc:.3f})')
    axes[0].plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label='Random Classifier')
    axes[0].set_xlim([0.0, 1.0])
    axes[0].set_ylim([0.0, 1.05])
    axes[0].set_xlabel('False Positive Rate')
    axes[0].set_ylabel('True Positive Rate')
    axes[0].set_title('ROC Curve (All CV Folds)')
    axes[0].legend(loc="lower right")
    axes[0].grid(alpha=0.3)
    
    # PR Curve
    precision, recall, _ = precision_recall_curve(labels, scores)
    pr_auc = auc(recall, precision)
    
    axes[1].plot(recall, precision, color='green', lw=2, label=f'PR curve (AUC = {pr_auc:.3f})')
    axes[1].axhline(y=np.mean(labels), color='gray', linestyle='--', label=f'Baseline ({np.mean(labels):.3f})')
    axes[1].set_xlim([0.0, 1.0])
    axes[1].set_ylim([0.0, 1.05])
    axes[1].set_xlabel('Recall')
    axes[1].set_ylabel('Precision')
    axes[1].set_title('Precision-Recall Curve (All CV Folds)')
    axes[1].legend(loc="upper right")
    axes[1].grid(alpha=0.3)
    
    fig.tight_layout()
    roc_pr_path = outdir / "roc_pr_curves.png"
    plt.savefig(roc_pr_path, dpi=300, bbox_inches='tight')
    print(f"✓ ROC/PR curves saved to {roc_pr_path}")
    plt.close()
    
    return {'roc_auc': roc_auc, 'pr_auc': pr_auc}


def plot_cv_metrics(fold_results_df, outdir):
    """Plot performance metrics across CV folds."""
    outdir = Path(outdir)
    
    metrics_to_plot = ['roc_auc', 'pr_auc', 'f1', 'accuracy', 'precision', 'recall']
    available_metrics = [m for m in metrics_to_plot if m in fold_results_df.columns]
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    for idx, metric in enumerate(available_metrics):
        ax = axes[idx]
        values = fold_results_df[metric].values
        
        ax.bar(range(len(values)), values, alpha=0.7, color='steelblue', edgecolor='black')
        ax.axhline(y=values.mean(), color='red', linestyle='--', lw=2, label=f'Mean: {values.mean():.3f}')
        ax.fill_between(range(len(values)), values.mean() - values.std(), 
                         values.mean() + values.std(), alpha=0.2, color='red', label=f'±1 Std: {values.std():.3f}')
        ax.set_xlabel('Fold')
        ax.set_ylabel(metric.upper())
        ax.set_title(f'{metric.upper()} Across Folds')
        ax.legend()
        ax.grid(axis='y', alpha=0.3)
    
    # Remove extra subplots
    for idx in range(len(available_metrics), len(axes)):
        fig.delaxes(axes[idx])
    
    fig.tight_layout()
    metrics_path = outdir / "cv_metrics.png"
    plt.savefig(metrics_path, dpi=300, bbox_inches='tight')
    print(f"✓ CV metrics plot saved to {metrics_path}")
    plt.close()


def extract_top_genes(top_n=50):
    """Extract top genes from novel predictions or DE genes."""
    results_dir = Path("results")
    
    # Try to load predictions
    pred_path = results_dir / "novel_predictions.csv"
    if pred_path.exists():
        df = pd.read_csv(pred_path)
        top_genes = df.nlargest(top_n, 'Score')['Target'].unique().tolist()
        print(f"✓ Extracted {len(top_genes)} top prediction targets")
        return top_genes
    
    # Fallback: load top DE genes
    mrna_path = results_dir / "DE_mRNA_logFC_abs_gt1.tsv"
    lncrna_path = results_dir / "DE_lncRNA_logFC_abs_gt1.tsv"
    
    top_genes = []
    if mrna_path.exists():
        df_mrna = pd.read_csv(mrna_path, sep='\t')
        top_genes.extend(df_mrna.nlargest(top_n // 2, 'log2FoldChange')['ensembl_id'].tolist())
    
    if lncrna_path.exists():
        df_lncrna = pd.read_csv(lncrna_path, sep='\t')
        top_genes.extend(df_lncrna.nlargest(top_n // 2, 'log2FoldChange')['ensembl_id'].tolist())
    
    print(f"✓ Extracted {len(top_genes)} top DE genes")
    return top_genes[:top_n]


def perform_enrichment(gene_list, outdir):
    """Perform GO and KEGG enrichment analysis."""
    outdir = Path(outdir)
    
    if not GSEAPY_AVAILABLE:
        print("⚠ gseapy not installed. Skipping enrichment analysis.")
        print("  Install with: pip install gseapy")
        return {}
    
    print("Running enrichment analysis...")
    results = {}
    
    # GO Biological Process
    try:
        go_bp = enrichr(gene_list, gene_sets='GO_Biological_Process_2021', outdir=str(outdir / 'go_bp'))
        if go_bp is not None and len(go_bp) > 0:
            go_bp_path = outdir / "go_bp_enrichment.csv"
            go_bp.to_csv(go_bp_path)
            print(f"✓ GO BP enrichment saved to {go_bp_path}")
            results['go_bp'] = go_bp.head(10)
    except Exception as e:
        print(f"⚠ GO BP enrichment failed: {e}")
    
    # KEGG
    try:
        kegg = enrichr(gene_list, gene_sets='KEGG_2019', outdir=str(outdir / 'kegg'))
        if kegg is not None and len(kegg) > 0:
            kegg_path = outdir / "kegg_enrichment.csv"
            kegg.to_csv(kegg_path)
            print(f"✓ KEGG enrichment saved to {kegg_path}")
            results['kegg'] = kegg.head(10)
    except Exception as e:
        print(f"⚠ KEGG enrichment failed: {e}")
    
    return results


def plot_enrichment(enrichment_results, outdir):
    """Visualize top enrichment results."""
    outdir = Path(outdir)
    
    if not enrichment_results:
        return
    
    fig, axes = plt.subplots(len(enrichment_results), 1, figsize=(12, 5 * len(enrichment_results)))
    if len(enrichment_results) == 1:
        axes = [axes]
    
    for idx, (analysis_type, df) in enumerate(enrichment_results.items()):
        if df is None or len(df) == 0:
            continue
        
        ax = axes[idx]
        top_10 = df.head(10).sort_values('Adjusted P-value')
        
        ax.barh(range(len(top_10)), -np.log10(top_10['Adjusted P-value']), alpha=0.7, color='steelblue')
        ax.set_yticks(range(len(top_10)))
        ax.set_yticklabels([t[:50] for t in top_10['Term']], fontsize=9)
        ax.set_xlabel('-log10(Adjusted P-value)')
        ax.set_title(f'{analysis_type.upper()} Enrichment')
        ax.grid(axis='x', alpha=0.3)
    
    fig.tight_layout()
    enrich_path = outdir / "enrichment_results.png"
    plt.savefig(enrich_path, dpi=300, bbox_inches='tight')
    print(f"✓ Enrichment plot saved to {enrich_path}")
    plt.close()


def generate_summary_report(fold_results_df, scores, labels, roc_metrics, outdir):
    """Generate a text summary report."""
    outdir = Path(outdir)
    
    report = []
    report.append("=" * 80)
    report.append("GBM LINK PREDICTION MODEL - CROSS-VALIDATION ANALYSIS REPORT")
    report.append("=" * 80)
    report.append("")
    
    # Overall metrics
    report.append("OVERALL PERFORMANCE (ALL CV FOLDS COMBINED)")
    report.append("-" * 80)
    report.append(f"  ROC-AUC:              {roc_metrics['roc_auc']:.4f}")
    report.append(f"  PR-AUC:               {roc_metrics['pr_auc']:.4f}")
    report.append("")
    
    # Per-metric statistics
    report.append("PER-FOLD METRIC STATISTICS")
    report.append("-" * 80)
    for col in fold_results_df.columns:
        mean_val = fold_results_df[col].mean()
        std_val = fold_results_df[col].std()
        min_val = fold_results_df[col].min()
        max_val = fold_results_df[col].max()
        report.append(f"  {col:15s}: {mean_val:.4f} ± {std_val:.4f} (min: {min_val:.4f}, max: {max_val:.4f})")
    
    report.append("")
    report.append("FILES GENERATED")
    report.append("-" * 80)
    report.append(f"  ROC/PR Curves:        roc_pr_curves.png")
    report.append(f"  CV Metrics Plot:      cv_metrics.png")
    report.append(f"  Fold Results CSV:     fold_results.csv")
    report.append(f"  Enrichment Plot:      enrichment_results.png")
    
    report.append("")
    report.append("RECOMMENDATIONS")
    report.append("-" * 80)
    if roc_metrics['roc_auc'] > 0.85:
        report.append("  [OK] Model shows strong discriminative ability (ROC-AUC > 0.85)")
    elif roc_metrics['roc_auc'] > 0.75:
        report.append("  [OK] Model shows fair discriminative ability (ROC-AUC > 0.75)")
    else:
        report.append("  [!] Model performance may be limited (ROC-AUC < 0.75)")
        report.append("    Consider: increasing hidden dimensions, more training epochs, data augmentation")
    
    report.append("")
    report.append("=" * 80)
    
    report_text = "\n".join(report)
    report_path = outdir / "analysis_report.txt"
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(report_text)
    
    print(f"\n{report_text}")
    print(f"\n✓ Report saved to {report_path}")


def main():
    parser = argparse.ArgumentParser(description="Analyze CV results and generate visualizations")
    parser.add_argument("--cv-results-dir", type=str, default="results/cv_results",
                        help="Directory containing CV results")
    parser.add_argument("--outdir", type=str, default="results/analysis",
                        help="Output directory for plots and reports")
    
    args = parser.parse_args()
    
    # Create output dir
    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    
    # Load results
    print(f"Loading CV results from {args.cv_results_dir}...")
    fold_results_df, scores, labels = load_cv_results(args.cv_results_dir)
    
    # Generate plots
    print("\nGenerating visualizations...")
    roc_metrics = plot_roc_pr_curves(scores, labels, args.outdir)
    plot_cv_metrics(fold_results_df, args.outdir)
    
    # Extract genes and run enrichment
    print("\nExtracting top genes for enrichment analysis...")
    top_genes = extract_top_genes(top_n=50)
    
    print("Performing GO/KEGG enrichment analysis...")
    enrichment_results = perform_enrichment(top_genes, args.outdir)
    plot_enrichment(enrichment_results, args.outdir)
    
    # Generate report
    print("\nGenerating summary report...")
    generate_summary_report(fold_results_df, scores, labels, roc_metrics, args.outdir)
    
    print("\n✓ Analysis complete!")
    print(f"  All results saved to {args.outdir}")
    
    return 0


if __name__ == "__main__":
    exit(main())
