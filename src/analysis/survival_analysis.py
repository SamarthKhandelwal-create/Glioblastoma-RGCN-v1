#!/usr/bin/env python3
"""
Survival Analysis for Top Biomarker Genes

Performs Kaplan-Meier survival analysis and Cox proportional hazards regression
on the top 3 genes (by model confidence score) from the biomarker discovery pipeline.

Supports two modes:
  1. Online: Fetches data from TCGA/Xena hubs (requires network)
  2. Offline: Uses cached/synthetic data for demonstration

Output:
  - results/survival/km_curves.png          : Kaplan-Meier curves for each gene
  - results/survival/cox_results.csv        : Cox regression coefficients and p-values
  - results/survival/survival_report.txt    : Full analysis report
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
import argparse

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# --- Configuration ---
RESULTS_DIR = Path(os.environ.get("RESULTS_DIR", "results"))
BIOMARKERS_FILE = RESULTS_DIR / "biomarkers" / "top_biomarkers.csv"
OUTPUT_DIR = RESULTS_DIR / "survival"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Cached data paths
CACHED_CLINICAL_FILE = RESULTS_DIR / "survival" / "cached_clinical.csv"
CACHED_EXPRESSION_FILE = RESULTS_DIR / "survival" / "cached_expression.csv"

# Xena hub settings for clinical + expression data
TOIL_HOST = "https://toil.xenahubs.net"
TOIL_COUNTS = "TcgaTargetGtex_gene_expected_count"
LEGACY_HOST = "https://tcga.xenahubs.net"
CLINICAL_DATASET = "TCGA.GBM.sampleMap/GBM_clinicalMatrix"

# Number of top genes to analyze (default: 3)
TOP_N_GENES = 3


def load_top_genes(biomarkers_path: Path, top_n: int = 3) -> pd.DataFrame:
    """Load top N mRNA genes from biomarker results (skip miRNAs for survival)."""
    if not biomarkers_path.exists():
        print(f"Error: Biomarkers file not found: {biomarkers_path}")
        print("Run find_biomarkers.py first to generate top_biomarkers.csv")
        sys.exit(1)
    
    df = pd.read_csv(biomarkers_path)
    
    # Filter to mRNA only (survival analysis typically uses protein-coding genes)
    mrna_df = df[df['type'] == 'mRNA'].copy()
    
    if len(mrna_df) < top_n:
        print(f"Warning: Only {len(mrna_df)} mRNA genes available, using all of them")
        top_n = len(mrna_df)
    
    top_genes = mrna_df.head(top_n)
    print(f"Top {top_n} mRNA genes by model confidence:")
    for i, (_, row) in enumerate(top_genes.iterrows(), 1):
        print(f"  {i}. {row['symbol']} ({row['gene_id']}) | Log2FC: {row['log2FC']:.3f}")
    
    return top_genes


def generate_synthetic_data(gene_ids: list, gene_symbols: list, n_samples: int = 150) -> tuple:
    """
    Generate synthetic survival and expression data for demonstration.
    
    Creates realistic GBM-like survival distributions with gene expression
    correlated to survival outcomes based on the log2FC direction.
    """
    print("\n[OFFLINE MODE] Generating synthetic TCGA-GBM-like data...")
    
    np.random.seed(42)
    
    # Generate sample IDs
    samples = [f"TCGA-GBM-{i:04d}" for i in range(n_samples)]
    
    # Generate survival times (Weibull distribution, median ~400 days for GBM)
    # GBM has poor prognosis: median OS ~14-16 months
    shape, scale = 1.5, 450
    survival_times = np.random.weibull(shape, n_samples) * scale
    survival_times = np.clip(survival_times, 30, 2000)  # Clip to realistic range
    
    # Generate censoring (about 40% censored in typical GBM cohort)
    max_followup = 1500
    censoring_times = np.random.uniform(200, max_followup, n_samples)
    
    # Observed time is min of survival and censoring
    observed_times = np.minimum(survival_times, censoring_times)
    events = (survival_times <= censoring_times).astype(int)
    
    # Clinical DataFrame
    clinical_df = pd.DataFrame({
        'sample': samples,
        'time': observed_times,
        'event': events
    })
    
    print(f"  Generated {n_samples} synthetic samples")
    print(f"  Events (deaths): {events.sum()} / {n_samples}")
    print(f"  Median follow-up: {np.median(observed_times):.0f} days")
    
    # Generate expression data correlated with survival
    # Higher expression of upregulated genes → worse survival (for GBM oncogenes)
    expr_data = {}
    for gene_id, gene_symbol in zip(gene_ids, gene_symbols):
        # Base expression (log2 scale, centered around 8)
        base_expr = np.random.normal(8, 2, n_samples)
        
        # Add correlation with survival (negative = high expr → short survival)
        # This simulates oncogene behavior in GBM
        survival_effect = -0.3 * (survival_times - np.mean(survival_times)) / np.std(survival_times)
        expr_data[gene_id] = base_expr + survival_effect + np.random.normal(0, 0.5, n_samples)
    
    expr_df = pd.DataFrame(expr_data, index=samples)
    
    print(f"  Generated expression for {len(gene_ids)} genes")
    
    return clinical_df, expr_df


def fetch_clinical_data(offline: bool = False) -> pd.DataFrame:
    """Fetch clinical data from TCGA GBM cohort via Xena or use cached data."""
    
    # Check for cached data first
    if CACHED_CLINICAL_FILE.exists() and not offline:
        print("\nLoading cached clinical data...")
        df = pd.read_csv(CACHED_CLINICAL_FILE)
        print(f"  Loaded {len(df)} samples from cache")
        return df
    
    if offline:
        return None  # Will use synthetic data
    
    try:
        import xenaPython as xp
    except ImportError:
        print("xenaPython not installed. Using offline mode.")
        return None
    
    print("\nFetching clinical data from TCGA-GBM...")
    
    try:
        # Get all samples
        samples = xp.dataset_samples(LEGACY_HOST, CLINICAL_DATASET, None)
        
        # Fields needed for survival analysis
        fields = ['sampleID', 'days_to_death', 'days_to_last_followup', 'vital_status']
        
        clinical_data = {}
        for field in fields:
            try:
                values = xp.dataset_probe_values(LEGACY_HOST, CLINICAL_DATASET, samples, [field])
                if values and len(values) > 0:
                    clinical_data[field] = values[0]
            except Exception as e:
                print(f"  Warning: Could not fetch {field}: {e}")
        
        # Build DataFrame
        df = pd.DataFrame({'sample': samples})
        for field in fields:
            if field in clinical_data:
                df[field] = clinical_data[field]
        
        # Calculate survival time and event
        df['time'] = df.apply(lambda row: _calculate_survival_time(row), axis=1)
        df['event'] = df['vital_status'].apply(lambda x: 1 if str(x).upper() == 'DECEASED' else 0)
        
        # Filter valid samples
        df = df.dropna(subset=['time'])
        df = df[df['time'] > 0]
        
        print(f"  Loaded {len(df)} samples with valid survival data")
        print(f"  Events (deaths): {df['event'].sum()} / {len(df)}")
        
        # Cache for future use
        if len(df) > 0:
            df.to_csv(CACHED_CLINICAL_FILE, index=False)
            print(f"  Cached to {CACHED_CLINICAL_FILE}")
        
        return df
    
    except Exception as e:
        print(f"  Error fetching from Xena: {e}")
        print("  Falling back to offline mode...")
        return None


def _calculate_survival_time(row) -> float:
    """Calculate survival time from clinical row."""
    try:
        if pd.notna(row.get('days_to_death')) and float(row['days_to_death']) > 0:
            return float(row['days_to_death'])
        elif pd.notna(row.get('days_to_last_followup')) and float(row['days_to_last_followup']) > 0:
            return float(row['days_to_last_followup'])
    except (ValueError, TypeError):
        pass
    return np.nan


def fetch_gene_expression(samples: list, gene_ids: list, offline: bool = False) -> pd.DataFrame:
    """Fetch gene expression for specified genes and samples."""
    
    if offline:
        return None  # Will use synthetic data
    
    try:
        import xenaPython as xp
    except ImportError:
        return None
    
    print(f"\nFetching expression data for {len(gene_ids)} genes...")
    
    expr_data = {}
    for gene_id in gene_ids:
        try:
            # Query by Ensembl ID (strip version for matching)
            base_id = gene_id.split('.')[0]
            values = xp.dataset_probe_values(TOIL_HOST, TOIL_COUNTS, samples, [base_id])
            if values and len(values) > 0:
                expr_data[gene_id] = values[0]
                print(f"  ✓ {gene_id}")
            else:
                print(f"  ✗ {gene_id} (no data)")
        except Exception as e:
            print(f"  ✗ {gene_id}: {e}")
    
    df = pd.DataFrame(expr_data, index=samples)
    return df


def run_kaplan_meier(clinical_df: pd.DataFrame, gene_col: str, gene_symbol: str, ax=None):
    """
    Run Kaplan-Meier analysis for a single gene.
    
    Stratifies patients into high/low expression groups by median split
    and performs log-rank test.
    """
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test
    
    # Median split
    median_expr = clinical_df[gene_col].median()
    high_mask = clinical_df[gene_col] >= median_expr
    low_mask = ~high_mask
    
    # KM fitters
    kmf_high = KaplanMeierFitter()
    kmf_low = KaplanMeierFitter()
    
    high_df = clinical_df[high_mask]
    low_df = clinical_df[low_mask]
    
    kmf_high.fit(high_df['time'], high_df['event'], label=f'{gene_symbol} High (n={len(high_df)})')
    kmf_low.fit(low_df['time'], low_df['event'], label=f'{gene_symbol} Low (n={len(low_df)})')
    
    # Log-rank test
    results = logrank_test(high_df['time'], low_df['time'], high_df['event'], low_df['event'])
    p_value = results.p_value
    
    # Plot
    if ax is not None:
        kmf_high.plot_survival_function(ax=ax, ci_show=True, color='red')
        kmf_low.plot_survival_function(ax=ax, ci_show=True, color='blue')
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('Survival Probability')
        ax.set_title(f'{gene_symbol}\nLog-rank p = {p_value:.4f}')
        ax.legend(loc='lower left')
    
    return {
        'gene': gene_symbol,
        'n_high': len(high_df),
        'n_low': len(low_df),
        'median_high': kmf_high.median_survival_time_,
        'median_low': kmf_low.median_survival_time_,
        'logrank_p': p_value
    }


def run_cox_regression(clinical_df: pd.DataFrame, gene_cols: list, gene_symbols: list) -> pd.DataFrame:
    """
    Run Cox proportional hazards regression with all genes as covariates.
    
    Returns DataFrame with hazard ratios, confidence intervals, and p-values.
    """
    from lifelines import CoxPHFitter
    
    print("\nRunning Cox Proportional Hazards regression...")
    
    # Prepare data for Cox model
    cox_df = clinical_df[['time', 'event'] + gene_cols].copy()
    
    # Rename columns to gene symbols for readability
    rename_map = {gene_cols[i]: gene_symbols[i] for i in range(len(gene_cols))}
    cox_df = cox_df.rename(columns=rename_map)
    
    # Z-score normalize expression values
    for symbol in gene_symbols:
        cox_df[symbol] = (cox_df[symbol] - cox_df[symbol].mean()) / cox_df[symbol].std()
    
    # Fit Cox model
    cph = CoxPHFitter()
    cph.fit(cox_df, duration_col='time', event_col='event')
    
    # Extract results
    summary = cph.summary.copy()
    summary['HR'] = np.exp(summary['coef'])
    summary['HR_lower'] = np.exp(summary['coef'] - 1.96 * summary['se(coef)'])
    summary['HR_upper'] = np.exp(summary['coef'] + 1.96 * summary['se(coef)'])
    
    print("\nCox Regression Results:")
    print("-" * 70)
    print(f"{'Gene':<15} {'HR':>10} {'95% CI':>20} {'p-value':>12}")
    print("-" * 70)
    for gene in gene_symbols:
        row = summary.loc[gene]
        ci = f"({row['HR_lower']:.3f} - {row['HR_upper']:.3f})"
        print(f"{gene:<15} {row['HR']:>10.3f} {ci:>20} {row['p']:>12.4f}")
    
    return summary


def generate_report(km_results: list, cox_summary: pd.DataFrame, output_path: Path):
    """Generate comprehensive survival analysis report."""
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("=" * 80 + "\n")
        f.write("GBM SURVIVAL ANALYSIS - TOP MODEL-IDENTIFIED GENES\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("OVERVIEW\n")
        f.write("-" * 80 + "\n")
        f.write("This analysis examines the prognostic value of the top genes identified\n")
        f.write("by the RGCN link prediction model, using TCGA-GBM patient data.\n\n")
        
        f.write("KAPLAN-MEIER ANALYSIS\n")
        f.write("-" * 80 + "\n")
        f.write("Patients stratified by median expression (high vs. low)\n\n")
        
        f.write(f"{'Gene':<15} {'N_High':>8} {'N_Low':>8} {'Med_High':>12} {'Med_Low':>12} {'LogRank_p':>12}\n")
        f.write("-" * 80 + "\n")
        
        for result in km_results:
            med_high = f"{result['median_high']:.1f}" if not np.isinf(result['median_high']) else "NR"
            med_low = f"{result['median_low']:.1f}" if not np.isinf(result['median_low']) else "NR"
            f.write(f"{result['gene']:<15} {result['n_high']:>8} {result['n_low']:>8} "
                    f"{med_high:>12} {med_low:>12} {result['logrank_p']:>12.4f}\n")
        
        f.write("\n* NR = Not Reached (median survival not reached in follow-up period)\n")
        
        f.write("\n\nCOX PROPORTIONAL HAZARDS REGRESSION\n")
        f.write("-" * 80 + "\n")
        f.write("Multivariate analysis with all genes as covariates (z-score normalized)\n\n")
        
        f.write(f"{'Gene':<15} {'Hazard Ratio':>15} {'95% CI':>25} {'p-value':>12}\n")
        f.write("-" * 80 + "\n")
        
        for gene in cox_summary.index:
            row = cox_summary.loc[gene]
            ci = f"({row['HR_lower']:.3f} - {row['HR_upper']:.3f})"
            sig = "*" if row['p'] < 0.05 else ""
            f.write(f"{gene:<15} {row['HR']:>15.3f} {ci:>25} {row['p']:>11.4f}{sig}\n")
        
        f.write("\n* p < 0.05\n")
        
        f.write("\n\nINTERPRETATION\n")
        f.write("-" * 80 + "\n")
        f.write("Hazard Ratio (HR) interpretation:\n")
        f.write("  - HR > 1: Higher expression associated with WORSE survival (risk factor)\n")
        f.write("  - HR < 1: Higher expression associated with BETTER survival (protective)\n")
        f.write("  - HR = 1: No association with survival\n\n")
        
        # Identify significant genes
        sig_genes = cox_summary[cox_summary['p'] < 0.05]
        if len(sig_genes) > 0:
            f.write("SIGNIFICANT FINDINGS:\n")
            for gene in sig_genes.index:
                row = sig_genes.loc[gene]
                direction = "poor prognosis" if row['HR'] > 1 else "better prognosis"
                f.write(f"  - {gene}: HR = {row['HR']:.3f}, p = {row['p']:.4f} ({direction})\n")
        else:
            f.write("No genes reached statistical significance (p < 0.05).\n")
            f.write("This may be due to small sample size or lack of true association.\n")
        
        f.write("\n\nFILES GENERATED\n")
        f.write("-" * 80 + "\n")
        f.write(f"  - {OUTPUT_DIR / 'km_curves.png'} : Kaplan-Meier survival curves\n")
        f.write(f"  - {OUTPUT_DIR / 'cox_results.csv'} : Cox regression results\n")
        f.write(f"  - {OUTPUT_DIR / 'survival_report.txt'} : This report\n")
    
    print(f"\n✓ Report saved to {output_path}")


def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Survival analysis on top model-identified genes")
    parser.add_argument('--offline', action='store_true', 
                        help='Use synthetic data (no network required)')
    parser.add_argument('--top-n', type=int, default=TOP_N_GENES,
                        help=f'Number of top genes to analyze (default: {TOP_N_GENES})')
    args = parser.parse_args()
    
    print("=" * 80)
    print("GBM SURVIVAL ANALYSIS")
    print("Kaplan-Meier and Cox Regression on Top Model-Identified Genes")
    print("=" * 80)
    
    # Check for required packages
    try:
        import lifelines
        import matplotlib.pyplot as plt
    except ImportError as e:
        print(f"\nError: Missing required package: {e}")
        print("Install with: pip install lifelines matplotlib")
        sys.exit(1)
    
    # Load top genes
    top_genes = load_top_genes(BIOMARKERS_FILE, args.top_n)
    gene_ids = top_genes['gene_id'].tolist()
    gene_symbols = top_genes['symbol'].tolist()
    
    # Try to fetch real data, fall back to synthetic if needed
    clinical_df = fetch_clinical_data(offline=args.offline)
    use_synthetic = clinical_df is None or len(clinical_df) == 0
    
    if use_synthetic:
        # Generate synthetic data
        clinical_df, expr_df = generate_synthetic_data(gene_ids, gene_symbols)
        merged_df = clinical_df.set_index('sample').join(expr_df, how='inner')
    else:
        # Fetch expression data
        expr_df = fetch_gene_expression(clinical_df['sample'].tolist(), gene_ids, offline=args.offline)
        
        if expr_df is None or len(expr_df.columns) == 0:
            print("\nCould not fetch expression data. Using synthetic data...")
            clinical_df, expr_df = generate_synthetic_data(gene_ids, gene_symbols)
            merged_df = clinical_df.set_index('sample').join(expr_df, how='inner')
            use_synthetic = True
        else:
            # Merge clinical and expression
            merged_df = clinical_df.set_index('sample').join(expr_df, how='inner')
            merged_df = merged_df.dropna(subset=gene_ids)
    
    if use_synthetic:
        print("\n⚠️  Using SYNTHETIC data for demonstration")
        print("   Run with network access for real TCGA-GBM analysis")
    
    print(f"\n{len(merged_df)} samples with complete clinical + expression data")
    
    if len(merged_df) < 20:
        print("Warning: Very few samples with complete data. Results may not be reliable.")
    
    # --- Kaplan-Meier Analysis ---
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(1, len(gene_ids), figsize=(5 * len(gene_ids), 4))
    if len(gene_ids) == 1:
        axes = [axes]
    
    km_results = []
    for i, (gene_id, gene_symbol) in enumerate(zip(gene_ids, gene_symbols)):
        result = run_kaplan_meier(merged_df, gene_id, gene_symbol, ax=axes[i])
        km_results.append(result)
    
    plt.tight_layout()
    km_path = OUTPUT_DIR / "km_curves.png"
    plt.savefig(km_path, dpi=150, bbox_inches='tight')
    print(f"\n✓ Kaplan-Meier curves saved to {km_path}")
    plt.close()
    
    # --- Cox Regression ---
    cox_summary = run_cox_regression(merged_df, gene_ids, gene_symbols)
    
    # Save Cox results
    cox_path = OUTPUT_DIR / "cox_results.csv"
    cox_summary.to_csv(cox_path)
    print(f"✓ Cox results saved to {cox_path}")
    
    # --- Generate Report ---
    report_path = OUTPUT_DIR / "survival_report.txt"
    generate_report(km_results, cox_summary, report_path)
    
    print("\n" + "=" * 80)
    print("SURVIVAL ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"\nOutput files in: {OUTPUT_DIR}/")
    print(f"  1. km_curves.png       - Kaplan-Meier survival curves")
    print(f"  2. cox_results.csv     - Cox regression coefficients")
    print(f"  3. survival_report.txt - Full analysis report")
    
    if use_synthetic:
        print("\n⚠️  Note: Results based on synthetic data for demonstration purposes")


if __name__ == "__main__":
    main()
