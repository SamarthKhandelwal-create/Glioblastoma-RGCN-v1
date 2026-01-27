#!/usr/bin/env python3
"""
GO and KEGG Enrichment Analysis for GBM Biomarkers

Performs comprehensive functional enrichment analysis on high-confidence
biomarker candidates identified from the link prediction model.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print("=" * 80)
print("GBM BIOMARKER ENRICHMENT ANALYSIS")
print("=" * 80)

# Paths
RESULTS_DIR = Path("results")
BIOMARKER_DIR = RESULTS_DIR / "biomarkers"
OUTPUT_DIR = BIOMARKER_DIR / "enrichment"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"\nLoading biomarkers...")

# Load biomarkers
biomarker_path = BIOMARKER_DIR / "top_biomarkers.csv"
if not biomarker_path.exists():
    print(f"ERROR: {biomarker_path} not found")
    exit(1)

biomarkers = pd.read_csv(biomarker_path)
print(f"✓ Loaded {len(biomarkers)} biomarkers")

# Filter for mRNA (suitable for GO/KEGG)
mrna_biomarkers = biomarkers[biomarkers['type'] == 'mRNA'].copy()
gene_symbols = mrna_biomarkers['symbol'].tolist()
gene_symbols = [g.strip() for g in gene_symbols if g.strip()]

print(f"✓ Extracted {len(gene_symbols)} mRNA genes for enrichment")

# === Try gseapy first ===
print(f"\nAttempting GO/KEGG enrichment analysis...")

enrichment_results = {}

try:
    from gseapy import enrichr
    
    print("  Using gseapy enrichr...")
    
    # GO Biological Process
    print("  [1/2] GO Biological Process 2021...")
    try:
        go_bp = enrichr(gene_symbols, 
                       gene_sets=['GO_Biological_Process_2021'],
                       outdir=str(OUTPUT_DIR / 'go_bp'))
        if go_bp is not None:
            # Save results
            go_bp_path = OUTPUT_DIR / "GO_BP_enrichment.csv"
            go_bp.to_csv(go_bp_path, index=False)
            enrichment_results['GO_BP'] = go_bp
            print(f"      ✓ Found {len(go_bp)} terms")
    except Exception as e:
        print(f"      ⚠ GO failed: {str(e)[:60]}")
    
    # KEGG
    print("  [2/2] KEGG 2021...")
    try:
        kegg = enrichr(gene_symbols,
                      gene_sets=['KEGG_2021_Human'],
                      outdir=str(OUTPUT_DIR / 'kegg'))
        if kegg is not None:
            # Save results
            kegg_path = OUTPUT_DIR / "KEGG_enrichment.csv"
            kegg.to_csv(kegg_path, index=False)
            enrichment_results['KEGG'] = kegg
            print(f"      ✓ Found {len(kegg)} pathways")
    except Exception as e:
        print(f"      ⚠ KEGG failed: {str(e)[:60]}")

except ImportError:
    print("  ⚠ gseapy not installed")

# === Generate Visualization if we have results ===
print(f"\nGenerating visualizations...")

if enrichment_results:
    # Plot GO results if available
    if 'GO_BP' in enrichment_results:
        go_df = enrichment_results['GO_BP']
        if len(go_df) > 0:
            # Top 15 by p-value
            top_go = go_df.nsmallest(15, 'Adjusted P-value').copy()
            
            fig, ax = plt.subplots(figsize=(12, 8))
            
            # Sort by p-value for plotting
            top_go = top_go.sort_values('Adjusted P-value', ascending=True)
            
            # Create bar plot with -log10(p-value)
            y_pos = np.arange(len(top_go))
            p_values = -np.log10(top_go['Adjusted P-value'].values + 1e-300)
            
            bars = ax.barh(y_pos, p_values, alpha=0.8, color='steelblue', edgecolor='black')
            
            # Color code by significance
            for i, bar in enumerate(bars):
                if p_values[i] > -np.log10(0.001):
                    bar.set_color('darkgreen')
                elif p_values[i] > -np.log10(0.01):
                    bar.set_color('steelblue')
                else:
                    bar.set_color('lightblue')
            
            ax.set_yticks(y_pos)
            ax.set_yticklabels([t[:60] for t in top_go['Term'].values], fontsize=9)
            ax.set_xlabel('-log10(Adjusted P-value)', fontsize=11, fontweight='bold')
            ax.set_title('Top 15 GO Biological Process Terms', fontsize=12, fontweight='bold')
            ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', lw=2, label='p=0.05')
            ax.axvline(x=-np.log10(0.001), color='darkred', linestyle='--', lw=2, label='p=0.001')
            ax.legend()
            ax.grid(axis='x', alpha=0.3)
            
            fig.tight_layout()
            go_plot = OUTPUT_DIR / "GO_BP_top15.jpg"
            plt.savefig(go_plot, dpi=150, format='jpg')
            print(f"  ✓ GO visualization saved to {go_plot}")
            plt.close()
    
    # Plot KEGG results if available
    if 'KEGG' in enrichment_results:
        kegg_df = enrichment_results['KEGG']
        if len(kegg_df) > 0:
            # Top 10 by p-value
            top_kegg = kegg_df.nsmallest(10, 'Adjusted P-value').copy()
            
            fig, ax = plt.subplots(figsize=(12, 6))
            
            # Sort by p-value for plotting
            top_kegg = top_kegg.sort_values('Adjusted P-value', ascending=True)
            
            # Create bar plot with -log10(p-value)
            y_pos = np.arange(len(top_kegg))
            p_values = -np.log10(top_kegg['Adjusted P-value'].values + 1e-300)
            
            bars = ax.barh(y_pos, p_values, alpha=0.8, color='coral', edgecolor='black')
            
            # Color code by significance
            for i, bar in enumerate(bars):
                if p_values[i] > -np.log10(0.001):
                    bar.set_color('darkred')
                elif p_values[i] > -np.log10(0.01):
                    bar.set_color('coral')
                else:
                    bar.set_color('lightsalmon')
            
            ax.set_yticks(y_pos)
            ax.set_yticklabels([t[:70] for t in top_kegg['Term'].values], fontsize=9)
            ax.set_xlabel('-log10(Adjusted P-value)', fontsize=11, fontweight='bold')
            ax.set_title('Top 10 KEGG Pathways', fontsize=12, fontweight='bold')
            ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', lw=2, label='p=0.05')
            ax.axvline(x=-np.log10(0.001), color='darkred', linestyle='--', lw=2, label='p=0.001')
            ax.legend()
            ax.grid(axis='x', alpha=0.3)
            
            fig.tight_layout()
            kegg_plot = OUTPUT_DIR / "KEGG_top10.jpg"
            plt.savefig(kegg_plot, dpi=150, format='jpg')
            print(f"  ✓ KEGG visualization saved to {kegg_plot}")
            plt.close()

# === Generate Summary Report ===
print(f"\nGenerating enrichment summary report...")

report_path = OUTPUT_DIR / "enrichment_summary.txt"
with open(report_path, 'w', encoding='utf-8') as f:
    f.write("=" * 80 + "\n")
    f.write("GO AND KEGG ENRICHMENT ANALYSIS - GBM BIOMARKERS\n")
    f.write("=" * 80 + "\n\n")
    
    f.write("ANALYSIS SUMMARY\n")
    f.write("-" * 80 + "\n")
    f.write(f"Input genes: {len(gene_symbols)} mRNA biomarkers\n")
    f.write(f"Analysis date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    
    if 'GO_BP' in enrichment_results and len(enrichment_results['GO_BP']) > 0:
        go_df = enrichment_results['GO_BP']
        f.write(f"GO BIOLOGICAL PROCESS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total significant terms: {len(go_df[go_df['Adjusted P-value'] < 0.05])}\n\n")
        
        f.write("TOP 15 TERMS\n")
        f.write(f"{'Rank':<6} {'Term':<50} {'P-value':<15} {'Genes':<10}\n")
        f.write("-" * 80 + "\n")
        
        top_go = go_df.nsmallest(15, 'Adjusted P-value')
        for idx, (_, row) in enumerate(top_go.iterrows(), 1):
            term = str(row.get('Term', ''))[:48]
            pval = row.get('Adjusted P-value', 1.0)
            genes = len(str(row.get('Genes', '')).split(';')) if 'Genes' in row else 0
            f.write(f"{idx:<6} {term:<50} {pval:<15.2e} {genes:<10}\n")
        
        f.write("\n")
    else:
        f.write("GO BIOLOGICAL PROCESS\n")
        f.write("-" * 80 + "\n")
        f.write("No significant terms found or analysis not run.\n\n")
    
    if 'KEGG' in enrichment_results and len(enrichment_results['KEGG']) > 0:
        kegg_df = enrichment_results['KEGG']
        f.write(f"KEGG PATHWAYS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total significant pathways: {len(kegg_df[kegg_df['Adjusted P-value'] < 0.05])}\n\n")
        
        f.write("TOP 10 PATHWAYS\n")
        f.write(f"{'Rank':<6} {'Pathway':<50} {'P-value':<15} {'Genes':<10}\n")
        f.write("-" * 80 + "\n")
        
        top_kegg = kegg_df.nsmallest(10, 'Adjusted P-value')
        for idx, (_, row) in enumerate(top_kegg.iterrows(), 1):
            pathway = str(row.get('Term', ''))[:48]
            pval = row.get('Adjusted P-value', 1.0)
            genes = len(str(row.get('Genes', '')).split(';')) if 'Genes' in row else 0
            f.write(f"{idx:<6} {pathway:<50} {pval:<15.2e} {genes:<10}\n")
        
        f.write("\n")
    else:
        f.write("KEGG PATHWAYS\n")
        f.write("-" * 80 + "\n")
        f.write("No significant pathways found or analysis not run.\n\n")
    
    f.write("INTERPRETATION\n")
    f.write("-" * 80 + "\n")
    f.write("GO terms and KEGG pathways indicate the functional roles and biological\n")
    f.write("processes associated with the identified GBM biomarkers.\n\n")
    
    f.write("Significant enrichment (p < 0.05) suggests that the biomarkers are\n")
    f.write("collectively involved in specific biological processes or pathways.\n\n")
    
    f.write("FILES GENERATED\n")
    f.write("-" * 80 + "\n")
    f.write(f"Location: {OUTPUT_DIR}/\n\n")
    
    if (OUTPUT_DIR / "GO_BP_enrichment.csv").exists():
        f.write("- GO_BP_enrichment.csv (full GO analysis results)\n")
    if (OUTPUT_DIR / "GO_BP_top15.jpg").exists():
        f.write("- GO_BP_top15.jpg (visualization)\n")
    if (OUTPUT_DIR / "KEGG_enrichment.csv").exists():
        f.write("- KEGG_enrichment.csv (full KEGG analysis results)\n")
    if (OUTPUT_DIR / "KEGG_top10.jpg").exists():
        f.write("- KEGG_top10.jpg (visualization)\n")
    
    f.write("\nNEXT STEPS\n")
    f.write("-" * 80 + "\n")
    f.write("1. Review GO terms to understand biological functions\n")
    f.write("2. Check KEGG pathways for known GBM mechanisms\n")
    f.write("3. Validate findings against published literature\n")
    f.write("4. Identify druggable targets in enriched pathways\n")
    f.write("5. Consider pathway-based therapeutic strategies\n")

print(f"✓ Report saved to {report_path}")

# === Create Quick Reference ===
quick_ref = BIOMARKER_DIR / "ENRICHMENT_QUICK_REFERENCE.txt"
with open(quick_ref, 'w', encoding='utf-8') as f:
    f.write("GBM BIOMARKER ENRICHMENT - QUICK REFERENCE\n")
    f.write("=" * 70 + "\n\n")
    
    if 'GO_BP' in enrichment_results and len(enrichment_results['GO_BP']) > 0:
        f.write("TOP GO BIOLOGICAL PROCESSES\n")
        f.write("-" * 70 + "\n")
        top_go = enrichment_results['GO_BP'].nsmallest(5, 'Adjusted P-value')
        for idx, (_, row) in enumerate(top_go.iterrows(), 1):
            f.write(f"{idx}. {row['Term']}\n")
            f.write(f"   p-value: {row['Adjusted P-value']:.2e}\n\n")
    
    if 'KEGG' in enrichment_results and len(enrichment_results['KEGG']) > 0:
        f.write("TOP KEGG PATHWAYS\n")
        f.write("-" * 70 + "\n")
        top_kegg = enrichment_results['KEGG'].nsmallest(5, 'Adjusted P-value')
        for idx, (_, row) in enumerate(top_kegg.iterrows(), 1):
            f.write(f"{idx}. {row['Term']}\n")
            f.write(f"   p-value: {row['Adjusted P-value']:.2e}\n\n")

print(f"✓ Quick reference saved to {quick_ref}")

print("\n" + "=" * 80)
print("ENRICHMENT ANALYSIS COMPLETE")
print("=" * 80)
print(f"\nOutput directory: {OUTPUT_DIR}/")
print(f"\nGenerated files:")
if (OUTPUT_DIR / "GO_BP_enrichment.csv").exists():
    print(f"  ✓ GO_BP_enrichment.csv")
if (OUTPUT_DIR / "GO_BP_top15.jpg").exists():
    print(f"  ✓ GO_BP_top15.jpg")
if (OUTPUT_DIR / "KEGG_enrichment.csv").exists():
    print(f"  ✓ KEGG_enrichment.csv")
if (OUTPUT_DIR / "KEGG_top10.jpg").exists():
    print(f"  ✓ KEGG_top10.jpg")
print(f"  ✓ enrichment_summary.txt")
print(f"  ✓ ENRICHMENT_QUICK_REFERENCE.txt")
