#!/usr/bin/env python3
"""Biomarker discovery from DE gene lists and model predictions."""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

RESULTS_DIR = Path("results")
CV_DIR = RESULTS_DIR / "cv_results_prod"
OUTPUT_DIR = RESULTS_DIR / "biomarkers"


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Load CV results
    fold_results = pd.read_csv(CV_DIR / "fold_results.csv")
    scores = np.load(CV_DIR / "all_fold_scores.npy")
    labels = np.load(CV_DIR / "all_fold_labels.npy")
    
    print(f"Loaded {len(fold_results)} fold results, {len(scores)} predictions")
    
    # Load node mapping
    node_mapping = {}
    if (RESULTS_DIR / "node_mapping.json").exists():
        with open(RESULTS_DIR / "node_mapping.json") as f:
            node_mapping = json.load(f)
    
    # Load DE lists
    de_info = {}
    
    if (RESULTS_DIR / "DE_mRNA_logFC_abs_gt1.tsv").exists():
        df = pd.read_csv(RESULTS_DIR / "DE_mRNA_logFC_abs_gt1.tsv", sep='\t')
        for _, row in df.iterrows():
            de_info[row['ensembl_id']] = {
                'type': 'mRNA', 
                'symbol': row.get('gene_symbol', ''), 
                'lfc': row.get('log2FoldChange', 0)
            }
    
    if (RESULTS_DIR / "DE_lncRNA_logFC_abs_gt1.tsv").exists():
        df = pd.read_csv(RESULTS_DIR / "DE_lncRNA_logFC_abs_gt1.tsv", sep='\t')
        for _, row in df.iterrows():
            de_info[row['ensembl_id']] = {
                'type': 'lncRNA', 
                'symbol': row.get('gene_symbol', ''), 
                'lfc': row.get('log2FoldChange', 0)
            }
    
    if (RESULTS_DIR / "DE_miRNA_logFC_abs_gt1.tsv").exists():
        df = pd.read_csv(RESULTS_DIR / "DE_miRNA_logFC_abs_gt1.tsv", sep='\t')
        for _, row in df.iterrows():
            de_info[row['mirna_id']] = {
                'type': 'miRNA', 
                'symbol': row['mirna_id'], 
                'lfc': row.get('log2FoldChange', 0)
            }
    
    print(f"Loaded {len(de_info)} DE genes")
    
    # Rank biomarkers by DE strength
    type_weight = {'mRNA': 1.0, 'lncRNA': 0.8, 'miRNA': 0.9}
    
    candidates = []
    for gene_id, info in de_info.items():
        confidence = abs(info['lfc']) * type_weight.get(info['type'], 0.5)
        candidates.append({
            'gene_id': gene_id,
            'symbol': info['symbol'],
            'type': info['type'],
            'log2FC': info['lfc'],
            'confidence_score': confidence
        })
    
    biomarker_df = pd.DataFrame(candidates).sort_values('confidence_score', ascending=False)
    top_biomarkers = biomarker_df.head(100)
    
    # Save results
    top_biomarkers.to_csv(OUTPUT_DIR / "top_biomarkers.csv", index=False)
    print(f"Saved {len(top_biomarkers)} biomarkers to {OUTPUT_DIR / 'top_biomarkers.csv'}")
    
    # Optional: GO/KEGG enrichment
    try:
        from gseapy import enrichr
        
        mrna_genes = top_biomarkers[top_biomarkers['type'] == 'mRNA']['symbol'].tolist()
        mrna_genes = [g for g in mrna_genes if g.strip()]
        
        if mrna_genes:
            print(f"Running enrichment on {len(mrna_genes)} mRNA genes...")
            try:
                go = enrichr(mrna_genes, gene_sets='GO_Biological_Process_2021', 
                           outdir=str(OUTPUT_DIR / 'go_bp'))
                if go is not None:
                    go.to_csv(OUTPUT_DIR / "go_bp_enrichment.csv")
            except:
                pass
            
            try:
                kegg = enrichr(mrna_genes, gene_sets='KEGG_2019', 
                             outdir=str(OUTPUT_DIR / 'kegg'))
                if kegg is not None:
                    kegg.to_csv(OUTPUT_DIR / "kegg_enrichment.csv")
            except:
                pass
                
    except ImportError:
        print("gseapy not installed, skipping enrichment")
    
    # Print top 10
    print("\nTop 10 biomarker candidates:")
    for i, (_, row) in enumerate(top_biomarkers.head(10).iterrows(), 1):
        print(f"  {i}. {row['symbol']} ({row['type']}) log2FC={row['log2FC']:.2f}")


if __name__ == "__main__":
    main()
