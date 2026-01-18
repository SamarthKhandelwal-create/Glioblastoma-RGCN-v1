# -*- coding: utf-8 -*-
"""
build_node_features.py

Generates a 10-dimensional feature matrix for the Genes in the GBM ceRNA network.
Feature 10 (Tumor vs Normal) is loaded from previous results.
Features 1-9 are calculated by running PyDESeq2 on TCGA-GBM subsets based on clinical metadata.

Output: results/node_features_matrix.csv
"""

import os
import sys
import numpy as np
import pandas as pd
import xenaPython as xp
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pathlib import Path

# --- Constants ---
RESULTS_DIR = Path("results")
DE_GENES_FILE = RESULTS_DIR / "DE_genes_significant.tsv"
OUTPUT_FILE = RESULTS_DIR / "node_features_matrix.csv"

TOIL_HOST = "https://toil.xenahubs.net"
TOIL_COUNTS = "TcgaTargetGtex_gene_expected_count"
LEGACY_HOST = "https://tcga.xenahubs.net"
CLINICAL_DATASET = "TCGA.GBM.sampleMap/GBM_clinicalMatrix"

# Feature definitions
# (Name, XenaField, Comparison_Logic)
# Logic: (Group1_HighRisk, Group0_LowRisk)
FEATURES = [
    {
        "id": 1, "name": "Age", "field": "age_at_initial_pathologic_diagnosis",
        "type": "numeric_binary", "cutoff": 60, "label_1": "Old (>60)", "label_0": "Young (<=60)"
    },
    {
        "id": 2, "name": "Sex", "field": "gender",
        "type": "categorical", "val_1": "MALE", "val_0": "FEMALE" # check caps
    },
    {
        "id": 3, "name": "Vital", "field": "vital_status",
        "type": "categorical", "val_1": "DECEASED", "val_0": "LIVING"
    },
    {
        "id": 4, "name": "Survival", "field": "days_to_death",
        "type": "numeric_median", "label_1": "Short (<Median)", "label_0": "Long (>=Median)"
    },
    {
        "id": 5, "name": "IDH", "field": "G_CIMP_STATUS",
        "type": "categorical", "val_1": "G-CIMP", "val_0": "non-G-CIMP" # G-CIMP is Mutant (better?), actually G-CIMP is usually better prognosis, but we just need a diff. Let's stick to G-CIMP vs Non.
        # Wait, usually WT is aggressive. G-CIMP is hypermethylated (mutant-like).
        # Let's define 1 as Aggressive (Non-G-CIMP) and 0 as G-CIMP if we want consistent risk direction?
        # Actually user said: IDH Wildtype (Aggressive) vs IDH Mutant.
        # G-CIMP ~= IDH Mutant. so Non-G-CIMP ~= WT.
        # So val_1 = "non-G-CIMP", val_0 = "G-CIMP"
    },
    {
        "id": 6, "name": "Subtype", "field": "GeneExp_Subtype",
        "type": "categorical", "val_1": "Mesenchymal", "val_0": "Proneural" # Or just vs Rest? Let's do Mesenchymal vs Not.
        # Special handler needed for "vs Others"
    },
    {
        "id": 7, "name": "Radiation", "field": "radiation_therapy",
        "type": "categorical", "val_1": "YES", "val_0": "NO"
    },
    {
        "id": 8, "name": "KPS", "field": "karnofsky_performance_score",
        "type": "numeric_binary", "cutoff": 80, "label_1": "Low (<80)", "label_0": "High (>=80)",
        "invert": True # cutoff logic is > is good. we want 1=Bad usually? or doesn't matter for GCN.
        # Let's just follow split.
    },
    {
        "id": 9, "name": "Chemo", "field": "chemo_therapy",
        "type": "categorical", "val_1": "YES", "val_0": "NO"
    }
]

def load_master_gene_list():
    m_path = RESULTS_DIR / "DE_mRNA_logFC_abs_gt1.tsv"
    l_path = RESULTS_DIR / "DE_lncRNA_logFC_abs_gt1.tsv"
    mir_path = RESULTS_DIR / "DE_miRNA_logFC_abs_gt1.tsv"
    
    if not m_path.exists() or not l_path.exists() or not mir_path.exists():
        print(f"Error: gt1 filtered files (mRNA, lncRNA, or miRNA) not found.")
        sys.exit(1)
        
    m_df = pd.read_csv(m_path, sep="\t")
    l_df = pd.read_csv(l_path, sep="\t")
    mir_df = pd.read_csv(mir_path, sep="\t")
    
    # Standardize ID columns
    # mRNA/lncRNA have 'ensembl_id', miRNA has 'mirna_id'
    # We will rename 'mirna_id' to 'ensembl_id' for consistency in the matrix index
    mir_df = mir_df.rename(columns={'mirna_id': 'ensembl_id'})
    
    combined = pd.concat([m_df, l_df, mir_df], axis=0)
    
    # Return list of IDs and the mixed DF
    return combined["ensembl_id"].tolist(), combined

def fetch_counts(samples, genes):
    import xenaPython as xq
    # Fetch in chunks
    print(f"Fetching counts for {len(genes)} genes across {len(samples)} samples...")
    chunk_size = 500
    all_df = []
    
    # Xena gene ID matching can be tricky. TOIL uses Ensembl without version usually?
    # Our ensembl_ids might have version.
    # Let's try to strip versions for query if needed, but TOIL usually has versions format "ENSG.v". 
    # Let's check a few.
    # If the previous script worked, genes are correct.
    
    # Filter for Ensembl IDs only for Xena (skip miRNAs which are 'hsa-mir-...')
    # miRNAs won't be found in the gene expression dataset anyway.
    ens_genes = [g for g in genes if str(g).startswith("ENSG")]
    
    print(f"Fetching counts for {len(ens_genes)} ENSG genes (skipping {len(genes)-len(ens_genes)} non-ENSG) across {len(samples)} samples...")
    
    # Only fetch for valid ENSG
    genes_to_fetch = ens_genes
    
    for i in range(0, len(genes_to_fetch), chunk_size):
        chunk = genes_to_fetch[i:i+chunk_size]
        try:
             # probeMap for Toil is usually just the ID
             # Xena returns [metadata, [values]] structure sometimes
             raw = xq.dataset_probe_values(TOIL_HOST, TOIL_COUNTS, samples, chunk)
             
             values = []
             if raw and isinstance(raw, list):
                 if len(raw) == len(chunk) and isinstance(raw[0], list): # Standard list of lists?
                     values = raw
                 elif len(raw) > 1 and isinstance(raw[1], list) and len(raw[1]) == len(chunk): # Wrapped [meta, data]
                     values = raw[1]
                 elif len(raw) > 0 and isinstance(raw[0], list) and len(raw[0]) == len(chunk): # Maybe wrapped [data]
                     values = raw[0]
                     
             if not values:
                 print(f"    [Warning] Empty values for chunk {i}")
                 continue
                 
             df = pd.DataFrame(values, index=chunk, columns=samples)
             all_df.append(df)
        except Exception as e:
             print(f"    [Error] fetch_counts chunk {i}: {e}")
             pass
             
    if not all_df:
        return pd.DataFrame() # Returns empty if no ENSG or failure
        
    result = pd.concat(all_df, axis=0)
    
    # Reindex to full gene list (filling missing/miRNAs with NaN -> 0 later)
    # This ensures the dataframe matches the 'genes' argument order if needed, 
    # but strictly we just need the IDs to exist.
    # Actually, main() loop expects counts_df to contain the genes it wants to process.
    # We will return the subset. The DESeq loop will handle missing genes by skipping them?
    # No, run_deseq_features uses counts_df[common]. 
    # So if miRNAs are missing from counts_df, their logFC will be NaN?
    # run_deseq_feature returns "None" if samples < 6.
    # Wait, we need to return a DF that can be queried.
    return result

def run_deseq_feature(counts_df, metadata, feature_col, feat_def):
    """
    Run DESeq2 comparing metadata[feature_col] == 1 vs 0.
    Returns Series of log2FoldChange for the genes.
    """
    # 1. Align
    common = counts_df.columns.intersection(metadata.index)
    counts_sub = counts_df[common].T # samples x genes
    meta_sub = metadata.loc[common]
    
    # 2. Filter for valid 0/1 (drop NaNs or -1)
    valid = meta_sub[meta_sub[feature_col].isin([0, 1])].index
    if len(valid) < 6: # Need min samples
        print(f"  [Skip] Not enough samples for {feat_def['name']} (n={len(valid)})")
        return None
        
    counts_run = counts_sub.loc[valid]
    meta_run = meta_sub.loc[valid]
    
    # Round counts for DESeq2 input
    counts_run = counts_run.fillna(0).clip(lower=0).round().astype(int)
    
    # 3. Design
    # FIX: Ensure design factor is string, not float
    try:
        # Convert 1.0 -> 1 -> "1"
        meta_run[feature_col] = meta_run[feature_col].astype(int).astype(str)
    except:
        meta_run[feature_col] = meta_run[feature_col].astype(str)
        
    meta_run[feature_col] = meta_run[feature_col].astype("category")
    
    print(f"  [DESeq2] {feat_def['name']}: {len(valid)} samples. 1-vs-0.")
    try:
        dds = DeseqDataSet(
            counts=counts_run,
            metadata=meta_run,
            design_factors=feature_col,
            quiet=True,
            n_cpus=2 # Assuming we have some cores
        )
        dds.deseq2()
        stat = DeseqStats(dds, contrast=(feature_col, "1", "0"), quiet=True, n_cpus=2)
        stat.summary()
        
        res = stat.results_df
        return res["log2FoldChange"]
    except Exception as e:
        print(f"  [Error] DESeq2 failed for {feat_def['name']}: {e}")
        return None

def fetch_field_data(samples, field):
    """
    Fetch data for a single field, handling weird Xena API structure 
    ([None, [[values]]]) and decoding codes.
    """
    import xenaPython as xp
    try:
        # Fetch values
        # The API seems to return [None, [[values...]]] or similar
        raw = xp.dataset_probe_values(LEGACY_HOST, CLINICAL_DATASET, samples, [field])
        
        # Find the actual data list
        data_list = None
        if raw and isinstance(raw, list):
            # Look for inner list of length samples
            for item in raw:
                if item and isinstance(item, list) and len(item) == 1 and isinstance(item[0], list) and len(item[0]) == len(samples):
                     data_list = item[0]
                     break
                # Backup: sometimes it might be just [ [values] ]
                if isinstance(item, list) and len(item) == len(samples) and not isinstance(item[0], list):
                     data_list = item
                     break
                     
        if data_list is None:
            # Fallback: maybe raw[1][0]
            if len(raw) > 1 and isinstance(raw[1], list) and len(raw[1]) > 0:
                 data_list = raw[1][0]
            elif len(raw) > 0 and isinstance(raw[0], list):
                 data_list = raw[0]
            else:
                 print(f"    [Error] Could not find data in response for {field}")
                 return pd.Series([np.nan]*len(samples), index=samples)

        # Handle 'NaN' string -> np.nan
        data_series = pd.Series(data_list, index=samples)
        data_series = data_series.replace('NaN', np.nan)
        
        # Fetch codes if available
        codes_res = xp.field_codes(LEGACY_HOST, CLINICAL_DATASET, [field])
        if codes_res and codes_res[0].get("code"):
            code_str = codes_res[0]["code"]
            # codes are tab separated. index = value
            code_map = code_str.split('\t')
            print(f"    [Decoding] {field} with {len(code_map)} codes (e.g. 0={code_map[0]})")
            
            def decode(x):
                try:
                    i = int(float(x))
                    if 0 <= i < len(code_map):
                        return code_map[i]
                except:
                    pass
                return x # keep as is if not int
                
            data_series = data_series.apply(decode)
            
        return data_series
        
    except Exception as e:
        print(f"    [Exception] fetching {field}: {e}")
        return pd.Series([np.nan]*len(samples), index=samples)

def main():
    if not RESULTS_DIR.exists():
        RESULTS_DIR.mkdir()

    # 1. Genes
    ensembl_ids, master_df = load_master_gene_list()
    print(f"Loaded {len(ensembl_ids)} genes.")
    
    # 2. Samples
    print("Fetching clinical metadata samples...")
    clinical_samples = xp.dataset_samples(LEGACY_HOST, CLINICAL_DATASET, None)
    
    # 3. Fetch each field and build clean DF
    clin_df = pd.DataFrame(index=clinical_samples)
    print(f"Querying fields one by one for {len(clinical_samples)} samples...")
    
    for feat in FEATURES:
        f = feat["field"]
        print(f"  Fetching {f}...")
        clin_df[f] = fetch_field_data(clinical_samples, f)
        print(f"    Unique values: {clin_df[f].dropna().unique()[:5]}") # debug
    
    # 3. Clean Metadata -> Create design matrix
    design_meta = pd.DataFrame(index=clin_df.index)
    
    for feat in FEATURES:
        col = feat["field"]
        name = feat["name"]
        
        # Helper to numeric
        s = pd.to_numeric(clin_df[col], errors='coerce')
        
        if feat["type"] == "numeric_binary":
            # 1 if > cutoff (or inverted)
            if feat.get("invert"):
                mask1 = s < feat["cutoff"]
                mask0 = s >= feat["cutoff"]
            else:
                mask1 = s > feat["cutoff"]
                mask0 = s <= feat["cutoff"]
            design_meta.loc[mask1, name] = 1
            design_meta.loc[mask0, name] = 0
            
        elif feat["type"] == "numeric_median":
            median = s.median()
            print(f"  {name} median: {median}")
            mask1 = s < median # Short survival is "High Risk" usually? 
                               # Paper: Short Survival (< Median) vs Long.
                               # Let's say 1 = Short (<), 0 = Long (>=)
            mask0 = s >= median
            design_meta.loc[mask1, name] = 1
            design_meta.loc[mask0, name] = 0
            
        elif feat["type"] == "categorical":
            clean_s = clin_df[col].astype(str).str.upper()
            
            # Custom logic for IDH (Non-G-CIMP vs G-CIMP)
            if name == "IDH":
                 # G-CIMP is 0 (Good), Non-G-CIMP is 1 (Bad)
                 # Check exact values with relaxed matching (space vs hyphen)
                 # "G-CIMP", "NON-G-CIMP" vs "NON G-CIMP"
                 design_meta.loc[clean_s.str.contains("NON"), name] = 1
                 design_meta.loc[clean_s.str.contains("G-CIMP") & (~clean_s.str.contains("NON")), name] = 0
                 
            # Custom for Subtype
            elif name == "Subtype":
                 design_meta.loc[clean_s == feat["val_1"].upper(), name] = 1
                 design_meta.loc[(clean_s != "NAN") & (clean_s != feat["val_1"].upper()), name] = 0
                 
            else:
                val1 = feat["val_1"].upper()
                val0 = feat["val_0"].upper()
                design_meta.loc[clean_s == val1, name] = 1
                design_meta.loc[clean_s == val0, name] = 0
                
    # 4. Fetch Counts for these clinical samples
    # Intersection of clinical samples and Toil samples
    toil_samples = xp.dataset_samples(TOIL_HOST, TOIL_COUNTS, None)
    valid_samples = list(set(design_meta.index) & set(toil_samples))
    print(f"Aligned samples: {len(valid_samples)} (Clinical + Toil)")
    
    if len(valid_samples) < 10:
        print("Error: Too few overlapping samples.")
        sys.exit(1)
        
    counts_df = fetch_counts(valid_samples, ensembl_ids)
    print(f"Counts DF shape: {counts_df.shape}")
    print(f"Counts DF total NaNs: {counts_df.isna().sum().sum()}")
    print(f"Counts DF sample values: {counts_df.iloc[:5, :5].values}")

    # 5. Run Loop
    results_matrix = pd.DataFrame(index=ensembl_ids)
    # ... Health feature logic ...
    if "log2FoldChange" in master_df.columns:
        results_matrix["Health"] = master_df.set_index("ensembl_id")["log2FoldChange"]
    else:
        results_matrix["Health"] = 0.0
    
    # DEBUG: Limit features
    DEBUG_FEATURES = ["Age", "Sex"] 
    
    for feat in FEATURES:
        name = feat["name"]
        # if name not in DEBUG_FEATURES: continue # Skip for debug
        
        print(f"\n--- Processing Feature: {name} ---")
        print(f"  Design META distribution:\n{design_meta[name].value_counts(dropna=False)}")
        
        logfc = run_deseq_feature(counts_df, design_meta, name, feat)
        
        if logfc is not None:
             # This assignment matches index. 
             # Use .reindex to force alignment with results_matrix (indexes are genes)
             # logfc only contains ENSG genes found in counts_df.
             aligned_logfc = logfc.reindex(results_matrix.index, fill_value=0.0)
             results_matrix[name] = aligned_logfc
             # print(f"  [Success] Got {len(logfc)} LogFC values.")
        else:
             print(f"  [Failure] LogFC is None.")
             results_matrix[name] = 0.0 # Fill neutral if failed
             
    # 6. Final Polish
    # Fill NaNs (genes that didn't converge or low count in subset)
    results_matrix = results_matrix.fillna(0.0)
    
    # Z-score normalization
    # (x - mean) / std
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    scaled_vals = scaler.fit_transform(results_matrix)
    final_df = pd.DataFrame(scaled_vals, index=results_matrix.index, columns=results_matrix.columns)
    
    # Save
    final_df.to_csv(OUTPUT_FILE)
    print(f"\nSuccess! Saved feature matrix to {OUTPUT_FILE}")
    print(final_df.head())

if __name__ == "__main__":
    main()
