import pandas as pd
from pathlib import Path
import os


# Results dir
RESULTS_DIR = Path(__import__("os").environ.get("RESULTS_DIR", "results"))
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

M2G_NAME = "mir2gene.csv"
M2L_NAME = "mir2lnc.csv"

# Determine data/raw location relative to repo root
if Path("data/raw").exists() and (Path("data/raw") / M2G_NAME).exists():
    DATA_DIR = Path("data/raw")
elif (Path(__file__).parent / "../../data/raw").resolve().exists():
    DATA_DIR = (Path(__file__).parent / "../../data/raw").resolve()
else:
    DATA_DIR = Path("data/raw")

MIR2GENE_PATH = DATA_DIR / M2G_NAME
MIR2LNC_PATH = DATA_DIR / M2L_NAME
INTERACTIONS_PATH = RESULTS_DIR / "interactions.csv"

DE_MRNA_PATH = RESULTS_DIR / "DE_mRNA_logFC_abs_gt1.tsv"
DE_LNCRNA_PATH = RESULTS_DIR / "DE_lncRNA_logFC_abs_gt1.tsv"
DE_MIRNA_PATH = RESULTS_DIR / "DE_miRNA_logFC_abs_gt1.tsv"


def normalize_mir_id(s):
    if pd.isna(s):
        return ""
    return str(s).strip().lower()


MIN_REF_COUNT = 1


def count_refs(lit_str):
    if pd.isna(lit_str) or str(lit_str) == "nan":
        return 0
    return len(str(lit_str).split("|"))


def main():
    print("--- Integrating ENCORI (Manual Downloads) ---")
    print(f"Filtering for interactions with >= {MIN_REF_COUNT} literature references.")

    # Basic checks
    for p in (DE_MRNA_PATH, DE_LNCRNA_PATH, DE_MIRNA_PATH):
        if not p.exists():
            print(f"Missing required DE file: {p}")
            return

    # 1. Load Mappings
    print("Loading DE lists...")
    de_mrna = pd.read_csv(DE_MRNA_PATH, sep='\t')
    de_lncrna = pd.read_csv(DE_LNCRNA_PATH, sep='\t')
    de_mirna = pd.read_csv(DE_MIRNA_PATH, sep='\t')

    mrna_map = dict(zip(de_mrna['gene_symbol'], de_mrna['ensembl_id']))
    lncrna_map = dict(zip(de_lncrna['gene_symbol'], de_lncrna['ensembl_id']))

    my_mirna_list = de_mirna['mirna_id'].astype(str).tolist()
    my_mirnas = set([normalize_mir_id(x) for x in my_mirna_list])


    def resolve_mirna(encori_id):
        encori_norm = normalize_mir_id(encori_id)
        if not encori_norm:
            return None
        if encori_norm in my_mirnas:
            # return DE-format ID (preserve original case from list)
            for m in my_mirna_list:
                if normalize_mir_id(m) == encori_norm:
                    return m
        # substring heuristics (match stems)
        for m in my_mirna_list:
            nm = normalize_mir_id(m)
            if nm == encori_norm:
                return m
            if nm in encori_norm and encori_norm.startswith(nm):
                return m
            if encori_norm in nm:
                return m
        return None


    new_edges = []

    # Process mir2gene
    if MIR2GENE_PATH.exists():
        print(f"Processing {MIR2GENE_PATH}...")
        df_genes = pd.read_csv(MIR2GENE_PATH, on_bad_lines='skip')
        count = 0
        for _, row in df_genes.iterrows():
            mir_raw = row.get('ID') or row.get('miRNA') or ''
            gene_sym = row.get('Target') or row.get('Gene') or ''
            my_mir = resolve_mirna(mir_raw)
            if not my_mir:
                continue
            my_gene = mrna_map.get(gene_sym)
            if not my_gene:
                continue
            refs = count_refs(row.get('Literature'))
            if refs < MIN_REF_COUNT:
                continue
            new_edges.append([my_mir, my_gene])
            count += 1
        print(f"  Added {count} miRNA-mRNA interactions.")
    else:
        print(f"  {MIR2GENE_PATH} not found.")

    # Process mir2lnc
    if MIR2LNC_PATH.exists():
        print(f"Processing {MIR2LNC_PATH}...")
        df_lnc = pd.read_csv(MIR2LNC_PATH, on_bad_lines='skip')
        count = 0
        for _, row in df_lnc.iterrows():
            mir_raw = row.get('ID') or row.get('miRNA') or ''
            target_sym = row.get('Target') or row.get('Gene') or ''
            my_mir = resolve_mirna(mir_raw)
            if not my_mir:
                continue
            my_lnc = lncrna_map.get(target_sym)
            if not my_lnc:
                continue
            new_edges.append([my_lnc, my_mir])
            count += 1
        print(f"  Added {count} lncRNA-miRNA interactions.")
    else:
        print(f"  {MIR2LNC_PATH} not found.")

    # Save merged interactions
    if new_edges:
        print(f"Merging {len(new_edges)} total new edges...")
        new_df = pd.DataFrame(new_edges, columns=['gene_A', 'gene_B'])
        if INTERACTIONS_PATH.exists():
            old_df = pd.read_csv(INTERACTIONS_PATH)
            combined = pd.concat([old_df, new_df], ignore_index=True)
        else:
            combined = new_df
        final_df = combined.drop_duplicates()
        final_df.to_csv(INTERACTIONS_PATH, index=False)
        print(f"Saved to {INTERACTIONS_PATH}. Total interactions: {len(final_df)}")
    else:
        print("No new interactions found matching our DE lists.")




if __name__ == "__main__":
    main()
