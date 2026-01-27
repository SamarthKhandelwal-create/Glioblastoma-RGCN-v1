import pandas as pd
from pathlib import Path

RESULTS_DIR = Path(__import__("os").environ.get("RESULTS_DIR", "results"))
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

MIRNET_PATH = RESULTS_DIR / "mirnet_interactions.csv"
INTERACTIONS_PATH = RESULTS_DIR / "interactions.csv"

DE_MRNA_PATH = RESULTS_DIR / "DE_mRNA_logFC_abs_gt1.tsv"
DE_LNCRNA_PATH = RESULTS_DIR / "DE_lncRNA_logFC_abs_gt1.tsv"
DE_MIRNA_PATH = RESULTS_DIR / "DE_miRNA_logFC_abs_gt1.tsv"


def strip_ensembl(x):
    if pd.isna(x):
        return x
    s = str(x).strip()
    if s.startswith('ENSG'):
        return s.split('.')[0]
    return s


def main():
    print("--- Integrating miRNet ---")
    if not MIRNET_PATH.exists():
        print(f"Error: {MIRNET_PATH} not found.")
        return

    df = pd.read_csv(MIRNET_PATH)
    print(f"Loaded miRNet data: {len(df)} rows")

    # Load mappings with checks
    if not (DE_MRNA_PATH.exists() and DE_MIRNA_PATH.exists() and DE_LNCRNA_PATH.exists()):
        print("Missing DE files in RESULTS_DIR. Aborting.")
        return

    de_mrna = pd.read_csv(DE_MRNA_PATH, sep='\t')
    de_lncrna = pd.read_csv(DE_LNCRNA_PATH, sep='\t')
    de_mirna = pd.read_csv(DE_MIRNA_PATH, sep='\t')

    mrna_map = {str(k).strip(): strip_ensembl(v) for k, v in zip(de_mrna['gene_symbol'], de_mrna['ensembl_id'])}
    lncrna_map = {str(k).strip(): strip_ensembl(v) for k, v in zip(de_lncrna['gene_symbol'], de_lncrna['ensembl_id'])}

    my_mirna_list = de_mirna['mirna_id'].astype(str).tolist()
    my_mirnas_lower = {m.lower(): m for m in my_mirna_list}

    new_edges = []

    col_source = 'Source'
    col_target = 'Target'
    if col_source not in df.columns:
        col_source = df.columns[0]
        col_target = df.columns[1]
        print(f"Guessing columns: Source={col_source}, Target={col_target}")

    for _, row in df.iterrows():
        src = str(row[col_source]).strip()
        dst = str(row[col_target]).strip()

        mir_id = None
        if src in my_mirna_list:
            mir_id = src
        else:
            mir_id = my_mirnas_lower.get(src.lower())

        if not mir_id:
            # substring heuristic
            for mlow, morig in my_mirnas_lower.items():
                if mlow in src.lower() or src.lower() in mlow:
                    mir_id = morig
                    break

        if mir_id:
            gene_id = mrna_map.get(dst) or lncrna_map.get(dst) or mrna_map.get(dst.upper())
            if gene_id:
                new_edges.append([mir_id, gene_id])

    print(f"Found {len(new_edges)} valid edges from miRNet.")
    if new_edges:
        new_df = pd.DataFrame(new_edges, columns=['gene_A', 'gene_B'])
        if INTERACTIONS_PATH.exists():
            old = pd.read_csv(INTERACTIONS_PATH)
            combined = pd.concat([old, new_df], ignore_index=True).drop_duplicates()
        else:
            combined = new_df.drop_duplicates()
        combined.to_csv(INTERACTIONS_PATH, index=False)
        print(f"Saved merged: {len(combined)} total rows.")


if __name__ == "__main__":
    main()
