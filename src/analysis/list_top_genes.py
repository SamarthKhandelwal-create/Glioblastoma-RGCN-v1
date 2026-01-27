import pandas as pd
from pathlib import Path

RESULTS_DIR = Path(__import__("os").environ.get("RESULTS_DIR", "results"))
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

PREDS_PATH = RESULTS_DIR / "novel_predictions.csv"
DE_MRNA = RESULTS_DIR / "DE_mRNA_logFC_abs_gt1.tsv"
DE_LNCRNA = RESULTS_DIR / "DE_lncRNA_logFC_abs_gt1.tsv"


def strip_version(ens_id):
    """Strip Ensembl version numbers."""
    s = str(ens_id).strip()
    if s.startswith('ENSG'):
        return s.split('.')[0]
    return s


def main():
    if not PREDS_PATH.exists():
        print(f"Error: {PREDS_PATH} not found. Run predictions first.")
        return

    df = pd.read_csv(PREDS_PATH)

    # Build Mapping (Ensembl -> Symbol) with version stripping
    symbol_map = {}

    if DE_MRNA.exists():
        m_df = pd.read_csv(DE_MRNA, sep='\t')
        for _, row in m_df.iterrows():
            ens_id = strip_version(row['ensembl_id'])
            symbol_map[ens_id] = row['gene_symbol']
            symbol_map[str(row['ensembl_id'])] = row['gene_symbol']

    if DE_LNCRNA.exists():
        l_df = pd.read_csv(DE_LNCRNA, sep='\t')
        for _, row in l_df.iterrows():
            ens_id = strip_version(row['ensembl_id'])
            symbol_map[ens_id] = row['gene_symbol']
            symbol_map[str(row['ensembl_id'])] = row['gene_symbol']

    def get_symbol(target_id):
        # Try exact match and version-stripped match
        if target_id in symbol_map:
            return symbol_map[target_id]
        stripped = strip_version(target_id)
        if stripped in symbol_map:
            return symbol_map[stripped]
        return target_id

    df['Gene Symbol'] = df['Target'].apply(get_symbol)

    # Filter for unique genes (best score for each gene)
    unique_genes = df.sort_values('Score', ascending=False).drop_duplicates(subset='Target')

    # Display and save
    print("--- Top 10 High-Confidence Gene Targets ---")
    top_10 = unique_genes[['Gene Symbol', 'miRNA', 'Score']].head(10)
    print(top_10.to_string(index=False))

    # Save top genes
    top_genes_path = RESULTS_DIR / "top_genes.csv"
    unique_genes.to_csv(top_genes_path, index=False)
    print(f"\nAll genes saved to {top_genes_path}")


if __name__ == "__main__":
    main()
