import pandas as pd
from pathlib import Path

RESULTS_DIR = Path(__import__("os").environ.get("RESULTS_DIR", "results"))
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DE_MIRNA_PATH = RESULTS_DIR / "DE_miRNA_logFC_abs_gt1.tsv"
OUTPUT_TXT_PATH = RESULTS_DIR / "mirna_list_for_mirnet.txt"


def main():
    if not DE_MIRNA_PATH.exists():
        print(f"Error: Could not find {DE_MIRNA_PATH}")
        return

    print(f"Reading miRNAs from {DE_MIRNA_PATH}...")
    df = pd.read_csv(DE_MIRNA_PATH, sep='\t')

    mirnas = df['mirna_id'].dropna().astype(str).unique().tolist()
    mirnas = sorted([m.strip() for m in mirnas])

    print(f"Found {len(mirnas)} unique miRNAs.")

    with open(OUTPUT_TXT_PATH, 'w') as f:
        for m in mirnas:
            f.write(m + "\n")

    print(f"Saved list to {OUTPUT_TXT_PATH}")
    print("You can upload this file directly to https://www.mirnet.ca/")


if __name__ == "__main__":
    main()
