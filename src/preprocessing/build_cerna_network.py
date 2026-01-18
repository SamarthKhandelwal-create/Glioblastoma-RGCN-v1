# -*- coding: utf-8 -*-
"""
build_cerna_network.py

GBM ceRNA prep pipeline (DE + enrichment):
1) Obtain expression data from UCSC Xena hubs:
   - Gene-level expected counts (TCGA + GTEx) from the Toil recompute hub
   - miRNA expression from the GDC hub (TCGA-GBM only on Xena)
2) Differential expression analysis with PyDESeq2:
   - Genes: TCGA-GBM tumor vs GTEx normal brain
   - miRNA: TCGA-GBM tumor vs TCGA-GBM solid tissue normal (since GTEx miRNA is
     not available on UCSC Xena public hubs)
3) Create DElncRNA, DEmRNA, DEmiRNA lists (padj < 0.05 by default) and save log2FC
4) Functional enrichment analysis (KEGG + GO BP) using GSEApy Enrichr wrapper
5) Output KEGG & GO enrichment graphs with the top 10 enriched terms

Outputs (default: ./results):
- DE_mRNA.tsv
- DE_lncRNA.tsv
- DE_miRNA.tsv
- enrichr_kegg.tsv, enrichr_go_bp.tsv
- KEGG_top10.png
- GO_BP_top10.png

Notes:
- The Toil hub provides RSEM expected counts (non-integer). For DESeq2, we round
  to integers.
- The TCGA miRNA Xena dataset appears log-scaled; we convert to pseudo-counts
  via inverse log2 transform for DESeq2 compatibility.
"""

from __future__ import annotations

import argparse
import os
import sys
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Iterable, List, Optional, Sequence, Tuple

if TYPE_CHECKING:  # pragma: no cover
import pandas as pd


# --- Xena configuration ---
TOIL_HOST = "https://toil.xenahubs.net"
TOIL_GENE_COUNTS_DATASET = "TcgaTargetGtex_gene_expected_count"
TOIL_PHENO_DATASET = "TcgaTargetGTEX_phenotype.txt"

GDC_HOST = "https://gdc.xenahubs.net"
GDC_MIRNA_DATASET = "TCGA-GBM.mirna.tsv"


def install_dependencies() -> None:
    """Best-effort dependency installer (useful for Colab / fresh environments)."""
    packages = [
        "numpy",
        "pandas",
        "xenaPython",
        "pydeseq2",
        "gseapy",
        "mygene",
        "matplotlib",
    ]
    for package in packages:
        try:
            __import__(package)
        except ImportError:
            print(f"[deps] Installing {package} ...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])


# -----------------------------
# Xena helpers
# -----------------------------

def _xena_dataset_samples(host: str, dataset: str) -> List[str]:
    import xenaPython.xenaQuery as xq

    return xq.dataset_samples(host, dataset)


def _xena_dataset_fields(host: str, dataset: str) -> List[str]:
    import xenaPython.xenaQuery as xq

    return xq.dataset_field(host, dataset)


def _xena_field_code_list(host: str, dataset: str, field: str) -> List[str]:
    """
    Return the code list for an enumerated Xena field.

    Xena stores many phenotype columns as integer codes; `xenaPython.field_codes`
    returns tab-separated label lists.
    """
    import xenaPython as xp

    rows = xp.field_codes(host, dataset, [field])
    for row in rows:
        if row.get("name") == field:
            return str(row.get("code", "")).split("\t")
    raise KeyError(f"Field codes not found for {dataset}:{field}")


def _xena_code_for_label(host: str, dataset: str, field: str, label: str) -> int:
    codes = _xena_field_code_list(host, dataset, field)
    try:
        return codes.index(label)
    except ValueError as e:
        raise ValueError(f"Label {label!r} not found for {dataset}:{field}") from e


def xena_fetch_matrix(
    host: str,
    dataset: str,
    samples: Sequence[str],
    fields: Sequence[str],
    *,
    chunk_size: int = 2000,
) -> "pd.DataFrame":
    """
    Fetch a Xena genomicMatrix (or phenotype matrix) as a DataFrame of shape:
      rows = fields, cols = samples
    """
    import xenaPython.xenaQuery as xq
    import pandas as pd

    frames: List[pd.DataFrame] = []
    fields_list = list(fields)
    samples_list = list(samples)

    for i in range(0, len(fields_list), chunk_size):
        chunk = fields_list[i : i + chunk_size]
        values = xq.dataset_probe_values(host, dataset, samples_list, chunk)  # rows=fields
        frames.append(pd.DataFrame(values, index=chunk, columns=samples_list))

    if not frames:
        return pd.DataFrame(index=list(fields_list), columns=list(samples_list))
    return pd.concat(frames, axis=0)


# -----------------------------
# Sample selection
# -----------------------------

def select_tcga_gbm_tumor_and_gtex_brain_normals(
    *,
    max_normal_samples: Optional[int] = None,
    random_seed: int = 0,
    gtex_brain_categories: Optional[Sequence[str]] = None,
) -> Tuple[List[str], List[str]]:
    """
    Select sample IDs for:
    - TCGA GBM tumor (primary tumor)
    - GTEx normal brain
    using the Toil phenotype dataset.
    """
    import xenaPython.xenaQuery as xq
    import pandas as pd

    pheno_samples = _xena_dataset_samples(TOIL_HOST, TOIL_PHENO_DATASET)
    fields = ["_study", "_primary_site", "_sample_type", "primary disease or tissue", "detailed_category"]
    pheno = xena_fetch_matrix(TOIL_HOST, TOIL_PHENO_DATASET, pheno_samples, fields, chunk_size=len(fields)).T
    # pheno is now samples x fields (numeric codes)
    pheno = pheno.apply(pd.to_numeric, errors="coerce")

    code_tcga = _xena_code_for_label(TOIL_HOST, TOIL_PHENO_DATASET, "_study", "TCGA")
    code_gtex = _xena_code_for_label(TOIL_HOST, TOIL_PHENO_DATASET, "_study", "GTEX")
    code_brain = _xena_code_for_label(TOIL_HOST, TOIL_PHENO_DATASET, "_primary_site", "Brain")
    code_gbm = _xena_code_for_label(
        TOIL_HOST, TOIL_PHENO_DATASET, "primary disease or tissue", "Glioblastoma Multiforme"
    )

    # "Primary Tumor" and "Primary Solid Tumor" both appear depending on encoding
    code_primary_tumor = _xena_code_for_label(TOIL_HOST, TOIL_PHENO_DATASET, "_sample_type", "Primary Tumor")
    code_primary_solid_tumor = _xena_code_for_label(
        TOIL_HOST, TOIL_PHENO_DATASET, "_sample_type", "Primary Solid Tumor"
    )
    primary_tumor_codes = {code_primary_tumor, code_primary_solid_tumor}

    mask_tumor = (
        (pheno["_study"] == code_tcga)
        & (pheno["primary disease or tissue"] == code_gbm)
        & (pheno["_sample_type"].isin(primary_tumor_codes))
    )
    mask_normal = (pheno["_study"] == code_gtex) & (pheno["_primary_site"] == code_brain)

    if gtex_brain_categories:
        detail_codes = _xena_field_code_list(TOIL_HOST, TOIL_PHENO_DATASET, "detailed_category")
        wanted: set[int] = set()
        for cat in gtex_brain_categories:
            try:
                wanted.add(detail_codes.index(cat))
            except ValueError:
                raise ValueError(f"Unknown GTEx brain category: {cat!r}")
        mask_normal = mask_normal & pheno["detailed_category"].isin(wanted)

    tumor_samples = pheno.index[mask_tumor].tolist()
    normal_samples = pheno.index[mask_normal].tolist()

    if max_normal_samples is not None and len(normal_samples) > max_normal_samples:
        normal_samples = (
            pd.Series(normal_samples).sample(n=max_normal_samples, random_state=random_seed).tolist()
        )

    # Intersect with the actual expression dataset sample IDs
    expr_samples = set(_xena_dataset_samples(TOIL_HOST, TOIL_GENE_COUNTS_DATASET))
    tumor_samples = [s for s in tumor_samples if s in expr_samples]
    normal_samples = [s for s in normal_samples if s in expr_samples]

    return tumor_samples, normal_samples


def split_tcga_tumor_normal(samples: Sequence[str]) -> Tuple[List[str], List[str]]:
    """Split TCGA sample IDs into tumor (01*) and normal (11*) using barcode sample-type."""
    tumor: List[str] = []
    normal: List[str] = []
    for s in samples:
        parts = s.split("-")
        if len(parts) < 4:
            continue
        code = parts[3]
        if code.startswith("01"):
            tumor.append(s)
        elif code.startswith("11"):
            normal.append(s)
    return tumor, normal


# -----------------------------
# Differential expression (PyDESeq2)
# -----------------------------

def run_pydeseq2(
    counts_df: "pd.DataFrame",
    metadata_df: "pd.DataFrame",
    *,
    contrast: Tuple[str, str, str] = ("condition", "tumor", "normal"),
    design_factor: str = "condition",
    n_cpus: int = 1,
    min_total_count: int = 10,
) -> "pd.DataFrame":
    """
    Run DESeq2 and return results_df (index = feature IDs).
    counts_df: samples x features
    metadata_df: samples x covariates (must align on index)
    """
    import numpy as np
    import pandas as pd
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    # Align indices
    metadata_df = metadata_df.loc[counts_df.index].copy()
    metadata_df[design_factor] = pd.Categorical(metadata_df[design_factor], categories=["normal", "tumor"])

    # Prefilter low-count features to speed up
    keep = counts_df.sum(axis=0) >= int(min_total_count)
    counts_df = counts_df.loc[:, keep]

    # Ensure non-negative integers
    counts_df = counts_df.fillna(0)
    counts_df = counts_df.clip(lower=0)
    if not np.issubdtype(counts_df.dtypes.iloc[0], np.integer):
        counts_df = counts_df.round().astype(int)

    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata_df,
        design_factors=design_factor,
        refit_cooks=True,
        n_cpus=int(n_cpus),
    )
    dds.deseq2()

    stat = DeseqStats(dds, contrast=contrast, n_cpus=int(n_cpus))
    stat.summary()

    # LFC shrink (best-effort; API can differ by version)
    try:
        stat.lfc_shrink()
    except Exception:
        pass

    return stat.results_df.copy()


# -----------------------------
# Annotation / splitting
# -----------------------------

def annotate_ensembl(
    ensembl_ids: Sequence[str],
    *,
    species: str = "human",
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Map Ensembl gene IDs -> symbol and -> type_of_gene using MyGene."""
    import mygene

    mg = mygene.MyGeneInfo()
    clean = [x.split(".")[0] for x in ensembl_ids]

    symbol: Dict[str, str] = {}
    gtype: Dict[str, str] = {}

    # Query in chunks to be gentle on the API
    chunk_size = 1000
    for i in range(0, len(clean), chunk_size):
        part = clean[i : i + chunk_size]
        res = mg.querymany(
            part,
            scopes="ensembl.gene",
            fields="symbol,type_of_gene",
            species=species,
            verbose=False,
        )
        for item in res:
            q = item.get("query")
            if not q:
                continue
            if "symbol" in item and item["symbol"] is not None:
                symbol[q] = item["symbol"]
            if "type_of_gene" in item and item["type_of_gene"] is not None:
                gtype[q] = item["type_of_gene"]

    return symbol, gtype


def classify_gene_type(type_of_gene: Optional[str]) -> str:
    # MyGene returns NaN (float) for missing values when loaded into pandas.
    if not isinstance(type_of_gene, str):
        return "other"
    s = type_of_gene.strip().lower()
    if s in {"protein-coding", "protein_coding", "protein coding"}:
        return "mRNA"
    # MyGene often uses broad labels like "ncRNA" for many long non-coding RNAs.
    # Treat ncRNA / lnc* as lncRNA, but keep small RNAs out.
    small_rna = {"mirna", "snorna", "snrna", "rrna", "trna", "scrna"}
    if s in small_rna:
        return "other"
    if "lnc" in s or s in {"ncrna", "nc-rna", "non-coding", "noncoding", "non_coding"}:
        return "lncRNA"
    return "other"


# -----------------------------
# Enrichment plotting
# -----------------------------

def plot_top10_enrichr(
    enrichr_df: "pd.DataFrame",
    *,
    title: str,
    outpath: Path,
    top_n: int = 10,
) -> None:
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    if enrichr_df is None or enrichr_df.empty:
        print(f"[enrich] No results for {title}; skipping plot.")
        return

    df = enrichr_df.copy()
    pcol = "Adjusted P-value" if "Adjusted P-value" in df.columns else "P-value"
    df = df.sort_values(pcol, ascending=True).head(int(top_n))
    df["-log10(p)"] = -np.log10(df[pcol].astype(float).clip(lower=1e-300))

    # Make plot height scale with term count
    h = max(4.0, 0.45 * len(df))
    plt.figure(figsize=(10, h))
    plt.barh(df["Term"][::-1], df["-log10(p)"][::-1])
    plt.xlabel(f"-log10({pcol})")
    plt.title(title)
    plt.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outpath, dpi=300)
    plt.close()


# -----------------------------
# Main pipeline
# -----------------------------

def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="GBM (TCGA) vs normal brain (GTEx) DE + enrichment pipeline.")
    parser.add_argument("--outdir", default="results", help="Output directory (default: results)")
    parser.add_argument("--pval-cutoff", type=float, default=0.05, help="Significance cutoff (default: 0.05)")
    parser.add_argument(
        "--pval-type",
        choices=["padj", "pvalue"],
        default="padj",
        help="Which p-value to filter on (default: padj)",
    )
    parser.add_argument(
        "--max-normal-samples",
        type=int,
        default=None,
        help="Optionally cap GTEx brain normal sample count (for speed).",
    )
    parser.add_argument(
        "--gtex-brain-category",
        action="append",
        default=None,
        help="Optionally restrict GTEx normals to a specific brain detailed category (repeatable).",
    )
    parser.add_argument("--min-total-count", type=int, default=10, help="Prefilter features by total count (default: 10)")
    parser.add_argument("--n-cpus", type=int, default=1, help="CPUs to use in PyDESeq2 (default: 1)")
    args = parser.parse_args(list(argv) if argv is not None else None)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=== Phase 0: Setup ===")
    install_dependencies()

    import numpy as np
    import pandas as pd
    import gseapy as gp

    # -------------------------
    # Genes: TCGA-GBM vs GTEx brain
    # -------------------------
    print("\n=== Phase 1: Fetch gene counts (TCGA-GBM tumor vs GTEx brain) ===")
    tumor_samples, normal_samples = select_tcga_gbm_tumor_and_gtex_brain_normals(
        max_normal_samples=args.max_normal_samples,
        random_seed=0,
        gtex_brain_categories=args.gtex_brain_category,
    )
    print(f"[samples] TCGA-GBM tumor: {len(tumor_samples)}")
    print(f"[samples] GTEx brain normal: {len(normal_samples)}")

    if len(tumor_samples) < 3 or len(normal_samples) < 3:
        raise RuntimeError("Not enough samples for DESeq2 (need >=3 per group).")

    all_samples = tumor_samples + normal_samples
    gene_fields = _xena_dataset_fields(TOIL_HOST, TOIL_GENE_COUNTS_DATASET)
    gene_fields = [g for g in gene_fields if g != "sampleID"]

    print(f"[xena] Fetching gene matrix: {len(gene_fields)} genes x {len(all_samples)} samples (chunked)...")
    gene_mat = xena_fetch_matrix(TOIL_HOST, TOIL_GENE_COUNTS_DATASET, all_samples, gene_fields, chunk_size=2000)
    # gene_mat: genes x samples -> transpose for PyDESeq2
    gene_counts = gene_mat.T
    gene_counts = gene_counts.apply(pd.to_numeric, errors="coerce").fillna(0)
    gene_counts = gene_counts.round().astype(int)

    meta_genes = pd.DataFrame(index=all_samples)
    meta_genes["condition"] = ["tumor"] * len(tumor_samples) + ["normal"] * len(normal_samples)

    print("\n=== Phase 2: Differential expression (genes) with PyDESeq2 ===")
    res_genes = run_pydeseq2(
        gene_counts,
        meta_genes,
        contrast=("condition", "tumor", "normal"),
        design_factor="condition",
        n_cpus=args.n_cpus,
        min_total_count=args.min_total_count,
    )

    # Filter significant
    pcol = "padj" if args.pval_type == "padj" else "pvalue"
    sig_genes = res_genes[(res_genes[pcol].notna()) & (res_genes[pcol] < float(args.pval_cutoff))].copy()
    print(f"[DE genes] significant ({pcol}<{args.pval_cutoff}): {len(sig_genes)}")

    # Annotate + split
    print("\n=== Phase 3: Annotate genes and write DE_mRNA / DE_lncRNA ===")
    ensembl_ids = sig_genes.index.tolist()
    symbol_map, type_map = annotate_ensembl(ensembl_ids)

    sig_genes = sig_genes.reset_index().rename(columns={"index": "ensembl_id"})
    sig_genes["ensembl_clean"] = sig_genes["ensembl_id"].astype(str).str.split(".").str[0]
    sig_genes["gene_symbol"] = sig_genes["ensembl_clean"].map(symbol_map)
    sig_genes["type_of_gene"] = sig_genes["ensembl_clean"].map(type_map)
    sig_genes["rna_class"] = sig_genes["type_of_gene"].apply(classify_gene_type)

    # Save full significant gene table (useful for debugging)
    sig_genes.to_csv(outdir / "DE_genes_significant.tsv", sep="\t", index=False)

    demrna = sig_genes[sig_genes["rna_class"] == "mRNA"].copy()
    delncrna = sig_genes[sig_genes["rna_class"] == "lncRNA"].copy()

    demrna.to_csv(outdir / "DE_mRNA.tsv", sep="\t", index=False)
    delncrna.to_csv(outdir / "DE_lncRNA.tsv", sep="\t", index=False)
    print(f"[export] DEmRNA: {len(demrna)} -> {outdir / 'DE_mRNA.tsv'}")
    print(f"[export] DElncRNA: {len(delncrna)} -> {outdir / 'DE_lncRNA.tsv'}")

    # -------------------------
    # miRNA: TCGA-GBM tumor vs TCGA-GBM normal
    # -------------------------
    print("\n=== Phase 4: Fetch miRNA matrix (TCGA-GBM) ===")
    mir_samples = _xena_dataset_samples(GDC_HOST, GDC_MIRNA_DATASET)
    mir_fields = _xena_dataset_fields(GDC_HOST, GDC_MIRNA_DATASET)
    mir_fields = [m for m in mir_fields if m != "sampleID"]

    tumor_mir, normal_mir = split_tcga_tumor_normal(mir_samples)
    print(f"[miRNA samples] tumor: {len(tumor_mir)} | normal: {len(normal_mir)} | total: {len(mir_samples)}")
    if len(tumor_mir) < 3 or len(normal_mir) < 3:
        print("[miRNA] WARNING: not enough tumor/normal miRNA samples for robust DESeq2; proceeding anyway.")

    mir_all = tumor_mir + normal_mir
    mir_mat = xena_fetch_matrix(GDC_HOST, GDC_MIRNA_DATASET, mir_all, mir_fields, chunk_size=2000)
    mir_df = mir_mat.T
    mir_df = mir_df.apply(pd.to_numeric, errors="coerce").fillna(0)

    # Heuristic: miRNA matrix is likely log2(x+1) scale on Xena; convert to pseudo-counts.
    vmax = float(np.nanmax(mir_df.values)) if mir_df.size else 0.0
    if vmax < 50.0:
        mir_counts = (np.power(2.0, mir_df) - 1.0).clip(lower=0)
    else:
        mir_counts = mir_df.clip(lower=0)
    mir_counts = mir_counts.round().astype(int)

    meta_mir = pd.DataFrame(index=mir_all)
    meta_mir["condition"] = ["tumor"] * len(tumor_mir) + ["normal"] * len(normal_mir)

    print("\n=== Phase 5: Differential expression (miRNA) with PyDESeq2 ===")
    res_mir = run_pydeseq2(
        mir_counts,
        meta_mir,
        contrast=("condition", "tumor", "normal"),
        design_factor="condition",
        n_cpus=args.n_cpus,
        min_total_count=max(1, args.min_total_count),
    )
    sig_mir = res_mir[(res_mir[pcol].notna()) & (res_mir[pcol] < float(args.pval_cutoff))].copy()
    sig_mir = sig_mir.reset_index().rename(columns={"index": "mirna_id"})
    sig_mir.to_csv(outdir / "DE_miRNA.tsv", sep="\t", index=False)
    print(f"[export] DEmiRNA: {len(sig_mir)} -> {outdir / 'DE_miRNA.tsv'}")

    # -------------------------
    # Enrichment (Enrichr via GSEApy)
    # -------------------------
    print("\n=== Phase 6: Enrichment (KEGG & GO BP) via Enrichr ===")
    gene_list = demrna["gene_symbol"].dropna().astype(str).unique().tolist()
    print(f"[enrich] Using {len(gene_list)} mRNA gene symbols.")

    if len(gene_list) == 0:
        print("[enrich] No mRNA genes available for enrichment; skipping.")
        return 0

    # KEGG
    enr_kegg = gp.enrichr(
        gene_list=gene_list,
        gene_sets="KEGG_2021_Human",
        organism="Human",
        outdir=None,
    )
    kegg_df = enr_kegg.results if enr_kegg is not None else None
    if kegg_df is not None and not kegg_df.empty:
        kegg_df.to_csv(outdir / "enrichr_kegg.tsv", sep="\t", index=False)
        plot_top10_enrichr(kegg_df, title="KEGG (top 10)", outpath=outdir / "KEGG_top10.png", top_n=10)
        print(f"[plot] {outdir / 'KEGG_top10.png'}")
    else:
        print("[enrich] No KEGG results.")

    # GO BP (try 2023 first; fallback to 2021 if missing)
    go_lib = "GO_Biological_Process_2023"
    try:
        enr_go = gp.enrichr(gene_list=gene_list, gene_sets=go_lib, organism="Human", outdir=None)
    except Exception:
        go_lib = "GO_Biological_Process_2021"
        enr_go = gp.enrichr(gene_list=gene_list, gene_sets=go_lib, organism="Human", outdir=None)

    go_df = enr_go.results if enr_go is not None else None
    if go_df is not None and not go_df.empty:
        go_df.to_csv(outdir / "enrichr_go_bp.tsv", sep="\t", index=False)
        plot_top10_enrichr(go_df, title=f"{go_lib} (top 10)", outpath=outdir / "GO_BP_top10.png", top_n=10)
        print(f"[plot] {outdir / 'GO_BP_top10.png'}")
            else:
        print("[enrich] No GO BP results.")
            
    print("\n=== Done ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


