import json
from pathlib import Path
from collections import defaultdict
import os
import sys

import torch
import pandas as pd
import numpy as np
from torch_geometric.data import Data

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.id_mapping import get_mapper

RESULTS_DIR = Path(os.environ.get("RESULTS_DIR", "results"))
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

FEATURE_MATRIX_PATH = RESULTS_DIR / "node_features_matrix.csv"
INTERACTIONS_PATH = RESULTS_DIR / "interactions.csv"
OUTPUT_GRAPH_PATH = RESULTS_DIR / "hetero_graph_GBM.pt"
OUTPUT_MAPPING_PATH = RESULTS_DIR / "node_mapping.json"

DE_MRNA_PATH = RESULTS_DIR / "DE_mRNA_logFC_abs_gt1.tsv"
DE_LNCRNA_PATH = RESULTS_DIR / "DE_lncRNA_logFC_abs_gt1.tsv"
DE_MIRNA_PATH = RESULTS_DIR / "DE_miRNA_logFC_abs_gt1.tsv"


def normalize_id(gene_id):
    s = str(gene_id).strip()
    if s.upper().startswith("ENSG"):
        return s.split(".")[0]
    return s.lower()


def load_external_interactions(path, mirna_col='mirna', gene_col='gene_symbol', score_col='score', default_score=0.8):
    if not path.exists():
        return {}
    try:
        df = pd.read_csv(path)
        interactions = {}
        for _, row in df.iterrows():
            mirna = normalize_id(row[mirna_col])
            gene = normalize_id(row.get(gene_col, '') or row.get('gene_ensembl', ''))
            if gene:
                interactions[(mirna, gene)] = float(row.get(score_col, default_score))
        print(f"Loaded {len(interactions)} interactions from {path.name}")
        return interactions
    except Exception as e:
        print(f"Warning: Failed to load {path}: {e}")
        return {}


def main():
    print("Building heterogeneous graph...")
    
    if not FEATURE_MATRIX_PATH.exists():
        raise FileNotFoundError(f"{FEATURE_MATRIX_PATH} not found. Run preprocessing first.")

    node_features_df = pd.read_csv(FEATURE_MATRIX_PATH, index_col=0)
    node_features_df.index = [normalize_id(g) for g in node_features_df.index]
    print(f"Loaded feature matrix: {node_features_df.shape}")
    
    # Optional miRNA injection
    if os.environ.get("INJECT_MIRNA", "").lower() in ("1", "true", "yes"):
        if INTERACTIONS_PATH.exists():
            interactions_tmp = pd.read_csv(INTERACTIONS_PATH)
            mi_keys = set()
            for col in ['gene_A', 'gene_B', 'A', 'B']:
                if col in interactions_tmp.columns:
                    mi_keys.update(interactions_tmp[col].dropna().astype(str).apply(normalize_id))
            mirna_keys = {k for k in mi_keys if 'mir' in k.lower() or 'let' in k.lower()}
            added = sum(1 for m in mirna_keys if m not in node_features_df.index 
                       and not node_features_df.loc.__setitem__((m, slice(None)), 0))
            for m in mirna_keys:
                if m not in node_features_df.index:
                    node_features_df.loc[m] = 0
            print(f"Injected {len([m for m in mirna_keys if m not in node_features_df.index]) or 0} miRNA rows")

    # Load DE gene lists
    de_mrna = pd.read_csv(DE_MRNA_PATH, sep='\t')
    de_lncrna = pd.read_csv(DE_LNCRNA_PATH, sep='\t')
    de_mirna = pd.read_csv(DE_MIRNA_PATH, sep='\t')

    mrna_set = set(de_mrna['ensembl_id'].astype(str).apply(normalize_id))
    lncrna_set = set(de_lncrna['ensembl_id'].astype(str).apply(normalize_id))
    mirna_set = set(de_mirna['mirna_id'].astype(str).apply(normalize_id))

    master_node_list = node_features_df.index.tolist()
    gene_to_index = {gene: idx for idx, gene in enumerate(master_node_list)}
    id_map = {gene: gene for gene in master_node_list}
    for gene in master_node_list:
        id_map[gene.lower()] = gene

    with open(OUTPUT_MAPPING_PATH, 'w') as f:
        json.dump(gene_to_index, f, indent=2)

    x = torch.tensor(node_features_df.values, dtype=torch.float)
    
    # Assign node types: 0=mRNA, 1=lncRNA, 2=miRNA
    gene_type_map = {}
    for gene in master_node_list:
        if gene in mrna_set:
            gene_type_map[gene] = 0
        elif gene in lncrna_set:
            gene_type_map[gene] = 1
        elif gene in mirna_set:
            gene_type_map[gene] = 2
        else:
            gene_type_map[gene] = -1

    # Load interactions
    interactions_df = pd.read_csv(INTERACTIONS_PATH)
    targetscan = load_external_interactions(RESULTS_DIR / "targetscan_predictions.csv", default_score=0.5)
    encori = load_external_interactions(RESULTS_DIR / "encori_interactions.csv", default_score=0.85)
    
    id_mapper = get_mapper()
    edge_list, edge_types, edge_sources, edge_confidence = [], [], [], []
    unique_edges = set()

    def resolve_id(raw_id):
        key = normalize_id(raw_id)
        if key in id_map:
            return id_map[key]
        if key.lower() in id_map:
            return id_map[key.lower()]
        for k, v in id_map.items():
            if key.lower() in k.lower():
                return v
        return None

    def add_edge(src_idx, dst_idx, rel_type, source='curated', conf=1.0):
        key = (src_idx, dst_idx, rel_type)
        if key not in unique_edges:
            unique_edges.add(key)
            edge_list.append([src_idx, dst_idx])
            edge_types.append(rel_type)
            edge_sources.append(source)
            edge_confidence.append(conf)
            return True
        return False

    # Process curated interactions
    for _, row in interactions_df.iterrows():
        raw_a = str(row.get('gene_A', row.get('A', ''))).strip()
        raw_b = str(row.get('gene_B', row.get('B', ''))).strip()
        
        gene_a, gene_b = resolve_id(raw_a), resolve_id(raw_b)
        if not gene_a or not gene_b:
            continue
        if gene_a not in gene_to_index or gene_b not in gene_to_index:
            continue
            
        idx_a, idx_b = gene_to_index[gene_a], gene_to_index[gene_b]
        type_a, type_b = gene_type_map[gene_a], gene_type_map[gene_b]
        
        # miRNA -> mRNA (type 0)
        if type_a == 2 and type_b == 0:
            add_edge(idx_a, idx_b, 0)
        elif type_b == 2 and type_a == 0:
            add_edge(idx_b, idx_a, 0)
        # lncRNA -> miRNA (type 1)
        elif type_a == 1 and type_b == 2:
            add_edge(idx_a, idx_b, 1)
        elif type_b == 1 and type_a == 2:
            add_edge(idx_b, idx_a, 1)

    print(f"Curated edges: {len(edge_list)}")

    # Infer ceRNA edges from shared miRNA targeting
    mirna_targets = defaultdict(set)
    for i, (src, dst) in enumerate(edge_list):
        if edge_types[i] == 0:
            src_gene = master_node_list[src]
            if gene_type_map.get(src_gene) == 2:
                mirna_targets[src].add(dst)

    cerna_added = 0
    for mirna_idx, targets in mirna_targets.items():
        target_list = list(targets)
        for i, t1 in enumerate(target_list):
            for t2 in target_list[i+1:]:
                if add_edge(min(t1, t2), max(t1, t2), 2, 'ceRNA'):
                    cerna_added += 1

    print(f"ceRNA edges: {cerna_added}")

    # Add TargetScan predictions
    ts_added = 0
    for (mirna, gene), score in targetscan.items():
        if score < 0.6:
            continue
        gene_a, gene_b = resolve_id(mirna), resolve_id(gene)
        if not gene_a or not gene_b:
            if gene.upper() in id_mapper.symbol_to_ensembl:
                gene_b = id_mapper.symbol_to_ensembl[gene.upper()]
                if gene_b not in gene_to_index:
                    gene_b = None
        if gene_a and gene_b and gene_a in gene_to_index and gene_b in gene_to_index:
            idx_a, idx_b = gene_to_index[gene_a], gene_to_index[gene_b]
            if gene_type_map.get(gene_a) == 2 and gene_type_map.get(gene_b) == 0:
                if add_edge(idx_a, idx_b, 0, 'TargetScan', score):
                    ts_added += 1

    print(f"TargetScan edges: {ts_added}")

    # Add ENCORI interactions
    enc_added = 0
    for (mirna, gene), score in encori.items():
        gene_a, gene_b = resolve_id(mirna), resolve_id(gene)
        if not gene_b and gene.upper() in id_mapper.symbol_to_ensembl:
            gene_b = id_mapper.symbol_to_ensembl[gene.upper()]
        if gene_a and gene_b and gene_a in gene_to_index and gene_b in gene_to_index:
            idx_a, idx_b = gene_to_index[gene_a], gene_to_index[gene_b]
            if gene_type_map.get(gene_a) == 2 and gene_type_map.get(gene_b) == 0:
                if add_edge(idx_a, idx_b, 0, 'ENCORI', score):
                    enc_added += 1

    print(f"ENCORI edges: {enc_added}")
    print(f"Total edges: {len(edge_list)}")

    # Build and save graph
    if edge_list:
        edge_index = torch.tensor(edge_list, dtype=torch.long).t().contiguous()
        edge_attr = torch.tensor(edge_types, dtype=torch.long)
    else:
        edge_index = torch.empty((2, 0), dtype=torch.long)
        edge_attr = torch.empty((0,), dtype=torch.long)

    data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)
    torch.save(data, OUTPUT_GRAPH_PATH)
    print(f"Saved graph to {OUTPUT_GRAPH_PATH}")
    print(data)

    # Save edge metadata
    if edge_sources:
        pd.DataFrame({
            'source': [e[0] for e in edge_list],
            'target': [e[1] for e in edge_list],
            'relation': edge_types,
            'provenance': edge_sources,
            'confidence': edge_confidence
        }).to_csv(RESULTS_DIR / "edge_metadata.csv", index=False)


if __name__ == "__main__":
    main()
