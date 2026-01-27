import json
from pathlib import Path
import os
import sys

import torch
import pandas as pd
import numpy as np
from torch_geometric.data import Data

# Add utils to path
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.id_mapping import get_mapper

# Results directory via environment or default 'results'
RESULTS_DIR = Path(__import__("os").environ.get("RESULTS_DIR", "results"))
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

FEATURE_MATRIX_PATH = RESULTS_DIR / "node_features_matrix.csv"
INTERACTIONS_PATH = RESULTS_DIR / "interactions.csv"
OUTPUT_GRAPH_PATH = RESULTS_DIR / "hetero_graph_GBM.pt"
OUTPUT_MAPPING_PATH = RESULTS_DIR / "node_mapping.json"

# DE Files for Node Type Assignment
DE_MRNA_PATH = RESULTS_DIR / "DE_mRNA_logFC_abs_gt1.tsv"
DE_LNCRNA_PATH = RESULTS_DIR / "DE_lncRNA_logFC_abs_gt1.tsv"
DE_MIRNA_PATH = RESULTS_DIR / "DE_miRNA_logFC_abs_gt1.tsv"

def normalize_id(gene_id):
    """Strips Ensembl version numbers and lowercases miRNAs.

    Keeps Ensembl base ID (no version) and normalizes miRNA names to lower-case.
    """
    s = str(gene_id).strip()
    if s.upper().startswith("ENSG"):
        # Strip version: ENSG00000000001.15 -> ENSG00000000001
        base = s.split(".")[0]
        return base
    # Lowercase miRNAs: hsa-mir-X -> hsa-mir-x
    return s.lower()

def load_targetscan_predictions(targetscan_path):
    """Load TargetScan predictions if file exists.
    
    Args:
        targetscan_path: Path to targetscan_predictions.csv
    
    Returns:
        dict mapping (mirna, gene) -> confidence score, or empty dict if file missing
    """
    if not targetscan_path.exists():
        return {}
    
    try:
        df = pd.read_csv(targetscan_path)
        predictions = {}
        for _, row in df.iterrows():
            mirna = normalize_id(row['mirna'])
            gene = normalize_id(row['gene_symbol'])
            score = float(row.get('score', 0.5))
            predictions[(mirna, gene)] = score
        print(f"Loaded {len(predictions)} TargetScan predictions from {targetscan_path}")
        return predictions
    except Exception as e:
        print(f"Warning: Failed to load TargetScan predictions: {e}")
        return {}

def load_encori_interactions(encori_path):
    """Load ENCORI ceRNA interactions if file exists.
    
    Args:
        encori_path: Path to encori_interactions.csv
    
    Returns:
        dict mapping (mirna, gene) -> confidence score, or empty dict if file missing
    """
    if not encori_path.exists():
        return {}
    
    try:
        df = pd.read_csv(encori_path)
        interactions = {}
        for _, row in df.iterrows():
            mirna = normalize_id(row['mirna'])
            # Try gene_symbol first, then gene_ensembl
            gene = normalize_id(row.get('gene_symbol', '') or row.get('gene_ensembl', ''))
            if not gene:
                continue
            score = float(row.get('score', 0.8))  # ENCORI data is typically high confidence
            interactions[(mirna, gene)] = score
        print(f"Loaded {len(interactions)} ENCORI interactions from {encori_path}")
        return interactions
    except Exception as e:
        print(f"Warning: Failed to load ENCORI interactions: {e}")
        return {}



def main():
    print("--- Starting Phase 2: Heterogeneous Graph Construction (Robust ID Matching) ---")
    print("[build_graph] miRNA merge/injection available")
    
    # 1. Load Feature Matrix
    print(f"Loading feature matrix from {FEATURE_MATRIX_PATH}...")
    if not FEATURE_MATRIX_PATH.exists():
        raise FileNotFoundError(f"{FEATURE_MATRIX_PATH} not found. Run Phase 1 first.")

    node_features_df = pd.read_csv(FEATURE_MATRIX_PATH, index_col=0)
    print(f"Loaded Feature Matrix: Shape {node_features_df.shape}")
    
    # CRITICAL: Normalize all node IDs immediately after loading
    print("\n--- CRITICAL: Normalizing all node IDs ---")
    print(f"Before normalization: {len(node_features_df)} rows")
    
    # Count versioned IDs before normalization
    versioned_before = sum(1 for idx in node_features_df.index if '.' in str(idx) and 'ENSG' in str(idx))
    print(f"  Ensembl IDs with versions: {versioned_before}")
    
    # Apply normalization to index
    normalized_index = [normalize_id(gene_id) for gene_id in node_features_df.index]
    node_features_df.index = normalized_index
    
    # Verify normalization
    versioned_after = sum(1 for idx in node_features_df.index if '.' in str(idx))
    print(f"After normalization: {len(node_features_df)} rows, {versioned_after} with dots")
    print(f"  All IDs normalized: {versioned_after == 0}")

    # Optional: Merge miRNA-specific feature file if present
    mirna_features_path = RESULTS_DIR / "miRNA_features.csv"
    if mirna_features_path.exists():
        try:
            mir_df = pd.read_csv(mirna_features_path, index_col=0)
            print(f"Loaded miRNA features from {mirna_features_path}: Shape {mir_df.shape}")
            mir_df.index = mir_df.index.map(lambda s: normalize_id(s))
            if not mir_df.columns.equals(node_features_df.columns):
                mir_df = mir_df.reindex(columns=node_features_df.columns, fill_value=0)
            for idx, row in mir_df.iterrows():
                if idx in node_features_df.index:
                    continue
                node_features_df.loc[idx] = row
            print(f"After merging miRNA features, feature matrix shape: {node_features_df.shape}")
        except Exception as e:
            print(f"Warning: Failed to merge miRNA features: {e}")
    else:
        inject_flag = os.environ.get("INJECT_MIRNA", "false").lower() in ("1", "true", "yes")
        if inject_flag:
            print("INJECT_MIRNA enabled: will inject zero-valued feature rows for miRNAs found in interactions file if missing.")
            try:
                if INTERACTIONS_PATH.exists():
                    interactions_tmp = pd.read_csv(INTERACTIONS_PATH)
                    mi_candidate_keys = set()
                    for col in ['gene_A', 'gene_B', 'A', 'B']:
                        if col in interactions_tmp.columns:
                            mi_candidate_keys.update(interactions_tmp[col].dropna().astype(str).apply(normalize_id).tolist())
                    mirna_keys = {k for k in mi_candidate_keys if ('mir' in k.lower() or 'let' in k.lower())}
                    added = 0
                    for m in mirna_keys:
                        if m not in node_features_df.index:
                            node_features_df.loc[m] = np.zeros(node_features_df.shape[1])
                            added += 1
                    print(f"Injected {added} zero-feature miRNA rows.")
                else:
                    print(f"Interactions file {INTERACTIONS_PATH} not found; cannot infer miRNAs to inject.")
            except Exception as e:
                print(f"Warning during miRNA zero-injection: {e}")
    
    # 2. Load Gene Lists for Type Assignment & ID Mapping
    try:
        de_mrna = pd.read_csv(DE_MRNA_PATH, sep='\t')
        de_lncrna = pd.read_csv(DE_LNCRNA_PATH, sep='\t')
        de_mirna = pd.read_csv(DE_MIRNA_PATH, sep='\t')
    except FileNotFoundError as e:
        print(f"Critical Error: DE files not found. {e}")
        return

    # Create Sets and Version Maps
    matrix_ids = set(node_features_df.index)
    id_map = {}
    for mid in matrix_ids:
        clean = normalize_id(mid)
        id_map[clean] = mid
        id_map[mid] = mid

    print(f"Built ID Map with {len(id_map)} entries from {len(matrix_ids)} matrix nodes.")

    # Determine Types based on DE membership
    mrna_col = 'ensembl_id' if 'ensembl_id' in de_mrna.columns else de_mrna.columns[0]
    lncrna_col = 'ensembl_id' if 'ensembl_id' in de_lncrna.columns else de_lncrna.columns[0]
    mirna_col = 'mirna_id' if 'mirna_id' in de_mirna.columns else de_mirna.columns[0]

    mrna_clean_set = set(de_mrna[mrna_col].astype(str).apply(normalize_id))
    lncrna_clean_set = set(de_lncrna[lncrna_col].astype(str).apply(normalize_id))
    mirna_clean_set = set(de_mirna[mirna_col].astype(str).apply(lambda s: normalize_id(s).lower()))

    for k in list(id_map.keys()):
        id_map[k.lower()] = id_map[k]

# 3. Align Matrix - NORMALIZE ALL IDs
    # Strip Ensembl versions and lowercase miRNAs for consistency
    print("\n--- Normalizing node IDs ---")
    normalized_index = [normalize_id(gene) for gene in node_features_df.index]
    node_features_df.index = normalized_index
    
    master_node_list = node_features_df.index.tolist()
    gene_to_index = {gene: idx for idx, gene in enumerate(master_node_list)}
    
    # Rebuild id_map with normalized keys
    id_map.clear()
    for gene in master_node_list:
        id_map[gene] = gene  # Normalized to itself
    
    print(f"After normalization: {len(set(master_node_list))} unique nodes")

    with open(OUTPUT_MAPPING_PATH, 'w') as f:
        json.dump(gene_to_index, f, indent=4)
    
    x = torch.tensor(node_features_df.values, dtype=torch.float)
    print(f"Feature Tensor x shape: {x.shape}")
    
    # 4. Assign Node Types
    gene_type_map = {}
    type_counts = {0: 0, 1: 0, 2: 0, -1: 0}
    
    for gene in master_node_list:
        clean = normalize_id(gene)
        g_type = -1
        
        if clean in mrna_clean_set:
            g_type = 0
        elif clean in lncrna_clean_set:
            g_type = 1
        elif clean in mirna_clean_set:
            g_type = 2
        
        gene_type_map[gene] = g_type
        type_counts[g_type] += 1
        
    print(f"Node Type Assignment: mRNA={type_counts[0]}, lncRNA={type_counts[1]}, miRNA={type_counts[2]}, Unknown={type_counts[-1]}")

    # 5. Process Interactions with Fuzzy Matching
    print(f"Loading interactions from {INTERACTIONS_PATH}...")
    if not INTERACTIONS_PATH.exists():
        raise FileNotFoundError(f"Interactions file {INTERACTIONS_PATH} not found.")
    interactions_df = pd.read_csv(INTERACTIONS_PATH)
    
    # Load TargetScan predictions for confidence annotation
    targetscan_path = RESULTS_DIR / "targetscan_predictions.csv"
    targetscan_predictions = load_targetscan_predictions(targetscan_path)
    print(f"TargetScan predictions loaded: {len(targetscan_predictions)} entries")
    
    # Load ENCORI interactions (experimentally validated)
    encori_path = RESULTS_DIR / "encori_interactions.csv"
    encori_interactions = load_encori_interactions(encori_path)
    print(f"ENCORI interactions loaded: {len(encori_interactions)} entries")
    
    # Initialize ID mapper for enhanced resolution
    id_mapper = get_mapper()
    
    edge_list = []
    edge_types = []
    edge_sources = []  # Track provenance: 'curated', 'TargetScan', 'mixed'
    edge_confidence = []  # Track confidence scores
    unique_edges = set()
    
    processed = 0
    accepted = 0
    
    for _, row in interactions_df.iterrows():
        raw_a = str(row.get('gene_A', row.get('A', ''))).strip()
        raw_b = str(row.get('gene_B', row.get('B', ''))).strip()
        processed += 1
        
        a_key = normalize_id(raw_a)
        b_key = normalize_id(raw_b)
        gene_a = id_map.get(a_key) or id_map.get(a_key.lower())
        gene_b = id_map.get(b_key) or id_map.get(b_key.lower())

        if gene_a is None:
            for k, v in id_map.items():
                if a_key.lower() in k.lower():
                    gene_a = v
                    break
        if gene_b is None:
            for k, v in id_map.items():
                if b_key.lower() in k.lower():
                    gene_b = v
                    break
        
        if not gene_a or not gene_b:
            if processed < 5:
                print(f"Skipping edge {raw_a}-{raw_b}: resolved ({gene_a}, {gene_b})")
            continue
            
        if gene_a not in gene_to_index or gene_b not in gene_to_index:
            continue
            
        idx_a = gene_to_index[gene_a]
        idx_b = gene_to_index[gene_b]
        
        type_a = gene_type_map[gene_a]
        type_b = gene_type_map[gene_b]
        
        src, dst, rel = None, None, None
        
        # 0: miRNA -> mRNA (Silencing)
        if type_a == 2 and type_b == 0:
            src, dst, rel = idx_a, idx_b, 0
        elif type_b == 2 and type_a == 0:
            src, dst, rel = idx_b, idx_a, 0
            
        # 1: lncRNA -> miRNA (Sequestering)
        elif type_a == 1 and type_b == 2:
            src, dst, rel = idx_a, idx_b, 1
        elif type_b == 1 and type_a == 2:
            src, dst, rel = idx_b, idx_a, 1
            
        if src is not None:
            if (src, dst, rel) not in unique_edges:
                unique_edges.add((src, dst, rel))
                edge_list.append([src, dst])
                edge_types.append(rel)
                
                # Track provenance and confidence
                # Check if this edge appears in TargetScan predictions
                ts_key = (normalize_id(gene_a), normalize_id(gene_b))
                if ts_key in targetscan_predictions:
                    edge_sources.append('curated+TargetScan')
                    edge_confidence.append(targetscan_predictions[ts_key])
                else:
                    edge_sources.append('curated')
                    edge_confidence.append(1.0)  # Default confidence for curated interactions
                
                accepted += 1
                
    print(f"Processed {processed} rows. Accepted {accepted} valid edges.")
    
    # 5.5 INFER ceRNA EDGES from shared miRNA co-targeting
    print("\n--- Inferring ceRNA Edges (Shared miRNA Co-targeting) ---")
    from collections import defaultdict
    
    # Build map: miRNA index -> set of (target_index, target_type)
    mirna_targets = defaultdict(set)
    
    for edge_idx, (src, dst) in enumerate(edge_list):
        if edge_types[edge_idx] == 0:  # miRNA -> target edge (type 0)
            src_gene = master_node_list[src]
            dst_gene = master_node_list[dst]
            src_type = gene_type_map.get(src_gene, -1)
            dst_type = gene_type_map.get(dst_gene, -1)
            
            # Source must be miRNA (type 2)
            if src_type == 2:
                mirna_targets[src].add((dst, dst_type))
    
    print(f"  Found {len(mirna_targets)} miRNAs with targets")
    
    # Find genes sharing the same miRNA regulator
    cerna_pairs = defaultdict(set)  # (gene_a, gene_b) -> set of shared miRNAs
    
    for mirna_idx, targets in mirna_targets.items():
        target_list = list(targets)
        
        # Find all pairs of targets sharing this miRNA
        for i, (tar1_idx, tar1_type) in enumerate(target_list):
            for tar2_idx, tar2_type in target_list[i+1:]:
                pair_key = tuple(sorted([tar1_idx, tar2_idx]))
                cerna_pairs[pair_key].add(mirna_idx)
    
    print(f"  Found {len(cerna_pairs)} candidate ceRNA pairs")
    
    # Add ceRNA edges (type 2: co-regulation via shared miRNA)
    cerna_added = 0
    for (gene_a, gene_b), shared_mirnas in cerna_pairs.items():
        # Add ceRNA edge if pair shares at least one miRNA
        edge_tuple = (gene_a, gene_b, 2)
        if edge_tuple not in unique_edges:
            unique_edges.add(edge_tuple)
            edge_list.append([gene_a, gene_b])
            edge_types.append(2)
            cerna_added += 1
    
    print(f"  Added {cerna_added} ceRNA edges")
    print(f"  Total edges after ceRNA inference: {len(edge_list)}")
    
    # Update edge_sources and edge_confidence for ceRNA edges
    for _ in range(cerna_added):
        edge_sources.append("ceRNA")
        edge_confidence.append(1.0)
    
    # 5.5 INFER ceRNA EDGES (Shared miRNA co-targeting)
    print("\n--- Inferring ceRNA Edges (Shared miRNA Co-targeting) ---")
    
    # Build miRNA-target map from existing edges
    from collections import defaultdict
    mirna_targets = defaultdict(set)  # mirna_idx -> set of (target_idx, target_type)
    
    for i, (src_idx, dst_idx) in enumerate(edge_list):
        if edge_types[i] == 0:  # miRNA -> target edge
            src_gene = master_node_list[src_idx]
            dst_gene = master_node_list[dst_idx]
            src_type = gene_type_map.get(src_gene, -1)
            dst_type = gene_type_map.get(dst_gene, -1)
            
            # If source is miRNA
            if src_type == 2:
                mirna_targets[src_idx].add((dst_idx, dst_type))
    
    print(f"Found {len(mirna_targets)} miRNAs with targets")
    
    # Find co-regulated genes (genes sharing the same miRNA regulator)
    cerna_candidates = defaultdict(set)  # (gene_a, gene_b) -> set of shared miRNAs
    
    for mirna_idx, targets in mirna_targets.items():
        # Get all mRNAs (type 0)
        mrna_targets = [t_idx for t_idx, t_type in targets if t_type == 0]
        # Get all lncRNAs (type 1)
        lncrna_targets = [t_idx for t_idx, t_type in targets if t_type == 1]
        
        # mRNA-mRNA pairs sharing this miRNA
        for i, mrna1 in enumerate(mrna_targets):
            for mrna2 in mrna_targets[i+1:]:
                key = tuple(sorted([mrna1, mrna2]))
                cerna_candidates[key].add(mirna_idx)
        
        # lncRNA-lncRNA pairs
        for i, lncrna1 in enumerate(lncrna_targets):
            for lncrna2 in lncrna_targets[i+1:]:
                key = tuple(sorted([lncrna1, lncrna2]))
                cerna_candidates[key].add(mirna_idx)
        
        # lncRNA-mRNA pairs (main ceRNA topology)
        for lncrna_idx in lncrna_targets:
            for mrna_idx in mrna_targets:
                key = tuple(sorted([lncrna_idx, mrna_idx]))
                cerna_candidates[key].add(mirna_idx)
    
    # Add ceRNA edges (type 2)
    cerna_added = 0
    for (gene_a_idx, gene_b_idx), shared_mirnas in cerna_candidates.items():
        if len(shared_mirnas) > 0:
            if (gene_a_idx, gene_b_idx, 2) not in unique_edges and \
               (gene_b_idx, gene_a_idx, 2) not in unique_edges:
                unique_edges.add((gene_a_idx, gene_b_idx, 2))
                edge_list.append([gene_a_idx, gene_b_idx])
                edge_types.append(2)
                cerna_added += 1
    
    print(f"Inferred {cerna_added} ceRNA edges from shared miRNA co-targeting")
    print(f"Total edges after ceRNA inference: {len(edge_list)}")
    
    # 6. Add high-confidence TargetScan predictions not in curated interactions
    print("Adding high-confidence TargetScan predictions...")
    ts_only_added = 0
    ts_skipped = 0
    for (ts_mirna, ts_gene), ts_score in targetscan_predictions.items():
        if ts_score < 0.6:  # Lowered threshold to 0.6 to add more edges
            continue
        
        # Try to resolve to matrix IDs (multiple strategies)
        gene_a = None
        gene_b = None
        
        # Strategy 1: Direct ID map lookup
        gene_a = id_map.get(ts_mirna) or id_map.get(ts_mirna.lower())
        gene_b = id_map.get(ts_gene) or id_map.get(ts_gene.lower())
        
        # Strategy 2: Use ID mapper for gene symbol -> Ensembl conversion
        if gene_b is None and ts_gene.upper() in id_mapper.symbol_to_ensembl:
            resolved_ensembl = id_mapper.symbol_to_ensembl[ts_gene.upper()]
            if resolved_ensembl in gene_to_index:
                gene_b = resolved_ensembl
        
        # Strategy 3: Fuzzy substring matching
        if gene_a is None:
            for k, v in id_map.items():
                if ts_mirna.lower() in k.lower() and ('mir' in k.lower() or 'let' in k.lower()):
                    gene_a = v
                    break
        if gene_b is None:
            for k, v in id_map.items():
                if ts_gene.lower() in k.lower():
                    gene_b = v
                    break
        
        if not gene_a or not gene_b or gene_a not in gene_to_index or gene_b not in gene_to_index:
            ts_skipped += 1
            continue
        
        idx_a = gene_to_index[gene_a]
        idx_b = gene_to_index[gene_b]
        
        type_a = gene_type_map[gene_a]
        type_b = gene_type_map[gene_b]
        
        # Only create miRNA→mRNA edges from TargetScan
        if type_a == 2 and type_b == 0:
            edge_key = (idx_a, idx_b, 0)
            if edge_key not in unique_edges:
                unique_edges.add(edge_key)
                edge_list.append([idx_a, idx_b])
                edge_types.append(0)
                edge_sources.append('TargetScan')
                edge_confidence.append(ts_score)
                ts_only_added += 1
    
    print(f"Added {ts_only_added} TargetScan-only edges with high confidence (score >= 0.6).")
    print(f"Skipped {ts_skipped} TargetScan predictions (ID resolution failed).")
    
    # 6b. Add ENCORI validated interactions not in curated interactions
    print("Adding ENCORI-validated interactions...")
    encori_only_added = 0
    encori_skipped = 0
    for (enc_mirna, enc_gene), enc_score in encori_interactions.items():
        # Try to resolve to matrix IDs
        gene_a = None
        gene_b = None
        
        gene_a = id_map.get(enc_mirna) or id_map.get(enc_mirna.lower())
        gene_b = id_map.get(enc_gene) or id_map.get(enc_gene.lower())
        
        # Use ID mapper
        if gene_b is None and enc_gene.upper() in id_mapper.symbol_to_ensembl:
            resolved_ensembl = id_mapper.symbol_to_ensembl[enc_gene.upper()]
            if resolved_ensembl in gene_to_index:
                gene_b = resolved_ensembl
        
        # Fuzzy matching
        if gene_a is None:
            for k, v in id_map.items():
                if enc_mirna.lower() in k.lower() and ('mir' in k.lower() or 'let' in k.lower()):
                    gene_a = v
                    break
        if gene_b is None:
            for k, v in id_map.items():
                if enc_gene.lower() in k.lower():
                    gene_b = v
                    break
        
        if not gene_a or not gene_b or gene_a not in gene_to_index or gene_b not in gene_to_index:
            encori_skipped += 1
            continue
        
        idx_a = gene_to_index[gene_a]
        idx_b = gene_to_index[gene_b]
        
        type_a = gene_type_map[gene_a]
        type_b = gene_type_map[gene_b]
        
        # Only create miRNA→mRNA edges from ENCORI
        if type_a == 2 and type_b == 0:
            edge_key = (idx_a, idx_b, 0)
            if edge_key not in unique_edges:
                unique_edges.add(edge_key)
                edge_list.append([idx_a, idx_b])
                edge_types.append(0)
                edge_sources.append('ENCORI')
                edge_confidence.append(enc_score)
                encori_only_added += 1
            elif (idx_a, idx_b, 0) in unique_edges:
                # Edge already exists; mark as mixed source
                edge_idx = next((i for i, e in enumerate(edge_list) if e == [idx_a, idx_b]), None)
                if edge_idx is not None and edge_sources[edge_idx] != 'ENCORI+curated':
                    edge_sources[edge_idx] = 'ENCORI+' + edge_sources[edge_idx]
                    edge_confidence[edge_idx] = max(edge_confidence[edge_idx], enc_score)
    
    print(f"Added {encori_only_added} ENCORI-only edges.")
    print(f"Skipped {encori_skipped} ENCORI predictions (ID resolution failed).")
    
    # 7. Save Graph
    if edge_list:
        edge_index = torch.tensor(edge_list, dtype=torch.long).t().contiguous()
        edge_attr = torch.tensor(edge_types, dtype=torch.long)
    else:
        edge_index = torch.empty((2, 0), dtype=torch.long)
        edge_attr = torch.empty((0,), dtype=torch.long)

    data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)
    torch.save(data, OUTPUT_GRAPH_PATH)
    print(f"Saved Heterogeneous Graph to {OUTPUT_GRAPH_PATH}")
    print(data)
    
    # Save edge metadata (provenance and confidence)
    edge_metadata_path = RESULTS_DIR / "edge_metadata.csv"
    if edge_sources:
        metadata_df = pd.DataFrame({
            'edge_idx': range(len(edge_sources)),
            'source_node': [edge_list[i][0] for i in range(len(edge_list))],
            'target_node': [edge_list[i][1] for i in range(len(edge_list))],
            'relation_type': edge_types,
            'provenance': edge_sources,
            'confidence': edge_confidence
        })
        metadata_df.to_csv(edge_metadata_path, index=False)
        print(f"Saved edge metadata to {edge_metadata_path}")
        
        # Print provenance summary
        prov_counts = pd.Series(edge_sources).value_counts()
        print("\nEdge Provenance Summary:")
        for source, count in prov_counts.items():
            print(f"  {source}: {count} edges")

if __name__ == "__main__":
    main()
