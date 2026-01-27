#!/usr/bin/env python3
"""
Extract Top GBM Biomarker Candidates from Trained RGCN Model

This script:
1. Loads the trained model and graph
2. Computes node embeddings
3. Ranks genes by predicted regulatory importance
4. Outputs top biomarker candidates

Usage:
    python extract_biomarkers.py --model-path results/model.pt --top 10
"""

import argparse
import json
import sys
from pathlib import Path

import torch
import pandas as pd
import numpy as np

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from model import RGCN_LinkPredictor


def load_model_and_data(model_path, graph_path, device='cpu'):
    """Load trained model and graph data."""
    # Load graph
    data = torch.load(graph_path, weights_only=False)
    print(f"Loaded graph: {data}")
    
    # Determine number of relations
    num_relations = int(data.edge_attr.max().item()) + 1 if data.edge_attr is not None else 3
    
    # Initialize model
    model = RGCN_LinkPredictor(
        num_nodes=data.x.shape[0],
        in_channels=data.x.shape[1],
        hidden_channels=64,
        out_channels=32,
        num_relations=num_relations
    ).to(device)
    
    # Load trained weights
    if Path(model_path).exists():
        model.load_state_dict(torch.load(model_path, map_location=device))
        print(f"Loaded model weights from {model_path}")
    else:
        print(f"WARNING: Model weights not found at {model_path}")
        print("Using untrained model - results will be based on graph structure only")
    
    model.eval()
    return model, data


def compute_node_importance(model, data, device='cpu'):
    """Compute node importance scores based on embeddings and connectivity."""
    
    with torch.no_grad():
        # Get node embeddings
        z = model.encode(
            data.x.to(device),
            data.edge_index.to(device),
            data.edge_attr.to(device)
        )
    
    # Compute importance metrics
    embeddings = z.cpu().numpy()
    edge_index = data.edge_index.cpu().numpy()
    edge_attr = data.edge_attr.cpu().numpy()
    
    num_nodes = embeddings.shape[0]
    
    # 1. Embedding magnitude (learned importance)
    embedding_magnitude = np.linalg.norm(embeddings, axis=1)
    
    # 2. Degree centrality (connectivity)
    out_degree = np.zeros(num_nodes)
    in_degree = np.zeros(num_nodes)
    for i in range(edge_index.shape[1]):
        src, dst = edge_index[0, i], edge_index[1, i]
        out_degree[src] += 1
        in_degree[dst] += 1
    total_degree = out_degree + in_degree
    
    # 3. ceRNA participation (edges with type 2)
    cerna_participation = np.zeros(num_nodes)
    for i in range(edge_index.shape[1]):
        if edge_attr[i] == 2:  # ceRNA edge
            src, dst = edge_index[0, i], edge_index[1, i]
            cerna_participation[src] += 1
            cerna_participation[dst] += 1
    
    # 4. miRNA targeting (being targeted by many miRNAs)
    mirna_targeting = np.zeros(num_nodes)
    for i in range(edge_index.shape[1]):
        if edge_attr[i] == 0:  # miRNA->target edge
            dst = edge_index[1, i]
            mirna_targeting[dst] += 1
    
    # Combined importance score (weighted)
    importance = (
        0.3 * (embedding_magnitude / embedding_magnitude.max()) +
        0.2 * (total_degree / max(total_degree.max(), 1)) +
        0.3 * (cerna_participation / max(cerna_participation.max(), 1)) +
        0.2 * (mirna_targeting / max(mirna_targeting.max(), 1))
    )
    
    return {
        'importance': importance,
        'embedding_magnitude': embedding_magnitude,
        'total_degree': total_degree,
        'cerna_participation': cerna_participation,
        'mirna_targeting': mirna_targeting,
        'embeddings': embeddings
    }


def get_node_info(results_dir):
    """Load node mapping and type information."""
    
    # Load node mapping
    mapping_path = results_dir / 'node_mapping.json'
    with open(mapping_path, 'r') as f:
        gene_to_idx = json.load(f)
    
    idx_to_gene = {v: k for k, v in gene_to_idx.items()}
    
    # Load edge metadata to determine node types
    metadata_path = results_dir / 'edge_metadata.csv'
    if metadata_path.exists():
        metadata_df = pd.read_csv(metadata_path)
    else:
        metadata_df = None
    
    # Load feature matrix to get node types
    features_path = results_dir / 'node_features_matrix.csv'
    features_df = pd.read_csv(features_path, index_col=0)
    
    # Determine node types
    node_types = {}
    for gene_id in gene_to_idx.keys():
        if 'hsa-' in gene_id.lower():
            node_types[gene_id] = 'miRNA'
        elif gene_id.startswith('ENSG'):
            # Check if it's in lncRNA or mRNA based on features or external annotation
            # For now, we'll mark as 'gene' and refine later
            node_types[gene_id] = 'gene'
        else:
            node_types[gene_id] = 'unknown'
    
    return idx_to_gene, node_types


def rank_biomarker_candidates(scores, idx_to_gene, node_types, top_n=10, gene_type_filter=None):
    """Rank and return top biomarker candidates."""
    
    results = []
    for idx, gene_id in idx_to_gene.items():
        node_type = node_types.get(gene_id, 'unknown')
        
        # Filter by gene type if specified
        if gene_type_filter and node_type != gene_type_filter:
            continue
        
        results.append({
            'rank': 0,
            'gene_id': gene_id,
            'node_type': node_type,
            'importance_score': scores['importance'][idx],
            'embedding_magnitude': scores['embedding_magnitude'][idx],
            'total_degree': int(scores['total_degree'][idx]),
            'cerna_participation': int(scores['cerna_participation'][idx]),
            'mirna_targeting': int(scores['mirna_targeting'][idx])
        })
    
    # Sort by importance score
    results.sort(key=lambda x: x['importance_score'], reverse=True)
    
    # Assign ranks
    for i, r in enumerate(results):
        r['rank'] = i + 1
    
    return results[:top_n]


def main():
    parser = argparse.ArgumentParser(description='Extract top GBM biomarker candidates')
    parser.add_argument('--model-path', type=str, default='results/model.pt',
                        help='Path to trained model weights')
    parser.add_argument('--graph-path', type=str, default='results/hetero_graph_GBM.pt',
                        help='Path to graph file')
    parser.add_argument('--results-dir', type=str, default='results',
                        help='Results directory')
    parser.add_argument('--top', type=int, default=10,
                        help='Number of top candidates to show')
    parser.add_argument('--filter', type=str, default=None,
                        choices=['gene', 'miRNA', None],
                        help='Filter by node type')
    parser.add_argument('--output', type=str, default=None,
                        help='Output CSV file path')
    args = parser.parse_args()
    
    results_dir = Path(args.results_dir)
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    
    print("=" * 70)
    print("  GBM BIOMARKER CANDIDATE EXTRACTION")
    print("=" * 70)
    print()
    
    # Load model and data
    print("Loading model and graph...")
    model, data = load_model_and_data(args.model_path, args.graph_path, device)
    
    # Compute importance scores
    print("Computing node importance scores...")
    scores = compute_node_importance(model, data, device)
    
    # Get node info
    print("Loading node information...")
    idx_to_gene, node_types = get_node_info(results_dir)
    
    # Rank candidates
    print(f"\nRanking top {args.top} biomarker candidates...")
    if args.filter:
        print(f"  Filtering by node type: {args.filter}")
    
    top_candidates = rank_biomarker_candidates(
        scores, idx_to_gene, node_types, 
        top_n=args.top, 
        gene_type_filter=args.filter
    )
    
    # Display results
    print("\n" + "=" * 70)
    print("  TOP GBM BIOMARKER CANDIDATES")
    print("=" * 70)
    print()
    print(f"{'Rank':<6}{'Gene ID':<20}{'Type':<10}{'Score':<10}{'Degree':<8}{'ceRNA':<8}{'miRNA':<8}")
    print("-" * 70)
    
    for c in top_candidates:
        print(f"{c['rank']:<6}{c['gene_id']:<20}{c['node_type']:<10}"
              f"{c['importance_score']:.4f}    {c['total_degree']:<8}"
              f"{c['cerna_participation']:<8}{c['mirna_targeting']:<8}")
    
    print("-" * 70)
    print()
    
    # Interpretation
    print("INTERPRETATION:")
    print("-" * 70)
    print("• Importance Score: Combined metric (embedding + connectivity + ceRNA + targeting)")
    print("• Degree: Total number of regulatory connections")
    print("• ceRNA: Participation in competing endogenous RNA relationships")
    print("• miRNA: Number of miRNAs targeting this gene")
    print()
    print("High-scoring genes are central to the ceRNA regulatory network")
    print("and may represent key drivers or therapeutic targets in GBM.")
    print()
    
    # Save to CSV if requested
    if args.output:
        df = pd.DataFrame(top_candidates)
        df.to_csv(args.output, index=False)
        print(f"Results saved to: {args.output}")
    
    # Also save full rankings
    all_candidates = rank_biomarker_candidates(
        scores, idx_to_gene, node_types, 
        top_n=len(idx_to_gene)
    )
    full_output = results_dir / 'biomarker_rankings.csv'
    pd.DataFrame(all_candidates).to_csv(full_output, index=False)
    print(f"Full rankings saved to: {full_output}")
    
    return top_candidates


if __name__ == '__main__':
    main()
