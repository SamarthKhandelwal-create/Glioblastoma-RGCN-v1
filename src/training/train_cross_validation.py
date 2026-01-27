"""
Cross-validation training pipeline for graph-based link prediction.

Supports both homogeneous (Data) and heterogeneous (HeteroData) graphs.
For homogeneous: 10-fold CV on positive edges with per-fold negative sampling.
For heterogeneous: relation-aware CV with per-relation negative sampling.

Metrics computed: ROC-AUC, PR-AUC, accuracy, precision, recall, F1, specificity, NPV.
Early stopping with model checkpointing per fold.
"""

import argparse
import csv
import os
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F
import random
from sklearn.metrics import (
    auc,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
)
from sklearn.model_selection import KFold
from torch_geometric.data import Data, HeteroData
from tqdm import tqdm

# Add repo root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from src.model import RGCN_LinkPredictor, HeteroRGCN_LinkPredictor


def compute_metrics(y_true, y_pred, y_scores):
    """Compute all validation metrics.
    
    Args:
        y_true: Ground truth binary labels
        y_pred: Predicted binary labels (threshold at 0.5)
        y_scores: Prediction scores (probabilities)
    
    Returns:
        dict with ROC-AUC, PR-AUC, accuracy, precision, recall, F1, specificity, NPV
    """
    metrics = {}
    
    # ROC-AUC
    if len(np.unique(y_true)) > 1:
        metrics["roc_auc"] = roc_auc_score(y_true, y_scores)
    else:
        metrics["roc_auc"] = 0.0
    
    # PR-AUC
    precision, recall, _ = precision_recall_curve(y_true, y_scores)
    metrics["pr_auc"] = auc(recall, precision)
    
    # Accuracy
    metrics["accuracy"] = accuracy_score(y_true, y_pred)
    
    # Precision
    metrics["precision"] = precision_score(y_true, y_pred, zero_division=0)
    
    # Recall
    metrics["recall"] = recall_score(y_true, y_pred, zero_division=0)
    
    # F1
    metrics["f1"] = f1_score(y_true, y_pred, zero_division=0)
    
    # Specificity: TN / (TN + FP)
    tn = np.sum((y_true == 0) & (y_pred == 0))
    fp = np.sum((y_true == 0) & (y_pred == 1))
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    metrics["specificity"] = specificity
    
    # NPV: TN / (TN + FN)
    fn = np.sum((y_true == 1) & (y_pred == 0))
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0.0
    metrics["npv"] = npv
    
    return metrics


def negative_sample(edge_index, num_nodes, num_samples=None, global_pos_set=None):
    """Generate negative samples (non-edges) for training/evaluation.
    
    Args:
        edge_index: [2, num_edges] edge indices
        num_nodes: Total number of nodes
        num_samples: Number of negatives to sample (default: num_edges)
    
    Returns:
        [2, num_samples] negative edge indices
    """
    if num_samples is None:
        num_samples = edge_index.shape[1]
    
    # Set of positive edges for fast lookup (local + optional global)
    pos_edges = set(map(tuple, edge_index.t().cpu().numpy()))
    if global_pos_set is not None:
        # union to ensure negatives avoid ANY known positive
        pos_edges = pos_edges.union(global_pos_set)
    
    neg_samples = []
    attempts = 0
    max_attempts = num_samples * 100
    
    while len(neg_samples) < num_samples and attempts < max_attempts:
        src = np.random.randint(0, num_nodes)
        dst = np.random.randint(0, num_nodes)
        if src != dst and (src, dst) not in pos_edges:
            neg_samples.append([src, dst])
        attempts += 1
    
    return torch.tensor(neg_samples, dtype=torch.long).t().contiguous()


def train_epoch_homogeneous(model, data, train_edge_index, train_edge_label, optimizer):
    """Train one epoch on homogeneous graph."""
    model.train()
    optimizer.zero_grad()
    
    # Encode using provided graph for message passing (avoid leakage)
    encode_eidx = getattr(data, 'encode_edge_index', None)
    encode_eattr = getattr(data, 'encode_edge_attr', None)
    if encode_eidx is None:
        z = model.encode(data.x, data.edge_index, data.edge_attr)
    else:
        z = model.encode(data.x, encode_eidx, encode_eattr)
    
    # Decode (both positive and negative edges with labels)
    logits = model.decode(z, train_edge_index)
    loss = F.binary_cross_entropy_with_logits(logits, train_edge_label.float())
    
    loss.backward()
    optimizer.step()
    
    return loss.item()


def eval_epoch_homogeneous(model, data, edge_index, edge_label):
    """Evaluate on homogeneous graph."""
    model.eval()
    with torch.no_grad():
        encode_eidx = getattr(data, 'encode_edge_index', None)
        encode_eattr = getattr(data, 'encode_edge_attr', None)
        if encode_eidx is None:
            z = model.encode(data.x, data.edge_index, data.edge_attr)
        else:
            z = model.encode(data.x, encode_eidx, encode_eattr)
        logits = model.decode(z, edge_index)
        scores = torch.sigmoid(logits).cpu().numpy()
        preds = (scores > 0.5).astype(int)
        
    return scores, preds


def train_cross_validation_homogeneous(
    graph_path, outdir, folds=10, epochs=100, seed=42, device="cpu"
):
    """10-fold cross-validation for homogeneous graphs.
    
    Splits positive edges, performs per-fold negative sampling,
    uses RGCN_LinkPredictor, computes all metrics, early stopping.
    """
    print("=" * 70)
    print("HOMOGENEOUS GRAPH CROSS-VALIDATION")
    print("=" * 70)
    
    # Load graph
    data = torch.load(graph_path, weights_only=False)
    data = data.to(device)
    
    print(f"Graph loaded: {data}")
    print(f"Number of nodes: {data.num_nodes}")
    print(f"Number of edges: {data.num_edges}")
    
    np.random.seed(seed)
    torch.manual_seed(seed)
    random.seed(seed)
    try:
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
    except Exception:
        pass
    random.seed(seed)
    # Deterministic flags for reproducibility (may impact performance)
    try:
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
    except Exception:
        pass
    
    # Get positive edges
    pos_edges = data.edge_index
    # Global positive set (avoid sampling true positives across folds)
    global_pos_set = set(map(tuple, pos_edges.t().cpu().numpy()))
    num_pos = pos_edges.shape[1]
    
    print(f"Positive edges: {num_pos}")
    
    # Split positive edges for CV
    kfold = KFold(n_splits=folds, shuffle=True, random_state=seed)
    all_indices = np.arange(num_pos)
    
    fold_results = []
    all_fold_scores = []
    all_fold_labels = []
    
    for fold_idx, (train_idx, test_idx) in enumerate(kfold.split(all_indices)):
        print(f"\n--- Fold {fold_idx + 1}/{folds} ---")
        
        # Split edges
        train_pos_edges = pos_edges[:, train_idx]
        test_pos_edges = pos_edges[:, test_idx]
        
        # Generate negatives (avoid ANY global positive edges)
        train_neg_edges = negative_sample(train_pos_edges, data.num_nodes, len(train_idx), global_pos_set)
        test_neg_edges = negative_sample(test_pos_edges, data.num_nodes, len(test_idx), global_pos_set)
        
        # Combine and create labels
        train_edge_index = torch.cat([train_pos_edges, train_neg_edges], dim=1)
        train_edge_label = torch.cat([
            torch.ones(train_pos_edges.shape[1]),
            torch.zeros(train_neg_edges.shape[1])
        ])
        
        test_edge_index = torch.cat([test_pos_edges, test_neg_edges], dim=1)
        test_edge_label = torch.cat([
            torch.ones(test_pos_edges.shape[1]),
            torch.zeros(test_neg_edges.shape[1])
        ])
        
        # Permute
        perm = torch.randperm(train_edge_index.shape[1])
        train_edge_index = train_edge_index[:, perm]
        train_edge_label = train_edge_label[perm]
        
        train_edge_index = train_edge_index.to(device)
        train_edge_label = train_edge_label.to(device)
        test_edge_index = test_edge_index.to(device)
        test_edge_label = test_edge_label.to(device)

        # Prepare graph used for encoding (mask out test positive edges)
        train_graph_edge_index = pos_edges[:, train_idx].to(device)
        train_graph_edge_attr = None
        if getattr(data, 'edge_attr', None) is not None:
            train_graph_edge_attr = data.edge_attr[train_idx].to(device)
        # Attach masked encode edges to data temporarily
        data.encode_edge_index = train_graph_edge_index
        data.encode_edge_attr = train_graph_edge_attr
        
        # Model
        model = RGCN_LinkPredictor(
            num_nodes=data.num_nodes,
            num_features=data.x.shape[1],
            num_relations=int(data.edge_attr.max().item()) + 1 if data.edge_attr is not None else 1,
            hidden_channels=64,
            out_channels=64
        ).to(device)
        
        optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
        
        # Training loop with early stopping
        best_val_auc = 0
        patience = 10
        patience_counter = 0
        
        for epoch in range(epochs):
            loss = train_epoch_homogeneous(model, data, train_edge_index, train_edge_label, optimizer)
            
            if epoch % 10 == 0:
                val_scores, val_preds = eval_epoch_homogeneous(model, data, test_edge_index, test_edge_label)
                val_auc = roc_auc_score(test_edge_label.cpu().numpy(), val_scores)
                print(f"  Epoch {epoch:3d} | Loss: {loss:.4f} | Val AUC: {val_auc:.4f}")
                
                if val_auc > best_val_auc:
                    best_val_auc = val_auc
                    patience_counter = 0
                    # Save best model
                    model_path = Path(outdir) / f"fold_{fold_idx}_best_model.pt"
                    torch.save(model.state_dict(), model_path)
                else:
                    patience_counter += 1
                    if patience_counter >= patience:
                        print(f"  Early stopping at epoch {epoch}")
                        break
        
        # Evaluate on test set
        val_scores, val_preds = eval_epoch_homogeneous(model, data, test_edge_index, test_edge_label)
        test_label_np = test_edge_label.cpu().numpy()
        
        metrics = compute_metrics(test_label_np, val_preds, val_scores)
        fold_results.append(metrics)
        
        all_fold_scores.extend(val_scores)
        all_fold_labels.extend(test_label_np)

        # Cleanup temporary encode attributes
        if hasattr(data, 'encode_edge_index'):
            delattr(data, 'encode_edge_index') if hasattr(data, 'encode_edge_index') else None
        if hasattr(data, 'encode_edge_attr'):
            delattr(data, 'encode_edge_attr') if hasattr(data, 'encode_edge_attr') else None
        
        print(f"  ROC-AUC: {metrics['roc_auc']:.4f} | PR-AUC: {metrics['pr_auc']:.4f} | F1: {metrics['f1']:.4f}")
    
    # Average metrics
    avg_metrics = {k: np.mean([fold[k] for fold in fold_results]) for k in fold_results[0].keys()}
    std_metrics = {k: np.std([fold[k] for fold in fold_results]) for k in fold_results[0].keys()}
    
    print("\n" + "=" * 70)
    print("CROSS-VALIDATION SUMMARY")
    print("=" * 70)
    for metric_name in sorted(avg_metrics.keys()):
        print(f"{metric_name:15s}: {avg_metrics[metric_name]:.4f} ± {std_metrics[metric_name]:.4f}")
    
    # Save results
    results_csv = Path(outdir) / "fold_results.csv"
    with open(results_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fold_results[0].keys())
        writer.writeheader()
        writer.writerows(fold_results)
    
    print(f"\nResults saved to {results_csv}")
    
    return fold_results, avg_metrics, all_fold_scores, all_fold_labels


def train_cross_validation_heterogeneous(
    graph_path, outdir, folds=10, epochs=100, seed=42, device="cpu"
):
    """10-fold CV for heterogeneous graphs with relation-aware sampling.
    
    Collects positive edges per relation, skips reverse relations,
    performs relation-specific negative sampling, uses HeteroRGCN_LinkPredictor.
    """
    print("=" * 70)
    print("HETEROGENEOUS GRAPH CROSS-VALIDATION")
    print("=" * 70)
    
    # Load graph
    data = torch.load(graph_path, weights_only=False)
    if device != "cpu":
        data = data.to(device)
    
    print(f"Graph loaded: {data}")
    print(f"Node types: {data.metadata()[0]}")
    print(f"Edge types: {data.metadata()[1]}")
    
    np.random.seed(seed)
    torch.manual_seed(seed)
    
    # Collect positive edges per relation
    edge_dict = {}
    for edge_type, edge_index in data.edge_index_dict.items():
        num_edges = edge_index.shape[1]
        # Skip reverse relations (simple heuristic: skip if reverse already seen)
        reverse_edge = (edge_type[2], "reverse_" + edge_type[1], edge_type[0])
        if reverse_edge not in edge_dict:
            edge_dict[edge_type] = edge_index
            print(f"  {edge_type}: {num_edges} edges")

    # Global pos sets per relation to avoid sampling true positives
    global_pos_sets = {et: set(map(tuple, e.t().cpu().numpy())) for et, e in edge_dict.items()}
    
    print(f"Total relations (non-reverse): {len(edge_dict)}")
    
    # Simple: perform CV on all edges combined (relation-aware scoring)
    kfold = KFold(n_splits=folds, shuffle=True, random_state=seed)
    
    fold_results = []
    all_fold_scores = []
    all_fold_labels = []
    
    # Collect all edges with their types
    all_edges_list = []
    for edge_type, edge_index in edge_dict.items():
        for i in range(edge_index.shape[1]):
            all_edges_list.append((edge_type, i, edge_index[:, i]))
    
    all_indices = np.arange(len(all_edges_list))
    
    for fold_idx, (train_idx, test_idx) in enumerate(kfold.split(all_indices)):
        print(f"\n--- Fold {fold_idx + 1}/{folds} ---")
        
        # Build train/test edge dicts
        train_edge_index_dict = {}
        test_edge_index_dict = {}
        
        for edge_type in edge_dict:
            train_edge_index_dict[edge_type] = []
            test_edge_index_dict[edge_type] = []
        
        for idx in train_idx:
            edge_type, _, edge = all_edges_list[idx]
            train_edge_index_dict[edge_type].append(edge.unsqueeze(1))
        
        for idx in test_idx:
            edge_type, _, edge = all_edges_list[idx]
            test_edge_index_dict[edge_type].append(edge.unsqueeze(1))
        
        # Concatenate per type
        for edge_type in edge_dict:
            if train_edge_index_dict[edge_type]:
                train_edge_index_dict[edge_type] = torch.cat(train_edge_index_dict[edge_type], dim=1)
            else:
                train_edge_index_dict[edge_type] = torch.empty((2, 0), dtype=torch.long)
            
            if test_edge_index_dict[edge_type]:
                test_edge_index_dict[edge_type] = torch.cat(test_edge_index_dict[edge_type], dim=1)
            else:
                test_edge_index_dict[edge_type] = torch.empty((2, 0), dtype=torch.long)
        
        # Generate negative samples per relation
        train_neg_dict = {}
        test_neg_dict = {}
        
        for edge_type, edge_index in train_edge_index_dict.items():
            src_type, _, dst_type = edge_type
            src_nodes = data[src_type].num_nodes
            dst_nodes = data[dst_type].num_nodes
            num_negs = edge_index.shape[1] if edge_index.shape[1] > 0 else 10
            
            # Simple neg sampling (avoid any global positives for this relation)
            neg_edges = []
            pos_set = global_pos_sets.get(edge_type, set())

            while len(neg_edges) < num_negs:
                src = np.random.randint(0, src_nodes)
                dst = np.random.randint(0, dst_nodes)
                if (src, dst) not in pos_set:
                    neg_edges.append([src, dst])

            train_neg_dict[edge_type] = torch.tensor(neg_edges, dtype=torch.long).t().contiguous()
        
        for edge_type, edge_index in test_edge_index_dict.items():
            src_type, _, dst_type = edge_type
            src_nodes = data[src_type].num_nodes
            dst_nodes = data[dst_type].num_nodes
            num_negs = edge_index.shape[1] if edge_index.shape[1] > 0 else 10
            
            neg_edges = []
            pos_set = global_pos_sets.get(edge_type, set())

            while len(neg_edges) < num_negs:
                src = np.random.randint(0, src_nodes)
                dst = np.random.randint(0, dst_nodes)
                if (src, dst) not in pos_set:
                    neg_edges.append([src, dst])

            test_neg_dict[edge_type] = torch.tensor(neg_edges, dtype=torch.long).t().contiguous()
        
        # Model
        in_channels_dict = {ntype: data[ntype].x.shape[1] for ntype in data.node_types}
        model = HeteroRGCN_LinkPredictor(
            metadata=data.metadata(),
            in_channels_dict=in_channels_dict,
            hidden_channels=64,
            out_channels=64,
            num_layers=2
        )
        if device != "cpu":
            model = model.to(device)
        
        optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
        
        # Training
        best_val_auc = 0
        patience = 10
        patience_counter = 0
        
        for epoch in range(epochs):
            model.train()
            optimizer.zero_grad()
            
            # Forward using only training edges for message passing (avoid leakage)
            x_dict = {ntype: data[ntype].x.to(device) for ntype in data.node_types}
            train_edge_index_dict_device = {k: v.to(device) for k, v in train_edge_index_dict.items()}
            out_dict = model(x_dict, train_edge_index_dict_device)
            
            # Compute loss on training edges
            loss = 0.0
            for edge_type, edge_index in train_edge_index_dict.items():
                if edge_index.shape[1] == 0:
                    continue
                src_type, _, dst_type = edge_type
                src_idx, dst_idx = edge_index
                
                src_emb = out_dict[src_type][src_idx]
                dst_emb = out_dict[dst_type][dst_idx]
                pos_scores = (src_emb * dst_emb).sum(dim=1)
                
                # Negative scores
                neg_idx = train_neg_dict[edge_type]
                neg_src_idx, neg_dst_idx = neg_idx
                neg_src_emb = out_dict[src_type][neg_src_idx]
                neg_dst_emb = out_dict[dst_type][neg_dst_idx]
                neg_scores = (neg_src_emb * neg_dst_emb).sum(dim=1)
                
                # BCE loss
                pos_labels = torch.ones(pos_scores.shape[0])
                neg_labels = torch.zeros(neg_scores.shape[0])
                all_scores = torch.cat([pos_scores, neg_scores])
                all_labels = torch.cat([pos_labels, neg_labels]).to(all_scores.device)
                
                loss += F.binary_cross_entropy_with_logits(all_scores, all_labels)
            
            loss.backward()
            optimizer.step()
            
            # Evaluation
            if epoch % 10 == 0:
                model.eval()
                with torch.no_grad():
                    x_dict = {ntype: data[ntype].x.to(device) for ntype in data.node_types}
                    out_dict = model(x_dict, train_edge_index_dict_device)
                    
                    all_val_scores = []
                    all_val_labels = []
                    
                    for edge_type, edge_index in test_edge_index_dict.items():
                        if edge_index.shape[1] == 0:
                            continue
                        src_type, _, dst_type = edge_type
                        src_idx, dst_idx = edge_index
                        
                        src_emb = out_dict[src_type][src_idx]
                        dst_emb = out_dict[dst_type][dst_idx]
                        pos_scores = torch.sigmoid((src_emb * dst_emb).sum(dim=1)).cpu().numpy()
                        
                        neg_idx = test_neg_dict[edge_type]
                        neg_src_idx, neg_dst_idx = neg_idx
                        neg_src_emb = out_dict[src_type][neg_src_idx]
                        neg_dst_emb = out_dict[dst_type][neg_dst_idx]
                        neg_scores = torch.sigmoid((neg_src_emb * neg_dst_emb).sum(dim=1)).cpu().numpy()
                        
                        all_val_scores.extend(list(pos_scores) + list(neg_scores))
                        all_val_labels.extend([1] * len(pos_scores) + [0] * len(neg_scores))
                    
                    if len(all_val_labels) > 0 and len(np.unique(all_val_labels)) > 1:
                        val_auc = roc_auc_score(all_val_labels, all_val_scores)
                        print(f"  Epoch {epoch:3d} | Loss: {loss:.4f} | Val AUC: {val_auc:.4f}")
                        
                        if val_auc > best_val_auc:
                            best_val_auc = val_auc
                            patience_counter = 0
                            model_path = Path(outdir) / f"fold_{fold_idx}_best_model.pt"
                            torch.save(model.state_dict(), model_path)
                        else:
                            patience_counter += 1
                            if patience_counter >= patience:
                                print(f"  Early stopping at epoch {epoch}")
                                break
        
        # Final evaluation
        model.eval()
        with torch.no_grad():
            x_dict = {ntype: data[ntype].x.to(device) for ntype in data.node_types}
            out_dict = model(x_dict, train_edge_index_dict_device)
            
            all_val_scores = []
            all_val_labels = []
            
            for edge_type, edge_index in test_edge_index_dict.items():
                if edge_index.shape[1] == 0:
                    continue
                src_type, _, dst_type = edge_type
                src_idx, dst_idx = edge_index
                
                src_emb = out_dict[src_type][src_idx]
                dst_emb = out_dict[dst_type][dst_idx]
                pos_scores = torch.sigmoid((src_emb * dst_emb).sum(dim=1)).cpu().numpy()
                
                neg_idx = test_neg_dict[edge_type]
                neg_src_idx, neg_dst_idx = neg_idx
                neg_src_emb = out_dict[src_type][neg_src_idx]
                neg_dst_emb = out_dict[dst_type][neg_dst_idx]
                neg_scores = torch.sigmoid((neg_src_emb * neg_dst_emb).sum(dim=1)).cpu().numpy()
                
                all_val_scores.extend(list(pos_scores) + list(neg_scores))
                all_val_labels.extend([1] * len(pos_scores) + [0] * len(neg_scores))
            
            all_val_scores = np.array(all_val_scores)
            all_val_labels = np.array(all_val_labels)
            all_val_preds = (all_val_scores > 0.5).astype(int)
        
        metrics = compute_metrics(all_val_labels, all_val_preds, all_val_scores)
        fold_results.append(metrics)
        
        all_fold_scores.extend(all_val_scores)
        all_fold_labels.extend(all_val_labels)
        
        print(f"  ROC-AUC: {metrics['roc_auc']:.4f} | PR-AUC: {metrics['pr_auc']:.4f} | F1: {metrics['f1']:.4f}")
    
    # Summary
    avg_metrics = {k: np.mean([fold[k] for fold in fold_results]) for k in fold_results[0].keys()}
    std_metrics = {k: np.std([fold[k] for fold in fold_results]) for k in fold_results[0].keys()}
    
    print("\n" + "=" * 70)
    print("CROSS-VALIDATION SUMMARY")
    print("=" * 70)
    for metric_name in sorted(avg_metrics.keys()):
        print(f"{metric_name:15s}: {avg_metrics[metric_name]:.4f} ± {std_metrics[metric_name]:.4f}")
    
    # Save
    results_csv = Path(outdir) / "fold_results.csv"
    with open(results_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fold_results[0].keys())
        writer.writeheader()
        writer.writerows(fold_results)
    
    print(f"\nResults saved to {results_csv}")
    
    return fold_results, avg_metrics, all_fold_scores, all_fold_labels


def main():
    parser = argparse.ArgumentParser(description="Cross-validation training for link prediction")
    parser.add_argument("--graph-path", type=str, default="results/hetero_graph_GBM.pt",
                        help="Path to graph file (.pt)")
    parser.add_argument("--outdir", type=str, default="results/cv_results",
                        help="Output directory for results")
    parser.add_argument("--folds", type=int, default=10,
                        help="Number of CV folds")
    parser.add_argument("--epochs", type=int, default=100,
                        help="Epochs per fold")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed")
    parser.add_argument("--device", type=str, default="cpu",
                        help="Device: cpu or cuda")
    
    args = parser.parse_args()
    
    # Check graph exists
    if not Path(args.graph_path).exists():
        print(f"ERROR: Graph file not found at {args.graph_path}")
        print("\nTo generate the graph, run:")
        print(f"  cd {Path(args.graph_path).parent.parent.parent}")
        print(f"  python src/graph/build_graph.py")
        return 1
    
    # Create output dir
    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    
    # Load and detect graph type
    data = torch.load(args.graph_path, weights_only=False)
    
    if isinstance(data, HeteroData):
        print("Detected HeteroData graph. Using HeteroRGCN_LinkPredictor.")
        fold_results, avg_metrics, all_scores, all_labels = train_cross_validation_heterogeneous(
            args.graph_path, args.outdir, args.folds, args.epochs, args.seed, args.device
        )
    elif isinstance(data, Data):
        print("Detected homogeneous Data graph. Using RGCN_LinkPredictor.")
        fold_results, avg_metrics, all_scores, all_labels = train_cross_validation_homogeneous(
            args.graph_path, args.outdir, args.folds, args.epochs, args.seed, args.device
        )
    else:
        print(f"ERROR: Unknown graph type {type(data)}")
        return 1
    
    # Save scores for visualization
    scores_file = Path(args.outdir) / "all_fold_scores.npy"
    labels_file = Path(args.outdir) / "all_fold_labels.npy"
    np.save(scores_file, np.array(all_scores))
    np.save(labels_file, np.array(all_labels))
    print(f"Scores saved to {scores_file}")
    print(f"Labels saved to {labels_file}")
    
    return 0


if __name__ == "__main__":
    exit(main())
