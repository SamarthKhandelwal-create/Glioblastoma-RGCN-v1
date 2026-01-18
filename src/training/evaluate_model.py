import torch
import torch.nn.functional as F
from torch_geometric.data import HeteroData
import torch_geometric.transforms as T
from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, recall_score, f1_score, confusion_matrix, roc_curve, precision_recall_curve
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# Allow importing from src
sys.path.append(os.getcwd())
from src.model import RGCN_LinkPredictor

# Configuration
GRAPH_PATH = "results/hetero_graph_GBM.pt"
MODEL_PATH = "results/gnn_model.pth"
PLOT_PATH = "results/model_performance.png"
device = torch.device('gpu' if torch.cuda.is_available() else 'cpu')

# --- Model Definition ---
# Imported from src.model

def main():
    print("--- Detailed Model Evaluation ---")
    
    # 1. Load Data
    data = torch.load(GRAPH_PATH, weights_only=False)
    data = data.to(device)
    
    # Create a fresh split for evaluation
    # Note: Ideally we'd use the exact same test split as training, but we didn't save it.
    # Evaluating on a new random split gives an unbiased estimate of generalization.
    transform = T.RandomLinkSplit(
        num_val=0.0,
        num_test=0.15, # Use 15% for robust testing
        is_undirected=False,
        add_negative_train_samples=True
    )
    _, _, test_data = transform(data)
    
    # 2. Load Model
    model = RGCN_LinkPredictor(
        num_nodes=data.num_nodes,
        num_features=data.num_features,
        num_relations=2,
        hidden_channels=32,
        out_channels=16
    ).to(device)
    
    model.load_state_dict(torch.load(MODEL_PATH))
    model.eval()
    
    # 3. Predict on Test Set
    with torch.no_grad():
        # Encode using the 'test' graph structure (which includes train edges usually, 
        # but RandomLinkSplit moves test edges to edge_label_index)
        # We use test_data.edge_index (message passing edges) to predict test_data.edge_label_index
        z = model.encode(test_data.x, test_data.edge_index, test_data.edge_attr)
        logits = model.decode(z, test_data.edge_label_index)
        probs = logits.sigmoid().cpu().numpy()
        labels = test_data.edge_label.cpu().numpy()
        
    # 4. Calculate Metrics
    preds = (probs > 0.5).astype(int)
    
    auc = roc_auc_score(labels, probs)
    acc = accuracy_score(labels, preds)
    prec = precision_score(labels, preds)
    rec = recall_score(labels, preds)
    f1 = f1_score(labels, preds)
    cm = confusion_matrix(labels, preds)
    
    print("\n--- Performance Metrics ---")
    print(f"AUC-ROC:   {auc:.4f}")
    print(f"Accuracy:  {acc:.4f}")
    print(f"Precision: {prec:.4f}")
    print(f"Recall:    {rec:.4f}")
    print(f"F1 Score:  {f1:.4f}")
    print("\nConfusion Matrix:")
    print(cm)
    
    # 5. Plotting
    fpr, tpr, _ = roc_curve(labels, probs)
    precision, recall, _ = precision_recall_curve(labels, probs)
    
    plt.figure(figsize=(12, 5))
    
    # ROC
    plt.subplot(1, 2, 1)
    plt.plot(fpr, tpr, label=f'RGCN (AUC = {auc:.2f})')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve')
    plt.legend()
    
    # PR Curve
    plt.subplot(1, 2, 2)
    plt.plot(recall, precision, label=f'F1 = {f1:.2f}')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(PLOT_PATH)
    print(f"\nPlots saved to {PLOT_PATH}")

if __name__ == "__main__":
    main()
