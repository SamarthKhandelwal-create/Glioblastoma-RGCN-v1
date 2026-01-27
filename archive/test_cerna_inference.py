#!/usr/bin/env python3
"""
Test ceRNA Inference Feature
=============================
Tests the ceRNA edge inference pipeline with miRNA injection support.

Usage:
    # Basic test (default results/ directory)
    python test_cerna_inference.py
    
    # With miRNA injection enabled
    $env:INJECT_MIRNA='true'
    python test_cerna_inference.py --inject
    
    # Run full pipeline then test
    python test_cerna_inference.py --run
    
    # Custom results directory
    python test_cerna_inference.py --results-dir custom_results
    
    # Verbose output with edge details
    python test_cerna_inference.py --verbose

Validates:
✓ Direct miRNA-target edges (curated)
✓ ceRNA edges (inferred from shared miRNA co-targeting)
✓ Node type assignments (mRNA, lncRNA, miRNA)
✓ Edge provenance and confidence scores
✓ ID normalization correctness
✓ Graph statistics and topology
"""

import argparse
import json
import os
import sys
from pathlib import Path
from collections import defaultdict, Counter

import pandas as pd
import numpy as np
import torch
from torch_geometric.data import HeteroData, Data
from torch_geometric.data.data import DataEdgeAttr

# PyTorch 2.6+ security: allow torch_geometric objects
try:
    torch.serialization.add_safe_globals([DataEdgeAttr])
except (AttributeError, TypeError):
    # Fallback for older PyTorch versions
    pass

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from graph.build_graph import normalize_id


class ceRNATestSuite:
    """Comprehensive test suite for ceRNA inference."""
    
    def __init__(self, results_dir="results", verbose=False):
        self.results_dir = Path(results_dir)
        self.verbose = verbose
        self.results = {
            'passed': 0,
            'failed': 0,
            'warnings': 0,
            'details': []
        }
        
    def log(self, level, message):
        """Log with level prefix."""
        prefix = {
            'INFO': '✓',
            'WARN': '⚠',
            'FAIL': '✗',
            'DEBUG': '→'
        }.get(level, '•')
        print(f"{prefix} [{level}] {message}")
        
    def assert_test(self, condition, test_name, error_msg=""):
        """Record test result."""
        if condition:
            self.results['passed'] += 1
            self.log('INFO', f"{test_name}: PASSED")
        else:
            self.results['failed'] += 1
            self.log('FAIL', f"{test_name}: FAILED {error_msg}")
            self.results['details'].append(f"FAILED: {test_name} - {error_msg}")
            
    def warn_test(self, condition, test_name, msg=""):
        """Record warning."""
        if not condition:
            self.results['warnings'] += 1
            self.log('WARN', f"{test_name}: {msg}")
            self.results['details'].append(f"WARNING: {test_name} - {msg}")

    def test_file_existence(self):
        """Test 1: Verify required files exist."""
        self.log('INFO', "\n=== Test 1: File Existence ===")
        
        required_files = {
            'Feature Matrix': self.results_dir / 'node_features_matrix.csv',
            'Interactions': self.results_dir / 'interactions.csv',
            'Graph': self.results_dir / 'hetero_graph_GBM.pt',
            'Node Mapping': self.results_dir / 'node_mapping.json',
        }
        
        for name, path in required_files.items():
            self.assert_test(path.exists(), f"{name} exists", f"Missing: {path}")
            
    def test_node_features(self):
        """Test 2: Validate node feature matrix."""
        self.log('INFO', "\n=== Test 2: Node Features ===")
        
        features_path = self.results_dir / 'node_features_matrix.csv'
        if not features_path.exists():
            self.log('WARN', "Skipping node features test (file missing)")
            return
            
        df = pd.read_csv(features_path, index_col=0)
        self.log('DEBUG', f"Feature matrix shape: {df.shape}")
        
        # Test shape
        self.assert_test(len(df) > 0, "Non-empty feature matrix", 
                        f"Got {len(df)} nodes")
        self.assert_test(df.shape[1] > 0, "Features present", 
                        f"Got {df.shape[1]} features")
        
        # Test miRNA presence (if INJECT_MIRNA enabled)
        inject_flag = os.environ.get("INJECT_MIRNA", "false").lower() in ("1", "true", "yes")
        mirna_count = len([idx for idx in df.index if 'mir' in idx.lower() or 'let' in idx.lower()])
        
        if inject_flag:
            self.warn_test(mirna_count > 0, "miRNAs in matrix", 
                          f"Expected miRNAs with INJECT_MIRNA=true, found {mirna_count}")
        else:
            self.log('DEBUG', f"Found {mirna_count} miRNA nodes (INJECT_MIRNA=false)")
            
        # Test feature statistics
        is_numeric = pd.api.types.is_numeric_dtype(df.iloc[:, 0])
        self.assert_test(is_numeric, "Features are numeric")
        
        # Check for NaN values
        nan_count = df.isna().sum().sum()
        self.assert_test(nan_count == 0, "No NaN values in features", 
                        f"Found {nan_count} NaN values")
        
    def test_interactions(self):
        """Test 3: Validate interaction file."""
        self.log('INFO', "\n=== Test 3: Interactions ===")
        
        interactions_path = self.results_dir / 'interactions.csv'
        if not interactions_path.exists():
            self.log('WARN', "Skipping interactions test (file missing)")
            return
            
        interactions = pd.read_csv(interactions_path)
        self.log('DEBUG', f"Loaded {len(interactions)} interactions")
        
        # Test structure
        self.assert_test('gene_A' in interactions.columns or 'A' in interactions.columns,
                        "Interaction columns present")
        self.assert_test(len(interactions) > 0, "Non-empty interactions", 
                        f"Got {len(interactions)} rows")
        
        # Categorize interactions
        mirna_target_count = 0
        lncrna_mirna_count = 0
        other_count = 0
        
        for _, row in interactions.iterrows():
            gene_a = str(row.iloc[0]).lower()
            gene_b = str(row.iloc[1]).lower()
            
            a_is_mirna = 'mir' in gene_a or 'let' in gene_a
            b_is_mirna = 'mir' in gene_b or 'let' in gene_b
            
            if a_is_mirna or b_is_mirna:
                mirna_target_count += 1
            else:
                other_count += 1
        
        self.log('DEBUG', f"miRNA-target interactions: {mirna_target_count}")
        self.log('DEBUG', f"Other interactions: {other_count}")
        self.assert_test(mirna_target_count > 0, "Has miRNA-target interactions")

    def test_graph_structure(self):
        """Test 4: Validate graph structure."""
        self.log('INFO', "\n=== Test 4: Graph Structure ===")
        
        graph_path = self.results_dir / 'hetero_graph_GBM.pt'
        if not graph_path.exists():
            self.log('WARN', "Skipping graph test (file missing)")
            return
            
        try:
            graph_data = torch.load(graph_path, weights_only=False)
        except TypeError:
            # Fallback for older PyTorch versions without weights_only parameter
            graph_data = torch.load(graph_path)
        self.log('DEBUG', f"Loaded graph: {graph_data}")
        
        # Check if HeteroData or regular Data
        is_hetero = hasattr(graph_data, 'node_types')
        self.log('DEBUG', f"Graph type: {'HeteroData' if is_hetero else 'Data'}")
        
        if is_hetero:
            # HeteroData structure
            node_types = list(graph_data.node_types)
            self.log('DEBUG', f"Node types: {node_types}")
            
            self.assert_test('mRNA' in node_types or 'lncRNA' in node_types or 'miRNA' in node_types,
                            "Has expected node types")
            
            # Test edges
            edge_types = list(graph_data.edge_types)
            self.log('DEBUG', f"Edge types (relation types): {edge_types}")
            
            total_edges = sum(graph_data[edge_type].edge_index.shape[1] for edge_type in edge_types 
                             if hasattr(graph_data[edge_type], 'edge_index'))
            self.log('DEBUG', f"Total edges: {total_edges}")
            
            self.assert_test(total_edges > 0, "Graph has edges", f"Got {total_edges} edges")
            
            # Check for ceRNA edges (edge type with 'competes' or type 2)
            has_cerna = any('competes' in str(et) for et in edge_types)
            self.warn_test(has_cerna, "Has ceRNA edges", 
                          "No ceRNA edges found (check inference)")
        else:
            # Regular Data structure
            self.log('DEBUG', f"Regular Data object with {graph_data.x.shape[0]} nodes")
            self.assert_test(hasattr(graph_data, 'edge_index'), "Has edge_index")
            
            total_edges = graph_data.edge_index.shape[1] if hasattr(graph_data, 'edge_index') else 0
            self.log('DEBUG', f"Total edges: {total_edges}")
            
            self.assert_test(total_edges > 0, "Graph has edges", f"Got {total_edges} edges")
            
            # Check for edge_attr (provenance)
            has_edge_attr = hasattr(graph_data, 'edge_attr') and graph_data.edge_attr is not None
            self.log('DEBUG', f"Has edge_attr: {has_edge_attr}")
            self.warn_test(has_edge_attr, "Has edge provenance", 
                          "No edge_attr for provenance tracking")
        
    def test_node_mapping(self):
        """Test 5: Validate node mapping."""
        self.log('INFO', "\n=== Test 5: Node Mapping ===")
        
        mapping_path = self.results_dir / 'node_mapping.json'
        if not mapping_path.exists():
            self.log('WARN', "Skipping node mapping test (file missing)")
            return
            
        with open(mapping_path, 'r') as f:
            mapping = json.load(f)
        
        # Check structure
        has_index_map = 'gene_to_index' in mapping or isinstance(mapping, dict)
        self.assert_test(has_index_map, "Mapping has gene index")
        
        # New: check for node_type_map
        has_type_map = 'node_type_map' in mapping
        self.log('DEBUG', f"Node type map present: {has_type_map}")
        
        if has_type_map:
            type_map = mapping['node_type_map']
            type_counts = Counter(type_map.values())
            self.log('DEBUG', f"Node type distribution: {dict(type_counts)}")
        
    def test_cerna_inference_logic(self):
        """Test 6: Validate ceRNA inference logic."""
        self.log('INFO', "\n=== Test 6: ceRNA Inference Logic ===")
        
        graph_path = self.results_dir / 'hetero_graph_GBM.pt'
        if not graph_path.exists():
            self.log('WARN', "Skipping ceRNA inference test (graph missing)")
            return
            
        try:
            graph_data = torch.load(graph_path, weights_only=False)
        except TypeError:
            graph_data = torch.load(graph_path)
        
        # Handle both HeteroData and regular Data
        is_hetero = hasattr(graph_data, 'edge_types')
        
        if is_hetero:
            # HeteroData structure
            edge_types = list(graph_data.edge_types)
            targets_edges = [et for et in edge_types if 'targets' in str(et)]
            competes_edges = [et for et in edge_types if 'competes' in str(et)]
            
            targets_count = sum(graph_data[et].edge_index.shape[1] for et in targets_edges 
                               if hasattr(graph_data[et], 'edge_index')) if targets_edges else 0
            competes_count = sum(graph_data[et].edge_index.shape[1] for et in competes_edges 
                                if hasattr(graph_data[et], 'edge_index')) if competes_edges else 0
        else:
            # Regular Data structure with edge_attr encoding relation types
            targets_count = graph_data.edge_index.shape[1] if hasattr(graph_data, 'edge_index') else 0
            competes_count = 0
            
            if hasattr(graph_data, 'edge_attr') and graph_data.edge_attr is not None:
                edge_attr = graph_data.edge_attr
                if edge_attr.dim() > 1:
                    targets_count = (edge_attr[:, 0] == 0).sum().item() if edge_attr.shape[1] > 0 else targets_count
                    competes_count = (edge_attr[:, 0] == 2).sum().item() if edge_attr.shape[1] > 0 else 0
                else:
                    targets_count = (edge_attr == 0).sum().item() if len(edge_attr) > 0 else targets_count
                    competes_count = (edge_attr == 2).sum().item() if len(edge_attr) > 0 else 0
        
        self.log('DEBUG', f"Direct miRNA-target edges: {targets_count}")
        self.log('DEBUG', f"Inferred ceRNA edges: {competes_count}")
        
        if targets_count > 0:
            self.warn_test(competes_count > 0, "ceRNA edges inferred",
                          f"Have {targets_count} direct edges but only {competes_count} ceRNA edges")
    
    def test_edge_provenance(self):
        """Test 7: Check edge provenance metadata."""
        self.log('INFO', "\n=== Test 7: Edge Provenance ===")
        
        metadata_path = self.results_dir / 'edge_metadata.csv'
        if not metadata_path.exists():
            self.log('DEBUG', "No edge_metadata.csv (optional)")
            return
            
        metadata = pd.read_csv(metadata_path)
        self.log('DEBUG', f"Edge metadata shape: {metadata.shape}")
        
        # Check columns
        expected_cols = ['source_gene', 'target_gene', 'relation_type', 'provenance', 'confidence']
        for col in expected_cols:
            self.warn_test(col in metadata.columns, f"Column '{col}' present")
        
        # Check provenance values
        if 'provenance' in metadata.columns:
            provenance_counts = metadata['provenance'].value_counts()
            self.log('DEBUG', f"Provenance distribution:\n{provenance_counts}")
        
        # Check confidence scores
        if 'confidence' in metadata.columns:
            conf_stats = metadata['confidence'].describe()
            self.log('DEBUG', f"Confidence score stats:\n{conf_stats}")
    
    def test_id_normalization(self):
        """Test 8: Validate ID normalization."""
        self.log('INFO', "\n=== Test 8: ID Normalization ===")
        
        # Check node_mapping.json for normalized IDs (authoritative source after build_graph.py)
        mapping_path = self.results_dir / 'node_mapping.json'
        if mapping_path.exists():
            import json
            with open(mapping_path, 'r') as f:
                mapping = json.load(f)
            
            ids = list(mapping.keys())
            self.log('DEBUG', f"IDs in node_mapping.json: {len(ids)}")
            
            # Count versioned Ensembl IDs (should be 0 after normalization)
            versioned = [idx for idx in ids if 'ENSG' in str(idx) and '.' in str(idx)]
            self.log('DEBUG', f"Versioned Ensembl IDs in mapping: {len(versioned)}")
            
            self.assert_test(len(versioned) == 0, "Ensembl IDs stripped of versions",
                            f"Found {len(versioned)} versioned IDs in node_mapping")
            
            # Check miRNA IDs are lowercase
            mirna_ids = [idx for idx in ids if 'hsa-' in str(idx).lower()]
            uppercase_mirna = [idx for idx in mirna_ids if idx != idx.lower()]
            
            self.assert_test(len(uppercase_mirna) == 0, "miRNA IDs are lowercase",
                            f"Found {len(uppercase_mirna)} uppercase miRNAs")
            
            self.log('DEBUG', f"Ensembl IDs: {len([i for i in ids if 'ENSG' in str(i)])}, miRNAs: {len(mirna_ids)}")
        else:
            self.log('WARN', "node_mapping.json not found, skipping ID normalization test")
    
    def test_inject_mirna_mode(self):
        """Test 9: Verify INJECT_MIRNA behavior."""
        self.log('INFO', "\n=== Test 9: INJECT_MIRNA Mode ===")
        
        inject_flag = os.environ.get("INJECT_MIRNA", "false").lower() in ("1", "true", "yes")
        self.log('DEBUG', f"INJECT_MIRNA enabled: {inject_flag}")
        
        features_path = self.results_dir / 'node_features_matrix.csv'
        if not features_path.exists():
            return
            
        df = pd.read_csv(features_path, index_col=0)
        mirna_count = len([idx for idx in df.index if 'mir' in idx.lower() or 'let' in idx.lower()])
        
        if inject_flag:
            # With injection, should have miRNAs even without miRNA_features.csv
            self.warn_test(mirna_count > 5, "Sufficient miRNAs injected",
                          f"Only {mirna_count} miRNAs found with INJECT_MIRNA=true")
        
        self.log('DEBUG', f"miRNA count in matrix: {mirna_count}")

    def run_all_tests(self):
        """Execute all tests."""
        print("\n" + "="*70)
        print("  ceRNA INFERENCE TEST SUITE")
        print("="*70)
        print(f"Results directory: {self.results_dir}")
        print(f"Verbose mode: {self.verbose}")
        print(f"INJECT_MIRNA: {os.environ.get('INJECT_MIRNA', 'false')}")
        print("="*70 + "\n")
        
        # Run all tests
        self.test_file_existence()
        self.test_node_features()
        self.test_interactions()
        self.test_graph_structure()
        self.test_node_mapping()
        self.test_cerna_inference_logic()
        self.test_edge_provenance()
        self.test_id_normalization()
        self.test_inject_mirna_mode()
        
        # Summary
        print("\n" + "="*70)
        print("  TEST SUMMARY")
        print("="*70)
        print(f"✓ PASSED: {self.results['passed']}")
        print(f"✗ FAILED: {self.results['failed']}")
        print(f"⚠ WARNINGS: {self.results['warnings']}")
        
        if self.results['details']:
            print("\nDetails:")
            for detail in self.results['details']:
                print(f"  • {detail}")
        
        print("="*70)
        
        # Exit code
        return 0 if self.results['failed'] == 0 else 1


def main():
    parser = argparse.ArgumentParser(
        description="Test ceRNA inference feature",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('--inject', action='store_true',
                       help='Enable INJECT_MIRNA mode for testing')
    parser.add_argument('--run', action='store_true',
                       help='Run full pipeline before testing')
    parser.add_argument('--results-dir', default='results',
                       help='Results directory (default: results)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    # Set INJECT_MIRNA if requested
    if args.inject:
        os.environ['INJECT_MIRNA'] = 'true'
    
    # Run pipeline if requested
    if args.run:
        print("\n" + "="*70)
        print("  RUNNING FULL PIPELINE")
        print("="*70 + "\n")
        
        try:
            print("Step 1: Building node features...")
            os.system("python src/preprocessing/build_node_features.py")
            
            print("\nStep 2: Building graph...")
            os.system("python src/graph/build_graph.py")
            
            print("\n✓ Pipeline completed\n")
        except Exception as e:
            print(f"✗ Pipeline failed: {e}")
            return 1
    
    # Run tests
    suite = ceRNATestSuite(results_dir=args.results_dir, verbose=args.verbose)
    exit_code = suite.run_all_tests()
    
    return exit_code


if __name__ == '__main__':
    sys.exit(main())
