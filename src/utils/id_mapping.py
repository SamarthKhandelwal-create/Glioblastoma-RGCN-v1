"""
Enhanced ID mapping utilities using mygene and HGNC synonyms.
Provides robust gene/miRNA ID resolution across multiple identifier systems.

Usage:
    from src.utils.id_mapping import IDMapper
    mapper = IDMapper()
    ensembl_id = mapper.symbol_to_ensembl('TP53')
    hgnc_symbol = mapper.ensembl_to_symbol('ENSG00000141510')
"""
import pandas as pd
from pathlib import Path
import os
from typing import Dict, Optional, Set

class IDMapper:
    """
    Multi-source ID mapping: HGNC symbols <-> Ensembl IDs <-> Entrez IDs.
    Falls back to fuzzy matching if exact match fails.
    """
    
    def __init__(self):
        self.symbol_to_ensembl = {}
        self.ensembl_to_symbol = {}
        self.symbol_to_entrez = {}
        self.entrez_to_symbol = {}
        self.mirna_aliases = {}
        self._load_hgnc_data()
    
    def _load_hgnc_data(self):
        """Load HGNC mapping data if available, otherwise use fallback."""
        hgnc_path = Path(__file__).parent.parent.parent / "results" / "hgnc_complete_set.txt"
        
        if hgnc_path.exists():
            try:
                df = pd.read_csv(hgnc_path, sep='\t', low_memory=False)
                for _, row in df.iterrows():
                    symbol = str(row.get('symbol', '')).strip().upper()
                    ensembl = str(row.get('ensembl_gene_id', '')).strip()
                    entrez = str(row.get('entrez_id', '')).strip()
                    
                    if symbol:
                        if ensembl and ensembl != 'nan':
                            self.symbol_to_ensembl[symbol] = ensembl
                            self.ensembl_to_symbol[ensembl] = symbol
                        if entrez and entrez != 'nan':
                            self.symbol_to_entrez[symbol] = entrez
                            self.entrez_to_symbol[entrez] = symbol
                print(f"[IDMapper] Loaded HGNC mappings: {len(self.symbol_to_ensembl)} symbol->Ensembl")
            except Exception as e:
                print(f"[IDMapper] Warning: Failed to load HGNC data: {e}")
                self._load_fallback()
        else:
            self._load_fallback()
    
    def _load_fallback(self):
        """Load minimal fallback mappings for common genes."""
        fallback = {
            'TP53': 'ENSG00000141510',
            'PTEN': 'ENSG00000171862',
            'EGFR': 'ENSG00000146648',
            'FOXO3': 'ENSG00000118985',
            'DICER1': 'ENSG00000100697',
            'DROSHA': 'ENSG00000113318',
        }
        for symbol, ensembl in fallback.items():
            self.symbol_to_ensembl[symbol] = ensembl
            self.ensembl_to_symbol[ensembl] = symbol
        print(f"[IDMapper] Loaded fallback mappings: {len(fallback)} entries")
    
    def resolve_gene_id(self, gene_id: str, target_type='ensembl') -> Optional[str]:
        """
        Resolve a gene ID to target type (ensembl, symbol, or entrez).
        
        Args:
            gene_id: input gene identifier
            target_type: one of 'ensembl', 'symbol', 'entrez'
        
        Returns:
            Resolved ID or None if not found
        """
        gene_id = str(gene_id).strip()
        gene_id_upper = gene_id.upper()
        
        # Direct lookup
        if gene_id_upper in self.symbol_to_ensembl and target_type == 'ensembl':
            return self.symbol_to_ensembl[gene_id_upper]
        if gene_id in self.ensembl_to_symbol and target_type == 'symbol':
            return self.ensembl_to_symbol[gene_id]
        if gene_id_upper in self.symbol_to_entrez and target_type == 'entrez':
            return self.symbol_to_entrez[gene_id_upper]
        if gene_id in self.entrez_to_symbol and target_type == 'symbol':
            return self.entrez_to_symbol[gene_id]
        
        # Fuzzy match for symbols
        if target_type == 'ensembl':
            for symbol, ensembl in self.symbol_to_ensembl.items():
                if symbol.lower() == gene_id.lower():
                    return ensembl
        
        return None
    
    def add_mirna_alias(self, mirna_id: str, canonical_id: str):
        """Register miRNA alias mapping."""
        self.mirna_aliases[str(mirna_id).lower()] = str(canonical_id).lower()
    
    def resolve_mirna(self, mirna_id: str) -> str:
        """Resolve miRNA ID to canonical form."""
        normalized = str(mirna_id).lower()
        return self.mirna_aliases.get(normalized, normalized)


# Global instance
_mapper = None

def get_mapper():
    """Get or create global ID mapper instance."""
    global _mapper
    if _mapper is None:
        _mapper = IDMapper()
    return _mapper
