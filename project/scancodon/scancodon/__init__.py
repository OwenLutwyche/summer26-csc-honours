"""
Scancodon - High-performance Codon port of Scanpy
"""
import sys
import os
import warnings
import numpy as np
from anndata import AnnData
from scipy import sparse as sp_sparse
import scanpy as sc


# 1. NATIVE EXTENSION IMPORT
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

try:
    import scancodon_native
    CODON_AVAILABLE = True
except Exception as e:
    print(f"SCANCODON load failed ({type(e).__name__}): {e}")
    CODON_AVAILABLE = False

# 3. PREPROCESSING
class Preprocessing:
    def _get_x(self, data):
        return data.X if isinstance(data, AnnData) else data

    def scrublet(
        self, 
        data, 
        sim_doublet_ratio=2.0, 
        expected_doublet_rate=0.1, 
        stdev_doublet_rate=0.02, 
        random_state=0, 
        inplace=True,
        **kwargs
    ):
        adata = data if inplace else data.copy()
        X = self._get_x(adata)
        print("[INFO] Running native scancodon Scrublet")
        
        # check sparse matrix
        if not sp_sparse.issparse(X):
            X = sp_sparse.csr_matrix(X)
        elif not isinstance(X, sp_sparse.csr_matrix):
            X = X.tocsr()

        if len(X.indptr) != X.shape[0] + 1:
            print(f"CRITICAL ERROR: Matrix has {X.shape[0]} rows but indptr has {len(X.indptr)} elements.")
            # Attempt to fix it
            X = X.tocsr() # Force a clean rebuild
            if len(X.indptr) != X.shape[0] + 1:
                raise ValueError("Matrix indptr is structurally invalid.")
            
        # force 64-bit strict types for the Codon bridge, they will be converted to a codon CSRMatrix later
        data_64 = np.asarray(X.data, dtype=np.float64)
        indices_64 = np.asarray(X.indices, dtype=np.int64)
        indptr_64 = np.asarray(X.indptr, dtype=np.int64)
        
        # native dispatcher call
        scores = scancodon_native.scrublet(
            data_64, 
            indices_64, 
            indptr_64, 
            X.shape[0], 
            X.shape[1],
            float(sim_doublet_ratio),
            float(expected_doublet_rate),
            float(stdev_doublet_rate),
            int(random_state)
        )

        return None if inplace else adata

# 5. EXPORT
pp = Preprocessing()

sys.modules[__name__ + '.pp'] = pp

__all__ = ['pp', 'AnnData']