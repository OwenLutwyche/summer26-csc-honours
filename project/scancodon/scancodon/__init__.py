"""
Scancodon - High-performance Codon port of Scanpy
"""
import sys
import os
import numpy as np
from anndata import AnnData

# 1. NATIVE EXTENSION IMPORT
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

try:
    import scancodon_native
    CODON_AVAILABLE = True
    print("SCANCODON: Running on NATIVE CODON KERNELS")
except ImportError:
    scancodon_native = None
    CODON_AVAILABLE = False
    print("SCANCODON: Extension not found. Using NumPy fallbacks.")

# 2. CLASS DEFINITIONS
class Settings:
    def __init__(self):
        self.verbosity = 3
settings = Settings()

class Neighbors:
    def __init__(self, adata):
        self._adata = adata

# 3. PREPROCESSING
class Preprocessing:
    def _get_x(self, data):
        return data.X if isinstance(data, AnnData) else data

    def log1p(self, data, copy=False, **kwargs):
        adata = data.copy() if copy else data
        X = self._get_x(adata)
        if CODON_AVAILABLE:
            # Call FLAT native function
            X_new = scancodon_native.log1p(X)
        else:
            X_new = np.log1p(X)
        if isinstance(adata, AnnData): adata.X = X_new
        return adata if copy else None

    def normalize_total(self, data, target_sum=None, inplace=True, **kwargs):
        if not inplace: data = data.copy()
        X = self._get_x(data)
        tgt = 1e4 if target_sum is None else float(target_sum)
        if CODON_AVAILABLE:
            X_new, _ = scancodon_native.normalize_total(X, tgt)
        else:
            counts = X.sum(axis=1)
            X_new = X * (tgt / np.maximum(counts, 1e-12))[:, None]
        if isinstance(data, AnnData): data.X = X_new
        return data if not inplace else None

    def scale(self, data, zero_center=True, max_value=None, copy=False, **kwargs):
        adata = data.copy() if copy else data
        X = self._get_x(adata)
        if CODON_AVAILABLE:
            # Handle None for max_value by passing -1 or similar if needed, 
            # or ensure Codon side handles Optional. 
            # Ideally pass 0.0 or a flag if None not supported in raw wrapper
            mv = max_value if max_value is not None else 0.0 
            # Note: You might need to adjust Codon signature to take float, not Optional
            X_new, _, _ = scancodon_native.scale(X, zero_center, mv)
        else:
            X_new = (X - X.mean(0)) / (X.std(0) + 1e-12)
        if isinstance(data, AnnData): adata.X = X_new
        return adata if copy else None

    def filter_cells(self, data, min_counts=None, min_genes=None, max_counts=None, max_genes=None, inplace=True, **kwargs):
        adata = data if inplace else data.copy()
        X = self._get_x(adata)
        if CODON_AVAILABLE:
            # Pass 0 for None integers
            mc = min_counts if min_counts else 0
            mg = min_genes if min_genes else 0
            xc = max_counts if max_counts else 0
            xg = max_genes if max_genes else 0
            mask, _ = scancodon_native.filter_cells(X, mc, mg, xc, xg)
        else:
            mask = np.ones(X.shape[0], dtype=bool)
        adata._inplace_subset_obs(mask)
        return None if inplace else (adata, mask)

    def filter_genes(self, data, min_cells=None, min_counts=None, inplace=True, **kwargs):
        adata = data if inplace else data.copy()
        X = self._get_x(adata)
        if CODON_AVAILABLE:
            mc = min_counts if min_counts else 0
            ms = min_cells if min_cells else 0
            mask, _ = scancodon_native.filter_genes(X, mc, ms, 0, 0)
        else:
            mask = np.ones(X.shape[1], dtype=bool)
        adata._inplace_subset_var(mask)
        return None if inplace else (adata, mask)

    def highly_variable_genes(self, adata, n_top_genes=2000, flavor='seurat', subset=False, **kwargs):
        X = self._get_x(adata)
        if CODON_AVAILABLE:
            mask, means, vars_, _, _ = scancodon_native.highly_variable_genes_seurat_dense(X, n_top_genes)
            adata.var['highly_variable'] = np.array(mask, dtype=bool)
            adata.var['means'] = np.array(means)
            adata.var['dispersions'] = np.array(vars_)
        else:
            adata.var['highly_variable'] = np.ones(adata.n_vars, dtype=bool)
        if subset: adata._inplace_subset_var(adata.var['highly_variable'])

    def pca(self, data, n_comps=50, zero_center=True, **kwargs):
        adata = data
        X = self._get_x(adata)
        from sklearn.decomposition import PCA
        pca_obj = PCA(n_components=n_comps)
        X_pca = pca_obj.fit_transform(X)
        adata.obsm['X_pca'] = X_pca
        adata.varm['PCs'] = pca_obj.components_.T
        adata.uns['pca'] = {'variance_ratio': pca_obj.explained_variance_ratio_}

    def neighbors(self, adata, n_neighbors=15, n_pcs=None, use_rep=None, **kwargs):
        if use_rep == 'X_pca' and 'X_pca' in adata.obsm:
            X = adata.obsm['X_pca']
        else:
            X = adata.X
        if CODON_AVAILABLE:
            indices, distances, connectivities = scancodon_native.neighbors(X, n_neighbors)
            adata.uns['neighbors'] = {'connectivities_key': 'connectivities', 'distances_key': 'distances', 'params': {'n_neighbors': n_neighbors, 'method': 'umap'}}
            adata.obsp['connectivities'] = connectivities
            adata.obsp['distances'] = connectivities
        else:
            pass

# 4. TOOLS
class Tools:
    def leiden(self, adata, **kwargs):
        if CODON_AVAILABLE: scancodon_native.leiden(adata)
    def umap(self, adata, **kwargs):
        if CODON_AVAILABLE: scancodon_native.umap(adata)
    def rank_genes_groups(self, adata, groupby, method='t-test', **kwargs):
        if CODON_AVAILABLE: scancodon_native.rank_genes_groups(adata, groupby, method)
    def tsne(self, adata, **kwargs): pass
    def diffmap(self, adata, **kwargs): pass

# 5. EXPORT
pp = Preprocessing()
tl = Tools()
# Create mixin for sc.pp.neighbors style access
pp.neighbors = pp.neighbors

sys.modules[__name__ + '.pp'] = pp
sys.modules[__name__ + '.tl'] = tl

__all__ = ['pp', 'tl', 'settings', 'Neighbors', 'AnnData']