"""
Scancodon - High-performance Codon port of Scanpy
"""
import sys
import os
import numpy as np
from anndata import AnnData
from scipy import sparse as sp_sparse

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

    def _log1p_numpy_inplace(self, target, base):
        if sp_sparse.issparse(target):
            np.log1p(target.data, out=target.data)
            if base is not None:
                target.data /= np.log(base)
            return target
        arr = np.asarray(target)
        np.log1p(arr, out=arr)
        if base is not None:
            arr /= np.log(base)
        if arr is not target and hasattr(target, "__setitem__"):
            target[...] = arr
        return arr

    def _log1p_chunked_numpy(self, target, base, chunk_size):
        n_obs = target.shape[0]
        step = chunk_size or min(1000, n_obs) or 1
        for start in range(0, n_obs, step):
            stop = min(n_obs, start + step)
            if isinstance(target, np.ndarray):
                block = target[start:stop]
                self._log1p_numpy_inplace(block, base)
            else:
                block = np.asarray(target[start:stop])
                self._log1p_numpy_inplace(block, base)
                target[start:stop] = block

    def log1p(self, data, copy=False, chunked=False, chunk_size=None, base=None, **kwargs):
        adata = data.copy() if copy else data
        X = self._get_x(adata)

        if sp_sparse.issparse(X):
            self._log1p_numpy_inplace(X, base)
            return adata if copy else None

        is_backed = isinstance(adata, AnnData) and getattr(adata, "isbacked", False)
        require_chunked = chunked or chunk_size is not None or is_backed

        if require_chunked:
            self._log1p_chunked_numpy(X, base, chunk_size)
            return adata if copy else None

        use_native = CODON_AVAILABLE and isinstance(X, np.ndarray)

        if use_native:
            X_new = scancodon_native.log1p(X, base)
        else:
            X_new = self._log1p_numpy_inplace(X, base)

        if isinstance(adata, AnnData):
            adata.X = X_new if use_native else X
            return adata if copy else None
        return X_new if use_native else X

    def normalize_total(self, data, target_sum=None, inplace=True, **kwargs):
        if not inplace:
            data = data.copy()
        X = self._get_x(data)
        tgt = 1e4 if target_sum is None else float(target_sum)

        is_sparse = sp_sparse.issparse(X)
        use_native = CODON_AVAILABLE and isinstance(X, np.ndarray)

        if use_native:
            result, _ = scancodon_native.normalize_total(X, tgt)
        elif is_sparse:
            counts = np.asarray(X.sum(axis=1)).flatten()
            scales = tgt / np.maximum(counts, 1e-12)
            result = sp_sparse.diags(scales).dot(X)
        else:
            arr = np.asarray(X)
            counts = arr.sum(axis=1)
            scales = tgt / np.maximum(counts, 1e-12)
            arr = arr * scales[:, None]
            if arr is not X and hasattr(X, "__setitem__"):
                X[...] = arr
                result = X
            else:
                result = arr

        if isinstance(data, AnnData):
            data.X = result
            return data if not inplace else None
        return result

    def _scale_numpy(self, X, zero_center, max_value):
        arr = np.asarray(X, dtype=np.float64)
        if zero_center:
            arr = arr - arr.mean(axis=0)
        std = arr.std(axis=0, ddof=1)
        std[std == 0] = 1.0
        arr = arr / std
        if max_value is not None:
            if zero_center:
                arr = np.clip(arr, -max_value, max_value)
            else:
                arr = np.minimum(arr, max_value)
        if arr is not X and hasattr(X, "__setitem__"):
            X[...] = arr
            return X
        return arr

    def scale(self, data, zero_center=True, max_value=None, copy=False, **kwargs):
        adata = data.copy() if copy else data
        X = self._get_x(adata)
        use_native = CODON_AVAILABLE and isinstance(X, np.ndarray) and not sp_sparse.issparse(X)
        if use_native:
            X_new, _, _ = scancodon_native.scale(X, zero_center, max_value)
        else:
            X_new = self._scale_numpy(X, zero_center, max_value)
        if isinstance(adata, AnnData):
            adata.X = X_new
            return adata if copy else None
        return X_new

    def _filter_cells_numpy(self, X, min_counts, min_genes, max_counts, max_genes):
        dense = X.toarray() if sp_sparse.issparse(X) else np.asarray(X)
        if min_genes is not None or max_genes is not None:
            stats = (dense > 0).sum(axis=1)
        else:
            stats = dense.sum(axis=1)
        if min_counts is not None:
            mask = stats >= min_counts
        elif min_genes is not None:
            mask = stats >= min_genes
        elif max_counts is not None:
            mask = stats <= max_counts
        elif max_genes is not None:
            mask = stats <= max_genes
        else:
            mask = np.ones(dense.shape[0], dtype=bool)
        return np.asarray(mask, dtype=bool)

    def filter_cells(self, data, min_counts=None, min_genes=None, max_counts=None, max_genes=None, inplace=True, **kwargs):
        adata = data if inplace else data.copy()
        X = self._get_x(adata)
        use_native = CODON_AVAILABLE and isinstance(X, np.ndarray) and not sp_sparse.issparse(X)
        if use_native:
            mask, _ = scancodon_native.filter_cells(X, min_counts, min_genes, max_counts, max_genes)
        else:
            mask = self._filter_cells_numpy(X, min_counts, min_genes, max_counts, max_genes)
        adata._inplace_subset_obs(np.asarray(mask, dtype=bool))
        return None if inplace else (adata, mask)

    def _filter_genes_numpy(self, X, min_cells, min_counts, max_cells, max_counts):
        dense = X.toarray() if sp_sparse.issparse(X) else np.asarray(X)
        if min_cells is not None or max_cells is not None:
            stats = (dense > 0).sum(axis=0)
        else:
            stats = dense.sum(axis=0)
        if min_counts is not None:
            mask = stats >= min_counts
        elif min_cells is not None:
            mask = stats >= min_cells
        elif max_counts is not None:
            mask = stats <= max_counts
        elif max_cells is not None:
            mask = stats <= max_cells
        else:
            mask = np.ones(dense.shape[1], dtype=bool)
        return np.asarray(mask, dtype=bool)

    def _dist_matrix_from_knn(self, indices, distances, n_obs):
        rows: list[int] = []
        cols: list[int] = []
        data_vals: list[float] = []
        for i in range(n_obs):
            for j, idx in enumerate(indices[i]):
                rows.append(i)
                cols.append(int(idx))
                data_vals.append(float(distances[i, j]))
        return sp_sparse.csr_matrix((data_vals, (rows, cols)), shape=(n_obs, n_obs))

    def filter_genes(self, data, min_cells=None, min_counts=None, max_cells=None, max_counts=None, inplace=True, **kwargs):
        adata = data if inplace else data.copy()
        X = self._get_x(adata)
        use_native = CODON_AVAILABLE and isinstance(X, np.ndarray) and not sp_sparse.issparse(X)
        if use_native:
            mask, _ = scancodon_native.filter_genes(X, min_counts, min_cells, max_counts, max_cells)
        else:
            mask = self._filter_genes_numpy(X, min_cells, min_counts, max_cells, max_counts)
        adata._inplace_subset_var(np.asarray(mask, dtype=bool))
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

        if sp_sparse.issparse(X):
            data_matrix = X.toarray()
        else:
            data_matrix = np.asarray(X)

        use_native = CODON_AVAILABLE and isinstance(data_matrix, np.ndarray)

        if use_native:
            indices, distances, connectivities = scancodon_native.neighbors(data_matrix, n_neighbors)
            distances_matrix = self._dist_matrix_from_knn(indices, distances, data_matrix.shape[0])
        else:
            from sklearn.neighbors import NearestNeighbors
            nn = NearestNeighbors(n_neighbors=n_neighbors)
            nn.fit(data_matrix)
            distances_matrix = nn.kneighbors_graph(data_matrix, mode='distance')
            connectivities = nn.kneighbors_graph(data_matrix, mode='connectivity')
            indices = None
            distances = None

        adata.uns['neighbors'] = {
            'connectivities_key': 'connectivities',
            'distances_key': 'distances',
            'params': {
                'n_neighbors': n_neighbors,
                'method': kwargs.get('method', 'umap'),
                'n_pcs': n_pcs,
                'use_rep': use_rep,
            },
        }
        adata.obsp['connectivities'] = connectivities
        adata.obsp['distances'] = distances_matrix

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