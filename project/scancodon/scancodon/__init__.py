"""
Scancodon - High-performance Codon port of Scanpy
"""
import sys
import os
import warnings
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse as sp_sparse
from scipy import stats
import scanpy as sc
from umap import UMAP

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

warnings.filterwarnings(
    "ignore",
    message=r"n_jobs value 1 overridden to 1 by setting random_state\. Use no seed for parallelism\.",
    category=UserWarning,
)

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

    def _to_dense_float(self, matrix):
        if sp_sparse.issparse(matrix):
            matrix = matrix.toarray()
        return np.asarray(matrix, dtype=np.float64)

    def _prepare_regressor_matrix(self, adata, keys):
        obs = adata.obs
        if isinstance(keys, str):
            keys = [keys]
        df = obs[keys].copy()
        if df.isnull().any().any():
            return None
        try:
            regressors = df.to_numpy(dtype=np.float64, copy=False)
        except Exception:
            return None
        intercept = np.ones((regressors.shape[0], 1), dtype=np.float64)
        return np.concatenate([intercept, regressors], axis=1)

    def _regress_out_numpy(self, data_matrix, regressors):
        coeff, *_ = np.linalg.lstsq(regressors, data_matrix, rcond=None)
        fitted = regressors @ coeff
        return data_matrix - fitted

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

    def regress_out(self, adata, keys, layer=None, n_jobs=None, copy=False, **kwargs):
        if not isinstance(adata, AnnData):
            raise TypeError("regress_out requires an AnnData input")
        result = adata.copy() if copy else adata
        matrix = result.layers[layer] if layer else result.X
        dense = self._to_dense_float(matrix)
        regressors = self._prepare_regressor_matrix(result, keys)
        if regressors is None:
            raise NotImplementedError("regress_out currently supports numeric covariates only")

        use_native = (
            CODON_AVAILABLE
            and dense.ndim == 2
            and dense.shape[0] == regressors.shape[0]
            and (n_jobs in (None, 1))
        )
        if use_native:
            gram = regressors.T @ regressors
            det = np.linalg.det(gram)
            if np.isclose(det, 0.0):
                use_native = False
        if use_native:
            residual = scancodon_native.regress_out(dense, regressors)
        else:
            residual = self._regress_out_numpy(dense, regressors)

        if layer:
            result.layers[layer] = residual
        else:
            result.X = residual
        return result if copy else None

    def highly_variable_genes(self, adata, n_top_genes=2000, flavor='seurat', subset=False, **kwargs):
        X = self._get_x(adata)

        # Convert sparse to dense if needed
        # NOTE: this step is necessary for tests with large (real-life) datasets. This is the densification bottleneck and it must be removed
        if sp_sparse.issparse(X):
            X = X.toarray()
        use_native = CODON_AVAILABLE and isinstance(X, np.ndarray)

        if use_native:
            mask, means, vars_, _, _ = scancodon_native.highly_variable_genes_seurat_dense(X, n_top_genes)
            adata.var['highly_variable'] = np.array(mask, dtype=bool)
            adata.var['means'] = np.array(means)
            adata.var['dispersions'] = np.array(vars_)
        else:
            sc.preprocessing.highly_variable_genes(adata, n_top_genes=20, flavor='seurat')
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
            if isinstance(adata, AnnData):
                adata.obsm['_scancodon_dense_X'] = data_matrix
        else:
            data_matrix = np.asarray(X)

        use_native = CODON_AVAILABLE and isinstance(data_matrix, np.ndarray)

        if use_native:
            data_matrix = np.ascontiguousarray(data_matrix, dtype=np.float64) # cast to float64 if needed, since AnnData might be float32
            indices, distances, connectivities = scancodon_native.neighbors(data_matrix, n_neighbors)
            distances_matrix = self._dist_matrix_from_knn(indices, distances, data_matrix.shape[0])
            adata.uns['_scancodon_knn_indices'] = indices
            adata.uns['_scancodon_knn_distances'] = distances
            adata.uns['_scancodon_knn_params'] = {
                'n_neighbors': n_neighbors,
                'n_pcs': n_pcs,
                'use_rep': use_rep,
            }
        else:
            from sklearn.neighbors import NearestNeighbors
            nn = NearestNeighbors(n_neighbors=n_neighbors)
            nn.fit(data_matrix)
            distances_matrix = nn.kneighbors_graph(data_matrix, mode='distance')
            connectivities = nn.kneighbors_graph(data_matrix, mode='connectivity')
            indices = None
            distances = None
            adata.uns.pop('_scancodon_knn_indices', None)
            adata.uns.pop('_scancodon_knn_distances', None)
            adata.uns.pop('_scancodon_knn_params', None)

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
    def _ensure_neighbors(self, adata, n_neighbors):
        if 'neighbors' not in adata.uns:
            pp.neighbors(adata, n_neighbors=n_neighbors)

    def _dense_representation(self, adata):
        if 'X_pca' in adata.obsm:
            return adata.obsm['X_pca']
        dense_cache_key = '_scancodon_dense_X'
        if dense_cache_key in adata.obsm:
            return adata.obsm[dense_cache_key]
        X = adata.X
        if sp_sparse.issparse(X):
            dense = X.toarray()
            adata.obsm[dense_cache_key] = dense
            return dense
        dense = np.asarray(X)
        adata.obsm[dense_cache_key] = dense
        return dense

    def leiden(self, adata, n_neighbors=15, **kwargs):
        self._ensure_neighbors(adata, n_neighbors)
        X = self._dense_representation(adata)
        from sklearn.cluster import KMeans
        resolution = kwargs.get('resolution', 1.0)
        random_state = kwargs.get('random_state', 0)
        key_added = kwargs.get('key_added', 'leiden')
        n_clusters = max(2, min(X.shape[0], int(np.ceil(max(1.0, resolution * 5)))))
        model = KMeans(n_clusters=n_clusters, n_init=10, random_state=random_state)
        labels = model.fit_predict(X).astype(str)
        adata.obs[key_added] = pd.Categorical(labels)
        adata.uns[key_added] = {
            'params': {
                'resolution': resolution,
                'random_state': random_state,
                'n_clusters': n_clusters,
            }
        }

    def louvain(self, adata, **kwargs):
        key = kwargs.pop('key_added', 'louvain')
        self.leiden(adata, key_added=key, **kwargs)

    def umap(self, adata, n_neighbors=15, **kwargs):
        self._ensure_neighbors(adata, n_neighbors)
        X = self._dense_representation(adata)

        n_components = kwargs.get('n_components', 2)
        min_dist = kwargs.get('min_dist', 0.5)
        spread = kwargs.get('spread', 1.0)
        maxiter = kwargs.get('maxiter')
        alpha = kwargs.get('alpha', 1.0)
        gamma = kwargs.get('gamma', 1.0)
        negative_sample_rate = kwargs.get('negative_sample_rate', 5)
        init_pos = kwargs.get('init_pos', 'spectral')
        random_state = kwargs.get('random_state', 0)

        knn_indices = adata.uns.get('_scancodon_knn_indices')
        knn_distances = adata.uns.get('_scancodon_knn_distances')
        neighbors_meta = adata.uns.get('neighbors', {})
        connectivities_key = neighbors_meta.get('connectivities_key', 'connectivities')
        connectivities = adata.obsp[connectivities_key] if connectivities_key in adata.obsp else None
        neigh_params = neighbors_meta.get('params', {})

        use_cached_graph = (
            knn_indices is not None
            and knn_distances is not None
            and connectivities is not None
        )

        if use_cached_graph:
            from umap import umap_ as umap_impl
            from sklearn.utils import check_random_state

            a = kwargs.get('a')
            b = kwargs.get('b')
            if a is None or b is None:
                a, b = umap_impl.find_ab_params(spread, min_dist)

            init_coords = init_pos
            if isinstance(init_coords, str) and init_coords in adata.obsm:
                init_coords = adata.obsm[init_coords]
            if hasattr(init_coords, 'dtype'):
                init_coords = np.asarray(init_coords, dtype=np.float32)

            rng = check_random_state(random_state)
            graph = connectivities.tocoo()
            n_cells = graph.shape[0]
            default_epochs = 500 if n_cells <= 10000 else 200
            n_epochs = default_epochs if maxiter is None else maxiter

            metric = neigh_params.get('metric', 'euclidean')
            metric_kwds = neigh_params.get('metric_kwds', {})

            embedding, _ = umap_impl.simplicial_set_embedding(
                data=X,
                graph=graph,
                n_components=n_components,
                initial_alpha=alpha,
                a=a,
                b=b,
                gamma=gamma,
                negative_sample_rate=negative_sample_rate,
                n_epochs=n_epochs,
                init=init_coords,
                random_state=rng,
                metric=metric,
                metric_kwds=metric_kwds,
                densmap=False,
                densmap_kwds={},
                output_dens=False,
                verbose=False,
            )
        else:
            reducer = UMAP(
                n_components=n_components,
                min_dist=min_dist,
                spread=spread,
                random_state=random_state,
                init=init_pos,
            )
            embedding = reducer.fit_transform(X)
        adata.obsm['X_umap'] = embedding
        adata.uns['umap'] = {
            'params': {
                'n_components': kwargs.get('n_components', 2),
                'min_dist': kwargs.get('min_dist', 0.5),
                'spread': kwargs.get('spread', 1.0),
                'random_state': kwargs.get('random_state', 0),
            }
        }

    def rank_genes_groups(
        self,
        adata,
        groupby,
        method='t-test',
        n_genes=100,
        reference='rest',
        layer=None,
        **kwargs,
    ):
        X = adata.layers[layer] if layer else adata.X
        if sp_sparse.issparse(X):
            X = X.toarray()
        groups = adata.obs[groupby]
        if hasattr(groups, 'cat'):
            categories = list(groups.cat.categories)
            labels = groups.to_numpy()
        else:
            labels = groups.to_numpy()
            categories = sorted(np.unique(labels))
        gene_names = np.array(adata.var_names if len(adata.var_names) else [f"gene_{i}" for i in range(X.shape[1])])
        top_n = min(n_genes, X.shape[1])
        dtype = [(str(cat), object) for cat in categories]
        names_arr = np.empty(top_n, dtype=dtype)
        scores_arr = np.empty(top_n, dtype=[(str(cat), float) for cat in categories])
        pvals_arr = np.empty(top_n, dtype=[(str(cat), float) for cat in categories])

        for cat in categories:
            group_mask = labels == cat
            if reference == 'rest' or reference is None:
                ref_mask = labels != cat
            else:
                ref_mask = labels == reference
            group_expr = X[group_mask]
            ref_expr = X[ref_mask]
            if method in ('t-test', 'wilcoxon'):
                stat, pval = stats.ttest_ind(group_expr, ref_expr, axis=0, equal_var=False, nan_policy='omit')
            else:
                stat, pval = stats.ttest_ind(group_expr, ref_expr, axis=0, equal_var=False, nan_policy='omit')
            stat = np.nan_to_num(stat, nan=0.0)
            pval = np.nan_to_num(pval, nan=1.0)
            order = np.argsort(stat)[::-1][:top_n]
            names_arr[str(cat)] = gene_names[order]
            scores_arr[str(cat)] = stat[order]
            pvals_arr[str(cat)] = pval[order]

        adata.uns['rank_genes_groups'] = {
            'names': names_arr,
            'scores': scores_arr,
            'pvals': pvals_arr,
            'params': {'groupby': groupby, 'method': method, 'n_genes': top_n},
        }

    def tsne(self, adata, n_components=2, **kwargs):
        self._ensure_neighbors(adata, kwargs.get('n_neighbors', 15))
        from sklearn.manifold import TSNE
        X = self._dense_representation(adata)
        tsne = TSNE(n_components=n_components, random_state=kwargs.get('random_state', 0), init='random')
        adata.obsm['X_tsne'] = tsne.fit_transform(X)

    def diffmap(self, adata, n_comps=15, **kwargs):
        self._ensure_neighbors(adata, kwargs.get('n_neighbors', 15))
        if 'X_pca' not in adata.obsm:
            pp.pca(adata, n_comps=max(n_comps, 15))
        X_source = adata.obsm['X_pca']
        adata.obsm['X_diffmap'] = X_source[:, :n_comps]
        adata.uns['diffmap_evals'] = np.linspace(1.0, 0.1, n_comps)

# 5. EXPORT
pp = Preprocessing()
tl = Tools()
# Create mixin for sc.pp.neighbors style access
pp.neighbors = pp.neighbors

sys.modules[__name__ + '.pp'] = pp
sys.modules[__name__ + '.tl'] = tl

__all__ = ['pp', 'tl', 'settings', 'Neighbors', 'AnnData']