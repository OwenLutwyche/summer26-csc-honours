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
import time

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
            X_native = np.ascontiguousarray(X, dtype=np.float64)
            X_new = scancodon_native.log1p(X_native, base)
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
        tgt = 1e4 if target_sum is None else float(target_sum) # should be median of nonzero rows, overwritten once we get counts

        is_sparse = sp_sparse.issparse(X)
        use_native = CODON_AVAILABLE and isinstance(X, np.ndarray)
        print("normalize_total")
        if use_native:
            print("using native")
            X_native = np.ascontiguousarray(X, dtype=np.float64)
            result, _ = scancodon_native.normalize_total(X_native, tgt)
        elif is_sparse:
            print("sparse case")
            if CODON_AVAILABLE:
                print("using native sparse total normalization")
                if not isinstance(X, sp_sparse.csr_matrix):
                    X = X.tocsr()
                elif not inplace:
                    X = X.copy()
                
                # 1. Force 64-bit strict types for the Codon bridge.
                # Reassigning ensures the in-place Codon mutation modifies the actual object.
                X.data = np.asarray(X.data, dtype=np.float64)
                X.indices = np.asarray(X.indices, dtype=np.int64)
                X.indptr = np.asarray(X.indptr, dtype=np.int64)

                if target_sum is None:
                    counts = np.asarray(X.sum(axis=1)).flatten()
                    tgt = float(np.median(counts[counts > 0]))
                else:
                    tgt = float(target_sum)
                
                # 2. Unpack the shape into two separate integers
                scancodon_native.normalize_total_sparse(
                    X.data, X.indices, X.indptr, X.shape[0], X.shape[1], tgt
                )
                result = X
            else:
                print("sparse, using sp_sparse.diags fallback")
                counts = np.asarray(X.sum(axis=1)).flatten()
                tgt = float(np.median(counts[counts > 0])) if target_sum is None else float(target_sum)
                scales = tgt / np.maximum(counts, 1e-12)
                result = sp_sparse.diags(scales).dot(X)
        else:
            ("else case, solution is all python.numpy")
            arr = np.asarray(X)
            counts = arr.sum(axis=1)
            tgt = float(np.median(counts[counts > 0])) if target_sum is None else float(target_sum)
            scales = tgt / np.maximum(counts, 1e-12)
            arr = arr * scales[:, None]
            if arr is not X and hasattr(X, "__setitem__"):
                X[...] = arr
                result = X
            else:
                result = arr

        if isinstance(data, AnnData):
            print("anndata case")
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
        
        # Densify if necessary before the Codon check
        if sp_sparse.issparse(X):
            X = X.toarray()
            if isinstance(adata, AnnData):
                adata.X = X
                
        use_native = CODON_AVAILABLE and isinstance(X, np.ndarray)
        if use_native:
            X_native = np.ascontiguousarray(X, dtype=np.float64)
            X_new, _, _ = scancodon_native.scale(X_native, zero_center, max_value)
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

        if sp_sparse.issparse(X):
            if CODON_AVAILABLE:
                print("filter cells sparse codon")
                mc = -1.0 if min_counts is None else float(min_counts)
                mg = -1.0 if min_genes is None else float(min_genes)
                xc = -1.0 if max_counts is None else float(max_counts)
                xg = -1.0 if max_genes is None else float(max_genes)
                
                # We only need data and indptr for row filtering!
                data_64 = np.asarray(X.data, dtype=np.float64)
                indptr_64 = np.asarray(X.indptr, dtype=np.int64)
                
                mask = scancodon_native.filter_cells_sparse(
                    data_64, indptr_64, X.shape[0], mc, mg, xc, xg
                )
        else:
            # Sparse-native path: scipy sparse reductions never materialise a dense matrix.
            # .astype(bool).sum() counts nnz per row directly from indptr — no toarray() needed.
            if min_genes is not None or max_genes is not None:
                number_per_cell = np.asarray(X.astype(bool).sum(axis=1)).flatten().astype(np.float64)
            else:
                number_per_cell = np.asarray(X.sum(axis=1)).flatten().astype(np.float64)
            if min_counts is not None:
                mask = number_per_cell >= float(min_counts)
            elif min_genes is not None:
                mask = number_per_cell >= float(min_genes)
            elif max_counts is not None:
                mask = number_per_cell <= float(max_counts)
            else:  # max_genes
                mask = number_per_cell <= float(max_genes)
            # Already dense — Codon kernel handles the reduction natively
            X_native = np.ascontiguousarray(X, dtype=np.float64)
            mask, _ = scancodon_native.filter_cells(X_native, min_counts, min_genes, max_counts, max_genes)
            #mask = self._filter_cells_numpy(X, min_counts, min_genes, max_counts, max_genes)

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

    def _ndarray_to_csr(self, data, row, col, n_obs):
        """Convert ndarrays returned by native Codon connectivity kernels
        (e.g. gauss_connectivity) into a scipy CSR matrix.

        Codon kernels return flat ndarray triples rather than scipy objects
        to avoid bridging overhead inside the compiled library. This helper
        sits at the Python receiver layer and does the final construction.

        Usage:
            data, row, col = scancodon_native.gauss_connectivity(
                indices, distances, n_obs)
            connectivities = self._ndarray_to_csr(data, row, col, n_obs)
        """
        return sp_sparse.csr_matrix(
            (np.array(data), (np.array(row), np.array(col))),
            shape=(n_obs, n_obs),
        )

    def filter_genes(self, data, min_cells=None, min_counts=None, max_cells=None, max_counts=None, inplace=True, **kwargs):
        adata = data if inplace else data.copy()
        X = self._get_x(adata)

        if sp_sparse.issparse(X):
            if CODON_AVAILABLE:
                print("filter genes sparse codon")
                mc = -1.0 if min_counts is None else float(min_counts)
                mcell = -1.0 if min_cells is None else float(min_cells)
                xc = -1.0 if max_counts is None else float(max_counts)
                xcell = -1.0 if max_cells is None else float(max_cells)
                
                # We only need data and indices for column filtering!
                data_64 = np.asarray(X.data, dtype=np.float64)
                indices_64 = np.asarray(X.indices, dtype=np.int64)
                
                mask = scancodon_native.filter_genes_sparse(
                    data_64, indices_64, X.shape[1], mc, mcell, xc, xcell
                )
        else:
            # Sparse-native path: column reductions over CSC/CSR without materialising dense matrix.
            # .astype(bool).sum() counts nnz per column directly — no toarray() needed.
            if min_cells is not None or max_cells is not None:
                number_per_gene = np.asarray(X.astype(bool).sum(axis=0)).flatten().astype(np.float64)
            else:
                number_per_gene = np.asarray(X.sum(axis=0)).flatten().astype(np.float64)
            if min_counts is not None:
                mask = number_per_gene >= float(min_counts)
            elif min_cells is not None:
                mask = number_per_gene >= float(min_cells)
            elif max_counts is not None:
                mask = number_per_gene <= float(max_counts)
            else:  # max_cells
                mask = number_per_gene <= float(max_cells)
        # elif CODON_AVAILABLE:
        #     # Already dense — Codon kernel handles the reduction natively
        #     X_native = np.ascontiguousarray(X, dtype=np.float64)
        #     mask, _ = scancodon_native.filter_genes(X_native, min_counts, min_cells, max_counts, max_cells)
        # else:
        #     mask = self._filter_genes_numpy(X, min_cells, min_counts, max_cells, max_counts)

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
            X_native = np.ascontiguousarray(X, dtype=np.float64)
            mask, means, vars_, dispersions, dispersions_norm = scancodon_native.highly_variable_genes_seurat_dense(X_native, n_top_genes)
            adata.var['highly_variable'] = np.array(mask, dtype=bool)
            adata.var['means'] = np.array(means)
            adata.var['dispersions'] = np.array(dispersions)
            adata.var['dispersions_norm'] = np.array(dispersions_norm)
        else:
            sc.preprocessing.highly_variable_genes(adata, n_top_genes=20, flavor='seurat')
        if subset: adata._inplace_subset_var(adata.var['highly_variable'])

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

        if CODON_AVAILABLE:
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
            
            # write back to anndata
            adata.obs['doublet_score'] = scores
            # TODO add predicted_doublets_ and uns
            
        else:
            print("[INFO] Falling back to scanpy.pp.scrublet")
            sc.pp.scrublet(
                adata, 
                sim_doublet_ratio=sim_doublet_ratio,
                expected_doublet_rate=expected_doublet_rate,
                stdev_doublet_rate=stdev_doublet_rate,
                random_state=random_state,
                **kwargs
            )

        return None if inplace else adata
    




    def pca(self, data, n_comps=50, zero_center=True, **kwargs):
        t0 = time.time()
        adata = data
        layer           = kwargs.get('layer', None)
        key_added       = kwargs.get('key_added', None)
        # scanpy's real default for use_highly_variable is None ("auto-detect":
        # use HVGs if the column exists, ignore otherwise) — NOT False.
        # Defaulting to False here meant PCA ran on the full gene set even
        # when highly_variable_genes() had already been run, causing a large
        # variance inflation relative to scanpy (verified: ~2.86x on 3k PBMC,
        # traced to 13714 genes used here vs ~2000 HVGs used by scanpy).
        use_highly_variable_kw = kwargs.get('use_highly_variable', None)
        if use_highly_variable_kw is None:
            use_highly_variable = 'highly_variable' in adata.var.columns
        else:
            use_highly_variable = bool(use_highly_variable_kw)
        standardize     = bool(kwargs.get('standardize', False))
        dtype           = kwargs.get('dtype', 'float32')

        # If use_highly_variable and the column exists, use highly variable genes only.
        # This follows the same pattern as scanpy's own PCA implementation.
        if use_highly_variable and 'highly_variable' in adata.var.columns:
            mask_var = adata.var['highly_variable'].values
            adata_comp = adata[:, mask_var]
        else:
            mask_var = None
            adata_comp = adata

        #X_raw = adata_comp.layers[layer] if layer is not None else adata_comp.X
        #X_raw = np.ascontiguousarray(adata_comp.X, dtype=np.float64)
        # Extract the matrix
        X_matrix = adata_comp.X
        
        # Explicitly densify if it's a sparse matrix
        if sp_sparse.issparse(X_matrix):
            X_matrix = X_matrix.toarray()
            
        # NOW it is safe to cast and align memory
        X_raw = np.ascontiguousarray(X_matrix, dtype=np.float64)
        
        if CODON_AVAILABLE:
            # Pass plain ndarray[float,2]
            X_pca, components, variance_ratio, variance = scancodon_native.pca(
                X_raw,
                n_comps=int(n_comps),
                zero_center=bool(zero_center),
                standardize=standardize,
            )
        else:
            from sklearn.decomposition import PCA
            pca_obj = PCA(n_components=n_comps)
            # densify before passing to sklearn PCA
            X_dense = X_raw.toarray() if sp_sparse.issparse(X_raw) else X_raw
            X_native = np.ascontiguousarray(X_dense, dtype=np.float64)
            X_pca = pca_obj.fit_transform(X_native)
            components   = pca_obj.components_
            variance_ratio = pca_obj.explained_variance_ratio_
            variance     = pca_obj.explained_variance_

        # write results back into the AnnData
        key_obsm = key_added if key_added else 'X_pca'
        key_varm = key_added if key_added else 'PCs'
        key_uns  = key_added if key_added else 'pca'

        adata.obsm[key_obsm] = np.array(X_pca, dtype=dtype)

        # If we subsetted to HVGs, zero-pad varm back to full gene width
        if mask_var is not None:
            n_genes = adata.n_vars
            PCs_full = np.zeros((n_genes, n_comps), dtype=dtype)
            PCs_full[mask_var] = np.array(components, dtype=dtype).T
            adata.varm[key_varm] = PCs_full
        else:
            adata.varm[key_varm] = np.array(components, dtype=dtype).T

        adata.uns[key_uns] = {
            'params': {
                'zero_center': zero_center,
                'use_highly_variable': use_highly_variable,
            },
            'variance_ratio': np.array(variance_ratio),
            'variance':       np.array(variance),
        }

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
    
import scanpy as sc
import numpy as np
import scipy.sparse as sp_sparse
from anndata import AnnData


def calculate_qc_metrics(
    adata: AnnData, 
    expr_type: str = "counts", 
    var_type: str = "genes", 
    qc_vars: dict = None, 
    percent_top: list[int] = (50, 100, 200, 500), 
    layer: str = None, 
    use_raw: bool = False, 
    inplace: bool = True, 
    log1p: bool = True,
    parallel: bool = True
):
    # 1. Prep data (Keep it sparse!)
    X = adata.raw.X if use_raw else adata.X
    if layer is not None:
        X = adata.layers[layer]
    
    # Ensure it's CSR for Codon
    X_csr = X.tocsr() if not sp_sparse.isspmatrix_csr(X) else X

    # 2. Check for native execution
    use_native = CODON_AVAILABLE and sp_sparse.isspmatrix_csr(X_csr)

    if use_native:
        # Prepare inputs
        qc_names = list(qc_vars.keys()) if qc_vars else []
        qc_masks = [qc_vars[name].values for name in qc_names]
        
        # Call Codon Bridge (passing data buffers, NOT the matrix object)
        obs_tuple, var_tuple = scancodon_native.calculate_qc_metrics(
            X_csr.data, X_csr.indices, X_csr.indptr, 
            X_csr.shape, qc_masks, percent_top, log1p
        )
        
        # Unpack obs results
        nnz, log1p_nnz, totals, log1p_totals, props, qc_totals, log1p_qc_totals, qc_pcts = obs_tuple
        
        # Stitch back to obs
        adata.obs[f"n_genes_by_{expr_type}"] = nnz
        adata.obs[f"log1p_n_genes_by_{expr_type}"] = log1p_nnz
        adata.obs[f"total_{expr_type}"] = totals
        adata.obs[f"log1p_total_{expr_type}"] = log1p_totals
        
        for i, n in enumerate(percent_top):
            adata.obs[f"pct_counts_in_top_{n}_{var_type}"] = props[:, i]
            
        for i, name in enumerate(qc_names):
            adata.obs[f"total_{expr_type}_{name}"] = qc_totals[:, i]
            adata.obs[f"log1p_total_{expr_type}_{name}"] = log1p_qc_totals[:, i]
            adata.obs[f"pct_counts_{name}"] = qc_pcts[:, i]

        # Unpack var results (nnz, means, log1p_means, dropout, totals, log1p_totals)
        v_nnz, v_means, v_log_means, v_dropout, v_totals, v_log_totals = var_tuple
        adata.var[f"n_cells_by_{expr_type}"] = v_nnz
        adata.var[f"mean_{expr_type}"] = v_means
        adata.var[f"log1p_mean_{expr_type}"] = v_log_means
        adata.var[f"pct_dropout_by_{expr_type}"] = v_dropout
        adata.var[f"total_{expr_type}"] = v_totals
        adata.var[f"log1p_total_{expr_type}"] = v_log_totals
        
        return None if inplace else adata

    else:
        # 3. Fallback
        return sc.pp.calculate_qc_metrics(
            adata, expr_type=expr_type, var_type=var_type, 
            qc_vars=qc_vars, percent_top=percent_top, 
            layer=layer, use_raw=use_raw, inplace=inplace, log1p=log1p
        )

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

        # if codon supports this function call, we pass it to the rank_genes_groups_dispatcher
        use_native = CODON_AVAILABLE and method in ('t-test', 't-test_overestim_var', 'wilcoxon')

        if use_native:
            X = np.ascontiguousarray(X, dtype=np.float64)
            
            # if codon is available, 
            # build groups_masks array (shape n_groups, n_cells), dtype = bool
            # resolve ireference
            # call scancodon_native.rank_genes_groups_dispatcher(X, groups_masks, method, ireference)
            # finally, unpack returned list 

            # groups_masks: 2d booleran array of shape (n_groups, n_cells) each row i is true if belongs to categories[i]
            groups_masks = np.array([labels == cat for cat in categories], dtype = bool)
            # index of reference group into categories
            ireference = None if reference == 'rest' else categories.index(reference)
            # track how long the dispatcher call takes:
            #start_time = time.time()
            if ireference is None:
                results = scancodon_native.rank_genes_groups_dispatcher(X, groups_masks, method)
            else:
                results = scancodon_native.rank_genes_groups_dispatcher(X, groups_masks, method, ireference)
            # tools is not exposed as part of the scancodon library, should be dispatched from top level instead
            #end_time = time.time()
            #ttest_time = end_time - start_time
            #print(f"ttest time: {ttest_time:.4f}")

            # now we unpack it
            # results is a list of group_idx, scores, pvals
            for group_idx, scores, pvals in results:
                cat = categories[group_idx]
                order = np.argsort(scores)[::-1][:top_n]
                names_arr[str(cat)] = gene_names[order]
                scores_arr[str(cat)] = scores[order]
                pvals_arr[str(cat)] = pvals[order]

        else:
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
# Create aliases to match scanpy API (scanpy has it located in both pp and tl)
pp.neighbors = pp.neighbors
tl.pca = pp.pca  # pca is available in both pp and tl in scanpy

sys.modules[__name__ + '.pp'] = pp
sys.modules[__name__ + '.tl'] = tl

__all__ = ['pp', 'tl', 'settings', 'Neighbors', 'AnnData']