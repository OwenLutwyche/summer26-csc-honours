"""
Standalone preprocessing tests without pytest dependencies.
These tests mirror the original scanpy tests but can run independently.
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
from anndata import AnnData
from scipy import sparse

import scanpy as sc


def test_log1p():
    """Test log1p transformation."""
    with tempfile.TemporaryDirectory() as tmp_path:
        tmp_path = Path(tmp_path)
        
        a = np.random.rand(200, 10).astype(np.float32)
        a_log = np.log1p(a)
        ad = AnnData(a.copy())
        ad2 = AnnData(a.copy())
        ad3 = AnnData(a.copy())
        ad3.filename = tmp_path / "test.h5ad"
        
        sc.pp.log1p(ad)
        assert np.allclose(ad.X, a_log), "log1p basic transform failed"
        
        sc.pp.log1p(ad2, chunked=True)
        assert np.allclose(ad2.X, ad.X), "log1p chunked transform failed"
        
        sc.pp.log1p(ad3, chunked=True)
        assert np.allclose(ad3.X, ad.X), "log1p backed transform failed"

        # Test base
        ad4 = AnnData(a)
        sc.pp.log1p(ad4, base=2)
        assert np.allclose(ad4.X, a_log / np.log(2)), "log1p with base=2 failed"


def test_normalize_total():
    """Test normalize_total (replacement for normalize_per_cell)."""
    x = np.array([[1, 0], [3, 0], [5, 6]], dtype=np.float32)
    adata = AnnData(x.copy())
    
    sc.pp.normalize_total(adata, target_sum=1)
    sums = adata.X.sum(axis=1).tolist()
    assert np.allclose(sums, [1.0, 1.0, 1.0]), f"normalize_total failed: {sums}"
    
    # Test with sparse matrix
    adata_sparse = AnnData(sparse.csr_matrix(x.copy()))
    sc.pp.normalize_total(adata_sparse, target_sum=1)
    sums_sparse = np.asarray(adata_sparse.X.sum(axis=1)).flatten().tolist()
    assert np.allclose(sums_sparse, [1.0, 1.0, 1.0]), f"normalize_total sparse failed: {sums_sparse}"


def test_filter_cells():
    """Test filter_cells function."""
    x = np.array([
        [1, 0, 0, 0],  # 1 gene expressed
        [1, 1, 1, 0],  # 3 genes expressed
        [1, 1, 1, 1],  # 4 genes expressed
        [0, 0, 0, 0],  # 0 genes expressed
    ], dtype=np.float32)
    
    adata = AnnData(x.copy())
    sc.pp.filter_cells(adata, min_genes=2)
    
    assert adata.n_obs == 2, f"Expected 2 cells, got {adata.n_obs}"
    assert np.allclose(adata.X[0], [1, 1, 1, 0]), "Wrong cell filtered"
    assert np.allclose(adata.X[1], [1, 1, 1, 1]), "Wrong cell filtered"


def test_filter_genes():
    """Test filter_genes function."""
    x = np.array([
        [1, 0, 1, 1],
        [1, 0, 1, 0],
        [1, 0, 1, 0],
    ], dtype=np.float32)
    
    adata = AnnData(x.copy())
    sc.pp.filter_genes(adata, min_cells=2)
    
    assert adata.n_vars == 2, f"Expected 2 genes, got {adata.n_vars}"


def test_scale():
    """Test scale function."""
    np.random.seed(42)
    x = np.random.rand(100, 50).astype(np.float32)
    adata = AnnData(x.copy())
    
    sc.pp.scale(adata, zero_center=True)
    
    # Check that mean is ~0 and std is ~1
    means = adata.X.mean(axis=0)
    stds = adata.X.std(axis=0)
    
    assert np.allclose(means, 0, atol=1e-5), f"Means not zero-centered: {means[:5]}"
    assert np.allclose(stds, 1, atol=0.1), f"Stds not ~1: {stds[:5]}"


def test_scale_no_zero_center():
    """Test scale function without zero centering."""
    np.random.seed(42)
    x = np.random.rand(100, 50).astype(np.float32)
    adata = AnnData(x.copy())
    
    sc.pp.scale(adata, zero_center=False)
    
    # Check that std is ~1 (mean may not be 0)
    stds = adata.X.std(axis=0)
    assert np.allclose(stds, 1, atol=0.1), f"Stds not ~1: {stds[:5]}"


def test_highly_variable_genes():
    """Test highly_variable_genes function."""
    np.random.seed(42)
    # Create data with some highly variable genes
    n_cells, n_genes = 200, 100
    x = np.random.poisson(5, (n_cells, n_genes)).astype(np.float32)
    
    # Make some genes more variable
    x[:, :10] = np.random.poisson(20, (n_cells, 10))
    
    adata = AnnData(x)
    sc.pp.highly_variable_genes(adata, n_top_genes=20, flavor='seurat')
    
    assert 'highly_variable' in adata.var.columns, "highly_variable column not created"
    n_hvg = adata.var['highly_variable'].sum()
    assert n_hvg == 20, f"Expected 20 highly variable genes, got {n_hvg}"


def test_pca():
    """Test PCA function."""
    np.random.seed(42)
    x = np.random.rand(100, 50).astype(np.float32)
    adata = AnnData(x.copy())
    
    sc.pp.pca(adata, n_comps=10)
    print(adata.obsm['X_pca'])
    assert 'X_pca' in adata.obsm, "X_pca not created"
    assert adata.obsm['X_pca'].shape == (100, 10), f"Wrong PCA shape: {adata.obsm['X_pca'].shape}"
    assert 'pca' in adata.uns, "pca params not stored"
    assert 'PCs' in adata.varm, "PCs not stored in varm"

def test_pca_correctness():
    """Verify PCA output quality without assuming sign convention or random seed.
 
    PCA components are unique only up to sign (and rotation within degenerate
    subspaces), so direct value comparison between implementations is wrong.
    Instead we check three properties that must hold for any correct PCA:
 
    1. Explained variance ratio sums to <= 1 and is descending.
    2. The projection X_pca has the correct shape and is zero-mean per component.
    3. The subspace spanned by our components matches sklearn's — measured by
       the principal angles between the two k-dimensional subspaces. If both
       are correct, all principal angles should be ~0.
    """
    from sklearn.decomposition import PCA as SklearnPCA
 
    np.random.seed(42)
    n_obs, n_vars, n_comps = 100, 50, 10
    x = np.random.rand(n_obs, n_vars).astype(np.float32)
 
    # --- scancodon result ---
    adata = AnnData(x.copy())
    sc.pp.pca(adata, n_comps=n_comps)
    X_pca_ours = adata.obsm["X_pca"]
    vr_ours    = adata.uns["pca"]["variance_ratio"]
 
    # --- sklearn ground truth ---
    pca_sk   = SklearnPCA(n_components=n_comps)
    X_pca_sk = pca_sk.fit_transform(x.astype(np.float64))
    vr_sk    = pca_sk.explained_variance_ratio_
 
    # 1. Variance ratio sanity
    assert float(np.sum(vr_ours)) <= 1.0 + 1e-6, \
        f"variance_ratio sums to {np.sum(vr_ours):.4f} > 1"
    assert all(vr_ours[i] >= vr_ours[i+1] - 1e-6 for i in range(len(vr_ours)-1)), \
        "variance_ratio not descending"
    assert np.allclose(vr_ours, vr_sk, atol=0.01), \
        f"variance_ratio differs from sklearn:\n  ours={vr_ours}\n  sk  ={vr_sk}"
 
    # 2. X_pca shape and approximate zero mean
    assert X_pca_ours.shape == (n_obs, n_comps), f"Wrong shape: {X_pca_ours.shape}"
    col_means = np.abs(X_pca_ours.mean(axis=0))
    assert np.all(col_means < 0.05), f"X_pca columns not zero-mean: {col_means}"
 
    # 3. Subspace agreement via principal angles
    # QR-orthogonalise both projection matrices to get orthonormal bases
    Q_ours, _ = np.linalg.qr(X_pca_ours)
    Q_sk,   _ = np.linalg.qr(X_pca_sk)
    # Singular values of Q_ours.T @ Q_sk are cosines of principal angles.
    # All should be ~1 (angles ~0) if the subspaces match.
    cosines   = np.linalg.svd(Q_ours.T @ Q_sk, compute_uv=False)
    cosines   = np.clip(cosines, -1.0, 1.0)
    angles_deg = np.degrees(np.arccos(cosines))
    assert np.all(angles_deg < 5.0), \
        f"Subspace mismatch — principal angles (deg): {angles_deg}"
    print(f"  max principal angle vs sklearn: {angles_deg.max():.4f} deg  [OK]")
 

def test_pca_sparse():
    """Test PCA with sparse input."""
    np.random.seed(42)
    x = sparse.random(100, 50, density=0.3, format='csr', dtype=np.float32)
    adata = AnnData(x)
    
    sc.pp.pca(adata, n_comps=10)
    
    assert 'X_pca' in adata.obsm, "X_pca not created for sparse"
    assert adata.obsm['X_pca'].shape == (100, 10), f"Wrong PCA shape: {adata.obsm['X_pca'].shape}"


# Registry of all tests
TESTS = [
    ("test_log1p", test_log1p),
    ("test_normalize_total", test_normalize_total),
    ("test_filter_cells", test_filter_cells),
    ("test_filter_genes", test_filter_genes),
    ("test_scale", test_scale),
    ("test_scale_no_zero_center", test_scale_no_zero_center),
    ("test_highly_variable_genes", test_highly_variable_genes),
    ("test_pca", test_pca),
    ("test_pca_correctness", test_pca_correctness),
    ("test_pca_sparse", test_pca_sparse),
]


def run_all():
    """Run all tests and return results."""
    results = []
    for name, func in TESTS:
        try:
            func()
            results.append((name, True, "PASSED"))
        except AssertionError as e:
            results.append((name, False, f"FAILED: {e}"))
        except Exception as e:
            results.append((name, False, f"ERROR: {type(e).__name__}: {e}"))
    return results


if __name__ == "__main__":
    results = run_all()
    for name, success, msg in results:
        status = "✅" if success else "❌"
        print(f"{status} {name}: {msg}")
