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
    
    assert 'X_pca' in adata.obsm, "X_pca not created"
    assert adata.obsm['X_pca'].shape == (100, 10), f"Wrong PCA shape: {adata.obsm['X_pca'].shape}"
    assert 'pca' in adata.uns, "pca params not stored"
    assert 'PCs' in adata.varm, "PCs not stored in varm"


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
