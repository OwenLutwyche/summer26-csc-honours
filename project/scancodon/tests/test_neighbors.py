"""
Standalone neighbors tests without pytest dependencies.
These tests mirror the original scanpy tests but can run independently.
"""
from __future__ import annotations

import numpy as np
from anndata import AnnData
from scipy import sparse

import scanpy as sc
from scanpy import Neighbors


# the input data
X = [[1, 0], [3, 0], [5, 6], [0, 4]]
n_neighbors = 3  # includes data points themselves

# Expected distances (Euclidean)
distances_euclidean = [
    [0.0, 2.0, 0.0, 4.123105525970459],
    [2.0, 0.0, 0.0, 5.0],
    [0.0, 6.324555397033691, 0.0, 5.385164737701416],
    [4.123105525970459, 5.0, 0.0, 0.0],
]


def get_neighbors() -> Neighbors:
    return Neighbors(AnnData(np.array(X)))


def test_neighbors_basic():
    """Test basic neighbors computation."""
    adata = AnnData(np.array(X, dtype=np.float32))
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X')
    
    assert 'neighbors' in adata.uns, "neighbors not stored in uns"
    assert 'distances' in adata.obsp, "distances not stored in obsp"
    assert 'connectivities' in adata.obsp, "connectivities not stored in obsp"


def test_neighbors_random_data():
    """Test neighbors with random data."""
    np.random.seed(42)
    adata = AnnData(np.random.rand(100, 50).astype(np.float32))
    
    sc.pp.pca(adata, n_comps=20)
    sc.pp.neighbors(adata, n_neighbors=15)
    
    assert 'neighbors' in adata.uns, "neighbors not stored"
    assert adata.obsp['distances'].shape == (100, 100), "Wrong distances shape"
    assert adata.obsp['connectivities'].shape == (100, 100), "Wrong connectivities shape"


def test_neighbors_methods():
    """Test different neighbor computation methods."""
    np.random.seed(42)
    adata = AnnData(np.random.rand(100, 50).astype(np.float32))
    sc.pp.pca(adata, n_comps=20)
    
    # Test with umap method (default)
    sc.pp.neighbors(adata, n_neighbors=15, method='umap')
    assert 'neighbors' in adata.uns, "umap method failed"
    
    # Test with gauss method
    adata2 = adata.copy()
    sc.pp.neighbors(adata2, n_neighbors=15, method='gauss')
    assert 'neighbors' in adata2.uns, "gauss method failed"


def test_neighbors_sparse():
    """Test neighbors with sparse input after PCA."""
    np.random.seed(42)
    x = sparse.random(100, 50, density=0.3, format='csr', dtype=np.float32)
    adata = AnnData(x)
    
    sc.pp.pca(adata, n_comps=20)
    sc.pp.neighbors(adata, n_neighbors=15)
    
    assert 'neighbors' in adata.uns, "neighbors not stored for sparse"


def test_neighbors_n_pcs():
    """Test neighbors with different n_pcs values."""
    np.random.seed(42)
    adata = AnnData(np.random.rand(100, 50).astype(np.float32))
    sc.pp.pca(adata, n_comps=30)
    
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=10)
    assert adata.uns['neighbors']['params']['n_pcs'] == 10, "n_pcs not stored correctly"


# Registry of all tests
TESTS = [
    ("test_neighbors_basic", test_neighbors_basic),
    ("test_neighbors_random_data", test_neighbors_random_data),
    ("test_neighbors_methods", test_neighbors_methods),
    ("test_neighbors_sparse", test_neighbors_sparse),
    ("test_neighbors_n_pcs", test_neighbors_n_pcs),
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
