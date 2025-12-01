"""
Standalone clustering tests without pytest dependencies.
These tests mirror the original scanpy tests but can run independently.
"""
from __future__ import annotations

import warnings
# Suppress FutureWarning about leiden default backend change
warnings.filterwarnings('ignore', category=FutureWarning, message='.*default backend for leiden.*')
warnings.filterwarnings('ignore', category=FutureWarning, message='.*Use `scanpy.tl.leiden`.*')

import numpy as np
from anndata import AnnData

import scanpy as sc


def get_adata_with_neighbors():
    """Create test AnnData with neighbors computed."""
    np.random.seed(42)
    adata = AnnData(np.random.rand(200, 50).astype(np.float32))
    sc.pp.pca(adata, n_comps=20)
    sc.pp.neighbors(adata, n_neighbors=15)
    return adata


def test_leiden_basic():
    """Test basic Leiden clustering."""
    adata = get_adata_with_neighbors()
    
    sc.tl.leiden(adata)
    
    assert 'leiden' in adata.obs.columns, "leiden column not created"
    assert adata.obs['leiden'].nunique() > 1, "Only one cluster found"
    assert 'leiden' in adata.uns, "leiden params not stored"


def test_leiden_resolution():
    """Test Leiden with different resolutions."""
    adata = get_adata_with_neighbors()
    
    # Lower resolution = fewer clusters
    sc.tl.leiden(adata, resolution=0.5, key_added='leiden_low')
    n_clusters_low = adata.obs['leiden_low'].nunique()
    
    # Higher resolution = more clusters
    sc.tl.leiden(adata, resolution=2.0, key_added='leiden_high')
    n_clusters_high = adata.obs['leiden_high'].nunique()
    
    assert n_clusters_high >= n_clusters_low, \
        f"Higher resolution should give more clusters: {n_clusters_high} vs {n_clusters_low}"


def test_leiden_random_state():
    """Test Leiden reproducibility with random_state."""
    adata = get_adata_with_neighbors()
    
    # Same random state should give same results
    adata1 = adata.copy()
    sc.tl.leiden(adata1, random_state=42)
    
    adata2 = adata.copy()
    sc.tl.leiden(adata2, random_state=42)
    
    assert (adata1.obs['leiden'] == adata2.obs['leiden']).all(), \
        "Same random_state should give same clusters"


def test_leiden_key_added():
    """Test Leiden with custom key_added."""
    adata = get_adata_with_neighbors()
    
    sc.tl.leiden(adata, key_added='my_clusters')
    
    assert 'my_clusters' in adata.obs.columns, "Custom key_added not used"
    assert 'leiden' not in adata.obs.columns, "Default key shouldn't exist"


def test_louvain_basic():
    """Test basic Louvain clustering (if available)."""
    adata = get_adata_with_neighbors()
    
    try:
        sc.tl.louvain(adata)
        assert 'louvain' in adata.obs.columns, "louvain column not created"
        assert adata.obs['louvain'].nunique() > 1, "Only one cluster found"
    except ImportError:
        # louvain package not installed, skip
        pass


# Registry of all tests
TESTS = [
    ("test_leiden_basic", test_leiden_basic),
    ("test_leiden_resolution", test_leiden_resolution),
    ("test_leiden_random_state", test_leiden_random_state),
    ("test_leiden_key_added", test_leiden_key_added),
    ("test_louvain_basic", test_louvain_basic),
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
