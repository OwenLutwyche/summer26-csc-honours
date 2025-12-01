"""
Standalone embedding tests without pytest dependencies.
These tests mirror the original scanpy tests but can run independently.
"""
from __future__ import annotations

import numpy as np
from anndata import AnnData

import scanpy as sc


def get_adata_with_neighbors():
    """Create test AnnData with neighbors computed."""
    np.random.seed(42)
    adata = AnnData(np.random.rand(100, 50).astype(np.float32))
    sc.pp.pca(adata, n_comps=20)
    sc.pp.neighbors(adata, n_neighbors=15)
    return adata


def test_umap_basic():
    """Test basic UMAP embedding."""
    adata = get_adata_with_neighbors()
    
    sc.tl.umap(adata)
    
    assert 'X_umap' in adata.obsm, "X_umap not created"
    assert adata.obsm['X_umap'].shape == (100, 2), f"Wrong UMAP shape: {adata.obsm['X_umap'].shape}"
    assert 'umap' in adata.uns, "umap params not stored"


def test_umap_n_components():
    """Test UMAP with different n_components."""
    adata = get_adata_with_neighbors()
    
    sc.tl.umap(adata, n_components=3)
    
    assert adata.obsm['X_umap'].shape == (100, 3), f"Wrong UMAP shape: {adata.obsm['X_umap'].shape}"


def test_umap_min_dist():
    """Test UMAP with different min_dist values."""
    adata = get_adata_with_neighbors()
    
    sc.tl.umap(adata, min_dist=0.1)
    umap1 = adata.obsm['X_umap'].copy()
    
    sc.tl.umap(adata, min_dist=0.5)
    umap2 = adata.obsm['X_umap'].copy()
    
    # Different min_dist should give different embeddings
    assert not np.allclose(umap1, umap2), "Different min_dist should give different embeddings"


def test_umap_random_state():
    """Test UMAP reproducibility with random_state."""
    adata1 = get_adata_with_neighbors()
    adata2 = adata1.copy()
    
    sc.tl.umap(adata1, random_state=42)
    sc.tl.umap(adata2, random_state=42)
    
    assert np.allclose(adata1.obsm['X_umap'], adata2.obsm['X_umap']), \
        "Same random_state should give same UMAP"


def test_umap_init_pos():
    """Test UMAP with custom initial positions."""
    adata = get_adata_with_neighbors()
    
    # Use PCA as initial positions
    init_pos = adata.obsm['X_pca'][:, :2].astype(np.float32)
    sc.tl.umap(adata, init_pos=init_pos)
    
    assert 'X_umap' in adata.obsm, "X_umap not created with init_pos"


def test_tsne_basic():
    """Test basic t-SNE embedding."""
    adata = get_adata_with_neighbors()
    
    sc.tl.tsne(adata)
    
    assert 'X_tsne' in adata.obsm, "X_tsne not created"
    assert adata.obsm['X_tsne'].shape == (100, 2), f"Wrong t-SNE shape: {adata.obsm['X_tsne'].shape}"


def test_diffmap_basic():
    """Test basic diffusion map."""
    adata = get_adata_with_neighbors()
    
    sc.tl.diffmap(adata)
    
    assert 'X_diffmap' in adata.obsm, "X_diffmap not created"
    assert 'diffmap_evals' in adata.uns, "diffmap_evals not stored"


def test_diffmap_reproducibility():
    """Test diffusion map reproducibility."""
    adata = get_adata_with_neighbors()
    
    sc.tl.diffmap(adata)
    d1 = adata.obsm['X_diffmap'].copy()
    
    sc.tl.diffmap(adata)
    d2 = adata.obsm['X_diffmap'].copy()
    
    assert np.allclose(d1, d2), "Diffmap should be reproducible"


# Registry of all tests
TESTS = [
    ("test_umap_basic", test_umap_basic),
    ("test_umap_n_components", test_umap_n_components),
    ("test_umap_min_dist", test_umap_min_dist),
    ("test_umap_random_state", test_umap_random_state),
    ("test_umap_init_pos", test_umap_init_pos),
    ("test_tsne_basic", test_tsne_basic),
    ("test_diffmap_basic", test_diffmap_basic),
    ("test_diffmap_reproducibility", test_diffmap_reproducibility),
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
