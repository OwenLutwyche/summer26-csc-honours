"""
Standalone rank_genes_groups tests without pytest dependencies.
These tests mirror the original scanpy tests but can run independently.
"""
from __future__ import annotations

import warnings
# Suppress the "raw count data" warning - we intentionally test on raw data
# to match the original scanpy tests. Must be set before importing scanpy.
warnings.filterwarnings('ignore', message='.*rank_genes_groups on the raw count data.*')

import numpy as np
import pandas as pd
from anndata import AnnData
from numpy.random import binomial, negative_binomial, seed
from scipy import sparse

import scanpy as sc

# Also suppress via scanpy's logging (it uses logg.warning, not warnings.warn)
sc.settings.verbosity = 0  # Suppress info/warning messages


def get_example_data(use_sparse=False):
    """Create test AnnData for rank_genes_groups."""
    seed(1234)
    
    # Create test object
    adata = AnnData(
        np.multiply(binomial(1, 0.15, (100, 20)), negative_binomial(2, 0.25, (100, 20))).astype(np.float32)
    )
    
    # Adapt marker genes for cluster (so as to have some form of reasonable input)
    adata.X[0:10, 0:5] = np.multiply(
        binomial(1, 0.9, (10, 5)), negative_binomial(1, 0.5, (10, 5))
    ).astype(np.float32)

    if use_sparse:
        adata.X = sparse.csr_matrix(adata.X)

    # Create cluster according to groups
    adata.obs["true_groups"] = pd.Categorical(
        np.concatenate((np.zeros((10,), dtype=int), np.ones((90,), dtype=int)))
    )

    return adata


def test_rank_genes_groups_ttest():
    """Test rank_genes_groups with t-test method."""
    adata = get_example_data()
    
    sc.tl.rank_genes_groups(adata, "true_groups", method="t-test", n_genes=20)
    
    assert 'rank_genes_groups' in adata.uns, "rank_genes_groups not stored"
    assert 'names' in adata.uns['rank_genes_groups'], "names not stored"
    assert 'scores' in adata.uns['rank_genes_groups'], "scores not stored"
    assert 'pvals' in adata.uns['rank_genes_groups'], "pvals not stored"
    

def test_rank_genes_groups_wilcoxon():
    """Test rank_genes_groups with Wilcoxon method."""
    adata = get_example_data()
    
    sc.tl.rank_genes_groups(adata, "true_groups", method="wilcoxon", n_genes=20)
    
    assert 'rank_genes_groups' in adata.uns, "rank_genes_groups not stored"
    assert 'names' in adata.uns['rank_genes_groups'], "names not stored"


def test_rank_genes_groups_sparse():
    """Test rank_genes_groups with sparse matrix."""
    adata = get_example_data(use_sparse=True)
    
    sc.tl.rank_genes_groups(adata, "true_groups", method="t-test", n_genes=20)
    
    assert 'rank_genes_groups' in adata.uns, "rank_genes_groups not stored for sparse"


def test_rank_genes_groups_multiple_groups():
    """Test rank_genes_groups with multiple groups."""
    seed(42)
    n_cells, n_genes = 150, 30
    adata = AnnData(np.random.rand(n_cells, n_genes).astype(np.float32))
    
    # Create 3 groups
    groups = np.array(['A'] * 50 + ['B'] * 50 + ['C'] * 50)
    adata.obs['group'] = pd.Categorical(groups)
    
    # Make group A have high expression of genes 0-9
    adata.X[:50, :10] += 5
    # Make group B have high expression of genes 10-19
    adata.X[50:100, 10:20] += 5
    
    sc.tl.rank_genes_groups(adata, 'group', method='t-test')
    
    assert 'rank_genes_groups' in adata.uns, "rank_genes_groups not stored"
    # Check that we have results for all groups
    names = adata.uns['rank_genes_groups']['names']
    assert names.dtype.names is not None, "No group names in results"


def test_rank_genes_groups_reference():
    """Test rank_genes_groups with a reference group."""
    adata = get_example_data()
    
    # Use group 1 (integer) as reference
    sc.tl.rank_genes_groups(adata, "true_groups", reference=1, method="t-test")
    
    assert 'rank_genes_groups' in adata.uns, "rank_genes_groups not stored"


def test_rank_genes_groups_layer():
    """Test rank_genes_groups using a layer."""
    adata = get_example_data()
    adata.layers['test_layer'] = adata.X.copy()
    
    # Modify X to be different
    adata.X = np.zeros_like(adata.X)
    
    sc.tl.rank_genes_groups(adata, "true_groups", method="t-test", layer='test_layer')
    
    assert 'rank_genes_groups' in adata.uns, "rank_genes_groups not stored"
    # Should have non-trivial results since we used the layer
    scores = adata.uns['rank_genes_groups']['scores']['0']
    assert not np.allclose(scores, 0), "Scores should not all be zero when using layer"


def test_rank_genes_groups_n_genes():
    """Test rank_genes_groups with different n_genes."""
    adata = get_example_data()
    
    sc.tl.rank_genes_groups(adata, "true_groups", method="t-test", n_genes=5)
    
    # Check that we get the right number of genes
    names = adata.uns['rank_genes_groups']['names']
    assert len(names) == 5, f"Expected 5 genes, got {len(names)}"


# Registry of all tests
TESTS = [
    ("test_rank_genes_groups_ttest", test_rank_genes_groups_ttest),
    ("test_rank_genes_groups_wilcoxon", test_rank_genes_groups_wilcoxon),
    ("test_rank_genes_groups_sparse", test_rank_genes_groups_sparse),
    ("test_rank_genes_groups_multiple_groups", test_rank_genes_groups_multiple_groups),
    ("test_rank_genes_groups_reference", test_rank_genes_groups_reference),
    ("test_rank_genes_groups_layer", test_rank_genes_groups_layer),
    ("test_rank_genes_groups_n_genes", test_rank_genes_groups_n_genes),
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
