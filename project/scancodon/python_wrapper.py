"""
Python wrapper for Codon preprocessing kernels.
This module follows the "Sandwich Pattern":
  - Python layer handles AnnData objects
  - Codon kernels do the computational work on numpy arrays
  - Python layer puts results back into AnnData
"""

import numpy as np
from typing import Optional, Union
from anndata import AnnData

# This will be the compiled Codon extension
# For now, we'll use Python scanpy as fallback
try:
    # When compiled as extension with: codon build --pyext --module scancodon_kernels
    import scancodon_kernels.preprocessing as codon_pp
    CODON_AVAILABLE = True
except ImportError:
    CODON_AVAILABLE = False
    import scanpy as fallback_scanpy


def filter_cells(
    data: Union[AnnData, np.ndarray],
    min_counts: Optional[int] = None,
    min_genes: Optional[int] = None,
    max_counts: Optional[int] = None,
    max_genes: Optional[int] = None,
    inplace: bool = True,
    copy: bool = False,
):
    """Filter cell outliers based on counts and numbers of genes expressed."""
    n_given = sum(x is not None for x in [min_counts, min_genes, max_counts, max_genes])
    if n_given != 1:
        raise ValueError(
            "Only provide one of the optional parameters `min_counts`, "
            "`min_genes`, `max_counts`, `max_genes` per call."
        )
    
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        
        # Extract matrix for Codon kernel
        x = adata.X
        if hasattr(x, "toarray"):
            x = x.toarray()
        
        # Call Codon kernel
        if CODON_AVAILABLE:
            cell_subset, number_per_cell = codon_pp.filter_cells_kernel(
                x, min_counts, min_genes, max_counts, max_genes
            )
        else:
            # Fallback to Python scanpy
            return fallback_scanpy.pp.filter_cells(
                data, min_counts=min_counts, min_genes=min_genes,
                max_counts=max_counts, max_genes=max_genes,
                inplace=inplace, copy=copy
            )
        
        if not inplace:
            return cell_subset, number_per_cell
        
        # Update AnnData
        if min_genes is None and max_genes is None:
            adata.obs["n_counts"] = number_per_cell
        else:
            adata.obs["n_genes"] = number_per_cell
        
        adata._inplace_subset_obs(cell_subset)
        return adata if copy else None
    else:
        # Direct matrix input
        if CODON_AVAILABLE:
            return codon_pp.filter_cells_kernel(
                data, min_counts, min_genes, max_counts, max_genes
            )
        else:
            return fallback_scanpy.pp.filter_cells(
                data, min_counts=min_counts, min_genes=min_genes,
                max_counts=max_counts, max_genes=max_genes,
                inplace=False
            )


def filter_genes(
    data: Union[AnnData, np.ndarray],
    min_counts: Optional[int] = None,
    min_cells: Optional[int] = None,
    max_counts: Optional[int] = None,
    max_cells: Optional[int] = None,
    inplace: bool = True,
    copy: bool = False,
):
    """Filter genes based on number of cells or counts."""
    n_given = sum(x is not None for x in [min_counts, min_cells, max_counts, max_cells])
    if n_given != 1:
        raise ValueError(
            "Only provide one of the optional parameters `min_counts`, "
            "`min_cells`, `max_counts`, `max_cells` per call."
        )
    
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        
        # Extract matrix for Codon kernel
        x = adata.X
        if hasattr(x, "toarray"):
            x = x.toarray()
        
        # Call Codon kernel
        if CODON_AVAILABLE:
            gene_subset, number_per_gene = codon_pp.filter_genes_kernel(
                x, min_counts, min_cells, max_counts, max_cells
            )
        else:
            return fallback_scanpy.pp.filter_genes(
                data, min_counts=min_counts, min_cells=min_cells,
                max_counts=max_counts, max_cells=max_cells,
                inplace=inplace, copy=copy
            )
        
        if not inplace:
            return gene_subset, number_per_gene
        
        # Update AnnData
        if min_cells is None and max_cells is None:
            adata.var["n_counts"] = number_per_gene
        else:
            adata.var["n_cells"] = number_per_gene
        
        adata._inplace_subset_var(gene_subset)
        return adata if copy else None
    else:
        # Direct matrix input
        if CODON_AVAILABLE:
            return codon_pp.filter_genes_kernel(
                data, min_counts, min_cells, max_counts, max_cells
            )
        else:
            return fallback_scanpy.pp.filter_genes(
                data, min_counts=min_counts, min_cells=min_cells,
                max_counts=max_counts, max_cells=max_cells,
                inplace=False
            )


def log1p(
    data: Union[AnnData, np.ndarray],
    base: Optional[float] = None,
    copy: bool = False,
    layer: Optional[str] = None,
    obsm: Optional[str] = None,
):
    """Logarithmize the data matrix. Computes X = log(X + 1)."""
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        
        # Extract matrix for Codon kernel
        if layer is not None:
            x = adata.layers[layer]
        elif obsm is not None:
            x = adata.obsm[obsm]
        else:
            x = adata.X
        
        if hasattr(x, "toarray"):
            x = x.toarray()
        
        # Call Codon kernel
        if CODON_AVAILABLE:
            x = codon_pp.log1p_kernel(x.copy() if not copy else x, base)
        else:
            return fallback_scanpy.pp.log1p(
                data, base=base, copy=copy, layer=layer, obsm=obsm
            )
        
        # Put back into AnnData
        if layer is not None:
            adata.layers[layer] = x
        elif obsm is not None:
            adata.obsm[obsm] = x
        else:
            adata.X = x
        
        adata.uns["log1p"] = {"base": base}
        return adata if copy else None
    else:
        # Direct matrix input
        if CODON_AVAILABLE:
            x = data.copy() if copy else data
            return codon_pp.log1p_kernel(x, base)
        else:
            return fallback_scanpy.pp.log1p(data, base=base, copy=copy)


def normalize_total(
    data: Union[AnnData, np.ndarray],
    target_sum: Optional[float] = None,
    exclude_highly_expressed: bool = False,
    max_fraction: float = 0.05,
    key_added: Optional[str] = None,
    layer: Optional[str] = None,
    inplace: bool = True,
    copy: bool = False,
):
    """Normalize counts per cell."""
    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        
        # Extract matrix for Codon kernel
        if layer is not None:
            x = adata.layers[layer]
        else:
            x = adata.X
        
        if hasattr(x, "toarray"):
            x = x.toarray()
        
        # Call Codon kernel
        if CODON_AVAILABLE:
            x, norm_factors = codon_pp.normalize_total_kernel(
                x.copy() if not copy else x,
                target_sum,
                exclude_highly_expressed,
                max_fraction
            )
        else:
            return fallback_scanpy.pp.normalize_total(
                data, target_sum=target_sum,
                exclude_highly_expressed=exclude_highly_expressed,
                max_fraction=max_fraction,
                key_added=key_added,
                layer=layer,
                inplace=inplace,
                copy=copy
            )
        
        # Put back into AnnData
        if layer is not None:
            adata.layers[layer] = x
        else:
            adata.X = x
        
        if key_added is not None:
            adata.obs[key_added] = norm_factors
        
        if not inplace:
            return {"X": x, "norm_factor": norm_factors}
        
        return adata if copy else None
    else:
        # Direct matrix input
        if CODON_AVAILABLE:
            x = data.copy() if copy else data
            result, _ = codon_pp.normalize_total_kernel(
                x, target_sum, exclude_highly_expressed, max_fraction
            )
            return result
        else:
            return fallback_scanpy.pp.normalize_total(
                data, target_sum=target_sum,
                exclude_highly_expressed=exclude_highly_expressed,
                max_fraction=max_fraction,
                inplace=False
            )
