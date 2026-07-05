'''
Test to run scancodon.scrublet
'''
import sys
import os
import time
import scanpy as sp

# setup: import both scanpy and scancodon

def setup_imports():
    """
    Set up both scanpy (Python) and scancodon (Codon) imports separately.
    Returns a tuple of (scanpy_module, scancodon_module).
    """

    
    # scancodon
    # Change to scancodon directory and import scancodon
    test_dir = os.path.dirname(os.path.abspath(__file__))
    scancodon_path = os.path.abspath(os.path.join(test_dir, "..", "scancodon"))
    original_cwd = os.getcwd()
    original_cwd = os.getcwd()
    
    try:
        os.chdir(scancodon_path)
        if scancodon_path not in sys.path:
            sys.path.insert(0, scancodon_path)
        
        import scancodon as sc
        print(f"[OK] Loaded Codon scancodon from: {os.path.dirname(sc.__file__)}")
        
        # Check if native Codon kernels are available
        if getattr(sc, "CODON_AVAILABLE", False):
            print("[STATUS] Running on NATIVE CODON KERNELS")
        else:
            print("[WARNING] Running on NUMPY FALLBACKS (Native extension not loaded)")
    except ImportError as e:
        print(f"[ERROR] Could not import scancodon: {e}")
        sys.exit(1)
    finally:
        os.chdir(original_cwd)
    
    return sc

def test_isolated_scrublet():
    """
    Isolated test for scrublet. Uses ONLY Python scanpy for dataset loading 
    and preliminary filtering, ensuring both scrublet implementations receive 
    the exact same starting CSR matrix.
    """
    import numpy as np
    
    print("=" * 80)
    print("--- ISOLATED BENCHMARK: Scrublet ---")
    print("=" * 80)

    sc = setup_imports()

    # 1. Load data
    print("[INFO] Loading 3k PBMC dataset (raw counts)...")
    adata = sp.datasets.pbmc3k()
    adata.var_names_make_unique()
    
    # 2. Preprocessing (Using ONLY scanpy to isolate the scancodon test)
    print("[INFO] Preprocessing with Scanpy (filtering)...")
    sp.pp.filter_cells(adata, min_genes=200)
    sp.pp.filter_genes(adata, min_cells=3)
    
    # Ensure it's CSR sparse before passing to Codon
    import scipy.sparse as sp_sparse
    if not sp_sparse.isspmatrix_csr(adata.X):
        adata.X = adata.X.tocsr()
    
    # Split into two independent copies
    adata_py = adata.copy()
    adata_cd = adata.copy()

    # 3. Run Python Scrublet
    print("\n[INFO] Running Python scanpy.pp.scrublet...")
    py_start = time.perf_counter()
    try:
        # random_state ensures reproducible doublet simulations
        sp.pp.scrublet(adata_py, random_state=0)
        py_time = time.perf_counter() - py_start
        print(f"  [OK] Python finished in {py_time:.3f}s")
    except Exception as e:
        print(f"  [ERROR] Python scrublet failed: {e}")
        py_time = None

    # 4. Run Codon Scrublet
    print("\n[INFO] Running Codon scancodon.pp.scrublet...")
    cd_start = time.perf_counter()
    try:
        sc.pp.scrublet(adata_cd, random_state=0)
        cd_time = time.perf_counter() - cd_start
        print(f"  [OK] Codon finished in {cd_time:.3f}s")
    except Exception as e:
        print(f"  [ERROR] Codon scrublet failed: {e}")
        cd_time = None

    # 5. Compare Results
    print("\n" + "=" * 80)
    print("COMPARISON: Scrublet Outputs")
    print("=" * 80)

    if py_time is not None and cd_time is not None:
        # Compare doublet scores
        py_scores = np.asarray(adata_py.obs.get('doublet_score', []))
        cd_scores = np.asarray(adata_cd.obs.get('doublet_score', []))
        
        if py_scores.size > 0 and cd_scores.size > 0:
            abs_diff = np.abs(py_scores - cd_scores)
            max_diff = np.nanmax(abs_diff)
            mean_diff = np.nanmean(abs_diff)
            
            print(f"{'Metric':<35} {'Value'}")
            print("-" * 50)
            print(f"{'Doublet Score Max Abs Diff:':<35} {max_diff:.6f}")
            print(f"{'Doublet Score Mean Abs Diff:':<35} {mean_diff:.6f}")
            
            if max_diff < 1e-4:
                print("  -> [PASS] Doublet scores match closely!")
            else:
                print("  -> [FAIL] Significant discrepancy in doublet scores.")
        else:
            print("  -> [ERROR] Doublet scores missing from output.")
            
        # Compare boolean predictions
        py_pred = np.asarray(adata_py.obs.get('predicted_doublet', []))
        cd_pred = np.asarray(adata_cd.obs.get('predicted_doublet', []))
        
        if py_pred.size > 0 and cd_pred.size > 0:
            mismatches = np.sum(py_pred != cd_pred)
            total = len(py_pred)
            print(f"\n{'Predicted Doublets Mismatches:':<35} {mismatches} / {total} cells")
            
            if mismatches == 0:
                print("  -> [PASS] Boolean predictions match exactly!")
            else:
                print("  -> [FAIL] Boolean predictions differ.")
        
        # Performance
        speedup = py_time / cd_time
        print(f"\nSpeedup (Python/Codon): {speedup:.2f}x")
        
    else:
        print("Comparison skipped due to execution failure (likely your loop bug crashing the Codon step).")



if __name__ == "__main__":
    test_isolated_scrublet()