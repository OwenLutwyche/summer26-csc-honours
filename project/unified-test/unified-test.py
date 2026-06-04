'''
Unified test that covers the fundamentals of scanpy, comparing performance between python and codon versions.
Imports both scanpy (Python) and scancodon (Codon) separately and runs tests side-by-side.
'''
import sys
import os
import importlib.util
import time
import anndata
import pooch

# setup: import both scanpy and scancodon

def setup_imports():
    """
    Set up both scanpy (Python) and scancodon (Codon) imports separately.
    Returns a tuple of (scanpy_module, scancodon_module).
    """
    
    # python scanpy
    # Add scanpy-main/src to path so we can import scanpy
    scanpy_path = "/Users/oweno/Desktop/honours/summer26-csc-honours/project/scanpy-main"
    scanpy_src = os.path.join(scanpy_path, "src")
    
    if scanpy_src not in sys.path:
        sys.path.insert(0, scanpy_src)
    
    os.environ["MPLBACKEND"] = "Agg"
    
    try:
        import scanpy as sp
        print(f"[OK] Loaded Python scanpy from: {os.path.dirname(sp.__file__)}")
    except ImportError as e:
        print(f"[ERROR] Could not import Python scanpy: {e}")
        sys.exit(1)
    
    # scancodon
    # Change to scancodon directory and import scancodon
    scancodon_path = "/Users/oweno/Desktop/honours/summer26-csc-honours/project/scancodon"
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
    
    return sp, sc

def load_test_module(filepath: str, module_prefix: str = "test_mod"):
    """Dynamically load a Python module from a file path."""
    module_name = module_prefix + "_" + os.path.basename(filepath).replace(".py", "")
    spec = importlib.util.spec_from_file_location(module_name, filepath)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load module from {filepath}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module

def scanpy_tutorial_test_suite():
    """Run the scanpy tutorial steps on both Python scanpy and Codon scancodon."""

    per_file_timings =[]
    python_timings = {}
    codon_timings = {}
    print("=" * 80)
    print("--- SCANPY TUTORIAL TEST SUITE ---")
    print("=" * 80)

    
    sp, sc = setup_imports()
    

    print("[INFO] Importing example data...")
    EXAMPLE_DATA = pooch.create(
    path=pooch.os_cache("scverse_tutorials"),
    base_url="doi:10.6084/m9.figshare.22716739.v1/",
    )
    EXAMPLE_DATA.load_registry_from_doi()
    samples = {
    "s1d1": "s1d1_filtered_feature_bc_matrix.h5",
    "s1d3": "s1d3_filtered_feature_bc_matrix.h5",
    }
    adatas = {}

    # marker gene set
    marker_genes = {
        "CD14+ Mono": ["FCN1", "CD14"],
        "CD16+ Mono": ["TCF7L2", "FCGR3A", "LYN"],
        # Note: DMXL2 should be negative
        "cDC2": ["CST3", "COTL1", "LYZ", "DMXL2", "CLEC10A", "FCER1A"],
        "Erythroblast": ["MKI67", "HBA1", "HBB"],
        # Note HBM and GYPA are negative markers
        "Proerythroblast": ["CDK6", "SYNGR1", "HBM", "GYPA"],
        "NK": ["GNLY", "NKG7", "CD247", "FCER1G", "TYROBP", "KLRG1", "FCGR3A"],
        "ILC": ["ID2", "PLCG2", "GNLY", "SYNE1"],
        "Naive CD20+ B": ["MS4A1", "IL4R", "IGHD", "FCRL1", "IGHM"],
        # Note IGHD and IGHM are negative markers
        "B cells": [
            "MS4A1",
            "ITGB1",
            "COL4A4",
            "PRDM1",
            "IRF4",
            "PAX5",
            "BCL11A",
            "BLK",
            "IGHD",
            "IGHM",
        ],
        "Plasma cells": ["MZB1", "HSP90B1", "FNDC3B", "PRDM1", "IGKC", "JCHAIN"],
        # Note PAX5 is a negative marker
        "Plasmablast": ["XBP1", "PRDM1", "PAX5"],
        "CD4+ T": ["CD4", "IL7R", "TRBC2"],
        "CD8+ T": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "GZMB", "GZMH", "GZMA"],
        "T naive": ["LEF1", "CCR7", "TCF7"],
        "pDC": ["GZMB", "IL3RA", "COBLL1", "TCF4"],
    }

    for sample_id, filename in samples.items():
        path = EXAMPLE_DATA.fetch(filename)
        sample_adata = sp.read_10x_h5(path)
        sample_adata.var_names_make_unique()
        adatas[sample_id] = sample_adata

    adata = anndata.concat(adatas, label="sample")
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    print(adata.obs["sample"].value_counts())
    adata


    # quality control
    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    # Saving count data for later steps
    adata.layers["counts"] = adata.X.copy()
    # Define test functions with their parameters
    # those that do not work in codon
    test_functions = [
        ("calculate_qc_metrics", lambda lib: lib.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)),
        ("filter_cells", lambda lib: lib.pp.filter_cells(adata, min_genes=100)),
        ("filter_genes", lambda lib: lib.pp.filter_genes(adata, min_cells=3)),
        ("scrublet", lambda lib: lib.pp.scrublet(adata, batch_key="sample")),
        ("normalize_total", lambda lib: lib.pp.normalize_total(adata)),
        ("log1p", lambda lib: lib.pp.log1p(adata)),
        ("highly_variable_genes", lambda lib: lib.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")),
        ("pca", lambda lib: lib.tl.pca(adata)),
        ("neighbors", lambda lib: lib.pp.neighbors(adata, n_neighbors=10, n_pcs=40)),
        ("umap", lambda lib: lib.tl.umap(adata)),
        ("leiden_clustering", lambda lib: lib.tl.leiden(adata, resolution=0.5)),
        ("rank_genes_groups", lambda lib: lib.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")), # should be leiden but not implemented
    ]

    # Run each test function for both libraries
    for func_name, func_lambda in test_functions:
        print(f"[INFO] Running {func_name}...")
        
        # Time Python version
        python_time_start = time.perf_counter()
        try:
            func_lambda(sp)
            python_time_end = time.perf_counter()
            python_timings[func_name] = python_time_end - python_time_start
        except Exception as e:
            print(f"  [ERROR] Python version failed: {e}")
            python_timings[func_name] = None
        
        # Time Codon version
        codon_time_start = time.perf_counter()
        try:
            func_lambda(sc)
            codon_time_end = time.perf_counter()
            codon_timings[func_name] = codon_time_end - codon_time_start
        except Exception as e:
            print(f"  [ERROR] Codon version failed: {e}")
            codon_timings[func_name] = None
    
    # Format results for print_comparison
    python_passed = sum(1 for t in python_timings.values() if t is not None)
    python_failed = sum(1 for t in python_timings.values() if t is None)
    python_total_time = sum(t for t in python_timings.values() if t is not None)
    python_per_file_timings = [(name, t) for name, t in python_timings.items() if t is not None]
    
    codon_passed = sum(1 for t in codon_timings.values() if t is not None)
    codon_failed = sum(1 for t in codon_timings.values() if t is None)
    codon_total_time = sum(t for t in codon_timings.values() if t is not None)
    codon_per_file_timings = [(name, t) for name, t in codon_timings.items() if t is not None]
    
    python_results = {
        "lib_name": "Python",
        "total_passed": python_passed,
        "total_failed": python_failed,
        "total_elapsed": python_total_time,
        "per_file_timings": python_per_file_timings,
    }
    
    codon_results = {
        "lib_name": "Codon",
        "total_passed": codon_passed,
        "total_failed": codon_failed,
        "total_elapsed": codon_total_time,
        "per_file_timings": codon_per_file_timings,
    }
    
    # Print comparison summary
    print()
    print("=" * 80)
    print("COMPARISON SUMMARY")
    print("=" * 80)
    print_comparison(python_results, codon_results)

def Preprocessing_and_clustering_3k_PBMCs():
    pass

def new_test_suite():
    print("=" * 80)
    print("--- UNIFIED BASIC TEST: Python Scanpy vs Codon Scancodon ---")
    print("=" * 80)
    print()
    
    # Import both libraries
    sp, sc = setup_imports()
    print()
    
    # Verify both imported successfully and separately
    print("Verification:")
    print(f"  Python scanpy module: {sp.__name__}")
    print(f"  Codon scancodon module: {sc.__name__}")
    print(f"  Are they different? {sp is not sc}")
    print()
    
    # Define test files for Python scanpy
    python_test_files = [
        "/Users/oweno/Desktop/honours/summer26-csc-honours/project/scanpy-main/new_tests/test_preprocessing.py",
        "/Users/oweno/Desktop/honours/summer26-csc-honours/project/scanpy-main/new_tests/test_neighbors.py",
        "/Users/oweno/Desktop/honours/summer26-csc-honours/project/scanpy-main/new_tests/test_clustering.py",
        "/Users/oweno/Desktop/honours/summer26-csc-honours/project/scanpy-main/new_tests/test_embedding.py",
        "/Users/oweno/Desktop/honours/summer26-csc-honours/project/scanpy-main/new_tests/test_rank_genes_groups.py",
    ]
    
    # Define test files for Codon scancodon
    codon_test_files = [
        "/Users/oweno/Desktop/honours/summer26-csc-honours/project/scancodon/tests/test_preprocessing.py",
        "/Users/oweno/Desktop/honours/summer26-csc-honours/project/scancodon/tests/test_neighbors.py",
        "/Users/oweno/Desktop/honours/summer26-csc-honours/project/scancodon/tests/test_clustering.py",
        "/Users/oweno/Desktop/honours/summer26-csc-honours/project/scancodon/tests/test_embedding.py",
        "/Users/oweno/Desktop/honours/summer26-csc-honours/project/scancodon/tests/test_rank_genes_groups.py",
    ]
    
    # Run tests for Python scanpy
    print("=" * 80)
    print("RUNNING TESTS: Python Scanpy")
    print("=" * 80)
    python_results = run_tests_for_library(python_test_files, "Python", sp)
    
    # Run tests for Codon scancodon
    print()
    print("=" * 80)
    print("RUNNING TESTS: Codon Scancodon")
    print("=" * 80)
    codon_results = run_tests_for_library(codon_test_files, "Codon", sc)
    
    # Print comparison summary
    print()
    print("=" * 80)
    print("COMPARISON SUMMARY")
    print("=" * 80)
    print_comparison(python_results, codon_results)


def run_tests_for_library(test_files, lib_name, lib_module):
    """Run all test files for a given library."""
    total_passed = 0
    total_failed = 0
    start_time = time.time()
    per_file_timings = []
    
    for test_file in test_files:
        test_name = os.path.basename(test_file)
        
        if not os.path.exists(test_file):
            print(f"\n[WARN] Test file not found: {test_file}")
            continue
        
        print(f"\n[INFO] Running {test_name}...")
        file_start = time.perf_counter()
        
        try:
            # Load the test module dynamically
            module = load_test_module(test_file, f"{lib_name}_test")
            
            if hasattr(module, "run_all"):
                results = module.run_all()
                for name, success, msg in results:
                    if success:
                        print(f"  [PASS] {name}")
                        total_passed += 1
                    else:
                        print(f"  [FAIL] {name}: {msg}")
                        total_failed += 1
            else:
                print("  [WARN] No run_all() function found in test file.")
        except Exception as e:
            print(f"  [CRASH] {e}")
            import traceback
            traceback.print_exc()
            total_failed += 1
        
        file_elapsed = time.perf_counter() - file_start
        per_file_timings.append((test_name, file_elapsed))
        print(f"  [TIME] {file_elapsed:.2f}s")
    
    total_elapsed = time.time() - start_time
    
    return {
        "lib_name": lib_name,
        "total_passed": total_passed,
        "total_failed": total_failed,
        "total_elapsed": total_elapsed,
        "per_file_timings": per_file_timings,
    }


def print_comparison(python_results, codon_results):
    """Print a side-by-side comparison of results."""
    print(f"\n{'Metric':<25} {'Python':<20} {'Codon':<20}")
    print("-" * 65)
    print(f"{'Passed Tests':<25} {python_results['total_passed']:<20} {codon_results['total_passed']:<20}")
    print(f"{'Failed Tests':<25} {python_results['total_failed']:<20} {codon_results['total_failed']:<20}")
    print(f"{'Total Time (s)':<25} {python_results['total_elapsed']:<20.2f} {codon_results['total_elapsed']:<20.2f}")
    
    if python_results["total_elapsed"] > 0:
        speedup = codon_results["total_elapsed"] / python_results["total_elapsed"]
        print(f"{'Speedup Factor':<25} {'1.0x':<20} {f'{speedup:.2f}x':<20}")
    
    print()
    print("Per-file timings:")
    print(f"{'Test File':<35} {'Python (s)':<15} {'Codon (s)':<15}")
    print("-" * 65)
    
    # Match up files by name for comparison
    python_times = {name: t for name, t in python_results["per_file_timings"]}
    codon_times = {name: t for name, t in codon_results["per_file_timings"]}
    
    all_files = set(python_times.keys()) | set(codon_times.keys())
    for test_name in sorted(all_files):
        python_t = python_times.get(test_name, 0)
        codon_t = codon_times.get(test_name, 0)
        print(f"{test_name:<35} {python_t:<15.2f} {codon_t:<15.2f}")


if __name__ == "__main__":
    scanpy_tutorial_test_suite()

