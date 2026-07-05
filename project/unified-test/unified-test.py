'''
Unified test that covers the fundamentals of scanpy, comparing performance between python and codon versions.
Imports both scanpy (Python) and scancodon (Codon) separately and runs tests side-by-side.
This file has undergone many iterations for debugging purposes.
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
        ("rank_genes_groups", lambda lib: lib.tl.rank_genes_groups(adata, "leiden", method="t-test")),
    ]

    # Run each test function for both libraries
    for func_name, func_lambda in test_functions:
        print(f"[INFO] Running {func_name}...")
        
        try:
            # Time Python version
            python_time_start = time.perf_counter()
            result_python = func_lambda(sp)
            python_time_end = time.perf_counter()
            python_timings[func_name] = python_time_end - python_time_start
        except Exception as e:
            print(f"  [ERROR] Python version failed: {e}")
            python_timings[func_name] = None
        
        try:
            # Time Codon version
            codon_time_start = time.perf_counter()
            result_codon = func_lambda(sc)
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

def debug_hvg_discrepancies(python_snapshots, codon_snapshots, var_names):
    """
    Prints a side-by-side comparison table of the worst offending genes 
    to pinpoint exactly where the Scanpy vs Scancodon calculation breaks down.
    """
    print("\n" + "="*90)
    print("  HIGHLY VARIABLE GENES: DETAILED DISCREPANCY REPORT  ")
    print("="*90)
    
    # Extract the internal dictionaries from the steps
    py_hvg = python_snapshots.get("highly_variable_genes", {})
    cd_hvg = codon_snapshots.get("highly_variable_genes", {})
    
    py_norm = py_hvg.get("var.dispersions_norm")
    cd_norm = cd_hvg.get("var.dispersions_norm")
    py_means = py_hvg.get("var.means")
    cd_means = cd_hvg.get("var.means")
    py_disps = py_hvg.get("var.dispersions")
    cd_disps = cd_hvg.get("var.dispersions")
    py_hv = py_hvg.get("var.highly_variable")
    cd_hv = cd_hvg.get("var.highly_variable")
    
    if py_norm is None or cd_norm is None:
        print("[ERROR] Could not find highly_variable_genes data in snapshots.")
        return

    import numpy as np
    
    # Calculate absolute difference on normalized dispersions
    abs_diff = np.abs(py_norm - cd_norm)
    worst_indices = np.argsort(abs_diff)[::-1]
    
    print(f"Total Genes Checked: {len(abs_diff)}")
    print(f"Max Absolute Error:  {np.nanmax(abs_diff):.5f}")
    print(f"Total Mismatched Masks (highly_variable): {np.sum(py_hv != cd_hv)} genes\n")
    
    print(f"{'Gene Name':<15} | {'Metric Array':<18} | {'Python (Scanpy)':<16} | {'Codon (Scancodon)':<18} | {'Abs Diff':<10}")
    print("-"*92)
    
    # Inspect the top 5 worst offenders
    for idx in worst_indices[:5]:
        gene_name = var_names[idx] if idx < len(var_names) else f"Idx_{idx}"
        
        print(f"{gene_name:<15} | {'dispersions_norm':<18} | {py_norm[idx]:<16.5f} | {cd_norm[idx]:<18.5f} | {abs_diff[idx]:<10.5f}")
        print(f"{'':<15} | {'means':<18} | {py_means[idx]:<16.5f} | {cd_means[idx]:<18.5f} | {abs(py_means[idx]-cd_means[idx]):<10.5f}")
        print(f"{'':<15} | {'dispersions':<18} | {py_disps[idx]:<16.5f} | {cd_disps[idx]:<18.5f} | {abs(py_disps[idx]-cd_disps[idx]):<10.5f}")
        print(f"{'':<15} | {'highly_variable':<18} | {str(py_hv[idx]):<16} | {str(cd_hv[idx]):<18} | {int(py_hv[idx]!=cd_hv[idx]):<10}")
        print("-"*92)

def debug_pca_discrepancies(python_snapshots, codon_snapshots):
    print("\n" + "="*85)
    print("  PCA: DETAILED DISCREPANCY REPORT  ")
    print("="*85)
 
    # python_snapshots / codon_snapshots are keyed by STEP LABEL (e.g. "pca"),
    # and snapshot() stores flat keys like "obsm.X_pca", not nested dicts.
    py_step = python_snapshots.get("pca", {})
    cd_step = codon_snapshots.get("pca", {})
 
    py_pca = py_step.get("obsm.X_pca")
    cd_pca = cd_step.get("obsm.X_pca")
 
    py_varm = py_step.get("varm.PCs")
    cd_varm = cd_step.get("varm.PCs")
 
    if py_pca is None or cd_pca is None:
        print("[ERROR] Could not find obsm['X_pca'] in snapshots.")
        return
 
    import numpy as np
 
    # 1. Check Variance (if available in uns) — uns.pca is a dict, so this
    # nested .get("variance") is correct as-is.
    py_var = (py_step.get("uns.pca") or {}).get("variance")
    cd_var = (cd_step.get("uns.pca") or {}).get("variance")
 
    if py_var is not None and cd_var is not None:
        print("1. VARIANCE COMPARISON (Eigenvalues)")
        print(f"{'Comp':<5} | {'Python':<15} | {'Codon':<15} | {'Abs Diff':<10}")
        print("-" * 55)
        for i in range(min(10, len(py_var))):
            print(f"{i:<5} | {py_var[i]:<15.4f} | {cd_var[i]:<15.4f} | {abs(py_var[i]-cd_var[i]):<10.6f}")
 
    # 2. Check Eigenvector Sign Flip (varm['PCs'])
    if py_varm is not None and cd_varm is not None:
        print("\n2. PRINCIPAL COMPONENTS (Checking for Sign Ambiguity)")
        print(f"{'Comp':<5} | {'Max Diff (Raw)':<16} | {'Max Diff (Flipped)':<20} | {'Is Flipped?':<12}")
        print("-" * 65)
 
        for i in range(min(10, py_varm.shape[1])):
            py_col = py_varm[:, i]
            cd_col = cd_varm[:, i]
 
            diff_raw = np.max(np.abs(py_col - cd_col))
            diff_flipped = np.max(np.abs(py_col - (-cd_col)))
 
            is_flipped = diff_flipped < diff_raw
            print(f"{i:<5} | {diff_raw:<16.4f} | {diff_flipped:<20.4f} | {str(is_flipped):<12}")
 
    # 3. Principal angle subspace check — independent of sign/permutation.
    # PCA components are only unique up to sign (and rotation within
    # degenerate subspaces), so per-component comparisons above can look bad
    # even when both implementations found the same subspace. This check
    # answers the real question: do X_pca (Python) and X_pca (Codon) span
    # the same k-dimensional subspace?
    #
    # Method: QR-orthonormalize both projection matrices, then take the SVD
    # of Q_py.T @ Q_cd — the resulting singular values are cosines of the
    # principal angles between the two subspaces. All ~1.0 (angle ~0 deg)
    # means the subspaces agree; values well below 1.0 mean a real
    # algorithmic discrepancy, not just sign/scale noise.
    print("\n3. SUBSPACE AGREEMENT (Principal Angles — sign/scale independent)")
    try:
        k = min(py_pca.shape[1], cd_pca.shape[1])
        Q_py, _ = np.linalg.qr(py_pca[:, :k])
        Q_cd, _ = np.linalg.qr(cd_pca[:, :k])
 
        cosines = np.linalg.svd(Q_py.T @ Q_cd, compute_uv=False)
        cosines = np.clip(cosines, -1.0, 1.0)
        angles_deg = np.degrees(np.arccos(cosines))
 
        print(f"{'Comp':<5} | {'Cosine':<10} | {'Angle (deg)':<12}")
        print("-" * 35)
        for i in range(len(angles_deg)):
            print(f"{i:<5} | {cosines[i]:<10.6f} | {angles_deg[i]:<12.4f}")
 
        print(f"\nMax principal angle: {angles_deg.max():.4f} deg")
        if angles_deg.max() < 2.0:
            print("  -> Subspaces agree closely. Discrepancies above are sign/scale artifacts.")
        elif angles_deg.max() < 15.0:
            print("  -> Mild subspace drift — check randomized SVD convergence (q, p params).")
        else:
            print("  -> SIGNIFICANT subspace mismatch. This points to an algorithmic bug")
            print("     (e.g. centering, SVD computation), not just sign ambiguity.")
    except Exception as exc:
        print(f"[ERROR] Could not compute principal angles: {exc}")
 
    print("="*85 + "\n")

def benchmark_3k_PBMCs():
    """Benchmark the 3k PBMC dataset from 10x Genomics."""
    import numpy as np

    print("=" * 80)
    print("--- BENCHMARK: 3k PBMCs ---")
    print("=" * 80)

    sp, sc = setup_imports()

    adata = sp.datasets.pbmc3k()
    adata.var_names_make_unique()

    # Each step is a (label, callable) pair so we can time them individually.
    def get_steps(lib):
        return [
            ("calculate_qc_metrics",  lambda a: lib.pp.calculate_qc_metrics(a)),
            ("filter_cells",          lambda a: lib.pp.filter_cells(a, min_genes=100)),
            ("filter_genes",          lambda a: lib.pp.filter_genes(a, min_cells=3)),
            ("scrublet",              lambda a: lib.pp.scrublet(a, random_state=0)),
            ("normalize_total",       lambda a: lib.pp.normalize_total(a)),
            ("log1p",                 lambda a: lib.pp.log1p(a)),
            ("highly_variable_genes", lambda a: lib.pp.highly_variable_genes(a, n_top_genes=2000)),
            ("scale",                 lambda a: lib.pp.scale(a, zero_center = False, max_value=10)),
            ("pca",                   lambda a: lib.tl.pca(a)),
            ("neighbors",             lambda a: lib.pp.neighbors(a, n_neighbors=10, n_pcs=40)),
            ("umap",                  lambda a: lib.tl.umap(a)),
            ("tsne",                  lambda a: lib.tl.tsne(a)),
            ("diffmap",               lambda a: lib.tl.diffmap(a)),
            ("leiden",                lambda a: lib.tl.leiden(a, resolution=0.5)),
            ("rank_genes_groups",     lambda a: lib.tl.rank_genes_groups(a, "leiden", method="t-test")),
        ]

    # ------------------------------------------------------------------
    # Keys written to adata by each step — used for correctness checks.
    # Each entry lists the (namespace, key) pairs to extract after the
    # step runs.  "obs" / "var" entries are column names; "obsm" / "varm"
    # / "uns" / "obsp" are dict keys.
    # ------------------------------------------------------------------
    STEP_OUTPUT_KEYS = {
        "filter_cells":          [],   # shape change only — checked via adata.shape
        "filter_genes":          [],   # shape change only
        "scrublet":              [("obs", "doublet_score"), ("obs", "predicted_doublet")],
        "normalize_total":       [("X", None)],
        "log1p":                 [("X", None)],
        "highly_variable_genes": [("var", "highly_variable"), ("var", "means"),
                                  ("var", "dispersions"), ("var", "dispersions_norm")],
        "scale":                 [("X", None)],
        "pca":                   [("obsm", "X_pca"), ("varm", "PCs"),
                                  ("uns",  "pca")],
        "neighbors":             [("obsp", "connectivities"), ("obsp", "distances")],
        "umap":                  [("obsm", "X_umap")],
        "tsne":                  [("obsm", "X_tsne")],
        "diffmap":               [("obsm", "X_diffmap"), ("uns", "diffmap_evals")],
        "leiden":                [("obs",  "leiden")],
        "rank_genes_groups":     [("uns",  "rank_genes_groups")],
    }

    def snapshot(a, step_label):
        """
        Capture the outputs relevant to *step_label* from AnnData *a*.
        Returns a dict of {key_description: value}.
        """
        result = {"shape": a.shape}
        for namespace, key in STEP_OUTPUT_KEYS.get(step_label, []):
            try:
                if namespace == "X":
                    val = a.X
                    result["X"] = val.toarray() if hasattr(val, "toarray") else np.array(val)
                elif namespace == "obs":
                    result[f"obs.{key}"] = a.obs[key].values.copy()
                elif namespace == "var":
                    result[f"var.{key}"] = a.var[key].values.copy()
                elif namespace == "obsm":
                    result[f"obsm.{key}"] = np.array(a.obsm[key])
                elif namespace == "varm":
                    result[f"varm.{key}"] = np.array(a.varm[key])
                elif namespace == "obsp":
                    val = a.obsp[key]
                    result[f"obsp.{key}"] = val.toarray() if hasattr(val, "toarray") else np.array(val)
                elif namespace == "uns":
                    result[f"uns.{key}"] = a.uns.get(key)
            except (KeyError, AttributeError) as exc:
                result[f"{namespace}.{key}"] = f"<missing: {exc}>"
        return result

    def compare_snapshots(step_label, py_snap, cd_snap, rtol=1e-4, atol=1e-6):
        """
        Compare two snapshots and return a list of deviation strings.
        Empty list means everything matched.
        """
        deviations = []

        # Shape
        if py_snap["shape"] != cd_snap["shape"]:
            deviations.append(
                f"shape mismatch: Python={py_snap['shape']}  Codon={cd_snap['shape']}"
            )

        all_keys = set(py_snap) | set(cd_snap)
        for k in sorted(all_keys):
            if k == "shape":
                continue

            if k not in py_snap:
                deviations.append(f"{k}: present in Codon but missing in Python")
                continue
            if k not in cd_snap:
                deviations.append(f"{k}: present in Python but missing in Codon")
                continue

            pv, cv = py_snap[k], cd_snap[k]

            # Both missing / error strings
            if isinstance(pv, str) and pv.startswith("<missing"):
                deviations.append(f"{k}: Python could not read value ({pv})")
                continue
            if isinstance(cv, str) and cv.startswith("<missing"):
                deviations.append(f"{k}: Codon could not read value ({cv})")
                continue

            # Numeric arrays
            if isinstance(pv, np.ndarray) and isinstance(cv, np.ndarray):
                if pv.shape != cv.shape:
                    deviations.append(
                        f"{k}: shape mismatch Python={pv.shape}  Codon={cv.shape}"
                    )
                else:
                    if not np.allclose(pv, cv, rtol=rtol, atol=atol, equal_nan=True):
                        max_diff = np.nanmax(np.abs(pv.astype(float) - cv.astype(float)))
                        deviations.append(
                            f"{k}: numeric deviation (max |diff|={max_diff:.3e}, "
                            f"rtol={rtol}, atol={atol})"
                        )
                continue

            # Categorical / string arrays (e.g. leiden labels)
            if isinstance(pv, np.ndarray) and pv.dtype.kind in ("U", "O"):
                if not np.array_equal(pv, cv):
                    n_diff = int(np.sum(pv != cv))
                    deviations.append(
                        f"{k}: {n_diff}/{len(pv)} label(s) differ"
                    )
                continue

            # uns dicts — shallow structural check
            if isinstance(pv, dict) and isinstance(cv, dict):
                py_keys = set(pv.keys())
                cd_keys = set(cv.keys())
                if py_keys != cd_keys:
                    deviations.append(
                        f"{k}: key sets differ  "
                        f"Python-only={py_keys - cd_keys}  "
                        f"Codon-only={cd_keys - py_keys}"
                    )
                continue

            # Fallback: direct equality
            try:
                if pv != cv:
                    deviations.append(f"{k}: values differ  Python={pv!r}  Codon={cv!r}")
            except Exception:
                pass  # comparison not meaningful; skip

        return deviations

    def run_pipeline_timed_with_snapshots(lib, a):
        """
        Run pipeline steps sequentially.
        Returns ({step: elapsed}, {step: snapshot}).
        """
        timings   = {}
        snapshots = {}
        for label, step_fn in get_steps(lib):
            t0 = time.perf_counter()
            try:
                step_fn(a)
                timings[label] = time.perf_counter() - t0
            except Exception as exc:
                timings[label] = None
                print(f"  [ERROR] {label} raised: {exc}")
            snapshots[label] = snapshot(a, label)
        return timings, snapshots

    # Run both pipelines on independent copies of the data
    print("[INFO] Running Python scanpy benchmark...")
    python_timings, python_snapshots = run_pipeline_timed_with_snapshots(sp, adata.copy())

    print("[INFO] Running Codon scancodon benchmark...")
    codon_timings, codon_snapshots = run_pipeline_timed_with_snapshots(sc, adata.copy())
    # Compare scale() output statistically — raw arrays are too large/noisy
    # to eyeball for a multiplicative scale discrepancy. We need mean, std,
    # max abs value, and total sum-of-squares (proportional to total variance,
    # which is exactly what feeds into PCA's eigenvalues).
    # print("\n" + "="*60)
    # print("  SCALE STEP: Python vs Codon — statistical comparison")
    # print("="*60)
    # cd_X = codon_snapshots.get("scale", {}).get("X")
    # py_X = python_snapshots.get("scale", {}).get("X")
    # if cd_X is not None and py_X is not None:
    #     cd_X = np.asarray(cd_X)
    #     py_X = np.asarray(py_X)
    #     print(f"{'Metric':<20} | {'Python':<15} | {'Codon':<15} | {'Ratio (cd/py)':<15}")
    #     print("-" * 70)
    #     py_mean = float(np.mean(py_X));      cd_mean = float(np.mean(cd_X))
    #     py_std  = float(np.std(py_X));       cd_std  = float(np.std(cd_X))
    #     py_max  = float(np.max(np.abs(py_X))); cd_max = float(np.max(np.abs(cd_X)))
    #     py_ss   = float(np.sum(py_X ** 2));  cd_ss   = float(np.sum(cd_X ** 2))
    #     for name, py_v, cd_v in [
    #         ("mean",          py_mean, cd_mean),
    #         ("std",           py_std,  cd_std),
    #         ("max abs value", py_max,  cd_max),
    #         ("sum of squares",py_ss,   cd_ss),
    #     ]:
    #         ratio = cd_v / py_v if py_v != 0 else float("nan")
    #         print(f"{name:<20} | {py_v:<15.4f} | {cd_v:<15.4f} | {ratio:<15.4f}")
    #     print(f"\nShapes — Python: {py_X.shape}, Codon: {cd_X.shape}")
    # else:
    #     print("[ERROR] Could not find 'scale' snapshot for one or both pipelines.")
    # print("="*60 + "\n")

    # var_names = adata.var_names.tolist()
    # debug_hvg_discrepancies(python_snapshots, codon_snapshots, var_names)
    # debug_pca_discrepancies(python_snapshots, codon_snapshots)
    # ------------------------------------------------------------------
    # Correctness report
    # ------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("CORRECTNESS REPORT")
    print("=" * 80)

    all_steps = [label for label, _ in get_steps(sp)]
    any_deviation = False

    for label in all_steps:
        py_snap = python_snapshots.get(label)
        cd_snap = codon_snapshots.get(label)

        py_failed = python_timings.get(label) is None
        cd_failed = codon_timings.get(label) is None

        if py_failed and cd_failed:
            print(f"  [SKIP ] {label:<25}  both versions failed — no comparison possible")
            continue
        if py_failed:
            print(f"  [SKIP ] {label:<25}  Python failed — cannot use as reference")
            continue
        if cd_failed:
            print(f"  [FAIL ] {label:<25}  Codon step raised an exception")
            any_deviation = True
            continue

        deviations = compare_snapshots(label, py_snap, cd_snap)
        if deviations:
            any_deviation = True
            print(f"  [FAIL ] {label:<25}  {len(deviations)} deviation(s):")
            for d in deviations:
                print(f"            • {d}")
        else:
            print(f"  [OK   ] {label:<25}  outputs match")

    if not any_deviation:
        print("\n  All steps produced matching outputs.")
    else:
        print("\n  One or more steps deviated — see details above.")

    # ------------------------------------------------------------------
    # Timing report
    # ------------------------------------------------------------------
    python_total = sum(t for t in python_timings.values() if t is not None)
    codon_total  = sum(t for t in codon_timings.values()  if t is not None)

    print("\n" + "=" * 80)
    print("BENCHMARK RESULTS")
    print("=" * 80)
    print(f"{'Step':<25} {'Python (s)':<14} {'Codon (s)':<14} {'Speedup':<10}")
    print("-" * 63)
    for label in all_steps:
        pt = python_timings.get(label)
        ct = codon_timings.get(label)
        pt_str  = f"{pt:.3f}" if pt is not None else "ERROR"
        ct_str  = f"{ct:.3f}" if ct is not None else "ERROR"
        speedup = f"{pt / ct:.2f}x" if (pt and ct) else "N/A"
        print(f"{label:<25} {pt_str:<14} {ct_str:<14} {speedup:<10}")
    print("-" * 63)
    overall_speedup = f"{python_total / codon_total:.2f}x" if codon_total else "N/A"
    print(f"{'TOTAL':<25} {python_total:<14.3f} {codon_total:<14.3f} {overall_speedup:<10}")
    
 
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
    benchmark_3k_PBMCs()