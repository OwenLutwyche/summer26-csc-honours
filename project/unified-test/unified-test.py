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
import numpy as np

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

def pipeline_benchmark_3k_PBMCs():
    """Benchmark the 3k PBMC dataset from 10x Genomics using case-by-case evaluation rules."""
    import numpy as np

    print("=" * 80)
    print("--- BENCHMARK: 3k PBMCs (Correctness Verification by Type of algorithm) ---")
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
            ("scale",                 lambda a: lib.pp.scale(a, max_value=10)),
            ("pca",                   lambda a: lib.tl.pca(a)),
            ("neighbors",             lambda a: lib.pp.neighbors(a, n_neighbors=10, n_pcs=40)),
            ("umap",                  lambda a: lib.tl.umap(a)),
            ("tsne",                  lambda a: lib.tl.tsne(a)),
            ("diffmap",               lambda a: lib.tl.diffmap(a)),
            ("leiden",                lambda a: lib.tl.leiden(a, resolution=0.5)),
            ("rank_genes_groups",     lambda a: lib.tl.rank_genes_groups(a, "leiden", method="t-test")),
        ]

    STEP_OUTPUT_KEYS = {
        "filter_cells":          [],   
        "filter_genes":          [],   
        "scrublet":              [("obs", "doublet_score"), ("obs", "predicted_doublet")],
        "normalize_total":       [("X", None)],
        "log1p":                 [("X", None)],
        "highly_variable_genes": [("var", "highly_variable"), ("var", "means"),
                                  ("var", "dispersions"), ("var", "dispersions_norm")],
        "scale":                 [("X", None)],
        "pca":                   [("obsm", "X_pca"), ("varm", "PCs"), ("uns",  "pca")],
        "neighbors":             [("obsp", "connectivities"), ("obsp", "distances")],
        "umap":                  [("obsm", "X_umap")],
        "tsne":                  [("obsm", "X_tsne")],
        "diffmap":               [("obsm", "X_diffmap"), ("uns", "diffmap_evals")],
        "leiden":                [("obs",  "leiden")],
        "rank_genes_groups":     [("uns",  "rank_genes_groups")],
    }

    def snapshot(a, step_label):
        """Capture the outputs relevant to *step_label* from AnnData *a*."""
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

    def compare_snapshots(step_label, py_snap, cd_snap, X_ref=None):
        """Compares snapshots routing components dynamically to their designated evaluation tier."""
        deviations = []

        if py_snap["shape"] != cd_snap["shape"]:
            deviations.append(f"shape mismatch: Python={py_snap['shape']}  Codon={cd_snap['shape']}")

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

            if isinstance(pv, str) and pv.startswith("<missing"):
                deviations.append(f"{k}: Python could not read value ({pv})")
                continue
            if isinstance(cv, str) and cv.startswith("<missing"):
                deviations.append(f"{k}: Codon could not read value ({cv})")
                continue

            err = None

            # Stochastic Manifolds (UMAP & TSNE)
            if "X_umap" in k or "X_tsne" in k:
                err = evaluate_stochastic_manifold(pv, cv, X_ref=X_ref, min_trustworthiness=0.82, max_disparity=0.50)
                if err:
                    deviations.append(f"{k} [Stochastic Manifold]: {err}")

            # Linear Subspaces (PCA Components & Diffusion Maps)
            elif "X_pca" in k or "PCs" in k or "X_diffmap" in k:
                err = evaluate_linear_subspace(pv, cv, max_disparity=1e-2)
                if err:
                    deviations.append(f"{k} [Linear Subspace]: {err}")

            # Graph Topology (Connectivities & Distances Matrices)
            elif "connectivities" in k or "distances" in k:
                err = evaluate_graph_topology(pv, cv, min_jaccard=0.85)
                if err:
                    deviations.append(f"{k} [3a Graph Topology]: {err}")

            # Clustering Groups (Leiden Cluster Maps)
            elif "leiden" in k:
                err = evaluate_clustering(pv, cv, min_ari=0.85)
                if err:
                    deviations.append(f"{k} [3b Clustering]: {err}")

            # Structured Array fallbacks
            elif isinstance(pv, np.ndarray) and isinstance(cv, np.ndarray):
                if pv.shape != cv.shape:
                    deviations.append(f"{k}: shape mismatch Python={pv.shape} Codon={cv.shape}")
                else:
                    if pv.dtype.kind in ("U", "O"):
                        if not np.array_equal(pv, cv):
                            n_diff = int(np.sum(pv != cv))
                            deviations.append(f"{k}: {n_diff}/{len(pv)} label mappings differ")
                    else:
                        err = evaluate_strict(pv, cv, rtol=1e-4, atol=1e-6)
                        if err:
                            deviations.append(f"{k} [Strict Parity]: {err}")

            # Deep verification logic for uns structured metrics (e.g., rank_genes_groups statistics arrays)
            elif isinstance(pv, dict) and isinstance(cv, dict):
                py_keys, cd_keys = set(pv.keys()), set(cv.keys())
                if py_keys != cd_keys:
                    deviations.append(f"{k}: internal dictionary keys differ. Python-only={py_keys - cd_keys} Codon-only={cd_keys - py_keys}")
                else:
                    for subk in sorted(py_keys):
                        p_sub, c_sub = pv[subk], cv[subk]
                        # Verify if fields represent native structured array layouts
                        if hasattr(p_sub, 'dtype') and hasattr(c_sub, 'dtype') and p_sub.dtype.names is not None:
                            for field in p_sub.dtype.names:
                                if field in c_sub.dtype.names:
                                    pf, cf = p_sub[field], c_sub[field]
                                    if pf.dtype.kind in ("U", "O"):
                                        if not np.array_equal(pf, cf):
                                            deviations.append(f"{k}.{subk}['{field}']: label metadata mismatch")
                                    else:
                                        # Yield differential expression stats minor allowance for precision drift
                                        sub_err = evaluate_strict(pf, cf, rtol=1e-2, atol=1e-4)
                                        if sub_err:
                                            deviations.append(f"{k}.{subk}['{field}'] [DE Array Metric Deviation]: {sub_err}")
                        elif isinstance(p_sub, np.ndarray) and isinstance(c_sub, np.ndarray):
                            if p_sub.dtype.kind not in ("U", "O"):
                                sub_err = evaluate_strict(p_sub, c_sub, rtol=1e-3, atol=1e-5)
                                if sub_err:
                                    deviations.append(f"{k}.{subk} [Array Deviation]: {sub_err}")
            else:
                try:
                    if pv != cv:
                        deviations.append(f"{k}: explicit static values differ. Python={pv!r} Codon={cv!r}")
                except Exception:
                    pass

        return deviations

    def run_pipeline_timed_with_snapshots(lib, a):
        """Run pipeline steps sequentially. Returns ({step: elapsed}, {step: snapshot})."""
        timings = {}
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

    # Run both pipelines on independent data frames
    print("[INFO] Running Python scanpy baseline pipeline...")
    python_timings, python_snapshots = run_pipeline_timed_with_snapshots(sp, adata.copy())

    print("[INFO] Running Codon scancodon native pipeline...")
    codon_timings, codon_snapshots = run_pipeline_timed_with_snapshots(sc, adata.copy())

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
            print(f" [SKIP ] {label:<25} both libraries failed execution.")
            continue
        if py_failed:
            print(f" [SKIP ] {label:<25} Baseline Python version errored out.")
            continue
        if cd_failed:
            print(f" [FAIL ] {label:<25} Native Codon layer broke during execution.")
            any_deviation = True
            continue

        # Extract X_pca coordinates from the python baseline to serve as high-dim neighborhood ground truths
        X_ref = python_snapshots.get("pca", {}).get("obsm.X_pca")
        
        deviations = compare_snapshots(label, py_snap, cd_snap, X_ref=X_ref)
        if not deviations:
            print(f" [PASS ] {label:<25} Verified aligned within safe structural bounds.")
        else:
            print(f" [FAIL ] {label:<25} Algorithmic drift uncovered:")
            for dev in deviations:
                print(f"   - {dev}")
            any_deviation = True

    # Format data for summary print structures
    py_passed = sum(1 for t in python_timings.values() if t is not None)
    py_failed = sum(1 for t in python_timings.values() if t is None)
    py_total = sum(t for t in python_timings.values() if t is not None)
    
    cd_passed = sum(1 for t in codon_timings.values() if t is not None)
    cd_failed = sum(1 for t in codon_timings.values() if t is None)
    cd_total = sum(t for t in codon_timings.values() if t is not None)

    py_results = {
        "total_passed": py_passed,
        "total_failed": py_failed,
        "total_elapsed": py_total,
        "per_file_timings": [(name, t) for name, t in python_timings.items() if t is not None]
    }
    cd_results = {
        "total_passed": cd_passed,
        "total_failed": cd_failed,
        "total_elapsed": cd_total,
        "per_file_timings": [(name, t) for name, t in codon_timings.items() if t is not None]
    }

    print("\n" + "=" * 80)
    print("PERFORMANCE BENCHMARK SUMMARY")
    print("=" * 80)
    print_comparison(py_results, cd_results)    


def evaluate_strict(pv, cv, rtol=1e-4, atol=1e-6):
    """
    Strict Determinism (Arithmetic exact matches)
    applies standard close tolerances for deterministic matrix operations, just a straight np.allclose
    should work for log1p, normalize_total, qc_metrics, filter_genes/cells, scale, hvg and scrublet (if it's seeded deterministic-style)
    """
    if pv.shape != cv.shape:
        return f"Shape mismatch: Python {pv.shape} vs Codon {cv.shape}"
    if not np.allclose(pv, cv, rtol=rtol, atol=atol, equal_nan=True):
        max_diff = np.nanmax(np.abs(pv.astype(float) - cv.astype(float)))
        return f"Numeric deviation exceeds tolerance (max |diff| = {max_diff:.3e})"
    return None


def evaluate_linear_subspace(pv, cv, max_disparity=1e-2):
    """
    Linear Subspaces (PCA / Diffmap structural alignment)
    uses procrustes to account for rotation and translation of eigenvectors in low-dimensional spaces.
    fairly evaluates pca and diffmap
    """
    if pv.shape != cv.shape:
        return f"Shape mismatch: Python {pv.shape} vs Codon {cv.shape}"
    
    from scipy.spatial import procrustes
    try:
        _, _, disparity = procrustes(pv, cv)
        if disparity > max_disparity:
            return f"Procrustes disparity {disparity:.4e} exceeds tolerance threshold {max_disparity:.4e}"
    except Exception as e:
        return f"Procrustes calculation failed: {e}"
    return None


def evaluate_graph_topology(pv, cv, min_jaccard=0.85):
    """
    graph tpology tests (Neighbors graph overlap)
    compute row-wise Jaccard similarity of neighborhood intersections to verify that the graph topology matches.
    essentially checks proportion of intersection, good for clustering neighbors and leiden
    """
    if pv.shape != cv.shape:
        return f"Shape mismatch: Python {pv.shape} vs Codon {cv.shape}"
    
    row_jaccards = []
    for i in range(pv.shape[0]):
        py_neighbors = set(np.where(pv[i] > 1e-5)[0])
        co_neighbors = set(np.where(cv[i] > 1e-5)[0])
        
        if not py_neighbors and not co_neighbors:
            row_jaccards.append(1.0)
            continue
            
        intersection = len(py_neighbors.intersection(co_neighbors))
        union = len(py_neighbors.union(co_neighbors))
        row_jaccards.append(intersection / union if union > 0 else 0.0)
        
    avg_jaccard = np.mean(row_jaccards)
    if avg_jaccard < min_jaccard:
        return f"Average row-wise Jaccard similarity {avg_jaccard:.4f} is below threshold {min_jaccard:.4f}"
    return None


def evaluate_clustering(pv, cv, min_ari=0.90):
    """
    clustering parity (leiden label assignments)
    Uses the Adjusted Rand Index (ARI) to determine if partition boundaries are functionally identical, independent of label permutation changes
    """
    import numpy as np
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    
    # 1. Check if the algorithm exploded into fragments or collapsed into a blob.
    # We allow a variance of +/- 1 cluster, as edge-case cells might form micro-clusters
    # depending on floating point drift during the Gumbel distribution roll.
    n_py = len(np.unique(pv))
    n_cd = len(np.unique(cv))
    
    if abs(n_py - n_cd) > 1:
        return f"Cluster count mismatch: Python={n_py}, Codon={n_cd}"
        
    # 2. Calculate structural parity
    ari = adjusted_rand_score(pv, cv)
    
    # 3. Check against threshold
    if ari < min_ari:
        # We calculate NMI only when it fails to provide richer diagnostic output
        nmi = normalized_mutual_info_score(pv, cv)
        return f"Structural divergence: ARI={ari:.4f} (threshold {min_ari:.4f}), NMI={nmi:.4f} | Clusters: Py={n_py}, Codon={n_cd}"
        
    return None


def evaluate_stochastic_manifold(pv, cv, X_ref=None, min_trustworthiness=0.80, max_disparity=0.5):
    """
    Stochastic Manifolds (UMAP / t-SNE neighborhood conservation)
    Combines Procrustes Analysis for global cluster macro-structures and sklearn.manifold.trustworthiness index to score local neighborhood integrity 
    against high-dim PCA space. 
    used for assessing umap and TSNE
    """
    if pv.shape != cv.shape:
        return f"Shape mismatch: Python {pv.shape} vs Codon {cv.shape}"
        
    errors = []
    
    #1. procrustes
    from scipy.spatial import procrustes
    try:
        _, _, disparity = procrustes(pv, cv)
        if disparity > max_disparity:
            errors.append(f"Global Procrustes disparity {disparity:.8f} exceeds threshold {max_disparity:.4f}")
    except Exception as e:
        errors.append(f"Procrustes failed: {e}")
        
    # 2. local neighbourhood conservation (trustworthiness)
    if X_ref is not None:
        from sklearn.manifold import trustworthiness
        try:
            t_score = trustworthiness(X_ref, cv, n_neighbors=15)
            if t_score < min_trustworthiness:
                errors.append(f"Local Trustworthiness score {t_score:.4f} below threshold {min_trustworthiness:.4f}")
        except Exception as e:
            errors.append(f"Trustworthiness benchmark failed: {e}")
            
    if errors:
        return " | ".join(errors)
    return None

def correctness_benchmark_3k_PBMCs():
    """Benchmark the 3k PBMC dataset isolating each step to prevent cascading errors."""
    import numpy as np

    print("=" * 80)
    print("--- CORRECTNESS BENCHMARK: ISOLATED STEPS (TIERED) ---")
    print("=" * 80)

    sp, sc = setup_imports()

    adata = sp.datasets.pbmc3k()
    adata.var_names_make_unique()

    def get_steps(lib):
        return [
            ("calculate_qc_metrics",  lambda a: lib.pp.calculate_qc_metrics(a)),
            ("filter_cells",          lambda a: lib.pp.filter_cells(a, min_genes=100)),
            ("filter_genes",          lambda a: lib.pp.filter_genes(a, min_cells=3)),
            #("scrublet",              lambda a: lib.pp.scrublet(a, random_state=0)), # consider skipping cause it takes too long
            ("normalize_total",       lambda a: lib.pp.normalize_total(a)),
            ("log1p",                 lambda a: lib.pp.log1p(a)),
            ("highly_variable_genes", lambda a: lib.pp.highly_variable_genes(a, n_top_genes=2000)),
            ("scale",                 lambda a: lib.pp.scale(a, max_value=10)),
            ("pca",                   lambda a: lib.tl.pca(a)),
            ("neighbors",             lambda a: lib.pp.neighbors(a, n_neighbors=10, n_pcs=40)),
            ("umap",                  lambda a: lib.tl.umap(a)),
            ("tsne",                  lambda a: lib.tl.tsne(a)),
            ("diffmap",               lambda a: lib.tl.diffmap(a)),
            ("leiden",                lambda a: lib.tl.leiden(a, resolution=0.5)),
            ("rank_genes_groups",     lambda a: lib.tl.rank_genes_groups(a, "leiden", method="t-test")),
        ]

    STEP_OUTPUT_KEYS = {
        "filter_cells":          [],   
        "filter_genes":          [],   
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

    # ==========================================
    # SNAPSHOT & COMPARISON LOGIC
    # ==========================================
    def snapshot(a, step_label):
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

    def compare_snapshots(step_label, py_snap, cd_snap, X_ref=None):
        deviations = []
        if py_snap["shape"] != cd_snap["shape"]:
            deviations.append(f"shape mismatch: Python={py_snap['shape']}  Codon={cd_snap['shape']}")

        all_keys = set(py_snap) | set(cd_snap)
        for k in sorted(all_keys):
            if k == "shape": continue
            if k not in py_snap: deviations.append(f"{k}: present in Codon but missing in Python"); continue
            if k not in cd_snap: deviations.append(f"{k}: present in Python but missing in Codon"); continue

            pv, cv = py_snap[k], cd_snap[k]
            
            # Catch missing values
            if isinstance(pv, str) and pv.startswith("<missing"):
                deviations.append(f"{k}: Python could not read value ({pv})"); continue
            if isinstance(cv, str) and cv.startswith("<missing"):
                deviations.append(f"{k}: Codon could not read value ({cv})"); continue

            err = None




            # -----------------------------------------------------------
            # The Multi-Tier Routing Logic
            # -----------------------------------------------------------
            
            #  Stochastic Manifolds
            if "X_umap" in k or "X_tsne" in k:
                err = evaluate_stochastic_manifold(pv, cv, X_ref=X_ref)
                if err: deviations.append(f"{k}: {err}")

            # Linear Subspaces
            elif "X_pca" in k or "PCs" in k or "X_diffmap" in k:
                err = evaluate_linear_subspace(pv, cv)
                if err: deviations.append(f"{k}: {err}")

            
            # Graph Topology
            elif "connectivities" in k or "distances" in k:
                err = evaluate_graph_topology(pv, cv)
                if err: deviations.append(f"{k}: {err}")
                


            # Clustering Groups
            elif "leiden" in k:
                err = evaluate_clustering(pv, cv)
                if err: deviations.append(f"{k}: {err}")
            

            # Structured Array fallbacks
            elif isinstance(pv, np.ndarray) and isinstance(cv, np.ndarray):
                if pv.dtype.kind in ("U", "O"):
                    if not np.array_equal(pv, cv):
                        n_diff = int(np.sum(pv != cv))
                        deviations.append(f"{k}: {n_diff}/{len(pv)} label mappings differ")
                else:
                    err = evaluate_strict(pv, cv)
                    if err: deviations.append(f"{k}: {err}")
                    
            # -----------------------------------------------------------

            # Structured dictionaries (like rank_genes_groups)
            elif isinstance(pv, dict) and isinstance(cv, dict):
                py_keys, cd_keys = set(pv.keys()), set(cv.keys())
                if py_keys != cd_keys:
                    deviations.append(f"{k}: dictionary keys differ Python-only={py_keys - cd_keys} Codon-only={cd_keys - py_keys}")
                else:
                    for subk in sorted(py_keys):
                        p_sub, c_sub = pv[subk], cv[subk]
                        if hasattr(p_sub, 'dtype') and hasattr(c_sub, 'dtype') and p_sub.dtype.names is not None:
                            for field in p_sub.dtype.names:
                                if field in c_sub.dtype.names:
                                    pf, cf = p_sub[field], c_sub[field]
                                    
                                    if pf.dtype.kind in ("U", "O"):
                                        # For names, just check if the sets of top genes are highly similar
                                        p_set, c_set = set(pf), set(cf)
                                        jaccard = len(p_set & c_set) / len(p_set | c_set) if p_set else 1.0
                                        if jaccard < 0.95:
                                            deviations.append(f"{k}.{subk}['{field}']: Gene set overlap too low (Jaccard={jaccard:.2f})")
                                    else:
                                        # --- THE ALIGNMENT FIX ---
                                        # Sort both the Python and Codon numeric arrays alphabetically by their respective gene names
                                        p_names = pv['names'][field]
                                        c_names = cv['names'][field]
                                        
                                        p_order = np.argsort(p_names)
                                        c_order = np.argsort(c_names)
                                        
                                        # compare genes with meaningful score
                                        mask = np.abs(pf) > 1e-4

                                        # Compare the aligned arrays
                                        sub_err = evaluate_strict(pf[p_order][mask], cf[c_order][mask], rtol=1e-2, atol=1e-4)
                                        if sub_err: 
                                            deviations.append(f"{k}.{subk}['{field}']: {sub_err}")
                        elif isinstance(p_sub, np.ndarray) and isinstance(c_sub, np.ndarray):
                            if p_sub.dtype.kind not in ("U", "O"):
                                sub_err = evaluate_strict(p_sub, c_sub, rtol=1e-3, atol=1e-5)
                                if sub_err: deviations.append(f"{k}.{subk}: {sub_err}")
                if deviations:
                    print(f"\n[DEBUG] Rank Genes Groups Mismatch detected. Probing Group '0':")
                    group_id = '0'
                    try:
                        # Look at the top gene in Python for group '0'
                        py_names = pv['names'][group_id]
                        top_gene = py_names[0]
                        
                        # Find index in Python (which is 0)
                        p_idx = 0
                        
                        # Find where THAT SAME GENE is in Codon
                        co_names = cv['names'][group_id]
                        c_idx = np.where(co_names == top_gene)[0][0]
                        
                        print(f"  Target Gene: '{top_gene}'")
                        print(f"  Python Rank: {p_idx} | Codon Rank: {c_idx}")
                        
                        # Compare Scores
                        py_score = pv['scores'][group_id][p_idx]
                        co_score = cv['scores'][group_id][c_idx]
                        print(f"  Scores: Py={py_score:.4f}, Co={co_score:.4f}")
                        
                        # Compare LogFoldChanges (this is where you were seeing 30.0 deviations)
                        py_lfc = pv['logfoldchanges'][group_id][p_idx]
                        co_lfc = cv['logfoldchanges'][group_id][c_idx]
                        print(f"  LFC:    Py={py_lfc:.4f}, Co={co_lfc:.4f}")
                        
                        # Compare P-values
                        py_pv = pv['pvals'][group_id][p_idx]
                        co_pv = cv['pvals'][group_id][c_idx]
                        print(f"  P-vals: Py={py_pv:.4e}, Co={co_pv:.4e}")
                        
                    except Exception as e:
                        print(f"  [!] Debug probe failed: {e}")
            
            else:
                try:
                    if pv != cv: deviations.append(f"{k}: values differ Python={pv!r} Codon={cv!r}")
                except Exception:
                    pass
            

        return deviations
    # ==========================================
    # EXECUTION PIPELINE
    # ==========================================
    adata_golden = adata.copy()
    python_timings = {}
    codon_timings = {}
    deviations_log = {}
    
    python_steps = get_steps(sp)
    codon_steps_dict = dict(get_steps(sc))

    print("[INFO] Running isolated correctness benchmark...")

    for label, step_fn_py in python_steps:
        step_fn_cd = codon_steps_dict[label]

        # 1. Isolate the current correct state for Codon to test on
        adata_codon_test = adata_golden.copy()

        # 2. Run Codon step on the isolated copy
        t0 = time.perf_counter()
        try:
            step_fn_cd(adata_codon_test)
            codon_timings[label] = time.perf_counter() - t0
            cd_snap = snapshot(adata_codon_test, label)
            cd_success = True
        except Exception as exc:
            codon_timings[label] = None
            print(f"  [ERROR] Codon {label} raised: {exc}")
            cd_success = False

        # 3. Run Python step on the golden object (advancing the pipeline truth)
        t0 = time.perf_counter()
        try:
            step_fn_py(adata_golden)
            python_timings[label] = time.perf_counter() - t0
            py_snap = snapshot(adata_golden, label)
            py_success = True
        except Exception as exc:
            python_timings[label] = None
            print(f"  [ERROR] Python {label} raised: {exc}")
            py_success = False

        if label == "neighbors":
            # After running both pipelines, compare:
            print("scanpy n_neighbors stored:", adata_golden.uns['neighbors']['params']['n_neighbors'])
            print("codon distances shape:", adata_codon_test.obsp['distances'].shape)
            print("codon distances nnz per row:", adata_codon_test.obsp['distances'].getnnz(axis=1)[:10])
            print("scanpy distances nnz per row:", adata_golden.obsp['distances'].getnnz(axis=1)[:10])

            # Check what X_pca looks like going in
            print("scanpy X_pca shape:", adata_golden.obsm['X_pca'].shape)
            print("codon X_pca shape:", adata_codon_test.obsm['X_pca'].shape)

            print("Scanpy row 0:", adata_golden.obsp['distances'][0])
            print("Codon  row 0:", adata_codon_test.obsp['distances'][0])

            print("Scanpy row 1:", adata_golden.obsp['distances'][1])
            print("Codon  row 1:", adata_codon_test.obsp['distances'][1])

        # 4. Compare the results (injecting X_ref for manifold checks)
        if py_success and cd_success:
            # Extract the true PCA matrix if it exists in the golden object yet
            X_ref = None
            if "X_pca" in adata_golden.obsm:
                X_ref = adata_golden.obsm["X_pca"]

            deviations = compare_snapshots(label, py_snap, cd_snap, X_ref=X_ref)
            deviations_log[label] = deviations
        else:
            deviations_log[label] = None

    # ------------------------------------------------------------------
    # Correctness report
    # ------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("ISOLATED CORRECTNESS REPORT (TIERED)")
    print("=" * 80)

    any_deviation = False
    for label, _ in python_steps:
        py_failed = python_timings.get(label) is None
        cd_failed = codon_timings.get(label) is None
        devs = deviations_log.get(label)

        if py_failed and cd_failed:
            print(f"  [SKIP ] {label:<25}  both versions failed")
            continue
        if py_failed:
            print(f"  [SKIP ] {label:<25}  Python failed — no reference")
            continue
        if cd_failed:
            print(f"  [FAIL ] {label:<25}  Codon raised an exception")
            any_deviation = True
            continue

        if devs:
            any_deviation = True
            print(f"  [FAIL ] {label:<25}  {len(devs)} deviation(s):")
            for d in devs:
                print(f"            • {d}")
        else:
            print(f"  [OK   ] {label:<25}  outputs verified")

    if not any_deviation:
        print("\n  All steps produced verified isolated outputs.")

    # ------------------------------------------------------------------
    # Timing report
    # ------------------------------------------------------------------
    python_total = sum(t for t in python_timings.values() if t is not None)
    codon_total  = sum(t for t in codon_timings.values()  if t is not None)

    print("\n" + "=" * 80)
    print("BENCHMARK RESULTS (ISOLATED EXECUTION)")
    print("=" * 80)
    print(f"{'Step':<25} {'Python (s)':<14} {'Codon (s)':<14} {'Speedup':<10}")
    print("-" * 63)
    for label, _ in python_steps:
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




def umap_noise_baseline_3k_PBMCs():
    """Debug-only run: UMAP disparity from RNG-driven negative sampling on PBMC3k."""
    sp, _ = setup_imports()

    adata = sp.datasets.pbmc3k()
    adata.var_names_make_unique()

    sp.pp.calculate_qc_metrics(adata)
    sp.pp.filter_cells(adata, min_genes=100)
    sp.pp.filter_genes(adata, min_cells=3)
    sp.pp.normalize_total(adata)
    sp.pp.log1p(adata)
    sp.pp.highly_variable_genes(adata, n_top_genes=2000)
    sp.pp.scale(adata, max_value=10)
    sp.tl.pca(adata)
    sp.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

    adata0 = adata.copy()
    adata1 = adata.copy()

    sp.tl.umap(adata0, random_state=0)
    sp.tl.umap(adata1, random_state=1)

    X0 = np.asarray(adata0.obsm["X_umap"])
    X1 = np.asarray(adata1.obsm["X_umap"])

    from scipy.spatial import procrustes
    _, _, disparity = procrustes(X0, X1)

    print("[INFO] UMAP noise baseline (random_state 0 vs 1) Procrustes disparity:", disparity)
    return disparity

def umap_native_noise_baseline_3k_PBMCs():
    """Debug-only run: native UMAP disparity from RNG-driven negative sampling, seed 0 vs seed 1."""
    sp, sc = setup_imports()

    adata = sp.datasets.pbmc3k()
    adata.var_names_make_unique()

    sp.pp.calculate_qc_metrics(adata)
    sp.pp.filter_cells(adata, min_genes=100)
    sp.pp.filter_genes(adata, min_cells=3)
    sp.pp.normalize_total(adata)
    sp.pp.log1p(adata)
    sp.pp.highly_variable_genes(adata, n_top_genes=2000)
    sp.pp.scale(adata, max_value=10)
    sp.tl.pca(adata)
    sp.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

    adata0 = adata.copy()
    adata1 = adata.copy()

    sc.tl.umap(adata0, random_state=0)
    sc.tl.umap(adata1, random_state=1)

    X0 = np.asarray(adata0.obsm["X_umap"])
    X1 = np.asarray(adata1.obsm["X_umap"])

    from scipy.spatial import procrustes
    _, _, disparity = procrustes(X0, X1)

    print("[INFO] NATIVE UMAP noise baseline (random_state 0 vs 1) Procrustes disparity:", disparity)
    return disparity

if __name__ == "__main__":
    #umap_noise_baseline_3k_PBMCs()
    #umap_native_noise_baseline_3k_PBMCs()
    correctness_benchmark_3k_PBMCs()