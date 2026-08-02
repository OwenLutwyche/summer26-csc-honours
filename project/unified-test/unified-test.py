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
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

SCANPY_PATH = os.path.join(PROJECT_ROOT, "scanpy-main")
SCANCODON_PATH = os.path.join(PROJECT_ROOT, "scancodon")

def setup_imports():
    """
    Set up both scanpy (Python) and scancodon (Codon) imports separately.
    Returns a tuple of (scanpy_module, scancodon_module).
    """
    # python scanpy
    scanpy_src = os.path.join(SCANPY_PATH, "src")
    
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
    original_cwd = os.getcwd()
    
    try:
        os.chdir(SCANCODON_PATH)
        if SCANCODON_PATH not in sys.path:
            sys.path.insert(0, SCANCODON_PATH)
        
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

def evaluate_strict(pv, cv, rtol=1e-4, atol=1e-6, debug_print = False):
    """
    Strict Determinism (Arithmetic exact matches)
    applies standard close tolerances for deterministic matrix operations, just a straight np.allclose
    should work for log1p, normalize_total, qc_metrics, filter_genes/cells, scale, hvg and scrublet (if it's seeded deterministic-style)
    """
    if pv.shape != cv.shape:
        return f"Shape mismatch: Python {pv.shape} vs Codon {cv.shape}"
    
    # Cast down to float32 to strip Codon's 64-bit precision tail and enable safe NaN checking
    pf = pv.astype(np.float32)
    cf = cv.astype(np.float32)
    
    # 1. Verify NaNs appear at the exact same indices in both arrays
    p_nans = np.isnan(pf)
    c_nans = np.isnan(cf)
    
    if not np.array_equal(p_nans, c_nans):
        return "NaN alignment mismatch: NaNs occur at different indices between Python and Codon."
        
    # 2. Create a mask of only the valid, non-NaN numbers
    valid_mask = ~p_nans
    
    # If the array is entirely NaNs, it's a perfect structural match
    if not np.any(valid_mask):
        return None
        
    # 3. Apply the mask to evaluate only the valid numeric values
    pf_valid = pf[valid_mask]
    cf_valid = cf[valid_mask]

    is_close = np.isclose(pf_valid, cf_valid, rtol=rtol, atol=atol)
    # 4. Compare the remaining valid numbers
    if not np.all(is_close):
        # Isolate the exact elements that failed
        mismatch_idx = np.where(~is_close)[0]
        p_dev = pf_valid[mismatch_idx]
        c_dev = cf_valid[mismatch_idx]
        
        # Calculate absolute differences
        diffs = np.abs(p_dev - c_dev)
        max_diff = np.max(diffs)
        
        if debug_print:
            print(f"\n[EVALUATE_STRICT] evaluate_strict failed. Found {len(mismatch_idx)} deviating elements.")
            print(f"{'Valid Array Index':<20} | {'Python Value':<22} | {'Codon Value':<22} | {'Absolute Diff'}")
            print("-" * 85)
            
            # Sort by largest difference to surface the worst offenders immediately
            sorted_args = np.argsort(diffs)[::-1]
            display_idx = sorted_args[:10]  # Show the top 10 worst deviations
            
            for i in display_idx:
                orig_i = mismatch_idx[i]
                print(f"{orig_i:<20} | {p_dev[i]:<22.8e} | {c_dev[i]:<22.8e} | {diffs[i]:.8e}")
                
            if len(mismatch_idx) > 10:
                print(f"... and {len(mismatch_idx) - 10} more mismatches hidden.")
                
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


def debug_and_evaluate_rank_genes_groups(pv, cv, k_prefix="uns.rank_genes_groups"):
    import numpy as np
    deviations = []
    
    metrics = ["names", "scores", "logfoldchanges", "pvals", "pvals_adj"]
    for metric in metrics:
        if metric not in pv or metric not in cv:
            deviations.append(f"{k_prefix}: Missing metric '{metric}'")
            continue
            
        p_rec = pv[metric]
        c_rec = cv[metric]
        
        if not (hasattr(p_rec, 'dtype') and p_rec.dtype.names):
            continue
            
        for group_id in p_rec.dtype.names:
            if group_id not in c_rec.dtype.names:
                continue
                
            p_vals = p_rec[group_id]
            c_vals = c_rec[group_id]
            
            # 1. Soft Evaluate Name Ordering
            if metric == "names":
                p_set, c_set = set(p_vals), set(c_vals)
                jaccard = len(p_set & c_set) / len(p_set | c_set) if p_set else 1.0
                if jaccard < 0.99:
                    deviations.append(f"{k_prefix}['names']['{group_id}']: Gene set composition mismatch (Jaccard={jaccard:.3f})")
                continue
            
            # 2. Dictionary Alignment
            p_names = pv["names"][group_id]
            c_names = cv["names"][group_id]
            
            p_map = dict(zip(p_names, p_vals))
            c_map = dict(zip(c_names, c_vals))
            
            common_genes = sorted(list(set(p_map.keys()) & set(c_map.keys())))
            if not common_genes:
                continue
                
            pf_aligned = np.array([p_map[g] for g in common_genes], dtype=np.float32)
            cf_aligned = np.array([c_map[g] for g in common_genes], dtype=np.float32)
            names_aligned = np.array(common_genes)
            
            # 3. Masking (Handle NaNs and Zeros)
            p_nans = np.isnan(pf_aligned)
            c_nans = np.isnan(cf_aligned)
            
            if not np.array_equal(p_nans, c_nans):
                deviations.append(f"{k_prefix}['{metric}']['{group_id}']: NaN alignment mismatch")
                continue
            
            valid_mask = ~p_nans
            
            if metric not in ("pvals", "pvals_adj"):
                valid_mask = valid_mask & (np.abs(pf_aligned) > 1e-4)
                
            if not np.any(valid_mask):
                continue
                
            pf_valid = pf_aligned[valid_mask]
            cf_valid = cf_aligned[valid_mask]
            
            # 4. Math Evaluation
            rtol = 1e-1 if metric in ("pvals", "pvals_adj") else 1e-2
            atol = 1e-3 if metric in ("pvals", "pvals_adj") else 1e-4
            
            is_close = np.isclose(pf_valid, cf_valid, rtol=rtol, atol=atol)
            
            if not np.all(is_close):
                mismatches = np.where(~is_close)[0]
                diffs = np.abs(pf_valid - cf_valid)
                max_diff = np.max(diffs)
                
                deviations.append(f"{k_prefix}['{metric}']['{group_id}']: {len(mismatches)} deviations (max |diff| = {max_diff:.3e})")
                
    return deviations

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
            
            if k == "uns.rank_genes_groups":
                devs = debug_and_evaluate_rank_genes_groups(pv, cv, k)
                deviations.extend(devs)
                continue


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
                    
            # Structured dictionaries (like rank_genes_groups)
            elif isinstance(pv, dict) and isinstance(cv, dict) and k != "uns.rank_genes_groups":
                py_keys, cd_keys = set(pv.keys()), set(cv.keys())
                if py_keys != cd_keys:
                    deviations.append(f"{k}: dictionary keys differ Python-only={py_keys - cd_keys} Codon-only={cd_keys - py_keys}")
                else:
                    for subk in sorted(py_keys):
                        if subk == 'names':
                            continue
                            
                        p_sub, c_sub = pv[subk], cv[subk]
                        if hasattr(p_sub, 'dtype') and hasattr(c_sub, 'dtype') and p_sub.dtype.names is not None:
                            for group_id in p_sub.dtype.names:
                                if group_id in c_sub.dtype.names:
                                    p_vals = p_sub[group_id]
                                    c_vals = c_sub[group_id]
                                    
                                    p_gene_names = pv['names'][group_id]
                                    c_gene_names = cv['names'][group_id]
                                    
                                    p_map = dict(zip(p_gene_names, p_vals))
                                    c_map = dict(zip(c_gene_names, c_vals))
                                    
                                    common_genes = sorted(list(set(p_map.keys()) & set(c_map.keys())))
                                    
                                    pf_aligned = np.array([p_map[g] for g in common_genes], dtype=float)
                                    cf_aligned = np.array([c_map[g] for g in common_genes], dtype=float)
                                    
                                    mask = np.abs(pf_aligned) > 1e-4
                                    
                                    if subk in ("pvals", "pvals_adj"):
                                        sub_err = evaluate_strict(pf_aligned[mask], cf_aligned[mask], rtol=1e-1, atol=1e-3)
                                    else:
                                        sub_err = evaluate_strict(pf_aligned[mask], cf_aligned[mask], rtol=1e-2, atol=1e-4)
                                        
                                    if sub_err: 
                                        deviations.append(f"{k}.{subk}['{group_id}']: {sub_err}")
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

        # 4. Compare the results (injecting X_ref for manifold checks)
        if py_success and cd_success:
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
    correctness_benchmark_3k_PBMCs()