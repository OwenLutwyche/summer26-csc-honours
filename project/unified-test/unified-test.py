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

# resolve paths dynamically
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

SCANPY_PATH = os.path.join(PROJECT_ROOT, "scanpy-main")
SCANCODON_PATH = os.path.join(PROJECT_ROOT, "scancodon")

def setup_imports():
    """
    Set up both scanpy (Python) and scancodon (Codon) imports separately.

    Returns: Tuple(scanpy_module, scancodon_module).
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
    """
    Dynamically load a Python module from a file path.
    used for executing indpendent unit tests without standard imports
    """
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
    Evaluate test output against a reference output using Strict Determinism (Arithmetic exact matches)
    """
    tol_str = f"strict(rtol={rtol}, atol={atol})"
    if pv.shape != cv.shape:
        return f"Shape mismatch: Python {pv.shape} vs Codon {cv.shape}", tol_str
    
    pf = pv.astype(np.float32)
    cf = cv.astype(np.float32)
    
    p_nans = np.isnan(pf)
    c_nans = np.isnan(cf)
    
    if not np.array_equal(p_nans, c_nans):
        return "NaN alignment mismatch: NaNs occur at different indices between Python and Codon.", tol_str
        
    valid_mask = ~p_nans
    
    if not np.any(valid_mask):
        return None, tol_str
        
    pf_valid = pf[valid_mask]
    cf_valid = cf[valid_mask]

    is_close = np.isclose(pf_valid, cf_valid, rtol=rtol, atol=atol)
    if not np.all(is_close):
        mismatch_idx = np.where(~is_close)[0]
        p_dev = pf_valid[mismatch_idx]
        c_dev = cf_valid[mismatch_idx]
        
        diffs = np.abs(p_dev - c_dev)
        max_diff = np.max(diffs)
        
        if debug_print:
            print(f"\n[EVALUATE_STRICT] evaluate_strict failed. Found {len(mismatch_idx)} deviating elements.")
            print(f"{'Valid Array Index':<20} | {'Python Value':<22} | {'Codon Value':<22} | {'Absolute Diff'}")
            print("-" * 85)
            
            sorted_args = np.argsort(diffs)[::-1]
            display_idx = sorted_args[:10]
            
            for i in display_idx:
                orig_i = mismatch_idx[i]
                print(f"{orig_i:<20} | {p_dev[i]:<22.8e} | {c_dev[i]:<22.8e} | {diffs[i]:.8e}")
                
            if len(mismatch_idx) > 10:
                print(f"... and {len(mismatch_idx) - 10} more mismatches hidden.")
                
        return f"Numeric deviation exceeds tolerance (max |diff| = {max_diff:.3e})", tol_str
        
    return None, tol_str


def evaluate_linear_subspace(pv, cv, max_disparity=1e-2):
    """
    Linear Subspaces (PCA / Diffmap structural alignment)
    """
    tol_str = f"procrustes(max_disp={max_disparity})"
    if pv.shape != cv.shape:
        return f"Shape mismatch: Python {pv.shape} vs Codon {cv.shape}", tol_str
    
    from scipy.spatial import procrustes
    try:
        _, _, disparity = procrustes(pv, cv)
        if disparity > max_disparity:
            return f"Procrustes disparity {disparity:.4e} exceeds tolerance threshold {max_disparity:.4e}", tol_str
    except Exception as e:
        return f"Procrustes calculation failed: {e}", tol_str
    return None, tol_str


def evaluate_graph_topology(pv, cv, min_jaccard=0.85):
    """
    Evaluate Graph Topology parity (neighbors graph overlap)
    """
    tol_str = f"jaccard(min_overlap={min_jaccard})"
    if pv.shape != cv.shape:
        return f"Shape mismatch: Python {pv.shape} vs Codon {cv.shape}", tol_str
    
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
        return f"Average row-wise Jaccard similarity {avg_jaccard:.4f} is below threshold {min_jaccard:.4f}", tol_str
    return None, tol_str


def evaluate_clustering(pv, cv, min_ari=0.90):
    """
    clustering parity (leiden label assignments)
    """
    tol_str = f"ARI(min_score={min_ari})"
    import numpy as np
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    
    n_py = len(np.unique(pv))
    n_cd = len(np.unique(cv))
    
    if abs(n_py - n_cd) > 1:
        return f"Cluster count mismatch: Python={n_py}, Codon={n_cd}", tol_str
        
    ari = adjusted_rand_score(pv, cv)
    
    if ari < min_ari:
        nmi = normalized_mutual_info_score(pv, cv)
        return f"Structural divergence: ARI={ari:.4f} (threshold {min_ari:.4f}), NMI={nmi:.4f} | Clusters: Py={n_py}, Codon={n_cd}", tol_str
        
    return None, tol_str


def evaluate_stochastic_manifold(pv, cv, X_ref=None, min_trustworthiness=0.80, max_disparity=0.5):
    """
    Stochastic Manifolds (UMAP / T-SNE neighborhood conservation)
    """
    tol_str = f"manifold(trust>={min_trustworthiness}, disp<={max_disparity})"
    if pv.shape != cv.shape:
        return f"Shape mismatch: Python {pv.shape} vs Codon {cv.shape}", tol_str
        
    errors = []
    
    from scipy.spatial import procrustes
    try:
        _, _, disparity = procrustes(pv, cv)
        if disparity > max_disparity:
            errors.append(f"Global Procrustes disparity {disparity:.8f} exceeds threshold {max_disparity:.4f}")
    except Exception as e:
        errors.append(f"Procrustes failed: {e}")
        
    if X_ref is not None:
        from sklearn.manifold import trustworthiness
        try:
            t_score = trustworthiness(X_ref, cv, n_neighbors=15)
            if t_score < min_trustworthiness:
                errors.append(f"Local Trustworthiness score {t_score:.4f} below threshold {min_trustworthiness:.4f}")
        except Exception as e:
            errors.append(f"Trustworthiness benchmark failed: {e}")
            
    if errors:
        return " | ".join(errors), tol_str
    return None, tol_str


def debug_and_evaluate_rank_genes_groups(pv, cv, k_prefix="uns.rank_genes_groups"):
    '''
    Closely assess differences in rank_genes_groups output
    '''
    import numpy as np
    deviations = []
    tol_str = "rank_genes(rtol=1e-1/1e-2, atol=1e-3/1e-4)"
    
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
            
            if metric == "names":
                p_set, c_set = set(p_vals), set(c_vals)
                jaccard = len(p_set & c_set) / len(p_set | c_set) if p_set else 1.0
                if jaccard < 0.99:
                    deviations.append(f"{k_prefix}['names']['{group_id}']: Gene set composition mismatch (Jaccard={jaccard:.3f})")
                continue
            
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
            
            rtol = 1e-1 if metric in ("pvals", "pvals_adj") else 1e-2
            atol = 1e-3 if metric in ("pvals", "pvals_adj") else 1e-4
            
            is_close = np.isclose(pf_valid, cf_valid, rtol=rtol, atol=atol)
            
            if not np.all(is_close):
                mismatches = np.where(~is_close)[0]
                diffs = np.abs(pf_valid - cf_valid)
                max_diff = np.max(diffs)
                
                deviations.append(f"{k_prefix}['{metric}']['{group_id}']: {len(mismatches)} deviations (max |diff| = {max_diff:.3e})")
                
    return deviations, tol_str

def run_isolated_correctness_benchmark(adata_loader, benchmark_label="benchmark"):
    """
    Execute a full single-cell pipeline correctness benchmark.

    adata_loader: callable(sp) -> AnnData. Receives the Python scanpy module so it
                  can use sp.datasets.*, sp.read_10x_h5, etc. to build/load the input.
    benchmark_label: short human-readable name, used only in log output.
    """
    import numpy as np

    print("=" * 80)
    print(f"--- CORRECTNESS BENCHMARK: ISOLATED STEPS (TIERED) [{benchmark_label}] ---")
    print("=" * 80)

    sp, sc = setup_imports()

    print(f"[INFO] Loading dataset for '{benchmark_label}'...")
    t_load = time.perf_counter()
    adata = adata_loader(sp)
    print(f"[INFO] Loaded AnnData shape={adata.shape} in {time.perf_counter() - t_load:.1f}s")
    adata.var_names_make_unique()

    def get_steps(lib):
        return [
            ("calculate_qc_metrics",  lambda a: lib.pp.calculate_qc_metrics(a)),
            ("filter_cells",          lambda a: lib.pp.filter_cells(a, min_genes=100)),
            ("filter_genes",          lambda a: lib.pp.filter_genes(a, min_cells=3)),
            #("scrublet",              lambda a: lib.pp.scrublet(a, random_state=0)),
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
        tolerances = set()
        
        if py_snap["shape"] != cd_snap["shape"]:
            deviations.append(f"shape mismatch: Python={py_snap['shape']}  Codon={cd_snap['shape']}")
            tolerances.add("exact dimensions")

        all_keys = set(py_snap) | set(cd_snap)
        for k in sorted(all_keys):
            if k == "shape": continue
            
            if k not in py_snap: deviations.append(f"{k}: present in Codon but missing in Python"); continue
            if k not in cd_snap: deviations.append(f"{k}: present in Python but missing in Codon"); continue

            pv, cv = py_snap[k], cd_snap[k]
            
            if isinstance(pv, str) and pv.startswith("<missing"):
                deviations.append(f"{k}: Python could not read value ({pv})"); continue
            if isinstance(cv, str) and cv.startswith("<missing"):
                deviations.append(f"{k}: Codon could not read value ({cv})"); continue

            err = None
            tol = None

            if k == "uns.rank_genes_groups":
                devs, tol = debug_and_evaluate_rank_genes_groups(pv, cv, k)
                deviations.extend(devs)
                tolerances.add(tol)
                continue

            if "X_umap" in k or "X_tsne" in k:
                err, tol = evaluate_stochastic_manifold(pv, cv, X_ref=X_ref)
                if err: deviations.append(f"{k}: {err}")
                tolerances.add(tol)

            elif "X_pca" in k or "PCs" in k or "X_diffmap" in k:
                err, tol = evaluate_linear_subspace(pv, cv)
                if err: deviations.append(f"{k}: {err}")
                tolerances.add(tol)

            elif "connectivities" in k or "distances" in k:
                err, tol = evaluate_graph_topology(pv, cv)
                if err: deviations.append(f"{k}: {err}")
                tolerances.add(tol)

            elif "leiden" in k:
                err, tol = evaluate_clustering(pv, cv)
                if err: deviations.append(f"{k}: {err}")
                tolerances.add(tol)

            elif isinstance(pv, np.ndarray) and isinstance(cv, np.ndarray):
                if pv.dtype.kind in ("U", "O"):
                    tolerances.add("exact string match")
                    if not np.array_equal(pv, cv):
                        n_diff = int(np.sum(pv != cv))
                        deviations.append(f"{k}: {n_diff}/{len(pv)} label mappings differ")
                else:
                    err, tol = evaluate_strict(pv, cv)
                    if err: deviations.append(f"{k}: {err}")
                    tolerances.add(tol)
                    
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
                                        sub_err, tol = evaluate_strict(pf_aligned[mask], cf_aligned[mask], rtol=1e-1, atol=1e-3)
                                    else:
                                        sub_err, tol = evaluate_strict(pf_aligned[mask], cf_aligned[mask], rtol=1e-2, atol=1e-4)
                                        
                                    if sub_err: 
                                        deviations.append(f"{k}.{subk}['{group_id}']: {sub_err}")
                                    tolerances.add(tol)
                        elif isinstance(p_sub, np.ndarray) and isinstance(c_sub, np.ndarray):
                            if p_sub.dtype.kind not in ("U", "O"):
                                sub_err, tol = evaluate_strict(p_sub, c_sub, rtol=1e-3, atol=1e-5)
                                if sub_err: deviations.append(f"{k}.{subk}: {sub_err}")
                                tolerances.add(tol)
                
                if deviations:
                    print(f"\n[DEBUG] Rank Genes Groups Mismatch detected. Probing Group '0':")
                    group_id = '0'
                    try:
                        py_names = pv['names'][group_id]
                        top_gene = py_names[0]
                        p_idx = 0
                        
                        co_names = cv['names'][group_id]
                        c_idx = np.where(co_names == top_gene)[0][0]
                        
                        print(f"  Target Gene: '{top_gene}'")
                        print(f"  Python Rank: {p_idx} | Codon Rank: {c_idx}")
                        
                        py_score = pv['scores'][group_id][p_idx]
                        co_score = cv['scores'][group_id][c_idx]
                        print(f"  Scores: Py={py_score:.4f}, Co={co_score:.4f}")
                        
                        py_lfc = pv['logfoldchanges'][group_id][p_idx]
                        co_lfc = cv['logfoldchanges'][group_id][c_idx]
                        print(f"  LFC:    Py={py_lfc:.4f}, Co={co_lfc:.4f}")
                        
                        py_pv = pv['pvals'][group_id][p_idx]
                        co_pv = cv['pvals'][group_id][c_idx]
                        print(f"  P-vals: Py={py_pv:.4e}, Co={co_pv:.4e}")
                        
                    except Exception as e:
                        print(f"  [!] Debug probe failed: {e}")
            
            else:
                tolerances.add("exact match")
                try:
                    if pv != cv: deviations.append(f"{k}: values differ Python={pv!r} Codon={cv!r}")
                except Exception:
                    pass
            
        return deviations, list(tolerances)

    # ==========================================
    # EXECUTION PIPELINE
    # ==========================================
    adata_golden = adata.copy()
    python_timings = {}
    codon_timings = {}
    deviations_log = {}
    tolerance_log = {}
    
    python_steps = get_steps(sp)
    codon_steps_dict = dict(get_steps(sc))

    print("[INFO] Running isolated correctness benchmark...")

    for label, step_fn_py in python_steps:
        step_fn_cd = codon_steps_dict[label]

        adata_codon_test = adata_golden.copy()

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

        if py_success and cd_success:
            X_ref = None
            if "X_pca" in adata_golden.obsm:
                X_ref = adata_golden.obsm["X_pca"]

            deviations, tols = compare_snapshots(label, py_snap, cd_snap, X_ref=X_ref)
            deviations_log[label] = deviations
            tolerance_log[label] = ", ".join(tols) if tols else "exact match"
        else:
            deviations_log[label] = None
            tolerance_log[label] = None

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
        tol = tolerance_log.get(label)

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
            print(f"  [OK   ] {label:<25}  outputs verified within tolerance: {tol}")

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
    


def correctness_benchmark_3k_PBMCs():
    """
    Isolated-step correctness/perf benchmark on the 3k PBMC dataset (small, fast baseline).
    """
    run_isolated_correctness_benchmark(lambda sp: sp.datasets.pbmc3k(), "3k PBMCs")


def run_tests_for_library(test_files, lib_name, lib_module):
    """
    Run all test files for a given library.
    Dynamically load and execute test_*.py files containing run_all() evaluation function
    Track success/failure
    """
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
    """
    Print a side-by-side comparison of results.
    Generate a table comparing unit test pass/failure rates
    """
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
    
    python_times = {name: t for name, t in python_results["per_file_timings"]}
    codon_times = {name: t for name, t in codon_results["per_file_timings"]}
    
    all_files = set(python_times.keys()) | set(codon_times.keys())
    for test_name in sorted(all_files):
        python_t = python_times.get(test_name, 0)
        codon_t = codon_times.get(test_name, 0)
        print(f"{test_name:<35} {python_t:<15.2f} {codon_t:<15.2f}")



def create_comprehensive_edge_case_adata():
    import pandas as pd
    import anndata as ad
    import numpy as np
    
    genes = [
        "Negative_Expr",     # 1. Scaled data simulation
        "Negative_ZeroVar",  # 2. Flatline negative data
        "Mixed_NaN_A",       # 3. Only ONE cell is NaN
        "Mixed_Inf_A",       # 4. Only ONE cell is Inf 
        "All_Inf",           # 5. All cells are Inf
        "Tiny_Negative",     # 6. Means right on -1e-9 boundary
        "Massive_T",         # 7. Huge effect size -> pushes t-stat high, risking pval underflow to 0.0
        "Zero_Diff",         # 8. Identical means & variance -> forces pval to 1.0
        "Extreme_Variance",  # 9. Wildly different group variances (stresses Welch-Satterthwaite DF)
        "Normal_Reference"   # 10. Baseline sanity check
    ]
    
    X = np.array([
        # NegExpr | NegZV | MixNaN | MixInf | AllInf | TinyNeg | MassiveT | ZeroDiff | ExtVar | Normal
        [ -2.5,     -5.0,   1.0,     1.0,     np.inf,  -1e-8,    100.0,     2.5,       0.1,     2.5 ], # A
        [ -2.3,     -5.0,   np.nan,  np.inf,  np.inf,  -1e-8,    100.0,     2.5,       10.0,    2.7 ], # A
        [ -2.6,     -5.0,   1.5,     1.5,     np.inf,  -1e-8,    100.0,     2.5,       0.1,     2.4 ], # A
        [ -2.1,     -5.0,   1.2,     1.2,     np.inf,  -1e-8,    100.0,     2.5,       10.0,    2.6 ], # A
        
        [ 0.5,      0.0,    1.0,     1.0,     np.inf,  0.0,      0.0,       2.5,       5.0,     1.2 ], # B
        [ 0.6,      0.0,    1.1,     1.1,     np.inf,  0.0,      0.0,       2.5,       5.0,     1.1 ], # B
        [ 0.4,      0.0,    1.0,     1.0,     np.inf,  0.0,      0.0,       2.5,       5.0,     1.3 ], # B
        [ 0.5,      0.0,    1.2,     1.2,     np.inf,  0.0,      0.0,       2.5,       5.0,     1.2 ], # B
        
        [ -1.0,     -2.0,   1.5,     1.5,     np.inf,  -1e-8,    50.0,      2.5,       5.0,     1.5 ], # C
        [ -1.1,     -2.0,   1.4,     1.4,     np.inf,  -1e-8,    50.0,      2.5,       5.0,     1.4 ], # C 
    ], dtype=np.float32)
    
    groups = pd.Categorical(
        ['A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'C', 'C'], 
        categories=['A', 'B', 'C']
    )

    additional_genes = [
        "Large_Offset_TinyVar",
        "Subnormal_Tiny",        # 12. Pushes past Tiny_Negative's scale (~1e-30) --
                                #     re-probes df_den underflowing to literal 0.0
        "DF_Boundary_Case",      # 13. Same shape as Tiny_Negative but at 1e-4 scale --
                                #     regression lock for the df_den fix
        "Tiny_Positive",     
    ]

    X_additional = np.array([
        # LargeOffset | Subnormal | DFBoundary | TinyPos
        [ 500.001,      -1e-30,     -1e-4,       1e-8 ], # A
        [ 500.002,      -1e-30,     -1e-4,       1e-8 ], # A
        [ 499.999,      -1e-30,     -1e-4,       1e-8 ], # A
        [ 500.000,      -1e-30,     -1e-4,       1e-8 ], # A

        [ 500.501,        0.0,       0.0,         0.0  ], # B
        [ 500.502,        0.0,       0.0,         0.0  ], # B
        [ 500.499,        0.0,       0.0,         0.0  ], # B
        [ 500.500,        0.0,       0.0,         0.0  ], # B

        [ 500.301,      -1e-30,     -1e-4,       1e-8 ], # C
        [ 500.300,      -1e-30,     -1e-4,       1e-8 ], # C
    ], dtype=np.float32)

    genes = genes + additional_genes
    X = np.hstack([X, X_additional])
    obs = pd.DataFrame({'group': groups})
    var = pd.DataFrame(index=genes)
    
    return ad.AnnData(X=X, obs=obs, var=var)

def evaluate_RGG_edge_cases():
    sp, sc = setup_imports()
    
    print("\n" + "=" * 80)
    print("--- DIFFERENTIAL EXPRESSION EDGE-CASE DIAGNOSTIC ---")
    print("=" * 80)
    
    adata_sp = create_comprehensive_edge_case_adata()
    adata_sc = adata_sp.copy()
    
    print("[INFO] Running Scanpy reference...")
    sp.tl.rank_genes_groups(adata_sp, groupby='group', method='t-test')
    
    print("[INFO] Running Scancodon implementation...")
    try:
        sc.tl.rank_genes_groups(adata_sc, groupby='group', method='t-test')
        cd_success = True
    except Exception as e:
        print(f"[ERROR] Scancodon rank_genes_groups failed: {e}")
        cd_success = False
        
    if not cd_success:
        return
        
    res_sp = adata_sp.uns['rank_genes_groups']
    res_sc = adata_sc.uns['rank_genes_groups']
    
    group_id = 'A'
    
    sp_names = res_sp['names'][group_id]
    sc_names = res_sc['names'][group_id]
    
    sp_scores = res_sp['scores'][group_id]
    sc_scores = res_sc['scores'][group_id]
    
    sp_pvals = res_sp['pvals'][group_id]
    sc_pvals = res_sc['pvals'][group_id]
    
    sp_lfc = res_sp['logfoldchanges'][group_id]
    sc_lfc = res_sc['logfoldchanges'][group_id]
    
    sp_dict = {n: (s, p, l) for n, s, p, l in zip(sp_names, sp_scores, sp_pvals, sp_lfc)}
    sc_dict = {n: (s, p, l) for n, s, p, l in zip(sc_names, sc_scores, sc_pvals, sc_lfc)}
    
    print(f"\nComparing Group '{group_id}' vs Rest:")
    print(f"{'Gene':<18} | {'Py Score':>10} | {'Co Score':>10} | {'Py P-val':>10} | {'Co P-val':>10} | {'Py LFC':>10} | {'Co LFC':>10}")
    print("-" * 105)
    
    genes = adata_sp.var_names
    for gene in genes:
        if gene not in sp_dict or gene not in sc_dict:
            print(f"{gene:<18} | Missing in one or both outputs")
            continue
            
        s_py, p_py, l_py = sp_dict[gene]
        s_co, p_co, l_co = sc_dict[gene]
        
        print(f"{gene:<18} | {s_py:>10.4f} | {s_co:>10.4f} | {p_py:>10.2e} | {p_co:>10.2e} | {l_py:>10.4f} | {l_co:>10.4f}")
    
    print("-" * 105)

if __name__ == "__main__":

    correctness_benchmark_3k_PBMCs()

    # Heavy: downloads ~4GB and may consume tens of GB of RAM / crash. Uncomment to run.
    #correctness_benchmark_1M_neurons()

    #evaluate_RGG_edge_cases()