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
import matplotlib.pyplot as plt

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

        applies standard close tolerances for deterministic matrix operations using np.allclose
    suitable for deterministic functions such as 
        - log1p, 
        - normalize_total, 
        - qc_metrics, 
        - filter_genes/cells, 
        - scale, 
        - highly variable genes 
        - and if seeded deterministically: scrublet
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

    Uses the Adjusted Rand Index (ARI) to determine if partition boundaries are functionally identical, independent of label permutation changes
    Suitable for Leiden clustering
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

        Combines Procrustes Analysis for global cluster macro-structures and sklearn.manifold.trustworthiness index 
    to score local neighborhood integrity against high-dimensional PCA space. 
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

    - Align unordered gene records by name
        - Check jaccard similarity
        - Gracefully mask NaNs
        - Apply varying tolerances across scores, log-fold changes, and p-values
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
        '''
        Extract output matrices and dicts corresponding to the given step_label
        '''
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
        '''
        Pass extracted snapshots to their appropriate evaluation function
        '''
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

            #  Stochastic Manifolds
            if "X_umap" in k or "X_tsne" in k:
                err, tol = evaluate_stochastic_manifold(pv, cv, X_ref=X_ref)
                if err: deviations.append(f"{k}: {err}")
                tolerances.add(tol)

            # Linear Subspaces
            elif "X_pca" in k or "PCs" in k or "X_diffmap" in k:
                err, tol = evaluate_linear_subspace(pv, cv)
                if err: deviations.append(f"{k}: {err}")
                tolerances.add(tol)

            # Graph Topology
            elif "connectivities" in k or "distances" in k:
                err, tol = evaluate_graph_topology(pv, cv)
                if err: deviations.append(f"{k}: {err}")
                tolerances.add(tol)

            # Clustering Groups
            elif "leiden" in k:
                err, tol = evaluate_clustering(pv, cv)
                if err: deviations.append(f"{k}: {err}")
                tolerances.add(tol)

            # Structured Array fallbacks
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
    py_snapshots = {}   # step_label -> snapshot dict, kept for post-hoc dot plots
    cd_snapshots = {}
    
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

        if py_success:
            py_snapshots[label] = py_snap
        if cd_success:
            cd_snapshots[label] = cd_snap

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
    # ------------------------------------------------------------------
    # Visual Comparison (From Isolated Snapshots)
    # ------------------------------------------------------------------
    import pandas as pd
    import anndata as ad
    from contextlib import contextmanager

    print("\n" + "=" * 80)
    print("VISUAL INSPECTION: SIDE-BY-SIDE DOT PLOTS")
    print("=" * 80)

    safe_label = benchmark_label.replace(" ", "_").lower()
    plot_dir = os.path.join(SCRIPT_DIR, "correctness_plots", safe_label)
    os.makedirs(plot_dir, exist_ok=True)

    # Leiden runs late in the pipeline, so its labels are only available
    # post-hoc here -- used to color every embedding plot by cluster.
    py_leiden = (py_snapshots.get("leiden") or {}).get("obs.leiden")
    cd_leiden = (cd_snapshots.get("leiden") or {}).get("obs.leiden")

    # Split stochastic manifolds from linear subspaces
    STOCHASTIC_EMBEDDING_STEPS = ("umap", "tsne")
    PROCRUSTES_STEPS = ("pca", "diffmap")
    
    # No capturable numeric output for these
    NO_PLOTTABLE_OUTPUT = ("calculate_qc_metrics", "filter_cells", "filter_genes")

    @contextmanager
    def side_by_side_figure(out_path, plot_label, figsize=(12, 5)):
        """Open a figure, yield its axes, then save+close on exit."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, sharey=True)
        try:
            yield ax1, ax2
            plt.tight_layout()
            plt.savefig(out_path, dpi=300)
            print(f"[INFO] Saved dot plot for '{plot_label}' -> {out_path}")
        except Exception as exc:
            print(f"[WARN] Failed to plot '{plot_label}': {exc}")
        finally:
            plt.close(fig)

    def flatten_for_plot(values, max_points=3000, seed=0):
        """Coerce a snapshot value to a 1D numeric array, subsampling large
        arrays (e.g. whole X matrices) so scatter plots stay legible/fast."""
        arr = np.asarray(values)
        if arr.dtype.kind == "b":
            arr = arr.astype(np.int8)
        elif arr.dtype.kind in ("U", "O"):
            arr = pd.Categorical(arr.ravel()).codes
        arr = arr.astype(float).ravel()
        if arr.size > max_points:
            rng = np.random.default_rng(seed)
            idx = np.sort(rng.choice(arr.size, size=max_points, replace=False))
            arr = arr[idx]
        return arr

    def scatter_dotplot(ax, values, title, color):
        arr = flatten_for_plot(values)
        ax.scatter(np.arange(arr.size), arr, s=6, alpha=0.5, color=color, edgecolors="none")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("sample index")

    for label, _ in python_steps:
        py_snap = py_snapshots.get(label)
        cd_snap = cd_snapshots.get(label)
        if py_snap is None or cd_snap is None:
            continue  # step failed on one or both sides -- nothing to plot
        if label in NO_PLOTTABLE_OUTPUT:
            continue

        # --- Linear Subspaces: Procrustes Before & After ---
        if label in PROCRUSTES_STEPS:
            obsm_key = f"X_{label}"
            snap_key = f"obsm.{obsm_key}"
            if snap_key not in py_snap or snap_key not in cd_snap:
                continue
                
            py_arr = py_snap[snap_key]
            cd_arr = cd_snap[snap_key]
            
            out_path = os.path.join(plot_dir, f"{label}_procrustes_alignment.png")
            fig, axes = plt.subplots(2, 2, figsize=(14, 12))
            
            # Use Leiden clusters for color if available, otherwise default colors
            c_py, c_cd = 'tab:blue', 'tab:orange'
            cmap = None
            if py_leiden is not None and len(py_leiden) == py_arr.shape[0]:
                c_py = pd.Categorical(py_leiden).codes
                cmap = 'tab20'
            if cd_leiden is not None and len(cd_leiden) == cd_arr.shape[0]:
                c_cd = pd.Categorical(cd_leiden).codes
                cmap = 'tab20'
            
            # [ROW 1] BEFORE PROCRUSTES (Raw Output)
            axes[0, 0].scatter(py_arr[:, 0], py_arr[:, 1], s=4, alpha=0.7, c=c_py, cmap=cmap, edgecolors="none")
            axes[0, 0].set_title(f"Scanpy (Python) - {label.upper()} [Raw Output]")
            axes[0, 0].set_xlabel(f"{label.upper()}1"); axes[0, 0].set_ylabel(f"{label.upper()}2")
            
            axes[0, 1].scatter(cd_arr[:, 0], cd_arr[:, 1], s=4, alpha=0.7, c=c_cd, cmap=cmap, edgecolors="none")
            axes[0, 1].set_title(f"Scancodon (Codon) - {label.upper()} [Raw Output]")
            axes[0, 1].set_xlabel(f"{label.upper()}1"); axes[0, 1].set_ylabel(f"{label.upper()}2")
            
            # [ROW 2] AFTER PROCRUSTES
            from scipy.spatial import procrustes
            try:
                # Procrustes applies optimal translation/rotation to align the matrices
                mtx1, mtx2, disp = procrustes(py_arr, cd_arr)
                
                axes[1, 0].scatter(mtx1[:, 0], mtx1[:, 1], s=4, alpha=0.7, c=c_py, cmap=cmap, edgecolors="none")
                axes[1, 0].set_title(f"Scanpy - [Procrustes Aligned]")
                axes[1, 0].set_xlabel("Aligned Dim 1"); axes[1, 0].set_ylabel("Aligned Dim 2")
                
                axes[1, 1].scatter(mtx2[:, 0], mtx2[:, 1], s=4, alpha=0.7, c=c_cd, cmap=cmap, edgecolors="none")
                axes[1, 1].set_title(f"Scancodon - [Procrustes Aligned]")
                axes[1, 1].set_xlabel("Aligned Dim 1"); axes[1, 1].set_ylabel("Aligned Dim 2")
                
                fig.suptitle(f"{label.upper()} Alignment Diagnostics (Procrustes Disparity: {disp:.4e})", fontsize=16)
            except Exception as e:
                axes[1, 0].text(0.5, 0.5, f"Procrustes failed:\n{e}", ha='center', va='center')
                axes[1, 1].text(0.5, 0.5, f"Procrustes failed:\n{e}", ha='center', va='center')
                fig.suptitle(f"{label.upper()} Alignment Diagnostics (Procrustes Failed)", fontsize=16)
                
            plt.tight_layout(rect=[0, 0.03, 1, 0.96]) # Adjust for suptitle
            plt.savefig(out_path, dpi=300)
            plt.close(fig)
            print(f"[INFO] Saved Procrustes alignment plot for '{label}' -> {out_path}")
            continue

        # --- Stochastic Embeddings: UMAP and t-SNE ---
        if label in STOCHASTIC_EMBEDDING_STEPS:
            obsm_key = f"X_{label}"
            snap_key = f"obsm.{obsm_key}"
            if snap_key not in py_snap or snap_key not in cd_snap:
                continue

            n_cells_py = py_snap["shape"][0]
            n_cells_cd = cd_snap["shape"][0]

            py_obs = pd.DataFrame(index=[str(i) for i in range(n_cells_py)])
            if py_leiden is not None and len(py_leiden) == n_cells_py:
                py_obs["leiden"] = pd.Categorical(py_leiden)

            cd_obs = pd.DataFrame(index=[str(i) for i in range(n_cells_cd)])
            if cd_leiden is not None and len(cd_leiden) == n_cells_cd:
                cd_obs["leiden"] = pd.Categorical(cd_leiden)

            dummy_py = ad.AnnData(X=np.empty((n_cells_py, 1), dtype=np.float32),
                                   obs=py_obs, obsm={obsm_key: py_snap[snap_key]})
            dummy_cd = ad.AnnData(X=np.empty((n_cells_cd, 1), dtype=np.float32),
                                   obs=cd_obs, obsm={obsm_key: cd_snap[snap_key]})
            color_arg = "leiden" if "leiden" in py_obs.columns else None

            out_path = os.path.join(plot_dir, f"{label}_embedding.png")
            with side_by_side_figure(out_path, f"{label} embedding") as (ax1, ax2):
                sp.pl.embedding(dummy_py, basis=label, color=color_arg, ax=ax1,
                                 show=False, title=f"Scanpy (Python) - {label.upper()}")
                sp.pl.embedding(dummy_cd, basis=label, color=color_arg, ax=ax2,
                                 show=False, title=f"Scancodon (Codon) - {label.upper()}")
            continue

        # --- rank_genes_groups: sorted per-gene score comparison ---
        if label == "rank_genes_groups":
            py_res = py_snap.get("uns.rank_genes_groups")
            cd_res = cd_snap.get("uns.rank_genes_groups")
            if not py_res or not cd_res or "names" not in py_res:
                continue
            group_id = py_res["names"].dtype.names[0]
            if group_id not in cd_res["names"].dtype.names:
                continue

            py_scores = np.sort(py_res["scores"][group_id])[::-1]
            cd_scores = np.sort(cd_res["scores"][group_id])[::-1]

            out_path = os.path.join(plot_dir, "rank_genes_groups_scores.png")
            with side_by_side_figure(out_path, f"rank_genes_groups (group '{group_id}')") as (ax1, ax2):
                scatter_dotplot(ax1, py_scores, f"Scanpy (Python) - scores ({group_id})", "tab:blue")
                scatter_dotplot(ax2, cd_scores, f"Scancodon (Codon) - scores ({group_id})", "tab:orange")
                ax1.set_ylabel("score (sorted desc.)")
            continue

        # --- neighbors: graph has no natural embedding, plot node degree instead ---
        if label == "neighbors":
            py_conn = py_snap.get("obsp.connectivities")
            cd_conn = cd_snap.get("obsp.connectivities")
            if py_conn is None or cd_conn is None:
                continue
            py_degree = np.count_nonzero(py_conn, axis=1)
            cd_degree = np.count_nonzero(cd_conn, axis=1)

            out_path = os.path.join(plot_dir, "neighbors_degree.png")
            with side_by_side_figure(out_path, "neighbors (node degree)") as (ax1, ax2):
                scatter_dotplot(ax1, py_degree, "Scanpy (Python) - node degree", "tab:blue")
                scatter_dotplot(ax2, cd_degree, "Scancodon (Codon) - node degree", "tab:orange")
                ax1.set_ylabel("# connections")
            continue

        # --- leiden: cluster id per cell (embeddings above get colored by this too) ---
        if label == "leiden":
            py_labels = py_snap.get("obs.leiden")
            cd_labels = cd_snap.get("obs.leiden")
            if py_labels is None or cd_labels is None:
                continue
            py_codes = pd.Categorical(py_labels).codes
            cd_codes = pd.Categorical(cd_labels).codes

            out_path = os.path.join(plot_dir, "leiden_clusters.png")
            with side_by_side_figure(out_path, "leiden (cluster assignment)") as (ax1, ax2):
                scatter_dotplot(ax1, py_codes, "Scanpy (Python) - cluster id", "tab:blue")
                scatter_dotplot(ax2, cd_codes, "Scancodon (Codon) - cluster id", "tab:orange")
                ax1.set_ylabel("cluster id")
            continue

        # --- everything else (scrublet, normalize_total, log1p, hvg, scale):
        #     one generic value dot plot per captured numeric snapshot key ---
        for key in sorted((set(py_snap) & set(cd_snap)) - {"shape"}):
            py_val, cd_val = py_snap[key], cd_snap[key]
            if isinstance(py_val, str) or isinstance(cd_val, str):
                continue  # "<missing: ...>" placeholders

            safe_key = key.replace(".", "_")
            out_path = os.path.join(plot_dir, f"{label}_{safe_key}.png")
            with side_by_side_figure(out_path, f"{label} / {key}") as (ax1, ax2):
                scatter_dotplot(ax1, py_val, f"Scanpy (Python) - {key}", "tab:blue")
                scatter_dotplot(ax2, cd_val, f"Scancodon (Codon) - {key}", "tab:orange")

    print(f"\n[INFO] All available dot plots saved under: {plot_dir}")
    


def correctness_benchmark_3k_PBMCs():
    """
    Isolated-step correctness/perf benchmark on the 3k PBMC dataset (small, fast baseline).
    Compartmentalized Execution model.
    Tests each Scancodon function using a golden AnnData object to prevent test inconsistencies from polluting later tests
    Compares output of each Scancodon function with the equivalent Scanpy function run on a copy of the golden AnnData object.
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

def run_pca_diagnostic_test():
    """
    Specifically isolates and diagnoses PCA alignment issues between Scanpy and Scancodon:
    1. Sign Ambiguity (Cosine similarity ~ -1.0)
    2. Clustered Eigenvectors (Subspace mixing due to tiny eigenvalue deltas)
    3. Gaussian Projections (Global subspace integrity via Procrustes)
    """
    import numpy as np

    print("\n" + "=" * 80)
    print("--- PCA DIAGNOSTIC: SIGN AMBIGUITY & SUBSPACE MIXING ---")
    print("=" * 80)
    sp, sc = setup_imports()
    
    print("[INFO] Loading and preprocessing dataset...")
    adata = sp.datasets.pbmc3k()
    adata.var_names_make_unique()
    
    # Preprocess so PCA has meaningful variance to capture
    sp.pp.filter_cells(adata, min_genes=200)
    sp.pp.filter_genes(adata, min_cells=3)
    sp.pp.normalize_total(adata, target_sum=1e4)
    sp.pp.log1p(adata)
    sp.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable].copy()
    sp.pp.scale(adata, max_value=10)
    
    adata_sp = adata.copy()
    adata_sc = adata.copy()
    
    print("[INFO] Executing PCA implementations...")
    # Run PCA (using default randomized SVD solver for both)
    sp.tl.pca(adata_sp, n_comps=20)
    sc.tl.pca(adata_sc, n_comps=20)
    
    # Extract Principal Components (eigenvectors) and Eigenvalues
    loadings_py = adata_sp.varm['PCs']
    loadings_cd = adata_sc.varm['PCs']
    evals_py = adata_sp.uns['pca']['variance']
    
    # ---------------------------------------------------------
    # 1. Sign Ambiguity & Component Alignment Check
    # ---------------------------------------------------------
    print("\n[1] Component Alignment Analysis")
    print(f"{'PC':<4} | {'Cosine Sim':<12} | {'Eigenvalue':<12} | {'Eval Delta (to next)':<22} | {'Status'}")
    print("-" * 80)
    
    for i in range(loadings_py.shape[1]):
        v_py = loadings_py[:, i]
        v_cd = loadings_cd[:, i]
        
        # Calculate Cosine Similarity to check alignment direction
        cos_sim = np.dot(v_py, v_cd) / (np.linalg.norm(v_py) * np.linalg.norm(v_cd))
        
        # Calculate the gap between current and next eigenvalue
        eval_curr = evals_py[i]
        eval_delta = eval_curr - evals_py[i+1] if i < len(evals_py)-1 else np.nan
        
        # Categorize the structural phenomenon
        status = "Aligned"
        if np.isclose(cos_sim, -1.0, atol=1e-2):
            status = "Sign Flipped"
        elif np.isclose(cos_sim, 1.0, atol=1e-2):
            status = "Exact Match"
        elif abs(cos_sim) < 0.99:
            # If the delta to the next eigenvalue is less than 5% of the current eigenvalue's magnitude,
            # the solver likely mixed the adjacent components.
            if not np.isnan(eval_delta) and eval_delta < (0.05 * eval_curr):
                status = "Mixed (Clustered Eigenvalues)"
            else:
                status = "Rotated (Stochastic Projection)"
                
        print(f"{i:<4} | {cos_sim:>12.4f} | {eval_curr:>12.4f} | {eval_delta:>22.4f} | {status}")

    # ---------------------------------------------------------
    # 2. Global Subspace Integrity Check (Procrustes)
    # ---------------------------------------------------------
    print("\n[2] Global Subspace Integrity")
    from scipy.spatial import procrustes
    try:
        _, _, disparity = procrustes(adata_sp.obsm['X_pca'], adata_sc.obsm['X_pca'])
        print(f"  Procrustes Disparity Score: {disparity:.6e}")
        if disparity < 1e-2:
            print("  Conclusion: The multidimensional manifolds are structurally identical.\n"
                  "              Any element-wise deviations in X_pca are strictly due to \n"
                  "              sign flipping or orthogonal rotations within the subspace.")
        else:
            print("  Conclusion: The multidimensional manifolds diverge structurally. \n"
                  "              The implementations are calculating fundamentally different matrices.")
    except Exception as e:
        print(f"  Procrustes calculation failed: {e}")
        
    print("=" * 80 + "\n")

def run_bridging_benchmark():
    print("Loading 3k PBMCs dataset...")
    sp, sc = setup_imports()
    # Using the standard 3k PBMC dataset for realistic biological data geometry
    adata = sp.datasets.pbmc3k()
    
    # Run the benchmark via the Tools class alias
    sc.tl.benchmark_bridge_comprehensive(adata)

def run_dense_sparse_kernel_benchmark():
    print("\nLoading 3k PBMCs dataset for Kernel Benchmarking...")
    sp, sc = setup_imports()
    # Using the standard 3k PBMC dataset for realistic biological data geometry
    adata = sp.datasets.pbmc3k()
    
    # Run the kernel benchmark via the Tools class alias
    sc.tl.benchmark_dense_sparse_kernels(adata)


def run_manifold_trustworthiness_test():
    """
    Evaluates the local topological preservation of UMAP and t-SNE embeddings
    by calculating the trustworthiness score against the high-dimensional PCA space.
    """
    import numpy as np
    from sklearn.manifold import trustworthiness

    print("\n" + "=" * 80)
    print("--- MANIFOLD TRUSTWORTHINESS DIAGNOSTIC: UMAP & T-SNE ---")
    print("=" * 80)
    sp, sc = setup_imports()
    
    print("[INFO] Loading and preprocessing dataset...")
    adata = sp.datasets.pbmc3k()
    adata.var_names_make_unique()
    
    # Standard preprocessing up to PCA
    sp.pp.filter_cells(adata, min_genes=200)
    sp.pp.filter_genes(adata, min_cells=3)
    sp.pp.normalize_total(adata, target_sum=1e4)
    sp.pp.log1p(adata)
    sp.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable].copy()
    sp.pp.scale(adata, max_value=10)
    
    print("[INFO] Computing High-Dimensional Reference (PCA)...")
    sp.tl.pca(adata, n_comps=40)
    X_ref = adata.obsm['X_pca']
    
    print("[INFO] Computing Neighborhood Graph...")
    sp.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
    
    adata_sp = adata.copy()
    adata_sc = adata.copy()
    
    print("[INFO] Executing UMAP and t-SNE implementations (Scanpy vs Scancodon)...")
    
    # Scanpy Manifolds
    t0 = time.perf_counter()
    sp.tl.umap(adata_sp)
    time_umap_sp = time.perf_counter() - t0
    
    t0 = time.perf_counter()
    sp.tl.tsne(adata_sp)
    time_tsne_sp = time.perf_counter() - t0
    
    # Scancodon Manifolds
    t0 = time.perf_counter()
    sc.tl.umap(adata_sc)
    time_umap_sc = time.perf_counter() - t0
    
    t0 = time.perf_counter()
    sc.tl.tsne(adata_sc)
    time_tsne_sc = time.perf_counter() - t0
    
    print("[INFO] Calculating Trustworthiness Scores (n_neighbors=15)...")
    # Trustworthiness scores how well the local 2D neighborhoods reflect the PCA space.
    
    tw_umap_sp = trustworthiness(X_ref, adata_sp.obsm['X_umap'], n_neighbors=15)
    tw_umap_sc = trustworthiness(X_ref, adata_sc.obsm['X_umap'], n_neighbors=15)
    
    tw_tsne_sp = trustworthiness(X_ref, adata_sp.obsm['X_tsne'], n_neighbors=15)
    tw_tsne_sc = trustworthiness(X_ref, adata_sc.obsm['X_tsne'], n_neighbors=15)
    
    print("\n[RESULTS] Local Topological Preservation")
    print(f"{'Algorithm':<12} | {'Implementation':<18} | {'Trustworthiness':<16} | {'Time (s)'}")
    print("-" * 65)
    print(f"{'UMAP':<12} | {'Scanpy (Python)':<18} | {tw_umap_sp:<16.4f} | {time_umap_sp:.2f}")
    print(f"{'UMAP':<12} | {'Scancodon (Codon)':<18} | {tw_umap_sc:<16.4f} | {time_umap_sc:.2f}")
    print("-" * 65)
    print(f"{'t-SNE':<12} | {'Scanpy (Python)':<18} | {tw_tsne_sp:<16.4f} | {time_tsne_sp:.2f}")
    print(f"{'t-SNE':<12} | {'Scancodon (Codon)':<18} | {tw_tsne_sc:<16.4f} | {time_tsne_sc:.2f}")
    print("=" * 80 + "\n")


if __name__ == "__main__":

    #correctness_benchmark_3k_PBMCs()

    # Heavy: downloads ~4GB and may consume tens of GB of RAM / crash. Uncomment to run.
    #correctness_benchmark_1M_neurons()

    #evaluate_RGG_edge_cases()

    #run_bridging_benchmark()
    run_dense_sparse_kernel_benchmark()
    #run_pca_diagnostic_test()
    #run_manifold_trustworthiness_test()