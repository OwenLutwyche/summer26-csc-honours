import sys
import os
import importlib.util
import time

# Ensure we are in the right directory (optional safety check)
required_path = "/Users/oweno/Desktop/honours/fall25-csc-bioinf/project/scanpy-main"
if os.getcwd() != required_path:
    print(f"Changing directory to {required_path}")
    try:
        os.chdir(required_path)
    except FileNotFoundError:
        print(f"ERROR: Could not find directory: {required_path}")
        sys.exit(1)

# Set environment variables
os.environ["PYTHONPATH"] = "src"
os.environ["MPLBACKEND"] = "Agg"

# Add src to path so scanpy can be imported
sys.path.insert(0, os.path.join(required_path, "src"))


def load_module_from_file(filepath: str):
    """Dynamically load a Python module from a file path."""
    module_name = os.path.basename(filepath).replace(".py", "")
    spec = importlib.util.spec_from_file_location(module_name, filepath)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load module from {filepath}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def run_evaluation():
    print("=" * 60)
    print("--- VERIFYING BASELINE SCANPY (PYTHON) ---")
    print("=" * 60)

    try:
        import scanpy as baseline_scanpy
        print(f"[OK] Loaded baseline 'scanpy' package from: {os.path.dirname(baseline_scanpy.__file__)}")
    except ImportError as exc:
        print("[ERROR] CRITICAL: Could not import baseline 'scanpy'.")
        print(f"   Error: {exc}")
        sys.exit(1)

    print("[INFO] Using upstream Python Scanpy implementation (no hijack)")
    print("-" * 60)

    test_files = [
        "new_tests/test_preprocessing.py",
        "new_tests/test_neighbors.py",
        "new_tests/test_clustering.py",
        "new_tests/test_embedding.py",
        "new_tests/test_rank_genes_groups.py",
    ]

    missing_files = [tf for tf in test_files if not os.path.exists(tf)]
    if missing_files:
        print("[ERROR] Missing test files:")
        for missing in missing_files:
            print(f"   - {missing}")
        sys.exit(1)

    total_passed = 0
    total_failed = 0
    start_time = time.time()
    per_file_timings = []

    for test_file in test_files:
        print(f"\n[INFO] Running {test_file}...")
        file_start = time.perf_counter()
        try:
            module = load_module_from_file(test_file)
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
        except Exception as exc:
            print(f"  [CRASH] {exc}")
            import traceback
            traceback.print_exc()
            total_failed += 1
        file_elapsed = time.perf_counter() - file_start
        per_file_timings.append((test_file, file_elapsed))
        print(f"  [TIME] {file_elapsed:.2f}s")

    total_elapsed = time.time() - start_time

    print("\n" + "=" * 60)
    print(f"SUMMARY: {total_passed} Passed, {total_failed} Failed")
    print(f"Time: {total_elapsed:.2f}s")
    print("=" * 60)

    if per_file_timings:
        print("Per-file timings (not added to total):")
        for test_file, duration in per_file_timings:
            print(f"  {test_file}: {duration:.2f}s")

    if total_failed > 0:
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    run_evaluation()