#!/usr/bin/env python3
"""
Evaluate Codon Port (scancodon) by hijacking standard tests.
"""
import sys
import os
import time
import importlib.util

# Ensure we are in the right directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(SCRIPT_DIR)

# Set matplotlib backend to non-interactive
os.environ["MPLBACKEND"] = "Agg"

def run_evaluation():
    print("=" * 60)
    print("--- VERIFYING CODON PORT (SCANCODON) ---")
    print("=" * 60)

    # 1. Import the wrapper package
    try:
        import scancodon
        print(f"[OK] Loaded 'scancodon' wrapper from: {os.path.dirname(scancodon.__file__)}")
    except ImportError as e:
        print(f"[ERROR] CRITICAL: Could not import 'scancodon' package.")
        print(f"   Error: {e}")
        print("   Ensure you have a folder named 'scancodon' with an '__init__.py'.")
        sys.exit(1)

    # 2. Verify Native Connection
    if getattr(scancodon, "CODON_AVAILABLE", False):
        print("[STATUS] Running on NATIVE CODON KERNELS")
    else:
        print("[WARNING] STATUS: Running on NUMPY FALLBACKS (Native extension not loaded)")
        print("   Did you run ./debug_build.sh successfully?")

    # 3. THE HIJACK
    sys.modules['scanpy'] = scancodon
    print("[OK] Hijacked 'scanpy' module -> pointing to 'scancodon'")
    print("-" * 60)

    # 4. Define Tests
    test_files = [
        "tests/test_preprocessing.py",
        "tests/test_neighbors.py",
        "tests/test_clustering.py",
        "tests/test_embedding.py",
        "tests/test_rank_genes_groups.py",
    ]

    total_passed = 0
    total_failed = 0
    start_time = time.time()

    # Helper to load test modules dynamically
    def load_test_module(filepath):
        module_name = "test_mod_" + os.path.basename(filepath).replace(".py", "")
        spec = importlib.util.spec_from_file_location(module_name, filepath)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    for f in test_files:
        print(f"\n[INFO] Running {f}...")
        if not os.path.exists(f):
            print(f"   [WARN] File not found: {f}")
            continue
            
        try:
            # Load the test file
            mod = load_test_module(f)
            
            if hasattr(mod, 'run_all'):
                results = mod.run_all()
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

    total_time = time.time() - start_time

    print("\n" + "=" * 60)
    print(f"SUMMARY: {total_passed} Passed, {total_failed} Failed")
    print(f"Time: {total_time:.2f}s")
    print("=" * 60)

    if total_failed > 0:
        sys.exit(1)

if __name__ == "__main__":
    run_evaluation()