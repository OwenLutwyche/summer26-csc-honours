import sys
import os
import importlib.util
import time

# Ensure we are in the right directory (optional safety check)
required_path = "/home/ryan/fall25-csc-bioinf/project/scanpy-main"
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


def run_baseline_tests():
    print("=" * 60)
    print("--- RUNNING BASELINE (PYTHON) SCANPY TESTS ---")
    print("=" * 60)
    
    # Target the new standalone test files
    test_files = [
        "new_tests/test_preprocessing.py",
        "new_tests/test_neighbors.py",
        "new_tests/test_clustering.py", 
        "new_tests/test_embedding.py",
        "new_tests/test_rank_genes_groups.py",
    ]

    # Verify all test files exist
    missing_files = []
    for test_file in test_files:
        if not os.path.exists(test_file):
            missing_files.append(test_file)
    
    if missing_files:
        print("WARNING: The following test files are missing:")
        for f in missing_files:
            print(f"   - {f}")
        print("\nAborting test run. Please check the file paths.")
        sys.exit(1)
    else:
        print(f"OK: All {len(test_files)} test files found.\n")
    
    # Stats
    total_passed = 0
    total_failed = 0
    total_errors = 0
    failed_tests = []
    
    total_start = time.time()
    
    # Load and run tests from each file
    for test_file in test_files:
        print(f"{'-' * 50}")
        print(f"FILE: {test_file}")
        print(f"{'-' * 50}")
        
        try:
            module = load_module_from_file(test_file)
            
            if hasattr(module, 'run_all'):
                results = module.run_all()
                for name, success, msg in results:
                    if success:
                        total_passed += 1
                        print(f"  PASS: {name}: {msg}")
                    else:
                        if "FAILED" in msg:
                            total_failed += 1
                        else:
                            total_errors += 1
                        failed_tests.append((test_file, name, msg))
                        print(f"  FAIL: {name}: {msg}")
            else:
                print("  WARNING: No run_all() function found")
                    
        except Exception as e:
            import traceback
            print(f"  ERROR loading module: {e}")
            traceback.print_exc()
            total_errors += 1
    
    total_elapsed = time.time() - total_start
    
    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Passed:  {total_passed}")
    print(f"  Failed:  {total_failed}")
    print(f"  Errors:  {total_errors}")
    print(f"  Time:    {total_elapsed:.2f}s")
    
    if failed_tests:
        print("\n" + "-" * 60)
        print("FAILED TESTS:")
        for file, name, msg in failed_tests:
            print(f"  {file}::{name}")
            print(f"    {msg}")
    
    print("=" * 60)
    
    # Exit with error code if any failures
    if total_failed > 0 or total_errors > 0:
        sys.exit(1)
    sys.exit(0)

if __name__ == "__main__":
    run_baseline_tests()