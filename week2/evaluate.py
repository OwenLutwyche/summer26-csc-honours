#!/usr/bin/env python3
import os
import subprocess
import sys


def run_test(description, command, cwd, env_vars=None):
    """Run a test command and report results"""
    print(f"\n{'='*60}")
    print(f"{description}")
    print(f"{'='*60}")
    
    env = os.environ.copy()
    if env_vars:
        env.update(env_vars)
    
    try:
        result = subprocess.run(
            command,
            shell=True,
            cwd=cwd,
            env=env,
            text=True,
            capture_output=True,
            executable='/bin/bash'
        )
        
        # Print output
        if result.stdout:
            print(result.stdout)
        
        if result.returncode != 0:
            print(f"\n[FAIL] {description} (exit code: {result.returncode})")
            if result.stderr:
                print(f"Error output:\n{result.stderr}")
            return False
        else:
            print(f"\n[PASS] {description}")
            return True
            
    except Exception as e:
        print(f"\n[FAIL] {description} with exception: {e}")
        return False


def main():
    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    results = []
    
    # Test 1: Run Python tests
    print("\n" + "="*60)
    print("WEEK 2 EVALUATION")
    print("="*60)
    
    success = run_test(
        "Python Test Suite",
        "python3 test.py",
        script_dir
    )
    results.append(("Python Tests", success))
    
    # Test 2: Run Codon tests
    # Set up environment variables for Codon
    codon_env = {}
    
    # Find codon binary
    codon_path = os.path.expanduser("~/.codon/bin")
    if os.path.exists(codon_path):
        path_env = os.environ.get('PATH', '')
        codon_env['PATH'] = f"{codon_path}:{path_env}"
    
    # Set CODON_PYTHON
    try:
        result = subprocess.run(
            "find_libpython",
            shell=True,
            capture_output=True,
            text=True
        )
        if result.returncode == 0 and result.stdout.strip():
            codon_env['CODON_PYTHON'] = result.stdout.strip()
    except:
        pass
    
    # Set PYTHON_PATH
    try:
        result = subprocess.run(
            "python3 -c \"import site; print(site.getsitepackages()[0])\"",
            shell=True,
            capture_output=True,
            text=True
        )
        if result.returncode == 0 and result.stdout.strip():
            codon_env['PYTHON_PATH'] = result.stdout.strip()
    except:
        pass
    
    success = run_test(
        "Codon Test Suite",
        "codon run test.py",
        script_dir,
        codon_env
    )
    results.append(("Codon Tests", success))
    
    # Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    for test_name, passed in results:
        status = "[PASS]" if passed else "[FAIL]"
        print(f"{test_name}: {status}")
    
    # Exit with error if any test failed
    if not all(passed for _, passed in results):
        print("\n[FAIL] Some tests failed")
        sys.exit(1)
    else:
        print("\n[PASS] All tests passed")
        sys.exit(0)


if __name__ == "__main__":
    main()
