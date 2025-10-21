#!/usr/bin/env python3
"""
Evaluation script for sequence alignment algorithms.
Tests both Python and Codon implementations.
Reports timing in milliseconds for each test.
"""

import os
import sys
import time
import subprocess


def read_fasta_sequences(filepath):
    """
    Read all sequences from a multi-FASTA file.
    
    Args:
        filepath: Path to the FASTA file
        
    Returns:
        list of tuples: [(header, sequence), ...]
    """
    sequences = []
    current_header = None
    current_seq = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_seq)))
                # Start new sequence
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_header is not None:
            sequences.append((current_header, ''.join(current_seq)))
    
    return sequences


def run_python_tests():
    """Run Python implementation tests"""
    # Add the week4 directory to the path
    week4_dir = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, week4_dir)
    
    from code.python import global_align, local_align, semiglobal_align, affine_align, read_fasta
    
    results = []
    data_dir = os.path.join(week4_dir, "data")
    
    # Test on q/t pairs
    q_sequences = read_fasta_sequences(os.path.join(data_dir, "q1.fa"))
    t_sequences = read_fasta_sequences(os.path.join(data_dir, "t1.fa"))
    
    for i, ((q_header, q_seq), (t_header, t_seq)) in enumerate(zip(q_sequences, t_sequences), 1):
        # Global alignment
        start = time.time()
        score, _, _ = global_align(t_seq, q_seq)
        runtime = max(1, int((time.time() - start) * 1000))
        results.append((f"global-q{i}", "python", runtime))
        
        # Local alignment
        start = time.time()
        score, _, _ = local_align(t_seq, q_seq)
        runtime = max(1, int((time.time() - start) * 1000))
        results.append((f"local-q{i}", "python", runtime))
        
        # Fitting alignment
        start = time.time()
        score, _, _ = semiglobal_align(t_seq, q_seq)
        runtime = max(1, int((time.time() - start) * 1000))
        results.append((f"fitting-q{i}", "python", runtime))
        
        # Affine alignment
        start = time.time()
        score, _, _ = affine_align(t_seq, q_seq, gap_open=-5, gap_extend=-1)
        runtime = max(1, int((time.time() - start) * 1000))
        results.append((f"affine-q{i}", "python", runtime))
    
    # Test on MT-human vs MT-orang
    human_seq = read_fasta(os.path.join(data_dir, "MT-human.fa"))
    orang_seq = read_fasta(os.path.join(data_dir, "MT-orang.fa"))
    
    # Global alignment
    start = time.time()
    score, _, _ = global_align(human_seq, orang_seq)
    runtime = int((time.time() - start) * 1000)
    results.append(("global-mt_human", "python", runtime))
    
    # Local alignment
    start = time.time()
    score, _, _ = local_align(human_seq, orang_seq)
    runtime = int((time.time() - start) * 1000)
    results.append(("local-mt_human", "python", runtime))
    
    # Fitting alignment
    start = time.time()
    score, _, _ = semiglobal_align(human_seq, orang_seq)
    runtime = int((time.time() - start) * 1000)
    results.append(("fitting-mt_human", "python", runtime))
    
    # Affine alignment
    start = time.time()
    score, _, _ = affine_align(human_seq, orang_seq, gap_open=-5, gap_extend=-1)
    runtime = int((time.time() - start) * 1000)
    results.append(("affine-mt_human", "python", runtime))
    
    return results


def run_codon_tests():
    """Run Codon implementation tests by creating and executing inline Codon code"""
    results = []
    week4_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Check if codon is available
    codon_path = "codon"
    if os.environ.get('GITHUB_ACTIONS'):
        potential_codon = os.path.expanduser("~/.codon/bin/codon")
        if os.path.exists(potential_codon):
            codon_path = potential_codon
    
    # Create inline Codon test script
    codon_test_code = '''#!/usr/bin/env codon
import time
from code.codon.utils import read_fasta
from code.codon.global_alignment import global_align
from code.codon.local_alignment import local_align
from code.codon.semiglobal_alignment import semiglobal_align
from code.codon.affine_alignment import affine_align

def read_fasta_sequences(filepath: str):
    sequences = []
    current_header = ""
    current_seq = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header != "":
                    sequences.append((current_header, ''.join(current_seq)))
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_header != "":
            sequences.append((current_header, ''.join(current_seq)))
    
    return sequences

data_dir = "data"

# Test on q/t pairs
q_sequences = read_fasta_sequences(f"{data_dir}/q1.fa")
t_sequences = read_fasta_sequences(f"{data_dir}/t1.fa")

for i in range(len(q_sequences)):
    q_header, q_seq = q_sequences[i]
    t_header, t_seq = t_sequences[i]
    test_num = i + 1
    
    # Global alignment
    start = time.time()
    score, _, _ = global_align(t_seq, q_seq)
    runtime = max(1, int((time.time() - start) * 1000.0))
    print(f"global-q{test_num}|codon|{runtime}")
    
    # Local alignment
    start = time.time()
    score, _, _ = local_align(t_seq, q_seq)
    runtime = max(1, int((time.time() - start) * 1000.0))
    print(f"local-q{test_num}|codon|{runtime}")
    
    # Fitting alignment
    start = time.time()
    score, _, _ = semiglobal_align(t_seq, q_seq)
    runtime = max(1, int((time.time() - start) * 1000.0))
    print(f"fitting-q{test_num}|codon|{runtime}")
    
    # Affine alignment
    start = time.time()
    score, _, _ = affine_align(t_seq, q_seq, gap_open=-5, gap_extend=-1)
    runtime = max(1, int((time.time() - start) * 1000.0))
    print(f"affine-q{test_num}|codon|{runtime}")

# Test on MT-human vs MT-orang
human_seq = read_fasta(f"{data_dir}/MT-human.fa")
orang_seq = read_fasta(f"{data_dir}/MT-orang.fa")

# Global alignment
start = time.time()
score, _, _ = global_align(human_seq, orang_seq)
runtime = int((time.time() - start) * 1000.0)
print(f"global-mt_human|codon|{runtime}")

# Local alignment
start = time.time()
score, _, _ = local_align(human_seq, orang_seq)
runtime = int((time.time() - start) * 1000.0)
print(f"local-mt_human|codon|{runtime}")

# Fitting alignment
start = time.time()
score, _, _ = semiglobal_align(human_seq, orang_seq)
runtime = int((time.time() - start) * 1000.0)
print(f"fitting-mt_human|codon|{runtime}")

# Affine alignment
start = time.time()
score, _, _ = affine_align(human_seq, orang_seq, gap_open=-5, gap_extend=-1)
runtime = int((time.time() - start) * 1000.0)
print(f"affine-mt_human|codon|{runtime}")
'''
    
    try:
        # Write the Codon code to a temporary file
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.codon', delete=False) as f:
            temp_codon_file = f.name
            f.write(codon_test_code)
        
        # Run the Codon test script
        result = subprocess.run(
            [codon_path, "run", temp_codon_file],
            cwd=week4_dir,
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )
        
        # Clean up temp file
        os.unlink(temp_codon_file)
        
        if result.returncode != 0:
            print(f"Warning: Codon tests failed with error:\n{result.stderr}", file=sys.stderr)
            return results
        
        # Parse the output - format: "method|language|runtime"
        lines = result.stdout.strip().split('\n')
        for line in lines:
            if '|' in line:
                parts = line.split('|')
                if len(parts) == 3:
                    method, language, runtime_str = parts
                    try:
                        runtime = int(runtime_str)
                        results.append((method, language, runtime))
                    except ValueError:
                        continue
        
    except subprocess.TimeoutExpired:
        print("Warning: Codon tests timed out", file=sys.stderr)
        if 'temp_codon_file' in locals():
            try:
                os.unlink(temp_codon_file)
            except:
                pass
    except FileNotFoundError:
        print(f"Warning: Codon not found at {codon_path}, skipping Codon tests", file=sys.stderr)
    except Exception as e:
        print(f"Warning: Error running Codon tests: {e}", file=sys.stderr)
        if 'temp_codon_file' in locals():
            try:
                os.unlink(temp_codon_file)
            except:
                pass
    
    return results


def main():
    """Main evaluation function"""
    try:
        # Run Python tests
        python_results = run_python_tests()
        
        # Run Codon tests
        codon_results = run_codon_tests()
        
        # Combine results
        all_results = python_results + codon_results
        
        # Print results in required format
        print(f"{'Method':<20} {'Language':<10} {'Runtime'}")
        print("-" * 38)
        
        for method, language, runtime in all_results:
            print(f"{method:<20} {language:<10} {runtime}ms")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
