#!/usr/bin/env python3
"""
Evaluation script for phylogenetic tree construction.
Tests both Python (biotite) and Codon implementations.
Reports timing in milliseconds for each test.
"""

import os
import subprocess
import time
import sys


def run_python_tests():
    """Run Python tests using biotite"""
    import numpy as np
    import biotite.sequence.phylo as phylo
    
    def get_distances():
        """Load test distances matrix from file."""
        week3_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(week3_dir, "data")
        distances = np.loadtxt(os.path.join(data_dir, "distances.txt"), dtype=float)
        return distances
    
    def get_upgma_newick():
        """Load the Newick notation for the UPGMA tree from file."""
        week3_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(week3_dir, "data")
        with open(os.path.join(data_dir, "newick_upgma.txt"), "r") as file:
            newick = file.read().strip()
        return newick
    
    def get_tree():
        """Create a phylogenetic tree using UPGMA"""
        distances = get_distances()
        return phylo.upgma(distances)
    
    def test_distances():
        """Test whether the distance_to() and topological distances work correctly."""
        tree = get_tree()
        
        # Tree is created via UPGMA
        # -> The distances to root should be equal for all leaf nodes
        dist = tree.root.distance_to(tree.leaves[0])
        for leaf in tree.leaves:
            assert leaf.distance_to(tree.root) == dist
        
        # Example topological distances
        assert tree.get_distance(0, 19, True) == 9
        assert tree.get_distance(4, 2, True) == 10
    
    def test_upgma():
        """Compare the results of upgma() with expected tree structure."""
        tree = get_tree()
        upgma_newick = get_upgma_newick()
        
        ref_tree = phylo.Tree.from_newick(upgma_newick)
        
        # Cannot apply direct tree equality assertion because the distance
        # might not be exactly equal due to floating point rounding errors
        for i in range(len(tree)):
            for j in range(len(tree)):
                # Check for equal distances and equal topologies
                dist_diff = abs(tree.get_distance(i, j) - ref_tree.get_distance(i, j))
                assert dist_diff < 1e-3, f"Distance mismatch at ({i},{j}): {dist_diff}"
                
                assert tree.get_distance(i, j, topological=True) == ref_tree.get_distance(
                    i, j, topological=True
                )
    
    def test_neighbor_joining():
        """Compare the results of neighbor_join() with a known tree."""
        dist = np.array([
            [ 0,  5,  4,  7,  6,  8],
            [ 5,  0,  7, 10,  9, 11],
            [ 4,  7,  0,  7,  6,  8],
            [ 7, 10,  7,  0,  5,  9],
            [ 6,  9,  6,  5,  0,  8],
            [ 8, 11,  8,  9,  8,  0],
        ])
        
        ref_tree = phylo.Tree(
            phylo.TreeNode(
                [
                    phylo.TreeNode(
                        [
                            phylo.TreeNode(
                                [
                                    phylo.TreeNode(index=0),
                                    phylo.TreeNode(index=1),
                                ],
                                [1, 4],
                            ),
                            phylo.TreeNode(index=2),
                        ],
                        [1, 2],
                    ),
                    phylo.TreeNode(
                        [
                            phylo.TreeNode(index=3),
                            phylo.TreeNode(index=4),
                        ],
                        [3, 2],
                    ),
                    phylo.TreeNode(index=5),
                ],
                [1, 1, 5],
            )
        )
        
        test_tree = phylo.neighbor_joining(dist)
        assert test_tree == ref_tree
    
    try:
        # Run all tests and time them
        start = time.time()
        test_distances()
        time_distances = (time.time() - start) * 1000
        
        start = time.time()
        test_upgma()
        time_upgma = (time.time() - start) * 1000
        
        start = time.time()
        test_neighbor_joining()
        time_neighbor = (time.time() - start) * 1000
        
        total_time = time_distances + time_upgma + time_neighbor
        return total_time, True
        
    except Exception as e:
        print(f"Python test error: {e}", file=sys.stderr)
        return None, False


def run_codon_tests():
    """Run Codon tests using biotite_codon"""
    # Get the week3 directory
    week3_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Check if codon is in PATH, otherwise use the expected GitHub Actions location
    codon_path = "codon"
    if os.environ.get('GITHUB_ACTIONS'):
        potential_codon = os.path.expanduser("~/.codon/bin/codon")
        if os.path.exists(potential_codon):
            codon_path = potential_codon
    
    # Create a test runner Codon file
    test_runner_path = os.path.join(week3_dir, "test_runner.codon")
    
    # Write the test runner
    with open(test_runner_path, "w") as f:
        f.write("""#!/usr/bin/env codon
\"\"\"
Codon test runner for phylogenetic tree construction.
\"\"\"

import time
from python import numpy as pnp
import numpy.pybridge
import numpy as np

# Import our Codon implementation
import code.biotite_codon as phylo


def get_distances():
    \"\"\"Load test distances matrix from file.\"\"\"
    distances: np.ndarray[int,2] = pnp.loadtxt("data/distances.txt", dtype=pnp.int64)
    return distances


def get_upgma_newick():
    \"\"\"Load the Newick notation for the UPGMA tree from file.\"\"\"
    with open("data/newick_upgma.txt", "r") as file:
        newick = file.read().strip()
    return newick


def get_tree():
    \"\"\"Create a phylogenetic tree using UPGMA\"\"\"
    distances = get_distances()
    return phylo.upgma(distances)


def test_distances():
    \"\"\"Test whether the distance_to() and topological distances work correctly.\"\"\"
    tree = get_tree()
    
    # Tree is created via UPGMA
    # -> The distances to root should be equal for all leaf nodes
    dist = tree.root.distance_to(tree.leaves[0])
    for leaf in tree.leaves:
        assert leaf.distance_to(tree.root) == dist, "Distance to root should be equal for all leaves"
    
    # Example topological distances
    assert tree.get_distance(0, 19, True) == 9, f"Expected topological distance 9, got {tree.get_distance(0, 19, True)}"
    assert tree.get_distance(4, 2, True) == 10, f"Expected topological distance 10, got {tree.get_distance(4, 2, True)}"
    
    print("test_distances passed")


def test_upgma():
    \"\"\"Compare the results of upgma() with expected tree structure.\"\"\"
    tree = get_tree()
    upgma_newick = get_upgma_newick()
    
    ref_tree = phylo.Tree.from_newick(upgma_newick)
    
    # Cannot apply direct tree equality assertion because the distance
    # might not be exactly equal due to floating point rounding errors
    for i in range(len(tree)):
        for j in range(len(tree)):
            # Check for equal distances and equal topologies
            dist_diff = abs(tree.get_distance(i, j) - ref_tree.get_distance(i, j))
            assert dist_diff < 1e-3, f"Distance mismatch at ({i},{j}): {dist_diff}"
            
            assert tree.get_distance(i, j, topological=True) == ref_tree.get_distance(
                i, j, topological=True
            ), f"Topological distance mismatch at ({i},{j})"
    
    print("test_upgma passed")


def test_neighbor_joining():
    \"\"\"Compare the results of neighbor_join() with a known tree.\"\"\"
    dist = np.array([
        [ 0.,  5.,  4.,  7.,  6.,  8.],
        [ 5.,  0.,  7., 10.,  9., 11.],
        [ 4.,  7.,  0.,  7.,  6.,  8.],
        [ 7., 10.,  7.,  0.,  5.,  9.],
        [ 6.,  9.,  6.,  5.,  0.,  8.],
        [ 8., 11.,  8.,  9.,  8.,  0.],
    ], dtype=np.float64)

    ref_tree = phylo.Tree(
        phylo.TreeNode(
            [
                phylo.TreeNode(
                    [
                        phylo.TreeNode(
                            [
                                phylo.TreeNode(index=0),
                                phylo.TreeNode(index=1),
                            ],
                            [1.0, 4.0],
                        ),
                        phylo.TreeNode(index=2),
                    ],
                    [1.0, 2.0],
                ),
                phylo.TreeNode(
                    [
                        phylo.TreeNode(index=3),
                        phylo.TreeNode(index=4),
                    ],
                    [3.0, 2.0],
                ),
                phylo.TreeNode(index=5),
            ],
            [1.0, 1.0, 5.0],
        )
    )

    test_tree = phylo.neighbor_joining(dist)

    # Check topological equivalence by comparing all pairwise distances
    # This is more robust than exact tree structure comparison
    for i in range(6):
        for j in range(6):
            assert test_tree.get_distance(i, j) == ref_tree.get_distance(i, j), \
                f"Distance mismatch at ({i},{j})"
            assert test_tree.get_distance(i, j, topological=True) == ref_tree.get_distance(i, j, topological=True), \
                f"Topological distance mismatch at ({i},{j})"
    
    print("test_neighbor_joining passed")


def run_tests():
    \"\"\"Run all tests and report timing\"\"\"
    print("Running Codon phylo tests...")
    print("=" + "=" * 39)
    
    # Test 1: test_distances
    start = time.time()
    test_distances()
    time_distances = (time.time() - start) * 1000
    
    # Test 2: test_upgma
    start = time.time()
    test_upgma()
    time_upgma = (time.time() - start) * 1000
    
    # Test 3: test_neighbor_joining
    start = time.time()
    test_neighbor_joining()
    time_neighbor = (time.time() - start) * 1000
    
    total_time = time_distances + time_upgma + time_neighbor
    
    print("=" + "=" * 39)
    print(f"test_distances:        {int(time_distances)}ms")
    print(f"test_upgma:            {int(time_upgma)}ms")
    print(f"test_neighbor_joining: {int(time_neighbor)}ms")
    print(f"Total:                 {int(total_time)}ms")
    print("=" + "=" * 39)
    
    return total_time


if __name__ == "__main__":
    try:
        total_time = run_tests()
        print(f"\\nAll tests passed! Total time: {int(total_time)}ms")
    except Exception as e:
        print(f"\\nTest failed: {e}")
        import sys
        sys.exit(1)
""")
    
    # Run the Codon test
    try:
        result = subprocess.run(
            [codon_path, "run", "-release", "test_runner.codon"],
            cwd=week3_dir,
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )
        
        if result.returncode == 0:
            # Parse the total time from the output
            # Looking for "Total: XXXms" in the output
            for line in result.stdout.split('\n'):
                if 'Total:' in line and 'ms' in line:
                    # Extract the number before 'ms'
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if 'ms' in part:
                            # Remove 'ms' and convert to float
                            time_str = part.replace('ms', '').strip()
                            try:
                                runtime_ms = float(time_str)
                                return runtime_ms, True
                            except ValueError:
                                pass
            # If we couldn't parse the time, return failure
            return None, False
        else:
            # Print error output for debugging
            if result.stderr:
                print(f"Codon error:\n{result.stderr}", file=sys.stderr)
            if result.stdout:
                print(f"Codon output:\n{result.stdout}", file=sys.stderr)
            return None, False
            
    except subprocess.TimeoutExpired:
        return None, False
    except Exception as e:
        return None, False
    finally:
        # Clean up the test runner file
        if os.path.exists(test_runner_path):
            os.remove(test_runner_path)


def main():
    """Main evaluation function"""
    # Run Python tests
    python_time, python_success = run_python_tests()
    
    # Run Codon tests
    codon_time, codon_success = run_codon_tests()
    
    # Print results in required format
    print("Language    Runtime")
    print("-------------------")
    if python_success and python_time is not None:
        print(f"python      {int(python_time)}ms")
    else:
        print("python      FAILED")
    
    if codon_success and codon_time is not None:
        print(f"codon       {int(codon_time)}ms")
    else:
        print("codon       FAILED")
    
    # Exit with error if any tests failed
    if not (python_success and codon_success):
        sys.exit(1)


if __name__ == "__main__":
    main()
