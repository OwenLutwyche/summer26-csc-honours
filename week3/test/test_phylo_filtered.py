"""
Filtered phylogenetic tree tests for biotite.sequence.phylo
Only includes: test_distances, test_upgma, test_neighbor_joining
"""

import numpy as np
import time
import os
import biotite.sequence.phylo as phylo


def get_distances():
    """
    Load test distances matrix from file.
    Based on the example "Dendrogram of the BLOSUM62 matrix"
    with the small modification M[i,j] += i+j
    """
    # Data folder is beside test folder, not inside it
    data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
    distances = np.loadtxt(os.path.join(data_dir, "distances.txt"), dtype=float)
    return distances


def get_upgma_newick():
    """
    Load the Newick notation for the UPGMA tree from file.
    """
    # Data folder is beside test folder, not inside it
    data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
    with open(os.path.join(data_dir, "newick_upgma.txt"), "r") as file:
        newick = file.read().strip()
    return newick


def get_tree():
    """Create a phylogenetic tree using UPGMA"""
    distances = get_distances()
    return phylo.upgma(distances)


def test_distances():
    """
    Test whether the `distance_to()` and topological distances work correctly.
    """
    tree = get_tree()
    
    # Tree is created via UPGMA
    # -> The distances to root should be equal for all leaf nodes
    dist = tree.root.distance_to(tree.leaves[0])
    for leaf in tree.leaves:
        assert leaf.distance_to(tree.root) == dist
    
    # Example topological distances
    assert tree.get_distance(0, 19, True) == 9
    assert tree.get_distance(4, 2, True) == 10
    
    print("test_distances passed")


def test_upgma():
    """
    Compare the results of `upgma()` with expected tree structure.
    """
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
    
    print("test_upgma passed")


def test_neighbor_joining():
    """
    Compare the results of `neighbor_join()` with a known tree.
    """
    dist = np.array([
        [ 0,  5,  4,  7,  6,  8],
        [ 5,  0,  7, 10,  9, 11],
        [ 4,  7,  0,  7,  6,  8],
        [ 7, 10,  7,  0,  5,  9],
        [ 6,  9,  6,  5,  0,  8],
        [ 8, 11,  8,  9,  8,  0],
    ])  # fmt: skip

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
    
    print("test_neighbor_joining passed")


def run_tests():
    """Run all tests and report timing"""
    print("Running Python phylo tests...")
    print("=" * 40)
    
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
    
    print("=" * 40)
    print(f"test_distances:        {time_distances:.2f}ms")
    print(f"test_upgma:            {time_upgma:.2f}ms")
    print(f"test_neighbor_joining: {time_neighbor:.2f}ms")
    print(f"Total:                 {total_time:.2f}ms")
    print("=" * 40)
    
    return total_time


if __name__ == "__main__":
    try:
        total_time = run_tests()
        print(f"\nAll tests passed! Total time: {total_time:.2f}ms")
    except Exception as e:
        print(f"\nTest failed: {e}")
        import traceback
        traceback.print_exc()
