# Deliverable 3

## Summary
This project ported biotite's phylogenetic tree construction module from Cython to Codon. The port includes four core files implementing UPGMA and neighbor-joining clustering algorithms, along with tree data structures. All three required tests pass in both Python (biotite) and Codon implementations, producing identical topological results with comparable performance (~3-7ms runtime for both).

Model Used: `Claude Sonnet 4.5`
Time Spent: Approximately 5-6 hours


## Cython-to-Codon Conversion

### Source Analysis
The original biotite.sequence.phylo module is written in Cython (.pyx files) for performance. Key components:
- **tree.pyx**: Defines `TreeNode` and `Tree` classes for representing phylogenetic trees
- **upgma.pyx**: Implements UPGMA (Unweighted Pair Group Method with Arithmetic mean) clustering
- **nj.pyx**: Implements neighbor-joining clustering algorithm

### Key Conversion Challenges:

1. **Class Ordering**
   - **Issue**: Codon requires classes to be defined before they're referenced in type annotations
   - **Solution**: Reordered `TreeNode` before `Tree` class definition

2. **Type System Adjustments**
   - **Field Declarations**: Added explicit type annotations for all class fields
   - **Optional Parameters**: Removed `Optional[int]` from `__init__` parameters (type inference handles None defaults)
   - **List Types**: Changed `_children` from tuple to `List[TreeNode]` (Codon tuples require compile-time known sizes)

3. **Collection Type Handling**
   - **Initial Approach**: Used `List[Optional[TreeNode]]` for node tracking
   - **Issue**: Type mismatch errors when assigning `None` to list elements
   - **Solution**: Used `List[TreeNode]` with separate boolean `is_clustered` array to track active nodes

4. **Set Hashing Extension**
   - **Challenge**: Codon doesn't have `frozenset` for immutable set hashing
   - **Solution**: Implemented custom `__hash__()` method for set class using `@extend` decorator with XOR-based hash algorithm

5. **Equality Comparison**
   - **Initial Approach**: Used set comparison for children nodes
   - **Issue**: Codon's `!=` operator unsupported for TreeNode objects in sets
   - **Solution**: Simplified `__eq__` to use list iteration with identity checks (`is not`)
   - **Hash Function**: Changed to XOR-based hash using object identity (`id(child)`)

6. **Float Precision**
   - Used `np.float64` instead of `np.float32` throughout
   - Added explicit `float()` casts in arithmetic operations
   - Changed to `MAX_FLOAT = np.finfo(np.float64).max`

7. **Format String Issues**
   - **Challenge**: Dynamic precision formatting (`f"{value:.{precision}f}"`) caused locale errors
   - **Solution**: Replaced with `round()` function: `rounded_dist = round(self._distance, round_distance)`

8. **NumPy Data Types**
   - Boolean arrays: Changed from `np.uint8` to `dtype=bool` for `is_clustered` arrays
   - Integer arrays: Used `np.int32` for cluster sizes and indices
   - All NumPy operations use explicit dtype specifications

9. **Exception Handling**
   - Replaced custom `TreeError` exceptions with standard `ValueError`
   - Simplified exception hierarchy for Codon compatibility

10. **Data Loading Workaround**
    - **Issue**: Codon's NumPy implementation has parser bug with `loadtxt`
    - **Solution**: Used Python's numpy via `from python import numpy as pnp`
    ```python
    distances: np.ndarray[int,2] = pnp.loadtxt("data/distances.txt", dtype=pnp.int64)
    ```

## Test Implementation

### Test Suite
Implemented three required tests from biotite's test suite:
- `test_distances`: Verifies distance calculations and topological distances
- `test_upgma`: Compares UPGMA results with expected tree structure
- `test_neighbor_joining`: Validates neighbor-joining algorithm output

## Evaluation Script

### Runtime Measurement Challenge
**Initial Problem**: First implementation showed extremely high runtimes:
- Python: ~134ms
- Codon: ~8429ms

**Root Cause**: Timing the entire subprocess call included:
- Python/Codon interpreter startup overhead
- Module imports and initialization
- For Codon: JIT compilation overhead
- Actual test execution (only ~3-5ms)

**Solution**: Parse internal timing from test output instead of timing subprocess
- Test files already report "Total: XXXms" with accurate internal timing
- Modified both `run_python_tests()` and `run_codon_tests()` to parse this output
- Extract timing by searching for "Total:" line and parsing the millisecond value

### Evaluation Results

#### Run 1
```
Language    Runtime
-------------------
python      7ms
codon       4ms
```

#### Run 2
```
Language    Runtime
-------------------
python      3ms
codon       3ms
```

#### Run 3
```
Language    Runtime
-------------------
python      6ms
codon       5ms
```

### Performance Analysis
Both implementations achieve similar performance (~3-7ms range), with slight variations due to:
- System load and other processes
- CPU cache warmth and frequency scaling
- Test timing precision for short-duration tests

The Codon implementation successfully matches or slightly exceeds Python/biotite performance while maintaining identical functional behavior.

## GitHub Actions Integration

### CI/CD Setup Challenge
**Issue**: GitHub Actions runs failed with "codon FAILED" error despite local success

**Root Cause**: The `CODON_PYTHON` environment variable was set in one step but not persisted across subsequent steps. Week3's Codon code requires this variable for Python bridge functionality (`from python import numpy as pnp`).

**Solution**: Updated `.github/workflows/actions.yml` to persist the environment variable:
```yaml
- name: Set up Codon Python bridge
  run: |
    CODON_PYTHON=$(find_libpython)
    echo "CODON_PYTHON=${CODON_PYTHON}" >> $GITHUB_ENV
    echo "Found Python at: ${CODON_PYTHON}"
```

This makes `CODON_PYTHON` available to all subsequent workflow steps, allowing the Codon Python bridge to function correctly.

## Hiccups and Solutions

1. **Class Forward References**: Fixed by reordering class definitions
2. **Optional Type Handling**: Switched to boolean tracking arrays
3. **Set Operations**: Implemented custom hash extension
4. **Format Strings**: Replaced dynamic formatting with round() function
5. **Runtime Measurement**: Parsed internal timing instead of subprocess timing
6. **GitHub Actions**: Persisted CODON_PYTHON environment variable
7. **NumPy Loading**: Used Python bridge for loadtxt due to Codon parser bug

All challenges were resolved through iterative debugging and applying instructor-provided hints where appropriate.
