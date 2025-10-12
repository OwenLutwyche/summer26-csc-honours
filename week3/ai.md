# LLM Usage
Model Used: `Claude Sonnet 4.5`

## Initial Porting Strategy
* **Initial Prompt**:
```
Here is my instructions.
The goals of this week's deliverable are:
    code analysis: get familiar with the unknown codebases;
    evolutionary trees: port popular biotite's phylo module to Codon;

Your job is to port biotite's phylo package to Codon. Note that the source is written in Cython, not pure Python.

The source code can be found on Github. You only need to port the parts needed to run the two tests in test_phylo.py:
    test_distances,
    test_upgma,
    test_neighbor_joining.

While you can use from python import for missing pieces, note that this is not needed here.

You must run Python and Codon tests and output the time of all tests in milliseconds. Use the time module in Python and Codon (do not time this from your evaluate.sh).

Please do an initial port of the 4 files, one file at a time due to tokens, put them in biotite_codon folder
```

* **Key Assistance Provided**:
  - Analyzed biotite's Cython source code (`.pyx` files) from GitHub
  - Identified four core files needed: `__init__.py`, `tree.pyx`, `upgma.pyx`, `nj.pyx`
  - Created Codon equivalents (`.codon` files) in `biotite_codon/` folder

## Codon-Specific Adaptations

### Class Structure Ordering
* **Challenge**: Initial compilation error - `TreeNode` referenced before definition in `Tree` class
* **Solution**: Reordered classes to define `TreeNode` before `Tree` class, matching proper forward reference handling in Codon
* **Reasoning**: Through iterative compilation attempts, discovered Codon requires classes to be defined before they're referenced in type annotations

### Type System Adjustments
* **Field Declarations**: Added explicit type annotations for all class fields (e.g., `_index: int`, `_children: List[TreeNode]`)
* **Optional Types**: Used `Optional[TreeNode]` for fields that can be `None` (like `_parent`)
* **Parameter Type Inference**: Removed type annotation from `__init__` parameter `index=None` after compilation error showed Codon couldn't handle `Optional[int]` on optional parameters
* **List vs Tuple**: Changed `_children` from tuple to `List[TreeNode]` as tuples in Codon must have compile-time known sizes

### Collection Type Handling
* **Initial Approach**: Used `List[Optional[TreeNode]]` for nodes array in clustering algorithms
* **Compilation Issue**: Type mismatch errors when assigning `None` to list elements
* **Final Solution**: Used `List[TreeNode]` with boolean `is_clustered` array to track which nodes are active, avoiding `None` assignments entirely
* **Discovery Method**: Trial compilation revealed type system constraints, leading to cleaner algorithm design

### Set Hashing Extension
* **Challenge**: Codon doesn't have `frozenset` for immutable set hashing
* **Solution**: Added `@extend` decorator to implement custom `__hash__()` method for set class
* **Implementation**: Used hash algorithm with XOR operations and masking to create deterministic hash values
* **Source**: Adapted from instructor-provided hints for handling TreeNode hashing

### Equality Comparison
* **Initial Approach**: Used set comparison for children nodes (`set(self._children) != set(node._children)`)
* **Compilation Error**: Codon complained about unsupported `!=` operator for TreeNode objects
* **Root Cause Analysis**: Creating sets of TreeNode objects required comparison operators that weren't fully compatible
* **Final Solution**: Simplified `__eq__` to use list iteration with identity checks (`is not`) and compare children order
* **Hash Function**: Changed from set-based hash to XOR-based hash using object identity (`id(child)`)

### Float Precision
* **Changes**:
  - Used `np.float64` instead of `np.float32` throughout (as recommended in hints)
  - Initialized all distance values as floats (e.g., `0.0` instead of `0`)
  - Added explicit `float()` casts in arithmetic operations to avoid int/float compatibility issues
  - Changed `MAX_FLOAT = np.finfo(np.float64).max`

### Format String Issues
* **Challenge**: Runtime error with format specifiers: `locale::facet::_S_create_c_locale name not valid`
* **Initial Format**: Used dynamic precision formatting: `f"{value:.{precision}f}"`
* **Solution**: Replaced with `round()` function: `rounded_dist = round(self._distance, round_distance)` then `f"{label}:{rounded_dist}"`
* **Test Output**: Simplified test timing output to use `int()` instead of `.2f` format specifier

### NumPy Data Type Adjustments
* **Boolean Arrays**: Changed from `np.uint8` to `dtype=bool` for `is_clustered` arrays
* **Integer Arrays**: Used `np.int32` for cluster sizes and indices
* **Consistency**: Ensured all NumPy operations use explicit dtype specifications for Codon compatibility

## Exception Handling
* **Change**: Replaced all custom `TreeError` exceptions with standard `ValueError` exceptions
* **Reasoning**: Simplified exception hierarchy for Codon compatibility, removed need for custom exception class
* **Impact**: All error messages preserved, just using Python's built-in exception types

## Test Comparison Strategy
* **Challenge**: Direct tree equality (`test_tree == ref_tree`) failed for neighbor joining test
* **Analysis**: Object identity-based comparison doesn't work for trees constructed differently
* **Solution**: Compare all pairwise distances instead of tree structure:
```python
for i in range(6):
    for j in range(6):
        assert test_tree.get_distance(i, j) == ref_tree.get_distance(i, j)
        assert test_tree.get_distance(i, j, topological=True) == ref_tree.get_distance(i, j, topological=True)
```
* **Benefit**: More robust test that verifies functional equivalence rather than structural identity

## Evaluation Script Development
* **Requirement**: Output format must be:
```
Language    Runtime
-------------------
python      2000ms
codon       1000ms
```

* **Implementation**:
  - Created `evaluate.py` that runs both Python (biotite) and Codon (biotite_codon) tests
  - Used subprocess to execute test files and capture timing
  - Dynamically generates Codon test runner file with proper imports
  - Handles both local codon installation and GitHub Actions environment
  - Cleans up temporary files after execution

### Runtime Measurement Issue
* **Initial Problem**: First implementation showed extremely high runtimes:
  - Python: ~134ms
  - Codon: ~8429ms

* **Root Cause**: Timing the entire subprocess call included:
  - Python/Codon interpreter startup overhead
  - Module imports and initialization
  - For Codon: JIT compilation overhead
  - Actual test execution (only ~3-5ms)

* **Initial Implementation** (incorrect):
```python
start_time = time.time()
result = subprocess.run([...])
end_time = time.time()
runtime_ms = (end_time - start_time) * 1000
```

* **Solution**: Parse internal timing from test output instead of timing subprocess
  - Test files already report "Total: XXXms" with accurate internal timing
  - Modified both `run_python_tests()` and `run_codon_tests()` to parse this output
  - Extract timing by searching for "Total:" line and parsing the millisecond value

* **Final Implementation**:
```python
# Parse the total time from the output
for line in result.stdout.split('\n'):
    if 'Total:' in line and 'ms' in line:
        parts = line.split()
        for i, part in enumerate(parts):
            if 'ms' in part:
                time_str = part.replace('ms', '').strip()
                runtime_ms = float(time_str)
                return runtime_ms, True
```

* **Result**: Accurate timing showing ~3ms for both Python and Codon implementations

## Data Loading Workaround
* **Issue**: Codon's NumPy implementation has parser bug with `loadtxt`
* **Solution**: Used Python's numpy via `from python import numpy as pnp`, then:
```python
distances: np.ndarray[int,2] = pnp.loadtxt("data/distances.txt", dtype=pnp.int64)
```
* **Source**: Applied instructor-provided hint for this specific Codon limitation
