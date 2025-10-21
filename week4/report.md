# Deliverable 4

## Summary
This project implements four sequence alignment algorithms in both Python and Codon: global alignment (Needleman-Wunsch), local alignment (Smith-Waterman), semi-global/fitting alignment, and affine gap penalty alignment. The implementations use dynamic programming with NumPy arrays for memory efficiency to handle large biological sequences like mitochondrial genomes (16k+ base pairs).

Model Used: `Claude Sonnet 4.5`
Time Spent: ~4 hours

## Algorithm Implementations

### Global Alignment (Needleman-Wunsch)
- **Purpose**: End-to-end alignment of two sequences
- **Approach**: Standard dynamic programming with gap penalties
- **Scoring**: Match +3, Mismatch -3, Gap -2
- **Complexity**: O(mn) time and space

### Local Alignment (Smith-Waterman)
- **Purpose**: Find best local match between sequences
- **Key Difference**: Allows 0 as minimum score (no negative alignments)
- **Traceback**: Starts from maximum score position, stops at 0
- **Use Case**: Finding conserved regions in sequences

### Semi-Global Alignment (Fitting)
- **Purpose**: Align sequences with free end gaps
- **Initialization**: First row and column set to 0 (no gap penalty at edges)
- **Traceback**: Starts from maximum in last row or column
- **Use Case**: Finding one sequence within another

### Affine Gap Penalty Alignment
- **Purpose**: More realistic gap modeling with different open/extend costs
- **Matrices**: Uses 3 DP matrices (M for match, I for insertion, D for deletion)
- **Scoring**: Gap open -5, Gap extend -1
- **Complexity**: O(mn) time, 3x space overhead

## Implementation Details

### Memory Optimization
All implementations use NumPy arrays with `dtype=np.int32` for memory efficiency:
```python
dp = np.zeros((m + 1, n + 1), dtype=np.int32)
```

**Memory Comparison for MT sequences (16,569 × 16,499 matrix)**:
- Python lists: ~28 bytes per cell = ~7.6GB for the DP matrix alone
- NumPy int32: 4 bytes per cell = ~1.1GB for the DP matrix
- **Improvement**: 85% memory reduction, enabling alignment of large genomes

The choice of NumPy was critical - attempting to use Python lists would exceed typical memory allocations for WSL or standard laptops.

### Traceback Strategy
All algorithms build alignments backwards from the optimal endpoint:
1. Identify starting position (varies by algorithm):
   - Global: Always (m, n)
   - Local: Maximum score position in matrix
   - Semi-global: Maximum in last row or column
   - Affine: Maximum among M[m,n], I[m,n], D[m,n]
2. Follow DP matrix to determine path taken using score comparisons
3. Build aligned sequences in reverse by appending characters
4. Reverse final alignment strings before returning

**Key insight**: Using explicit score recalculation during traceback (e.g., `dp[i-1,j-1] + match`) ensures we follow the correct path even when multiple paths yield the same optimal score.

### Python Implementation
**Files**: `/week4/code/python/*.py`
- `global_alignment.py` - Needleman-Wunsch
- `local_alignment.py` - Smith-Waterman  
- `semiglobal_alignment.py` - Fitting alignment
- `affine_alignment.py` - Affine gap penalty
- `utils.py` - FASTA I/O utilities
- `__init__.py` - Package exports

### Codon Implementation
**Files**: `/week4/code/codon/*.codon`
- Same structure as Python implementation
- Uses type annotations for Codon compatibility
- Explicit int casting for NumPy array indexing
- Modified string formatting to avoid f-string locale issues

## Codon Porting Challenges

### Type Compatibility Issues
Codon's strict type system required explicit casting when accessing NumPy arrays:
```python
# Codon requires explicit int() casting
diag_score = int(dp[i-1, j-1]) + match
```

**Problem**: NumPy array elements in Codon have type `Int[32]`, which cannot be directly mixed with Python's native `int` type in arithmetic operations. This differs from Python where NumPy handles type conversions automatically.

**Solution**: Wrap all NumPy array accesses with `int()` to convert to native type before arithmetic:
```python
# Before (causes type error in Codon)
dp[i, j] = dp[i-1, j-1] + match

# After (works in Codon)
dp[i, j] = int(dp[i-1, j-1]) + match
```

### String Formatting
F-string format specifiers (`{var:<20}`) caused locale errors in Codon:
```
ValueError: invalid format specifier: locale::facet::_S_create_c_locale name not valid
```

**Problem**: Codon doesn't support advanced f-string format specifiers like `:<width` for padding and alignment. These specifiers rely on locale settings that aren't available in Codon's runtime.

**Solution**: Manual padding instead of format specifiers:
```python
# Before (causes error in Codon)
print(f"{method:<20} {language:<10} {runtime}ms")

# After (works in Codon)
method_padded = method + " " * (20 - len(method))
language_padded = language + " " * (10 - len(language))
print(f"{method_padded} {language_padded} {runtime}ms")
```

### List vs Tuple Handling
Codon requires consistent collection types. Changed alignment construction to use lists throughout:
```python
# Build as lists, then join
aligned1 = []
aligned2 = []
# ... append characters ...
aligned1.reverse()
aligned2.reverse()
aligned1_str = ''.join(aligned1)
aligned2_str = ''.join(aligned2)
```

**Rationale**: Using `.reverse()` on lists is more explicit than `reversed()` which returns an iterator. Codon handles in-place list operations more efficiently.

## Evaluation Script

### Design
The `evaluate.py` script tests all implementations:
1. Runs Python tests by importing from `code.python`
2. Generates and executes inline Codon test code
3. Parses outputs and combines results
4. Displays formatted comparison table

### Inline Codon Execution
Instead of maintaining a separate test file, `evaluate.py` creates temporary Codon code:
```python
def run_codon_tests():
    # Create inline Codon test script
    codon_test_code = '''...'''
    
    # Write to temp file and execute
    with tempfile.NamedTemporaryFile(mode='w', suffix='.codon', delete=False) as f:
        f.write(codon_test_code)
    
    # Run with subprocess
    result = subprocess.run([codon_path, "run", temp_codon_file], ...)
```

### Test Coverage
- **Small sequences**: q1.fa and t1.fa (5 sequence pairs each)
- **Large sequences**: MT-human.fa and MT-orang.fa (16,569 bp)
- **Total tests**: 24 per implementation (48 total)

## Performance Results

### Complete Test Results

```
Method               Language   Runtime
--------------------------------------
global-q1            python     1ms
local-q1             python     1ms
fitting-q1           python     1ms
affine-q1            python     1ms
global-q2            python     1ms
local-q2             python     1ms
fitting-q2           python     1ms
affine-q2            python     1ms
global-q3            python     1ms
local-q3             python     1ms
fitting-q3           python     1ms
affine-q3            python     2ms
global-q4            python     1ms
local-q4             python     1ms
fitting-q4           python     1ms
affine-q4            python     1ms
global-q5            python     1ms
local-q5             python     1ms
fitting-q5           python     1ms
affine-q5            python     1ms
global-mt_human      python     73507ms
local-mt_human       python     94188ms
fitting-mt_human     python     73482ms
affine-mt_human      python     218176ms
global-q1            codon      1ms
local-q1             codon      1ms
fitting-q1           codon      1ms
affine-q1            codon      1ms
global-q2            codon      1ms
local-q2             codon      1ms
fitting-q2           codon      1ms
affine-q2            codon      1ms
global-q3            codon      1ms
local-q3             codon      1ms
fitting-q3           codon      2ms
affine-q3            codon      3ms
global-q4            codon      1ms
local-q4             codon      1ms
fitting-q4           codon      1ms
affine-q4            codon      1ms
global-q5            codon      1ms
local-q5             codon      1ms
fitting-q5           codon      1ms
affine-q5            codon      2ms
global-mt_human      codon      53866ms
local-mt_human       codon      64838ms
fitting-mt_human     codon      49705ms
affine-mt_human      codon      158556ms

Total tests: 48 (Python: 24, Codon: 24)
```

### Python Results (NumPy int32)
- **Small sequences (q1-q5)**: 1-2ms per alignment
- **MT Global**: 73,507ms (~73 seconds)
- **MT Local**: 94,188ms (~94 seconds)
- **MT Fitting**: 73,482ms (~73 seconds)
- **MT Affine**: 218,176ms (~218 seconds)

### Codon Results
- **Small sequences (q1-q5)**: 1-3ms per alignment
- **MT Global**: 53,866ms (~54 seconds) - **27% faster than Python**
- **MT Local**: 64,838ms (~65 seconds) - **31% faster than Python**
- **MT Fitting**: 49,705ms (~50 seconds) - **32% faster than Python**
- **MT Affine**: 158,556ms (~159 seconds) - **27% faster than Python**

### Performance Analysis
Codon consistently outperforms Python by 27-32% on large sequences due to:
- Ahead-of-time compilation to native code
- No interpreter overhead
- Optimized memory access patterns
- Type specialization during compilation

## Key Debugging Steps

1. **Memory Optimization**: Selected NumPy arrays from the start to avoid memory issues with Python lists. For MT sequences (16,569 × 16,499 matrix), Python lists would use ~7.6GB while NumPy int32 arrays use only ~1.1GB.

2. **Traceback Logic**: Added comprehensive fallback conditions and safety breaks to prevent infinite loops during alignment reconstruction. Key insight was ensuring every branch in the while loop either decrements i, j, or breaks.

3. **Codon Type Errors**: Discovered need for explicit `int()` casting when accessing NumPy array elements in Codon. This was the primary blocker for porting and required systematic replacement throughout all algorithms.

4. **Format String Errors**: Replaced f-string format specifiers with manual string padding for Codon compatibility. Tested padding logic to ensure proper alignment of output columns.

5. **Test Integration**: Consolidated testing into single `evaluate.py` file using inline Codon code generation. This eliminated the need for separate test files and simplified maintenance - all test logic now lives in one place.

6. **Alignment Construction**: Changed from using `reversed()` iterator to in-place `.reverse()` method for better Codon compatibility and clarity.

## Technical Insights

### Algorithm Complexity
- All algorithms: O(mn) time complexity
- Space: O(mn) for standard implementations
- Affine: 3x space overhead due to three DP matrices

### Memory Efficiency Impact
For MT sequences (16,569 × 16,499 matrix):
- Matrix cells: 273,476,931
- Python lists: 7.6GB (28 bytes/cell)
- NumPy int32: 1.1GB (4 bytes/cell)
- Peak usage with overhead: ~4.6GB (Python), 5.8-8.0GB (Codon)

**Note**: Codon uses more memory than Python despite being faster, likely due to compilation overhead and LLVM runtime. Affine alignment uses the most memory (up to 8.0GB in Codon) due to three DP matrices.

### Codon vs Python
- **Compilation**: Codon uses LLVM for native code generation
- **Type System**: Stricter than Python, requires explicit types
- **Performance**: 27-32% faster on compute-intensive operations
- **Memory**: Uses 70% more memory (5.8-8.0GB vs 4.6GB) but delivers faster execution
- **Development**: Requires more type annotations and casting

## Conclusion

Both Python and Codon implementations successfully handle small and large sequence alignments. The NumPy-based approach provides excellent memory efficiency, while Codon's compilation model delivers superior runtime performance. The key to successful porting was understanding Codon's type system requirements and working within its constraints for string formatting and type compatibility.
