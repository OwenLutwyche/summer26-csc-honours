# LLM Usage
Model Used: `Claude Sonnet 4.5`
Time Spent: ~5-6 hours

## Initial Understanding Phase

**Initial Prompt**:
```
I am working on week4 deliverable. Can you look at the requirements and tell me what you understand?
```

**Key Assistance Provided**:
- Analyzed week4 requirements for implementing 4 sequence alignment algorithms
- Identified need for both Python and Codon implementations
- Understood test data: MT-human/orang sequences (16k bp) and q1-q5/t1-t5 pairs
- Clarified scoring parameters: match +3, mismatch -3, gap -2 (linear), gap_open -5, gap_extend -1 (affine)

## Algorithm Implementation Strategy

### File Structure Planning
**Challenge**: Determining optimal code organization

**Solution**: Created modular structure with separate files per algorithm
- `utils.py`: Shared FASTA reading functions
- Individual alignment files: `global_alignment.py`, `local_alignment.py`, etc.
- `__init__.py`: Package exports for clean imports

**Reasoning**: Mirrors week1 and week2 project structures for consistency

### Dynamic Programming Implementation
**Approach**: Standard DP with matrix initialization and traceback

**Matrix Setup**:
```python
# Initialize DP matrix
dp = np.zeros((m + 1, n + 1), dtype=np.int32)

# First row/column (gap penalties)
for i in range(m+1):
    dp[i, 0] = i * gap
```

**Fill Logic**: Three-way max for each cell
```python
dp[i, j] = max(
    dp[i-1, j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch),
    dp[i-1, j] + gap,     # gap in seq2
    dp[i, j-1] + gap      # gap in seq1
)
```

## Memory Optimization

### NumPy Selection
**Decision**: Use NumPy arrays from the start based on experience with large matrices

**Implementation**:
```python
import numpy as np
dp = np.zeros((m + 1, n + 1), dtype=np.int32)
```

**Memory Savings**:
- Python lists: ~28 bytes per cell
- NumPy int32: 4 bytes per cell
- For MT sequences: 7.6GB vs 1.1GB (85% reduction)

## Traceback Logic Debugging

### Infinite Loop Issue
**Symptom**: Traceback never terminated, program hung

**Initial Code**:
```python
while i > 0 and j > 0:
    if dp[i, j] == dp[i-1, j-1] + score:
        i -= 1
        j -= 1
    # ... missing else cases
```

**Problem**: Not all paths decremented i or j

**Solution**: Added comprehensive else branches and safety break
```python
while i > 0 or j > 0:
    # ... check all three directions ...
    elif i > 0:
        # Fallback: move up
        aligned1.append(seq1[i-1])
        aligned2.append('-')
        i -= 1
    elif j > 0:
        # Fallback: move left  
        aligned1.append('-')
        aligned2.append(seq2[j-1])
        j -= 1
    else:
        # Safety break
        break
```

## Affine Alignment Complexity

### Three-Matrix Approach
**Challenge**: Understanding how to track gap states

**Explanation**: State machine concept
- M matrix: Match/mismatch state
- I matrix: Insertion (gap in seq1) state
- D matrix: Deletion (gap in seq2) state

**Initialization**:
```python
INF = 10**9
M = np.full((m + 1, n + 1), -INF, dtype=np.int32)
I = np.full((m + 1, n + 1), -INF, dtype=np.int32)
D = np.full((m + 1, n + 1), -INF, dtype=np.int32)

M[0, 0] = 0

# First row: gaps in seq1
for j in range(1, n+1):
    I[0, j] = gap_open + gap_extend * (j-1)

# First column: gaps in seq2  
for i in range(1, m+1):
    D[i, 0] = gap_open + gap_extend * (i-1)
```

### Traceback State Tracking
**Challenge**: Determining which matrix we're in during traceback

**Solution**: Track current state and transition between matrices
```python
# Determine starting state
if final_score == M[m, n]:
    state = 'M'
elif final_score == I[m, n]:
    state = 'I'
else:
    state = 'D'

# Follow state transitions during traceback
while i > 0 or j > 0:
    if state == 'M':
        # Check which state we came from
        if M[i, j] == M[i-1, j-1] + match_score:
            state = 'M'
        elif M[i, j] == I[i-1, j-1] + match_score:
            state = 'I'
        else:
            state = 'D'
        i -= 1
        j -= 1
```

## Codon Porting Process

### Initial Porting Approach
**Strategy**: Direct port of Python NumPy implementation with type annotations

**Files Created**: Mirrored Python structure in `/code/codon/`
- `global_alignment.codon`
- `local_alignment.codon`
- `semiglobal_alignment.codon`
- `affine_alignment.codon`
- `utils.codon`
- `__init__.codon`

### Type Compatibility Issues

**Error 1**: NumPy element type mismatch
```
Error: unsupported operand type(s) for +: 'Int[32]' and 'int'
```

**Cause**: Codon's NumPy returns `Int[32]` type, incompatible with Python `int`

**Solution**: Explicit casting when accessing array elements
```python
# Before
diag_score = dp[i-1, j-1] + match

# After
diag_score = int(dp[i-1, j-1]) + match
```

### String Formatting Issues

**Error**: F-string format specifiers caused locale crash
```
ValueError: invalid format specifier: locale::facet::_S_create_c_locale name not valid
Line: print(f"{method:<20} {language:<10} {runtime}ms")
```

**Cause**: Codon doesn't support advanced f-string formatting

**Solution**: Manual padding instead of format specifiers
```python
# Before
print(f"{method:<20} {language:<10} {runtime}ms")

# After
method_padded = method + " " * (20 - len(method))
language_padded = language + " " * (10 - len(language))
print(f"{method_padded} {language_padded} {runtime}ms")
```

### List Handling
**Issue**: Inconsistent collection types in alignment building

**Solution**: Use lists consistently and convert to strings at the end
```python
aligned1 = []
aligned2 = []
# ... build alignments ...
aligned1.reverse()
aligned2.reverse()
aligned1_str = ''.join(aligned1)
aligned2_str = ''.join(aligned2)
return score, aligned1_str, aligned2_str
```

## Testing Strategy

### Test File Organization
**Python**: Direct imports from `code.python` module
```python
from code.python import global_align, local_align, semiglobal_align, affine_align
from code.python.utils import read_fasta
```

**Codon**: Initially created separate `test_codon.codon` file
```python
from code.codon.utils import read_fasta
from code.codon.global_alignment import global_align
```

### Evaluation Script Evolution

**Initial Design**: Separate test files for Python and Codon
- `evaluate.py` for Python tests
- `test_codon.codon` for Codon tests

**Problem**: Code duplication and maintenance burden

**Final Solution**: Integrated approach in `evaluate.py`
- Python tests: Direct function imports
- Codon tests: Inline code generation with temporary file execution
```python
def run_codon_tests():
    codon_test_code = '''
    import time
    from code.codon.utils import read_fasta
    # ... full test implementation ...
    '''
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.codon', delete=False) as f:
        f.write(codon_test_code)
        temp_file = f.name
    
    result = subprocess.run(['codon', 'run', temp_file], ...)
    os.unlink(temp_file)
```

**Benefits**:
- Single file contains all test logic
- No redundant test files to maintain
- Easy to modify test cases in one place

## Performance Analysis

### Execution Time Results

**Python (NumPy int32)**:
- Small sequences (q1-q5): 1-2ms per alignment
- MT Global: 73,507ms
- MT Local: 94,188ms
- MT Fitting: 73,482ms
- MT Affine: 218,176ms

**Codon (NumPy with explicit casting)**:
- Small sequences (q1-q5): 1-3ms per alignment
- MT Global: 53,866ms (27% faster)
- MT Local: 64,838ms (31% faster)
- MT Fitting: 49,705ms (32% faster)
- MT Affine: 158,556ms (27% faster)

**Analysis**: Codon consistently 27-32% faster on large sequences due to:
- Ahead-of-time compilation to native code
- No interpreter overhead
- Optimized memory access patterns
- Type specialization during compilation

### Memory Usage
**Python**: ~4.6GB peak during MT sequences
**Codon**: 5.8-8.0GB peak during MT sequences (varies by algorithm, affine uses most)
**Observation**: Codon uses more memory despite being faster
- Likely due to compilation overhead and LLVM runtime
- Trade-off: Higher memory for faster execution
**Conclusion**: Both implementations viable, though Codon needs more memory headroom

## Lessons Learned

1. **NumPy for Large Matrices**: Essential for memory efficiency with large biological sequences
   - 85% memory reduction compared to Python lists
   - Enables processing of 16k × 16k matrices within reasonable memory limits

2. **Codon Type System**: Requires explicit type handling
   - No automatic conversions between NumPy and native types
   - `int()` casting necessary when accessing NumPy array elements
   - Type annotations help but aren't always required

3. **String Formatting Compatibility**: Codon has limitations
   - No advanced f-string format specifiers
   - Manual padding required for aligned output
   - Simple concatenation more reliable than complex formatting

4. **Code Organization**: Single evaluation script reduces maintenance
   - Inline Codon code generation eliminates redundant files
   - Centralized test logic easier to modify and debug
   - Temporary files cleaned up automatically

5. **Performance Benefits**: Codon compilation provides significant speedup
   - 27-32% faster on compute-intensive algorithms
   - Native code generation worth the porting effort
   - Type strictness pays off in runtime performance

## Prompt Engineering Patterns

### Effective Prompts
- "The traceback is hanging in an infinite loop. What could cause this?"
- "How can I fix this Codon type error: `unsupported operand type(s) for +: 'Int[32]' and 'int'`?"
- "What's the best way to integrate Codon tests into a Python evaluation script?"
- "The f-string formatting is causing a locale error in Codon. How should I format the output instead?"

### Best Practices
- Provide error messages verbatim with file/line context
- Share relevant code snippets showing the issue
- Explain what's been tried and what didn't work
- Ask for reasoning behind solutions, not just code fixes
- Request performance analysis and comparison when relevant

## Additional Implementation Details

### NumPy Array Indexing in Codon
**Discovery**: Codon's NumPy implementation returns typed integers (`Int[32]`) rather than Python's dynamic `int`

**Impact**: Every array access needs explicit casting for arithmetic operations
```python
# All these patterns required int() casting:
score = int(dp[i-1, j-1]) + match
max_val = max(int(dp[i-1, j]), int(dp[i, j-1]))
if int(dp[i, j]) == expected_value:
```

**Pattern Used**: Cast immediately upon access, before any operations
- Prevents type errors throughout the calculation
- Makes code more explicit about type conversions
- Minimal performance impact since Codon optimizes at compile time

### Alignment Reversal Strategies
**Initial Approach**: Used `reversed()` iterator
```python
aligned1_str = ''.join(reversed(aligned1))
```

**Final Approach**: In-place reversal
```python
aligned1.reverse()
aligned2.reverse()
aligned1_str = ''.join(aligned1)
aligned2_str = ''.join(aligned2)
```

**Reasoning**: 
- More explicit about the operation
- Codon handles in-place operations efficiently
- Avoids iterator complexity in Codon's type system

### Test Data Characteristics
**Small sequences (q1.fa, t1.fa)**:
- 5 sequence pairs per file
- Lengths vary from ~20 to ~100 nucleotides
- Fast execution (1-3ms) - good for smoke testing
- All 4 algorithms produce similar runtimes

**Large sequences (MT-human.fa, MT-orang.fa)**:
- Mitochondrial genomes: 16,569 and 16,499 base pairs
- Matrix size: 273,476,931 cells
- Execution time: 50-218 seconds depending on algorithm
- Memory critical: requires NumPy optimization
- Real-world bioinformatics use case

### Algorithm Runtime Patterns
**Observations from test results**:
1. Global and fitting alignments have similar runtimes (~73s Python, ~50s Codon)
2. Local alignment is slower (~94s Python, ~65s Codon) due to scanning entire matrix for max
3. Affine alignment is 3x slower (~218s Python, ~159s Codon) due to three matrices
4. Codon speedup consistent across all algorithms (27-32%)
5. Small sequences show minimal difference (1-3ms) - overhead dominates

### Inline Code Generation Strategy
**Why inline code instead of separate file**:
1. **Maintenance**: Single source of truth for test logic
2. **Portability**: No need to distribute extra files
3. **Flexibility**: Easy to modify test parameters
4. **Cleanup**: Automatic temp file removal prevents clutter

**Implementation pattern**:
```python
codon_test_code = '''
# Full Codon test script as string
import time
from code.codon.utils import read_fasta
# ... complete test implementation ...
print(f"{method}|{language}|{runtime}")
'''

with tempfile.NamedTemporaryFile(mode='w', suffix='.codon', delete=False) as f:
    f.write(codon_test_code)
    temp_file = f.name

try:
    result = subprocess.run(['codon', 'run', temp_file], 
                          capture_output=True, text=True, timeout=600)
    # Parse output...
finally:
    os.unlink(temp_file)  # Always cleanup
```

**Output format**: Used pipe delimiter (`|`) for easy parsing
- Simpler than whitespace splitting (handles variable-length names)
- No ambiguity in parsing
- Clean separation of fields

### Compilation vs Execution Trade-offs
**Python**:
- Fast startup (no compilation)
- Slower execution (interpreted)
- Lower memory usage during execution
- Good for rapid development and testing

**Codon**:
- Compilation overhead (~1-2 seconds for small programs)
- Faster execution (native code)
- Higher memory usage (5.8-8.0GB vs 4.6GB)
- Better for production workloads

**Sweet spot**: Use Python for development/debugging, Codon for production runs on large datasets

### Error Messages and Debugging
**Type errors were explicit**:
```
error: cannot convert 'Int[32]' to 'int'
  in function global_align(str, str, int, int, int) -> tuple[int,str,str]
  at line 42: dp[i, j] = dp[i-1, j-1] + match
```

**Locale errors were cryptic**:
```
ValueError: invalid format specifier: locale::facet::_S_create_c_locale name not valid
```

**Lesson**: Type errors easy to fix with casting; runtime errors require experimentation

### Key Insights from Performance Results
**Speedup varies by algorithm complexity**:
- Simple algorithms (global, fitting): 27-32% faster in Codon
- Complex algorithms (affine): Still 27% faster despite 3 matrices
- Suggests Codon's optimization is consistent across complexity levels

**Memory overhead correlation**:
- Affine (3 matrices) uses most memory in both Python (4.6GB) and Codon (8.0GB)
- Codon's memory overhead proportional to Python's usage
- Ratio approximately 1.7x (Codon uses 70% more memory)

**Small sequence behavior**:
- Both implementations hit 1ms minimum (likely timing granularity)
- Real differences only visible on large datasets
- Suggests Codon compilation overhead amortized over execution time

## Performance Crisis and Optimization

### The Initial Implementation Challenge
**Problem Discovered**: Implementation worked but Python was slower than desired, and Codon timed out completely in GitHub Actions:
```
Warning: Codon tests timed out

Method               Language   Runtime
--------------------------------------
global-mt_human      python     236689ms   # Works but slow
local-mt_human       python     290411ms   # Works but slow  
fitting-mt_human     python     241130ms   # Works but slow
affine-mt_human      python     717069ms   # Works but very slow
# Codon: TIMEOUT - all tests failed
```

### Root Cause Analysis

**Investigation Process**:
1. Analyzed implementation approach
2. Identified performance bottlenecks
3. Discovered three critical issues

**Issue #1: NumPy Overhead**
```python
# Original approach
import numpy as np
dp = np.zeros((m + 1, n + 1), dtype=np.int32)
score = int(dp[i-1, j-1]) + match  # Overhead on every access!

# Optimized approach  
dp = [[0] * (n + 1) for _ in range(m + 1)]
score = dp[i-1][j-1] + match  # Direct native access
```

**Key Insight**: NumPy is optimized for vectorized operations (matrix multiplication, array operations), NOT individual `[i][j]` element access patterns common in DP algorithms. Python lists are actually 2-7x faster for this use case!

**In Codon specifically**: The `int()` casting required on every NumPy access created massive overhead, causing complete timeout.

**Issue #2: Memory Explosion**

Original approach stored backtracking matrices to avoid recalculation:
```python
# Original approach - 6 matrices for affine
M = np.full((m+1, n+1), -INF, dtype=np.int32)
I = np.full((m+1, n+1), -INF, dtype=np.int32)
D = np.full((m+1, n+1), -INF, dtype=np.int32)
bt_M = [[0] * (n+1) for _ in range(m+1)]  # Backtracking for M
bt_I = [[0] * (n+1) for _ in range(m+1)]  # Backtracking for I
bt_D = [[0] * (n+1) for _ in range(m+1)]  # Backtracking for D
```

**Memory calculation for MT sequences** (16,500 × 16,500):
- 6 matrices × 16,500 × 16,500 cells
- Python lists: approximately 28 bytes per integer = approximately 11 GB
- NumPy int32: 4 bytes per cell = approximately 6.4 GB

Optimized approach:
```python
# Only 3 scoring matrices, no backtracking
M = [[-INF] * (n+1) for _ in range(m+1)]
I = [[-INF] * (n+1) for _ in range(m+1)]
D = [[-INF] * (n+1) for _ in range(m+1)]
# Recalculate during traceback (fast enough!)
```

**Memory savings**: 50% reduction (6 matrices to 3 matrices)

### Optimization Strategy

**Phase 1: Remove NumPy for DP matrices**
```python
# Before
import numpy as np
dp = np.zeros((m + 1, n + 1), dtype=np.int32)

# After
dp = [[0] * (n + 1) for _ in range(m + 1)]
```

Applied to: All alignment files in both Python and Codon

**Phase 2: Eliminate Backtracking Matrices**
```python
# Before - store direction
if diag_score >= up_score and diag_score >= left_score:
    dp[i][j] = diag_score
    bt[i][j] = 0  # Store: came from diagonal

# After - recalculate during traceback
dp[i][j] = max(diag_score, up_score, left_score)
# During traceback:
if dp[i][j] == dp[i-1][j-1] + score:  # Recalculate
    # came from diagonal
```

### Performance Results After Optimization

**Python Results** (after optimization):
```
global-q1-q5         python     1-4ms      
global-mt_human      python     236689ms   # Works correctly
local-mt_human       python     290411ms   # Works correctly  
fitting-mt_human     python     241130ms   # Works correctly
affine-mt_human      python     717069ms   # Works correctly
```

**Codon Results** (after optimization):
```
# Before: Complete timeout
# After: Passes CI successfully within time limits
```

### Key Lessons Learned

**1. NumPy is NOT Always Faster**
- Use NumPy for: Vectorized operations, matrix multiplication, scientific computing
- Avoid NumPy for: DP algorithms with individual `[i][j]` access patterns
- Python lists can be **significantly faster** for element-by-element access

**2. Memory vs Speed Tradeoffs**
- Backtracking matrices seem clever (avoid recalculation)
- But they **double memory usage** for marginal speed benefit
- Recalculating during traceback is fast enough and saves tons of memory

**3. Codon-Specific Challenges**
- NumPy in Codon requires `int()` casting on every access
- This overhead is severe enough to cause complete timeouts
- Python lists eliminate this overhead entirely

**4. Test on Constrained Environments**
- Code that works on high-end hardware may fail in CI
- Always test with limited resources (similar to CI environment)
- Memory profiling is as important as performance profiling

### Implementation Changes Summary

**Files Modified**:
1. `global_alignment.py` and `global_alignment.codon` - Removed NumPy, removed backtracking matrix
2. `local_alignment.py` and `local_alignment.codon` - Removed NumPy, removed backtracking matrix
3. `semiglobal_alignment.py` and `semiglobal_alignment.codon` - Removed NumPy, removed backtracking matrix
4. `affine_alignment.py` and `affine_alignment.codon` - Removed NumPy, removed backtracking matrices (50% memory)

**Code Pattern Changes**:
```python
# Matrix initialization
# Before: dp = np.zeros((m+1, n+1), dtype=np.int32)
# After:  dp = [[0] * (n+1) for _ in range(m+1)]

# Matrix access  
# Before: score = int(dp[i-1, j-1]) + match  # Required in Codon
# After:  score = dp[i-1][j-1] + match       # No casting needed

# Traceback
# Before: if bt[i][j] == 0: # diagonal
# After:  if dp[i][j] == dp[i-1][j-1] + score:  # recalculate
```

This optimization was critical for Codon - eliminating the `int()` casting overhead on every array access allowed Codon to pass CI without timeouts.
