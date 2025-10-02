# Deliverable 2

## Summary
This project ported the Bio.motifs package from Biopython to Codon. The port includes four core modules (__init__, matrix, minimal, thresholds) and a unified test suite that runs identically in both Python and Codon. All 8 tests pass in both environments, demonstrating functional equivalence between the implementations.

Model Used: `Claude Sonnet 4.5`

## Repository Setup
```
week2/
├── code/
│   ├── bio_codon/
│   │   └── motifs/
│   │       ├── __init__.codon
│   │       ├── matrix.codon
│   │       ├── minimal.codon
│   │       └── thresholds.codon
│   └── bio_python/
│       └── motifs/
│           ├── __init__.py
│           ├── matrix.py
│           ├── minimal.py
│           └── thresholds.py
├── data/
│   └── motifs/
│       ├── minimal_test.meme
│       └── minimal_test_rna.meme
├── test/
│   └── motifs/
│       ├── test_motifs.py
│       ├── test_motifs_online.py
│       └── test_motifs_FULL.py
├── test.py         # Unified test file
├── report.md       # This report
└── ai.md          # AI usage documentation
```

## Python-to-Codon Conversion

### Key Changes Required:
1. **Import Strategy**: Used `from python import` for external dependencies (urllib, Bio.Align, numpy) while keeping local imports standard
2. **C Extension Replacement**: Replaced `_pwm.calculate()` C extension in matrix.py with pure Python/NumPy implementation
3. **Type Annotations**: Added explicit type hints to class attributes for Codon's static analysis
4. **Property Decorators**: Converted `property()` function calls to `@property` decorators
5. **Exception Handling**: Replaced specific exceptions like `ImportError` with generic `Exception`
6. **Dead Code Removal**: Removed unused imports like `warnings` and `numbers`

### Compilation Results:
All four .codon files compile successfully:
- `__init__.codon`: Compiles after removing dead `warnings` import
- `matrix.codon`: Compiles after replacing `numbers.Integral` with `int`
- `minimal.codon`: Compiles successfully
- `thresholds.codon`: Compiles successfully

## Test Framework

### Unified test.py
Created a single test file that works with both Python and Codon using environment detection:
```python
if hasattr(str, 'memcpy'):
    # Codon
    import code.bio_codon.motifs as motifs
    from python import Bio as _Bio
    from python import numpy as np
else:
    # Python
    from Bio import motifs
    from Bio.Seq import Seq
    import numpy as np
```

### Test Results:
```
Dataset: Python
----------------------------------------------------------------------
[PASS] test_format_jaspar
[PASS] test_format_transfac
[PASS] test_relative_entropy_basic
[PASS] test_relative_entropy_background
[PASS] test_reverse_complement_forward
[PASS] test_reverse_complement_reverse
[PASS] test_minimal_meme_parser
[PASS] test_meme_parser_rna
Tests passed: 8/8

Dataset: Codon
----------------------------------------------------------------------
[PASS] test_format_jaspar
[PASS] test_format_transfac
[PASS] test_relative_entropy_basic
[PASS] test_relative_entropy_background
[PASS] test_reverse_complement_forward
[PASS] test_reverse_complement_reverse
[PASS] test_minimal_meme_parser
[PASS] test_meme_parser_rna
Tests passed: 8/8
```

### Running Tests:
**Python:**
```bash
python3 test.py
```

**Codon:**
```bash
export CODON_PYTHON=/usr/lib/x86_64-linux-gnu/libpython3.13.so
export PYTHON_PATH=/home/ryan/.local/lib/python3.13/site-packages
codon run test.py
```

## Hiccups

1. **Environment Detection**: Initial attempts using `try-except` with `__codon__` failed because Codon parses the entire file at compile time. Solution: `hasattr(str, 'memcpy')` for runtime detection.

2. **Import Errors**: Codon couldn't find `code.bio_codon` module. Solution: Create `__init__.codon` files in parent directories.

3. **NumPy Import**: Had `import numpy as np` instead of `from python import numpy as np` in `__init__.codon`. Fixed by adding `from python import`.

4. **Data File Paths**: Tests expected files in current directory but they were in `./data/motifs/`. Updated all file paths.

5. **Multi-line Strings**: Codon's parser struggled with triple-quoted strings. Converted to single-line strings with `\n` escapes.

6. **Python Interop Setup**: Needed to set `CODON_PYTHON` and `PYTHON_PATH` environment variables for Codon to access Biopython.

## Key Learnings

1. **Codon Differences**: Codon requires explicit type hints, doesn't support `property()` function, and has different exception handling
2. **Python Interop**: `from python import` provides access to Python libraries but adds runtime overhead
3. **Environment Detection**: `hasattr(str, 'memcpy')` is the reliable way to detect Codon vs Python
4. **Single File Testing**: Possible to create unified test files that work in both environments with proper conditional imports
5. **Package Structure**: Codon needs `__init__.codon` files to recognize package hierarchies
