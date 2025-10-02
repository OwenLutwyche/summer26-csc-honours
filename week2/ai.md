# LLM Usage
Model Used: `Claude Sonnet 4.5`

## Model Used
* **Model**: Claude 4.5 Sonnet
* **Purpose**: Assisted with Python-to-Codon conversion and test framework porting

## Python-to-Codon Conversion
* **File Context**: Provided all 4 Python files from Biopython (`__init__.py`, `matrix.py`, `minimal.py`, `thresholds.py`)
* **Initial Prompt**:
```
I am porting Bio.motifs package to Codon. Specifically only the following files from this repository: 
github.com/biopython/biopython Bio.motifs (init.py) Bio.motifs.matrix Bio.motifs.minimal Bio.motifs.thresholds

I have the Python source files. If I want to convert them to Codon, what are the key differences I need 
to handle? What imports need to change? Are there any C extensions or Python-specific features I need to replace?
```

* **Key Assistance Provided**:
  - Identified C extension (`_pwm.calculate`) that needed pure Python replacement
  - Recommended `from python import` for external dependencies (urllib, Bio.Align, numpy)
  - Added explicit type hints to class attributes
  - Module-by-module conversion strategy

## Compilation Debugging
* **Context**: All 4 .codon files after initial conversion
* **Issues Identified**:
  - `__init__.codon`: Dead `warnings` import (removed)
  - `matrix.codon`: Missing `numbers` module (replaced `numbers.Integral` with `int`)
  - Property syntax: `property()` function not supported in Codon
  - Exception types: `ImportError` not available (used generic `Exception`)
  - Type annotations: Removed `Optional`, `Dict`, `Tuple`, `List`

* **Fixes Applied**:
  - Converted `mask = property(__get_mask, __set_mask)` to `@property` decorators
  - Replaced Python-specific type hints with Codon equivalents
  - Changed `bytes` type hints to `str`
  - All 4 files now compile successfully

## Test Framework Creation
* **File Context**: Biopython test files `test_motifs.py` and `test_motifs_online.py`
* **Prompt**:
```
Using test_motifs.py and test_motifs_online.py I want you to put all the tests in test.py but 
you must make some changes because I am required to do this: Have a single file test.py that 
contains tests for both Codon and Python. Both codon test.py and python test.py should produce 
same results. The output should state whether or not the test passes. Do it like [PASS] or [FAIL]
```

* **Assistance Provided**:
  - Converted unittest methods to `@python` decorated functions
  - Restructured tests to return result tuples instead of using assertions
  - Fixed data file paths (./data/motifs/ instead of current directory)
  - Converted multi-line strings to single-line with \n escapes
  - Implemented environment detection using `hasattr(str, 'memcpy')`

## Environment Detection Solution
* **Problem**: Initial `try-except` with `__codon__` failed because Codon parses entire file at compile time
* **Solution**: 
```python
if hasattr(str, 'memcpy'):
    # Codon
    import code.bio_codon.motifs as motifs
    from python import Bio as _Bio
else:
    # Python  
    from Bio import motifs
    from Bio.Seq import Seq
```

## Codon Setup Requirements
* **Package Structure**:
  - Created `__init__.codon` in `code/` and `code/bio_codon/`
  - Fixed numpy import: `from python import numpy as np`

* **Environment Variables**:
```bash
export CODON_PYTHON=/usr/lib/x86_64-linux-gnu/libpython3.13.so
export PYTHON_PATH=/home/ryan/.local/lib/python3.13/site-packages
```

## Key Learnings with AI Assistance
1. **Codon Differences**: Stricter type checking, no `property()` function, different exception handling
2. **Runtime Detection**: `hasattr(str, 'memcpy')` works better than compile-time checks
3. **Python Interop**: `from python import` allows using Python libraries at runtime cost
4. **Dead Code**: Codon catches unused imports that Python ignores
5. **Single Test File**: Possible to create unified tests with conditional imports