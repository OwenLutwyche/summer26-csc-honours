# LLM Usage
Model Used: `Claude Sonnet 4.5`

## Test File Filtering and Setup
* **File Context**: Original biotite test file (`test_phylo.py`) with comprehensive test suite
* **Initial Prompt**:
```
Your job is to port biotite's phylo package to Codon. Note that the source is written in Cython, not pure Python.

The source code can be found on Github. You only need to port the parts needed to run the two tests in test_phylo.py:
    test_distances,
    test_upgma,
    test_neighbor_joining.

Using the original `test_phylo.py`, make a new `test_phylo_filtered.py` that only tests the parts I need for when I have to have to test my Codon implementation. You can remove all other tests and any code that is not needed for these three tests to run.
```

* **Key Assistance Provided**:
  - Analyzed original test file to identify required tests and their dependencies
  - Created filtered test file (`test_phylo_filtered.py`) with only the three required tests
  - Implemented timing functionality using Python's `time` module to measure test execution in milliseconds
  - Removed pytest dependency and converted fixtures to regular functions for standalone execution

## Data Management
* **Assistance Provided**:
  - Used `os.path` operations to ensure cross-platform compatibility

## Test Structure Modifications
* **Original Structure**: pytest-based with fixtures (`@pytest.fixture`)
* **Modified Structure**: Standalone Python script with:
  - Helper functions `get_distances()` and `get_upgma_newick()` to load data
  - Individual test functions for each required test
  - Custom `run_tests()` function with timing for each test
  - User-friendly output with pass/fail indicators and millisecond timing
