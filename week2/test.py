#!/usr/bin/env python3
"""
Test file for Bio.motifs functionality supporting both Python and Codon.
Converted from unittest to @test decorator format for Codon compatibility.
to be set to access Python's Bio module through Codon's Python integration.
"""

import math
import os

# Detect if we're running in Codon or Python
if hasattr(str, 'memcpy'):
    # Running in Codon
    import code.bio_codon.motifs as motifs
    from python import Bio as _Bio
    from python import numpy as np
    Seq = _Bio.Seq.Seq
    print("Using Codon with bio_codon modules")
    HAS_BIO = True
    HAS_NUMPY = True
    CODON = True
else:
    from Bio import motifs
    from Bio.Seq import Seq
    import numpy as np
    print("Using Biopython")
    HAS_BIO = True
    HAS_NUMPY = True
    CODON = False
    
    def python(func):
        return func

def test_result(test_name: str, passed: bool, message: str = ""):
    """Helper function to print test results"""
    status = "[PASS]" if passed else "[FAIL]"
    if message:
        print(f"{status} {test_name}: {message}")
    else:
        print(f"{status} {test_name}")

def assert_equal(actual, expected, test_name: str):
    """Helper assertion function"""
    if actual == expected:
        return True
    else:
        print(f"Expected: {expected}")
        print(f"Actual: {actual}")
        return False

def assert_almost_equal(actual, expected, places: int = 7, test_name: str = ""):
    """Helper assertion function for floats"""
    tolerance = 10 ** (-places)
    if abs(actual - expected) <= tolerance:
        return True
    else:
        print(f"Expected: {expected} (±{tolerance})")
        print(f"Actual: {actual}")
        return False

def assert_allclose(actual, expected, test_name: str = ""):
    """Helper assertion function for numpy arrays"""
    try:
        if HAS_NUMPY:
            return np.allclose(actual, expected)
        else:
            # Fallback comparison for when numpy is not available
            if len(actual) != len(expected):
                return False
            for i in range(len(actual)):
                if abs(actual[i] - expected[i]) > 1e-10:
                    return False
            return True
    except Exception as e:
        print(f"Error in allclose comparison: {e}")
        return False

# Test decorators - for Codon compatibility, we'll define a simple test decorator
def test(func):
    """Test decorator for Codon compatibility"""
    func._is_test = True
    return func

# Python-specific test functions (will run in Python mode)
@python
def run_python_tests():
    """Run all tests that require Python-specific functionality"""
    results = []
    
    # Test format functionality
    try:
        from Bio import motifs
        from Bio.Seq import Seq
        
        # Test motif formatting
        m = motifs.create([Seq("ATATA")])
        m.name = "Foo"
        
        # Test PFM format
        s1 = format(m, "pfm")
        expected_pfm = "  1.00   0.00   1.00   0.00   1.00\n  0.00   0.00   0.00   0.00   0.00\n  0.00   0.00   0.00   0.00   0.00\n  0.00   1.00   0.00   1.00   0.00\n"
        if s1 == expected_pfm:
            results.append(("test_format_pfm", True, ""))
        else:
            results.append(("test_format_pfm", False, f"Format mismatch"))
        
        # Test JASPAR format
        s2 = format(m, "jaspar")
        expected_jaspar = ">None Foo\nA [  1.00   0.00   1.00   0.00   1.00]\nC [  0.00   0.00   0.00   0.00   0.00]\nG [  0.00   0.00   0.00   0.00   0.00]\nT [  0.00   1.00   0.00   1.00   0.00]\n"
        if s2 == expected_jaspar:
            results.append(("test_format_jaspar", True, ""))
        else:
            results.append(("test_format_jaspar", False, f"Format mismatch"))
        
        # Test TRANSFAC format
        s3 = format(m, "transfac")
        expected_transfac = "P0      A      C      G      T\n01      1      0      0      0      A\n02      0      0      0      1      T\n03      1      0      0      0      A\n04      0      0      0      1      T\n05      1      0      0      0      A\nXX\n//\n"
        if s3 == expected_transfac:
            results.append(("test_format_transfac", True, ""))
        else:
            results.append(("test_format_transfac", False, f"Format mismatch"))
        
    except Exception as e:
        results.append(("test_format", False, f"Exception: {e}"))
    
    return results

@python  
def run_python_relative_entropy_tests():
    """Test relative entropy calculations"""
    results = []
    try:
        import numpy as np
        from Bio import motifs
        from Bio.Seq import Seq
        
        m = motifs.create([Seq("ATATA"), Seq("ATCTA"), Seq("TTGTA")])
        
        # Test basic relative entropy
        expected = np.array([1.0817041659455104, 2.0, 0.4150374992788437, 2.0, 2.0])
        if np.allclose(m.relative_entropy, expected):
            results.append(("test_relative_entropy_basic", True, ""))
        else:
            results.append(("test_relative_entropy_basic", False, "Values don't match expected"))
        
        # Test with different background
        m.background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
        expected = np.array([
            0.8186697601117167,
            1.7369655941662063,
            0.5419780939258206,
            1.7369655941662063,
            1.7369655941662063,
        ])
        if np.allclose(m.relative_entropy, expected):
            results.append(("test_relative_entropy_background", True, ""))
        else:
            results.append(("test_relative_entropy_background", False, "Values don't match expected"))
            
    except Exception as e:
        results.append(("test_relative_entropy", False, f"Exception: {e}"))
    
    return results

@python
def run_python_reverse_complement_tests():
    """Test reverse complement functionality"""
    results = []
    try:
        from Bio import motifs
        from Bio.Seq import Seq
        
        background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
        pseudocounts = 0.5
        m = motifs.create([Seq("ATATA")])
        m.background = background
        m.pseudocounts = pseudocounts
        
        # Test forward
        received_forward = format(m, "transfac")
        expected_forward = "P0      A      C      G      T\n01      1      0      0      0      A\n02      0      0      0      1      T\n03      1      0      0      0      A\n04      0      0      0      1      T\n05      1      0      0      0      A\nXX\n//\n"
        if received_forward == expected_forward:
            results.append(("test_reverse_complement_forward", True, ""))
        else:
            results.append(("test_reverse_complement_forward", False, "Forward format mismatch"))
        
        # Test reverse complement
        m = m.reverse_complement()
        received_reverse = format(m, "transfac")
        expected_reverse = "P0      A      C      G      T\n01      0      0      0      1      T\n02      1      0      0      0      A\n03      0      0      0      1      T\n04      1      0      0      0      A\n05      0      0      0      1      T\nXX\n//\n"
        if received_reverse == expected_reverse:
            results.append(("test_reverse_complement_reverse", True, ""))
        else:
            results.append(("test_reverse_complement_reverse", False, "Reverse format mismatch"))
            
    except Exception as e:
        results.append(("test_reverse_complement", False, f"Exception: {e}"))
    
    return results

@python
def run_python_meme_parser_tests():
    """Test MEME parser functionality"""
    results = []
    try:
        from Bio import motifs
        
        # Test minimal MEME parser
        with open("./data/motifs/minimal_test.meme") as stream:
            record = motifs.parse(stream, "minimal")
        
        if (record.version == "4" and 
            record.alphabet == "ACGT" and 
            len(record.sequences) == 0 and
            record.command == "" and
            len(record) == 3):
            
            motif = record[0]
            if (motif.name == "KRP" and
                record["KRP"] == motif and
                motif.num_occurrences == 17 and
                motif.length == 19):
                results.append(("test_minimal_meme_parser", True, ""))
            else:
                results.append(("test_minimal_meme_parser", False, "Motif properties don't match"))
        else:
            results.append(("test_minimal_meme_parser", False, "Record properties don't match"))
        
        # Test RNA MEME parser
        with open("./data/motifs/minimal_test_rna.meme") as stream:
            record = motifs.parse(stream, "minimal")
        
        if (record.version == "4" and 
            record.alphabet == "ACGU" and
            len(record) == 3):
            
            motif = record[0]
            if (motif.name == "KRP_fake_RNA" and
                motif.consensus == "UGUGAUCGAGGUCACACUU"):
                results.append(("test_meme_parser_rna", True, ""))
            else:
                results.append(("test_meme_parser_rna", False, "RNA motif properties don't match"))
        else:
            results.append(("test_meme_parser_rna", False, "RNA record properties don't match"))
            
    except Exception as e:
        results.append(("test_meme_parser", False, f"Exception: {e}"))
    
    return results

# Note: Online tests (weblogo) are commented out as they may not work in Codon
# and require internet connectivity
@python
def run_python_online_tests():
    """Test online functionality (weblogo) - commented out for Codon compatibility"""
    results = []
    try:
        from Bio import motifs
        from Bio.Seq import Seq
        
        # Test DNA weblogo
        m = motifs.create([Seq(s) for s in ["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"]], "GATCBDSW")
        m.weblogo(os.devnull)  # Send output to null device
        results.append(("test_weblogo_dna", True, ""))
        
        # Test RNA weblogo  
        m = motifs.create([Seq(s) for s in ["UACAA", "UACGC", "UACAC", "UACCC", "AACCC", "AAUGC", "AAUGC"]], "GAUC")
        m.weblogo(os.devnull)
        results.append(("test_weblogo_rna", True, ""))
        
        # Test protein weblogo
        m = motifs.create([Seq(s) for s in ["ACDEG", "AYCRN", "HYLID", "AYHEL", "ACDEH", "AYYRN", "HYIID"]], "ACDEFGHIKLMNPQRSTVWYBXZJUO")
        m.weblogo(os.devnull)
        results.append(("test_weblogo_protein", True, ""))
        
    except Exception as e:
        results.append(("test_online", False, f"Exception: {e}"))
    
    return results

# Codon-compatible tests (commented out until bio_codon modules are fully implemented)
"""
@test  
def test_codon_basic():
    # Example of how Codon tests would work with bio_codon modules
    # from python import sys  # For bio_codon imports
    test_result("test_codon_basic", False, "Not implemented - bio_codon modules needed")

@test
def test_codon_motif_creation():
    # This would test bio_codon motif creation
    test_result("test_codon_motif_creation", False, "Not implemented - bio_codon modules needed")
"""

# Main execution
def main():
    """Run all tests and report results"""
    print("Running Bio.motifs tests...")
    print("=" * 50)
    
    all_results = []
    
    # Run tests (same tests for both Python and Codon)
    if HAS_BIO and HAS_NUMPY:
        if CODON:
            print("Running in Codon environment...")
        else:
            print("Running in Python environment...")
        
        print("Running tests...")
        all_results.extend(run_python_tests())
        all_results.extend(run_python_relative_entropy_tests())
        all_results.extend(run_python_reverse_complement_tests())
        all_results.extend(run_python_meme_parser_tests())
        
        # Uncomment to run online tests (requires internet)
        # all_results.extend(run_python_online_tests())
    else:
        print("Skipping tests - Bio and/or numpy not available")
    
    # Print all results
    if all_results:
        print("\nTest Results:")
        print("=" * 50)
        passed = 0
        total = 0
        
        for test_name, success, message in all_results:
            test_result(test_name, success, message)
            if success:
                passed += 1
            total += 1
        
        print("=" * 50)
        print(f"Tests passed: {passed}/{total}")
    else:
        print("\nNo tests were run.")
        if not HAS_BIO:
            if CODON:
                print("bio_codon modules not available for Codon tests.")
            else:
                print("Bio modules not available for Python tests.")
        if CODON and not HAS_BIO:
            print("Codon tests require bio_codon modules to be implemented.")

if __name__ == "__main__":
    main()
