"""
Local Alignment (Smith-Waterman Algorithm)

Finds the best local alignment between two sequences using dynamic programming.
"""

import numpy as np

def local_align(seq1, seq2, match=3, mismatch=-3, gap=-2):
    """
    Perform local alignment of two sequences using Smith-Waterman algorithm.
    
    Args:
        seq1: First sequence (string)
        seq2: Second sequence (string)
        match: Score for matching characters (default: 3)
        mismatch: Score for mismatching characters (default: -3)
        gap: Score for gap (default: -2)
        
    Returns:
        tuple: (alignment_score, aligned_seq1, aligned_seq2)
    """
    m, n = len(seq1), len(seq2)
    
    # Initialize DP matrix using numpy for memory efficiency
    dp = np.zeros((m + 1, n + 1), dtype=np.int32)
    
    # Track maximum score and its position
    max_score = 0
    max_i, max_j = 0, 0
    
    # Fill DP matrix (no initialization of first row/column needed - they stay 0)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Score for match or mismatch
            if seq1[i-1] == seq2[j-1]:
                diag_score = dp[i-1, j-1] + match
            else:
                diag_score = dp[i-1, j-1] + mismatch
            
            # Score for gap in seq2
            up_score = dp[i-1, j] + gap
            
            # Score for gap in seq1
            left_score = dp[i, j-1] + gap
            
            # Take maximum, but at least 0 (this is the key difference from global)
            dp[i, j] = max(0, diag_score, up_score, left_score)
            
            # Track maximum score position
            if dp[i, j] > max_score:
                max_score = int(dp[i, j])
                max_i, max_j = i, j
    
    # Traceback from maximum score position until we hit 0
    aligned1, aligned2 = [], []
    i, j = max_i, max_j
    
    while i > 0 and j > 0 and dp[i, j] > 0:
        # Check if current cell came from diagonal
        if seq1[i-1] == seq2[j-1]:
            score_diag = dp[i-1, j-1] + match
        else:
            score_diag = dp[i-1, j-1] + mismatch
        
        if dp[i, j] == score_diag and score_diag > 0:
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and dp[i, j] == dp[i-1, j] + gap and dp[i-1, j] + gap > 0:
            # Came from above (gap in seq2)
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        elif j > 0 and dp[i, j] == dp[i, j-1] + gap and dp[i, j-1] + gap > 0:
            # Came from left (gap in seq1)
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1
        else:
            # Stop if we can't determine the path or hit 0
            break
    
    # Reverse alignments (we built them backwards)
    aligned1 = ''.join(reversed(aligned1))
    aligned2 = ''.join(reversed(aligned2))
    
    return max_score, aligned1, aligned2
