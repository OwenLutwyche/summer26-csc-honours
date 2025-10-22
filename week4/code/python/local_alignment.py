"""
Local Alignment (Smith-Waterman Algorithm)

Finds the best local alignment between two sequences using dynamic programming.
"""

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
    
    # Initialize DP matrix using Python lists
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Track maximum score and its position
    max_score = 0
    max_i, max_j = 0, 0
    
    # Fill DP matrix (no initialization of first row/column needed - they stay 0)
    for i in range(1, m + 1):
        curr_row = dp[i]
        prev_row = dp[i-1]
        for j in range(1, n + 1):
            # Cache character comparison
            if seq1[i-1] == seq2[j-1]:
                diag_score = prev_row[j-1] + match
            else:
                diag_score = prev_row[j-1] + mismatch
            
            # Score for gap in seq2
            up_score = prev_row[j] + gap
            
            # Score for gap in seq1
            left_score = curr_row[j-1] + gap
            
            # Take maximum, but at least 0 (inline comparison)
            temp = diag_score if diag_score > up_score else up_score
            temp = temp if temp > left_score else left_score
            curr_row[j] = temp if temp > 0 else 0
            
            # Track maximum score position
            if curr_row[j] > max_score:
                max_score = curr_row[j]
                max_i, max_j = i, j
    
    # Traceback from maximum score position until we hit 0
    aligned1, aligned2 = [], []
    i, j = max_i, max_j
    
    while i > 0 and j > 0 and dp[i][j] > 0:
        # Check if current cell came from diagonal
        if seq1[i-1] == seq2[j-1]:
            score_diag = dp[i-1][j-1] + match
        else:
            score_diag = dp[i-1][j-1] + mismatch
        
        if dp[i][j] == score_diag:
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i-1][j] + gap:
            # Came from above (gap in seq2)
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        elif j > 0 and dp[i][j] == dp[i][j-1] + gap:
            # Came from left (gap in seq1)
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1
        else:
            # Stop if we can't determine the path
            break
    
    # Reverse alignments (we built them backwards)
    aligned1 = ''.join(reversed(aligned1))
    aligned2 = ''.join(reversed(aligned2))
    
    return max_score, aligned1, aligned2
