"""
Semi-Global (Fitting) Alignment

Aligns two sequences allowing free gaps at the beginning and end.
Useful when one sequence might be contained within another.
"""

def semiglobal_align(seq1, seq2, match=3, mismatch=-3, gap=-2):
    """
    Perform semi-global (fitting) alignment of two sequences.
    
    This alignment allows free gaps at the beginning and end of both sequences.
    It's useful for finding where a shorter sequence fits within a longer one.
    
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
    
    # Initialize first row and column to 0 (free gaps at start)
    # Already done by list initialization
    
    # Fill DP matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Score for match or mismatch
            if seq1[i-1] == seq2[j-1]:
                diag_score = dp[i-1][j-1] + match
            else:
                diag_score = dp[i-1][j-1] + mismatch
            
            # Score for gap in seq2
            up_score = dp[i-1][j] + gap
            
            # Score for gap in seq1
            left_score = dp[i][j-1] + gap
            
            # Take maximum
            dp[i][j] = max(diag_score, up_score, left_score)
    
    # Find maximum score in last row or last column (free gaps at end)
    max_score = float('-inf')
    max_i, max_j = m, n
    
    # Check last row
    for j in range(n + 1):
        if dp[m][j] > max_score:
            max_score = dp[m][j]
            max_i, max_j = m, j
    
    # Check last column
    for i in range(m + 1):
        if dp[i][n] > max_score:
            max_score = dp[i][n]
            max_i, max_j = i, n
    
    # Traceback from the maximum score position
    aligned1, aligned2 = [], []
    i, j = max_i, max_j
    
    while i > 0 and j > 0:
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
            # Fallback
            break
    
    # Reverse alignments (we built them backwards)
    aligned1 = ''.join(reversed(aligned1))
    aligned2 = ''.join(reversed(aligned2))
    
    return max_score, aligned1, aligned2
