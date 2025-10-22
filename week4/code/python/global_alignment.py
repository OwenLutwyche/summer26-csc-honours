"""
Global Alignment (Needleman-Wunsch Algorithm)

Aligns two sequences end-to-end using dynamic programming.
"""

def global_align(seq1, seq2, match=3, mismatch=-3, gap=-2):
    """
    Perform global alignment of two sequences using Needleman-Wunsch algorithm.
    
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
    
    # Initialize DP matrix using Python lists (much faster than NumPy for this use case)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Initialize first row and column (gap penalties)
    for i in range(1, m + 1):
        dp[i][0] = i * gap
    for j in range(1, n + 1):
        dp[0][j] = j * gap
    
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
    
    # Traceback to get alignment
    aligned1, aligned2 = [], []
    i, j = m, n
    
    while i > 0 or j > 0:
        if i > 0 and j > 0:
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
                continue
        
        if i > 0 and dp[i][j] == dp[i-1][j] + gap:
            # Came from above (gap in seq2)
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        elif j > 0:
            # Came from left (gap in seq1)
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1
        else:
            # Only up moves remain
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
    
    # Reverse alignments (we built them backwards)
    aligned1 = ''.join(reversed(aligned1))
    aligned2 = ''.join(reversed(aligned2))
    
    return dp[m][n], aligned1, aligned2
