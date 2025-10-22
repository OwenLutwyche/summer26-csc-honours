"""
Affine Gap Penalty Global Alignment

Global alignment with affine gap penalties (different costs for opening vs extending gaps).
Uses three matrices: M (match), I (insertion/gap in seq1), D (deletion/gap in seq2).
"""

import numpy as np

def affine_align(seq1, seq2, match=3, mismatch=-3, gap_open=-5, gap_extend=-1):
    """
    Perform global alignment with affine gap penalties.
    
    Uses three DP matrices:
    - M[i][j]: best score ending with alignment of seq1[i-1] and seq2[j-1]
    - I[i][j]: best score ending with gap in seq1 (insertion in seq2)
    - D[i][j]: best score ending with gap in seq2 (deletion from seq1)
    
    Args:
        seq1: First sequence (string)
        seq2: Second sequence (string)
        match: Score for matching characters (default: 3)
        mismatch: Score for mismatching characters (default: -3)
        gap_open: Penalty for opening a gap (default: -5)
        gap_extend: Penalty for extending a gap (default: -1)
        
    Returns:
        tuple: (alignment_score, aligned_seq1, aligned_seq2)
    """
    m, n = len(seq1), len(seq2)
    INF = 10**9  # Use a large number instead of float('inf')
    
    # For large sequences, use NumPy for memory efficiency
    # For small sequences, Python lists are fine
    use_numpy = m * n > 100000  # Use NumPy for sequences creating >100k cells
    
    if use_numpy:
        # Initialize three DP matrices using NumPy (more memory efficient for large matrices)
        M = np.full((m + 1, n + 1), -INF, dtype=np.int32)  # Match/mismatch
        I = np.full((m + 1, n + 1), -INF, dtype=np.int32)  # Gap in seq1 (insertion)
        D = np.full((m + 1, n + 1), -INF, dtype=np.int32)  # Gap in seq2 (deletion)
    else:
        # Use Python lists for small sequences (faster access)
        M = [[-INF] * (n + 1) for _ in range(m + 1)]  # Match/mismatch
        I = [[-INF] * (n + 1) for _ in range(m + 1)]  # Gap in seq1 (insertion)
        D = [[-INF] * (n + 1) for _ in range(m + 1)]  # Gap in seq2 (deletion)
    
    # Base case: M[0][0] = 0
    M[0][0] = 0
    
    # Initialize first row (gaps in seq2)
    for i in range(1, m + 1):
        D[i][0] = gap_open + (i - 1) * gap_extend
    
    # Initialize first column (gaps in seq1)
    for j in range(1, n + 1):
        I[0][j] = gap_open + (j - 1) * gap_extend
    
    # Fill DP matrices
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # M[i][j]: match or mismatch
            if seq1[i-1] == seq2[j-1]:
                match_score = match
            else:
                match_score = mismatch
            
            M[i][j] = max(
                M[i-1][j-1] + match_score,
                I[i-1][j-1] + match_score,
                D[i-1][j-1] + match_score
            )
            
            # I[i][j]: gap in seq1
            I[i][j] = max(
                M[i][j-1] + gap_open + gap_extend,
                I[i][j-1] + gap_extend,
                D[i][j-1] + gap_open + gap_extend
            )
            
            # D[i][j]: gap in seq2
            D[i][j] = max(
                M[i-1][j] + gap_open + gap_extend,
                I[i-1][j] + gap_open + gap_extend,
                D[i-1][j] + gap_extend
            )
    
    # Find best final score
    final_score = int(max(M[m, n] if use_numpy else M[m][n], 
                         I[m, n] if use_numpy else I[m][n], 
                         D[m, n] if use_numpy else D[m][n]))
    
    # Traceback to reconstruct alignment (recalculate to determine path)
    aligned1, aligned2 = [], []
    i, j = m, n
    
    # Helper function for accessing matrices (handles both NumPy and lists)
    def get_val(mat, i, j):
        return int(mat[i, j]) if use_numpy else mat[i][j]
    
    # Determine which matrix we ended in
    m_val = get_val(M, m, n)
    i_val = get_val(I, m, n)
    d_val = get_val(D, m, n)
    
    if final_score == m_val:
        state = 'M'
    elif final_score == i_val:
        state = 'I'
    else:
        state = 'D'
    
    while i > 0 or j > 0:
        if state == 'M':
            if i == 0 or j == 0:
                break
            
            # We're in match state
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            
            # Determine which state we came from by recalculating
            if seq1[i-1] == seq2[j-1]:
                match_score = match
            else:
                match_score = mismatch
            
            curr = get_val(M, i, j)
            if curr == get_val(M, i-1, j-1) + match_score:
                state = 'M'
            elif curr == get_val(I, i-1, j-1) + match_score:
                state = 'I'
            else:
                state = 'D'
            
            i -= 1
            j -= 1
            
        elif state == 'I':
            if j == 0:
                break
            
            # We're in insertion state (gap in seq1)
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            
            # Determine which state we came from
            curr = get_val(I, i, j)
            if curr == get_val(M, i, j-1) + gap_open + gap_extend:
                state = 'M'
            elif curr == get_val(I, i, j-1) + gap_extend:
                state = 'I'
            else:
                state = 'D'
            
            j -= 1
            
        else:  # state == 'D'
            if i == 0:
                break
            
            # We're in deletion state (gap in seq2)
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            
            # Determine which state we came from
            curr = get_val(D, i, j)
            if curr == get_val(M, i-1, j) + gap_open + gap_extend:
                state = 'M'
            elif curr == get_val(I, i-1, j) + gap_open + gap_extend:
                state = 'I'
            else:
                state = 'D'
            
            i -= 1
    
    # Handle remaining characters at the beginning (if any)
    while i > 0:
        aligned1.append(seq1[i-1])
        aligned2.append('-')
        i -= 1
    
    while j > 0:
        aligned1.append('-')
        aligned2.append(seq2[j-1])
        j -= 1
    
    # Reverse alignments (we built them backwards)
    aligned1 = ''.join(reversed(aligned1))
    aligned2 = ''.join(reversed(aligned2))
    
    return final_score, aligned1, aligned2
