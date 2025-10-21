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
    INF = 10**9  # Use a large number instead of float('inf') for numpy
    
    # Initialize three DP matrices using numpy for memory efficiency
    M = np.full((m + 1, n + 1), -INF, dtype=np.int32)  # Match/mismatch
    I = np.full((m + 1, n + 1), -INF, dtype=np.int32)  # Gap in seq1 (insertion)
    D = np.full((m + 1, n + 1), -INF, dtype=np.int32)  # Gap in seq2 (deletion)
    
    # Base case: M[0][0] = 0
    M[0, 0] = 0
    
    # Initialize first row (gaps in seq2)
    for i in range(1, m + 1):
        D[i, 0] = gap_open + (i - 1) * gap_extend
    
    # Initialize first column (gaps in seq1)
    for j in range(1, n + 1):
        I[0, j] = gap_open + (j - 1) * gap_extend
    
    # Fill DP matrices
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # M[i][j]: match or mismatch
            if seq1[i-1] == seq2[j-1]:
                match_score = match
            else:
                match_score = mismatch
            
            M[i, j] = max(
                M[i-1, j-1] + match_score,  # From match state
                I[i-1, j-1] + match_score,  # From insertion state
                D[i-1, j-1] + match_score   # From deletion state
            )
            
            # I[i][j]: gap in seq1
            I[i, j] = max(
                M[i, j-1] + gap_open + gap_extend,  # Open new gap
                I[i, j-1] + gap_extend,              # Extend gap
                D[i, j-1] + gap_open + gap_extend   # Switch from D to I
            )
            
            # D[i][j]: gap in seq2
            D[i, j] = max(
                M[i-1, j] + gap_open + gap_extend,  # Open new gap
                I[i-1, j] + gap_open + gap_extend,  # Switch from I to D
                D[i-1, j] + gap_extend              # Extend gap
            )
    
    # Find best final score
    final_score = int(max(M[m, n], I[m, n], D[m, n]))
    
    # Traceback to reconstruct alignment
    aligned1, aligned2 = [], []
    i, j = m, n
    
    # Determine which matrix we ended in
    if final_score == M[m, n]:
        state = 'M'
    elif final_score == I[m, n]:
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
            
            if seq1[i-1] == seq2[j-1]:
                match_score = match
            else:
                match_score = mismatch
            
            # Find which state we came from
            if M[i, j] == M[i-1, j-1] + match_score:
                state = 'M'
            elif M[i, j] == I[i-1, j-1] + match_score:
                state = 'I'
            else:  # M[i, j] == D[i-1, j-1] + match_score
                state = 'D'
            
            i -= 1
            j -= 1
            
        elif state == 'I':
            if j == 0:
                break
            
            # We're in insertion state (gap in seq1)
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            
            # Find which state we came from
            if I[i, j] == M[i, j-1] + gap_open + gap_extend:
                state = 'M'
            elif I[i, j] == I[i, j-1] + gap_extend:
                state = 'I'
            else:  # I[i, j] == D[i, j-1] + gap_open + gap_extend
                state = 'D'
            
            j -= 1
            
        else:  # state == 'D'
            if i == 0:
                break
            
            # We're in deletion state (gap in seq2)
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            
            # Find which state we came from
            if D[i, j] == M[i-1, j] + gap_open + gap_extend:
                state = 'M'
            elif D[i, j] == I[i-1, j] + gap_open + gap_extend:
                state = 'I'
            else:  # D[i, j] == D[i-1, j] + gap_extend
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
