"""
Week 4: Sequence Alignment Algorithms
Dynamic Programming Implementation

This module implements various sequence alignment algorithms:
- Global alignment (Needleman-Wunsch)
- Local alignment (Smith-Waterman)
- Semi-global/Fitting alignment
- Affine gap penalty global alignment
"""

from .global_alignment import global_align
from .local_alignment import local_align
from .semiglobal_alignment import semiglobal_align
from .affine_alignment import affine_align
from .utils import read_fasta

__all__ = [
    'global_align',
    'local_align',
    'semiglobal_align',
    'affine_align',
    'read_fasta',
]
