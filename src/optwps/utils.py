"""
Utility functions for optwps package.

This module provides helper functions for BAM file processing and file I/O operations.
"""

_open = open

import os
import sys
import pigz
import gzip

from contextlib import nullcontext


def is_soft_clipped(cigar):
    """
    Check if a read has soft clipping in its CIGAR string.

    Soft clipping (op=4) indicates that some bases at the start or end of the read
    are not aligned to the reference but are present in the sequence.

    Args:
        cigar (list): CIGAR tuples from pysam AlignedSegment.cigartuples
            Each tuple is (operation, length)

    Returns:
        bool: True if any soft clipping operation is present, False otherwise
    """
    return any(op == 4 for op, _ in cigar)


def ref_aln_length(cigar):
    """
    Calculate the length of alignment on the reference sequence from CIGAR.

    Computes the total length consumed on the reference by summing lengths of
    operations that consume reference bases: M(0), D(2), N(3), =(7), X(8).

    Args:
        cigar (list): CIGAR tuples from pysam AlignedSegment.cigartuples
            Each tuple is (operation, length)

    Returns:
        int: Total length on reference sequence
    """
    return sum(l for op, l in cigar if op in (0, 2, 3, 7, 8))


def exopen(fil: str, mode: str = "r", *args, use_pigz=True, njobs=-1, **kwargs):
    """
    Open a file with automatic gzip support and parallel compression.

    This function wraps the standard open() function with automatic detection
    and handling of gzipped files. When writing gzipped files, parallel compression
    (pigz) can be used for better performance on multi-core systems. Also supports
    writing to stdout when fil='stdout'.

    Args:
        fil (str): Path to the file to open, or 'stdout' for standard output
        mode (str, optional): File open mode ('r', 'w', 'rb', 'wb', etc.).
            Default: 'r'
        *args: Additional positional arguments passed to open function
        use_pigz (bool, optional): Whether to use pigz for parallel gzip compression.
            Falls back to standard gzip when pigz is unavailable or when handling
            multiple concurrent writers. Default: True
        njobs (int, optional): Number of parallel jobs for gzip compression.
            If -1, uses all available CPU cores. Default: -1
        **kwargs: Additional keyword arguments passed to open function

    Returns:
        file object: Opened file handle (stdout, standard file, or gzipped file)
    """
    if njobs == -1:
        njobs = os.cpu_count()
    if fil == "stdout":
        assert "r" not in mode, "Cannot open stdout in read mode"
        return nullcontext(sys.stdout)
    if fil.endswith(".gz"):
        open_func = pigz.open if use_pigz else gzip.open
        try:
            return open_func(
                fil, mode + "t" if not mode.endswith("b") else mode, *args, **kwargs
            )
        except BaseException:
            return open_func(fil, mode + "t" if not mode.endswith("b") else mode)

    return _open(fil, mode, *args, **kwargs)
