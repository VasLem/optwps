"""
Utility functions for fast_wps package.

This module provides helper functions for BAM file processing and file I/O operations.
"""

_open = open

import os
import pgzip


def isSoftClipped(cigar):
    """
    Check if a read has soft clipping in its CIGAR string.

    Soft clipping (op=4) indicates that some bases at the start or end of the read
    are not aligned to the reference but are present in the sequence.

    Args:
        cigar (list): CIGAR tuples from pysam AlignedSegment.cigartuples
            Each tuple is (operation, length)

    Returns:
        bool: True if any soft clipping operation is present, False otherwise

    Example:
        >>> cigar = [(0, 50), (4, 10)]  # 50M10S
        >>> isSoftClipped(cigar)
        True
        >>> cigar = [(0, 60)]  # 60M
        >>> isSoftClipped(cigar)
        False
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

    Example:
        >>> cigar = [(0, 50), (2, 5), (0, 45)]  # 50M5D45M
        >>> ref_aln_length(cigar)
        100

    Note:
        Operations included:
            - 0 (M): alignment match/mismatch
            - 2 (D): deletion from reference
            - 3 (N): skipped region from reference
            - 7 (=): sequence match
            - 8 (X): sequence mismatch
    """
    return sum(l for op, l in cigar if op in (0, 2, 3, 7, 8))


def exopen(fil: str, mode: str = "r", *args, njobs=-1, **kwargs):
    """
    Open a file with automatic gzip support and parallel compression.

    This function wraps the standard open() function with automatic detection
    and handling of gzipped files using parallel compression (pgzip) for better
    performance on multi-core systems.

    Args:
        fil (str): Path to the file to open
        mode (str, optional): File open mode ('r', 'w', 'rb', 'wb', etc.).
            Default: 'r'
        *args: Additional positional arguments passed to open function
        njobs (int, optional): Number of parallel jobs for gzip compression.
            If -1, uses all available CPU cores. Default: -1
        **kwargs: Additional keyword arguments passed to open function

    Returns:
        file object: Opened file handle (either standard or pgzip)

    Example:
        >>> # Reading a gzipped file
        >>> with exopen('data.txt.gz', 'r') as f:
        ...     content = f.read()

        >>> # Writing to a gzipped file with parallel compression
        >>> with exopen('output.txt.gz', 'w') as f:
        ...     f.write('Hello, world!')

        >>> # Regular file (no gzip)
        >>> with exopen('data.txt', 'r') as f:
        ...     content = f.read()

    Note:
        - Automatically detects gzipped files by .gz extension
        - Uses pgzip for parallel compression/decompression
        - Falls back to standard open() for non-gzipped files
        - Handles text and binary modes appropriately
    """

    if njobs == -1:
        njobs = os.cpu_count()
    if fil.endswith(".gz"):
        try:
            return pgzip.open(
                fil, mode + "t" if not mode.endswith("b") else mode, *args, **kwargs
            )
        except BaseException:
            return pgzip.open(fil, mode + "t" if not mode.endswith("b") else mode)
    return _open(fil, mode, *args, **kwargs)
