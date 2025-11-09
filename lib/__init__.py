"""
fast_wps - Fast Window Protection Score Calculator

A high-performance Python package for computing Window Protection Score (WPS) from BAM files,
designed for cell-free DNA (cfDNA) analysis.

Main Components:
    - WPS: Main class for Window Protection Score calculation
    - ROIGenerator: Generates regions of interest from BED files or genomes
    - Utility functions: File I/O and BAM processing helpers

Example:
    >>> from lib.fast_wps import WPS
    >>> wps = WPS(protection_size=120)
    >>> wps.run(bamfile='input.bam', out_filepath='output.tsv')

Version: 1.0.0
"""

__version__ = "1.0.0"
__author__ = "JorisVermeeschLab"

from .fast_wps import WPS, ROIGenerator
from .utils import exopen, isSoftClipped, ref_aln_length

__all__ = ["WPS", "ROIGenerator", "exopen", "isSoftClipped", "ref_aln_length"]
