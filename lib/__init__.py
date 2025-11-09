"""
optwps - Fast Window Protection Score Calculator

A high-performance Python package for computing Window Protection Score (WPS) from BAM files,
designed for cell-free DNA (cfDNA) analysis.

Example:
    >>> from lib.optwps import WPS
    >>> wps = WPS(protection_size=120)
    >>> wps.run(bamfile='input.bam', out_filepath='output.tsv')

Version: 1.0.0
"""

__version__ = "1.0.0"
__author__ = "VasLem"

from .optwps import WPS

__all__ = ["WPS"]
