#!/usr/bin/env python
"""
Fast Window Protection Score (WPS) Calculator for Cell-Free DNA Analysis.

This module provides efficient computation of Window Protection Scores from BAM files,
designed for analyzing cell-free DNA fragmentation patterns and nucleosome positioning.

The WPS algorithm identifies protected genomic regions by analyzing DNA fragment
coverage patterns. It counts fragments that span versus fragments whose endpoints
fall within a protection window around each genomic position.

Example:
    Basic usage::

        from lib.optwps import WPS

        wps_calc = WPS(protection_size=120)
        wps_calc.run(bamfile='input.bam', out_filepath='output.tsv')

Note:
    - Op BAM Description
    - M 0 alignment match (can be a sequence match or mismatch)
    - I 1 insertion to the reference
    - D 2 deletion from the reference
    - N 3 skipped region from the reference
    - S 4 soft clipping (clipped sequences present in SEQ)
    - H 5 hard clipping (clipped sequences NOT present in SEQ)
    - P 6 padding (silent deletion from padded reference)
    - = 7 sequence match
    - X 8 sequence mismatch
    - pysam always uses 0-based coordinates, converting from 1-based of the sam files!
"""

import os
import pysam
import random
import numpy as np

from lib.utils import exopen, is_soft_clipped, ref_aln_length
from tqdm.auto import tqdm


class ROIGenerator:
    """
    Generate regions of interest (ROI) for processing.

    This class generates genomic regions either from a BED file or from the entire
    genome in the BAM file. Regions are yielded in chunks for memory-efficient processing.

    Args:
        bed_file (str, optional): Path to BED file containing regions to process.
            If None, the entire genome will be processed. Default: None
        chunk_size (int, optional): Size of chunks in base pairs for processing.
            Default: 1e8 (100 megabases)
    """

    def __init__(self, bed_file=None, chunk_size=1e8):
        self.bed_file = bed_file
        self.chunk_size = chunk_size

    def regions(self, bam_file=None):
        """
        Generate regions for processing.

        Yields genomic regions either from a BED file or from the entire genome
        referenced in the BAM file. Regions are chunked according to chunk_size.

        Args:
            bam_file (str, optional): Path to BAM file. Required if no BED file
                is provided to determine genome regions. Default: None

        Yields:
            tuple: (chromosome, start, end) for each region chunk

        Raises:
            ValueError: If neither bed_file nor bam_file can provide regions
        """
        if (self.bed_file is None) or (not os.path.exists(self.bed_file)):
            input_file = pysam.Samfile(bam_file, "rb")
            nchunks = sum(
                (input_file.get_reference_length(chrom) - 1) // self.chunk_size + 1
                for chrom in input_file.references
            )
            iterator = tqdm(total=nchunks, desc="Processing genome regions")
            for chrom in input_file.references:
                chrom_length = input_file.get_reference_length(chrom)
                region_start = 0
                while region_start < chrom_length:
                    region_end = min(region_start + self.chunk_size, chrom_length)
                    yield chrom, region_start, region_end
                    region_start = region_end
                    iterator.update(1)
            iterator.close()
        else:
            # read number of lines in bed file
            nlines = sum(1 for _ in exopen(self.bed_file, "r"))
            with exopen(self.bed_file, "r") as bed:
                for line in tqdm(bed, total=nlines, desc="Processing BED regions"):
                    chrom, start, end = line.split()[:3]
                    chunk_start = int(start)
                    chunk_end = min(chunk_start + self.chunk_size, int(end))
                    while chunk_start < int(end):
                        yield chrom, chunk_start, chunk_end
                        chunk_start = chunk_end
                        chunk_end = min(chunk_start + self.chunk_size, int(end))


class WPS:
    """
    Window Protection Score (WPS) calculator for cell-free DNA analysis.

    This class computes Window Protection Scores from aligned sequencing reads in BAM format.
    WPS quantifies the protection of DNA fragments around each genomic position, useful for
    identifying nucleosome positioning and other protected regions in cell-free DNA.

    The algorithm calculates:
        - Outside score: fragments that completely span the protection window
        - Inside score: fragment endpoints falling within the protection window
        - WPS = outside - inside

    Args:
        bed_file (str, optional): Path to BED file with regions to process.
            If None, processes entire genome. Default: None
        protection_size (int, optional): Total protection window size in base pairs.
            This value is divided by 2 to get the window on each side of the position.
            Default: 120
        min_insert_size (int, optional): Minimum insert/fragment size to include.
            If None, no minimum filter applied. Default: None
        max_insert_size (int, optional): Maximum insert/fragment size to include.
            If None, no maximum filter applied. Default: None
        valid_chroms (set, optional): Set of valid chromosome names to process.
            Default: chromosomes 1-22, X, Y
        chunk_size (float, optional): Region chunk size for processing.
            Default: 1e8 (100 Mb)

    Attributes:
        bed_file (str): Path to BED file or None
        protection_size (int): Half of the protection window size
        valid_chroms (set): Set of valid chromosome names
        min_insert_size (int): Minimum fragment size filter
        max_insert_size (int): Maximum fragment size filter
        chunk_size (float): Chunk size for processing
        roi_generator (ROIGenerator): Region generator instance

    Example:
        >>> wps = WPS(protection_size=120, min_insert_size=120, max_insert_size=180)
        >>> wps.run(bamfile='sample.bam', out_filepath='wps_output.tsv')

    Note:
        - Automatically filters duplicate, QC-failed, and unmapped reads
        - Handles both paired-end and single-end sequencing data
        - Supports downsampling for high-coverage samples
    """

    def __init__(
        self,
        bed_file=None,
        protection_size=120,
        min_insert_size=None,
        max_insert_size=None,
        valid_chroms=set(map(str, list(range(1, 23)) + ["X", "Y"])),
        chunk_size=1e8,
    ):
        self.bed_file = bed_file
        self.protection_size = protection_size // 2
        if valid_chroms is not None:
            self.valid_chroms = [x.replace("chr", "") for x in valid_chroms]
        else:
            self.valid_chroms = None
        self.min_insert_size = min_insert_size
        self.max_insert_size = max_insert_size
        self.chunk_size = chunk_size
        self.roi_generator = ROIGenerator(bed_file=self.bed_file)

    def __call__(self, *args, **kwargs):
        return self.run(*args, **kwargs)

    def run(
        self,
        bamfile,
        out_filepath=None,
        downsample_ratio=None,
        verbose_output=False,
    ):
        """
        Calculate Window Protection Score for all regions and write to file.

        Processes the BAM file to compute WPS values for each genomic position
        in the specified regions (or entire genome). Results are written to a
        tab-separated output file.

        Args:
            bamfile (str): Path to input BAM file (must be sorted and indexed)
            out_filepath (str): Path to output TSV file, or stdout (default)
            downsample_ratio (float, optional): Fraction of reads to randomly keep
                (0.0 to 1.0). Useful for high-coverage samples. Default: None (no downsampling)

        Returns:
            None: Results are written directly to out_filepath

        Output Format:
            Tab-separated file with columns:
                1. chromosome: Chromosome name (without 'chr' prefix)
                2. start: Start position (0-based)
                3. end: End position (1-based, start + 1)
                4. outside: Count of fragments spanning the protection window (if verbose_output)
                5. inside: Count of fragment endpoints in protection window (if verbose_output)
                6. wps: Window Protection Score (outside - inside)

        Raises:
            FileNotFoundError: If bamfile does not exist
            ValueError: If downsample_ratio is not between 0 and 1

        Example:
            >>> wps = WPS()
            >>> wps.run(bamfile='input.bam', out_filepath='output.tsv')
            >>> # With downsampling
            >>> wps.run(bamfile='input.bam', out_filepath='output.tsv',
            ...         downsample_ratio=0.5)
        """
        if out_filepath is None:
            out_filepath = "stdout"
        input_file = pysam.Samfile(bamfile, "rb")
        prefix = (
            "chr" if any(r.startswith("chr") for r in input_file.references) else ""
        )
        with exopen(out_filepath, "w") as outfile:

            for chrom, start, end in self.roi_generator.regions(bam_file=bamfile):
                if "chr" in chrom:
                    chrom = chrom.replace("chr", "")
                if self.valid_chroms is not None and chrom not in self.valid_chroms:
                    continue
                try:
                    regionStart, regionEnd = int(start), int(end)
                except ValueError:
                    continue

                starts = []
                ends = []
                for read in input_file.fetch(
                    prefix + chrom,
                    max(0, regionStart - self.protection_size - 1),
                    regionEnd + self.protection_size + 1,
                ):
                    if read.is_duplicate or read.is_qcfail or read.is_unmapped:
                        continue
                    if is_soft_clipped(read.cigartuples):
                        continue

                    if read.is_paired:
                        if read.mate_is_unmapped:
                            continue
                        if read.rnext != read.tid:
                            continue
                        if read.is_read1 or (
                            read.is_read2
                            and read.pnext + read.qlen
                            < regionStart - self.protection_size - 1
                        ):
                            if read.isize == 0:
                                continue
                            if (
                                downsample_ratio is not None
                                and random.random() >= downsample_ratio
                            ):
                                continue
                            lseq = abs(read.isize)
                            if (
                                self.min_insert_size is not None
                                and lseq < self.min_insert_size
                            ):
                                continue
                            if (
                                self.max_insert_size is not None
                                and lseq > self.max_insert_size
                            ):
                                continue
                            rstart = min(read.pos, read.pnext)
                            rend = rstart + lseq - 1
                            starts.append(rstart)
                            ends.append(rend)
                    else:
                        if (
                            downsample_ratio is not None
                            and random.random() >= downsample_ratio
                        ):
                            continue
                        rstart = read.pos
                        lseq = ref_aln_length(read.cigartuples)
                        if self.min_insert_size is not None and (
                            (lseq < self.min_insert_size)
                            or (lseq > self.max_insert_size)
                        ):
                            continue
                        rend = rstart + lseq - 1  # end included
                        starts.append(rstart)
                        ends.append(rend)
                n = regionEnd - regionStart + 1
                if len(starts) > 0:
                    starts = np.array(starts)
                    ends = np.array(ends)

                    # Fragments fully spanning the window boundaries
                    span_start = starts + self.protection_size - regionStart
                    span_end = ends - self.protection_size - regionStart + 2
                    valid = span_end >= span_start
                    outside = np.zeros(n + 2, dtype=int)
                    np.add.at(
                        outside,
                        np.clip(span_start[valid] + 1, 0, n + 1),
                        1,
                    )
                    np.add.at(
                        outside,
                        np.clip(span_end[valid], 0, n + 1),
                        -1,
                    )
                    np.add.at(
                        outside,
                        np.clip(span_start[~valid] + 1, 0, n + 1),
                        -1,
                    )
                    np.add.at(
                        outside,
                        np.clip(span_end[~valid], 0, n + 1),
                        1,
                    )

                    outside_cum = np.cumsum(outside)[:-1]

                    # Fragments whose endpoints fall inside windows
                    all_ends = np.concatenate([starts, ends]) - regionStart
                    left = np.clip(all_ends - self.protection_size + 2, 0, n + 1)
                    right = np.clip(all_ends + self.protection_size + 1, 0, n + 1)
                    inside = np.zeros(n + 2, dtype=int)
                    np.add.at(inside, left, 1)
                    np.add.at(inside, right, -1)
                    inside_cum = np.cumsum(inside)[:-1]

                    wps = outside_cum - inside_cum
                else:
                    outside_cum = np.zeros(n, dtype=int)
                    inside_cum = np.zeros(n, dtype=int)
                    wps = np.zeros(n, dtype=int)
                for i in range(regionEnd - regionStart + 1):
                    if verbose_output:
                        outfile.write(
                            "%s\t%d\t%d\t%d\t%d\t%d\n"
                            % (
                                chrom,
                                regionStart + i,
                                regionStart + i + 1,  # add end coordinate
                                outside_cum[i],
                                inside_cum[i],
                                wps[i],
                            )
                        )
                    else:
                        outfile.write(
                            "%s\t%d\t%d\t%d\n"
                            % (
                                chrom,
                                regionStart + i,
                                regionStart + i + 1,  # add end coordinate
                                wps[i],
                            )
                        )
