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

        from optwps import WPS

        wps_calc = WPS(protection_size=120)
        wps_calc.run(bamfile='input.bam', out_filepath='output.tsv')

    Creating separate output files::

        # Per chromosome
        wps_calc.run(bamfile='input.bam', out_filepath='wps_{chrom}.tsv')

        # Per target region
        wps_calc.run(bamfile='input.bam', out_filepath='wps_{target}.tsv.gz')

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
import sys
import pysam
import random
import numpy as np
import pandas as pd
import joblib

from optwps.utils import exopen, is_soft_clipped, ref_aln_length
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
            tuple: (chromosome, chunk_start, chunk_end, region_id)
                where chunk_start/chunk_end define the current chunk being processed
                and region_id is either from the BED file or constructed from chrom/start/end.

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
                    yield chrom, region_start, region_end, chrom
                    region_start = region_end
                    iterator.update(1)
            iterator.close()
        else:
            # read number of lines in bed file
            nlines = sum(1 for _ in exopen(self.bed_file, "r"))
            with exopen(self.bed_file, "r") as bed:
                for line in tqdm(bed, total=nlines, desc="Processing BED regions"):
                    ret = line.strip().split("\t")
                    chrom, start, end = ret[:3]
                    chrom = chrom.replace("chr", "")
                    try:
                        region_id = ret[3]
                    except IndexError:
                        region_id = f"{chrom}_{start}_{end}"
                    chunk_start = int(start)
                    chunk_end = min(chunk_start + self.chunk_size, int(end))
                    while chunk_start < int(end):
                        yield chrom, chunk_start, chunk_end, region_id
                        chunk_start = chunk_end
                        chunk_end = min(chunk_start + self.chunk_size, int(end))


class CMWriterUnpacker:
    """
    Wrapper for file handles to provide consistent write and close interface.

    This class handles both context managers and direct file handles, providing
    a unified interface for writing data and closing files. It also prevents
    closing stdout to avoid unexpected behavior.

    Args:
        cm: A context manager or file handle object

    Attributes:
        cm: The original context manager
        handle: The unpacked file handle
    """

    def __init__(self, cm):
        self.cm = cm
        try:
            self.handle = self.cm.__enter__()
        except AttributeError:
            self.handle = self.cm

    def write(self, data):
        """Write data to the file handle."""
        self.handle.write(data)

    def close(self):
        """Close the file handle, avoiding closure of stdout."""
        try:
            if self.handle != sys.stdout:
                self.handle.close()
        except AttributeError:
            self.cm.__exit__(None, None, None)


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
        njobs (int, optional): Number of threads to use for input bam file decompression. If negative, uses
            (number of CPUs + njobs). Default: 1

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
        >>> # Process specific chromosomes
        >>> wps = WPS(valid_chroms={'1', '2', '3', 'X', 'Y'})
        >>> wps.run(bamfile='sample.bam', out_filepath='chr_wps.tsv')

    Note:
        - Automatically filters duplicate, QC-failed, and unmapped reads
        - Handles both paired-end and single-end sequencing data
        - Supports downsampling for high-coverage samples
        - Can write to separate files per chromosome or target region using placeholders
    """

    def __init__(
        self,
        bed_file=None,
        protection_size=120,
        min_insert_size=None,
        max_insert_size=None,
        valid_chroms=set(map(str, list(range(1, 23)) + ["X", "Y"])),
        chunk_size=1e8,
        njobs=1,
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
        self.njobs = njobs
        if self.njobs < 0:
            self.njobs = joblib.cpu_count() + self.njobs
        self.roi_generator = ROIGenerator(bed_file=self.bed_file)

    def __call__(self, *args, **kwargs):
        return self.run(*args, **kwargs)

    def run(
        self,
        bamfile,
        out_filepath=None,
        downsample_ratio=None,
        compute_coverage=False,
        verbose_output=False,
        add_header=False,
    ):
        """
        Calculate Window Protection Score for all regions and write to file.

        Processes the BAM file to compute WPS values for each genomic position
        in the specified regions (or entire genome). Results are written to a
        tab-separated output file.

        Args:
            bamfile (str): Path to input BAM file (must be sorted and indexed)
            out_filepath (str): Path to output TSV file, or stdout.
                If it contains a formatting substring {target} or {chrom}, it will be used to create per-target  or
                per-chromosome files. Default: None (stdout)
            downsample_ratio (float, optional): Fraction of reads to randomly keep
                (0.0 to 1.0). Useful for high-coverage samples. Default: None (no downsampling)
            compute_coverage (bool, optional): Whether to compute and include base coverage
            verbose_output (bool, optional): Whether to include detailed counts
            add_header (bool, optional): Whether to add header to the output

        Returns:
            None: Results are written directly to out_filepath

        Output Format:
            Tab-separated file with columns:
                - chromosome: Chromosome name (without 'chr' prefix)
                - start: Start position (0-based)
                - end: End position (1-based, start + 1)
                - base read coverage (if compute_coverage=True)
                - outside: Count of fragments spanning the protection window (if verbose_output=True)
                - inside: Count of fragment endpoints in protection window (if verbose_output=True)
                - wps: Window Protection Score (outside - inside)

        Raises:
            FileNotFoundError: If bamfile does not exist
            ValueError: If downsample_ratio is not between 0 and 1

        Example:
            >>> wps = WPS()
            >>> wps.run(bamfile='input.bam', out_filepath='output.tsv')
            >>> # With downsampling
            >>> wps.run(bamfile='input.bam', out_filepath='output.tsv',
            ...         downsample_ratio=0.5)
            >>> # Creating separate files per chromosome
            >>> wps.run(bamfile='input.bam', out_filepath='wps_{chrom}.tsv')
            >>> # Creating separate files per target region
            >>> wps.run(bamfile='input.bam', out_filepath='wps_{target}.tsv.gz')
        """
        if out_filepath is None:
            out_filepath = "stdout"
        input_file = pysam.Samfile(bamfile, "rb", threads=self.njobs)
        prefix = (
            "chr" if any(r.startswith("chr") for r in input_file.references) else ""
        )
        use_partial_writer = "{target}" in out_filepath or "{chrom}" in out_filepath
        partial_writers = dict()
        total_outfile = None
        header_added = set()
        if not use_partial_writer:
            total_outfile = CMWriterUnpacker(exopen(out_filepath, "w"))
            try:
                total_outfile = total_outfile.__enter__()
            except AttributeError:
                pass

        for chrom, start, end, region_id in self.roi_generator.regions(
            bam_file=bamfile
        ):
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
                        (lseq < self.min_insert_size) or (lseq > self.max_insert_size)
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

                outside_cum = np.cumsum(outside)[:-2]

                # Fragments whose endpoints fall inside windows
                all_ends = np.concatenate([starts, ends]) - regionStart
                left = np.clip(all_ends - self.protection_size + 2, 0, n + 1)
                right = np.clip(all_ends + self.protection_size + 1, 0, n + 1)
                inside = np.zeros(n + 2, dtype=int)
                np.add.at(inside, left, 1)
                np.add.at(inside, right, -1)
                inside_cum = np.cumsum(inside)[:-2]

                wps = outside_cum - inside_cum
                coverage = None
                if compute_coverage:
                    coverage = np.zeros(n + 1, dtype=int)
                    np.add.at(
                        coverage,
                        np.clip(starts - regionStart, 0, n),
                        1,
                    )
                    np.add.at(
                        coverage,
                        np.clip(ends - regionStart + 1, 0, n),
                        -1,
                    )
                    coverage = np.cumsum(coverage)[:-1]
            else:
                outside_cum = np.zeros(n, dtype=int)
                inside_cum = np.zeros(n, dtype=int)
                wps = np.zeros(n, dtype=int)
                coverage = None
                if compute_coverage:
                    coverage = np.zeros(n, dtype=int)
            partial_outfile = None
            if use_partial_writer:
                formatted_out_filepath = out_filepath.format(
                    target=region_id,
                    chrom=chrom,
                )
                partial_outfile = partial_writers.get(
                    formatted_out_filepath,
                    CMWriterUnpacker(
                        exopen(formatted_out_filepath, "w", use_pigz=False)
                    ),
                )
                partial_writers[formatted_out_filepath] = partial_outfile
                try:
                    partial_outfile = partial_outfile.__enter__()
                except AttributeError:
                    pass

            outfile = partial_outfile if use_partial_writer else total_outfile
            st = np.arange(regionStart, regionEnd + 1)
            en = st + 1  # add end coordinate
            df = pd.DataFrame({"#chrom": chrom, "start": st, "end": en})
            if compute_coverage:
                df["coverage"] = coverage
            if verbose_output:
                df["outside"] = outside_cum
                df["inside"] = inside_cum
            df["wps"] = wps
            header = False
            if add_header and outfile not in header_added:
                header_added.add(outfile)
                header = True
            df.to_csv(outfile, sep="\t", header=header, index=False)

        if use_partial_writer:
            for writer in partial_writers.values():
                writer.close()
        else:
            total_outfile.close()


def main():
    """Main entry point for the command-line interface."""
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        help="Input BAM file",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        help="The output file path for WPS results. If not provided, results will be printed to stdout.",
        required=False,
    )
    parser.add_argument(
        "-r",
        "--regions",
        dest="regions",
        help="BED file with regions of interest (default: whole genome)",
        default=None,
    )
    parser.add_argument(
        "-w",
        "--protection",
        dest="protection",
        help="Base pair protection window (default: 120)",
        default=120,
        type=int,
    )
    parser.add_argument(
        "--min-insert-size",
        dest="min_insert_size",
        help="Minimum read length threshold to consider (Optional)",
        default=None,
        type=int,
    )
    parser.add_argument(
        "--max-insert-size",
        dest="max_insert_size",
        help="Minimum read length threshold to consider (Optional)",
        default=None,
        type=int,
    )
    parser.add_argument(
        "--downsample",
        dest="downsample",
        help="Ratio to down sample reads (default OFF)",
        default=None,
        type=float,
    )
    parser.add_argument(
        "--chunk-size",
        dest="chunk_size",
        help="Chunk size for processing in pieces, in case of low memory (default 1e8)",
        default=1e8,
        type=int,
    )
    parser.add_argument(
        "--valid-chroms",
        dest="valid_chroms",
        help="Comma-separated list of valid chromosomes to include (e.g., '1,2,3,X,Y') or 'canonical' for chromosomes 1-22, X, Y. Optional",
        default=None,
        type=str,
    )
    parser.add_argument(
        "--compute-coverage",
        dest="compute_coverage",
        help="If provided, output will include base read coverage as the 4th column.",
        action="store_true",
    )
    parser.add_argument(
        "--verbose-output",
        dest="verbose_output",
        help="If provided, output will include separate counts for 'outside' and 'inside' along with WPS.",
        action="store_true",
    )
    parser.add_argument(
        "--add-header",
        dest="add_header",
        help="If provided, output files will include a header line.",
        action="store_true",
    )
    parser.add_argument(
        "--njobs",
        dest="njobs",
        help="Number of threads to use for input bam file decompression. If negative, uses (number of CPUs + njobs). Default: 1",
        default=1,
        type=int,
    )
    args = parser.parse_args()
    valid_chroms = None
    if args.valid_chroms == "canonical":
        valid_chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]
    else:
        valid_chroms = args.valid_chroms.split(",") if args.valid_chroms else None
    optwps = WPS(
        bed_file=args.regions,
        protection_size=args.protection,
        min_insert_size=args.min_insert_size,
        max_insert_size=args.max_insert_size,
        chunk_size=args.chunk_size,
        valid_chroms=valid_chroms,
        njobs=args.njobs,
    )
    optwps.run(
        bamfile=args.input,
        out_filepath=args.output,
        downsample_ratio=args.downsample,
        compute_coverage=args.compute_coverage,
        verbose_output=args.verbose_output,
        add_header=args.add_header,
    )


if __name__ == "__main__":
    main()
