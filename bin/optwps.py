#!/usr/bin/env python
"""
Command-line interface for optwps package.

This script provides a command-line interface for calculating Window Protection Scores
from BAM files. It can be used directly or installed as a console script.
"""

from lib.optwps import WPS


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
        "--outfile",
        dest="outfile",
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
        help="Chunk size for processing in pieces (default 1e6)",
        default=1e6,
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
        "--verbose-output",
        dest="verbose_output",
        help="If provided, output will include separate counts for 'outside' and 'inside' along with WPS.",
        action="store_true",
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
    )
    optwps.run(
        bamfile=args.input,
        out_filepath=args.outfile,
        downsample_ratio=args.downsample,
        verbose_output=args.verbose_output,
    )


if __name__ == "__main__":
    main()
