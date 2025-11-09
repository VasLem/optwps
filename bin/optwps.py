#!/usr/bin/env python
"""
Command-line interface for optwps package.

This script provides a command-line interface for calculating Window Protection Scores
from BAM files. It can be used directly or installed as a console script.
"""

import logging
from lib.optwps import WPS, LOGGER


def main():
    """Main entry point for the command-line interface."""
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        help="Use regions transcript file (def transcriptAnno.tsv)",
        default="transcriptAnno.tsv",
    )
    parser.add_argument(
        "-m",
        "--merged",
        dest="merged",
        help="Assume reads are merged (default Off)",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "-w",
        "--protection",
        dest="protection",
        help="Base pair protection assumed for elements (default 120)",
        default=120,
        type="int",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        dest="outfile",
        help="Outfile",
    )  # reserve atleast 6 digits
    parser.add_argument(
        "-e",
        "--empty",
        dest="empty",
        help="Keep files of empty blocks (def Off)",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--minInsert",
        dest="minInsSize",
        help="Minimum read length threshold to consider (def None)",
        default=-1,
        type="int",
    )
    parser.add_argument(
        "--maxInsert",
        dest="maxInsSize",
        help="Minimum read length threshold to consider (def None)",
        default=-1,
        type="int",
    )
    parser.add_argument(
        "--max_length",
        dest="max_length",
        help="Assumed maximum insert size (default 1000)",
        default=1000,
        type="int",
    )
    parser.add_argument(
        "--downsample",
        dest="downsample",
        help="Ratio to down sample reads (default OFF)",
        default=None,
        type="float",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        help="Turn debug output on",
        action="store_true",
    )
    args = parser.parse_args()
    if args.verbose:
        LOGGER.setLevel(logging.DEBUG)
    else:
        LOGGER.setLevel(logging.INFO)
    optwps = WPS(
        protection_size=args.protection,
        min_insert_size=args.minInsSize if args.minInsSize > 0 else None,
        max_insert_size=args.maxInsSize if args.maxInsSize > 0 else None,
    )
    optwps.run(
        bamfile=args.input, out_filepath=args.outfile, downsample_ratio=args.downsample
    )


if __name__ == "__main__":
    main()
