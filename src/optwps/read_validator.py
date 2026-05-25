import random

import pysam

from .utils import is_soft_clipped, ref_aln_length


class ReadValidator:
    """
    Validates reads based on criteria such as duplication, mapping quality,
    insert size, and mappability.
    This class encapsulates the logic for determining whether a read should be
    included in WPS calculations."""

    def __init__(
        self,
        min_insert_size=None,
        max_insert_size=None,
        mappability_file=None,
        min_mappability_threshold=0.9,
    ):
        self.mappability_file = None
        if mappability_file is not None:
            import pyBigWig

            self.mappability_file = pyBigWig.open(mappability_file)
            if not self.mappability_file.isBigWig():
                raise ValueError(
                    "Provided mappability file is not a valid BigWig file."
                )
        self.min_insert_size = min_insert_size
        self.max_insert_size = max_insert_size
        self.min_mappability_threshold = min_mappability_threshold

    def valid_read(
        self,
        read: pysam.AlignedSegment,
        upstream_limit: int = None,
        downsample_ratio: float = None,
    ):
        """Check if a read is valid for WPS calculation based on various criteria."""
        if read.is_duplicate or read.is_qcfail or read.is_unmapped:
            return False
        if is_soft_clipped(read.cigartuples):
            return False
        if self.mappability_file is not None:
            chrom = read.reference_name
            start = read.reference_start
            end = read.reference_end
            try:
                mappability = self.mappability_file.values(chrom, start, end)
                if not mappability:
                    return False
                if any(value != value for value in mappability):
                    return False
                if self.min_mappability_threshold is not None and (
                    sum(mappability) / len(mappability)
                    < self.min_mappability_threshold
                ):
                    return False  # Filter out reads in low mappability regions
            except RuntimeError:
                return (
                    False  # Filter out reads in regions not found in mappability file
                )
        if read.is_paired:
            if read.mate_is_unmapped:
                return False
            if read.rnext != read.tid:
                return False
            if not (
                read.is_read1
                or (
                    read.is_read2
                    and upstream_limit is not None
                    and read.pnext + read.qlen < upstream_limit
                )
            ):
                return False
            lseq = abs(read.isize)
            if lseq == 0:
                return False
            if downsample_ratio is not None and random.random() >= downsample_ratio:
                return False
            if self.min_insert_size is not None and lseq < self.min_insert_size:
                return False
            if self.max_insert_size is not None and lseq > self.max_insert_size:
                return False
            return True
        else:
            if downsample_ratio is not None and random.random() >= downsample_ratio:
                return False
            lseq = ref_aln_length(read.cigartuples)
            if self.min_insert_size is not None and lseq < self.min_insert_size:
                return False
            if self.max_insert_size is not None and lseq > self.max_insert_size:
                return False
        return True
