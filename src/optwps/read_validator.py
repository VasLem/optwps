import pysam

from .read_processing import valid_read_info


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
        return valid_read_info(
            read,
            min_insert_size=self.min_insert_size,
            max_insert_size=self.max_insert_size,
            mappability_file=self.mappability_file,
            min_mappability_threshold=self.min_mappability_threshold,
            upstream_limit=upstream_limit,
            downsample_ratio=downsample_ratio,
        )
