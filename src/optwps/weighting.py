import pysam
from typing import Union
import numpy as np

from .read_validator import ReadValidator
from .utils import ref_aln_length


class FragmentFeatureExtractor:
    """Updates bam file with fragment features such as GC content and mappability."""

    def fit(self, X: Union[str, pysam.AlignmentFile], y=None):
        return self

    def transform(self, read: pysam.AlignedSegment):
        fragment_length = (
            abs(read.template_length)
            if read.is_paired
            else ref_aln_length(read.cigartuples)
        )
        gc_content = self.compute_gc_content(read)
        return [fragment_length, gc_content]

    def compute_gc_content(self, read: pysam.AlignedSegment):
        sequence = read.query_sequence or ""
        read_length = len(sequence)
        if read_length == 0:
            return 0.0
        gc_count = sum(1 for base in sequence if base in "GCgc")
        return gc_count / read_length


class WeightsCalculator:
    """Calculates correction weights for WPS based on fragment features.
    Weights are computed based on the distribution of fragment features in the dataset,
    allowing for correction of biases in WPS calculations.
    The weights are determined by binning the fragment features and calculating the inverse
    of the frequency of fragments in each bin, which can then be applied to adjust WPS
    scores accordingly.
    """

    def __init__(
        self,
        mappability_file: str = None,
        nbins=10,
        subsample=0.05,
        min_insert_size=None,
        max_insert_size=None,
        min_mappability_threshold=0.9,
    ):
        self.subsample = subsample
        self.nbins = nbins
        self.extractor = FragmentFeatureExtractor()
        self.validator = ReadValidator(
            min_insert_size=min_insert_size,
            max_insert_size=max_insert_size,
            mappability_file=mappability_file,
            min_mappability_threshold=min_mappability_threshold,
        )
        self.bin_edges = None
        self.weights = None

    def fit(self, bam: Union[str, pysam.AlignmentFile], y=None):
        close_bam = False
        if isinstance(bam, str):
            bam = pysam.AlignmentFile(bam, "rb")
            close_bam = True
        fragment_features = []
        for read in bam.fetch():
            if not self.validator.valid_read(read, downsample_ratio=self.subsample):
                continue
            fragment_features.append(self.extractor.transform(read))
        if close_bam:
            bam.close()
        if not fragment_features:
            self.bin_edges = None
            self.weights = None
            return self
        fragment_features = np.array(fragment_features)
        histogram, bin_edges = np.histogramdd(fragment_features, bins=self.nbins)
        self.weights = np.ones_like(histogram, dtype=float)
        observed = histogram > 0
        self.weights[observed] = np.mean(histogram[observed]) / histogram[observed]
        self.bin_edges = bin_edges
        return self

    def transform(self, read: pysam.AlignedSegment):
        if self.weights is None:
            return 1.0
        features = self.extractor.transform(read)
        bin_indices = []
        for i in range(len(features)):
            bin_index = (
                np.searchsorted(self.bin_edges[i], features[i], side="right") - 1
            )
            bin_index = np.clip(bin_index, 0, self.weights.shape[i] - 1)
            bin_indices.append(bin_index)
        return self.weights[tuple(bin_indices)]
