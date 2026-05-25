import pysam
from typing import Union
import numpy as np
import joblib

from .read_processing import (
    collect_fragment_features,
    fragment_features,
    iter_pysam_reads,
    read_info,
    weight_from_features,
)


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
        prior_count=20.0,
        min_weight=0.2,
        max_weight=5.0,
        njobs=1,
        read_buffer_size=10000,
    ):
        self.subsample = subsample
        self.nbins = nbins
        self.mappability_file = mappability_file
        self.min_insert_size = min_insert_size
        self.max_insert_size = max_insert_size
        self.min_mappability_threshold = min_mappability_threshold
        self.prior_count = prior_count
        self.min_weight = min_weight
        self.max_weight = max_weight
        self.njobs = joblib.cpu_count() + njobs if njobs < 0 else njobs
        self.njobs = max(1, self.njobs)
        self.read_buffer_size = read_buffer_size
        self.bin_edges = None
        self.weights = None

    def fit(self, bam: Union[str, pysam.AlignmentFile]):
        close_bam = False
        if isinstance(bam, str):
            bam = (
                pysam.AlignmentFile(bam, "rb")
                if self.njobs == 1
                else pysam.AlignmentFile(bam, "rb", threads=self.njobs)
            )
            close_bam = True
        read_batches = []
        batch = []
        for read in iter_pysam_reads(bam.fetch()):
            batch.append(read_info(read))
            if len(batch) >= self.read_buffer_size:
                read_batches.append(batch)
                batch = []
        if batch:
            read_batches.append(batch)
        features = collect_fragment_features(
            read_batches,
            min_insert_size=self.min_insert_size,
            max_insert_size=self.max_insert_size,
            mappability_path=self.mappability_file,
            min_mappability_threshold=self.min_mappability_threshold,
            downsample_ratio=self.subsample,
            njobs=self.njobs,
            read_buffer_size=self.read_buffer_size,
        )
        if close_bam:
            bam.close()
        if not features:
            self.bin_edges = None
            self.weights = None
            return self
        histogram, bin_edges = np.histogramdd(np.array(features), bins=self.nbins)
        self.weights = np.ones_like(histogram, dtype=float)
        observed = histogram > 0
        mean_count = np.mean(histogram[observed])
        self.weights[observed] = (mean_count + self.prior_count) / (
            histogram[observed] + self.prior_count
        )
        if self.min_weight is not None or self.max_weight is not None:
            self.weights[observed] = np.clip(
                self.weights[observed],
                -np.inf if self.min_weight is None else self.min_weight,
                np.inf if self.max_weight is None else self.max_weight,
            )
        self.bin_edges = bin_edges
        return self

    def transform_features(self, features):
        return weight_from_features(features, self.bin_edges, self.weights)

    def transform(self, read: pysam.AlignedSegment):
        if self.weights is None:
            return 1.0
        return self.transform_features(fragment_features(read))
