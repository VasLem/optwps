import random
from dataclasses import dataclass

import joblib
import numpy as np

from .utils import is_soft_clipped, ref_aln_length


@dataclass
class ReadInfo:
    is_duplicate: bool
    is_qcfail: bool
    is_unmapped: bool
    cigartuples: tuple
    reference_name: str
    reference_start: int
    reference_end: int
    is_paired: bool
    mate_is_unmapped: bool
    rnext: int
    tid: int
    is_read1: bool
    is_read2: bool
    pnext: int
    qlen: int
    template_length: int
    pos: int
    query_sequence: str


def read_info(read):
    try:
        reference_name = read.reference_name
    except ValueError:
        reference_name = None
    return ReadInfo(
        read.is_duplicate,
        read.is_qcfail,
        read.is_unmapped,
        tuple(read.cigartuples or ()),
        reference_name,
        read.reference_start,
        read.reference_end,
        read.is_paired,
        read.mate_is_unmapped,
        read.rnext,
        read.tid,
        read.is_read1,
        read.is_read2,
        read.pnext,
        read.qlen,
        read.template_length,
        read.pos,
        read.query_sequence or "",
    )


def read_info_batches(reads, buffer_size):
    batch = []
    for read in reads:
        batch.append(read_info(read))
        if len(batch) >= buffer_size:
            yield batch
            batch = []
    if batch:
        yield batch


def iter_pysam_reads(reads):
    while True:
        try:
            yield next(reads)
        except StopIteration:
            return
        except ValueError as error:
            if "Firing event 10 with no exception set" not in str(error):
                raise
            return


def valid_read_info(
    read,
    min_insert_size=None,
    max_insert_size=None,
    mappability_file=None,
    min_mappability_threshold=0.9,
    upstream_limit=None,
    downsample_ratio=None,
):
    if read.is_duplicate or read.is_qcfail or read.is_unmapped:
        return False
    if is_soft_clipped(read.cigartuples):
        return False
    if mappability_file is not None:
        if read.reference_name is None or read.reference_end is None:
            return False
        try:
            mappability = mappability_file.values(
                read.reference_name, read.reference_start, read.reference_end
            )
            if not mappability or any(value != value for value in mappability):
                return False
            if min_mappability_threshold is not None and (
                sum(mappability) / len(mappability) < min_mappability_threshold
            ):
                return False
        except RuntimeError:
            return False
    if read.is_paired:
        if read.mate_is_unmapped or read.rnext != read.tid:
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
        lseq = abs(read.template_length)
        if lseq == 0:
            return False
    else:
        lseq = ref_aln_length(read.cigartuples)
    if downsample_ratio is not None and random.random() >= downsample_ratio:
        return False
    if min_insert_size is not None and lseq < min_insert_size:
        return False
    if max_insert_size is not None and lseq > max_insert_size:
        return False
    return True


def fragment_interval(read):
    if read.is_paired:
        lseq = abs(read.template_length)
        rstart = min(read.pos, read.pnext)
    else:
        lseq = ref_aln_length(read.cigartuples)
        rstart = read.pos
    return rstart, rstart + lseq - 1


def fragment_features(read):
    fragment_length = (
        abs(read.template_length)
        if read.is_paired
        else ref_aln_length(read.cigartuples)
    )
    sequence = read.query_sequence
    if not sequence:
        return [fragment_length, 0.0]
    gc_count = sum(1 for base in sequence if base in "GCgc")
    return [fragment_length, gc_count / len(sequence)]


def weight_from_features(features, bin_edges, weights):
    if weights is None:
        return 1.0
    bin_indices = []
    for i, feature in enumerate(features):
        bin_index = np.searchsorted(bin_edges[i], feature, side="right") - 1
        bin_indices.append(np.clip(bin_index, 0, weights.shape[i] - 1))
    return weights[tuple(bin_indices)]


def _open_mappability(path):
    if path is None:
        return None
    import pyBigWig

    return pyBigWig.open(path)


def valid_fragment_features_batch(
    reads,
    min_insert_size=None,
    max_insert_size=None,
    mappability_path=None,
    min_mappability_threshold=0.9,
    downsample_ratio=None,
):
    mappability_file = _open_mappability(mappability_path)
    try:
        return [
            fragment_features(read)
            for read in reads
            if valid_read_info(
                read,
                min_insert_size=min_insert_size,
                max_insert_size=max_insert_size,
                mappability_file=mappability_file,
                min_mappability_threshold=min_mappability_threshold,
                downsample_ratio=downsample_ratio,
            )
        ]
    finally:
        if mappability_file is not None:
            mappability_file.close()


def valid_fragment_intervals_batch(
    reads,
    min_insert_size=None,
    max_insert_size=None,
    mappability_path=None,
    min_mappability_threshold=0.9,
    upstream_limit=None,
    downsample_ratio=None,
    bin_edges=None,
    weight_values=None,
    use_weights=False,
):
    mappability_file = _open_mappability(mappability_path)
    starts, ends, weights = [], [], []
    try:
        for read in reads:
            if not valid_read_info(
                read,
                min_insert_size=min_insert_size,
                max_insert_size=max_insert_size,
                mappability_file=mappability_file,
                min_mappability_threshold=min_mappability_threshold,
                upstream_limit=upstream_limit,
                downsample_ratio=downsample_ratio,
            ):
                continue
            start, end = fragment_interval(read)
            starts.append(start)
            ends.append(end)
            if use_weights:
                features = fragment_features(read)
                weights.append(weight_from_features(features, bin_edges, weight_values))
    finally:
        if mappability_file is not None:
            mappability_file.close()
    return starts, ends, weights


def process_read_batches(batches, njobs, batch_func, max_queued_batches=None, **kwargs):
    if njobs == 1:
        for batch in batches:
            yield batch_func(batch, **kwargs)
    else:
        max_queued_batches = max_queued_batches or max(1, njobs * 2)
        queued = []
        for batch in batches:
            queued.append(batch)
            if len(queued) >= max_queued_batches:
                yield from joblib.Parallel(n_jobs=njobs)(
                    joblib.delayed(batch_func)(queued_batch, **kwargs)
                    for queued_batch in queued
                )
                queued = []
        if queued:
            yield from joblib.Parallel(n_jobs=njobs)(
                joblib.delayed(batch_func)(queued_batch, **kwargs)
                for queued_batch in queued
            )


def collect_fragment_features(
    batches,
    min_insert_size=None,
    max_insert_size=None,
    mappability_path=None,
    min_mappability_threshold=0.9,
    downsample_ratio=None,
    njobs=1,
):
    return [
        feature
        for batch in process_read_batches(
            batches,
            njobs,
            valid_fragment_features_batch,
            min_insert_size=min_insert_size,
            max_insert_size=max_insert_size,
            mappability_path=mappability_path,
            min_mappability_threshold=min_mappability_threshold,
            downsample_ratio=downsample_ratio,
        )
        for feature in batch
    ]


def collect_fragment_intervals(
    batches,
    min_insert_size=None,
    max_insert_size=None,
    mappability_path=None,
    min_mappability_threshold=0.9,
    upstream_limit=None,
    downsample_ratio=None,
    bin_edges=None,
    weight_values=None,
    use_weights=False,
    njobs=1,
):
    starts, ends, weights = [], [], []
    for read_starts, read_ends, read_weights in process_read_batches(
        batches,
        njobs,
        valid_fragment_intervals_batch,
        min_insert_size=min_insert_size,
        max_insert_size=max_insert_size,
        mappability_path=mappability_path,
        min_mappability_threshold=min_mappability_threshold,
        upstream_limit=upstream_limit,
        downsample_ratio=downsample_ratio,
        bin_edges=bin_edges,
        weight_values=weight_values,
        use_weights=use_weights,
    ):
        starts.extend(read_starts)
        ends.extend(read_ends)
        weights.extend(read_weights)
    return starts, ends, weights
