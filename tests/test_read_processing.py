import pytest
import numpy as np

from optwps.read_processing import (
    ReadInfo,
    collect_fragment_features,
    collect_fragment_intervals,
    iter_pysam_reads,
    read_info_batches,
    valid_read_info,
    weight_from_features,
)


def _read(**kwargs):
    defaults = dict(
        is_duplicate=False,
        is_qcfail=False,
        is_unmapped=False,
        cigartuples=((0, 100),),
        reference_name="chr1",
        reference_start=10,
        reference_end=110,
        is_paired=True,
        mate_is_unmapped=False,
        rnext=0,
        tid=0,
        is_read1=True,
        is_read2=False,
        pnext=150,
        qlen=100,
        template_length=240,
        pos=10,
        query_sequence="ACGT" * 25,
    )
    defaults.update(kwargs)
    return ReadInfo(**defaults)


def test_iter_pysam_reads_handles_coverage_iterator_artifact():
    class Reads:
        def __init__(self):
            self.values = iter([1, 2])

        def __next__(self):
            try:
                return next(self.values)
            except StopIteration:
                raise ValueError("Firing event 10 with no exception set")

    assert list(iter_pysam_reads(Reads())) == [1, 2]

    class BrokenReads:
        def __next__(self):
            raise ValueError("real failure")

    with pytest.raises(ValueError, match="real failure"):
        list(iter_pysam_reads(BrokenReads()))


def test_read_info_batches_splits_and_flushes_tail():
    batches = list(read_info_batches([_read(pos=i) for i in range(3)], 2))

    assert [[read.pos for read in batch] for batch in batches] == [[0, 1], [2]]


def test_valid_read_info_filter_branches(monkeypatch):
    assert not valid_read_info(_read(is_duplicate=True))
    assert not valid_read_info(_read(cigartuples=((4, 5), (0, 95))))
    assert not valid_read_info(_read(mate_is_unmapped=True))
    assert not valid_read_info(_read(rnext=1))
    assert not valid_read_info(_read(is_read1=False, is_read2=True))
    assert not valid_read_info(_read(template_length=0))
    assert not valid_read_info(_read(template_length=99), min_insert_size=100)
    assert not valid_read_info(_read(template_length=301), max_insert_size=300)

    monkeypatch.setattr("optwps.read_processing.random.random", lambda: 0.9)
    assert not valid_read_info(_read(), downsample_ratio=0.5)

    read2 = _read(is_read1=False, is_read2=True, pnext=10, qlen=50)
    assert valid_read_info(read2, upstream_limit=100)


def test_mappability_filter_branches():
    class Mappability:
        def __init__(self, values=None, raises=False):
            self._values = values
            self.raises = raises
            self.closed = False

        def values(self, *args):
            if self.raises:
                raise RuntimeError("missing chromosome")
            return self._values

    assert not valid_read_info(
        _read(reference_name=None), mappability_file=Mappability([1])
    )
    assert not valid_read_info(_read(), mappability_file=Mappability([]))
    assert not valid_read_info(_read(), mappability_file=Mappability([float("nan")]))
    assert not valid_read_info(_read(), mappability_file=Mappability([0.1]))
    assert not valid_read_info(_read(), mappability_file=Mappability(raises=True))
    assert valid_read_info(_read(), mappability_file=Mappability([1.0]))


def test_collectors_and_weight_lookup():
    read = _read(template_length=200, pos=10, pnext=150, query_sequence="GCGCAA")
    features = collect_fragment_features([[read]])
    starts, ends, weights = collect_fragment_intervals(
        [[read]],
        bin_edges=[[100, 300], [0.0, 1.0]],
        weight_values=np.array([[2.0]]),
        use_weights=True,
    )

    assert features[0][0] == 200
    assert features[0][1] == pytest.approx(4 / 6)
    assert (starts, ends, weights) == ([10], [209], [2.0])
    assert weight_from_features(features[0], None, None) == 1.0
