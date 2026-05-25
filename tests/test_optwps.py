import pysam
import pytest

from optwps import exopen, is_soft_clipped, ref_aln_length
from optwps.read_processing import iter_pysam_reads
from bx.intervals import Intersecter, Interval


def old_implementation(
    bed_file,
    bamfile,
    out_filepath,
    protection_size,
    valid_chroms,
    min_insert_size,
    max_insert_size,
    njobs=1,
):
    protection = protection_size // 2
    with exopen(bed_file, "r") as infile:
        for cnt, line in enumerate(infile):
            chrom, start, end = line.split(
                "\t"
            )  # positions should be 0-based and end non-inclusive
            chrom = chrom.replace("chr", "")
            if chrom not in valid_chroms:
                continue

            regionStart, regionEnd = int(start), int(end)

            filteredReads = Intersecter()

            input_file = (
                pysam.Samfile(bamfile, "rb")
                if njobs == 1
                else pysam.Samfile(bamfile, "rb", threads=njobs)
            )
            prefix = ""
            for tchrom in input_file.references:
                if tchrom.startswith("chr"):
                    prefix = "chr"
                    break

            for read in iter_pysam_reads(
                input_file.fetch(
                    prefix + chrom,
                    max(0, regionStart - protection - 1),
                    regionEnd + protection + 1,
                )
            ):
                if read.is_duplicate or read.is_qcfail or read.is_unmapped:
                    continue
                if is_soft_clipped(read.cigar):
                    continue

                if read.is_paired:
                    if read.mate_is_unmapped:
                        continue
                    if read.rnext != read.tid:
                        continue
                    if read.is_read1 or (
                        read.is_read2
                        and read.pnext + read.qlen < regionStart - protection - 1
                    ):
                        if read.isize == 0:
                            continue
                        rstart = min(read.pos, read.pnext) + 1  # 1-based
                        lseq = abs(read.isize)
                        rend = rstart + lseq - 1  # end included
                        if min_insert_size != None and (
                            (lseq < min_insert_size) or (lseq > max_insert_size)
                        ):
                            continue

                        filteredReads.add_interval(Interval(rstart, rend))
                else:
                    rstart = read.pos + 1  # 1-based
                    lseq = ref_aln_length(read.cigar)
                    rend = rstart + lseq - 1  # end included
                    if min_insert_size != None and (
                        (lseq < min_insert_size) or (lseq > max_insert_size)
                    ):
                        continue

                    filteredReads.add_interval(Interval(rstart, rend))
            with exopen(out_filepath, "a" if cnt > 0 else "w") as out:
                outLines = []
                for pos in range(regionStart, regionEnd + 1):
                    rstart, rend = pos - protection, pos + protection
                    gcount, bcount = 0, 0
                    for read in filteredReads.find(rstart, rend):
                        if (read.start > rstart) or (read.end < rend):
                            bcount += 1
                        else:
                            gcount += 1

                    outLines.append(
                        "%s\t%d\t%d\t%d\t%d\t%d\n"
                        % (chrom, pos, pos + 1, gcount, bcount, gcount - bcount)
                    )

                for line in outLines:
                    out.write(line)
                out.close()


def test_optwps_does_the_same_as_old_version_paired_end(
    make_test_bed_file, make_test_bam_file_paired, tmp_path
):
    from optwps import WPS

    maker = WPS(
        bed_file=str(make_test_bed_file),
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
    )
    tmp_new_output = tmp_path / "new_wps_output.tsv"
    tmp_new_output = str(tmp_new_output)
    maker.run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=tmp_new_output,
    )
    tmp_old_output = tmp_path / "old_wps_output.tsv"
    tmp_old_output = str(tmp_old_output)
    old_implementation(
        bed_file=str(make_test_bed_file),
        bamfile=str(make_test_bam_file_paired),
        out_filepath=tmp_old_output,
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
        min_insert_size=None,
        max_insert_size=None,
    )
    new_lines = open(tmp_new_output).readlines()
    old_lines = open(tmp_old_output).readlines()
    assert len(new_lines) == len(old_lines)
    for cnt, (nl, ol) in enumerate(zip(new_lines, old_lines)):
        if nl.split("\t")[-1] != ol.split("\t")[-1]:
            raise AssertionError(
                "Mismatch in line {}:\nNew: {}\nOld: {}".format(cnt, nl, ol)
            )


def test_optwps_does_the_same_as_old_version_single_end(
    make_test_bed_file, make_test_bam_file_single, tmp_path
):
    from optwps import WPS

    maker = WPS(
        bed_file=str(make_test_bed_file),
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
    )
    tmp_new_output = tmp_path / "new_wps_output_single.tsv"
    tmp_new_output = str(tmp_new_output)
    maker.run(
        bamfile=str(make_test_bam_file_single),
        out_filepath=tmp_new_output,
    )
    tmp_old_output = tmp_path / "old_wps_output_single.tsv"
    tmp_old_output = str(tmp_old_output)
    old_implementation(
        bed_file=str(make_test_bed_file),
        bamfile=str(make_test_bam_file_single),
        out_filepath=tmp_old_output,
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
        min_insert_size=None,
        max_insert_size=None,
    )
    new_lines = open(tmp_new_output).readlines()
    old_lines = open(tmp_old_output).readlines()
    assert len(new_lines) == len(old_lines)
    for cnt, (nl, ol) in enumerate(zip(new_lines, old_lines)):
        if nl.split("\t")[-1] != ol.split("\t")[-1]:
            raise AssertionError(
                "Mismatch in line {}:\nNew: {}\nOld: {}".format(cnt, nl, ol)
            )


def test_optwps_downsampling(
    make_test_bed_file, make_test_bam_file_paired, make_test_bam_file_single, tmp_path
):
    from optwps import WPS

    maker = WPS(
        bed_file=str(make_test_bed_file),
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
    )
    tmp_output = tmp_path / "wps_output_downsampled.tsv"
    tmp_output = str(tmp_output)
    maker.run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=tmp_output,
        downsample_ratio=0.5,
    )
    lines = open(tmp_output).readlines()
    assert len(lines) > 0  # Just check that some output is produced
    maker.run(
        bamfile=str(make_test_bam_file_single),
        out_filepath=tmp_output,
        downsample_ratio=0.5,
    )
    lines = open(tmp_output).readlines()
    assert len(lines) > 0  # Just check that some output is produced


def test_optwps_with_parallelism(
    make_test_bed_file, make_test_bam_file_paired, tmp_path
):
    from optwps import WPS

    serial = WPS(
        bed_file=str(make_test_bed_file),
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
    )
    parallel = WPS(
        bed_file=str(make_test_bed_file),
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
        njobs=2,
        read_buffer_size=2,
    )
    serial_output = str(tmp_path / "wps_output_serial.tsv")
    parallel_output = str(tmp_path / "wps_output_parallel.tsv")
    serial.run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=serial_output,
    )
    parallel.run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=parallel_output,
    )
    assert open(parallel_output).readlines() == open(serial_output).readlines()


def test_optwps_no_bed_file(make_test_bam_file_paired, tmp_path):
    from optwps import WPS

    maker = WPS(
        bed_file=None,
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
    )
    tmp_output = tmp_path / "wps_output_no_bed.tsv"
    tmp_output = str(tmp_output)
    maker.run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=tmp_output,
    )
    lines = open(tmp_output).readlines()
    assert len(lines) > 0  # Just check that some output is produced


def test_gzbedfile_handling(make_test_bed_file, make_test_bam_file_paired, tmp_path):
    import gzip

    # Create a gzipped version of the bed file
    gz_bed_path = tmp_path / "test_regions.bed.gz"
    with open(make_test_bed_file, "rb") as f_in:
        with gzip.open(gz_bed_path, "wb") as f_out:
            f_out.writelines(f_in)
    from optwps import WPS

    maker = WPS(
        bed_file=str(gz_bed_path),
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
    )
    tmp_output = tmp_path / "wps_output_gzbed.tsv.gz"
    tmp_output = str(tmp_output)
    maker.run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=tmp_output,
    )


def test_with_minsize_maxsize(
    make_test_bed_file, make_test_bam_file_paired, make_test_bam_file_single, tmp_path
):
    from optwps import WPS

    maker = WPS(
        bed_file=str(make_test_bed_file),
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
        min_insert_size=400,
        max_insert_size=800,
    )
    tmp_output = tmp_path / "wps_output_minsize_maxsize.tsv"
    tmp_output = str(tmp_output)
    maker.run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=tmp_output,
    )
    lines = open(tmp_output).readlines()
    assert len(lines) > 0  # Just check that some output is produced
    maker.run(
        bamfile=str(make_test_bam_file_single),
        out_filepath=tmp_output,
    )
    lines = open(tmp_output).readlines()
    assert len(lines) > 0  # Just check that some output is produced


def test_with_coverage_computation(
    make_test_bed_file, make_test_bam_file_paired, tmp_path
):
    from optwps import WPS

    maker = WPS(
        bed_file=str(make_test_bed_file),
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
    )
    tmp_output = tmp_path / "wps_output_with_coverage.tsv"
    tmp_output = str(tmp_output)
    maker.run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=tmp_output,
        compute_coverage=True,
    )
    lines = open(tmp_output).readlines()
    assert len(lines) > 0  # Just check that some output is produced
    # Check that coverage column is present
    assert len(lines[0].strip().split("\t")) == 5  # chrom, start, end, coverage, wps


def test_printed_to_stdout(
    make_test_bed_file, make_test_bam_file_paired, tmp_path, capsys
):
    from optwps import WPS

    maker = WPS(
        bed_file=str(make_test_bed_file),
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
    )
    maker.run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=None,
        verbose_output=True,
    )
    captured = capsys.readouterr()
    lines = captured.out.strip().split("\n")
    assert len(lines) > 0  # Just check that some output is produced
    # Check that verbose output columns are present
    for line in lines:
        cols = line.split("\t")
        assert len(cols) == 6  # chrom, pos_start, pos_end, gcount, bcount, wps


def test_with_header(make_test_bed_file, make_test_bam_file_paired, tmp_path, capsys):
    from optwps import WPS

    tmp_output = tmp_path / "wps_output_with_header.tsv"

    maker = WPS(
        bed_file=str(make_test_bed_file),
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
    )
    maker.run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=str(tmp_output),
        add_header=True,
    )
    contents = open(tmp_output, "r").read()
    assert contents.startswith("#chrom\tstart\tend\twps\n")


def _read_verbose_wps(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            chrom, start, end, outside, inside, wps = line.strip().split("\t")
            rows.append(
                (
                    chrom,
                    int(start),
                    int(end),
                    float(outside),
                    float(inside),
                    float(wps),
                )
            )
    return rows


def _make_mappability_bigwig(path, bamfile, low_chroms=None):
    import pyBigWig

    low_chroms = set(low_chroms or [])
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        chroms = list(bam.references)
        lengths = list(bam.lengths)

    bigwig = pyBigWig.open(str(path), "w")
    bigwig.addHeader(list(zip(chroms, lengths)))
    for chrom, length in zip(chroms, lengths):
        value = 0.0 if chrom in low_chroms else 1.0
        bigwig.addEntries([chrom], [0], ends=[length], values=[value])
    bigwig.close()
    return path


def test_mappability_bigwig_filters_low_mappability_reads(
    make_test_bed_file, make_test_bam_file_paired, tmp_path
):
    from optwps import WPS

    mappability_file = _make_mappability_bigwig(
        tmp_path / "mappability.bw",
        make_test_bam_file_paired,
        low_chroms={"chr1"},
    )
    output = str(tmp_path / "mappability_filtered.tsv")
    WPS(
        bed_file=str(make_test_bed_file),
        protection_size=120,
        valid_chroms=set(["1", "2"]),
        mappability_file=str(mappability_file),
        min_mappability=0.9,
    ).run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=output,
        verbose_output=True,
    )

    rows = _read_verbose_wps(output)
    chr1_rows = [row for row in rows if row[0] == "1"]
    chr2_rows = [row for row in rows if row[0] == "2"]

    assert chr1_rows
    assert all(row[3:] == (0.0, 0.0, 0.0) for row in chr1_rows)
    assert any(row[3:] != (0.0, 0.0, 0.0) for row in chr2_rows)


def test_bias_weights_are_applied_to_inside_and_outside_counts(
    make_test_bed_file, make_test_bam_file_paired, tmp_path
):
    from optwps import WPS

    class FixedWeights:
        def fit(self, bam):
            return self

        def transform(self, read):
            return 2.5

    valid_chroms = set(["1", "2", "X", "3", "4", "5"])
    unweighted = WPS(
        bed_file=str(make_test_bed_file),
        protection_size=120,
        valid_chroms=valid_chroms,
    )
    weighted = WPS(
        bed_file=str(make_test_bed_file),
        protection_size=120,
        valid_chroms=valid_chroms,
    )
    weighted.weights_calculator = FixedWeights()

    unweighted_output = str(tmp_path / "unweighted.tsv")
    weighted_output = str(tmp_path / "weighted.tsv")
    unweighted.run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=unweighted_output,
        verbose_output=True,
    )
    weighted.run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=weighted_output,
        verbose_output=True,
    )

    for plain, corrected in zip(
        _read_verbose_wps(unweighted_output), _read_verbose_wps(weighted_output)
    ):
        assert corrected[:3] == plain[:3]
        assert corrected[3] == pytest.approx(plain[3] * 2.5)
        assert corrected[4] == pytest.approx(plain[4] * 2.5)
        assert corrected[5] == pytest.approx(plain[5] * 2.5)


def test_bias_correction_runs_with_estimated_weights(
    make_test_bed_file, make_test_bam_file_paired, tmp_path
):
    from optwps import WPS

    output = str(tmp_path / "bias_corrected.tsv")
    WPS(
        bed_file=str(make_test_bed_file),
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
        correct_for_bias=True,
        bias_bins=2,
        bias_subsample=1.0,
    ).run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=output,
        verbose_output=True,
    )

    assert len(_read_verbose_wps(output)) > 0


def test_parallel_bias_correction_matches_serial(
    make_test_bed_file, make_test_bam_file_paired, tmp_path
):
    from optwps import WPS

    kwargs = dict(
        bed_file=str(make_test_bed_file),
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
        correct_for_bias=True,
        bias_bins=2,
        bias_subsample=1.0,
    )
    serial_output = str(tmp_path / "bias_corrected_serial.tsv")
    parallel_output = str(tmp_path / "bias_corrected_parallel.tsv")
    WPS(**kwargs).run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=serial_output,
        verbose_output=True,
    )
    WPS(**kwargs, njobs=2, read_buffer_size=2).run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=parallel_output,
        verbose_output=True,
    )

    assert _read_verbose_wps(parallel_output) == _read_verbose_wps(serial_output)
