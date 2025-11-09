import os
import tempfile

import pytest
from pytest import fixture
import pysam


from lib.utils import exopen, is_soft_clipped, ref_aln_length
from collections import defaultdict
from bx.intervals import Intersecter, Interval


@fixture
def make_test_bed_file(tmp_path):
    bed_content = """chr1\t1000\t2000
chr2\t1500\t2500
chr3\t0\t100
chr4\t50\t150
chr5\t3000\t4000
chrM\t200\t800
chrX\t500\t1500
chrUn_gl000220\t100\t300
"""
    bed_file = tmp_path / "test_regions.bed"
    bed_file.write_text(bed_content)
    return bed_file


def _create_read(name, ref_id, pos, length, flag, mate_pos=None, isize=None):
    """Helper function to create a properly configured read.

    For paired-end reads, provide mate_pos and isize.
    For single-end reads, leave mate_pos and isize as None.
    """
    read = pysam.AlignedSegment()
    read.query_name = name
    read.reference_id = ref_id
    read.reference_start = pos
    read.cigar = ((0, length),)
    read.mapping_quality = 60
    read.query_sequence = "A" * length
    read.query_qualities = pysam.qualitystring_to_array("I" * length)
    read.flag = flag

    # Set mate information only for paired-end reads
    if mate_pos is not None:
        read.next_reference_id = ref_id
        read.next_reference_start = mate_pos
        read.template_length = isize if isize is not None else 0

    return read


@fixture
def make_test_bam_file_paired(tmp_path, make_test_bed_file):
    bam_path = tmp_path / "test_reads.bam"
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [
            {"LN": 5000, "SN": "chr1"},
            {"LN": 5000, "SN": "chr2"},
            {"LN": 5000, "SN": "chrX"},
            {"LN": 1000, "SN": "chrM"},
            {"LN": 300, "SN": "chrUn_gl000220"},
            {"LN": 5000, "SN": "chr3"},
            {"LN": 5000, "SN": "chr4"},
            {"LN": 5000, "SN": "chr5"},
        ],
    }

    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        for line in make_test_bed_file.read_text().strip().split("\n"):
            chrom, start, end = line.split("\t")[:3]
            start, end = int(start), int(end)
            ref_id = outf.get_tid(chrom)

            # Pair 1: Completely spans the region
            pos1, pos2, length = max(1, start - 10), end + 5, 100
            isize = (pos2 + length) - pos1
            name = f"pair1_{chrom}_{start}"
            outf.write(_create_read(name, ref_id, pos1, length, 99, pos2, isize))
            outf.write(_create_read(name, ref_id, pos2, length, 147, pos1, -isize))

            # Pair 2: Partially outside the region
            pos1, pos2, length = max(1, end - 50), end + 30, 80
            isize = (pos2 + length) - pos1
            name = f"pair2_{chrom}_{start}"
            outf.write(_create_read(name, ref_id, pos1, length, 99, pos2, isize))
            outf.write(_create_read(name, ref_id, pos2, length, 147, pos1, -isize))

    # Sort and index the BAM file
    sorted_bam_path = str(tmp_path / "test_reads_sorted.bam")
    pysam.sort("-o", sorted_bam_path, str(bam_path))
    pysam.index(sorted_bam_path)
    return sorted_bam_path


@fixture
def make_test_bam_file_single(tmp_path, make_test_bed_file):
    """Create a single-end BAM file for testing."""
    bam_path = tmp_path / "test_reads_single.bam"
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [
            {"LN": 5000, "SN": "chr1"},
            {"LN": 5000, "SN": "chr2"},
            {"LN": 5000, "SN": "chrX"},
            {"LN": 1000, "SN": "chrM"},
            {"LN": 300, "SN": "chrUn_gl000220"},
            {"LN": 5000, "SN": "chr3"},
            {"LN": 5000, "SN": "chr4"},
            {"LN": 5000, "SN": "chr5"},
        ],
    }

    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        for line in make_test_bed_file.read_text().strip().split("\n"):
            chrom, start, end = line.split("\t")[:3]
            start, end = int(start), int(end)
            ref_id = outf.get_tid(chrom)

            # Read 1: Completely spans the region (forward strand)
            pos, length = max(1, start - 10), end - start + 20
            outf.write(
                _create_read(f"read1_{chrom}_{start}", ref_id, pos, length, flag=0)
            )

            # Read 2: Partially outside the region (reverse strand)
            pos, length = end - 30, 60
            outf.write(
                _create_read(f"read2_{chrom}_{start}", ref_id, pos, length, flag=16)
            )

    # Sort and index the BAM file
    sorted_bam_path = str(tmp_path / "test_reads_single_sorted.bam")
    pysam.sort("-o", sorted_bam_path, str(bam_path))
    pysam.index(sorted_bam_path)
    return sorted_bam_path


def old_implementation(
    bed_file,
    bamfile,
    out_filepath,
    protection_size,
    valid_chroms,
    min_insert_size,
    max_insert_size,
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

            input_file = pysam.Samfile(bamfile, "rb")
            prefix = ""
            for tchrom in input_file.references:
                if tchrom.startswith("chr"):
                    prefix = "chr"
                    break

            for read in input_file.fetch(
                prefix + chrom,
                max(0, regionStart - protection - 1),
                regionEnd + protection + 1,
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
    from lib.optwps import WPS

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
    from lib.optwps import WPS

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
    from lib.optwps import WPS

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


def test_optwps_no_bed_file(make_test_bam_file_paired, tmp_path):
    from lib.optwps import WPS

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
    from lib.utils import exopen

    # Create a gzipped version of the bed file
    gz_bed_path = tmp_path / "test_regions.bed.gz"
    with open(make_test_bed_file, "rb") as f_in:
        with gzip.open(gz_bed_path, "wb") as f_out:
            f_out.writelines(f_in)
    from lib.optwps import WPS

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
    from lib.optwps import WPS

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
