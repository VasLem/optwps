from pytest import fixture
import pysam


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
