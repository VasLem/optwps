from optwps import WPS


def test_multiple_outputs_with_chrom_placeholder(
    make_test_bed_file, make_test_bam_file_paired, tmp_path
):
    """Test that {chrom} placeholder creates separate files per chromosome."""
    wps = WPS(
        bed_file=str(make_test_bed_file),
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
    )

    output_template = str(tmp_path / "wps_{chrom}.tsv")
    wps.run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=output_template,
    )

    # Check that separate files were created for each chromosome in the bed file
    expected_chroms = ["1", "2", "3", "4", "5", "X"]
    for chrom in expected_chroms:
        output_file = tmp_path / f"wps_{chrom}.tsv"
        assert (
            output_file.exists()
        ), f"Expected output file for chromosome {chrom} not found"

        # Verify that the file contains only data for that chromosome
        with open(output_file) as f:
            lines = f.readlines()
            assert len(lines) > 0, f"Output file for chromosome {chrom} is empty"
            for line in lines:
                cols = line.strip().split("\t")
                assert (
                    cols[0] == chrom
                ), f"Found data for wrong chromosome in {output_file}"


def test_multiple_outputs_with_target_placeholder(
    make_test_bed_file, make_test_bam_file_paired, tmp_path
):
    """Test that {target} placeholder creates separate files per target region."""
    wps = WPS(
        bed_file=str(make_test_bed_file),
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
    )

    output_template = str(tmp_path / "wps_{target}.tsv")
    wps.run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=output_template,
    )

    # Read the bed file to know expected targets
    with open(make_test_bed_file) as f:
        bed_lines = f.readlines()

    # Check that files were created for each target
    created_files = list(tmp_path.glob("wps_*.tsv"))
    assert len(created_files) > 0, "No output files were created"

    # Verify each target file has correct region data
    for line in bed_lines:
        chrom, start, end = line.strip().split()[:3]
        chrom = chrom.replace("chr", "")
        if chrom not in ["1", "2", "X", "3", "4", "5"]:
            continue

        target_name = f"{chrom}_{start}_{end}"
        target_file = tmp_path / f"wps_{target_name}.tsv"
        assert target_file.exists(), f"Expected target file {target_file} not found"

        # Verify content is for the correct region
        with open(target_file) as f:
            file_lines = f.readlines()
            assert len(file_lines) > 0
            for file_line in file_lines:
                cols = file_line.strip().split("\t")
                file_chrom = cols[0]
                file_start = int(cols[1])
                file_end = int(cols[2])

                assert file_chrom == chrom
                assert file_start >= int(start)
                assert file_end <= int(end) + 1  # end is exclusive


def test_multiple_outputs_with_both_placeholders(
    make_test_bed_file, make_test_bam_file_paired, tmp_path
):
    """Test that both {chrom} and {target} placeholders can be used together."""
    wps = WPS(
        bed_file=str(make_test_bed_file),
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
    )

    output_template = str(tmp_path / "wps_{chrom}_{target}.tsv")
    wps.run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=output_template,
    )

    # Check that files were created
    created_files = list(tmp_path.glob("wps_*.tsv"))
    assert len(created_files) > 0, "No output files were created"

    # Verify file naming pattern
    for output_file in created_files:
        filename = output_file.name
        assert filename.startswith("wps_")
        assert filename.endswith(".tsv")
        # Should have format: wps_{chrom}_{chrom}_{start}_{end}.tsv
        parts = filename.replace("wps_", "").replace(".tsv", "").split("_")
        assert len(parts) >= 4, f"Unexpected filename format: {filename}"


def test_multiple_outputs_handles_chunking(
    make_test_bed_file, make_test_bam_file_paired, tmp_path
):
    """Test that multiple chunks of the same region are written to the same file."""
    # Use a very small chunk size to force multiple chunks per region
    wps = WPS(
        bed_file=str(make_test_bed_file),
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
        chunk_size=100,  # Very small chunk size
    )

    output_template = str(tmp_path / "wps_{chrom}.tsv")
    wps.run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=output_template,
    )

    # Verify files were created and have continuous data
    for chrom in ["1", "2", "3", "4", "5", "X"]:
        output_file = tmp_path / f"wps_{chrom}.tsv"
        if output_file.exists():
            with open(output_file) as f:
                lines = f.readlines()
                assert len(lines) > 0

                # Check that positions are continuous (accounting for potential gaps from bed regions)
                positions = [int(line.split("\t")[1]) for line in lines]
                # Positions should be sorted
                assert positions == sorted(
                    positions
                ), f"Positions are not sorted in {output_file}"


def test_single_output_still_works(
    make_test_bed_file, make_test_bam_file_paired, tmp_path
):
    """Test that normal single output mode still works when no placeholders are used."""
    wps = WPS(
        bed_file=str(make_test_bed_file),
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
    )

    output_file = str(tmp_path / "wps_single.tsv")
    wps.run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=output_file,
    )

    # Check that only one file was created
    output_files = list(tmp_path.glob("wps_*.tsv"))
    assert len(output_files) == 1, f"Expected 1 output file, found {len(output_files)}"

    # Verify the file has data
    with open(output_file) as f:
        lines = f.readlines()
        assert len(lines) > 0


def test_multiple_outputs_with_gzip(
    make_test_bed_file, make_test_bam_file_paired, tmp_path
):
    """Test that multiple outputs work with gzipped files."""
    wps = WPS(
        bed_file=str(make_test_bed_file),
        protection_size=120,
        valid_chroms=set(["1", "2", "X", "3", "4", "5"]),
    )

    output_template = str(tmp_path / "wps_{chrom}.tsv.gz")
    wps.run(
        bamfile=str(make_test_bam_file_paired),
        out_filepath=output_template,
    )

    # Check that gzipped files were created
    gz_files = list(tmp_path.glob("wps_*.tsv.gz"))
    assert len(gz_files) > 0, "No gzipped output files were created"

    # Verify we can read the gzipped files
    import gzip

    for gz_file in gz_files:
        with gzip.open(gz_file, "rt") as f:
            lines = f.readlines()
            assert len(lines) > 0
