# optwps

A high-performance Python package for computing Window Protection Score (WPS) from BAM files, designed for cell-free DNA (cfDNA) analysis.

## Overview

`optwps` is a fast and efficient tool for calculating Window Protection Scores from aligned sequencing reads. WPS is a metric used in cell-free DNA analysis to identify nucleosome positioning and protected regions by analyzing fragment coverage patterns.

## Features

- **Fast Processing**: Optimized numpy-based implementation for efficient WPS calculation
- **Flexible Input**: Supports both paired-end and single-end sequencing data
- **Region-Based Analysis**: Process specific genomic regions via BED files or entire genomes
- **Customizable Parameters**: Adjust protection size, insert size filters, and downsampling
- **Parallel Processing**: Automatic multi-core support for compressed file handling
- **Memory Efficient**: Chunked processing for large genomic regions

## Installation

### From Source

```bash
git clone <repository-url>
cd wps
pip install -r requirements.txt
```

### Dependencies

- Python >= 3.7
- pysam
- numpy
- pgzip
- tqdm
- bx-python

## Usage

### Command Line Interface

Basic usage:

```bash
python bin/make_wps.py -i input.bam -o output.tsv
```

With custom parameters:

```bash
python bin/make_wps.py \
    -i input.bam \
    -o output.tsv \
    -w 120 \
    --min_insert_size 120 \
    --max_insert_size 180 \
    --downsample 0.5
```

### Command Line Arguments

- `-i, --input`: Input BAM file (required)
- `-o, --outfile`: Output file path (required)
- `-w, --protection`: Base pair protection size assumed for elements (default: 120)
- `--min_insert_size`: Minimum insert size threshold to consider (optional)
- `--max_insert_size`: Maximum insert size threshold to consider (optional)
- `--downsample`: Ratio to downsample reads (0.0-1.0, optional)
- `-v, --verbose`: Enable debug output

### Python API

```python
from lib.optwps import WPS

# Initialize WPS calculator
wps_calculator = WPS(
    protection_size=120,
    min_insert_size=120,
    max_insert_size=180,
    valid_chroms=set(map(str, list(range(1, 23)) + ['X', 'Y']))
)

# Run WPS calculation
wps_calculator.run(
    bamfile='input.bam',
    out_filepath='output.tsv',
    downsample_ratio=0.5
)
```

### Using BED Files for Specific Regions

```python
from lib.optwps import WPS

# Process only specific regions
wps_calculator = WPS(
    bed_file='regions.bed',
    protection_size=120
)

wps_calculator.run(
    bamfile='input.bam',
    out_filepath='output.tsv'
)
```

## Output Format

The output is a tab-separated file with the following columns:

1. **chromosome**: Chromosome name
2. **start**: Start position (0-based)
3. **end**: End position (1-based)
4. **outside**: Count of fragments fully spanning the protection window
5. **inside**: Count of fragment endpoints falling inside the protection window
6. **wps**: Window Protection Score (outside - inside)

Example output:

```
1    1000    1001    15    3    12
1    1001    1002    16    2    14
1    1002    1003    14    4    10
```

## Algorithm

The WPS algorithm works as follows:

1. **Fragment Collection**: For each genomic position, collect all DNA fragments (paired-end reads or single reads) in the region

2. **Protection Window**: Define a protection window of size `protection_size` (default 120bp, or ±60bp from the center)

3. **Score Calculation**:
   - **Outside Score**: Count fragments that completely span the protection window
   - **Inside Score**: Count fragment endpoints that fall within the protection window
   - **WPS**: Subtract inside score from outside score: `WPS = outside - inside`

4. **Interpretation**: Positive WPS values indicate protected regions (likely nucleosome-bound), while negative values suggest accessible regions

## Performance Considerations

- **Chunk Size**: Processing is done in chunks (default 1MB) to manage memory usage
- **Downsampling**: Use `--downsample` to reduce processing time for high-coverage samples
- **Parallel Processing**: Compressed output files automatically use multiple cores via `pgzip`
- **Valid Chromosomes**: By default, only autosomes and sex chromosomes are processed

## Quality Control

The tool automatically filters out:
- Duplicate reads (`read.is_duplicate`)
- QC-failed reads (`read.is_qcfail`)
- Unmapped reads (`read.is_unmapped`)
- Soft-clipped reads
- Reads with unmapped mates (paired-end mode)
- Inter-chromosomal pairs (paired-end mode)

## Examples

### Example 1: Basic WPS Calculation

```bash
python bin/make_wps.py -i sample.bam -o sample_wps.tsv
```

### Example 2: Nucleosome-Sized Fragments Only

```bash
python bin/make_wps.py \
    -i sample.bam \
    -o nucleosome_wps.tsv \
    --min_insert_size 120 \
    --max_insert_size 180
```

### Example 3: Specific Regions with Downsampling

```bash
python bin/make_wps.py \
    -i high_coverage.bam \
    -o regions_wps.tsv \
    --bed regions_of_interest.bed \
    --downsample 0.3
```

## Testing

Run the test suite:

```bash
pytest tests/
```

Run with verbose output:

```bash
pytest tests/ -v
```

## Citation

If you use `optwps` in your research, please cite:

[Add citation information here]

## License

[Add license information here]

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Contact

[Add contact information here]

## Acknowledgments

This tool implements the Window Protection Score algorithm as described in the cell-free DNA literature for nucleosome positioning analysis.

## Troubleshooting

### Common Issues

**Issue**: "No module named 'lib.optwps'"
- **Solution**: Make sure you're running from the project root directory or install the package

**Issue**: Slow processing for large BAM files
- **Solution**: Use `--downsample` option or process specific regions with a BED file

**Issue**: High memory usage
- **Solution**: Reduce chunk_size parameter when initializing WPS class

**Issue**: Missing chromosomes in output
- **Solution**: Check that chromosome names match (with/without "chr" prefix) and verify `valid_chroms` parameter

## Version History

### Version 1.0.0
- Initial release
- Support for paired-end and single-end reads
- BED file region support
- Configurable protection size and insert size filters
- Downsampling support
