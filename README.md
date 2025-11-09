# optwps

A high-performance Python package for computing Window Protection Score (WPS) from BAM files, designed for cell-free DNA (cfDNA) analysis. It was built as a direct alternative of a script provided by the [Kircher Lab](https://github.com/kircherlab/cfDNA.git), and has been tested to replicate the exact numbers.

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

The Windowed Protection Score [![DOI](https://img.shields.io/badge/DOI-110.1016%2Fj.cell.2015.11.050-blue?style=flat-square)](https://doi.org/10.1016/j.cell.2015.11.050) algorithm has the following steps:

1. **Fragment Collection**: For each genomic position, collect all DNA fragments (paired-end reads or single reads) in the region

2. **Protection Window**: Define a protection window of size `protection_size` (default 120bp, or ±60bp from the center)

3. **Score Calculation**:
   - **Outside Score**: Count fragments that completely span the protection window
   - **Inside Score**: Count fragment endpoints that fall within the protection window (exclusive boundaries)
   - **WPS**: Subtract inside score from outside score: `WPS = outside - inside`

4. **Interpretation**: Positive WPS values indicate protected regions (likely nucleosome-bound), while negative values suggest accessible regions


## Examples

### Example 1: Basic WPS Calculation

```bash
python bin/make_wps.py -i sample.bam -o sample_wps.tsv
```

### Example 2: Providing a regions bed file, limiting the range of the size of the inserts considered, and printing to the terminal

```bash
python bin/make_wps.py \
    -i sample.bam \
    -r regions.tsv \
    --min_insert_size 120 \
    --max_insert_size 180
```

### Example 3: Specific Regions with Downsampling

```bash
python bin/make_wps.py \
    -i high_coverage.bam \
    -o regions_wps.tsv \
    -r regions_of_interest.bed \
    --downsample 0.3
```