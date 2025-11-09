# optwps

[![PyPI version](https://badge.fury.io/py/optwps.svg)](https://badge.fury.io/py/optwps)
[![codecov](https://codecov.io/gh/VasLem/optwps/branch/master/graph/badge.svg)](https://codecov.io/gh/VasLem/optwps)
[![DOI](https://zenodo.org/badge/1092793606.svg)](https://doi.org/10.5281/zenodo.17566994)

A high-performance Python package for computing Window Protection Score (WPS) from BAM files, designed for cell-free DNA (cfDNA) analysis. It was built as a direct alternative of a script provided by the [Kircher Lab](https://github.com/kircherlab/cfDNA.git), and has been tested to replicate the exact numbers.

## Overview

`optwps` is a fast and efficient tool for calculating Window Protection Scores from aligned sequencing reads. WPS is a metric used in cell-free DNA analysis to identify nucleosome positioning and protected regions by analyzing fragment coverage patterns.

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
optwps -i input.bam -o output.tsv
```

With custom parameters:

```bash
optwps \
    -i input.bam \
    -o output.tsv \
    -w 120 \
    --min_insert_size 120 \
    --max_insert_size 180 \
    --downsample 0.5
```

### Command Line Arguments

- `-i, --input`: Input BAM file (required)
- `-o, --outfile`: Output file path for WPS results. If not provided, results will be printed to stdout (optional)
- `-r, --regions`: BED file with regions of interest (default: whole genome, optional)
- `-w, --protection`: Base pair protection window (default: 120)
- `--min-insert-size`: Minimum read length threshold to consider (optional)
- `--max-insert-size`: Maximum read length threshold to consider (optional)
- `--downsample`: Ratio to downsample reads (default OFF, optional)
- `--chunk-size`: Chunk size for processing in pieces (default: 1e8)
- `--valid-chroms`: Comma-separated list of valid chromosomes to include (e.g., '1,2,3,X,Y') or 'canonical' for chromosomes 1-22, X, Y (optional)
- `--verbose-output`: If provided, output will include separate counts for 'outside' and 'inside' along with WPS

### Python API

```python
from optwps import WPS

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
4. **outside**: Count of fragments fully spanning the protection window (if `--verbose-output`)
5. **inside**: Count of fragment endpoints falling inside the protection window (if `--verbose-output`)
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
optwps -i sample.bam -o sample_wps.tsv
```

### Example 2: Providing a regions bed file, limiting the range of the size of the inserts considered, and printing to the terminal

```bash
optwps \
    -i sample.bam \
    -r regions.tsv \
    --min_insert_size 120 \
    --max_insert_size 180
```

### Example 3: Specific Regions with Downsampling

```bash
optwps \
    -i high_coverage.bam \
    -o regions_wps.tsv \
    -r regions_of_interest.bed \
    --downsample 0.3
```