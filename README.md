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
pip install optwps
```

### Dependencies

- Python >= 3.7
- pysam
- pandas
- pigz
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
- `-o, --output`: Output file path for WPS results. If not provided, results will be printed to stdout. Supports placeholders `{chrom}` and `{target}` for creating separate files per chromosome or region (optional)
- `-r, --regions`: BED file with regions of interest (default: whole genome, optional)
- `-w, --protection`: Base pair protection window (default: 120)
- `--min-insert-size`: Minimum read length threshold to consider (optional)
- `--max-insert-size`: Maximum read length threshold to consider (optional)
- `--downsample`: Ratio to downsample reads (default OFF, optional)
- `--chunk-size`: Chunk size for processing in pieces (default: 1e8)
- `--valid-chroms`: Comma-separated list of valid chromosomes to include (e.g., '1,2,3,X,Y') or 'canonical' for chromosomes 1-22, X, Y (optional)
- `--compute-coverage`: If provided, output will include base coverage
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

The output is a tab-separated no-header file with the following columns:

    - Chromosome name (without 'chr' prefix)
    - Start position (0-based)
    - End position (start + 1)
    - Base read coverage (if `--compute-coverage`)
    - Count of fragments spanning the protection window (if `--verbose-output`)
    - Count of fragment endpoints in protection window (if `--verbose-output`)
    - Window Protection Score (outside - inside)

Example output:

```
1\t1000\t1001\t12
1\t1001\t1002\t14
1\t1002\t1003\t10
```

With `--compute-coverage`
```
1\t1000\t1001\t20\t12
1\t1001\t1002\t20\t14
1\t1002\t1003\t19\t10
```

With `--verbose-output`:

```
1\t1000\t1001\t15\t3\t12
1\t1001\t1002\t16\t2\t14
1\t1002\t1003\t14\t4\t10
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
    -o wps.tsv \
    --downsample 0.3
```

### Example 4: Creating Separate Output Files per Chromosome

```bash
optwps \
    -i sample.bam \
    -o "wps_{chrom}.tsv"
```

### Example 5: Include coverage

```bash
optwps \
    -i sample.bam \
    --compute_coverage \
    -o "wps.tsv"
```