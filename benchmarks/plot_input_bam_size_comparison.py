#!/usr/bin/env python
"""Benchmark WPS runtime while varying BAM size for a fixed 2,000 bp target."""

import argparse
import sys
import tempfile
from pathlib import Path
from statistics import mean, median, stdev
from time import perf_counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT)]

from optwps import WPS  # noqa: E402
from tests.test_optwps import old_implementation  # noqa: E402


def create_read(name, ref_id, pos, length, flag, mate_pos, isize):
    read = pysam.AlignedSegment()
    read.query_name = name
    read.reference_id = ref_id
    read.reference_start = pos
    read.cigar = ((0, length),)
    read.mapping_quality = 60
    read.query_sequence = "A" * length
    read.query_qualities = pysam.qualitystring_to_array("I" * length)
    read.flag = flag
    read.next_reference_id = ref_id
    read.next_reference_start = mate_pos
    read.template_length = isize
    return read


def make_inputs(workdir, read_pairs, instance, target_size):
    start, end, read_len = 1_000, 1_000 + target_size - 1, 75
    bed = workdir / f"regions_{read_pairs}_{instance}.bed"
    bed.write_text(f"chr1\t{start}\t{end}\n")
    bam = workdir / f"reads_{read_pairs}_{instance}.bam"
    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": end + 1_000, "SN": "chr1"}]}
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for i in range(read_pairs):
            insert = target_size + 220 + (i + instance * 13) % 80
            pos1 = max(1, start - 120 - (i * 7 + instance * 11) % 80)
            pos2 = pos1 + insert - read_len
            name = f"pair_{read_pairs}_{instance}_{i}"
            out.write(create_read(name, 0, pos1, read_len, 99, pos2, insert))
            out.write(create_read(name, 0, pos2, read_len, 147, pos1, -insert))
    sorted_bam = workdir / f"reads_{read_pairs}_{instance}.sorted.bam"
    pysam.sort("-o", str(sorted_bam), str(bam))
    pysam.index(str(sorted_bam))
    return bed, sorted_bam


def timed(func, repeats):
    times = []
    for _ in range(repeats):
        start = perf_counter()
        func()
        times.append(perf_counter() - start)
    return median(times)


def output_lines(path):
    return Path(path).read_text().splitlines()


def benchmark_instance(workdir, read_pairs, instance, repeats, target_size):
    bed, bam = make_inputs(workdir, read_pairs, instance, target_size)
    new = WPS(bed_file=str(bed), protection_size=120, valid_chroms={"1"}, njobs=1)

    def run_new():
        new.run(
            str(bam),
            str(workdir / f"new_{read_pairs}_{instance}.tsv"),
            verbose_output=True,
        )

    def run_old():
        old_implementation(
            str(bed),
            str(bam),
            str(workdir / f"old_{read_pairs}_{instance}.tsv"),
            120,
            {"1"},
            None,
            None,
        )

    new_seconds = timed(run_new, repeats)
    old_seconds = timed(run_old, repeats)
    new_output = output_lines(workdir / f"new_{read_pairs}_{instance}.tsv")
    old_output = output_lines(workdir / f"old_{read_pairs}_{instance}.tsv")
    if new_output != old_output:
        if len(new_output) != len(old_output):
            raise AssertionError(
                f"Output length mismatch for read_pairs={read_pairs}, "
                f"instance={instance}: "
                f"new={len(new_output)}, old={len(old_output)}"
            )
        mismatch = next(
            i for i, (new, old) in enumerate(zip(new_output, old_output)) if new != old
        )
        raise AssertionError(
            f"Output mismatch for read_pairs={read_pairs}, instance={instance}, "
            f"line={mismatch + 1}:\nnew: {new_output[mismatch]}\n"
            f"old: {old_output[mismatch]}"
        )
    return new_seconds, old_seconds


def summarize(values):
    return mean(values), stdev(values) if len(values) > 1 else 0


def benchmark_size(workdir, read_pairs, instances, repeats, target_size):
    timings = [
        benchmark_instance(workdir, read_pairs, instance, repeats, target_size)
        for instance in range(instances)
    ]
    new_times = [new for new, _ in timings]
    old_times = [old for _, old in timings]
    return (*summarize(new_times), *summarize(old_times))


def plot(results, output, target_size):
    read_pair_counts = [read_pairs for read_pairs, *_ in results]
    new_means = [new_mean for _, new_mean, _, _, _ in results]
    new_errors = [new_error for _, _, new_error, _, _ in results]
    old_means = [old_mean for _, _, _, old_mean, _ in results]
    old_errors = [old_error for _, _, _, _, old_error in results]
    plt.figure(figsize=(7.2, 4.6), dpi=160)
    plt.errorbar(
        read_pair_counts,
        old_means,
        yerr=old_errors,
        fmt="o-",
        capsize=4,
        label="Old implementation (Kircher Lab)",
        color="#9b3d2f",
    )
    plt.errorbar(
        read_pair_counts,
        new_means,
        yerr=new_errors,
        fmt="o-",
        capsize=4,
        label="optwps",
        color="#1f6f8b",
    )
    plt.xlabel("Synthetic BAM size (read pairs)")
    plt.ylabel("Runtime (seconds)")
    plt.title(f"WPS runtime by BAM size for a fixed {target_size:,} bp target")
    plt.grid(True, axis="y", alpha=0.25)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(output)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", default=ROOT / "benchmarks" / "input_bam_size_comparison.svg"
    )
    parser.add_argument("--target-size", type=int, default=2000)
    parser.add_argument("--instances", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument(
        "--read-pairs",
        type=int,
        nargs="+",
        default=[25, 50, 100, 250, 500, 1000],
    )
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        results = [
            (
                read_pairs,
                *benchmark_size(
                    workdir,
                    read_pairs,
                    args.instances,
                    args.repeats,
                    args.target_size,
                ),
            )
            for read_pairs in args.read_pairs
        ]
    plot(results, args.output, args.target_size)
    for read_pairs, new_mean, new_error, old_mean, old_error in results:
        print(
            f"read_pairs={read_pairs}\t"
            f"target_size={args.target_size}\t"
            f"new={new_mean:.4f}+/-{new_error:.4f}s\t"
            f"old={old_mean:.4f}+/-{old_error:.4f}s\t"
            f"speedup={old_mean / new_mean:.1f}x"
        )


if __name__ == "__main__":
    main()
