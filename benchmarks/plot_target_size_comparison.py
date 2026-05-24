#!/usr/bin/env python
"""Benchmark WPS runtime while varying BED target size for a fixed BAM."""

import argparse
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam

from plot_input_bam_size_comparison import (
    ROOT,
    WPS,
    create_read,
    old_implementation,
    output_lines,
    summarize,
    timed,
)


def make_bam(workdir, instance, bam_size):
    start, read_len = 1_000, 75
    bam = workdir / f"target_size_reads_{instance}.bam"
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{"LN": start + bam_size + 1_000, "SN": "chr1"}],
    }
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for i in range(max(25, bam_size // 4)):
            insert = bam_size + 220 + (i + instance * 13) % 80
            pos1 = max(1, start - 120 - (i * 7 + instance * 11) % 80)
            pos2 = pos1 + insert - read_len
            name = f"pair_target_size_{instance}_{i}"
            out.write(create_read(name, 0, pos1, read_len, 99, pos2, insert))
            out.write(create_read(name, 0, pos2, read_len, 147, pos1, -insert))
    sorted_bam = workdir / f"target_size_reads_{instance}.sorted.bam"
    pysam.sort("-o", str(sorted_bam), str(bam))
    pysam.index(str(sorted_bam))
    return sorted_bam


def make_bed(workdir, instance, target_size):
    start = 1_000
    bed = workdir / f"target_size_regions_{instance}_{target_size}.bed"
    bed.write_text(f"chr1\t{start}\t{start + target_size - 1}\n")
    return bed


def compare_outputs(new_path, old_path, target_size, instance):
    new_output = output_lines(new_path)
    old_output = output_lines(old_path)
    if new_output == old_output:
        return
    if len(new_output) != len(old_output):
        raise AssertionError(
            f"Output length mismatch for size={target_size}, instance={instance}: "
            f"new={len(new_output)}, old={len(old_output)}"
        )
    mismatch = next(
        i for i, (new, old) in enumerate(zip(new_output, old_output)) if new != old
    )
    raise AssertionError(
        f"Output mismatch for size={target_size}, instance={instance}, "
        f"line={mismatch + 1}:\nnew: {new_output[mismatch]}\n"
        f"old: {old_output[mismatch]}"
    )


def benchmark_instance(workdir, bam, target_size, instance, repeats):
    bed = make_bed(workdir, instance, target_size)
    new_path = workdir / f"new_target_size_{target_size}_{instance}.tsv"
    old_path = workdir / f"old_target_size_{target_size}_{instance}.tsv"
    new = WPS(bed_file=str(bed), protection_size=120, valid_chroms={"1"})

    def run_new():
        new.run(str(bam), str(new_path), verbose_output=True)

    def run_old():
        old_implementation(str(bed), str(bam), str(old_path), 120, {"1"}, None, None)

    new_seconds = timed(run_new, repeats)
    old_seconds = timed(run_old, repeats)
    compare_outputs(new_path, old_path, target_size, instance)
    return new_seconds, old_seconds


def benchmark_target_size(workdir, target_size, instances, repeats, bam_size):
    timings = [
        benchmark_instance(
            workdir,
            make_bam(workdir, instance, bam_size),
            target_size,
            instance,
            repeats,
        )
        for instance in range(instances)
    ]
    new_times = [new for new, _ in timings]
    old_times = [old for _, old in timings]
    return (*summarize(new_times), *summarize(old_times))


def plot(results, output, bam_size):
    target_sizes = [size for size, *_ in results]
    new_means = [new_mean for _, new_mean, _, _, _ in results]
    new_errors = [new_error for _, _, new_error, _, _ in results]
    old_means = [old_mean for _, _, _, old_mean, _ in results]
    old_errors = [old_error for _, _, _, _, old_error in results]
    plt.figure(figsize=(7.2, 4.6), dpi=160)
    plt.errorbar(
        target_sizes,
        old_means,
        yerr=old_errors,
        fmt="o-",
        capsize=4,
        label="Old implementation (Kircher Lab)",
        color="#9b3d2f",
    )
    plt.errorbar(
        target_sizes,
        new_means,
        yerr=new_errors,
        fmt="o-",
        capsize=4,
        label="optwps",
        color="#1f6f8b",
    )
    plt.xlabel("BED target size (bp)")
    plt.ylabel("Runtime (seconds)")
    plt.title(f"WPS runtime by target size in a fixed {bam_size:,} bp BAM")
    plt.grid(True, axis="y", alpha=0.25)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(output)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", default=ROOT / "benchmarks" / "target_size_comparison.svg"
    )
    parser.add_argument("--bam-size", type=int, default=10000)
    parser.add_argument("--instances", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument(
        "--target-sizes",
        type=int,
        nargs="+",
        default=[250, 500, 1000, 2000, 5000, 10000],
    )
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        results = [
            (
                size,
                *benchmark_target_size(
                    workdir, size, args.instances, args.repeats, args.bam_size
                ),
            )
            for size in args.target_sizes
        ]
    plot(results, args.output, args.bam_size)
    for size, new_mean, new_error, old_mean, old_error in results:
        print(
            f"{size}\t"
            f"new={new_mean:.4f}+/-{new_error:.4f}s\t"
            f"old={old_mean:.4f}+/-{old_error:.4f}s\t"
            f"speedup={old_mean / new_mean:.1f}x"
        )


if __name__ == "__main__":
    main()
