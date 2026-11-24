#!/usr/bin/env python3
"""Compare uniform and biased regional depth with the same low-depth floor."""

import argparse
import csv
import math
import random
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from matrix import DNA, fasta, score

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def regional_recall(reference, contigs, segments):
    k = 31
    index = {}
    for i in range(len(reference) - k + 1):
        word = reference[i:i + k]
        reverse = word.translate(COMPLEMENT)[::-1]
        index.setdefault(min(word, reverse), []).append(i)
    recovered = bytearray(len(reference) - k + 1)
    for contig in contigs:
        for i in range(len(contig) - k + 1):
            word = contig[i:i + k]
            reverse = word.translate(COMPLEMENT)[::-1]
            for position in index.get(min(word, reverse), ()):
                recovered[position] = 1
    bounds = [len(reference) * i // segments for i in range(segments + 1)]
    return [
        sum(recovered[left:min(right, len(recovered))])
        / max(1, min(right, len(recovered)) - left)
        for left, right in zip(bounds, bounds[1:])
    ]


def make_profile_reads(path, reference, mean_length, low_depth, high_depth, error, rng):
    length = len(reference)
    segments = 2 if length < 2000 else 4
    bounds = [length * i // segments for i in range(segments + 1)]
    depths = (
        [low_depth, high_depth]
        if segments == 2
        else [high_depth, low_depth, low_depth, high_depth]
    )
    segment_weights = [
        (bounds[i + 1] - bounds[i]) * depth for i, depth in enumerate(depths)
    ]

    target_bases = round(sum((bounds[i + 1] - bounds[i]) * depth for i, depth in enumerate(depths)))
    covered = [0] * length
    bases = reads = 0
    log_survival = math.log1p(-error) if error else None
    with path.open("w") as out:
        while bases < target_bases:
            read_length = min(length, max(100, round(rng.gauss(mean_length, mean_length / 4))))
            segment = rng.choices(range(segments), weights=segment_weights, k=1)[0]
            center = rng.randrange(bounds[segment], bounds[segment + 1])
            start = center - read_length // 2
            end = start + read_length
            clipped_start, clipped_end = max(0, start), min(length, end)
            fragment = reference[clipped_start:clipped_end]
            for position in range(clipped_start, clipped_end):
                covered[position] += 1
            bases += len(fragment)
            sequence = []
            pos = 0
            if log_survival is not None:
                while True:
                    gap = int(math.log1p(-rng.random()) / log_survival)
                    event = pos + gap
                    if event >= len(fragment):
                        break
                    sequence.append(fragment[pos:event])
                    kind = rng.randrange(4)
                    if kind < 2:
                        original = fragment[event]
                        sequence.append(rng.choice(DNA.replace(original, "")))
                    elif kind == 2:
                        sequence.append(rng.choice(DNA) + fragment[event])
                    pos = event + 1
            sequence.append(fragment[pos:])
            read = "".join(sequence)
            if rng.randrange(2):
                read = read.translate(COMPLEMENT)[::-1]
            out.write(f">r{reads}\n{read}\n")
            reads += 1
    regional = []
    for i in range(segments):
        left, right = bounds[i], bounds[i + 1]
        trim = min(mean_length // 2, (right - left - 1) // 2)
        interior = covered[left + trim:right - trim]
        regional.append(sum(interior) / len(interior))
    return reads, bases, regional


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("assembler", type=Path)
    parser.add_argument("--sizes", type=int, nargs="+", default=[800, 3000, 10000])
    parser.add_argument("--means", type=int, nargs="+", default=[300, 500])
    parser.add_argument("--floor", type=int, default=20)
    parser.add_argument("--ratios", type=int, nargs="+", default=[1, 10, 20])
    parser.add_argument("--error", type=float, default=0.05)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--timeout", type=float, default=90)
    parser.add_argument("--seed", type=int, default=20261124)
    args = parser.parse_args()
    if min(args.sizes) < 100 or min(args.means) < 100 or args.floor <= 0 or min(args.ratios) < 1:
        parser.error("invalid size, read length, floor, or ratio")
    configs = [
        ("k21-count3-f04", [
            "--kmer", "21", "--min-count", "3", "--min-fraction", "0.4",
        ]),
        ("k21-count3-f02", [
            "--kmer", "21", "--min-count", "3", "--min-fraction", "0.2",
        ]),
        ("k21-count3-f01", [
            "--kmer", "21", "--min-count", "3", "--min-fraction", "0.1",
        ]),
        ("k21-count3-f00", [
            "--kmer", "21", "--min-count", "3", "--min-fraction", "0",
        ]),
        ("k21-count3-no-clean", [
            "--kmer", "21", "--min-count", "3", "--min-fraction", "0.4",
            "--no-clean-errors",
        ]),
        ("k21-count2-f04", [
            "--kmer", "21", "--min-count", "2", "--min-fraction", "0.4",
        ]),
        ("k21-count2-f01", [
            "--kmer", "21", "--min-count", "2", "--min-fraction", "0.1",
        ]),
        ("k15-count3-f04", [
            "--kmer", "15", "--min-count", "3", "--min-fraction", "0.4",
        ]),
        ("k15-count2-f04", [
            "--kmer", "15", "--min-count", "2", "--min-fraction", "0.4",
        ]),
        ("k21-count2-normalized", [
            "--kmer", "21", "--min-count", "2", "--min-fraction", "0.4",
            "--normalize-depth", "25",
        ]),
    ]
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow([
        "size", "mean_read_length", "coverage_ratio", "profile", "target_mean_depth",
        "reads", "total_bases",
        "regional_mean_coverage", "mode", "status", "seconds", "contigs",
        "total_contig_bp", "longest_contig_bp", "reference_31mer_recall",
        "contig_31mer_precision", "regional_31mer_recall",
    ])
    binary = str(args.assembler.resolve())
    with tempfile.TemporaryDirectory(prefix="viral-bias-") as directory:
        directory = Path(directory)
        reads_path, output_path = directory / "reads.fa", directory / "contigs.fa"
        for size in args.sizes:
            reference = "".join(random.Random(args.seed + size).choices(DNA, k=size))
            for mean in args.means:
                for ratio in args.ratios:
                    high_depth = args.floor * ratio
                    patterns = [("uniform", args.floor, args.floor)]
                    if ratio > 1:
                        matched_mean = round(args.floor * (ratio + 1) / 2)
                        patterns = [
                            ("uniform-mean", matched_mean, matched_mean),
                            ("biased", args.floor, high_depth),
                        ]
                    for profile, low, high in patterns:
                        reads, bases, regional = make_profile_reads(
                            reads_path, reference, mean, low, high, args.error,
                            random.Random(args.seed + size * 1000 + mean * 10 + ratio),
                        )
                        actual_depth = ";".join(f"{value:.1f}" for value in regional)
                        for mode, options in configs:
                            start = time.monotonic()
                            try:
                                with output_path.open("w") as out:
                                    result = subprocess.run(
                                        [binary, "--threads", str(args.threads), *options, str(reads_path)],
                                        stdout=out, stderr=subprocess.PIPE, check=False,
                                        timeout=args.timeout,
                                    )
                                elapsed = time.monotonic() - start
                                status = "ok" if result.returncode == 0 else f"exit_{result.returncode}"
                            except subprocess.TimeoutExpired:
                                elapsed = time.monotonic() - start
                                status = "timeout"
                            if status == "ok":
                                contigs = list(fasta(output_path))
                                recall, precision, _ = score(reference, contigs)
                                region_recall = regional_recall(
                                    reference, contigs, 2 if size < 2000 else 4,
                                )
                                metrics = [
                                    len(contigs), sum(map(len, contigs)),
                                    max(map(len, contigs), default=0),
                                    f"{recall:.4f}", f"{precision:.4f}",
                                    ";".join(f"{value:.4f}" for value in region_recall),
                                ]
                            else:
                                metrics = [""] * 6
                            writer.writerow([
                                size, mean, ratio, profile, f"{(low + high) / 2:.1f}",
                                reads, bases, actual_depth, mode, status,
                                f"{elapsed:.2f}", *metrics,
                            ])
                            sys.stdout.flush()


if __name__ == "__main__":
    main()
