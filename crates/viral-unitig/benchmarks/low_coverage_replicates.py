#!/usr/bin/env python3
"""Repeat promising low-coverage settings across independent random seeds."""

import argparse
import csv
import random
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from matrix import DNA, fasta, make_reads, score


CONFIGS = [
    ("k21-count3", ["--kmer", "21", "--min-count", "3"]),
    ("k21-count2", ["--kmer", "21", "--min-count", "2"]),
    ("k15-count3", ["--kmer", "15", "--min-count", "3"]),
    ("k15-count2", ["--kmer", "15", "--min-count", "2"]),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("assembler", type=Path)
    parser.add_argument("--sizes", type=int, nargs="+", default=[600, 3000, 10000])
    parser.add_argument("--depths", type=int, nargs="+", default=[5, 10])
    parser.add_argument("--seeds", type=int, nargs="+", default=[20261124, 20261125, 20261126, 20261127, 20261128])
    parser.add_argument("--mean", type=int, default=500)
    parser.add_argument("--error", type=float, default=0.05)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--timeout", type=float, default=90)
    args = parser.parse_args()
    if min(args.sizes) < 100 or min(args.depths) <= 0 or args.mean < 100 or not 0 <= args.error < 1:
        parser.error("invalid size, depth, mean, or error")

    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow([
        "seed", "size", "depth", "reads", "mode", "status", "seconds", "contigs",
        "total_contig_bp", "longest_contig_bp", "reference_31mer_recall",
        "contig_31mer_precision",
    ])
    binary = str(args.assembler.resolve())
    with tempfile.TemporaryDirectory(prefix="viral-replicates-") as directory:
        directory = Path(directory)
        reads_path, output_path = directory / "reads.fa", directory / "contigs.fa"
        for seed in args.seeds:
            for size in args.sizes:
                reference = "".join(random.Random(seed + size).choices(DNA, k=size))
                for depth in args.depths:
                    reads, _ = make_reads(
                        reads_path, reference, depth, args.mean, args.error,
                        random.Random(seed + size * 10000 + depth),
                    )
                    for mode, options in CONFIGS:
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
                            metrics = [
                                len(contigs), sum(map(len, contigs)),
                                max(map(len, contigs), default=0),
                                f"{recall:.4f}", f"{precision:.4f}",
                            ]
                        else:
                            metrics = [""] * 5
                        writer.writerow([
                            seed, size, depth, reads, mode, status,
                            f"{elapsed:.2f}", *metrics,
                        ])
                        sys.stdout.flush()


if __name__ == "__main__":
    main()
