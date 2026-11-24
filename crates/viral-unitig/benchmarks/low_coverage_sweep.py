#!/usr/bin/env python3
"""Parameter sweep for low-depth viral-unitig sensitivity and contiguity."""

import argparse
import csv
import itertools
import random
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from matrix import DNA, fasta, make_reads, score


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("assembler", type=Path)
    parser.add_argument("--sizes", type=int, nargs="+", default=[600, 3000, 10000])
    parser.add_argument("--depths", type=int, nargs="+", default=[5, 10])
    parser.add_argument("--mean", type=int, default=500)
    parser.add_argument("--error", type=float, default=0.05)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--timeout", type=float, default=90)
    parser.add_argument("--seed", type=int, default=20261124)
    parser.add_argument("--bubble-ablation", action="store_true")
    parser.add_argument("--support-ablation", action="store_true")
    args = parser.parse_args()
    if min(args.sizes) < 100 or min(args.depths) <= 0 or args.mean < 100 or not 0 <= args.error < 1:
        parser.error("invalid size, depth, mean, or error")

    if args.bubble_ablation and args.support_ablation:
        parser.error("choose only one ablation mode")
    if args.bubble_ablation:
        configurations = [
            (k, count, fraction, True, bubbles)
            for k, count, fraction, bubbles in itertools.product(
                (15, 21), (2, 3), (0.1, 0.4), (True, False),
            )
        ]
    elif args.support_ablation:
        configurations = []
        for k, count, fraction in itertools.product(
            (15, 21), (2, 3), (0.0, 0.02, 0.05, 0.1, 0.2, 0.4),
        ):
            for clean_errors in ((True,) if fraction == 0 else (True, False)):
                configurations.append((k, count, fraction, clean_errors, False))
    else:
        configurations = [
            (*config, False)
            for config in itertools.product(
                (15, 17, 19, 21),
                (1, 2, 3),
                (0.1, 0.2, 0.4),
                (True, False),
            )
        ]
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow([
        "size", "depth", "reads", "k", "min_count", "min_fraction", "clean_errors", "clean_bubbles",
        "status", "seconds", "contigs", "total_contig_bp", "longest_contig_bp",
        "reference_31mer_recall", "contig_31mer_precision",
    ])
    binary = str(args.assembler.resolve())
    with tempfile.TemporaryDirectory(prefix="viral-lowcov-") as directory:
        directory = Path(directory)
        reads_path = directory / "reads.fa"
        output_path = directory / "contigs.fa"
        for size in args.sizes:
            reference = "".join(random.Random(args.seed + size).choices(DNA, k=size))
            for depth in args.depths:
                reads, _ = make_reads(
                    reads_path, reference, depth, args.mean, args.error,
                    random.Random(args.seed + size * 10000 + depth),
                )
                for k, min_count, min_fraction, clean_errors, clean_bubbles in configurations:
                    options = [
                        "--threads", str(args.threads),
                        "--kmer", str(k),
                        "--min-count", str(min_count),
                        "--min-fraction", str(min_fraction),
                    ]
                    if not clean_errors:
                        options.append("--no-clean-errors")
                    if clean_bubbles:
                        options.append("--clean-bubbles")
                    start = time.monotonic()
                    try:
                        with output_path.open("w") as out:
                            result = subprocess.run(
                                [binary, *options, str(reads_path)],
                                stdout=out,
                                stderr=subprocess.PIPE,
                                check=False,
                                timeout=args.timeout,
                            )
                        elapsed = time.monotonic() - start
                        status = "ok" if result.returncode == 0 else f"exit_{result.returncode}"
                        if result.returncode:
                            print(
                                f"{size}bp {depth}x k={k} min-count={min_count} "
                                f"min-fraction={min_fraction} clean={clean_errors}: "
                                f"{result.stderr.decode(errors='replace')[-500:]}",
                                file=sys.stderr,
                                flush=True,
                            )
                    except subprocess.TimeoutExpired:
                        elapsed = time.monotonic() - start
                        status = "timeout"
                    if status == "ok":
                        contigs = list(fasta(output_path))
                        recall, precision, _ = score(reference, contigs)
                        metrics = [
                            len(contigs),
                            sum(map(len, contigs)),
                            max(map(len, contigs), default=0),
                            f"{recall:.4f}",
                            f"{precision:.4f}",
                        ]
                    else:
                        metrics = [""] * 5
                    writer.writerow([
                        size, depth, reads, k, min_count, min_fraction, clean_errors, clean_bubbles,
                        status, f"{elapsed:.2f}", *metrics,
                    ])
                    sys.stdout.flush()


if __name__ == "__main__":
    main()
