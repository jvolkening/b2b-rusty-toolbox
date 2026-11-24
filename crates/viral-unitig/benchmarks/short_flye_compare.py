#!/usr/bin/env python3
"""Compare viral-unitig settings with the short-ONT Flye branch."""

import argparse
import csv
import random
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from coverage_bias import make_profile_reads, regional_recall
from matrix import DNA, fasta, make_reads, score


UNITIG_MODES = [
    ("unitig-k21-count3", ["--kmer", "21", "--min-count", "3"]),
    ("unitig-k21-count2", ["--kmer", "21", "--min-count", "2"]),
    ("unitig-k15-count3", ["--kmer", "15", "--min-count", "3"]),
    ("unitig-k15-count2", ["--kmer", "15", "--min-count", "2"]),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("viral_unitig", type=Path)
    parser.add_argument("flye", type=Path)
    parser.add_argument("--sizes", type=int, nargs="+", default=[600, 3000, 10000])
    parser.add_argument("--depth", type=int, default=10)
    parser.add_argument("--profile", choices=["uniform", "biased"], default="uniform")
    parser.add_argument("--floor", type=int, default=20)
    parser.add_argument("--ratio", type=int, default=10)
    parser.add_argument("--mean", type=int, default=500)
    parser.add_argument("--error", type=float, default=0.05)
    parser.add_argument("--flye-iterations", type=int, default=2)
    parser.add_argument("--skip-extend-ends", action="store_true")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--timeout", type=float, default=900)
    parser.add_argument("--seed", type=int, default=20261124)
    parser.add_argument("--seeds", type=int, nargs="+")
    args = parser.parse_args()
    if (
        min(args.sizes) < 100 or args.depth <= 0 or args.mean < 100
        or args.floor <= 0 or args.ratio < 1 or args.flye_iterations <= 0
        or not 0 <= args.error < 1
    ):
        parser.error("invalid size, depth, mean, or error")

    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow([
        "seed", "size", "profile", "regional_mean_coverage", "reads", "assembler",
        "mode", "status", "seconds", "contigs", "total_contig_bp", "longest_contig_bp",
        "reference_31mer_recall", "contig_31mer_precision", "regional_31mer_recall",
    ])
    with tempfile.TemporaryDirectory(prefix="viral-flye-compare-") as directory:
        directory = Path(directory)
        reads_path = directory / "reads.fa"
        for seed in args.seeds or [args.seed]:
            for size in args.sizes:
                reference = "".join(random.Random(seed + size).choices(DNA, k=size))
                rng = random.Random(seed + size * 10000 + args.depth)
                if args.profile == "uniform":
                    reads, _ = make_reads(
                        reads_path, reference, args.depth, args.mean, args.error, rng,
                    )
                    regional_depth = ""
                else:
                    reads, _, regional = make_profile_reads(
                        reads_path, reference, args.mean, args.floor,
                        args.floor * args.ratio, args.error, rng,
                    )
                    regional_depth = ";".join(f"{value:.1f}" for value in regional)
                for mode, options in UNITIG_MODES:
                    output_path = directory / "unitigs.fa"
                    start = time.monotonic()
                    try:
                        with output_path.open("w") as out:
                            result = subprocess.run(
                                [str(args.viral_unitig.resolve()), "--threads", str(args.threads),
                                 *options, str(reads_path)],
                                stdout=out, stderr=subprocess.PIPE, check=False,
                                timeout=args.timeout,
                            )
                        elapsed = time.monotonic() - start
                        status = "ok" if result.returncode == 0 else f"exit_{result.returncode}"
                        sequences = list(fasta(output_path)) if status == "ok" else []
                    except subprocess.TimeoutExpired:
                        elapsed, status, sequences = time.monotonic() - start, "timeout", []
                    write_result(writer, seed, reference, reads, size, args.profile, regional_depth,
                                 "viral-unitig", mode, status, elapsed, sequences)

                extension_modes = (False,) if args.skip_extend_ends else (False, True)
                for extend in extension_modes:
                    mode = "short-ont" + ("-extend-ends" if extend else "")
                    out_dir = directory / f"flye-{seed}-{size}-{'extend' if extend else 'base'}"
                    command = [
                        str(args.flye.resolve()), "--nano-raw", str(reads_path),
                        "--short-reads", "--threads", str(args.threads),
                        "--iterations", str(args.flye_iterations), "--out-dir", str(out_dir),
                    ]
                    if extend:
                        command.append("--extend-ends")
                    start = time.monotonic()
                    try:
                        result = subprocess.run(
                            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, check=False, timeout=args.timeout,
                        )
                        elapsed = time.monotonic() - start
                        assembly = out_dir / "assembly.fasta"
                        status = "ok" if result.returncode == 0 and assembly.exists() else f"exit_{result.returncode}"
                        sequences = list(fasta(assembly)) if status == "ok" else []
                        if status != "ok":
                            print(
                                f"Flye {size} {mode}: {result.stderr[-1000:]}",
                                file=sys.stderr, flush=True,
                            )
                    except subprocess.TimeoutExpired:
                        elapsed, status, sequences = time.monotonic() - start, "timeout", []
                    write_result(writer, seed, reference, reads, size, args.profile, regional_depth,
                                 "Flye", mode, status, elapsed, sequences)
                sys.stdout.flush()


def write_result(
    writer, seed, reference, reads, size, profile, regional_depth, assembler,
    mode, status, elapsed, sequences,
):
    if status == "ok":
        recall, precision, _ = score(reference, sequences)
        regions = regional_recall(reference, sequences, 2 if size < 2000 else 4)
        metrics = [
            len(sequences), sum(map(len, sequences)), max(map(len, sequences), default=0),
            f"{recall:.4f}", f"{precision:.4f}",
            ";".join(f"{value:.4f}" for value in regions),
        ]
    else:
        metrics = [""] * 6
    writer.writerow([
        seed, size, profile, regional_depth, reads, assembler, mode, status,
        f"{elapsed:.2f}", *metrics,
    ])


if __name__ == "__main__":
    main()
