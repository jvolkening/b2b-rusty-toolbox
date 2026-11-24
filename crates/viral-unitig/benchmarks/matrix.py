#!/usr/bin/env python3
"""Reproducible, bounded single-genome benchmark; see README for limitations."""

import argparse
import csv
import math
import random
import subprocess
import sys
import tempfile
import time
from pathlib import Path

DNA = "ACGT"
COMP = str.maketrans("ACGT", "TGCA")
MODES = {
    "aggressive": ["--kmer", "21", "--min-count", "10", "--min-fraction", "0.4", "--normalize-depth", "0"],
    "sensitive": ["--kmer", "21", "--min-count", "3", "--min-fraction", "0.4", "--normalize-depth", "0"],
    "normalized": ["--kmer", "21", "--min-count", "3", "--min-fraction", "0.4", "--normalize-depth", "50"],
    "normalized-25": ["--kmer", "21", "--min-count", "3", "--min-fraction", "0.4", "--normalize-depth", "25"],
    "count2-normalized": ["--kmer", "21", "--min-count", "2", "--min-fraction", "0.4", "--normalize-depth", "25"],
    "normalized-15": ["--kmer", "21", "--min-count", "3", "--min-fraction", "0.4", "--normalize-depth", "15"],
    "short-k": ["--kmer", "15", "--min-count", "2", "--min-fraction", "0.4", "--normalize-depth", "50"],
    "fast": ["--kmer", "21", "--min-count", "3", "--min-fraction", "0.4", "--normalize-depth", "50", "--no-clean-errors"],
}


def reverse_complement(sequence):
    return sequence.translate(COMP)[::-1]


def make_reads(path, reference, depth, mean, error, rng):
    target = len(reference) * depth
    bases = reads = 0
    log_survival = math.log1p(-error) if error else None
    with path.open("w") as out:
        while bases < target:
            length = min(len(reference), max(100, round(rng.gauss(mean, mean / 4))))
            start = rng.randrange(-length + 1, len(reference))
            fragment = reference[max(start, 0):min(start + length, len(reference))]
            if len(fragment) < 50:
                continue
            bases += len(fragment)
            if rng.randrange(2):
                fragment = reverse_complement(fragment)
            seq = []
            pos = 0
            while log_survival is not None:
                gap = int(math.log1p(-rng.random()) / log_survival)
                event = pos + gap
                if event >= len(fragment):
                    break
                seq.append(fragment[pos:event])
                kind = rng.randrange(4)
                if kind < 2:
                    old = fragment[event]
                    seq.append(rng.choice(DNA.replace(old, "")))
                elif kind == 2:
                    seq.append(rng.choice(DNA) + fragment[event])
                pos = event + 1
            seq.append(fragment[pos:])
            out.write(f">r{reads}\n{''.join(seq)}\n")
            reads += 1
    return reads, bases


def fasta(path):
    seq = []
    with path.open() as inp:
        for line in inp:
            if line.startswith(">"):
                if seq:
                    yield "".join(seq)
                seq = []
            else:
                seq.append(line.strip())
    if seq:
        yield "".join(seq)


def score(reference, sequences, k=31):
    index = {}
    for i in range(len(reference) - k + 1):
        word = reference[i:i + k]
        index.setdefault(min(word, reverse_complement(word)), []).append(i)
    covered = bytearray(len(reference) - k + 1)
    matching = total = 0
    for seq in sequences:
        for i in range(len(seq) - k + 1):
            word = seq[i:i + k]
            positions = index.get(min(word, reverse_complement(word)), ())
            total += 1
            if positions:
                matching += 1
                for position in positions:
                    covered[position] = 1
    unique = sum(covered)
    return (
        unique / len(covered),
        matching / total if total else 0,
        matching / unique if unique else 0,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("assembler", type=Path)
    parser.add_argument("--sizes", type=int, nargs="+", default=[800, 3000, 10000, 30000, 100000])
    parser.add_argument("--depths", type=int, nargs="+", default=[5, 10, 100, 1000, 5000])
    parser.add_argument("--modes", choices=MODES, nargs="+", default=list(MODES))
    parser.add_argument("--mean", type=int, default=500)
    parser.add_argument("--error", type=float, default=0.05)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--seed", type=int, default=20261124)
    args = parser.parse_args()
    if min(args.sizes) < 100 or min(args.depths) <= 0 or args.mean < 100 or not 0 <= args.error < 1 or args.timeout <= 0 or args.threads < 1:
        parser.error("invalid size, depth, mean, error, timeout, or thread count")
    binary = str(args.assembler.resolve())
    output = csv.writer(sys.stdout, lineterminator="\n")
    output.writerow(["size", "depth", "reads", "bases", "mode", "status", "seconds",
                     "contigs", "longest", "reference_31mer_recall", "contig_31mer_precision",
                     "matching_31mer_redundancy"])
    with tempfile.TemporaryDirectory(prefix="viral-matrix-") as directory:
        reads_path = Path(directory) / "reads.fa"
        contigs_path = Path(directory) / "contigs.fa"
        for size in args.sizes:
            reference = "".join(random.Random(args.seed + size).choices(DNA, k=size))
            for depth in args.depths:
                rng = random.Random(args.seed + size * 10000 + depth)
                reads, bases = make_reads(reads_path, reference, depth, args.mean, args.error, rng)
                for mode in args.modes:
                    command = [binary, "--threads", str(args.threads), *MODES[mode], str(reads_path)]
                    start = time.monotonic()
                    try:
                        with contigs_path.open("w") as out:
                            proc = subprocess.run(command, stdout=out, stderr=subprocess.PIPE,
                                                  timeout=args.timeout, check=False)
                        elapsed = time.monotonic() - start
                        status = "ok" if proc.returncode == 0 else f"exit_{proc.returncode}"
                        if proc.returncode:
                            print(f"{size} {depth} {mode}: {proc.stderr.decode(errors='replace')[-1000:]}",
                                  file=sys.stderr, flush=True)
                    except subprocess.TimeoutExpired:
                        elapsed = time.monotonic() - start
                        status = "timeout"
                    if status == "ok":
                        seqs = list(fasta(contigs_path))
                        recall, precision, redundancy = score(reference, seqs)
                        metrics = [len(seqs), max(map(len, seqs), default=0),
                                   f"{recall:.4f}", f"{precision:.4f}", f"{redundancy:.4f}"]
                    else:
                        metrics = [""] * 5
                    output.writerow([size, depth, reads, bases, mode, status, f"{elapsed:.2f}", *metrics])
                    sys.stdout.flush()


if __name__ == "__main__":
    main()
