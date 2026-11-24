#!/usr/bin/env python3
"""Reproducible 10-kb mixed-strain smoke benchmark, not real nanopore data."""

import json
import random
import subprocess
import sys
import tempfile
import time
from pathlib import Path

if len(sys.argv) < 2:
    sys.exit("usage: python3 strains.py /path/to/viral-unitig [assembler options]")

binary = Path(sys.argv[1]).resolve()
options = sys.argv[2:]
if "--uncertain" in options:
    sys.exit("the benchmark manages --uncertain itself")
complement = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence):
    return sequence.translate(complement)[::-1]


def canonical(sequence):
    return min(sequence, reverse_complement(sequence))


def fasta_sequences(path):
    sequence = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if sequence:
                    yield "".join(sequence)
                    sequence.clear()
            else:
                sequence.append(line.strip())
        if sequence:
            yield "".join(sequence)


def run(divergence):
    reference_rng = random.Random(39173)
    major = "".join(reference_rng.choices("ACGT", k=10_000))
    mutation_rng = random.Random(4094)
    positions = mutation_rng.sample(range(10_000), 1000)[: 100 * divergence]
    minor = list(major)
    for pos in positions:
        minor[pos] = mutation_rng.choice([base for base in "ACGT" if base != major[pos]])
    minor = "".join(minor)
    reads_rng = random.Random(1717 + divergence)

    with tempfile.TemporaryDirectory() as directory:
        directory = Path(directory)
        reads, primary, uncertain = (directory / name for name in ("reads.fa", "primary.fa", "uncertain.fa"))
        with reads.open("w") as handle:
            for label, sequence, depth in [("major", major, 1000), ("minor", minor, 200)]:
                for index in range(depth * len(sequence) // 250):
                    start = reads_rng.randrange(len(sequence) - 250 + 1)
                    read = []
                    for base in sequence[start : start + 250]:
                        event = reads_rng.random()
                        if event < 0.01:
                            continue
                        if event < 0.02:
                            read.append(reads_rng.choice("ACGT"))
                        if reads_rng.random() < 0.02:
                            base = reads_rng.choice([other for other in "ACGT" if other != base])
                        read.append(base)
                    sequence_read = "".join(read)
                    if reads_rng.randrange(2):
                        sequence_read = reverse_complement(sequence_read)
                    handle.write(f">{label}_{index}\n{sequence_read}\n")

        start = time.monotonic()
        with primary.open("w") as output:
            result = subprocess.run(
                [str(binary), *options, "--uncertain", str(uncertain), str(reads)],
                stdout=output,
                stderr=subprocess.PIPE,
                text=True,
            )
        elapsed = time.monotonic() - start
        if result.returncode:
            sys.exit(f"assembler failed:\n{result.stderr}")

        contigs = list(fasta_sequences(primary))
        uncertain_count = sum(1 for _ in fasta_sequences(uncertain))
        k = 21
        words_major = [canonical(major[i : i + k]) for i in range(len(major) - k + 1)]
        words_minor = [canonical(minor[i : i + k]) for i in range(len(minor) - k + 1)]
        set_major, set_minor = set(words_major), set(words_minor)
        observed = {
            canonical(contig[i : i + k])
            for contig in contigs
            for i in range(len(contig) - k + 1)
        }

        def summarize(reference, reference_words, other_words):
            informative = [word for word in reference_words if word not in other_words]
            matched = sum(word in observed for word in informative)
            longest = max(
                (len(contig) for contig in contigs if contig in reference or reverse_complement(contig) in reference),
                default=0,
            )
            return {
                "informative_21mer_positions": len(informative),
                "matched_informative_21mer_positions": matched,
                "informative_21mer_recovery_percent": round(100 * matched / len(informative), 2),
                "longest_exact_unitig": longest,
            }

        return {
            "divergence_percent": divergence,
            "options": options,
            "elapsed_seconds": round(elapsed, 2),
            "primary_unitigs": len(contigs),
            "longest_unitigs": sorted((len(contig) for contig in contigs), reverse=True)[:10],
            "uncertain_unitigs": uncertain_count,
            "major": summarize(major, words_major, set_minor),
            "minor": summarize(minor, words_minor, set_major),
            "assembler_stderr": result.stderr.strip(),
        }


for divergence in (10, 5):
    print(json.dumps(run(divergence), sort_keys=True), flush=True)
