#!/usr/bin/env python3
"""Small synthetic mixture benchmark; not a model of real nanopore data."""

import random
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: python3 mixed.py /path/to/viral-unitig")
    rng = random.Random(9204)
    complement = str.maketrans("ACGT", "TGCA")
    first = "".join(rng.choices("ACGT", k=1000))
    unrelated = "".join(rng.choices("ACGT", k=1000))
    sites = [150, 300, 450, 600, 750, 850]
    related = list(first)
    for pos in sites:
        related[pos] = next(base for base in "ACGT" if base != related[pos])
    related = "".join(related)

    def reads(reference, depth, label):
        result = []
        for i in range(round(depth * len(reference) / 200)):
            start = rng.randrange(len(reference) - 200 + 1)
            sequence = []
            for base in reference[start : start + 200]:
                event = rng.random()
                if event < 0.005:
                    continue
                if event < 0.01:
                    sequence.append(rng.choice("ACGT"))
                if rng.random() < 0.01:
                    base = rng.choice([other for other in "ACGT" if other != base])
                sequence.append(base)
            result.extend((">" + label + str(i), "".join(sequence)))
        return result

    print("mixture,primary_unitigs,uncertain_unitigs,longest_exact_a,longest_exact_b,alternate_sites")
    for name, second, depth_a, depth_b in [
        ("unrelated", unrelated, 1000, 50),
        ("related", related, 600, 300),
    ]:
        input_data = ("\n".join(reads(first, depth_a, "a") + reads(second, depth_b, "b")) + "\n").encode()
        with tempfile.TemporaryDirectory() as directory:
            uncertain = Path(directory) / "uncertain.fa"
            result = subprocess.run(
                [sys.argv[1], "--uncertain", str(uncertain)], input=input_data,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True,
            )
            sequences = [
                "".join(record.splitlines()[1:])
                for record in result.stdout.decode("ascii").split(">")
                if record
            ]
            oriented = sequences + [seq.translate(complement)[::-1] for seq in sequences]
            longest = lambda ref: max((len(seq) for seq in oriented if seq in ref), default=0)
            alternate_sites = sum(
                any(second[pos - 10 : pos + 11] in seq for seq in oriented)
                for pos in sites
            )
            n_uncertain = sum(line.startswith(">") for line in uncertain.read_text().splitlines())
            print(f"{name},{len(sequences)},{n_uncertain},{longest(first)},{longest(second)},{alternate_sites}")


if __name__ == "__main__":
    main()
