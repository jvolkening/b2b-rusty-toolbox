#!/usr/bin/env python3
"""Small synthetic smoke benchmark; not a model of real nanopore data."""

import random
import subprocess
import sys


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: python3 synthetic.py /path/to/viral-unitig")
    random.seed(4096)
    reference = "".join(random.choices("ACGT", k=3000))
    complement = str.maketrans("ACGT", "TGCA")
    print("substitution_rate,indel_rate,nominal_depth,unitigs,longest_exact")
    for substitutions, indels, depth in [(0.01, 0.01, 100), (0.02, 0.02, 100), (0.02, 0.02, 1000)]:
        reads = []
        for i in range(round(depth * len(reference) / 250)):
            start = random.randrange(len(reference) - 250 + 1)
            sequence = []
            for base in reference[start : start + 250]:
                event = random.random()
                if event < indels / 2:
                    continue
                if event < indels:
                    sequence.append(random.choice("ACGT"))
                if random.random() < substitutions:
                    base = random.choice([other for other in "ACGT" if other != base])
                sequence.append(base)
            reads.extend((">r" + str(i), "".join(sequence)))
        result = subprocess.run(
            [sys.argv[1]], input=("\n".join(reads) + "\n").encode(),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True,
        )
        unitigs = [
            "".join(record.splitlines()[1:])
            for record in result.stdout.decode("ascii").split(">")
            if record
        ]
        exact = [
            len(seq) for seq in unitigs
            if seq in reference or seq.translate(complement)[::-1] in reference
        ]
        print(f"{substitutions},{indels},{depth},{len(unitigs)},{max(exact, default=0)}")


if __name__ == "__main__":
    main()
