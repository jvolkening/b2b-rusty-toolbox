import random
import tempfile
import unittest
from pathlib import Path

from matrix import fasta, make_reads, reverse_complement, score


class MatrixTest(unittest.TestCase):
    def test_scoring_both_strands_and_redundancy(self):
        rng = random.Random(12)
        reference = "".join(rng.choices("ACGT", k=200))
        self.assertEqual(score(reference, [reference]), (1.0, 1.0, 1.0))
        self.assertEqual(score(reference, [reverse_complement(reference)]), (1.0, 1.0, 1.0))
        self.assertEqual(score(reference, [reference, reference]), (1.0, 1.0, 2.0))
        self.assertEqual(score(reference, ["A" * 20]), (0.0, 0.0, 0.0))

    def test_generator_is_repeatable_and_accounts_for_template_bases(self):
        rng = random.Random(10)
        reference = "".join(rng.choices("ACGT", k=800))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.fa"
            first = make_reads(path, reference, 5, 500, 0, random.Random(5))
            reads = list(fasta(path))
            self.assertEqual(first[0], len(reads))
            self.assertEqual(first[1], sum(map(len, reads)))
            self.assertGreaterEqual(first[1], 5 * len(reference))
            self.assertTrue(all(read in reference or reverse_complement(read) in reference for read in reads))
            second = make_reads(path, reference, 5, 500, 0, random.Random(5))
            self.assertEqual(first, second)
            self.assertEqual(reads, list(fasta(path)))


if __name__ == "__main__":
    unittest.main()
