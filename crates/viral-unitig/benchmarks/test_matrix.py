import csv
import io
import random
import tempfile
import unittest
from pathlib import Path

from coverage_bias import make_profile_reads, regional_recall
from matrix import fasta, make_reads, reverse_complement, score
from short_flye_compare import write_result


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

    def test_biased_generator_retains_the_requested_regional_contrast(self):
        rng = random.Random(12)
        reference = "".join(rng.choices("ACGT", k=10_000))
        with tempfile.TemporaryDirectory() as directory:
            reads, bases, regional = make_profile_reads(
                Path(directory) / "reads.fa", reference, 500, 20, 200, 0,
                random.Random(42),
            )
            self.assertGreater(reads, 0)
            self.assertGreaterEqual(bases, 1_100_000)
            self.assertGreater(min(regional), 10)
            self.assertGreater(max(regional) / min(regional), 5)

    def test_regional_recall_scores_a_complete_assembly(self):
        rng = random.Random(12)
        reference = "".join(rng.choices("ACGT", k=10_000))
        self.assertEqual(regional_recall(reference, [reference], 4), [1.0] * 4)

    def test_short_flye_result_reports_regional_recall(self):
        rng = random.Random(12)
        reference = "".join(rng.choices("ACGT", k=1_000))
        output = io.StringIO()
        writer = csv.writer(output)
        write_result(
            writer, 1, reference, 10, len(reference), "biased", "20;200",
            "Flye", "short-ont", "ok", 1.0, [reference],
        )
        row = next(csv.reader(output.getvalue().splitlines()))
        self.assertEqual(row[-1], "1.0000;1.0000")


if __name__ == "__main__":
    unittest.main()
