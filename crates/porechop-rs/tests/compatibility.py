"""Differential tests against the original executable.

Set PORECHOP_REFERENCE_ROOT to a checkout of the original implementation
with its C++ component built. Run after `cargo build --release -p porechop-rs`:
python3 crates/porechop-rs/tests/compatibility.py
"""
import gzip
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(os.environ["PORECHOP_REFERENCE_ROOT"]).resolve()
ORIGINAL = [sys.executable, str(ROOT / "porechop-runner.py")]
RUST = [os.environ.get("PORECHOP_RS_BIN", str(
    Path(__file__).resolve().parents[3] / "target/release/porechop-rs"))]


def run(executable, *args):
    p = subprocess.run(executable + [str(x) for x in args],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if p.returncode:
        raise AssertionError("{} failed:\n{}".format(executable, p.stderr.decode()))
    return p.stdout


def contents(path):
    if path.suffix == ".gz":
        with gzip.open(str(path), "rb") as stream:
            return stream.read()
    return path.read_bytes()


class Compatibility(unittest.TestCase):
    def test_reads_and_options(self):
        examples = [
            ("test_one_adapter_set.fastq", []),
            ("test_two_adapter_sets.fastq", []),
            ("test_one_adapter_set.fastq", ["--no_split"]),
            ("test_one_adapter_set.fastq", ["--discard_middle"]),
            ("test_one_adapter_set.fastq", ["--check_reads", "1"]),
            ("test_one_adapter_set.fastq", ["--check_reads", "2"]),
            ("test_one_adapter_set.fastq", ["--end_size", "50"]),
            ("test_one_adapter_set.fastq", ["--end_size", "200"]),
            ("test_one_adapter_set.fastq", ["--extra_end_trim", "10"]),
            ("test_one_adapter_set.fastq", ["--extra_end_trim", "100"]),
            ("test_one_adapter_set.fastq", ["--adapter_threshold", "80"]),
            ("test_one_adapter_set.fastq", ["--scoring_scheme", "2,-3,-4,-1"]),
            ("test_one_adapter_set.fastq", ["--check_reads", "0", "--scoring_scheme=-1,-6,-5,-2",
                                            "--end_threshold", "-1"]),
            ("test_one_adapter_set.fastq", ["--min_trim_size", "6"]),
            ("test_one_adapter_set.fastq", ["--middle_threshold", "96"]),
            ("test_one_adapter_set.fastq", ["--middle_threshold", "97"]),
            ("test_one_adapter_set.fastq", ["--min_split_read_size", "1"]),
            ("test_one_adapter_set.fastq", ["--min_split_read_size", "429"]),
            ("test_one_adapter_set.fastq", ["--extra_middle_trim_good_side", "0",
                                            "--extra_middle_trim_bad_side", "0"]),
            ("test_one_adapter_set.fastq", ["--extra_middle_trim_good_side", "123",
                                            "--extra_middle_trim_bad_side", "456"]),
            ("test_format.fasta", ["--format", "fastq"]),
            ("test_format.fastq", ["--format", "fasta"]),
            ("test_format.fastq", ["--format", "fasta.gz"]),
            ("test_format.fasta.gz", []),
            ("test_format.fastq.gz", []),
        ]
        for fixture, extra in examples:
            with self.subTest(fixture=fixture, extra=extra):
                args = ["-i", ROOT / "test" / fixture, "-v", "0", "-t", "1"] + extra
                self.assertEqual(run(ORIGINAL, *args), run(RUST, *args))

    def test_barcode_binning(self):
        for fixture in ("test_barcodes.fastq", "test_albacore_directory"):
            for extra in ([], ["--untrimmed"], ["--require_two_barcodes"],
                          ["--discard_unassigned"], ["--format", "fasta.gz"]):
                with self.subTest(fixture=fixture, extra=extra), tempfile.TemporaryDirectory() as tmp:
                    old = Path(tmp) / "old"
                    new = Path(tmp) / "new"
                    args = ["-i", ROOT / "test" / fixture, "-v", "0", "-t", "1"] + extra
                    run(ORIGINAL, *(args + ["-b", old]))
                    run(RUST, *(args + ["-b", new]))
                    self.assertEqual(sorted(p.name for p in old.iterdir()),
                                     sorted(p.name for p in new.iterdir()))
                    for path in old.iterdir():
                        self.assertEqual(contents(path), contents(new / path.name))

    def test_compressed_output_and_custom_config(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            config = tmp / "adapters.json"
            run(RUST, "--write_default_adapters", config)
            original_config = tmp / "original-adapters.json"
            run(ORIGINAL, "--write_default_adapters", original_config)
            self.assertEqual(json.loads(config.read_text()), json.loads(original_config.read_text()))
            config_data = json.loads(config.read_text())
            config_data["adapters"] = [config_data["adapters"][0]]
            config.write_text(json.dumps(config_data))
            args = ["-i", ROOT / "test/test_one_adapter_set.fastq", "-v", "0", "-t", "1",
                    "--adapter_file", config]
            old, new = tmp / "old.fastq.gz", tmp / "new.fastq.gz"
            run(ORIGINAL, *(args + ["-o", old]))
            run(RUST, *(args + ["-o", new]))
            self.assertEqual(contents(old), contents(new))
            with self.assertRaises(AssertionError):
                run(RUST, "--write_default_adapters", config)

    def test_custom_barcode_and_contexts(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            config = tmp / "adapters.json"
            run(RUST, "--write_default_adapters", config)
            settings = json.loads(config.read_text())
            settings["adapters"] = [a for a in settings["adapters"] if a["name"] in (
                "SQK-NSK007", "Barcode 1 (reverse)", "Barcode 2 (reverse)")]
            barcode = next(a for a in settings["adapters"] if a["name"] == "Barcode 1 (reverse)")
            barcode["start"]["sequence"] = "GCTAGCTAGTCGATCGTAGC"
            barcode["end"]["sequence"] = "GCTACGATCGACTAGCTAGC"
            settings["barcode_contexts"]["native_start_prefix"] = "ACGTCGACTGACTGAC"
            settings["barcode_contexts"]["native_start_suffix"] = "TACGATCGA"
            settings["barcode_contexts"]["native_end_prefix"] = "TCGATCGTA"
            settings["barcode_contexts"]["native_end_suffix"] = "GTCAGTCAGTCGACGT"
            config.write_text(json.dumps(settings))
            kit = settings["adapters"][0]
            ctx = settings["barcode_contexts"]
            seq = (kit["start"]["sequence"] + ctx["native_start_prefix"]
                   + barcode["start"]["sequence"] + ctx["native_start_suffix"]
                   + "GC" * 500 + ctx["native_end_prefix"] + barcode["end"]["sequence"]
                   + ctx["native_end_suffix"] + kit["end"]["sequence"])
            reads = tmp / "reads.fastq"
            reads.write_text("@read\n{}\n+\n{}\n".format(seq, "I" * len(seq)))
            old, new = tmp / "old", tmp / "new"
            args = ["-i", reads, "-v", "0", "-t", "1", "--no_split", "--adapter_file", config]
            run(ORIGINAL, *(args + ["-b", old]))
            run(RUST, *(args + ["-b", new]))
            self.assertTrue((old / "BC01.fastq").exists())
            self.assertEqual(contents(old / "BC01.fastq"), contents(new / "BC01.fastq"))

    def test_short_reads_and_rna(self):
        rng = random.Random(7)
        adapter = "AATGTACTTCGTTCAGTTACGTATTGCT"
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "reads.fastq"
            records = []
            for index in range(60):
                length = rng.randrange(3, 100)
                sequence = "".join(rng.choice("ACGT") for _ in range(length))
                if index % 3 == 0:
                    sequence = adapter + sequence
                elif index % 3 == 1:
                    sequence += adapter
                if index % 7 == 0:
                    sequence = sequence.replace("T", "U")
                records.append("@read{} description\n{}\n+\n{}\n".format(
                    index, sequence, "I" * len(sequence)))
            path.write_text("".join(records))
            for options in (["--no_split"], ["--no_split", "--format", "fasta"]):
                with self.subTest(options=options):
                    args = ["-i", path, "-v", "0", "-t", "1", "--check_reads", "60"] + options
                    self.assertEqual(run(ORIGINAL, *args), run(RUST, *args))

    def test_verbose_diagnostics_and_input_errors(self):
        fixture = ROOT / "test/test_one_adapter_set.fastq"
        result = subprocess.run(RUST + ["-i", str(fixture), "-v", "3", "-t", "1"],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, run(ORIGINAL, "-i", fixture, "-v", "0", "-t", "1"))
        for expected in (b"start alignments:", b"end alignments:", b"read coords:"):
            self.assertIn(expected, result.stderr)
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / "bad.fastq"
            bad.write_text("@read\nACGT\n+\n")
            result = subprocess.run(RUST + ["-i", str(bad), "-v", "0"], stderr=subprocess.PIPE)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(b"truncated FASTQ", result.stderr)


if __name__ == "__main__":
    unittest.main()
