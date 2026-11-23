# b2b-rusty-toolbox

This repo is expected to eventually contain a miscellaneous set of bioinformatics tools written in Rust. These are not general-purpose tools -- there are plenty of toolkits (`seqtk`, `seqkit`, `bedtools`, etc) for this purpose. The utilities in this toolkit will perform specific, possibly unusual, tasks that the author has a need for and which need to be fast/efficient.

For now, consider this a personal Rust playground.

**NOTE:** Most (possibly all at this point) of the tools in this toolbox are AI-authored and human-reviewed/tested; USE AT YOUR OWN RISK.

## Tools

### `n-every`

Insert an `N` after every K bases in FASTA input. See `n-every --help`.

### `bg-window`

Shrink a bedgraph using a per-sequence adaptive sliding window (min/max/median/mean/sum) so each sequence yields a user-specified range of entries. See `bg-window --help`.

### `porechop-rs`

Trim Oxford Nanopore adapters, split reads with internal adapters, and demultiplex barcodes. This is a mostly drop-in Rust re-implementation of [Porechop](https://github.com/rrwick/Porechop) (thank you Ryan Wick!). It accepts FASTA/FASTQ files (including gzip) and directories of FASTQ files and preserves the original trimming and barcode options; progress messages differ slightly. End-trimming should be identical to the original Porechop (although no guarantees are made). The middle adapter search has been refactored for speed (see below), although results are *nearly* identical to the original behavior in our testing.

When middle adapter searching is performed, this implementation runs roughly 5x faster than the original Porechop in our hands; your mileage may vary.

```sh
porechop-rs -i reads.fastq.gz -o trimmed.fastq.gz -t 8
```
To customize adapter and barcode sequences, export the defaults to a new JSON file, edit it, and pass it with `--adapter_file`:

```sh
porechop-rs --write_default_adapters adapters.json
porechop-rs -i reads.fastq.gz -o trimmed.fastq.gz --adapter_file adapters.json
```

The file replaces the defaults. Each adapter has a unique `name` and at least one of `start` and `end`, each containing a `name` and DNA `sequence` (A, C, G, T or N). Barcode set names start with `Barcode ` and include `(forward)` or `(reverse)`; both ends are required, and reverse barcode start names contain `_rev`. The `barcode_contexts` entries specify flanking sequences for the original full native and rapid barcode adapters; their expansion uses the original kit and barcode naming conventions.

Middle-adapter searches use a four-exact-seed heuristic followed by Parasail alignment with default scoring and a middle identity threshold of at least 90%. The heuristic can miss adapters even above the threshold if all seeds are disrupted, and alignments near the threshold or split coordinates may differ from Porechop. Thresholds below 90%, custom scoring, and adapters shorter than nine bases use the original SeqAn alignment; adapter discovery and end trimming always use SeqAn. To avoid the heuristic, use a threshold below 90% or custom scoring.

This crate is **GPL-3.0-or-later**; see [`crates/porechop-rs/NOTICE.md`](crates/porechop-rs/NOTICE.md) for attribution and bundled third-party licenses. Include the notice and license files when distributing the executable.

### `ir-chimera`

Filter nanopore FASTA/FASTQ reads by inverted-repeat self-alignment. Reads with a plausible inverted seed are aligned to themselves using an in-process minimap2 `ava-ont` index, searching only the reverse strand. By default, reads without a qualifying alignment are written to stdout; `--chimeric` selects the others.

A qualifying alignment has two non-overlapping arms, each spanning at least 500 bases (by default), and at least 80% identity (by default) (exact A/C/G/T matches divided by alignment columns, including gaps). Set these thresholds with `--min-length` and `--min-identity` (percent). Alignment portions crossing the self-alignment diagonal are split into non-overlapping arms before scoring. Ambiguous bases do not count as matches. Reads without a non-overlapping inverted 15-mer pair are retained without alignment; this and minimap2's heuristic seeding favor retaining questionable reads. Detection is not exhaustive, and biological inverted repeats can also be detected. Alignments can extend a shorter inverted repeat into flanking sequence, including the intervening sequence aligned to itself, and thereby exceed the thresholds; tune cutoffs against representative data.

Optionally, `--min-fraction F` (0 < F <= 1) also accepts an alignment when the sum of bases spanned by its two non-overlapping arms is at least `F` times the full read length. Either this fractional threshold or the per-arm `--min-length` threshold can qualify a read, whichever is easier to meet; identity still applies. This allows for the graceful handling of short read distributions.

```sh
ir-chimera --chimeric --min-length 1000 --min-identity 85 reads.fasta > chimeric.fasta
  ```

Input can also be piped on stdin or supplied as multiple plain or gzip-compressed files of the same format. FASTQ must use four lines per record; output preserves the input record format, headers, sequences, and qualities. Reads are processed in bounded batches in input order; `--threads` controls worker count (by default, available CPUs up to 8). Building requires a C compiler and zlib development libraries.

Use `--log summary.tsv` to write a two-line tab-delimited summary after successful processing. Its columns are `n_total`, `n_kept`, and `fraction` (reads written divided by reads consumed, or 0.0 for empty input). Counts include all input files and reflect the selected output mode. The log file must not already exist, to avoid overwriting input or other files.

Use `--alignments alignments.txt` to write one qualifying alignment for every detected chimeric read, independently of whether `--chimeric` is selected. Normal FASTA/FASTQ output still goes to stdout. The report identifies reads by the sequence ID before the first whitespace and shows only the non-overlapping portion used for classification, in 60-column blocks with 1-based inclusive read coordinates. The target is shown in the forward orientation, and the query is reverse-complemented with decreasing coordinates. `|` marks exact A/C/G/T matches, spaces indicate mismatches or ambiguous bases, and `-` indicates a gap. Identity includes gaps and ambiguous bases. The report file must not already exist and may be partial if processing fails. Only the first qualifying alignment is reported, not necessarily the best or all possible alignments.
