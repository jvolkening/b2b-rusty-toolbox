# b2b-rusty-toolbox

This repo is expected to eventually contain a miscellaneous set of bioinformatics tools written in Rust. These are not general purpose tools -- there are plenty of toolkits (`seqtk`, `seqkit`, `bedtools`, etc) for this purpose. The utilities in this toolkit will perform specific, probably unusual, tasks that the author has a need for and which need to be fast/efficient.

More documentation may come at some point. For now, consider this a personal Rust playground.

## Tools

- `n-every`: insert an `N` after every K bases in FASTA input.
- `bg-window`: shrink a bedgraph using a per-sequence adaptive sliding window (min/max/median/mean/sum) so each sequence yields a user-specified range of entries. See `bg-window --help`.
- `ir-chimera`: filter nanopore FASTA/FASTQ reads by inverted-repeat self-alignment. Reads with a plausible inverted seed are aligned to themselves using an in-process minimap2 `ava-ont` index, searching only the reverse strand. By default, reads without a qualifying alignment are written to stdout; `--chimeric` selects the others.

  A qualifying alignment has two non-overlapping arms, each spanning at least 500 bases, and at least 80% identity (exact A/C/G/T matches divided by alignment columns, including gaps). Set these thresholds with `--min-length` and `--min-identity` (percent). Alignment portions crossing the self-alignment diagonal are split into non-overlapping arms before scoring. Ambiguous bases do not count as matches. Reads without a non-overlapping inverted 15-mer pair are retained without alignment; this and minimap2's heuristic seeding favor retaining questionable reads. This is not exhaustive, and biological inverted repeats can also be detected.

  ```sh
  cargo build --release -p ir-chimera
  target/release/ir-chimera reads.fastq.gz > nonchimeric.fastq
  target/release/ir-chimera --chimeric --min-length 1000 --min-identity 85 reads.fasta > chimeric.fasta
  ```

  Input can also be piped on stdin or supplied as multiple plain or gzip-compressed files of the same format. FASTQ must use four lines per record; output preserves the input record format, headers, sequences, and qualities. Reads are processed in bounded batches in input order; `--threads` controls worker count (by default, available CPUs up to 8). Building requires a C compiler and zlib development libraries.
