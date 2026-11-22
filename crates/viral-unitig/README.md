# viral-unitig

`viral-unitig` is an **experimental research baseline**, not a validated
metagenomic assembler. It counts canonical k-mers from FASTA or FASTQ reads,
removes k-mers below an absolute count and a local relative-support threshold,
and emits maximal non-branching de Bruijn graph paths (unitigs). It does not
require a virus database or a minimum contig length designed for bacteria.

```sh
cargo build --release -p viral-unitig
set -o pipefail
gzip -dc reads.fastq.gz | target/release/viral-unitig > unitigs.fa
# Or use one or more uncompressed FASTA/FASTQ files:
target/release/viral-unitig sample1.fq sample2.fq > unitigs.fa
```

The defaults are `k=21`, minimum count 3, and minimum branch support 0.1.
Only A/C/G/T participate in k-mers; IUPAC ambiguous bases break k-mers.
Reads shorter than k are accepted but contribute no k-mers. Both strands are
combined, and FASTQ qualities are checked for format but not used in assembly.
Input can be wrapped FASTA or FASTQ, but gzip must be decompressed externally.
For available settings, run `viral-unitig --help`.

The minimum branch support is the k-mer count divided by the highest k-mer
count at each of its two junctions, calculated before filtering. It removes
an edge if it falls below the threshold at either junction. This can suppress
low-frequency strains or real low-abundance organisms where their graph paths
overlap abundant ones. Lower `--min-fraction` to retain more alternatives,
at the cost of more sequencing-error branches. `--min-count` similarly trades
low-coverage recovery against error suppression. K-mer counts count
occurrences, not independent reads; repeats can inflate support.

Unitigs stop at branches. Alternative alleles are **not phased into complete
strain contigs**, and individual unitigs are not assigned to organisms.
Different strains, repeats, sequencing errors and related species can all
fragment assemblies or produce misleading paths. Circular paths are marked
`circular=true` and include the duplicated k-1 base overlap. There is no
consensus polishing, base-quality filtering, host-read removal, error model,
coverage normalization, graph export or long-read scaffolding.

This single-threaded in-memory prototype keeps even singleton k-mers until
input has been counted. `--max-kmers` limits distinct k-mers and fails instead
of silently discarding them; it is **not a memory limit**. Graph construction
uses additional memory, and the 5-billion-base / 64-GB target has not been
demonstrated. The defaults and algorithm have not been evaluated against
existing assemblers, real nanopore reads or mixed-infection truth sets.
Do not interpret missing unitigs as absence of an organism. Independent
benchmarking and validation are required before drawing biological conclusions.

## Synthetic smoke benchmark

From the repository root, run:

```sh
python3 crates/viral-unitig/benchmarks/synthetic.py target/release/viral-unitig
```

This simulates a single random 3 kb sequence, uniformly sampled 250-base
fragments, and independent uniform substitutions, insertions and deletions.
It reports nominal depth, unitig count and the longest unitig exactly matching
the reference on either strand. It does not model nanopore error distributions,
repeats, mixed infections, bacterial reads, or large datasets.

For example, one run at 100x nominal depth with 2% substitutions and 2% indel
events produced 449 unitigs and a longest exact unitig of 219 bases; at 1000x
it produced 16,992 unitigs and a longest exact unitig of 2,878 bases. This
illustrates both fragmentation and high-depth spurious output, not superiority
over other assemblers. Approximate-overlap and repeat-graph approaches, such
as Flye's, are possible future directions; they are not implemented here.
