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
target/release/viral-unitig --uncertain uncertain.fa sample1.fq sample2.fq > unitigs.fa
```

The defaults are `k=21`, minimum count 10, minimum branch support 0.4, and
coverage-relative error cleaning enabled. **These aggressive defaults can
discard low-coverage organisms and collapse or erase real minority strains.**
For a less aggressive investigation, use `--min-count 3 --min-fraction 0.2`,
at the cost of more error-derived fragments.
Only A/C/G/T participate in k-mers; IUPAC ambiguous bases break k-mers.
Reads shorter than k are accepted but contribute no k-mers. Both strands are
combined, and FASTQ qualities are checked for format but not used in assembly.
Input can be wrapped FASTA or FASTQ, but gzip must be decompressed externally.
For available settings, run `viral-unitig --help`.

The error-k-mer searches and short-unitig redundancy checks use worker
threads; reading and counting k-mers, graph construction, traversal, and
output remain serial. By default the worker pool uses the available CPU
count, capped at 8. Use `--threads N` (or `-t N`) to set a positive thread
count, including `--threads 1` for a single worker. More workers can increase
memory consumption, particularly when checking short unitigs. Speedup depends
on the dataset and on the time spent in the parallel stages.

The minimum branch support is the k-mer count divided by the highest k-mer
count at each of its two junctions, calculated before filtering. It removes
an edge if it falls below the threshold at either junction. This can suppress
low-frequency strains or real low-abundance organisms where their graph paths
overlap abundant ones. Lower `--min-fraction` to retain more alternatives,
at the cost of more sequencing-error branches. `--min-count` similarly trades
low-coverage recovery against error suppression. K-mer counts count
occurrences, not independent reads; repeats can inflate support.

The default cleaner compares low-count k-mers to much more abundant similar
k-mers, including substitution and indel-like neighbors. It removes weak
k-mers only when their extended local graph contexts are similar. A second
pass also considers two-substitution neighbors. After graph traversal, short
unitigs (`<= 2*k+2` bases) that closely match a much better-supported path are
treated as uncertain. Without `--uncertain`, they are **discarded** and their
count is reported on stderr. With `--uncertain FILE`, they are written to that
new file instead; the file must not already exist. These matches are not proof
that a fragment is an error. In complex metagenomes, a real virus or strain can
be filtered or routed to the uncertain file, including by chance similarity.
The uncertain file does not capture k-mers removed before unitig generation.

Every primary and uncertain FASTA header includes `mean_kmer_depth`, the
arithmetic mean of the input occurrence counts of the canonical k-mers along
the emitted unitig, rounded to two decimal places. For circular unitigs the
duplicated k-1 overlap does not add extra k-mer positions. Counts combine both
strands and all input files and saturate at 2^32-1 per distinct k-mer.
This is a rough relative-support signal, **not a count of reads or an unbiased
estimate of base coverage or organism abundance**: repeated sequences, shared
k-mers, read length, errors, uneven coverage and filtering can bias it.
Discarded or corrected error k-mers are not reassigned to the retained path.

Use `--no-clean-errors --min-fraction 0 --min-count 3` to disable both the
near-neighbor cleaner and relative branch filtering when investigating
low-frequency alternatives, while retaining a minimum count of 3. Cleaning
does not use base qualities or read linkage and can collapse true closely
related populations.
The short-unitig search checks at most 50,000 neighbors per unitig; on reaching
that limit it retains the unitig and reports a warning.

`--clean-bubbles` enables an additional **experimental, opt-in** cleanup.
Before edge filtering, it looks for short unbranched weak paths that diverge
and rejoin a stronger path, and removes the entire weak path if its maximum
k-mer count is below `--min-fraction` times the stronger path's minimum count.
The stronger path is followed greedily; paths are limited to `2*k+2` edges and
must differ by at most two edits. Counts of removed paths and k-mers are
reported on stderr. Incomplete bubbles and unanchored fragments are left alone.

This evidence cannot distinguish an error from a true low-frequency strain
or related virus: the option can eliminate real variation. Disabling it does
not disable the other cleanup and filtering steps.

Unitigs stop at branches. Alternative alleles are **not phased into complete
strain contigs**, and individual unitigs are not assigned to organisms.
Different strains, repeats, sequencing errors and related species can all
fragment assemblies or produce misleading paths. Circular paths are marked
`circular=true` and include the duplicated k-1 base overlap. There is no
consensus polishing, base-quality filtering, host-read removal, error model,
coverage normalization, graph export or long-read scaffolding.

This in-memory prototype keeps even singleton k-mers until
input has been counted. `--max-kmers` limits distinct k-mers and fails instead
of silently discarding them; it is **not a memory limit**. Graph construction
uses additional memory. Neighbor searches and buffering short unitigs can be
expensive. The 5-billion-base / 64-GB target has not been demonstrated, and
there has been no comparison against existing assemblers or validation on real
nanopore reads or mixed-infection truth sets.
Do not interpret missing unitigs as absence of an organism. Independent
benchmarking and validation are required before drawing biological conclusions.

## Synthetic smoke benchmark

From the repository root, run:

```sh
python3 crates/viral-unitig/benchmarks/synthetic.py target/release/viral-unitig
python3 crates/viral-unitig/benchmarks/mixed.py target/release/viral-unitig
python3 crates/viral-unitig/benchmarks/strains.py target/release/viral-unitig
# Compare the previous, less aggressive settings:
python3 crates/viral-unitig/benchmarks/strains.py target/release/viral-unitig --min-count 3 --min-fraction 0.2
```

This simulates a single random 3 kb sequence, uniformly sampled 250-base
fragments, and independent uniform substitutions, insertions and deletions.
It reports nominal depth, unitig count and the longest unitig exactly matching
the reference on either strand, and the number of fragments written to the
uncertain file. It does not model nanopore error distributions, repeats,
bacterial reads, or large datasets.

In the 1000x case with 2% substitutions and 2% indel events, the default
produced one 2,985-base exact primary unitig, while 1,398 short unitigs were
classified as uncertain. At 100x with the same errors, it produced five
primary unitigs, with a longest exact unitig of 2,876 bases.

The second benchmark simulates two 1 kb genomes using 200-base fragments and
1% substitutions plus 1% indel events. A 1000x genome mixed with an unrelated
50x genome produced three primary unitigs, with longest exact matches of
995 and 895 bases, respectively. Two related genomes at 600x and 300x produced
19 primary unitigs: all six alternate sites were represented, but the longest
exact match to either genome was only 149 bases. This does **not** demonstrate
strain phasing or superiority over other assemblers. Approximate-overlap and
repeat-graph approaches, such as Flye's, are not implemented.

The third benchmark simulates two 10 kb genomes differing at uniformly
sampled substitutions, with 1,000x and 200x nominal coverage, 250-base
fragments, and independent 2% substitutions and 2% insertion/deletion events.
It reports primary and uncertain unitig counts, longest exact matches and
recovery of strain-specific 21-mer positions. These k-mers can be scattered
across unitigs and do **not** establish complete or correctly phased genomes.

For these fixed synthetic seeds, the defaults trade minority recovery for
longer major-strain unitigs compared with `--min-count 3 --min-fraction 0.2`:

| Divergence | Settings | Primary unitigs | Longest exact major unitig | Minor-specific 21-mer recovery |
|---|---|---:|---:|---:|
| 10% | less aggressive | 528 | 578 bp | 93.38% |
| 10% | defaults | 162 | 9,990 bp | 62.07% |
| 5% | less aggressive | 674 | 267 bp | 86.70% |
| 5% | defaults | 78 | 9,987 bp | 25.81% |

Neither setting produces a complete minor-strain contig in this simulation.
Real errors and population structure may behave very differently.
