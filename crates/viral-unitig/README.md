# viral-unitig

`viral-unitig` is an **experimental research baseline**, not a validated
metagenomic assembler. It counts canonical k-mers from FASTA or FASTQ reads,
removes k-mers below an absolute count and a local relative-support threshold,
and emits maximal non-branching de Bruijn graph paths (unitigs). It does not
require a virus database or a minimum contig length designed for bacteria.
Unlike the short-read mode in the
[Flye fork](https://github.com/jvolkening/Flye/tree/short-ont-metagenomes),
it does not build an overlap/repeat graph or use mapped reads for end
extension. Those approaches can exploit linkage across noisy k-mers; this
prototype instead favors a smaller, coverage-filtered graph.

```sh
cargo build --release -p viral-unitig
set -o pipefail
gzip -dc reads.fastq.gz | target/release/viral-unitig > unitigs.fa
# Or use one or more uncompressed FASTA/FASTQ files:
target/release/viral-unitig --uncertain uncertain.fa sample1.fq sample2.fq > unitigs.fa
```

The defaults are `k=21`, minimum count 3, minimum branch support 0.4, and
coverage-relative error cleaning enabled. **These filters can
discard low-coverage organisms and collapse or erase real minority strains.**
For a less aggressive investigation, use `--min-fraction 0.2`; this retains
more alternatives but can add error-derived fragments and does not guarantee
better recall or contiguity.
For very low-depth samples, shorter k-mers can recover more sequence, but are
less specific and may increase graph branching, misjoins and error-derived
contigs. In a five-seed 10x sweep, `--kmer 15 --min-count 3` gave better
contiguity than `--kmer 15 --min-count 2`, with somewhat less recall; neither
setting guarantees recovery at 5x.
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

`--normalize-depth 50` enables experimental, order-dependent streaming
normalization. A read is skipped if at least 10% of its valid k-mer
occurrences already have count 50 or more in previously retained reads,
and no run of 200 consecutive valid k-mer occurrences falls below that count.
The threshold is an occurrence count, not read or base depth. Zero (the
default) disables normalization. The number of skipped reads is reported on
stderr. This can reduce high-depth time and memory, but shared sequences,
repeats and early abundant reads can cause reads from a real rare organism or
strain to be skipped, including entire later populations. Even a single
organism can lose short novel extensions when its overlapping reads are
ordered. The 200-k-mer guard protects only long undercovered runs, not all
new sequence; ambiguous bases do not reset this run. Counts and
`mean_kmer_depth` then describe retained reads only, not the original sample.
It is not a substitute for quality-aware correction or a guarantee of
uniform coverage. Keep the original reads and compare with normalization
disabled when sensitivity matters.

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
graph export or long-read scaffolding.

This in-memory prototype keeps even singleton k-mers until
input has been counted. `--max-kmers` limits distinct k-mers and fails instead
of silently discarding them; it is **not a memory limit**. Graph construction
uses additional memory. Neighbor searches and buffering short unitigs can be
expensive. The 5-billion-base / 64-GB target has not been demonstrated, and
there has been no comparison against existing assemblers or validation on real
nanopore reads or mixed-infection truth sets.
Do not interpret missing unitigs as absence of an organism. Independent
benchmarking and validation are required before drawing biological conclusions.

## Synthetic size and coverage matrix

```sh
cargo build --release -p viral-unitig
python3 crates/viral-unitig/benchmarks/matrix.py target/release/viral-unitig \
  --modes aggressive sensitive normalized --threads 2 --timeout 15 > results.csv
```

The default matrix uses one independently generated random linear genome per
size (800, 3,000, 10,000, 30,000 and 100,000 bp), each at nominal 5, 10,
100, 1,000 and 5,000x. Read lengths are sampled from a normal distribution
(mean 500 bp, standard deviation 125 bp, minimum 100 bp before end
clipping), and both strands are sampled. Per reference base, independent
events occur at 5%: half substitutions, one quarter insertions and one
quarter deletions. The seed is fixed; nominal depth is sampled reference bases
divided by genome length. Cases run separately; each mode sees identical reads.
`--sizes`, `--depths`, `--mean`, `--error`, `--seed`, `--timeout` and
`--threads` permit additional experiments. The benchmark includes
`aggressive` (`--min-count 10`), `sensitive` (`--min-count 3`), `normalized`
(`--min-count 3 --normalize-depth 50`) and `short-k` (`k=15`,
`--min-count 2 --normalize-depth 50`) modes, plus `fast` (the normalized
settings with near-neighbor cleaning disabled) and `normalized-25` (a
normalization threshold of 25).

The CSV reports contig count, longest contig, reference 31-mer recall, contig
31-mer precision and matching 31-mer redundancy (matching contig occurrences
divided by distinct reference positions recovered). A contig shorter than
31 bp contributes to count and length but not to those k-mer metrics.
These are exact-match proxies, not alignment-based identity, genome breadth,
misassembly or strain-phasing measurements. In particular, disconnected
fragments can yield high recall, erroneous regions are penalized over a
31-base window, and repeat positions can distort redundancy. Timeout and
nonzero exit rows have no assembly metrics; they are not zero-recovery results.
Generation and scoring are excluded from assembler elapsed time. The 15-second
limit is a benchmark budget, not the requested one-hour target.

The checked-in [matrix results](benchmarks/matrix-results.csv) and
[candidate results](benchmarks/candidate-results.csv) used an earlier 20%
normalization cutoff; the command above uses the current 10% cutoff. The
[10% results](benchmarks/normalization-10pct-results.csv) used 10% without
an undercoverage guard; [guard-100](benchmarks/guard100-results.csv) and
[guard-200](benchmarks/guard200-results.csv) use 100- and 200-k-mer guards,
respectively. These runs used:

```sh
python3 crates/viral-unitig/benchmarks/matrix.py target/release/viral-unitig \
  --sizes 800 10000 100000 --depths 100 1000 5000 \
  --modes normalized normalized-25 --threads 2 --timeout 30
```

These ran on a two-core Intel Xeon Platinum 8370C, with 7.5 GiB RAM and
Rust 1.100.0. A single random sequence
and seed per size cannot represent repeats, real nanopore systematic errors,
coverage bias, chimeras, contamination, mixed taxa, or related strains. A
comparison with the experimental short-ONT Flye branch is now included below,
but there is no comparison with MEGAHIT, real nanopore data, mixed-taxon truth
sets, or five million reads. Performance and biological utility in those
settings remain unverified.

In the earlier matrix, at 5x the aggressive settings recovered no reference 31-mers
for any size, while the normalized settings recovered 2.9% to 14.1% in
fragmented unitigs. At 100x, normalized settings recovered 99.97% or more of
reference 31-mers in one or two unitigs. At 5,000x and 10 kb, the aggressive
run exceeded the default 10-million-distinct-k-mer limit; normalization
completed in 6.07 seconds with one 10-kb contig and full 31-mer recall.
However, normalization also exceeded that limit at 100 kb and 1,000x and
5,000x, and the 30-kb 5,000x run timed out at 15 seconds. These results
demonstrate both the improvement and substantial remaining limitations;
they do not establish completeness or accuracy on real samples.

At 100 kb and 5,000x, the unguarded 10% rule with threshold 25 completed in
28.62 seconds with two contigs and 99.99% reference 31-mer recall; threshold
50 exceeded the distinct-k-mer limit. At 100 kb and 100x, threshold 25 gave
eight contigs, while threshold 50 gave one. Earlier 20% tests of `short-k`
improved 5x 31-mer recall but fragmented 100-kb high-depth assemblies, and
disabling error cleanup produced thousands of extra contigs at high depth.
These tradeoffs are not evidence of performance on real data.

The final 200-k-mer guard at threshold 25 produced one 100-kb contig with
full reference 31-mer recall at 5,000x in 28.68 seconds; threshold 50 exceeded
the k-mer limit. The stricter 100-k-mer guard exceeded the limit at 100 kb
and 1,000x and 5,000x. A full run with the final normalization rule is in
[final-matrix.csv](benchmarks/final-matrix.csv), and additional 600-bp and
10-kb cases with mean fragment lengths of
[300](benchmarks/mean300-results.csv) and
[800](benchmarks/mean800-results.csv) bp are included.

```sh
python3 crates/viral-unitig/benchmarks/matrix.py target/release/viral-unitig \
  --modes normalized-25 --threads 2 --timeout 40
python3 crates/viral-unitig/benchmarks/matrix.py target/release/viral-unitig \
  --sizes 600 10000 --depths 5 10 100 --mean 300 \
  --modes sensitive normalized-25 short-k --threads 2 --timeout 30
# Repeat with --mean 800.
```

## Low-depth and uneven-coverage investigations

The [low-depth sweep](benchmarks/low-coverage-sweep.csv) tested k=15, 17, 19,
and 21; minimum counts 1, 2, and 3; relative support 0.1, 0.2, and 0.4; and
error-cleaning on/off at 5x and 10x. Five-seed results for the more promising
settings are in [low-coverage-replicates.csv](benchmarks/low-coverage-replicates.csv).
At 10x, k=15/min-count=3 improves the mean exact 31-mer recall and median
longest unitig over the default k=21/min-count=3, but remains highly
fragmented:

| Genome | k=21/count=3 recall; median contigs; longest | k=15/count=3 recall; median contigs; longest |
|---:|---:|---:|
| 600 bp | 51.3%; 20; 110 bp | 72.4%; 10; 216 bp |
| 3,000 bp | 45.1%; 91; 282 bp | 62.7%; 76; 498 bp |
| 10,000 bp | 44.6%; 295; 255 bp | 60.5%; 230; 632 bp |

These are five-seed means for recall and medians for contigs/longest contig.
Minimum count 1 produced too many error-derived fragments to be useful;
lowering relative support below 0.4 did not improve the promising settings,
and disabling error cleanup increased false fragments without improving
reference recall. The tested weak-bubble cleanup had no effect at 10x.
Consequently, these experiments do not justify weakening the defaults
globally. `k=15, min-count=3` is a candidate sensitivity profile, not a
validated replacement.

The [regional-bias benchmark](benchmarks/coverage-bias-final.csv) compares
uniform controls with profiles having the same mean depth and a 10:1 or 20:1
high-to-low regional depth ratio. Reads have mean length 500 bp and 5%
independent substitution/insertion/deletion events. In the 10-kb, 10:1 case,
the default assembler recovered 92.28% of exact reference 31-mers but emitted
70 contigs (longest 2,926 bp), versus one 10-kb contig from the matched-mean
uniform control. The realized regional depths were 210.4, 16.8, 20.8, and
208.3x. At 20:1, the biased case emitted 67 contigs (longest 3,366 bp) with
92.64% recall; the matched-mean uniform control again emitted one contig.
In both cases the high-depth end regions had full recall, while lower-depth
regions lost recall and contiguity. These are single-seed examples, not a
claim that every coverage gradient fragments this way.

The [threshold/cleanup sweep](benchmarks/coverage-bias-threshold-sweep.csv)
tested whether those losses could be fixed with the existing graph filters.
For the 10-kb 10:1 case, lowering `--min-fraction` from 0.4 to 0.2 changed
70 contigs/92.28% recall to 78/91.98%; lowering it to 0.1 produced 319
contigs/82.53%, and disabling relative filtering produced 7,767
contigs/41.29%. Disabling error cleanup at 0.4 retained the same recall but
produced 3,181 contigs and 79.21% exact 31-mer precision. The tested
k=15/min-count=3 mode improved this case to 36 contigs, a 3,334-bp longest
unitig, and 95.92% recall, but did not make it contiguous. Lowering the
minimum count to 2 also increased fragmentation. The current filters are
therefore not a remedy for strong regional depth imbalance; a locally
adaptive error model or a read-overlap/linkage approach remains necessary.

As a direct overlap-assembly comparison, the
[short-Flye benchmark](benchmarks/short-flye-compare-5seeds.csv) used the
`short-ont-metagenomes` branch at
[`f9cb6d8`](https://github.com/jvolkening/Flye/commit/f9cb6d87086cc9e8f57e6208a71b9c37fdd369af).
On the same five independent, uniform 10x read sets (mean 500 bp, 5% errors),
that branch's short-read overlap mode emitted one contig for every 600-bp and
3-kb case, and one contig in four of five 10-kb cases. Mean exact 31-mer recall
was 98.94%, 99.75%, and 98.96%, respectively; median runtimes were 0.69, 0.77,
and 1.00 seconds using two threads. In this experiment, default
`viral-unitig` had median contig counts of 20, 91, and 295 and mean recall of
51.33%, 45.08%, and 44.58%. These results strongly favor the overlap-based
approach on this particular low-depth synthetic model; they do not establish
performance on real reads or high-depth mixtures. Optional Flye end extension
did not materially change the results.

On three 100-kb, 10x replicates (about 2,000 reads each), Flye produced
2-5 contigs, with 98.79-99.58% exact 31-mer recall, 98.58-99.50% precision,
and 3.14-3.63-second runtimes. In contrast, the default de Bruijn settings
produced about 3,000 fragments and 44-46% recall; k=15/count=3 improved this
to about 2,250 fragments and 60-61% recall. See the
[100-kb comparison](benchmarks/short-flye-100kb-10x-3seeds.csv).

At 5x, Flye was also much more sensitive than the de Bruijn baseline but was
not complete: it failed to emit contigs in one of five 600-bp runs; the other
four had one contig and mean 84.78% recall. For 3-kb genomes it produced a
median of two contigs and 86.09% mean recall; for 10-kb genomes, three contigs
and 78.32% mean recall. On the same seeds, default `viral-unitig` had 9.47%,
10.48%, and 11.29% mean recall and median contig counts of 15, 68, and 217.
Full results are in the [5x comparison](benchmarks/short-flye-compare-5x-5seeds.csv).

On three 10:1 biased-depth replicates with observed low regions of 16.7-30.4x
and high regions of 204-220x, Flye also emitted one contig in every 3-kb and
10-kb case. Exact 31-mer recall was 99.97% and 99.99%, respectively, with
median runtimes of 3.0 and 4.5 seconds. On the same reads, default
`viral-unitig` produced 2-11 contigs at 3 kb and 49-68 at 10 kb; k=15/count=3
reduced those counts to 2-7 and 28-35 but did not make the 10-kb assemblies
contiguous. The 600-bp profiles have substantial read-boundary smoothing, so
their nominal 10:1 target does not produce a 10:1 realized regional depth.
Flye still emitted one contig for each 600-bp genome, but the unextended contigs
recovered 84-99.65% of exact reference 31-mers; `--extend-ends` made no
material difference in these runs. See the
[three-seed regional comparison](benchmarks/short-flye-biased-ratio10-3seeds.csv).
Three 20:1 profiles were also tested, with observed low regions of 18.8-24.4x
and high regions of 409-428x. Flye again produced one contig in every 3-kb and
10-kb case, with 99.97% and 99.99% exact 31-mer recall and median runtimes of
5.2 and 8.7 seconds. Default `viral-unitig` produced 9-16 contigs at 3 kb and
52-80 at 10 kb; k=15/count=3 reduced these to 3-7 and 28-41. At 600 bp, Flye
produced one 599-bp contig in all three runs (99.82% recall). Results are in
the [20:1 regional comparison](benchmarks/short-flye-biased-ratio20-3seeds.csv).
In this Flye branch, short mode selects a minimum overlap from read N90
(clamped to 100-500 bp), while raw short reads use 11-mers and a 3-minimizer
window. This overlap-based read linkage is the strongest remaining algorithmic
avenue to pursue in Rust; lowering k-mer support thresholds alone did not
reproduce its contiguity.

The same overlap mode did not handle extreme depth: on 10,356 reads from a
600-bp genome at 5,000x it failed to assemble any disjointigs in 15.12 seconds.
On those same reads, unnormalized k=21/min-count=3 emitted 84 contigs, whereas
streaming normalization at depth 25 or 50 produced one 600-bp contig with
full exact 31-mer recall and precision in 0.13-0.18 seconds. The matching
[normalization results](benchmarks/short-flye-600-5000x-unitig-modes.csv)
illustrate a useful high-depth path, but are single-genome simulations:
normalization is order-dependent and can discard real low-abundance organisms
or variants.

The sweeps can be reproduced from the repository root:

```sh
python3 crates/viral-unitig/benchmarks/low_coverage_sweep.py \
  target/release/viral-unitig --threads 2 > low-coverage-sweep.csv
python3 crates/viral-unitig/benchmarks/low_coverage_replicates.py \
  target/release/viral-unitig --threads 2 > low-coverage-replicates.csv
python3 crates/viral-unitig/benchmarks/coverage_bias.py \
  target/release/viral-unitig --sizes 3000 10000 --means 500 \
  --floor 20 --ratios 1 10 20 --threads 2 > coverage-bias.csv
python3 crates/viral-unitig/benchmarks/matrix.py \
  target/release/viral-unitig --sizes 600 --depths 5000 \
  --modes sensitive normalized normalized-25 count2-normalized short-k \
  --threads 2 > small-high-depth.csv
python3 crates/viral-unitig/benchmarks/short_flye_compare.py \
  target/release/viral-unitig /path/to/Flye/bin/flye \
  --sizes 600 3000 10000 --seeds 20261124 20261125 20261126 20261127 20261128 \
  --depth 10 --mean 500 --threads 2 > short-flye-compare.csv
python3 crates/viral-unitig/benchmarks/short_flye_compare.py \
  target/release/viral-unitig /path/to/Flye/bin/flye \
  --sizes 600 3000 10000 --seeds 20261124 20261125 20261126 20261127 20261128 \
  --depth 5 --mean 500 --threads 2 > short-flye-5x.csv
python3 crates/viral-unitig/benchmarks/short_flye_compare.py \
  target/release/viral-unitig /path/to/Flye/bin/flye \
  --sizes 100000 --seeds 20261124 20261125 20261126 \
  --depth 10 --mean 500 --threads 2 --skip-extend-ends \
  > short-flye-100kb.csv
python3 crates/viral-unitig/benchmarks/short_flye_compare.py \
  target/release/viral-unitig /path/to/Flye/bin/flye \
  --sizes 600 3000 10000 --seeds 20261124 20261125 20261126 \
  --profile biased --floor 20 --ratio 10 --mean 500 --threads 2 \
  > short-flye-biased.csv
python3 crates/viral-unitig/benchmarks/short_flye_compare.py \
  target/release/viral-unitig /path/to/Flye/bin/flye \
  --sizes 600 3000 10000 --seeds 20261124 20261125 20261126 \
  --profile biased --floor 20 --ratio 20 --mean 500 --threads 2 \
  > short-flye-biased-20to1.csv
```

The Flye result is specific to the cited experimental branch and parameter
set; it is a useful direction for the Rust prototype to investigate, not an
endorsement of that branch as a production metagenomics assembler.

The [mixed benchmark](benchmarks/mixed-guard200-results.csv), with reads
ordered by genome, recovered long unitigs from an unrelated 50x genome
alongside a 1,000x genome, but recovered none of six alternate sites in a
related 300x genome alongside a 600x genome. In the separate 10-kb
[strain benchmark](benchmarks/strains-guard200-results.jsonl), which also
orders reads by genome, both 5% and 10% divergent mixtures were highly
fragmented. These small, independent-substitution simulations do not establish
strain reconstruction or the safety of normalization in metagenomes.

## Synthetic smoke benchmark

From the repository root, run:

```sh
python3 crates/viral-unitig/benchmarks/synthetic.py target/release/viral-unitig
python3 crates/viral-unitig/benchmarks/mixed.py target/release/viral-unitig
python3 crates/viral-unitig/benchmarks/strains.py target/release/viral-unitig --min-count 10
# Compare less aggressive settings:
python3 crates/viral-unitig/benchmarks/strains.py target/release/viral-unitig --min-count 3 --min-fraction 0.2
```

This simulates a single random 3 kb sequence, uniformly sampled 250-base
fragments, and independent uniform substitutions, insertions and deletions.
It reports nominal depth, unitig count and the longest unitig exactly matching
the reference on either strand, and the number of fragments written to the
uncertain file. It does not model nanopore error distributions, repeats,
bacterial reads, or large datasets.

In an earlier run using `--min-count 10`, in the 1000x case with 2%
substitutions and 2% indel events, the assembler
produced one 2,985-base exact primary unitig, while 1,398 short unitigs were
classified as uncertain. At 100x with the same errors, it produced five
primary unitigs, with a longest exact unitig of 2,876 bases.

The second benchmark simulates two 1 kb genomes using 200-base fragments and
1% substitutions plus 1% indel events. A 1000x genome mixed with an unrelated
50x genome in the earlier `--min-count 10` run produced three primary unitigs, with longest exact matches of
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

For these fixed synthetic seeds, `--min-count 10 --min-fraction 0.4` trades minority recovery for
longer major-strain unitigs compared with `--min-count 3 --min-fraction 0.2`:

| Divergence | Settings | Primary unitigs | Longest exact major unitig | Minor-specific 21-mer recovery |
|---|---|---:|---:|---:|
| 10% | less aggressive | 528 | 578 bp | 93.38% |
| 10% | aggressive | 162 | 9,990 bp | 62.07% |
| 5% | less aggressive | 674 | 267 bp | 86.70% |
| 5% | aggressive | 78 | 9,987 bp | 25.81% |

Neither setting produces a complete minor-strain contig in this simulation.
Real errors and population structure may behave very differently.
