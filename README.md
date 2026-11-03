# b2b-rusty-toolbox

This repo is expected to eventually contain a miscellaneous set of bioinformatics tools written in Rust. These are not general purpose tools -- there are plenty of toolkits (`seqtk`, `seqkit`, `bedtools`, etc) for this purpose. The utilities in this toolkit will perform specific, probably unusual, tasks that the author has a need for and which need to be fast/efficient.

More documentation may come at some point. For now, consider this a personal Rust playground.

## Tools

- `n-every`: insert an `N` after every K bases in FASTA input.
- `bg-window`: shrink a bedgraph using a per-sequence adaptive sliding window (min/max/median/mean/sum) so each sequence yields a user-specified range of entries. See `bg-window --help`.
