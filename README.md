# b2b-rusty-toolbox

This repo is expected to eventually contain a miscellaneous set of bioinformatics tools written in Rust. These are not general purpose tools -- there are plenty of toolkits (`seqtk`, `seqkit`, `bedtools`, etc) for this purpose. The utilities in this toolkit will perform specific, probably unusual, tasks that the author has a need for and which need to be fast/efficient.

More documentation may come at some point. For now, consider this a personal Rust playground.

## Tools

- `n-every`: see `crates/n-every`.
- `mega-blastn`: megablast-style search against BLAST v5 nucleotide databases, with blastn-compatible tabular output (`crates/mega-blastn/README.md`).
