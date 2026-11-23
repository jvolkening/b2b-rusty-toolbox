# Attribution and licenses

`porechop-rs` is a modified Rust implementation of
[Porechop](https://github.com/rrwick/Porechop), Copyright 2017 Ryan Wick,
distributed under GPL-3.0-or-later; see [LICENSE](LICENSE). The adapter
sequences, alignment scoring implementation and test fixtures are derived from
Porechop. The Rust implementation was moved from the `jvolkening/Porechop`
branch at commit `5250b96480735cef5316767c830eb5a090f7385c` and adapted for
this workspace.

The vendored SeqAn 2.3.1 headers are distributed under the 3-clause BSD license
in [LICENSE-SEQAN](LICENSE-SEQAN). They include portions covered by
[LICENSE-NVIDIA](LICENSE-NVIDIA), Boost's [LICENSE-BOOST](LICENSE-BOOST),
and samtools' [LICENSE-SAMTOOLS](LICENSE-SAMTOOLS). Their original notices
remain in the headers.

The vendored [Parasail](https://github.com/jeffdaily/parasail) 2.6.1 source
is Copyright 2015 Battelle Memorial Institute, distributed under the license
in [LICENSE-PARASAIL](LICENSE-PARASAIL).
Its traceback allocation handling has been modified to report allocation
failures rather than dereferencing null pointers.

Daily, J. (2016). Parasail: SIMD C library for global, semi-global, and local
pairwise sequence alignments. *BMC Bioinformatics*, 17, 81.
https://doi.org/10.1186/s12859-016-0930-z
