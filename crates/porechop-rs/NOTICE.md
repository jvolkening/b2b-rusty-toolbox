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
