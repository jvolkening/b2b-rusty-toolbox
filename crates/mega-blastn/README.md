# mega-blastn

A megablast-style nucleotide search (`blastn -task megablast`) that reads existing
BLAST **v4/v5** nucleotide databases and writes `-outfmt 6` / `-outfmt 7` tabular output
in the same format as NCBI `blastn`.

```
mega-blastn -db nt -query queries.fa -outfmt 6 -num_threads 16 > hits.tsv
mega-blastn -db mydb -query q.fa -outfmt "7 std qlen slen qcovs sstrand btop"
```

Run `mega-blastn -h` for all options.

## How it works

* The database volumes (`.nin`/`.nsq`/`.nhr`, plus `.nal` aliases for multi-volume
  databases) are memory-mapped. The packed 2-bit sequence data is streamed once per
  query batch, so memory use does not grow with database size (resident pages are
  just page cache) and no database-sized structure is ever built.
* Query k-mers (both strands, DUST-masked like blastn's default) go into a small
  cache-friendly prefilter + bloom filter + hash table. Subjects are scanned at byte
  boundaries every `stride` bases (16, 12, 8 or 4, chosen from the query size so the
  k-mer table stays sparse) using k = word_size - stride + 1. Any exact run of at least
  `-word_size` (default 28) bases contains an aligned k-mer, so no seed that blastn
  would find is missed. Benchmarks (2 cores): 500 MB DB 0.26 s vs 0.69 s for blastn;
  3 Mb x 3 Mb alignment 3 s vs 168 s (more RAM: ~300 MB vs 137 MB).
* Seeds are extended like megablast: ungapped x-drop, then gapped x-drop (X = 100 bits)
  with traceback. Scoring is blastn's megablast default (+1/-2, linear gap cost 2.5),
  and E-values, bit scores, effective lengths, HSP containment tests and subject
  ordering follow the BLAST+ algorithm (Karlin-Altschul parameters 1.28/0.46).
* The DUST implementation is a port of the symmetric DUST reference code in `sdust.c`
  (Heng Li, MIT license).

## Validation

Against BLAST+ 2.12.0 on synthetic databases: identical hit sets, scores, bit scores,
E-values, subject/HSP ordering and formatting for default settings; HSP end points and
gap placement can differ in cases where equal-scoring alternatives exist (greedy vs.
dynamic-programming extension). With smaller `-word_size` values it reports a superset
of blastn's hits.

## Limitations

* Only megablast scoring (reward 1, penalty -2) and `-outfmt 6`/`7` are supported.
* Database alias restrictions (GI/seqid/taxid lists, OID masks) are rejected, not applied.
* A single (very large) subject sequence is processed by one thread.
* `sseqid`/`sacc` follow blastn for parsed and non-parsed databases, but exotic id
  types are only lightly covered.
