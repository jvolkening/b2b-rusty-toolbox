use flate2::read::MultiGzDecoder;
use minimap2::{Aligner, Mapping, Strand, ffi};
use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::io::{BufRead, BufReader, Cursor, Read, Write};

pub type Result<T> = std::result::Result<T, Box<dyn Error + Send + Sync>>;

#[derive(Default)]
pub struct Counts {
    pub n_total: u64,
    pub n_kept: u64,
}

#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub min_length: usize,
    pub min_fraction: Option<f64>,
    pub min_identity: f64,
    pub chimeric: bool,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            min_length: 500,
            min_fraction: None,
            min_identity: 0.80,
            chimeric: false,
        }
    }
}

#[derive(Default, Debug)]
struct Arm {
    matches: usize,
    columns: usize,
    target_start: usize,
    target_end: usize,
    query_start: usize,
    query_end: usize,
    target_bases: usize,
    query_bases: usize,
}

impl Arm {
    fn add(&mut self, target: Option<usize>, query: Option<usize>, matches: bool) {
        self.columns += 1;
        self.matches += usize::from(matches);
        if let Some(t) = target {
            if self.target_bases == 0 {
                self.target_start = t;
                self.target_end = t + 1;
            } else {
                self.target_start = self.target_start.min(t);
                self.target_end = self.target_end.max(t + 1);
            }
            self.target_bases += 1;
        }
        if let Some(q) = query {
            if self.query_bases == 0 {
                self.query_start = q;
                self.query_end = q + 1;
            } else {
                self.query_start = self.query_start.min(q);
                self.query_end = self.query_end.max(q + 1);
            }
            self.query_bases += 1;
        }
    }

    fn passes(&self, options: Options, read_length: usize) -> bool {
        let length_passes = self.target_bases.min(self.query_bases) >= options.min_length
            || options.min_fraction.is_some_and(|fraction| {
                (self.target_bases + self.query_bases) as f64 >= fraction * read_length as f64
            });
        length_passes
            && self.target_bases > 0
            && self.query_bases > 0
            && (self.target_end <= self.query_start || self.query_end <= self.target_start)
            && self.matches as f64 >= options.min_identity * self.columns as f64
    }
}

fn complement(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(b'T'),
        b'C' => Some(b'G'),
        b'G' => Some(b'C'),
        b'T' => Some(b'A'),
        _ => None,
    }
}

fn has_inverted_seed(seq: &[u8]) -> bool {
    const K: usize = 15;
    const MASK: u64 = (1 << (2 * K)) - 1;
    let mut seen = HashMap::with_capacity(seq.len());
    let (mut forward, mut reverse, mut valid) = (0u64, 0u64, 0usize);
    for (i, base) in seq.iter().enumerate() {
        let code = match base.to_ascii_uppercase() {
            b'A' => 0u64,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => {
                valid = 0;
                forward = 0;
                reverse = 0;
                continue;
            }
        };
        forward = ((forward << 2) | code) & MASK;
        reverse = (reverse >> 2) | ((3 - code) << (2 * (K - 1)));
        valid += 1;
        if valid >= K {
            let start = i + 1 - K;
            if seen
                .get(&reverse)
                .is_some_and(|&earliest| start >= earliest + K)
            {
                return true;
            }
            seen.entry(forward).or_insert(start);
        }
    }
    false
}

fn alignment_columns(
    seq: &[u8],
    mapping: &Mapping,
    mut visit: impl FnMut(usize, Option<usize>, Option<usize>, bool),
) -> Result<()> {
    let cigar = mapping
        .alignment
        .as_ref()
        .and_then(|a| a.cigar.as_ref())
        .ok_or("minimap2 returned an alignment without a CIGAR")?;
    let mut t = mapping.target_start as usize;
    let mut q = mapping.query_end as usize;
    // On a reverse alignment the query runs backwards. Split at the point
    // where its coordinates cross the target to exclude overlapping arms.
    for &(len, op) in cigar {
        for _ in 0..len {
            let (target, query) = match op {
                0 | 7 | 8 => {
                    q -= 1;
                    let pair = (Some(t), Some(q));
                    t += 1;
                    pair
                }
                1 => {
                    q -= 1;
                    (None, Some(q))
                }
                2 => {
                    let pos = t;
                    t += 1;
                    (Some(pos), None)
                }
                _ => return Err(format!("unexpected CIGAR operation {op}").into()),
            };
            let tp = target.unwrap_or(t);
            let qp = query.unwrap_or(q.saturating_sub(1));
            let side = match tp.cmp(&qp) {
                std::cmp::Ordering::Less => 0,
                std::cmp::Ordering::Greater => 1,
                std::cmp::Ordering::Equal => continue,
            };
            let matches = match (target, query) {
                (Some(t), Some(q)) => complement(seq[q]) == Some(seq[t].to_ascii_uppercase()),
                _ => false,
            };
            visit(side, target, query, matches);
        }
    }
    Ok(())
}

fn passing_alignment(
    seq: &[u8],
    mapping: &Mapping,
    options: Options,
) -> Result<Option<(Arm, usize)>> {
    if mapping.strand != Strand::Reverse {
        return Ok(None);
    }
    let mut arms = [Arm::default(), Arm::default()];
    alignment_columns(seq, mapping, |side, target, query, matches| {
        arms[side].add(target, query, matches);
    })?;
    Ok(arms
        .into_iter()
        .enumerate()
        .find(|(_, arm)| arm.passes(options, seq.len()))
        .map(|(side, arm)| (arm, side)))
}

struct Hit {
    mapping: Option<Mapping>,
    arm: Arm,
    side: usize,
}

fn find_chimera(seq: &[u8], options: Options, keep_alignment: bool) -> Result<Option<Hit>> {
    if seq.len() > i32::MAX as usize {
        return Err("sequence exceeds minimap2's length limit".into());
    }
    if !seq.iter().all(|b| {
        matches!(
            b.to_ascii_uppercase(),
            b'A' | b'C'
                | b'G'
                | b'T'
                | b'N'
                | b'R'
                | b'Y'
                | b'S'
                | b'W'
                | b'K'
                | b'M'
                | b'B'
                | b'D'
                | b'H'
                | b'V'
        )
    }) {
        return Err("sequence contains a non-IUPAC DNA character".into());
    }
    if options.min_length == 0
        || options
            .min_fraction
            .is_some_and(|f| !f.is_finite() || f <= 0.0 || f > 1.0)
        || !options.min_identity.is_finite()
        || !(0.0..=1.0).contains(&options.min_identity)
    {
        return Err("invalid alignment thresholds".into());
    }
    if (options.min_fraction.is_none() && seq.len() / 2 < options.min_length)
        || !has_inverted_seed(seq)
    {
        return Ok(None);
    }
    let mut builder = Aligner::builder().ava_ont().with_cigar();
    builder.mapopt.flag |= ffi::MM_F_REV_ONLY as i64;
    builder.idxopt.bucket_bits = 10;
    let aligner = builder
        .with_seq(seq)
        .map_err(|e| format!("indexing: {e}"))?;
    let hits = aligner
        .map(seq, false, false, None, None, None)
        .map_err(|e| format!("mapping: {e}"))?;
    for mapping in hits {
        if let Some((arm, side)) = passing_alignment(seq, &mapping, options)? {
            return Ok(Some(Hit {
                mapping: keep_alignment.then_some(mapping),
                arm,
                side,
            }));
        }
    }
    Ok(None)
}

pub fn is_chimeric(seq: &[u8], options: Options) -> Result<bool> {
    Ok(find_chimera(seq, options, false)?.is_some())
}

fn display_complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        base => base,
    }
}

fn write_alignment(writer: &mut dyn Write, id: &[u8], seq: &[u8], hit: &Hit) -> Result<()> {
    let arm = &hit.arm;
    writer.write_all(b"Read: ")?;
    writer.write_all(id)?;
    writeln!(
        writer,
        "\nLength: {}\nTarget (+): {}-{}\nQuery (-): {}-{}\nIdentity: {}/{} ({:.2}%)\n",
        seq.len(),
        arm.target_start + 1,
        arm.target_end,
        arm.query_end,
        arm.query_start + 1,
        arm.matches,
        arm.columns,
        100.0 * arm.matches as f64 / arm.columns as f64
    )?;
    let mut columns = Vec::with_capacity(arm.columns);
    let mapping = hit.mapping.as_ref().ok_or("missing alignment for report")?;
    alignment_columns(seq, mapping, |side, target, query, matches| {
        if side == hit.side {
            columns.push((target, query, matches));
        }
    })?;
    let width = seq.len().to_string().len();
    for chunk in columns.chunks(60) {
        let mut target = String::with_capacity(chunk.len());
        let mut query = String::with_capacity(chunk.len());
        let mut markers = String::with_capacity(chunk.len());
        let (mut t_first, mut t_last, mut q_first, mut q_last) = (None, None, None, None);
        for &(t, q, matches) in chunk {
            target.push(t.map_or('-', |pos| seq[pos].to_ascii_uppercase() as char));
            query.push(q.map_or('-', |pos| display_complement(seq[pos]) as char));
            markers.push(if matches { '|' } else { ' ' });
            if let Some(pos) = t {
                t_first.get_or_insert(pos + 1);
                t_last = Some(pos + 1);
            }
            if let Some(pos) = q {
                q_first.get_or_insert(pos + 1);
                q_last = Some(pos + 1);
            }
        }
        let coord = |value: Option<usize>| value.map_or_else(|| "-".to_owned(), |n| n.to_string());
        writeln!(
            writer,
            "Target {:>width$}  {}  {}\n       {:>width$}  {}\nQuery  {:>width$}  {}  {}\n",
            coord(t_first),
            target,
            coord(t_last),
            "",
            markers,
            coord(q_first),
            query,
            coord(q_last),
        )?;
    }
    Ok(())
}

struct Record {
    sequence: Vec<u8>,
    raw: Vec<u8>,
    id: Vec<u8>,
}

fn write_batch<W: Write>(
    batch: &mut Vec<Record>,
    writer: &mut W,
    options: Options,
    pool: &rayon::ThreadPool,
    counts: &mut Counts,
    report: &mut Option<&mut dyn Write>,
) -> Result<()> {
    let keep_alignment = report.is_some();
    let results: Vec<_> = pool.install(|| {
        batch
            .par_iter()
            .map(|record| find_chimera(&record.sequence, options, keep_alignment))
            .collect()
    });
    for (record, result) in batch.iter().zip(results) {
        let hit =
            result.map_err(|e| format!("record {}: {e}", String::from_utf8_lossy(&record.id)))?;
        if let (Some(hit), Some(report)) = (&hit, report.as_deref_mut()) {
            write_alignment(report, &record.id, &record.sequence, hit)?;
        }
        if hit.is_some() == options.chimeric {
            writer.write_all(&record.raw)?;
            writer.write_all(b"\n")?;
            counts.n_kept += 1;
        }
    }
    batch.clear();
    Ok(())
}

pub fn filter<R: BufRead + Send, W: Write>(
    mut reader: R,
    writer: &mut W,
    options: Options,
    pool: &rayon::ThreadPool,
    format: &mut Option<needletail::parser::Format>,
    counts: &mut Counts,
    report: &mut Option<&mut dyn Write>,
) -> Result<()> {
    if reader.fill_buf()?.is_empty() {
        return Ok(());
    }
    let mut prefix = [0; 2];
    reader.read_exact(&mut prefix)?;
    let input = Cursor::new(prefix).chain(reader);
    let mut input: Box<dyn BufRead + Send + '_> = if prefix == [0x1f, 0x8b] {
        Box::new(BufReader::new(MultiGzDecoder::new(input)))
    } else {
        Box::new(BufReader::new(input))
    };
    if input.fill_buf()?.is_empty() {
        return Ok(());
    }
    let mut records = needletail::parse_fastx_reader(input)?;
    let batch_size = pool.current_num_threads().saturating_mul(8).clamp(1, 256);
    let mut batch = Vec::with_capacity(batch_size);
    let mut bytes = 0usize;
    while let Some(record) = records.next() {
        let record = record?;
        if format.is_some_and(|previous| previous != record.format()) {
            return Err("cannot mix FASTA and FASTQ input files".into());
        }
        *format = Some(record.format());
        counts.n_total += 1;
        let seq = record.seq().into_owned();
        let raw = record.all().to_vec();
        bytes = bytes.saturating_add(seq.len()).saturating_add(raw.len());
        batch.push(Record {
            sequence: seq,
            raw,
            id: record.id().to_vec(),
        });
        if batch.len() >= batch_size || bytes >= 16 * 1024 * 1024 {
            write_batch(&mut batch, writer, options, pool, counts, report)?;
            bytes = 0;
        }
    }
    write_batch(&mut batch, writer, options, pool, counts, report)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use minimap2::Alignment;

    #[test]
    fn splits_full_hairpin_at_diagonal_and_scores_length_and_identity() {
        let mut seq = vec![b'A'; 600];
        seq.extend(vec![b'T'; 600]);
        let mapping = Mapping {
            strand: Strand::Reverse,
            query_end: 1200,
            alignment: Some(Alignment {
                nm: 0,
                cigar: Some(vec![(1200, 0)]),
                cigar_str: None,
                md: None,
                cs: None,
                alignment_score: None,
            }),
            ..Mapping::default()
        };
        let mut options = Options {
            min_length: 600,
            min_fraction: None,
            min_identity: 1.0,
            chimeric: true,
        };
        assert!(
            passing_alignment(&seq, &mapping, options)
                .unwrap()
                .is_some()
        );
        options.min_length = 601;
        assert!(
            passing_alignment(&seq, &mapping, options)
                .unwrap()
                .is_none()
        );
        options.min_length = 600;
        seq[0] = b'C';
        assert!(
            passing_alignment(&seq, &mapping, options)
                .unwrap()
                .is_none()
        );
        options.min_identity = 0.99;
        assert!(
            passing_alignment(&seq, &mapping, options)
                .unwrap()
                .is_some()
        );
    }

    #[test]
    fn fractional_threshold_uses_combined_arms_and_either_length_rule() {
        let arm = Arm {
            matches: 490,
            columns: 500,
            target_start: 0,
            target_end: 500,
            query_start: 900,
            query_end: 1390,
            target_bases: 500,
            query_bases: 490,
        };
        let mut options = Options {
            min_length: 500,
            min_fraction: Some(0.495),
            min_identity: 0.98,
            chimeric: false,
        };
        assert!(arm.passes(options, 2000));
        options.min_fraction = Some(0.4951);
        assert!(!arm.passes(options, 2000));
        options.min_length = 490;
        options.min_fraction = Some(0.9);
        assert!(arm.passes(options, 2000));
        options.min_fraction = None;
        options.min_length = 500;
        assert!(!arm.passes(options, 2000));
        options.min_fraction = Some(0.495);
        options.min_identity = 0.99;
        assert!(!arm.passes(options, 2000));
    }

    #[test]
    fn inverted_seed_requires_distinct_positions_and_handles_ambiguous_bases() {
        assert!(!has_inverted_seed(b"ACGTACGTACGTACG"));
        assert!(!has_inverted_seed(b"AAAAAAAAAAAAAAANNN"));
        assert!(has_inverted_seed(b"AAAAAAAAAAAAAAANNTTTTTTTTTTTTTTT"));
        assert!(has_inverted_seed(b"aaaaaaaaaaaaaaannTTTTTTTTTTTTTTT"));
    }

    #[test]
    fn inverted_seed_matches_brute_force_search() {
        fn brute(seq: &[u8]) -> bool {
            const K: usize = 15;
            seq.windows(K).enumerate().any(|(i, a)| {
                seq.windows(K).skip(i + K).any(|b| {
                    a.iter()
                        .zip(b.iter().rev())
                        .all(|(&x, &y)| complement(x).is_some_and(|c| c == y.to_ascii_uppercase()))
                })
            })
        }

        let mut state = 0x12345678u64;
        for case in 0..300 {
            let mut seq = Vec::new();
            for _ in 0..case % 120 {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                seq.push(b"ACGTNacgt"[(state as usize) % 9]);
            }
            if case % 3 == 0 {
                let seed = b"ACGTTAGCCGATACG";
                seq.extend(seed);
                seq.extend(seed.iter().rev().map(|&b| complement(b).unwrap()));
            }
            assert_eq!(has_inverted_seed(&seq), brute(&seq), "case {case}");
        }
    }

    #[test]
    fn identity_includes_insertions_deletions_and_ambiguous_bases() {
        let mut seq = b"AAAACCCCGGNNNNNNNNNNCCGGGCTTTT".to_vec();
        let mapping = Mapping {
            strand: Strand::Reverse,
            query_end: seq.len() as i32,
            alignment: Some(Alignment {
                nm: 2,
                cigar: Some(vec![(4, 0), (1, 1), (3, 0), (1, 2), (2, 0)]),
                cigar_str: None,
                md: None,
                cs: None,
                alignment_score: None,
            }),
            ..Mapping::default()
        };
        let mut options = Options {
            min_length: 10,
            min_fraction: None,
            min_identity: 9.0 / 11.0,
            chimeric: true,
        };
        assert!(
            passing_alignment(&seq, &mapping, options)
                .unwrap()
                .is_some()
        );
        options.min_identity = 9.01 / 11.0;
        assert!(
            passing_alignment(&seq, &mapping, options)
                .unwrap()
                .is_none()
        );
        options.min_identity = 9.0 / 11.0;
        seq[0] = b'N';
        assert!(
            passing_alignment(&seq, &mapping, options)
                .unwrap()
                .is_none()
        );
        options.min_identity = 8.0 / 11.0;
        assert!(
            passing_alignment(&seq, &mapping, options)
                .unwrap()
                .is_some()
        );
        options.min_length = 11;
        assert!(
            passing_alignment(&seq, &mapping, options)
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn report_shows_gaps_mismatches_and_reverse_coordinates() {
        let seq = b"AAAACCCCGGNNNNNNNNNNCCGGGCTTTT";
        let mapping = Mapping {
            strand: Strand::Reverse,
            query_end: seq.len() as i32,
            alignment: Some(Alignment {
                nm: 2,
                cigar: Some(vec![(4, 0), (1, 1), (3, 0), (1, 2), (2, 0)]),
                cigar_str: None,
                md: None,
                cs: None,
                alignment_score: None,
            }),
            ..Mapping::default()
        };
        let options = Options {
            min_length: 10,
            min_identity: 0.8,
            ..Options::default()
        };
        let (arm, side) = passing_alignment(seq, &mapping, options).unwrap().unwrap();
        let mut out = Vec::new();
        write_alignment(
            &mut out,
            b"example",
            seq,
            &Hit {
                mapping: Some(mapping),
                arm,
                side,
            },
        )
        .unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "Read: example\nLength: 30\nTarget (+): 1-10\nQuery (-): 30-21\nIdentity: 9/11 (81.82%)\n\n\
             Target  1  AAAA-CCCCGG  10\n           |||| ||| ||\nQuery  30  AAAAGCCC-GG  21\n\n"
        );
    }

    #[test]
    fn report_splits_diagonal_and_wraps_at_sixty_columns() {
        let mut seq = vec![b'A'; 61];
        seq.extend(vec![b'T'; 61]);
        let mapping = Mapping {
            strand: Strand::Reverse,
            query_end: 122,
            alignment: Some(Alignment {
                nm: 0,
                cigar: Some(vec![(122, 0)]),
                cigar_str: None,
                md: None,
                cs: None,
                alignment_score: None,
            }),
            ..Mapping::default()
        };
        let options = Options {
            min_length: 61,
            ..Options::default()
        };
        let (arm, side) = passing_alignment(&seq, &mapping, options).unwrap().unwrap();
        let mut out = Vec::new();
        write_alignment(
            &mut out,
            b"hairpin",
            &seq,
            &Hit {
                mapping: Some(mapping),
                arm,
                side,
            },
        )
        .unwrap();
        let report = String::from_utf8(out).unwrap();
        assert!(report.contains("Target (+): 1-61\nQuery (-): 122-62\nIdentity: 61/61 (100.00%)"));
        assert!(report.contains(&format!("Target   1  {}  60", "A".repeat(60))));
        assert!(report.contains(&format!("Query  122  {}  63", "A".repeat(60))));
        assert!(report.contains("Target  61  A  61\n            |\nQuery   62  A  62\n"));
    }

    #[test]
    fn report_displays_iupac_reverse_complements() {
        assert_eq!(
            b"ACGTRYSWKMBDHVN"
                .iter()
                .map(|&b| display_complement(b))
                .collect::<Vec<_>>(),
            b"TGCAYRSWMKVHDBN"
        );
    }

    #[test]
    fn report_handles_the_other_diagonal_side_and_write_errors() {
        let mut seq = vec![b'A'; 61];
        seq.extend(vec![b'T'; 61]);
        let mapping = Mapping {
            strand: Strand::Reverse,
            query_end: 122,
            alignment: Some(Alignment {
                nm: 0,
                cigar: Some(vec![(122, 0)]),
                cigar_str: None,
                md: None,
                cs: None,
                alignment_score: None,
            }),
            ..Mapping::default()
        };
        let mut arm = Arm::default();
        alignment_columns(&seq, &mapping, |side, target, query, matches| {
            if side == 1 {
                arm.add(target, query, matches);
            }
        })
        .unwrap();
        let hit = Hit {
            mapping: Some(mapping),
            arm,
            side: 1,
        };
        let mut out = Vec::new();
        write_alignment(&mut out, b"example", &seq, &hit).unwrap();
        let report = String::from_utf8(out).unwrap();
        assert!(report.contains("Target (+): 62-122\nQuery (-): 61-1\nIdentity: 61/61 (100.00%)"));
        assert!(report.contains("Target  62  "));
        assert!(report.contains("Query   61  "));

        struct BrokenWriter;
        impl Write for BrokenWriter {
            fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
                Err(std::io::Error::other("report write failed"))
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }
        let error = write_alignment(&mut BrokenWriter, b"example", &seq, &hit).unwrap_err();
        assert!(error.to_string().contains("report write failed"));
    }
}
