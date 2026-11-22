use minimap2::{Aligner, Mapping, Strand, ffi};
use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::io::{BufRead, Write};

pub type Result<T> = std::result::Result<T, Box<dyn Error + Send + Sync>>;

#[derive(Default)]
pub struct Counts {
    pub n_total: u64,
    pub n_kept: u64,
}

#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub min_length: usize,
    pub min_identity: f64,
    pub chimeric: bool,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            min_length: 500,
            min_identity: 0.80,
            chimeric: false,
        }
    }
}

#[derive(Default)]
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

    fn passes(&self, options: Options) -> bool {
        self.target_bases.min(self.query_bases) >= options.min_length
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

fn passing_alignment(seq: &[u8], mapping: &Mapping, options: Options) -> Result<bool> {
    if mapping.strand != Strand::Reverse {
        return Ok(false);
    }
    let cigar = mapping
        .alignment
        .as_ref()
        .and_then(|a| a.cigar.as_ref())
        .ok_or("minimap2 returned an alignment without a CIGAR")?;
    let mut t = mapping.target_start as usize;
    let mut q = mapping.query_end as usize;
    let mut arms = [Arm::default(), Arm::default()];

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
            arms[side].add(target, query, matches);
        }
    }
    Ok(arms.iter().any(|arm| arm.passes(options)))
}

pub fn is_chimeric(seq: &[u8], options: Options) -> Result<bool> {
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
        || !options.min_identity.is_finite()
        || !(0.0..=1.0).contains(&options.min_identity)
    {
        return Err("invalid alignment thresholds".into());
    }
    if seq.len() / 2 < options.min_length || !has_inverted_seed(seq) {
        return Ok(false);
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
    for hit in &hits {
        if passing_alignment(seq, hit, options)? {
            return Ok(true);
        }
    }
    Ok(false)
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
) -> Result<()> {
    let results: Vec<_> = pool.install(|| {
        batch
            .par_iter()
            .map(|record| is_chimeric(&record.sequence, options))
            .collect()
    });
    for (record, result) in batch.iter().zip(results) {
        let chimeric =
            result.map_err(|e| format!("record {}: {e}", String::from_utf8_lossy(&record.id)))?;
        if chimeric == options.chimeric {
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
) -> Result<()> {
    if reader.fill_buf()?.is_empty() {
        return Ok(());
    }
    let mut records = needletail::parse_fastx_reader(reader)?;
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
            write_batch(&mut batch, writer, options, pool, counts)?;
            bytes = 0;
        }
    }
    write_batch(&mut batch, writer, options, pool, counts)?;
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
            min_identity: 1.0,
            chimeric: true,
        };
        assert!(passing_alignment(&seq, &mapping, options).unwrap());
        options.min_length = 601;
        assert!(!passing_alignment(&seq, &mapping, options).unwrap());
        options.min_length = 600;
        seq[0] = b'C';
        assert!(!passing_alignment(&seq, &mapping, options).unwrap());
        options.min_identity = 0.99;
        assert!(passing_alignment(&seq, &mapping, options).unwrap());
    }

    #[test]
    fn inverted_seed_requires_distinct_positions_and_handles_ambiguous_bases() {
        assert!(!has_inverted_seed(b"ACGTACGTACGTACG"));
        assert!(!has_inverted_seed(b"AAAAAAAAAAAAAAANNN"));
        assert!(has_inverted_seed(b"AAAAAAAAAAAAAAANNTTTTTTTTTTTTTTT"));
        assert!(has_inverted_seed(b"aaaaaaaaaaaaaaannTTTTTTTTTTTTTTT"));
    }
}
