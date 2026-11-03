use std::collections::{HashMap, HashSet};
use std::io::{self, BufRead, Write};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Stat {
    Min,
    Max,
    Median,
    Mean,
    Sum,
}

impl Stat {
    pub fn parse(s: &str) -> Option<Stat> {
        match s {
            "min" => Some(Stat::Min),
            "max" => Some(Stat::Max),
            "median" => Some(Stat::Median),
            "mean" => Some(Stat::Mean),
            "sum" => Some(Stat::Sum),
            _ => None,
        }
    }
}

/// Half-open interval [start, end) with a depth value.
pub type Interval = (u64, u64, f64);

fn err(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

/// Returns (name, length) for each FASTA record, in file order.
pub fn read_fasta_lengths<R: BufRead>(reader: R) -> io::Result<Vec<(String, u64)>> {
    let mut out: Vec<(String, u64)> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(h) = line.strip_prefix('>') {
            let name = h.split_whitespace().next().unwrap_or("").to_string();
            out.push((name, 0));
        } else if let Some(last) = out.last_mut() {
            last.1 += line.trim_end().len() as u64;
        } else if !line.trim().is_empty() {
            return Err(err("FASTA sequence data before first header".into()));
        }
    }
    Ok(out)
}

/// Reads a bedgraph into per-sequence interval lists, sorted by start.
pub fn read_bedgraph<R: BufRead>(reader: R) -> io::Result<HashMap<String, Vec<Interval>>> {
    let mut map: HashMap<String, Vec<Interval>> = HashMap::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let l = line.trim();
        if l.is_empty() || l.starts_with('#') || l.starts_with("track") || l.starts_with("browser") {
            continue;
        }
        let f: Vec<&str> = l.split_whitespace().collect();
        let bad = || err(format!("Malformed bedgraph line {}: {}", i + 1, l));
        if f.len() < 4 {
            return Err(bad());
        }
        let s: u64 = f[1].parse().map_err(|_| bad())?;
        let e: u64 = f[2].parse().map_err(|_| bad())?;
        let v: f64 = f[3].parse().map_err(|_| bad())?;
        if e <= s || !v.is_finite() {
            return Err(bad());
        }
        map.entry(f[0].to_string()).or_default().push((s, e, v));
    }
    for v in map.values_mut() {
        v.sort_by_key(|x| (x.0, x.1));
        if v.windows(2).any(|w| w[0].1 > w[1].0) {
            return Err(err("Overlapping bedgraph intervals".into()));
        }
    }
    Ok(map)
}

/// Expands intervals to a gap-free tiling of [0, len) with uncovered bases as depth 0.
fn dense(intervals: &[Interval], len: u64) -> Vec<Interval> {
    let mut out = Vec::with_capacity(intervals.len() * 2 + 1);
    let mut pos = 0;
    for &(s, e, v) in intervals {
        let e = e.min(len);
        if s >= e {
            continue;
        }
        if s > pos {
            out.push((pos, s, 0.0));
        }
        out.push((s, e, v));
        pos = e;
    }
    if pos < len {
        out.push((pos, len, 0.0));
    }
    out
}

fn summarize(segs: &mut [(f64, u64)], stat: Stat) -> f64 {
    let total: u64 = segs.iter().map(|s| s.1).sum();
    match stat {
        Stat::Min => segs.iter().map(|s| s.0).fold(f64::INFINITY, f64::min),
        Stat::Max => segs.iter().map(|s| s.0).fold(f64::NEG_INFINITY, f64::max),
        Stat::Sum => segs.iter().map(|s| s.0 * s.1 as f64).sum(),
        Stat::Mean => segs.iter().map(|s| s.0 * s.1 as f64).sum::<f64>() / total as f64,
        Stat::Median => {
            segs.sort_by(|a, b| a.0.total_cmp(&b.0));
            let kth = |k: u64| {
                let mut acc = 0;
                for s in segs.iter() {
                    acc += s.1;
                    if k < acc {
                        return s.0;
                    }
                }
                segs.last().unwrap().0
            };
            if total % 2 == 1 {
                kth(total / 2)
            } else {
                (kth(total / 2 - 1) + kth(total / 2)) / 2.0
            }
        }
    }
}

/// Summarizes `dense` (tiling of [0, len)) in consecutive windows of size `w`,
/// merging adjacent windows with identical values.
fn window_pass(dense: &[Interval], len: u64, w: u64, stat: Stat) -> Vec<Interval> {
    let mut out: Vec<Interval> = Vec::new();
    let mut idx = 0;
    let mut start = 0;
    let mut segs: Vec<(f64, u64)> = Vec::new();
    while start < len {
        let end = (start + w).min(len);
        segs.clear();
        while idx < dense.len() && dense[idx].0 < end {
            let (s, e, v) = dense[idx];
            let (cs, ce) = (s.max(start), e.min(end));
            segs.push((v, ce - cs));
            if e <= end {
                idx += 1;
            } else {
                break;
            }
        }
        let val = summarize(&mut segs, stat);
        match out.last_mut() {
            Some(last) if last.2 == val => last.1 = end,
            _ => out.push((start, end, val)),
        }
        start = end;
    }
    out
}

/// Finds a window size so the entry count for this sequence falls in [min, max].
/// If impossible (e.g. sequence too short), returns the closest result found.
pub fn adaptive_windows(
    intervals: &[Interval],
    len: u64,
    min: usize,
    max: usize,
    stat: Stat,
) -> Vec<Interval> {
    if len == 0 {
        return Vec::new();
    }
    let d = dense(intervals, len);
    let target = ((min + max) / 2).max(1) as u64;
    let mut w = len.div_ceil(target).max(1);
    let mut seen = HashSet::new();
    let mut best: Option<(usize, Vec<Interval>)> = None;
    for _ in 0..64 {
        if seen.contains(&w) {
            let Some(next) = (1..=len).find(|candidate| !seen.contains(candidate)) else {
                break;
            };
            w = next;
        }
        if !seen.insert(w) {
            break;
        }
        let res = window_pass(&d, len, w, stat);
        let n = res.len();
        let dist = if n < min { min - n } else { n.saturating_sub(max) };
        if dist == 0 {
            return res;
        }
        if best.as_ref().is_none_or(|b| dist < b.0) {
            best = Some((dist, res));
        }
        let scaled = (w as f64 * n as f64 / target as f64) as u64;
        w = if n > max { scaled.max(w + 1) } else { scaled.min(w.saturating_sub(1)) };
        if w == 0 || w > len {
            break;
        }
    }
    best.map(|b| b.1).unwrap_or_default()
}

pub fn run<W: Write>(
    fasta: &[(String, u64)],
    bedgraph: &HashMap<String, Vec<Interval>>,
    min: usize,
    max: usize,
    stat: Stat,
    mut out: W,
) -> io::Result<()> {
    let known: HashSet<&str> = fasta.iter().map(|f| f.0.as_str()).collect();
    if let Some(k) = bedgraph.keys().find(|k| !known.contains(k.as_str())) {
        return Err(err(format!("Bedgraph sequence '{k}' not found in FASTA")));
    }
    for (name, len) in fasta {
        let iv = bedgraph.get(name).map(Vec::as_slice).unwrap_or(&[]);
        if let Some(last) = iv.last() {
            if last.1 > *len {
                return Err(err(format!("Bedgraph interval exceeds length of '{name}'")));
            }
        }
        for (s, e, v) in adaptive_windows(iv, *len, min, max, stat) {
            writeln!(out, "{name}\t{s}\t{e}\t{v}")?;
        }
    }
    Ok(())
}
