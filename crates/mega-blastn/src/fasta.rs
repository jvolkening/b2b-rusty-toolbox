//! Query FASTA parsing, sequence encoding and DUST low-complexity masking.

use std::io::{self, BufRead};

/// Base codes: 0-3 = A,C,G,T; >=4 = ambiguity (4 + NCBI4na mask).
pub fn encode_base(c: u8) -> Option<u8> {
    let m = match c.to_ascii_uppercase() {
        b'A' => return Some(0),
        b'C' => return Some(1),
        b'G' => return Some(2),
        b'T' | b'U' => return Some(3),
        b'M' => 3,
        b'R' => 5,
        b'S' => 6,
        b'V' => 7,
        b'W' => 9,
        b'Y' => 10,
        b'H' => 11,
        b'K' => 12,
        b'D' => 13,
        b'B' => 14,
        b'N' => 15,
        b'X' => 15,
        _ => return None,
    };
    Some(4 + m)
}

pub fn decode_base(code: u8) -> u8 {
    match code {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        c => match c - 4 {
            3 => b'M',
            5 => b'R',
            6 => b'S',
            7 => b'V',
            9 => b'W',
            10 => b'Y',
            11 => b'H',
            12 => b'K',
            13 => b'D',
            14 => b'B',
            _ => b'N',
        },
    }
}

pub fn complement(code: u8) -> u8 {
    if code < 4 {
        3 - code
    } else {
        let n = code - 4;
        let c = ((n & 1) << 3) | ((n & 8) >> 3) | ((n & 2) << 1) | ((n & 4) >> 1);
        4 + c
    }
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&c| complement(c)).collect()
}

pub struct Query {
    pub id: String,
    pub title: String,
    pub fwd: Vec<u8>,
    pub rc: Vec<u8>,
    /// Half-open DUST-masked intervals on the forward strand.
    pub masked: Vec<(usize, usize)>,
}

impl Query {
    pub fn len(&self) -> usize {
        self.fwd.len()
    }

    pub fn is_empty(&self) -> bool {
        self.fwd.is_empty()
    }
}

pub struct QueryReader<R: BufRead> {
    inner: R,
    pending: Option<String>,
    dust: bool,
    line: String,
}

impl<R: BufRead> QueryReader<R> {
    pub fn new(inner: R, dust: bool) -> Self {
        QueryReader { inner, pending: None, dust, line: String::new() }
    }

    fn next_record(&mut self) -> io::Result<Option<Query>> {
        let header = loop {
            if let Some(h) = self.pending.take() {
                break h;
            }
            self.line.clear();
            if self.inner.read_line(&mut self.line)? == 0 {
                return Ok(None);
            }
            let l = self.line.trim_end();
            if let Some(h) = l.strip_prefix('>') {
                break h.to_string();
            }
            // text before the first header is ignored
        };
        let mut seq = Vec::new();
        loop {
            self.line.clear();
            if self.inner.read_line(&mut self.line)? == 0 {
                break;
            }
            let l = self.line.trim_end();
            if let Some(h) = l.strip_prefix('>') {
                self.pending = Some(h.to_string());
                break;
            }
            seq.extend(l.bytes().filter_map(encode_base));
        }
        let id = header.split_whitespace().next().unwrap_or("").to_string();
        let masked = if self.dust { dust_mask(&seq, 20, 64) } else { Vec::new() };
        let rc = reverse_complement(&seq);
        Ok(Some(Query { id, title: header, fwd: seq, rc, masked }))
    }

    /// Read queries until roughly `max_bases` residues have been collected.
    pub fn next_batch(&mut self, max_bases: usize) -> io::Result<Vec<Query>> {
        let mut out = Vec::new();
        let mut total = 0usize;
        while total < max_bases || out.is_empty() {
            match self.next_record()? {
                Some(q) => {
                    total += q.len();
                    out.push(q);
                }
                None => break,
            }
        }
        Ok(out)
    }
}

/// Symmetric DUST (Morgulis et al. 2006), ported from the reference
/// implementation in `sdust.c` by Heng Li (minimap2, MIT license). Returns
/// merged half-open masked intervals.
pub fn dust_mask(seq: &[u8], t: i32, w: i32) -> Vec<(usize, usize)> {
    const WLEN: i32 = 3;
    const WTOT: usize = 64;
    #[derive(Clone, Copy)]
    struct Perf {
        start: i32,
        finish: i32,
        r: i32,
        l: i32,
    }
    let mut p: Vec<Perf> = Vec::new(); // sorted by descending start
    let mut res: Vec<(i32, i32)> = Vec::new();
    let mut win: std::collections::VecDeque<usize> = std::collections::VecDeque::new();
    let (mut rv, mut rw, mut ll) = (0i32, 0i32, 0i32);
    let mut cv = [0i32; WTOT];
    let mut cw = [0i32; WTOT];

    let save = |res: &mut Vec<(i32, i32)>, p: &mut Vec<Perf>, start: i32| {
        if p.is_empty() || p[p.len() - 1].start >= start {
            return;
        }
        let last = p[p.len() - 1];
        let mut saved = false;
        if let Some(prev) = res.last_mut()
            && last.start <= prev.1 {
                saved = true;
                prev.1 = prev.1.max(last.finish);
            }
        if !saved {
            res.push((last.start, last.finish));
        }
        let mut i = p.len() as i32 - 1;
        while i >= 0 && p[i as usize].start < start {
            i -= 1;
        }
        p.truncate((i + 1) as usize);
    };

    let n = seq.len();
    let (mut l, mut tword) = (0i32, 0usize);
    for i in 0..=n {
        let b = if i < n && seq[i] < 4 { seq[i] as i32 } else { 4 };
        if b < 4 {
            l += 1;
            tword = ((tword << 2) | b as usize) & (WTOT - 1);
            if l >= WLEN {
                let start = (if l - w > 0 { l - w } else { 0 }) + (i as i32 + 1 - l);
                save(&mut res, &mut p, start);
                // shift_window
                if win.len() as i32 > w - WLEN {
                    let s = win.pop_front().unwrap();
                    cw[s] -= 1;
                    rw -= cw[s];
                    if ll > win.len() as i32 {
                        ll -= 1;
                        cv[s] -= 1;
                        rv -= cv[s];
                    }
                }
                win.push_back(tword);
                ll += 1;
                rw += cw[tword];
                cw[tword] += 1;
                rv += cv[tword];
                cv[tword] += 1;
                if cv[tword] * 10 > t << 1 {
                    loop {
                        let s = win[win.len() - ll as usize];
                        cv[s] -= 1;
                        rv -= cv[s];
                        ll -= 1;
                        if s == tword {
                            break;
                        }
                    }
                }
                if rw * 10 > ll * t {
                    // find_perfect
                    let mut c = cv;
                    let mut r = rv;
                    let (mut max_r, mut max_l) = (0i32, 0i32);
                    let mut i2 = win.len() as i32 - ll - 1;
                    while i2 >= 0 {
                        let tt = win[i2 as usize];
                        r += c[tt];
                        c[tt] += 1;
                        let new_r = r;
                        let new_l = win.len() as i32 - i2 - 1;
                        if new_r * 10 > t * new_l {
                            let mut j = 0usize;
                            while j < p.len() && p[j].start >= i2 + start {
                                if max_r == 0 || p[j].r * max_l > max_r * p[j].l {
                                    max_r = p[j].r;
                                    max_l = p[j].l;
                                }
                                j += 1;
                            }
                            if max_r == 0 || new_r * max_l >= max_r * new_l {
                                max_r = new_r;
                                max_l = new_l;
                                p.insert(
                                    j,
                                    Perf {
                                        start: i2 + start,
                                        finish: win.len() as i32 + (WLEN - 1) + start,
                                        r: new_r,
                                        l: new_l,
                                    },
                                );
                            }
                        }
                        i2 -= 1;
                    }
                }
            }
        } else {
            let mut start = (if l - w + 1 > 0 { l - w + 1 } else { 0 }) + (i as i32 + 1 - l);
            while !p.is_empty() {
                save(&mut res, &mut p, start);
                start += 1;
            }
            l = 0;
            tword = 0;
            // the window restarts after an ambiguous base
            win.clear();
            rv = 0;
            rw = 0;
            ll = 0;
            cv = [0; WTOT];
            cw = [0; WTOT];
        }
    }
    res.into_iter().map(|(a, b)| (a as usize, b as usize)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rc_roundtrip() {
        let s: Vec<u8> = b"ACGTNRY".iter().map(|&c| encode_base(c).unwrap()).collect();
        assert_eq!(reverse_complement(&reverse_complement(&s)), s);
    }

    #[test]
    fn dust_masks_homopolymer() {
        let mut s: Vec<u8> = Vec::new();
        let mut x = 12345u64;
        for _ in 0..100 {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            s.push(((x >> 33) & 3) as u8);
        }
        s.extend(std::iter::repeat_n(0u8, 80));
        for _ in 0..100 {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            s.push(((x >> 33) & 3) as u8);
        }
        let m = dust_mask(&s, 20, 64);
        assert!(!m.is_empty());
        assert!(m.iter().any(|&(a, b)| a >= 90 && b <= 190 && b - a >= 60));
    }
}
