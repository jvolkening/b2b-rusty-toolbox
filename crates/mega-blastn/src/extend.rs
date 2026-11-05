//! Ungapped and gapped x-drop extension.
//!
//! Gapped scoring is done in half-units (match +2, mismatch -4, gap -5 per
//! position), which is exactly megablast's greedy scoring (reward 1,
//! penalty -2, gap cost 2.5 per position) in integer arithmetic.

use crate::blastdb::Subject;

pub const MATCH: i32 = 2;
pub const MISMATCH: i32 = -4;
pub const GAP: i32 = -5;
const NEG: i32 = i32::MIN / 2;

#[inline(always)]
fn sc(a: u8, b: u8) -> i32 {
    if a == b && a < 4 { MATCH } else { MISMATCH }
}

/// Edit operation: 0 = aligned pair, 1 = subject-only (gap in query),
/// 2 = query-only (gap in subject).
pub type Op = (u8, u32);

#[derive(Clone, Debug)]
pub struct Hsp {
    pub qs: usize,
    pub qe: usize,
    pub ss: usize,
    pub se: usize,
    /// Raw (integer) BLAST score.
    pub score: i32,
    pub ops: Vec<Op>,
    pub length: u32,
    pub nident: u32,
    pub gaps: u32,
    pub gap_opens: u32,
}

impl Hsp {
    /// BLAST's HSP containment test (`s_HSPIsContained`): the start and end
    /// points of the input lie inside this HSP's rectangle and the starts or
    /// the ends are on nearby diagonals. `score` is the input's score.
    pub fn envelops(&self, qs: usize, qe: usize, ss: usize, se: usize, score: i32) -> bool {
        const MIN_DIAG_SEP: i64 = 6;
        let inside = |q: usize, s: usize| q >= self.qs && q <= self.qe && s >= self.ss && s <= self.se;
        let close = |q1: usize, s1: usize, q2: usize, s2: usize| {
            ((q1 as i64 - s1 as i64) - (q2 as i64 - s2 as i64)).abs() < MIN_DIAG_SEP
        };
        score <= self.score
            && inside(qs, ss)
            && inside(qe, se)
            && (close(self.qs, self.ss, qs, ss) || close(self.qe, self.se, qe, se))
    }

    /// True if the alignment pairs query position `qpos` with subject
    /// position `spos`.
    pub fn on_path(&self, qpos: usize, spos: usize) -> bool {
        if qpos < self.qs || qpos >= self.qe {
            return false;
        }
        let (mut qi, mut si) = (self.qs, self.ss);
        for &(op, n) in &self.ops {
            let n = n as usize;
            match op {
                0 => {
                    if qpos < qi + n {
                        return spos == si + (qpos - qi);
                    }
                    qi += n;
                    si += n;
                }
                1 => si += n,
                _ => {
                    if qpos < qi + n {
                        return false;
                    }
                    qi += n;
                }
            }
        }
        false
    }

    pub fn mismatches(&self) -> u32 {
        self.length - self.nident - self.gaps
    }

    pub fn compute_stats(&mut self, q: &[u8], s: &Subject) {
        let (mut qi, mut si) = (self.qs, self.ss);
        let (mut len, mut id, mut gaps, mut opens) = (0u32, 0u32, 0u32, 0u32);
        for &(op, n) in &self.ops {
            len += n;
            match op {
                0 => {
                    for _ in 0..n {
                        let (a, b) = (q[qi], s.get(si));
                        if a == b && a < 4 {
                            id += 1;
                        }
                        qi += 1;
                        si += 1;
                    }
                }
                1 => {
                    si += n as usize;
                    gaps += n;
                    opens += 1;
                }
                _ => {
                    qi += n as usize;
                    gaps += n;
                    opens += 1;
                }
            }
        }
        self.length = len;
        self.nident = id;
        self.gaps = gaps;
        self.gap_opens = opens;
    }
}

pub fn push_op(ops: &mut Vec<Op>, op: u8, n: u32) {
    if let Some(l) = ops.last_mut()
        && l.0 == op {
            l.1 += n;
            return;
        }
    ops.push((op, n));
}

/// Ungapped x-drop extension of an exact seed. `x` is in raw units
/// (match +1, mismatch -2). Returns (qs, qe, ss, se, score).
pub fn ungapped(
    q: &[u8],
    s: &Subject,
    qpos: usize,
    spos: usize,
    seed_len: usize,
    x: i32,
) -> (usize, usize, usize, usize, i32) {
    let mut score = seed_len as i32;
    // right
    let (mut best, mut cur) = (score, score);
    let mut bi = 0usize;
    let maxr = (q.len() - qpos - seed_len).min(s.len() - spos - seed_len);
    for i in 0..maxr {
        let (a, b) = (q[qpos + seed_len + i], s.get(spos + seed_len + i));
        cur += if a == b && a < 4 { 1 } else { -2 };
        if cur > best {
            best = cur;
            bi = i + 1;
        } else if best - cur > x {
            break;
        }
    }
    let right = bi;
    score = best;
    // left
    let (mut best, mut cur) = (score, score);
    let mut bl = 0usize;
    let maxl = qpos.min(spos);
    for i in 0..maxl {
        let (a, b) = (q[qpos - 1 - i], s.get(spos - 1 - i));
        cur += if a == b && a < 4 { 1 } else { -2 };
        if cur > best {
            best = cur;
            bl = i + 1;
        } else if best - cur > x {
            break;
        }
    }
    (qpos - bl, qpos + seed_len + right, spos - bl, spos + seed_len + right, best)
}

/// Center of the best-scoring window of `win` positions inside an ungapped
/// segment; returns offset from the segment start.
pub fn gapped_start(q: &[u8], s: &Subject, qs: usize, ss: usize, len: usize, win: usize) -> usize {
    if len <= win {
        return len / 2;
    }
    let m = |i: usize| -> i32 {
        let (a, b) = (q[qs + i], s.get(ss + i));
        if a == b && a < 4 { 1 } else { -2 }
    };
    let mut cur: i32 = (0..win).map(m).sum();
    let (mut best, mut bi) = (cur, 0usize);
    for i in 1..=(len - win) {
        cur += m(i + win - 1) - m(i - 1);
        if cur > best {
            best = cur;
            bi = i;
        }
    }
    bi + win / 2
}

/// Traceback storage at 2 bits per DP cell.
struct PackedOps {
    data: Vec<u8>,
    len: usize,
}

impl PackedOps {
    fn push(&mut self, v: u8) {
        if self.len % 4 == 0 {
            self.data.push(0);
        }
        let sh = 2 * (self.len % 4);
        let b = &mut self.data[self.len / 4];
        *b = (*b & !(3 << sh)) | (v << sh);
        self.len += 1;
    }

    fn truncate(&mut self, n: usize) {
        self.len = n;
        self.data.truncate(n.div_ceil(4));
    }

    #[inline]
    fn get(&self, i: usize) -> u8 {
        (self.data[i / 4] >> (2 * (i % 4))) & 3
    }
}

/// One-directional gapped x-drop extension. Returns (best half-score, query
/// consumed, subject consumed, ops from the far end back toward the origin).
pub fn xdrop_extend<FQ: Fn(usize) -> u8, FS: Fn(usize) -> u8>(
    fq: FQ,
    qlen: usize,
    fs: FS,
    slen: usize,
    x: i32,
) -> (i32, usize, usize, Vec<u8>) {
    let mut best = 0i32;
    let (mut bi, mut bj) = (0usize, 0usize);
    let mut ops = PackedOps { data: Vec::new(), len: 0 };
    let mut row_start: Vec<usize> = vec![0];
    let mut row_off: Vec<usize> = vec![0];
    let mut prev: Vec<i32> = vec![0];
    ops.push(3);
    let mut j = 1usize;
    while j <= slen && (GAP * j as i32) >= -x {
        prev.push(GAP * j as i32);
        ops.push(1);
        j += 1;
    }
    let mut plo = 0usize;
    let mut cur: Vec<i32> = Vec::new();
    for i in 1..=qlen {
        let phi = plo + prev.len();
        let qc = fq(i - 1);
        cur.clear();
        let off = ops.len;
        let thr = best - x;
        let mut j = plo;
        while j <= slen {
            let mut v = NEG;
            let mut op = 0u8;
            if j > plo && j - 1 < phi {
                v = prev[j - 1 - plo] + sc(qc, fs(j - 1));
            }
            if j < phi {
                let u = prev[j - plo] + GAP;
                if u > v {
                    v = u;
                    op = 2;
                }
            }
            if let Some(&l) = cur.last() {
                let l = l + GAP;
                if l > v {
                    v = l;
                    op = 1;
                }
            }
            if j >= phi && v < thr {
                break;
            }
            if v < thr {
                v = NEG;
            }
            cur.push(v);
            ops.push(op);
            j += 1;
        }
        row_start.push(plo);
        row_off.push(off);
        let mut first = usize::MAX;
        let mut last = 0usize;
        for (k, &v) in cur.iter().enumerate() {
            if v > NEG {
                if first == usize::MAX {
                    first = k;
                }
                last = k;
                if v >= best {
                    best = v;
                    bi = i;
                    bj = plo + k;
                }
            }
        }
        if first == usize::MAX {
            break;
        }
        ops.truncate(off + last + 1);
        cur.truncate(last + 1);
        plo += first;
        prev.clear();
        prev.extend_from_slice(&cur[first..]);
    }
    let mut out = Vec::new();
    let (mut i, mut j) = (bi, bj);
    while i > 0 || j > 0 {
        let op = ops.get(row_off[i] + (j - row_start[i]));
        match op {
            0 => {
                i -= 1;
                j -= 1;
            }
            2 => i -= 1,
            _ => j -= 1,
        }
        out.push(op);
    }
    (best, bi, bj, out)
}

/// Gapped extension around (qpos, spos). `x` is in half-units.
pub fn gapped(q: &[u8], s: &Subject, qpos: usize, spos: usize, x: i32) -> Hsp {
    let (rs, rq, rsj, rops) = xdrop_extend(
        |i| q[qpos + i],
        q.len() - qpos,
        |j| s.get(spos + j),
        s.len() - spos,
        x,
    );
    let (ls, lq, lsj, lops) = xdrop_extend(
        |i| q[qpos - 1 - i],
        qpos,
        |j| s.get(spos - 1 - j),
        spos,
        x,
    );
    let mut ops: Vec<Op> = Vec::new();
    for &o in &lops {
        push_op(&mut ops, o, 1);
    }
    for &o in rops.iter().rev() {
        push_op(&mut ops, o, 1);
    }
    Hsp {
        qs: qpos - lq,
        qe: qpos + rq,
        ss: spos - lsj,
        se: spos + rsj,
        score: (ls + rs).div_euclid(2),
        ops,
        length: 0,
        nident: 0,
        gaps: 0,
        gap_opens: 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extend_identical() {
        let q: Vec<u8> = (0..50).map(|i| ((i * 7 + i / 3) % 4) as u8).collect();
        let f = |i: usize| q[i];
        let (s, bi, bj, ops) = xdrop_extend(f, 50, f, 50, 108);
        assert_eq!((s, bi, bj), (100, 50, 50));
        assert!(ops.iter().all(|&o| o == 0));
    }

    #[test]
    fn extend_with_gap() {
        let a: Vec<u8> = (0..60).map(|i| (((i * 13) ^ (i >> 2)) % 4) as u8).collect();
        let mut b = a.clone();
        b.remove(30);
        let (s, bi, bj, _) = xdrop_extend(|i| a[i], 60, |j| b[j], 59, 108);
        assert_eq!((bi, bj), (60, 59));
        assert_eq!(s, 2 * 59 - 5);
    }
}
