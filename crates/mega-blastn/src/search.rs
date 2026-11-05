//! The database scan: seed detection, extension and hit bookkeeping.

use crate::blastdb::{Database, Subject};
use crate::extend::{self, Hsp};
use crate::fasta::Query;
use crate::index::Index;
use crate::stats::SearchSpace;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};
use std::sync::Mutex;

/// Raw-unit x-drop for the ungapped extension (megablast default: 11).
const UNGAPPED_XDROP: i32 = 11;
/// Gapped x-drop in half-units: 100 bits -> 54 raw -> 108.
const GAPPED_XDROP: i32 = 108;
const START_WINDOW: usize = 11;
const CHUNK: u64 = 16;

pub struct Params {
    pub word_size: usize,
    pub evalue: f64,
    pub max_target_seqs: usize,
    pub max_hsps: Option<usize>,
    pub perc_identity: f64,
    pub threads: usize,
    pub plus: bool,
    pub minus: bool,
}

#[derive(Clone, Debug)]
pub struct QHsp {
    pub h: Hsp,
    pub minus: bool,
    pub evalue: f64,
}

#[derive(Clone, Debug)]
pub struct SubjectHit {
    pub gid: u64,
    pub subject_len: usize,
    pub hsps: Vec<QHsp>,
}

impl SubjectHit {
    fn best_evalue(&self) -> f64 {
        self.hsps.iter().map(|h| h.evalue).fold(f64::INFINITY, f64::min)
    }

    fn top_score(&self) -> i32 {
        self.hsps.iter().map(|h| h.h.score).max().unwrap_or(0)
    }
}

/// BLAST's subject ordering: e-value, then score, then OID.
fn hit_order(a: &SubjectHit, b: &SubjectHit) -> Ordering {
    let (ea, eb) = (a.best_evalue(), b.best_evalue());
    if !(ea < 1.0e-180 && eb < 1.0e-180) && ea != eb {
        return ea.partial_cmp(&eb).unwrap_or(Ordering::Equal);
    }
    b.top_score().cmp(&a.top_score()).then(b.gid.cmp(&a.gid))
}

#[derive(Default)]
struct IdHasher(u64);
impl Hasher for IdHasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0 ^ b as u64).wrapping_mul(0x100_0000_01b3);
        }
    }
    fn write_u64(&mut self, v: u64) {
        self.0 = (self.0 ^ v).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        self.0 ^= self.0 >> 29;
    }
    fn write_i64(&mut self, v: i64) {
        self.write_u64(v as u64);
    }
}
type FastMap<K, V> = HashMap<K, V, BuildHasherDefault<IdHasher>>;

struct Seg {
    ctx: u32,
    qs: usize,
    qe: usize,
    ss: usize,
    se: usize,
    score: i32,
}

/// Per-query lookup context shared by all threads.
pub struct QueryCtx {
    /// Masked intervals per strand, sorted, in that strand's coordinates.
    masks: [Vec<(usize, usize)>; 2],
    space: SearchSpace,
}

pub fn make_contexts(queries: &[Query], db: &Database) -> Vec<QueryCtx> {
    queries
        .iter()
        .map(|q| {
            let n = q.len();
            let mut rc: Vec<(usize, usize)> = q.masked.iter().map(|&(a, b)| (n - b, n - a)).collect();
            rc.reverse();
            QueryCtx {
                masks: [q.masked.clone(), rc],
                space: SearchSpace::new(n, db.total_len, db.total_seqs),
            }
        })
        .collect()
}

#[inline(always)]
fn load_kmer(packed: &[u8], bi: usize) -> u64 {
    if bi >= 7 {
        u64::from_be_bytes(packed[bi - 7..=bi].try_into().unwrap())
    } else {
        packed[..=bi].iter().fold(0u64, |a, &b| (a << 8) | b as u64)
    }
}

fn mask_limits(masks: &[(usize, usize)], p: usize, qlen: usize) -> (usize, usize) {
    let i = masks.partition_point(|iv| iv.1 <= p);
    let lo = if i > 0 { masks[i - 1].1 } else { 0 };
    let hi = masks.get(i).map(|iv| iv.0).unwrap_or(qlen);
    (lo, hi)
}

struct Worker<'a> {
    queries: &'a [Query],
    ctxs: &'a [QueryCtx],
    index: &'a Index,
    params: &'a Params,
    last_end: FastMap<(u32, i64), usize>,
    segs: Vec<Seg>,
}

impl<'a> Worker<'a> {
    fn qseq(&self, ctx: u32) -> &'a [u8] {
        let q = &self.queries[(ctx >> 1) as usize];
        if ctx & 1 == 0 { &q.fwd } else { &q.rc }
    }

    /// Probe the index every `stride` bases with a k-mer ending on a byte
    /// boundary. Any exact run of at least `word_size` bases contains such a
    /// window because k + stride - 1 <= word_size.
    fn scan(&mut self, sub: &Subject) {
        let k = self.index.k;
        let step = self.index.stride / 4;
        let packed = sub.packed();
        let nfull = sub.len() / 4;
        let mask = self.index.kmer_mask;
        let mut bi = step - 1;
        while bi < nfull {
            let end = 4 * (bi + 1);
            if end >= k {
                let kmer = load_kmer(packed, bi) & mask;
                if self.index.maybe_contains(kmer) {
                    self.seed_hits(sub, kmer, end - k);
                }
            }
            bi += step;
        }
        // a window ending exactly at the end of a sequence whose length is
        // not a multiple of 4 is not byte-aligned
        let len = sub.len();
        if len % 4 != 0 && len >= k && len.is_multiple_of(self.index.stride) {
            let mut kmer = 0u64;
            for i in len - k..len {
                let c = sub.get(i);
                if c >= 4 {
                    return;
                }
                kmer = (kmer << 2) | c as u64;
            }
            if self.index.maybe_contains(kmer) {
                self.seed_hits(sub, kmer, len - k);
            }
        }
    }

    fn seed_hits(&mut self, sub: &Subject, kmer: u64, spos: usize) {
        let k = self.index.k;
        let w = self.params.word_size;
        let amb = sub.has_ambiguity();
        {
            for &(_, pos) in self.index.lookup(kmer) {
                let ctx = (pos >> 32) as u32;
                if (ctx & 1 == 0 && !self.params.plus) || (ctx & 1 == 1 && !self.params.minus) {
                    continue;
                }
                let qoff = (pos & 0xffff_ffff) as usize;
                let diag = spos as i64 - qoff as i64;
                if let Some(&e) = self.last_end.get(&(ctx, diag))
                    && spos < e {
                        continue;
                    }
                let q = self.qseq(ctx);
                let qlen = q.len();
                // exact run around the seed
                let mut a = 0usize;
                if amb {
                    // verify the seed itself against true residues
                    if (0..k).any(|i| sub.get(spos + i) != q[qoff + i]) {
                        continue;
                    }
                }
                while a < qoff.min(spos) && {
                    let (x, y) = (q[qoff - a - 1], sub.get(spos - a - 1));
                    x == y && x < 4
                } {
                    a += 1;
                }
                let mut b = k;
                while qoff + b < qlen && spos + b < sub.len() && {
                    let (x, y) = (q[qoff + b], sub.get(spos + b));
                    x == y && x < 4
                } {
                    b += 1;
                }
                let strand = (ctx & 1) as usize;
                let masks = &self.ctxs[(ctx >> 1) as usize].masks[strand];
                let (lo, hi) = mask_limits(masks, qoff, qlen);
                let ra = (qoff - a).max(lo);
                let rb = (qoff + b).min(hi);
                if rb < ra + w {
                    continue;
                }
                let sa = spos + ra - qoff;
                let (qs, qe, ss, se, score) =
                    extend::ungapped(q, sub, ra, sa, rb - ra, UNGAPPED_XDROP);
                self.last_end.insert((ctx, diag), se);
                self.segs.push(Seg { ctx, qs, qe, ss, se, score });
            }
        }
    }

    fn extend_segments(&mut self, sub: &Subject, out: &mut Vec<(usize, Vec<QHsp>)>) {
        let mut segs = std::mem::take(&mut self.segs);
        segs.sort_by(|a, b| {
            a.ctx.cmp(&b.ctx).then(b.score.cmp(&a.score)).then(a.qs.cmp(&b.qs)).then(a.ss.cmp(&b.ss))
        });
        let mut per_query: Vec<QHsp> = Vec::new();
        let mut cur_ctx = u32::MAX;
        let mut group: Vec<Hsp> = Vec::new();
        let mut cur_query = usize::MAX;
        let flush = |per_query: &mut Vec<QHsp>, group: &mut Vec<Hsp>, ctx: u32, ctxs: &[QueryCtx]| {
            if ctx == u32::MAX {
                return;
            }
            let space = &ctxs[(ctx >> 1) as usize].space;
            for h in group.drain(..) {
                let evalue = space.evalue(h.score);
                per_query.push(QHsp { h, minus: ctx & 1 == 1, evalue });
            }
        };
        for seg in &segs {
            if seg.ctx != cur_ctx {
                flush(&mut per_query, &mut group, cur_ctx, self.ctxs);
                cur_ctx = seg.ctx;
                if (seg.ctx >> 1) as usize != cur_query {
                    if cur_query != usize::MAX {
                        out.push((cur_query, std::mem::take(&mut per_query)));
                    }
                    cur_query = (seg.ctx >> 1) as usize;
                }
            }
            let mid = (seg.qe - seg.qs) / 2;
            let contained = group.iter().any(|h| {
                h.on_path(seg.qs + mid, seg.ss + mid)
                    || h.envelops(seg.qs, seg.qe, seg.ss, seg.se, seg.score)
            });
            if contained {
                continue;
            }
            let q = self.qseq(seg.ctx);
            let off = extend::gapped_start(q, sub, seg.qs, seg.ss, seg.qe - seg.qs, START_WINDOW);
            let mut h = extend::gapped(q, sub, seg.qs + off, seg.ss + off, GAPPED_XDROP);
            if h.qe == h.qs {
                continue;
            }
            h.compute_stats(q, sub);
            let space = &self.ctxs[(seg.ctx >> 1) as usize].space;
            if space.evalue(h.score) > self.params.evalue {
                continue;
            }
            group.push(h);
        }
        flush(&mut per_query, &mut group, cur_ctx, self.ctxs);
        if cur_query != usize::MAX {
            out.push((cur_query, per_query));
        }
        self.segs = segs;
        self.segs.clear();
    }

    fn process(&mut self, sub: &Subject) -> Vec<(usize, Vec<QHsp>)> {
        self.last_end.clear();
        self.segs.clear();
        self.scan(sub);
        let mut out = Vec::new();
        if !self.segs.is_empty() {
            self.extend_segments(sub, &mut out);
        }
        for (qi, hs) in out.iter_mut() {
            let qlen = self.queries[*qi].len();
            finalize_hsps(hs, qlen, self.params);
        }
        out.retain(|(_, h)| !h.is_empty());
        out
    }
}

fn hsp_order(qlen: usize) -> impl Fn(&QHsp, &QHsp) -> Ordering {
    move |a, b| {
        let qa = a.h.qs + if a.minus { qlen + 1 } else { 0 };
        let qb = b.h.qs + if b.minus { qlen + 1 } else { 0 };
        let ea = a.h.qe + if a.minus { qlen + 1 } else { 0 };
        let eb = b.h.qe + if b.minus { qlen + 1 } else { 0 };
        b.h.score
            .cmp(&a.h.score)
            .then(a.h.ss.cmp(&b.h.ss))
            .then(b.h.se.cmp(&a.h.se))
            .then(qa.cmp(&qb))
            .then(eb.cmp(&ea))
    }
}

/// Order, cull contained HSPs, apply identity/HSP-count limits.
fn finalize_hsps(hs: &mut Vec<QHsp>, qlen: usize, p: &Params) {
    hs.sort_by(hsp_order(qlen));
    let mut keep: Vec<QHsp> = Vec::with_capacity(hs.len());
    for h in hs.drain(..) {
        let dup = keep.iter().any(|k| {
            k.minus == h.minus && k.h.envelops(h.h.qs, h.h.qe, h.h.ss, h.h.se, h.h.score)
        });
        if !dup {
            keep.push(h);
        }
    }
    if p.perc_identity > 0.0 {
        keep.retain(|h| 100.0 * h.h.nident as f64 / h.h.length as f64 >= p.perc_identity);
    }
    if let Some(m) = p.max_hsps {
        keep.truncate(m);
    }
    *hs = keep;
}

fn prune(v: &mut Vec<SubjectHit>, n: usize) {
    v.sort_by(hit_order);
    v.truncate(n);
}

/// Search all queries against the database and return, per query, the
/// ranked subject hits.
pub fn search(
    db: &Database,
    queries: &[Query],
    ctxs: &[QueryCtx],
    index: &Index,
    params: &Params,
) -> Vec<Vec<SubjectHit>> {
    let total = db.num_oids();
    let counter = AtomicU64::new(0);
    let results: Mutex<Vec<Vec<SubjectHit>>> = Mutex::new(vec![Vec::new(); queries.len()]);
    let nthreads = params.threads.max(1);
    std::thread::scope(|sc| {
        for _ in 0..nthreads {
            sc.spawn(|| {
                let mut w = Worker {
                    queries,
                    ctxs,
                    index,
                    params,
                    last_end: FastMap::default(),
                    segs: Vec::new(),
                };
                let mut local: Vec<Vec<SubjectHit>> = vec![Vec::new(); queries.len()];
                loop {
                    let start = counter.fetch_add(CHUNK, AtomicOrdering::Relaxed);
                    if start >= total {
                        break;
                    }
                    for gid in start..(start + CHUNK).min(total) {
                        let (v, oid) = db.locate(gid);
                        let sub = db.volumes[v].subject(oid);
                        if sub.len() < params.word_size {
                            continue;
                        }
                        for (qi, hsps) in w.process(&sub) {
                            let l = &mut local[qi];
                            l.push(SubjectHit { gid, subject_len: sub.len(), hsps });
                            if l.len() > 4 * params.max_target_seqs.max(16) {
                                prune(l, params.max_target_seqs);
                            }
                        }
                    }
                }
                let mut g = results.lock().unwrap();
                for (qi, mut l) in local.into_iter().enumerate() {
                    g[qi].append(&mut l);
                    if g[qi].len() > 4 * params.max_target_seqs.max(16) {
                        prune(&mut g[qi], params.max_target_seqs);
                    }
                }
            });
        }
    });
    let mut res = results.into_inner().unwrap();
    for v in res.iter_mut() {
        prune(v, params.max_target_seqs);
    }
    res
}
