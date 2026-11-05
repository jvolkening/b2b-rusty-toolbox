//! K-mer index of the query batch: a bloom filter in front of an
//! open-addressing hash table.

use crate::fasta::Query;

pub struct Index {
    pub k: usize,
    /// Subject scan stride in bases (a multiple of 4).
    pub stride: usize,
    /// Cache-resident prefilter over the low 12 bases of each k-mer.
    pre: Vec<u64>,
    bloom: Vec<u64>,
    bloom_shift: u32,
    /// Open-addressing table of indexes into `pairs` (first entry per k-mer).
    slots: Vec<u32>,
    mask: usize,
    /// (k-mer, ctx << 32 | offset), sorted.
    pairs: Vec<(u64, u64)>,
    pub kmer_mask: u64,
}

const EMPTY: u32 = u32::MAX;

#[inline(always)]
fn mix(key: u64) -> u64 {
    key.wrapping_mul(0x9E37_79B9_7F4A_7C15)
}

impl Index {
    /// Build from all unmasked, unambiguous k-mers of both strands. Context
    /// `2*i` is the forward strand of query `i` and `2*i+1` its reverse
    /// complement.
    pub fn build(queries: &[Query], k: usize, stride: usize) -> Index {
        let kmer_mask = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
        let total: usize = queries.iter().map(|q| 2 * q.len()).sum();
        let mut pairs: Vec<(u64, u64)> = Vec::with_capacity(total);
        for (qi, q) in queries.iter().enumerate() {
            for strand in 0..2 {
                let seq = if strand == 0 { &q.fwd } else { &q.rc };
                let n = seq.len();
                let ctx = (qi as u64) * 2 + strand as u64;
                let mut bad = vec![false; n];
                for &(a, b) in &q.masked {
                    let (a, b) = if strand == 0 { (a, b) } else { (n - b, n - a) };
                    bad[a..b].iter_mut().for_each(|x| *x = true);
                }
                let (mut kmer, mut valid) = (0u64, 0usize);
                for p in 0..n {
                    let c = seq[p];
                    if c < 4 && !bad[p] {
                        kmer = ((kmer << 2) | c as u64) & kmer_mask;
                        valid += 1;
                        if valid >= k {
                            pairs.push((kmer, (ctx << 32) | (p + 1 - k) as u64));
                        }
                    } else {
                        valid = 0;
                        kmer = 0;
                    }
                }
            }
        }
        pairs.sort_unstable();
        pairs.shrink_to_fit();
        let uniq = pairs.windows(2).filter(|w| w[0].0 != w[1].0).count() + !pairs.is_empty() as usize;
        let cap = (uniq * 2).max(16).next_power_of_two();
        let mut slots = vec![EMPTY; cap];
        let mask = cap - 1;
        let bloom_bits = ((uniq as u64 * 32).max(1 << 16)).next_power_of_two().min(1 << 30);
        let bloom_shift = 64 - bloom_bits.trailing_zeros();
        let mut bloom = vec![0u64; (bloom_bits / 64) as usize];
        let mut pre = vec![0u64; 1 << 18];
        for i in 0..pairs.len() {
            let key = pairs[i].0;
            if i > 0 && pairs[i - 1].0 == key {
                continue;
            }
            let h = mix(key);
            let mut slot = (h >> 20) as usize & mask;
            while slots[slot] != EMPTY {
                slot = (slot + 1) & mask;
            }
            slots[slot] = i as u32;
            let bit = (h >> bloom_shift) as usize;
            bloom[bit >> 6] |= 1 << (bit & 63);
            let p = (key & 0xff_ffff) as usize;
            pre[p >> 6] |= 1 << (p & 63);
        }
        Index { k, stride, pre, bloom, bloom_shift, slots, mask, pairs, kmer_mask }
    }

    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    #[inline(always)]
    pub fn maybe_contains(&self, key: u64) -> bool {
        let p = (key & 0xff_ffff) as usize;
        if self.pre[p >> 6] >> (p & 63) & 1 == 0 {
            return false;
        }
        let bit = (mix(key) >> self.bloom_shift) as usize;
        self.bloom[bit >> 6] >> (bit & 63) & 1 != 0
    }

    /// All (k-mer, ctx << 32 | offset) entries of a k-mer.
    #[inline]
    pub fn lookup(&self, key: u64) -> &[(u64, u64)] {
        let mut slot = (mix(key) >> 20) as usize & self.mask;
        loop {
            let i = self.slots[slot];
            if i == EMPTY {
                return &[];
            }
            let i = i as usize;
            if self.pairs[i].0 == key {
                let n = self.pairs[i..].iter().take_while(|e| e.0 == key).count();
                return &self.pairs[i..i + n];
            }
            slot = (slot + 1) & self.mask;
        }
    }
}

/// Choose the k-mer length and subject stride for a word size: the largest
/// stride whose k-mers (k = word_size - stride + 1) stay selective enough
/// that most probes miss.
pub fn choose_seed(word_size: usize, nkmers: usize) -> (usize, usize) {
    for stride in [16usize, 12, 8] {
        if word_size + 1 <= stride {
            continue;
        }
        let k = (word_size + 1 - stride).min(32);
        if k >= 10 && (4f64).powi(k as i32) >= 16.0 * nkmers as f64 {
            return (k, stride);
        }
    }
    (word_size.min(32) - 3, 4)
}
