use minimap2::{Aligner, Mapping, Strand};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::OpenOptions;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

static TEMP_ID: AtomicU64 = AtomicU64::new(0);

#[derive(Clone, Debug)]
pub struct Read {
    pub id: usize,
    pub sequence: Vec<u8>,
}

#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub min_overlap: usize,
    pub min_identity: f64,
    pub adaptive_identity: bool,
    pub max_secondary: i32,
    pub threads: usize,
    pub polish: bool,
}

pub struct AssemblyOutput {
    pub contigs: Vec<Vec<u8>>,
    pub input_reads: usize,
    pub retained_reads: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            min_overlap: 100,
            min_identity: 0.80,
            adaptive_identity: true,
            max_secondary: 20,
            threads: std::thread::available_parallelism()
                .map_or(1, usize::from)
                .min(8),
            polish: true,
        }
    }
}

#[derive(Clone, Debug)]
struct Edge {
    from: usize,
    to: usize,
    trim: usize,
    overlap: usize,
    identity: f64,
}

struct TempFasta(PathBuf);

impl TempFasta {
    fn create(records: impl IntoIterator<Item = (String, Vec<u8>)>) -> io::Result<Self> {
        let path = std::env::temp_dir().join(format!(
            "viral-overlap-{}-{}.fa",
            std::process::id(),
            TEMP_ID.fetch_add(1, Ordering::Relaxed)
        ));
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)?;
        let mut writer = BufWriter::new(file);
        for (name, sequence) in records {
            writeln!(writer, ">{name}")?;
            for chunk in sequence.chunks(80) {
                writer.write_all(chunk)?;
                writer.write_all(b"\n")?;
            }
        }
        writer.flush()?;
        Ok(Self(path))
    }
}

impl Drop for TempFasta {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

pub fn auto_min_overlap(reads: &[Read]) -> usize {
    if reads.is_empty() {
        return 100;
    }
    let mut lengths: Vec<_> = reads.iter().map(|read| read.sequence.len()).collect();
    lengths.sort_unstable();
    let n90_index = ((lengths.len() as f64 * 0.10).floor() as usize).min(lengths.len() - 1);
    (lengths[n90_index] / 3).clamp(100, 500)
}

pub fn assemble(reads: &[Read], options: Options) -> io::Result<Vec<Vec<u8>>> {
    assemble_with_stats(reads, options).map(|output| output.contigs)
}

pub fn assemble_with_stats(reads: &[Read], options: Options) -> io::Result<AssemblyOutput> {
    if reads.is_empty() {
        return Ok(AssemblyOutput {
            contigs: Vec::new(),
            input_reads: 0,
            retained_reads: 0,
        });
    }
    if options.min_overlap == 0
        || !options.min_identity.is_finite()
        || !(0.0..=1.0).contains(&options.min_identity)
        || options.max_secondary < 1
        || options.threads == 0
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid overlap assembler options",
        ));
    }
    let normalized_reads = normalize_reads(reads);
    let retained_reads = normalized_reads.len();
    let mut assembly_options = options;
    assembly_options.min_identity =
        effective_min_identity(options, reads.len(), normalized_reads.len());
    let contigs = assemble_normalized(&normalized_reads, assembly_options)?;
    Ok(AssemblyOutput {
        contigs,
        input_reads: reads.len(),
        retained_reads,
    })
}

fn effective_min_identity(options: Options, input_reads: usize, retained_reads: usize) -> f64 {
    const SENSITIVE_IDENTITY: f64 = 0.85;
    const MIN_RETAINED_READ_FRACTION: f64 = 0.10;
    if options.adaptive_identity
        && input_reads > 0
        && retained_reads as f64 / input_reads as f64 >= MIN_RETAINED_READ_FRACTION
    {
        options.min_identity.max(SENSITIVE_IDENTITY)
    } else {
        options.min_identity
    }
}

fn assemble_normalized(reads: &[Read], options: Options) -> io::Result<Vec<Vec<u8>>> {
    let edges = find_overlaps(reads, options)?;
    let mut contigs = layout(reads, edges);
    contigs = deduplicate_reverse_complements(contigs);
    let mut depths = if options.polish {
        polish(reads, &mut contigs, options.threads)?
    } else {
        vec![0; contigs.len()]
    };
    for _ in 0..5 {
        if contigs.len() < 2 {
            break;
        }
        let intermediate: Vec<_> = contigs
            .iter()
            .enumerate()
            .map(|(id, sequence)| Read {
                id,
                sequence: sequence.clone(),
            })
            .collect();
        let merged = deduplicate_reverse_complements(layout(
            &intermediate,
            find_overlaps(&intermediate, options)?,
        ));
        if merged.len() >= contigs.len() {
            break;
        }
        contigs = merged;
        if options.polish {
            depths = polish(reads, &mut contigs, options.threads)?;
        } else {
            depths = vec![0; contigs.len()];
        }
    }
    if options.polish {
        contigs = remove_redundant_contigs(contigs, depths, options.min_overlap);
    } else {
        contigs.retain(|contig| contig.len() >= options.min_overlap);
    }
    contigs.sort_by_key(|contig| std::cmp::Reverse(contig.len()));
    Ok(contigs)
}

fn normalize_reads(reads: &[Read]) -> Vec<Read> {
    const K: usize = 11;
    const TARGET_KMER_DEPTH: u32 = 50;
    const MIN_UNDERCOVERED_PERCENT: usize = 60;

    let mut order: Vec<_> = (0..reads.len()).collect();
    let mut random = 0x9e3779b97f4a7c15u64;
    for index in (1..order.len()).rev() {
        random ^= random << 13;
        random ^= random >> 7;
        random ^= random << 17;
        order.swap(index, random as usize % (index + 1));
    }

    let mut selected_depth = HashMap::<u64, u32>::new();
    let mut selected = vec![false; reads.len()];
    for index in order {
        let kmers = canonical_kmers(&reads[index].sequence, K);
        if kmers.is_empty() {
            selected[index] = true;
            continue;
        }
        let undercovered = kmers
            .iter()
            .filter(|kmer| selected_depth.get(kmer).copied().unwrap_or(0) < TARGET_KMER_DEPTH)
            .count();
        if undercovered * 100 < kmers.len() * MIN_UNDERCOVERED_PERCENT {
            continue;
        }
        selected[index] = true;
        for kmer in kmers {
            let depth = selected_depth.entry(kmer).or_default();
            *depth = depth.saturating_add(1);
        }
    }
    let mut normalized = Vec::new();
    for (read, selected) in reads.iter().zip(selected) {
        if selected {
            normalized.push(Read {
                id: normalized.len(),
                sequence: read.sequence.clone(),
            });
        }
    }
    normalized
}

fn find_overlaps(reads: &[Read], options: Options) -> io::Result<Vec<Edge>> {
    let read_fasta = TempFasta::create(
        reads
            .iter()
            .map(|read| (format!("r{}", read.id), read.sequence.clone())),
    )?;
    let mut builder = Aligner::builder().ava_ont().with_cigar();
    builder.idxopt.k = 11;
    builder.mapopt.best_n = options.max_secondary;
    builder.mapopt.mid_occ = 1000;
    let aligner = builder
        .with_index_threads(options.threads)
        .with_index(&read_fasta.0, None)
        .map_err(|error| io::Error::other(format!("cannot index reads: {error}")))?;
    let mappings: Vec<io::Result<Vec<Mapping>>> = reads
        .par_iter()
        .map(|read| {
            aligner
                .map(
                    &read.sequence,
                    false,
                    false,
                    None,
                    None,
                    Some(format!("r{}", read.id).as_bytes()),
                )
                .map_err(|error| {
                    io::Error::other(format!("read {} mapping failed: {error}", read.id))
                })
        })
        .collect();
    let mut edges = Vec::new();
    for (read, result) in reads.iter().zip(mappings) {
        for mapping in result? {
            add_mapping_edges(read, &mapping, reads, options, &mut edges)?;
        }
    }
    Ok(edges)
}

fn remove_redundant_contigs(
    contigs: Vec<Vec<u8>>,
    depths: Vec<u32>,
    min_contig_length: usize,
) -> Vec<Vec<u8>> {
    const HIGH_DEPTH_THRESHOLD: u32 = 30;
    const HIGH_DEPTH_K: usize = 15;
    const LOW_DEPTH_K: usize = 11;
    const FRAGMENT_K: usize = 9;
    const FRAGMENT_LENGTH_PERCENT: usize = 80;
    const FRAGMENT_DEPTH_RATIO: u32 = 4;
    const HIGH_DEPTH_MATCH_PERCENT: usize = 40;
    const LOW_DEPTH_MATCH_PERCENT: usize = 20;
    let mut contigs: Vec<_> = contigs.into_iter().zip(depths).collect();
    contigs.retain(|(sequence, _)| sequence.len() >= min_contig_length);
    contigs.sort_by_key(|(sequence, _)| std::cmp::Reverse(sequence.len()));
    let mut high_depth_index = std::collections::HashSet::new();
    let mut low_depth_index = std::collections::HashSet::new();
    let mut fragment_index = std::collections::HashSet::new();
    let mut retained = Vec::new();
    let mut low_depth = Vec::new();
    let mut longest_retained_length = 0;
    let mut longest_retained_depth = 0;
    for (sequence, depth) in contigs {
        let high_depth_kmers = canonical_kmers(&sequence, HIGH_DEPTH_K);
        let shared = high_depth_kmers
            .iter()
            .filter(|kmer| high_depth_index.contains(*kmer))
            .count();
        if depth >= HIGH_DEPTH_THRESHOLD {
            if high_depth_kmers.is_empty()
                || shared * 100 < high_depth_kmers.len() * HIGH_DEPTH_MATCH_PERCENT
            {
                high_depth_index.extend(high_depth_kmers);
                low_depth_index.extend(canonical_kmers(&sequence, LOW_DEPTH_K));
                fragment_index.extend(canonical_kmers(&sequence, FRAGMENT_K));
                if sequence.len() > longest_retained_length {
                    longest_retained_length = sequence.len();
                    longest_retained_depth = depth;
                }
                retained.push(sequence);
            }
        } else {
            low_depth.push((sequence, depth));
        }
    }
    for (sequence, depth) in low_depth {
        let is_low_support_fragment = longest_retained_length > 0
            && sequence.len() * 100 < longest_retained_length * FRAGMENT_LENGTH_PERCENT
            && depth.saturating_mul(FRAGMENT_DEPTH_RATIO) < longest_retained_depth;
        let k = if is_low_support_fragment {
            FRAGMENT_K
        } else {
            LOW_DEPTH_K
        };
        let kmers = canonical_kmers(&sequence, k);
        let index = if is_low_support_fragment {
            &fragment_index
        } else {
            &low_depth_index
        };
        let shared = kmers.iter().filter(|kmer| index.contains(*kmer)).count();
        if kmers.is_empty() || shared * 100 < kmers.len() * LOW_DEPTH_MATCH_PERCENT {
            if sequence.len() > longest_retained_length {
                longest_retained_length = sequence.len();
                longest_retained_depth = depth;
            }
            low_depth_index.extend(canonical_kmers(&sequence, LOW_DEPTH_K));
            fragment_index.extend(canonical_kmers(&sequence, FRAGMENT_K));
            retained.push(sequence);
        }
    }
    retained
}

fn canonical_kmers(sequence: &[u8], k: usize) -> Vec<u64> {
    let mask = (1u64 << (2 * k)) - 1;
    let reverse_shift = 2 * (k - 1);
    let (mut forward, mut reverse, mut length) = (0u64, 0u64, 0usize);
    let mut kmers = Vec::new();
    for &base in sequence {
        let Some(base) = base_index(base) else {
            forward = 0;
            reverse = 0;
            length = 0;
            continue;
        };
        forward = ((forward << 2) | base as u64) & mask;
        reverse = (reverse >> 2) | ((3 - base as u64) << reverse_shift);
        length += 1;
        if length >= k {
            kmers.push(forward.min(reverse));
        }
    }
    kmers
}

fn add_mapping_edges(
    query: &Read,
    mapping: &Mapping,
    reads: &[Read],
    options: Options,
    edges: &mut Vec<Edge>,
) -> io::Result<()> {
    let target_name = mapping
        .target_name
        .as_ref()
        .ok_or_else(|| io::Error::other("minimap2 returned a mapping without a target name"))?;
    let Some(target_id) = target_name
        .strip_prefix('r')
        .and_then(|id| id.parse::<usize>().ok())
    else {
        return Err(io::Error::other(format!(
            "unexpected minimap2 target name {target_name}"
        )));
    };
    if target_id >= reads.len() {
        return Err(io::Error::other(format!(
            "minimap2 target id {target_id} is out of range"
        )));
    }
    if target_id == query.id || mapping.block_len <= 0 {
        return Ok(());
    }
    let target = &reads[target_id];
    let target_start = usize::try_from(mapping.target_start)
        .map_err(|_| io::Error::other("negative target start from minimap2"))?;
    let target_end = usize::try_from(mapping.target_end)
        .map_err(|_| io::Error::other("negative target end from minimap2"))?;
    let query_len = query.sequence.len();
    let raw_query_start = usize::try_from(mapping.query_start)
        .map_err(|_| io::Error::other("negative query start from minimap2"))?;
    let raw_query_end = usize::try_from(mapping.query_end)
        .map_err(|_| io::Error::other("negative query end from minimap2"))?;
    let (query_start, query_end, query_reverse) = if mapping.strand == Strand::Reverse {
        (query_len - raw_query_end, query_len - raw_query_start, true)
    } else {
        (raw_query_start, raw_query_end, false)
    };
    let identity = mapping.match_len as f64 / mapping.block_len as f64;
    if identity < options.min_identity
        || target_end > target.sequence.len()
        || query_end > query_len
        || target_start >= target_end
        || query_start >= query_end
    {
        return Ok(());
    }
    let target_span = target_end - target_start;
    let query_span = query_end - query_start;
    if target_span.abs_diff(query_span) > target_span.min(query_span) / 2 {
        return Ok(());
    }
    let overlap = target_span.min(query_span);
    if overlap < options.min_overlap {
        return Ok(());
    }
    let end_slop = (options.min_overlap / 4).clamp(20, 75);
    let target_prefix = target_start <= end_slop;
    let target_suffix = target.sequence.len() - target_end <= end_slop;
    let query_prefix = query_start <= end_slop;
    let query_suffix = query_len - query_end <= end_slop;
    let query_vertex = query.id * 2 + usize::from(query_reverse);
    let target_vertex = target.id * 2;

    if target_suffix && query_prefix {
        push_edge(
            edges,
            target_vertex,
            query_vertex,
            query_end,
            target.sequence.len() - target_start,
            identity,
        );
        push_edge(
            edges,
            query_vertex ^ 1,
            target_vertex ^ 1,
            target.sequence.len() - target_start,
            target.sequence.len() - target_start,
            identity,
        );
    }
    if target_prefix && query_suffix {
        push_edge(
            edges,
            query_vertex,
            target_vertex,
            target_end,
            target_end - target_start,
            identity,
        );
        push_edge(
            edges,
            target_vertex ^ 1,
            query_vertex ^ 1,
            query_len - query_start,
            query_end - query_start,
            identity,
        );
    }
    Ok(())
}

fn push_edge(
    edges: &mut Vec<Edge>,
    from: usize,
    to: usize,
    trim: usize,
    overlap: usize,
    identity: f64,
) {
    if from != to {
        edges.push(Edge {
            from,
            to,
            trim,
            overlap,
            identity,
        });
    }
}

fn oriented_sequence(read: &Read, reverse: bool) -> Vec<u8> {
    if reverse {
        reverse_complement(&read.sequence)
    } else {
        read.sequence.clone()
    }
}

fn vertex_sequence(reads: &[Read], vertex: usize) -> Vec<u8> {
    oriented_sequence(&reads[vertex / 2], vertex % 2 == 1)
}

fn layout(reads: &[Read], mut edges: Vec<Edge>) -> Vec<Vec<u8>> {
    edges.sort_by(|left, right| {
        let left_extension = reads[left.to / 2].sequence.len().saturating_sub(left.trim);
        let right_extension = reads[right.to / 2]
            .sequence
            .len()
            .saturating_sub(right.trim);
        right
            .identity
            .total_cmp(&left.identity)
            .then_with(|| right.overlap.cmp(&left.overlap))
            .then_with(|| right_extension.cmp(&left_extension))
            .then_with(|| left.from.cmp(&right.from))
            .then_with(|| left.to.cmp(&right.to))
    });
    edges.dedup_by(|right, left| right.from == left.from && right.to == left.to);

    let mut successor = vec![None; reads.len() * 2];
    let mut predecessor = vec![false; reads.len() * 2];
    let mut parent: Vec<_> = (0..reads.len() * 2).collect();
    let mut rank = vec![0u8; reads.len() * 2];
    for edge in edges {
        if successor[edge.from].is_some() || predecessor[edge.to] {
            continue;
        }
        let from_root = find(&mut parent, edge.from);
        let to_root = find(&mut parent, edge.to);
        if from_root == to_root {
            continue;
        }
        union(&mut parent, &mut rank, from_root, to_root);
        predecessor[edge.to] = true;
        let from = edge.from;
        successor[from] = Some(edge);
    }

    let mut used_reads = vec![false; reads.len()];
    let mut contigs = Vec::new();
    for start in 0..reads.len() * 2 {
        if predecessor[start] || used_reads[start / 2] {
            continue;
        }
        let mut sequence = vertex_sequence(reads, start);
        let mut current = start;
        loop {
            used_reads[current / 2] = true;
            let Some(edge) = successor[current].as_ref() else {
                break;
            };
            if used_reads[edge.to / 2] {
                break;
            }
            let next = vertex_sequence(reads, edge.to);
            if edge.trim >= next.len() {
                break;
            }
            sequence.extend_from_slice(&next[edge.trim..]);
            current = edge.to;
        }
        contigs.push(sequence);
    }
    contigs
}

fn find(parent: &mut [usize], item: usize) -> usize {
    if parent[item] != item {
        parent[item] = find(parent, parent[item]);
    }
    parent[item]
}

fn union(parent: &mut [usize], rank: &mut [u8], left: usize, right: usize) {
    let (mut left, mut right) = (left, right);
    if rank[left] < rank[right] {
        std::mem::swap(&mut left, &mut right);
    }
    parent[right] = left;
    if rank[left] == rank[right] {
        rank[left] += 1;
    }
}

fn deduplicate_reverse_complements(contigs: Vec<Vec<u8>>) -> Vec<Vec<u8>> {
    let mut seen = std::collections::HashSet::new();
    contigs
        .into_iter()
        .filter_map(|sequence| {
            let reverse = reverse_complement(&sequence);
            let key = if sequence <= reverse {
                sequence.clone()
            } else {
                reverse
            };
            seen.insert(key.clone()).then_some(key)
        })
        .collect()
}

#[derive(Default)]
struct Pileup {
    bases: Vec<[u32; 4]>,
    coverage: Vec<u32>,
    deletions: Vec<u32>,
    boundary_coverage: Vec<u32>,
    insertions: Vec<HashMap<Vec<u8>, u32>>,
}

impl Pileup {
    fn new(length: usize) -> Self {
        Self {
            bases: vec![[0; 4]; length],
            coverage: vec![0; length],
            deletions: vec![0; length],
            boundary_coverage: vec![0; length + 1],
            insertions: vec![HashMap::new(); length + 1],
        }
    }

    fn add_mapping(&mut self, sequence: &[u8], mapping: &Mapping) {
        let Some(alignment) = mapping.alignment.as_ref() else {
            return;
        };
        let Some(cigar) = alignment.cigar.as_ref() else {
            return;
        };
        let reverse = mapping.strand == Strand::Reverse;
        let mut target = mapping.target_start.max(0) as usize;
        let mut query = if reverse {
            mapping.query_end.max(0) as usize
        } else {
            mapping.query_start.max(0) as usize
        };
        self.boundary_coverage[target] += 1;
        let mut insertion = Vec::new();
        for &(length, operation) in cigar {
            match operation {
                0 | 7 | 8 => {
                    for _ in 0..length {
                        if target >= self.bases.len() {
                            return;
                        }
                        let Some(base) = query_base(sequence, query, reverse) else {
                            return;
                        };
                        let base = if reverse { complement(base) } else { base };
                        if let Some(index) = base_index(base) {
                            self.bases[target][index] += 1;
                            self.coverage[target] += 1;
                        }
                        target += 1;
                        query = advance_query(query, reverse, 1);
                        self.boundary_coverage[target] += 1;
                    }
                }
                1 => {
                    for _ in 0..length {
                        if let Some(base) = query_base(sequence, query, reverse) {
                            insertion.push(if reverse {
                                complement(base)
                            } else {
                                base.to_ascii_uppercase()
                            });
                        }
                        query = advance_query(query, reverse, 1);
                    }
                    if !insertion.is_empty() && target < self.insertions.len() {
                        *self.insertions[target]
                            .entry(std::mem::take(&mut insertion))
                            .or_default() += 1;
                    }
                }
                2 => {
                    for _ in 0..length {
                        if target >= self.bases.len() {
                            return;
                        }
                        self.coverage[target] += 1;
                        self.deletions[target] += 1;
                        target += 1;
                        self.boundary_coverage[target] += 1;
                    }
                }
                4 => query = advance_query(query, reverse, length as usize),
                5 | 6 => {}
                _ => return,
            }
        }
    }

    fn consensus(self, reference: &[u8]) -> Vec<u8> {
        let mut output = Vec::with_capacity(reference.len());
        for boundary in 0..=reference.len() {
            if let Some((insertion, count)) = self.insertions[boundary]
                .iter()
                .max_by_key(|(_, count)| *count)
                && *count > self.boundary_coverage[boundary] / 2
            {
                output.extend_from_slice(insertion);
            }
            if boundary == reference.len() {
                break;
            }
            if self.coverage[boundary] > 0 && self.deletions[boundary] * 2 > self.coverage[boundary]
            {
                continue;
            }
            let base = self.bases[boundary]
                .iter()
                .enumerate()
                .max_by_key(|&(_, count)| count)
                .filter(|(_, count)| **count > 0)
                .map(|(index, _)| b"ACGT"[index])
                .unwrap_or(reference[boundary]);
            output.push(base);
        }
        output
    }
}

fn polish(reads: &[Read], contigs: &mut [Vec<u8>], threads: usize) -> io::Result<Vec<u32>> {
    if contigs.is_empty() {
        return Ok(Vec::new());
    }
    let contig_fasta = TempFasta::create(
        contigs
            .iter()
            .enumerate()
            .map(|(index, sequence)| (format!("c{index}"), sequence.clone())),
    )?;
    let mut builder = Aligner::builder().map_ont().with_cigar();
    builder.idxopt.k = 11;
    let aligner = builder
        .with_index_threads(threads)
        .with_index(&contig_fasta.0, None)
        .map_err(|error| io::Error::other(format!("cannot index contigs: {error}")))?;
    let mut pileups: Vec<_> = contigs
        .iter()
        .map(|contig| Pileup::new(contig.len()))
        .collect();
    for read in reads {
        let mappings = aligner
            .map(
                &read.sequence,
                false,
                false,
                None,
                None,
                Some(format!("r{}", read.id).as_bytes()),
            )
            .map_err(|error| io::Error::other(format!("read polishing failed: {error}")))?;
        if let Some(mapping) = mappings.iter().find(|mapping| mapping.is_primary) {
            let Some(name) = mapping.target_name.as_ref() else {
                return Err(io::Error::other(
                    "minimap2 returned a polishing mapping without a target name",
                ));
            };
            let Some(index) = name
                .strip_prefix('c')
                .and_then(|index| index.parse::<usize>().ok())
            else {
                return Err(io::Error::other(format!(
                    "unexpected contig name from minimap2: {name}"
                )));
            };
            if let Some(pileup) = pileups.get_mut(index) {
                pileup.add_mapping(&read.sequence, mapping);
            }
        }
    }
    let mut depths = Vec::with_capacity(contigs.len());
    for (contig, pileup) in contigs.iter_mut().zip(pileups) {
        let mut covered_depths: Vec<_> = pileup
            .coverage
            .iter()
            .copied()
            .filter(|depth| *depth > 0)
            .collect();
        covered_depths.sort_unstable();
        let median_depth = covered_depths
            .get(covered_depths.len() / 2)
            .copied()
            .unwrap_or(0);
        depths.push(median_depth);
        let old = std::mem::take(contig);
        *contig = pileup.consensus(&old);
    }
    Ok(depths)
}

fn query_base(sequence: &[u8], position: usize, reverse: bool) -> Option<u8> {
    let index = if reverse {
        position.checked_sub(1)?
    } else {
        position
    };
    sequence.get(index).copied()
}

fn advance_query(position: usize, reverse: bool, length: usize) -> usize {
    if reverse {
        position.saturating_sub(length)
    } else {
        position.saturating_add(length)
    }
}

fn base_index(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|&base| complement(base))
        .collect()
}

pub fn estimate_depths(
    reads: &[Read],
    contigs: &[Vec<u8>],
    threads: usize,
) -> io::Result<Vec<u64>> {
    if threads == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "threads must be positive",
        ));
    }
    if contigs.is_empty() {
        return Ok(Vec::new());
    }
    if contigs.iter().any(Vec::is_empty) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "contigs must not be empty",
        ));
    }
    let contig_fasta = TempFasta::create(
        contigs
            .iter()
            .enumerate()
            .map(|(index, sequence)| (format!("c{index}"), sequence.clone())),
    )?;
    let mut builder = Aligner::builder().map_ont().with_cigar();
    builder.idxopt.k = 11;
    let aligner = builder
        .with_index_threads(threads)
        .with_index(&contig_fasta.0, None)
        .map_err(|error| io::Error::other(format!("cannot index contigs for depth: {error}")))?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|error| io::Error::other(format!("cannot create depth worker pool: {error}")))?;
    let mappings: Vec<io::Result<Option<(usize, u64)>>> = pool.install(|| {
        reads
            .par_iter()
            .map(|read| {
                let mappings = aligner
                    .map(&read.sequence, false, false, None, None, None)
                    .map_err(|error| {
                        io::Error::other(format!("read {} depth mapping failed: {error}", read.id))
                    })?;
                let Some(mapping) = mappings.iter().find(|mapping| mapping.is_primary) else {
                    return Ok(None);
                };
                let name = mapping.target_name.as_deref().ok_or_else(|| {
                    io::Error::other("minimap2 returned a depth mapping without a target name")
                })?;
                let index = name
                    .strip_prefix('c')
                    .and_then(|id| id.parse::<usize>().ok())
                    .filter(|&id| id < contigs.len())
                    .ok_or_else(|| {
                        io::Error::other(format!("unexpected depth target name {name}"))
                    })?;
                let cigar = mapping
                    .alignment
                    .as_ref()
                    .and_then(|alignment| alignment.cigar.as_ref())
                    .ok_or_else(|| {
                        io::Error::other("minimap2 returned a depth mapping without a CIGAR")
                    })?;
                let aligned_bases = cigar
                    .iter()
                    .filter(|(_, operation)| matches!(operation, 0 | 7 | 8))
                    .map(|(length, _)| u64::from(*length))
                    .sum();
                Ok(Some((index, aligned_bases)))
            })
            .collect()
    });
    let mut aligned_bases = vec![0u64; contigs.len()];
    for mapping in mappings {
        if let Some((index, bases)) = mapping? {
            aligned_bases[index] = aligned_bases[index]
                .checked_add(bases)
                .ok_or_else(|| io::Error::other("depth aligned-base count overflow"))?;
        }
    }
    Ok(contigs
        .iter()
        .zip(aligned_bases)
        .map(|(contig, bases)| {
            let length = contig.len() as u64;
            bases / length + u64::from(bases % length >= length.div_ceil(2))
        })
        .collect())
}

pub fn write_fasta<W: Write>(contigs: &[Vec<u8>], depths: &[u64], mut writer: W) -> io::Result<()> {
    if contigs.len() != depths.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "contig and depth counts differ",
        ));
    }
    for (index, (contig, depth)) in contigs.iter().zip(depths).enumerate() {
        writeln!(
            writer,
            ">ctg_{} LN:i:{} DP:i:{}",
            index + 1,
            contig.len(),
            depth
        )?;
        for chunk in contig.chunks(80) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
    }
    Ok(())
}

pub fn parse_reads(paths: &[PathBuf]) -> io::Result<Vec<Read>> {
    let mut reads = Vec::new();
    for path in paths {
        let mut parser = needletail::parse_fastx_file(path).map_err(|error| {
            io::Error::other(format!("cannot read {}: {error}", path.display()))
        })?;
        while let Some(record) = parser.next() {
            let record = record.map_err(|error| {
                io::Error::other(format!("invalid reads in {}: {error}", path.display()))
            })?;
            let sequence = record.seq().into_owned();
            if sequence.is_empty() {
                continue;
            }
            reads.push(Read {
                id: reads.len(),
                sequence,
            });
        }
    }
    Ok(reads)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fasta_headers_include_length_and_depth_for_each_contig() {
        let contigs = vec![vec![b'A'; 81], b"ACGT".to_vec()];
        let mut output = Vec::new();
        write_fasta(&contigs, &[12, 0], &mut output).unwrap();
        assert_eq!(
            String::from_utf8(output).unwrap(),
            format!(
                ">ctg_1 LN:i:81 DP:i:12\n{}\nA\n>ctg_2 LN:i:4 DP:i:0\nACGT\n",
                "A".repeat(80)
            )
        );
        assert!(write_fasta(&contigs, &[12], Vec::new()).is_err());
    }

    #[test]
    fn depth_counts_original_mapped_bases_across_full_contigs() {
        let mut state = 17u64;
        let mut sequence = || {
            (0..800)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    b"ACGT"[(state >> 62) as usize]
                })
                .collect::<Vec<_>>()
        };
        let first = sequence();
        let second = sequence();
        let uncovered = sequence();
        let reads = vec![
            Read {
                id: 0,
                sequence: first.clone(),
            },
            Read {
                id: 1,
                sequence: reverse_complement(&first),
            },
            Read {
                id: 2,
                sequence: first[200..600].to_vec(),
            },
            Read {
                id: 3,
                sequence: second.clone(),
            },
        ];
        assert_eq!(
            estimate_depths(&reads, &[first, second, uncovered], 2).unwrap(),
            vec![3, 1, 0]
        );
    }

    #[test]
    fn depth_uses_reads_discarded_by_normalization() {
        let mut state = 47u64;
        let contig: Vec<_> = (0..500)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state >> 62) as usize]
            })
            .collect();
        let reads: Vec<_> = (0..120)
            .map(|id| Read {
                id,
                sequence: contig.clone(),
            })
            .collect();
        assert!(normalize_reads(&reads).len() < reads.len());
        assert_eq!(estimate_depths(&reads, &[contig], 2).unwrap(), vec![120]);
    }

    #[test]
    fn auto_overlap_clamps_to_flye_short_read_range() {
        let reads = vec![
            Read {
                id: 0,
                sequence: vec![b'A'; 180],
            },
            Read {
                id: 1,
                sequence: vec![b'C'; 360],
            },
        ];
        assert_eq!(auto_min_overlap(&reads), 100);
        let long = vec![Read {
            id: 0,
            sequence: vec![b'A'; 2400],
        }];
        assert_eq!(auto_min_overlap(&long), 500);
    }

    #[test]
    fn adaptive_identity_uses_normalization_to_protect_high_depth() {
        let options = Options::default();
        assert_eq!(effective_min_identity(options, 1000, 1000), 0.85);
        assert_eq!(effective_min_identity(options, 1000, 99), 0.80);
        let fixed = Options {
            adaptive_identity: false,
            ..options
        };
        assert_eq!(effective_min_identity(fixed, 1000, 1000), 0.80);
    }

    #[test]
    fn reverse_complement_is_involutive() {
        let sequence = b"ACGTTAGC";
        assert_eq!(reverse_complement(&reverse_complement(sequence)), sequence);
    }

    #[test]
    fn removes_near_identical_contained_contigs() {
        let mut state = 17u64;
        let sequence: Vec<_> = (0..1000)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(state >> 32) as usize % 4]
            })
            .collect();
        let mut contained = sequence[100..900].to_vec();
        let mutation = base_index(contained[400]).unwrap();
        contained[400] = b"ACGT"[(mutation + 1) % 4];
        assert_eq!(
            remove_redundant_contigs(vec![sequence.clone(), contained.clone()], vec![100, 5], 100,),
            vec![sequence.clone()]
        );
        assert_eq!(
            remove_redundant_contigs(vec![sequence.clone(), contained], vec![10, 5], 100),
            vec![sequence]
        );
    }

    #[test]
    fn retains_low_depth_contigs_with_divergent_kmers() {
        let mut state = 17u64;
        let sequence: Vec<_> = (0..1000)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(state >> 32) as usize % 4]
            })
            .collect();
        let mut variant = sequence.clone();
        for position in (0..variant.len()).step_by(7) {
            let current = base_index(variant[position]).unwrap();
            variant[position] = b"ACGT"[(current + 1) % 4];
        }
        assert_eq!(
            remove_redundant_contigs(vec![sequence, variant], vec![5, 5], 100).len(),
            2
        );
    }

    #[test]
    fn retains_supported_divergent_contigs_below_a_longer_contig() {
        let mut state = 17u64;
        let sequence: Vec<_> = (0..1000)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(state >> 32) as usize % 4]
            })
            .collect();
        let mut variant = sequence[100..800].to_vec();
        for base in &mut variant {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            if state % 7 == 0 {
                let current = base_index(*base).unwrap();
                *base = b"ACGT"[(current + 1) % 4];
            }
        }
        assert_eq!(
            remove_redundant_contigs(vec![sequence, variant], vec![10, 10], 100).len(),
            2
        );
    }

    #[test]
    fn normalizes_high_coverage_but_retains_low_coverage_reads() {
        let mut state = 17u64;
        let reference: Vec<_> = (0..500)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(state >> 32) as usize % 4]
            })
            .collect();
        let low_coverage: Vec<_> = (0..10)
            .map(|id| Read {
                id,
                sequence: reference.clone(),
            })
            .collect();
        assert_eq!(normalize_reads(&low_coverage).len(), low_coverage.len());

        let mut high_coverage = Vec::new();
        for id in 0..400 {
            let mut sequence = reference.clone();
            for base in &mut sequence {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                if state.is_multiple_of(20) {
                    let current = base_index(*base).unwrap();
                    *base = b"ACGT"[(current + 1) % 4];
                }
            }
            high_coverage.push(Read { id, sequence });
        }
        let normalized = normalize_reads(&high_coverage);
        let retained = normalized.len();
        assert!(
            retained > 10 && retained < high_coverage.len(),
            "{retained}"
        );
        assert!(
            normalized
                .iter()
                .enumerate()
                .all(|(index, read)| read.id == index)
        );
    }

    #[test]
    fn pileup_applies_majority_substitution_and_deletion() {
        let reference = b"ACGT";
        let mut pileup = Pileup::new(reference.len());
        pileup.bases[0][1] = 3;
        pileup.coverage[0] = 3;
        pileup.coverage[1] = 3;
        pileup.deletions[1] = 2;
        pileup.coverage[2] = 3;
        pileup.bases[2][0] = 2;
        pileup.bases[2][2] = 1;
        pileup.coverage[3] = 3;
        assert_eq!(pileup.consensus(reference), b"CAT");
    }

    #[test]
    fn minimap2_pileup_corrects_substitutions_on_both_strands() {
        let mut state = 17u64;
        let reference: Vec<u8> = (0..1200)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state >> 62) as usize]
            })
            .collect();
        let mut reads = Vec::new();
        for start in (0..800).step_by(40) {
            let end = (start + 400).min(reference.len());
            let mut sequence = reference[start..end].to_vec();
            for base in &mut sequence {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                if state % 17 == 0 {
                    *base = match *base {
                        b'A' => b'C',
                        _ => b'A',
                    };
                }
            }
            if reads.len().is_multiple_of(2) {
                sequence = reverse_complement(&sequence);
            }
            reads.push(Read {
                id: reads.len(),
                sequence,
            });
        }
        let mut contigs = vec![reference.clone()];
        polish(&reads, &mut contigs, 2).unwrap();
        let differences = contigs[0]
            .iter()
            .zip(&reference)
            .filter(|(left, right)| left != right)
            .count()
            + contigs[0].len().abs_diff(reference.len());
        assert!(
            differences <= 10,
            "consensus length {} vs reference {}, differences {differences}",
            contigs[0].len(),
            reference.len()
        );
    }

    #[test]
    fn reverse_complement_read_votes_in_reference_orientation() {
        let mut state = 17u64;
        let reference: Vec<u8> = (0..1200)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state >> 62) as usize]
            })
            .collect();
        let reads = vec![Read {
            id: 0,
            sequence: reverse_complement(&reference),
        }];
        let mut contigs = vec![reference.clone()];
        polish(&reads, &mut contigs, 1).unwrap();
        assert_eq!(contigs[0], reference);
    }
}
