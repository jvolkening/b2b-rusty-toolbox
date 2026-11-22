use std::collections::{HashMap, HashSet, hash_map::Entry};
use std::io::{self, BufRead, Write};

#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub k: usize,
    pub min_count: u32,
    pub min_fraction: f64,
    pub min_length: usize,
    pub max_kmers: usize,
    pub clean_bubbles: bool,
    pub clean_errors: bool,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            k: 21,
            min_count: 3,
            min_fraction: 0.2,
            min_length: 0,
            max_kmers: 10_000_000,
            clean_bubbles: false,
            clean_errors: true,
        }
    }
}

#[derive(Default)]
struct Edge {
    count: u32,
    visited: bool,
}

#[derive(Default)]
struct Node {
    max_in: u32,
    max_out: u32,
    incoming: u8,
    outgoing: u8,
}

#[derive(Default, Debug)]
pub struct AssemblyStats {
    pub unitigs: usize,
    pub removed_bubble_paths: usize,
    pub removed_bubble_kmers: usize,
    pub removed_error_kmers: usize,
    pub separated_uncertain_unitigs: usize,
    pub discarded_uncertain_unitigs: usize,
    pub redundancy_search_limit: usize,
}

pub struct Assembler {
    options: Options,
    edges: HashMap<u64, Edge>,
    pub reads: u64,
    pub bases: u64,
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

fn reverse_complement(word: u64, k: usize) -> u64 {
    let swapped = ((word & 0x5555_5555_5555_5555) << 1) | ((word & 0xaaaa_aaaa_aaaa_aaaa) >> 1);
    (swapped.reverse_bits() >> (64 - 2 * k)) ^ ((1u64 << (2 * k)) - 1)
}

fn canonical(word: u64, k: usize) -> u64 {
    word.min(reverse_complement(word, k))
}

fn nodes(word: u64, k: usize) -> (u64, u64) {
    (word >> 2, word & ((1u64 << (2 * (k - 1))) - 1))
}

fn decode(word: u64, len: usize, out: &mut Vec<u8>) {
    for i in (0..len).rev() {
        out.push(b"ACGT"[((word >> (2 * i)) & 3) as usize]);
    }
}

fn encode(sequence: &[u8]) -> u64 {
    sequence.iter().fold(0, |word, base| {
        (word << 2)
            | match base {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => unreachable!("unitigs contain only A/C/G/T"),
            }
    })
}

fn one_error_neighbors(word: u64, k: usize) -> Vec<u64> {
    let mut neighbors = Vec::with_capacity(20 * k);
    for_each_one_error_neighbor(word, k, |candidate| neighbors.push(candidate));
    neighbors
}

fn for_each_one_error_neighbor(word: u64, k: usize, mut visit: impl FnMut(u64)) {
    let short_mask = (1u64 << (2 * (k - 1))) - 1;
    let without_first = word & short_mask;
    let without_last = word >> 2;
    for pos in 0..k {
        let shift = 2 * pos;
        let low_mask = (1u64 << shift) - 1;
        let deleted = ((word >> (shift + 2)) << shift) | (word & low_mask);
        for base in 0..4 {
            visit((word & !(3 << shift)) | (base << shift));
            visit((deleted << 2) | base);
            visit((base << (2 * (k - 1))) | deleted);
            for shorter in [without_first, without_last] {
                visit(((shorter >> shift) << (shift + 2)) | (base << shift) | (shorter & low_mask));
            }
        }
    }
}

fn write_record<W: Write>(
    writer: &mut W,
    sequence: &[u8],
    circular: bool,
    prefix: &str,
    count: &mut usize,
) -> io::Result<()> {
    *count += 1;
    writeln!(
        writer,
        ">{prefix}_{} length={} circular={circular}",
        *count,
        sequence.len()
    )?;
    for chunk in sequence.chunks(80) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

fn within_two_edits(a: &[u8], b: &[u8]) -> bool {
    if a.len().abs_diff(b.len()) > 2 {
        return false;
    }
    let mut previous: Vec<usize> = (0..=b.len()).collect();
    let mut current = vec![0; b.len() + 1];
    for (i, &left) in a.iter().enumerate() {
        current[0] = i + 1;
        for (j, &right) in b.iter().enumerate() {
            current[j + 1] = (previous[j + 1] + 1)
                .min(current[j] + 1)
                .min(previous[j] + usize::from(left != right));
        }
        std::mem::swap(&mut previous, &mut current);
    }
    previous[b.len()] <= 2
}

struct Path {
    words: Vec<u64>,
    sequence: Vec<u8>,
    min_count: u32,
    max_count: u32,
    end: u64,
}

impl Assembler {
    pub fn new(options: Options) -> io::Result<Self> {
        if !(15..=31).contains(&options.k) || options.k.is_multiple_of(2) {
            return Err(invalid("k must be an odd number between 15 and 31"));
        }
        if options.min_count == 0 || options.max_kmers == 0 {
            return Err(invalid("min-count and max-kmers must be positive"));
        }
        if !options.min_fraction.is_finite() || !(0.0..1.0).contains(&options.min_fraction) {
            return Err(invalid("min-fraction must be at least 0 and less than 1"));
        }
        Ok(Self {
            options,
            edges: HashMap::new(),
            reads: 0,
            bases: 0,
        })
    }

    pub fn add_read(&mut self, sequence: &[u8]) -> io::Result<()> {
        if sequence.is_empty() {
            return Err(invalid("empty sequence"));
        }
        let k = self.options.k;
        let mask = (1u64 << (2 * k)) - 1;
        let (mut forward, mut reverse, mut run) = (0u64, 0u64, 0usize);
        for &base in sequence {
            let value = match base.to_ascii_uppercase() {
                b'A' => 0u64,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                b'N' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' => {
                    forward = 0;
                    reverse = 0;
                    run = 0;
                    continue;
                }
                _ => return Err(invalid(format!("invalid sequence character: {base:?}"))),
            };
            forward = ((forward << 2) | value) & mask;
            reverse = (reverse >> 2) | ((3 - value) << (2 * (k - 1)));
            run += 1;
            if run >= k {
                let at_limit = self.edges.len() >= self.options.max_kmers;
                match self.edges.entry(forward.min(reverse)) {
                    Entry::Occupied(mut entry) => {
                        entry.get_mut().count = entry.get().count.saturating_add(1);
                    }
                    Entry::Vacant(entry) => {
                        if at_limit {
                            return Err(invalid(format!(
                                "distinct k-mer limit ({}) exceeded; increase --max-kmers if memory permits",
                                self.options.max_kmers
                            )));
                        }
                        entry.insert(Edge {
                            count: 1,
                            visited: false,
                        });
                    }
                }
            }
        }
        self.reads = self.reads.saturating_add(1);
        self.bases = self.bases.saturating_add(sequence.len() as u64);
        Ok(())
    }

    pub fn read_fastx<R: BufRead>(&mut self, mut reader: R) -> io::Result<()> {
        let mut line = Vec::new();
        if read_line(&mut reader, &mut line)? == 0 {
            return Err(invalid("empty input"));
        }
        match line.first() {
            Some(b'>') => {
                let mut sequence = Vec::new();
                loop {
                    if read_line(&mut reader, &mut line)? == 0 {
                        self.add_read(&sequence)?;
                        break;
                    }
                    if line.starts_with(b">") {
                        self.add_read(&sequence)?;
                        sequence.clear();
                    } else {
                        sequence.extend_from_slice(&line);
                    }
                }
            }
            Some(b'@') => loop {
                let mut sequence = Vec::new();
                loop {
                    if read_line(&mut reader, &mut line)? == 0 {
                        return Err(invalid("truncated FASTQ sequence"));
                    }
                    if line.starts_with(b"+") {
                        break;
                    }
                    sequence.extend_from_slice(&line);
                }
                if sequence.is_empty() {
                    return Err(invalid("empty FASTQ sequence"));
                }
                let mut quality_len = 0usize;
                while quality_len < sequence.len() {
                    if read_line(&mut reader, &mut line)? == 0 {
                        return Err(invalid("truncated FASTQ quality"));
                    }
                    if line.iter().any(|b| !(33..=126).contains(b)) {
                        return Err(invalid("invalid FASTQ quality character"));
                    }
                    quality_len += line.len();
                }
                if quality_len != sequence.len() {
                    return Err(invalid("FASTQ sequence and quality lengths differ"));
                }
                self.add_read(&sequence)?;
                if read_line(&mut reader, &mut line)? == 0 {
                    break;
                }
                if !line.starts_with(b"@") {
                    return Err(invalid("expected FASTQ header"));
                }
            },
            _ => return Err(invalid("input must be uncompressed FASTA or FASTQ")),
        }
        Ok(())
    }

    pub fn write_unitigs<W: Write>(self, writer: W) -> io::Result<AssemblyStats> {
        self.write_unitigs_with_uncertain(writer, None)
    }

    pub fn write_unitigs_with_uncertain<W: Write>(
        mut self,
        mut writer: W,
        mut uncertain: Option<&mut dyn Write>,
    ) -> io::Result<AssemblyStats> {
        let k = self.options.k;
        self.edges
            .retain(|_, edge| edge.count >= self.options.min_count);

        let fraction = self.options.min_fraction;
        let mut stats = AssemblyStats::default();
        if self.options.clean_errors && fraction > 0.0 {
            stats.removed_error_kmers = self.clean_error_kmers(false);
        }
        let mut graph: HashMap<u64, Node> = HashMap::new();
        for (&word, edge) in &self.edges {
            for oriented in [word, reverse_complement(word, k)] {
                let (from, to) = nodes(oriented, k);
                let source = graph.entry(from).or_default();
                source.max_out = source.max_out.max(edge.count);
                source.outgoing += 1;
                let target = graph.entry(to).or_default();
                target.max_in = target.max_in.max(edge.count);
                target.incoming += 1;
            }
        }

        if self.options.clean_bubbles && fraction > 0.0 {
            let (removed, paths) = self.weak_bubbles(&graph);
            stats.removed_bubble_paths = paths;
            stats.removed_bubble_kmers = removed.len();
            self.edges.retain(|word, _| !removed.contains(word));
        }
        self.edges.retain(|&word, edge| {
            let (from, to) = nodes(word, k);
            (edge.count as f64) >= fraction * (graph[&from].max_out as f64)
                && (edge.count as f64) >= fraction * (graph[&to].max_in as f64)
        });
        if self.options.clean_errors && fraction > 0.0 {
            stats.removed_error_kmers += self.clean_error_kmers(true);
        }
        for node in graph.values_mut() {
            node.incoming = 0;
            node.outgoing = 0;
        }
        for &word in self.edges.keys() {
            for oriented in [word, reverse_complement(word, k)] {
                let (from, to) = nodes(oriented, k);
                graph.get_mut(&from).unwrap().outgoing += 1;
                graph.get_mut(&to).unwrap().incoming += 1;
            }
        }
        let mut words: Vec<u64> = self.edges.keys().copied().collect();
        words.sort_unstable();
        let mut short_unitigs = Vec::new();

        // Start at junctions first, then consume isolated cycles.
        for cycles in [false, true] {
            for &word in &words {
                for oriented in [word, reverse_complement(word, k)] {
                    let (from, _) = nodes(oriented, k);
                    let source = &graph[&from];
                    if !cycles && source.incoming == 1 && source.outgoing == 1 {
                        continue;
                    }
                    if self.edges[&word].visited {
                        continue;
                    }
                    let (sequence, circular) = self.walk(oriented, &graph);
                    if sequence.len() < self.options.min_length {
                        continue;
                    }
                    if self.options.clean_errors && !circular && sequence.len() <= 2 * k + 2 {
                        short_unitigs.push(sequence);
                    } else {
                        write_record(
                            &mut writer,
                            &sequence,
                            circular,
                            "unitig",
                            &mut stats.unitigs,
                        )?;
                    }
                }
            }
        }
        for sequence in short_unitigs {
            let (redundant, exhausted) = self.redundant_unitig(&sequence);
            stats.redundancy_search_limit += usize::from(exhausted);
            if redundant {
                if let Some(output) = uncertain.as_mut() {
                    write_record(
                        output,
                        &sequence,
                        false,
                        "uncertain",
                        &mut stats.separated_uncertain_unitigs,
                    )?;
                } else {
                    stats.discarded_uncertain_unitigs += 1;
                }
            } else {
                write_record(&mut writer, &sequence, false, "unitig", &mut stats.unitigs)?;
            }
        }
        writer.flush()?;
        if let Some(uncertain) = uncertain {
            uncertain.flush()?;
        }
        Ok(stats)
    }

    fn redundant_unitig(&self, sequence: &[u8]) -> (bool, bool) {
        const MAX_CANDIDATES: usize = 50_000;
        let k = self.options.k;
        if self.options.min_fraction == 0.0 {
            return (false, false);
        }
        let weak_max = sequence
            .windows(k)
            .map(|w| self.edges[&canonical(encode(w), k)].count)
            .max()
            .unwrap_or(0);
        let seed = encode(&sequence[..k]);
        let mut seen = HashSet::new();
        let first = one_error_neighbors(seed, k);
        for candidate in first
            .iter()
            .copied()
            .chain(first.iter().flat_map(|&word| one_error_neighbors(word, k)))
        {
            if !seen.insert(candidate) {
                continue;
            }
            if seen.len() > MAX_CANDIDATES {
                return (false, true);
            }
            let key = canonical(candidate, k);
            let Some(edge) = self.edges.get(&key) else {
                continue;
            };
            if weak_max as f64 >= self.options.min_fraction * edge.count as f64 {
                continue;
            }
            let context = self.strong_context(candidate, sequence.len(), weak_max);
            if context.len() <= sequence.len() {
                continue;
            }
            for start in 0..=context
                .len()
                .saturating_sub(sequence.len().saturating_sub(2))
            {
                for length in sequence.len().saturating_sub(2)..=sequence.len() + 2 {
                    if start + length <= context.len()
                        && within_two_edits(sequence, &context[start..start + length])
                    {
                        return (true, false);
                    }
                }
            }
        }
        (false, false)
    }

    fn strong_context(&self, seed: u64, length: usize, weak_max: u32) -> Vec<u8> {
        let k = self.options.k;
        let mut result = Vec::new();
        decode(seed, k, &mut result);
        let mut used = HashSet::from([canonical(seed, k)]);
        let mut word = seed;
        for _ in 0..4 {
            let node = word >> 2;
            let next = (0..4)
                .map(|base| (base << (2 * (k - 1))) | node)
                .filter_map(|candidate| {
                    self.edges
                        .get(&canonical(candidate, k))
                        .map(|edge| (edge.count, candidate))
                })
                .max();
            let Some((count, candidate)) = next else {
                break;
            };
            if weak_max as f64 >= self.options.min_fraction * count as f64
                || !used.insert(canonical(candidate, k))
            {
                break;
            }
            result.insert(0, b"ACGT"[(candidate >> (2 * (k - 1))) as usize]);
            word = candidate;
        }
        word = seed;
        for _ in 0..length.saturating_sub(k) + 4 {
            let node = nodes(word, k).1;
            let next = (0..4)
                .map(|base| (node << 2) | base)
                .filter_map(|candidate| {
                    self.edges
                        .get(&canonical(candidate, k))
                        .map(|edge| (edge.count, candidate))
                })
                .max();
            let Some((count, candidate)) = next else {
                break;
            };
            if weak_max as f64 >= self.options.min_fraction * count as f64
                || !used.insert(canonical(candidate, k))
            {
                break;
            }
            result.push(b"ACGT"[(candidate & 3) as usize]);
            word = candidate;
        }
        result
    }

    fn comparable_context(&self, weak: u64, strong: u64, strong_count: u32) -> bool {
        let k = self.options.k;
        let mut sequence = Vec::new();
        decode(weak, k, &mut sequence);
        let mut weak_max = self.edges[&canonical(weak, k)].count;
        let mut used = HashSet::from([canonical(weak, k)]);
        let mut word = weak;
        for _ in 0..4 {
            let node = word >> 2;
            let next = (0..4)
                .map(|base| (base << (2 * (k - 1))) | node)
                .filter_map(|candidate| {
                    self.edges
                        .get(&canonical(candidate, k))
                        .filter(|edge| {
                            (edge.count as f64) < self.options.min_fraction * strong_count as f64
                        })
                        .map(|edge| (edge.count, candidate))
                })
                .max();
            let Some((count, candidate)) = next else {
                return false;
            };
            if !used.insert(canonical(candidate, k)) {
                return false;
            }
            weak_max = weak_max.max(count);
            sequence.insert(0, b"ACGT"[(candidate >> (2 * (k - 1))) as usize]);
            word = candidate;
        }
        word = weak;
        for _ in 0..4 {
            let node = nodes(word, k).1;
            let next = (0..4)
                .map(|base| (node << 2) | base)
                .filter_map(|candidate| {
                    self.edges
                        .get(&canonical(candidate, k))
                        .filter(|edge| {
                            (edge.count as f64) < self.options.min_fraction * strong_count as f64
                        })
                        .map(|edge| (edge.count, candidate))
                })
                .max();
            let Some((count, candidate)) = next else {
                return false;
            };
            if !used.insert(canonical(candidate, k)) {
                return false;
            }
            weak_max = weak_max.max(count);
            sequence.push(b"ACGT"[(candidate & 3) as usize]);
            word = candidate;
        }
        let context = self.strong_context(strong, sequence.len(), weak_max);
        for start in 0..=context
            .len()
            .saturating_sub(sequence.len().saturating_sub(2))
        {
            for length in sequence.len().saturating_sub(2)..=sequence.len() + 2 {
                if start + length <= context.len()
                    && within_two_edits(&sequence, &context[start..start + length])
                {
                    return true;
                }
            }
        }
        false
    }

    fn clean_error_kmers(&mut self, two_substitutions: bool) -> usize {
        let max_count = self
            .edges
            .values()
            .map(|edge| edge.count)
            .max()
            .unwrap_or(0);
        let mut weak = HashMap::new();
        for (&word, edge) in &self.edges {
            if edge.count as f64 >= self.options.min_fraction * max_count as f64 {
                continue;
            }
            let neighbor = if two_substitutions {
                self.strongest_two_substitutions(word, edge.count)
            } else {
                self.strongest_similar_kmer(word, edge.count)
            };
            if let Some((neighbor, _)) = neighbor {
                weak.insert(word, neighbor);
            }
        }
        let removed: HashSet<u64> = weak
            .iter()
            .filter_map(|(&word, neighbor)| (!weak.contains_key(neighbor)).then_some(word))
            .collect();
        self.edges.retain(|word, _| !removed.contains(word));
        removed.len()
    }

    fn strongest_two_substitutions(&self, word: u64, count: u32) -> Option<(u64, u32)> {
        let k = self.options.k;
        let mut best = None;
        for i in 0..k {
            let shift_i = 2 * i;
            let old_i = (word >> shift_i) & 3;
            for j in i + 1..k {
                let shift_j = 2 * j;
                let old_j = (word >> shift_j) & 3;
                let remainder = word & !(3 << shift_i) & !(3 << shift_j);
                for a in 0..4 {
                    if a == old_i {
                        continue;
                    }
                    for b in 0..4 {
                        if b == old_j {
                            continue;
                        }
                        let candidate = remainder | (a << shift_i) | (b << shift_j);
                        let key = canonical(candidate, k);
                        if key == word {
                            continue;
                        }
                        if let Some(edge) = self.edges.get(&key)
                            && (count as f64) < self.options.min_fraction * edge.count as f64
                            && best.is_none_or(|(_, best_count)| edge.count > best_count)
                            && self.comparable_context(word, candidate, edge.count)
                        {
                            best = Some((key, edge.count));
                        }
                    }
                }
            }
        }
        best
    }

    fn strongest_similar_kmer(&self, word: u64, count: u32) -> Option<(u64, u32)> {
        let k = self.options.k;
        let mut best = None;
        let mut consider = |candidate: u64| {
            let key = canonical(candidate, k);
            if key == word {
                return;
            }
            if let Some(edge) = self.edges.get(&key)
                && (count as f64) < self.options.min_fraction * (edge.count as f64)
                && best.is_none_or(|(_, best_count)| edge.count > best_count)
                && self.comparable_context(word, candidate, edge.count)
            {
                best = Some((key, edge.count));
            }
        };
        for_each_one_error_neighbor(word, k, &mut consider);
        best
    }

    fn weak_bubbles(&self, graph: &HashMap<u64, Node>) -> (HashSet<u64>, usize) {
        let k = self.options.k;
        let mut removed = HashSet::new();
        let mut protected = HashSet::new();
        let mut paths = 0;
        let mut words: Vec<u64> = self.edges.keys().copied().collect();
        words.sort_unstable();
        for word in words {
            let edge = &self.edges[&word];
            for first in [word, reverse_complement(word, k)] {
                let (start, _) = nodes(first, k);
                let source = &graph[&start];
                if source.outgoing < 2
                    || edge.count as f64 >= self.options.min_fraction * source.max_out as f64
                {
                    continue;
                }
                let Some(weak) = self.trace_weak_path(first, graph) else {
                    continue;
                };
                if graph[&weak.end].incoming < 2 {
                    continue;
                }
                let Some(strong) = self.trace_strong_path(start, weak.end, &weak.words) else {
                    continue;
                };
                if weak.max_count as f64 >= self.options.min_fraction * strong.min_count as f64
                    || !within_two_edits(&weak.sequence, &strong.sequence)
                {
                    continue;
                }
                // Never remove a path used as evidence for another removal.
                if weak
                    .words
                    .iter()
                    .any(|key| removed.contains(key) || protected.contains(key))
                    || strong.words.iter().any(|key| {
                        let (from, to) = nodes(*key, k);
                        let count = self.edges[key].count as f64;
                        removed.contains(key)
                            || count < self.options.min_fraction * graph[&from].max_out as f64
                            || count < self.options.min_fraction * graph[&to].max_in as f64
                    })
                {
                    continue;
                }
                protected.extend(strong.words);
                removed.extend(weak.words);
                paths += 1;
            }
        }
        (removed, paths)
    }

    fn trace_weak_path(&self, first: u64, graph: &HashMap<u64, Node>) -> Option<Path> {
        let k = self.options.k;
        let mut path = Path {
            words: Vec::new(),
            sequence: Vec::new(),
            min_count: u32::MAX,
            max_count: 0,
            end: 0,
        };
        let mut word = first;
        loop {
            let key = canonical(word, k);
            if path.words.contains(&key) || path.words.len() >= 2 * k + 2 {
                return None;
            }
            let count = self.edges[&key].count;
            path.words.push(key);
            path.sequence.push(b"ACGT"[(word & 3) as usize]);
            path.min_count = path.min_count.min(count);
            path.max_count = path.max_count.max(count);
            path.end = nodes(word, k).1;
            let node = &graph[&path.end];
            if node.incoming != 1 || node.outgoing != 1 {
                return Some(path);
            }
            word = (0..4)
                .map(|base| (path.end << 2) | base)
                .find(|&candidate| self.edges.contains_key(&canonical(candidate, k)))?;
        }
    }

    fn trace_strong_path(&self, start: u64, end: u64, weak: &[u64]) -> Option<Path> {
        let k = self.options.k;
        let mut path = Path {
            words: Vec::new(),
            sequence: Vec::new(),
            min_count: u32::MAX,
            max_count: 0,
            end: start,
        };
        while path.end != end || path.words.is_empty() {
            if path.words.len() >= 2 * k + 2 {
                return None;
            }
            let word = (0..4)
                .map(|base| (path.end << 2) | base)
                .filter_map(|candidate| {
                    let key = canonical(candidate, k);
                    self.edges
                        .get(&key)
                        .map(|edge| (edge.count, candidate, key))
                })
                .max()?;
            if weak.contains(&word.2) || path.words.contains(&word.2) {
                return None;
            }
            path.words.push(word.2);
            path.sequence.push(b"ACGT"[(word.1 & 3) as usize]);
            path.min_count = path.min_count.min(word.0);
            path.max_count = path.max_count.max(word.0);
            path.end = nodes(word.1, k).1;
        }
        Some(path)
    }

    fn walk(&mut self, first: u64, graph: &HashMap<u64, Node>) -> (Vec<u8>, bool) {
        let k = self.options.k;
        let (start, _) = nodes(first, k);
        let mut word = first;
        let mut sequence = Vec::new();
        decode(word, k, &mut sequence);
        self.edges.get_mut(&canonical(word, k)).unwrap().visited = true;
        loop {
            let (_, to) = nodes(word, k);
            let node = &graph[&to];
            if node.incoming != 1 || node.outgoing != 1 {
                return (sequence, false);
            }
            let next = (0..4)
                .map(|base| (to << 2) | base)
                .find(|&candidate| self.edges.contains_key(&canonical(candidate, k)))
                .expect("graph degree must match edges");
            let edge = self.edges.get_mut(&canonical(next, k)).unwrap();
            if edge.visited {
                return (sequence, to == start);
            }
            edge.visited = true;
            sequence.push(b"ACGT"[(next & 3) as usize]);
            word = next;
        }
    }
}

fn read_line<R: BufRead>(reader: &mut R, line: &mut Vec<u8>) -> io::Result<usize> {
    line.clear();
    let n = reader.read_until(b'\n', line)?;
    if line.last() == Some(&b'\n') {
        line.pop();
        if line.last() == Some(&b'\r') {
            line.pop();
        }
    }
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packed_reverse_complement_matches_basewise_version() {
        let mut state = 0x1234_5678_9abc_def0u64;
        for k in (15..=31).step_by(2) {
            for _ in 0..1000 {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                let word = state & ((1u64 << (2 * k)) - 1);
                let mut slow = 0;
                let mut remaining = word;
                for _ in 0..k {
                    slow = (slow << 2) | (3 - (remaining & 3));
                    remaining >>= 2;
                }
                assert_eq!(reverse_complement(word, k), slow);
            }
        }
    }
}
