use std::collections::{HashMap, hash_map::Entry};
use std::io::{self, BufRead, Write};

#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub k: usize,
    pub min_count: u32,
    pub min_fraction: f64,
    pub min_length: usize,
    pub max_kmers: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            k: 21,
            min_count: 3,
            min_fraction: 0.1,
            min_length: 0,
            max_kmers: 10_000_000,
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

pub struct Assembler {
    options: Options,
    edges: HashMap<u64, Edge>,
    pub reads: u64,
    pub bases: u64,
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

fn reverse_complement(mut word: u64, k: usize) -> u64 {
    let mut result = 0;
    for _ in 0..k {
        result = (result << 2) | (3 - (word & 3));
        word >>= 2;
    }
    result
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

    pub fn write_unitigs<W: Write>(mut self, mut writer: W) -> io::Result<usize> {
        let k = self.options.k;
        self.edges
            .retain(|_, edge| edge.count >= self.options.min_count);

        let mut graph: HashMap<u64, Node> = HashMap::new();
        for (&word, edge) in &self.edges {
            for oriented in [word, reverse_complement(word, k)] {
                let (from, to) = nodes(oriented, k);
                let source = graph.entry(from).or_default();
                source.max_out = source.max_out.max(edge.count);
                let target = graph.entry(to).or_default();
                target.max_in = target.max_in.max(edge.count);
            }
        }

        let fraction = self.options.min_fraction;
        self.edges.retain(|&word, edge| {
            let (from, to) = nodes(word, k);
            (edge.count as f64) >= fraction * (graph[&from].max_out as f64)
                && (edge.count as f64) >= fraction * (graph[&to].max_in as f64)
        });
        for &word in self.edges.keys() {
            for oriented in [word, reverse_complement(word, k)] {
                let (from, to) = nodes(oriented, k);
                graph.get_mut(&from).unwrap().outgoing += 1;
                graph.get_mut(&to).unwrap().incoming += 1;
            }
        }
        let mut words: Vec<u64> = self.edges.keys().copied().collect();
        words.sort_unstable();
        let mut emitted = 0;

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
                    emitted += 1;
                    writeln!(
                        writer,
                        ">unitig_{emitted} length={} circular={circular}",
                        sequence.len()
                    )?;
                    for chunk in sequence.chunks(80) {
                        writer.write_all(chunk)?;
                        writer.write_all(b"\n")?;
                    }
                }
            }
        }
        writer.flush()?;
        Ok(emitted)
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
