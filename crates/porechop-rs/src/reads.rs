use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use walkdir::WalkDir;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Fasta,
    Fastq,
}

#[derive(Clone)]
pub struct ReadRecord {
    pub name: String,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub rna: bool,
    pub albacore: Option<String>,
    pub start_trim: usize,
    pub end_trim: usize,
    pub middle: Vec<(usize, usize)>,
    pub barcode: String,
    pub start_hits: Vec<(String, crate::align::Alignment)>,
    pub end_hits: Vec<(String, crate::align::Alignment)>,
    pub middle_hits: Vec<(String, crate::align::Alignment)>,
    pub start_scores: Vec<(String, f64)>,
    pub end_scores: Vec<(String, f64)>,
}

impl ReadRecord {
    fn new(name: String, seq: String, qual: String, albacore: Option<String>) -> Self {
        let mut seq = seq.to_ascii_uppercase().into_bytes();
        let rna =
            seq.iter().filter(|&&b| b == b'U').count() > seq.iter().filter(|&&b| b == b'T').count();
        if rna {
            for b in &mut seq {
                if *b == b'U' {
                    *b = b'T';
                }
            }
        }
        let mut qual = qual.into_bytes();
        qual.resize(qual.len().max(seq.len()), b'+');
        Self {
            name,
            seq,
            qual,
            rna,
            albacore,
            start_trim: 0,
            end_trim: 0,
            middle: Vec::new(),
            barcode: "none".into(),
            start_hits: Vec::new(),
            end_hits: Vec::new(),
            middle_hits: Vec::new(),
            start_scores: Vec::new(),
            end_scores: Vec::new(),
        }
    }

    pub fn trimmed(&self) -> &[u8] {
        let end = python_end(self.seq.len(), self.end_trim);
        &self.seq[self.start_trim.min(end)..end]
    }

    fn render(
        &self,
        format: Format,
        min_size: usize,
        discard_middle: bool,
        untrimmed: bool,
    ) -> Vec<u8> {
        if !self.middle.is_empty() && discard_middle {
            return Vec::new();
        }
        let end = python_end(self.seq.len(), self.end_trim);
        let start = self.start_trim.min(end);
        let (seq, qual) = if untrimmed && self.middle.is_empty() {
            (&self.seq[..], &self.qual[..])
        } else {
            let qend = python_end(self.qual.len(), self.end_trim);
            (&self.seq[start..end], &self.qual[start.min(qend)..qend])
        };
        let mut ranges = Vec::new();
        if self.middle.is_empty() {
            if !seq.is_empty() {
                ranges.push((0, seq.len()));
            }
        } else {
            let mut mask = vec![false; seq.len()];
            for &(a, b) in &self.middle {
                mask[a.min(seq.len())..b.min(seq.len())].fill(true);
            }
            let mut i = 0;
            while i < seq.len() {
                while i < seq.len() && mask[i] {
                    i += 1;
                }
                let begin = i;
                while i < seq.len() && !mask[i] {
                    i += 1;
                }
                if i - begin >= min_size && i > begin {
                    ranges.push((begin, i));
                }
            }
        }
        let mut out = Vec::new();
        for (index, (a, b)) in ranges.into_iter().enumerate() {
            let name = if self.middle.is_empty() {
                self.name.clone()
            } else if let Some((first, rest)) = self.name.split_once(' ') {
                format!("{first}_{} {rest}", index + 1)
            } else {
                format!("{}_{}", self.name, index + 1)
            };
            out.push(if format == Format::Fasta { b'>' } else { b'@' });
            out.extend_from_slice(name.as_bytes());
            out.push(b'\n');
            let mut part = seq[a..b].to_vec();
            if self.rna {
                for base in &mut part {
                    if *base == b'T' {
                        *base = b'U';
                    }
                }
            }
            if format == Format::Fasta {
                for chunk in part.chunks(70) {
                    out.extend_from_slice(chunk);
                    out.push(b'\n');
                }
            } else {
                out.extend_from_slice(&part);
                out.extend_from_slice(b"\n+\n");
                if self.middle.is_empty() {
                    out.extend_from_slice(qual);
                } else {
                    out.extend_from_slice(&qual[a..b]);
                }
                out.push(b'\n');
            }
        }
        out
    }
}

fn python_end(length: usize, trim: usize) -> usize {
    if trim <= length {
        length - trim
    } else {
        length.saturating_sub(trim - length)
    }
}

fn input_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut magic = [0; 2];
    let n = file.read(&mut magic)?;
    drop(file);
    let file = File::open(path)?;
    if n == 2 && magic == [0x1f, 0x8b] {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn load_file(path: &Path, albacore: Option<String>) -> Result<(Vec<ReadRecord>, Format)> {
    let reader = input_reader(path)?;
    let mut lines = reader.lines();
    let first = lines.next().context("empty sequence file")??;
    let mut reads = Vec::new();
    let format = if let Some(first_name) = first.strip_prefix('>') {
        let mut name = first_name.trim().to_owned();
        let mut seq = String::new();
        for line in lines {
            let line = line?;
            let line = line.trim();
            if let Some(header) = line.strip_prefix('>') {
                if !name.is_empty() {
                    reads.push(ReadRecord::new(name, seq, String::new(), albacore.clone()));
                }
                name = header.to_owned();
                seq = String::new();
            } else {
                seq.push_str(line);
            }
        }
        if !name.is_empty() {
            reads.push(ReadRecord::new(name, seq, String::new(), albacore));
        }
        Format::Fasta
    } else if first.starts_with('@') {
        let mut header = first;
        loop {
            if !header.starts_with('@') {
                bail!("invalid FASTQ header in {}", path.display());
            }
            let seq = lines.next().context("truncated FASTQ sequence")??;
            let plus = lines.next().context("truncated FASTQ separator")??;
            let qual = lines.next().context("truncated FASTQ qualities")??;
            if !plus.starts_with('+') {
                bail!("invalid FASTQ separator in {}", path.display());
            }
            reads.push(ReadRecord::new(
                header[1..].trim().to_owned(),
                seq.trim().to_owned(),
                qual.trim().to_owned(),
                albacore.clone(),
            ));
            match lines.next() {
                Some(line) => header = line?,
                None => break,
            }
        }
        Format::Fastq
    } else {
        bail!("{} is neither FASTA nor FASTQ", path.display());
    };
    Ok((reads, format))
}

fn albacore(path: &Path) -> Option<String> {
    let mut barcode = None;
    for component in path.components() {
        let text = component.as_os_str().to_string_lossy();
        if text == "unclassified" {
            return Some("none".into());
        }
        if let Some(suffix) = text.strip_prefix("barcode")
            && suffix.len() == 2
            && suffix.bytes().all(|c| c.is_ascii_digit())
        {
            barcode = Some(format!("BC{suffix}"));
        }
    }
    barcode
}

pub fn load(input: &Path, check_count: usize) -> Result<(Vec<ReadRecord>, Vec<usize>, Format)> {
    if input.is_file() {
        let (reads, format) = load_file(input, None)?;
        let check = (0..reads.len().min(check_count)).collect();
        return Ok((reads, check, format));
    }
    if !input.is_dir() {
        bail!("could not find {}", input.display());
    }
    let mut paths: Vec<PathBuf> = WalkDir::new(input)
        .into_iter()
        .map(|entry| entry.map(|e| e.into_path()))
        .collect::<std::result::Result<_, _>>()?;
    paths.retain(|p| {
        p.is_file() && {
            let p = p.to_string_lossy().to_lowercase();
            p.ends_with(".fastq") || p.ends_with(".fastq.gz")
        }
    });
    paths.sort();
    if paths.is_empty() {
        bail!("could not find fastq files in {}", input.display());
    }
    let per_file = (check_count as f64 / paths.len() as f64).round_ties_even() as usize;
    let mut reads = Vec::new();
    let mut check = Vec::new();
    for path in paths {
        let (batch, format) = load_file(&path, albacore(&path))?;
        if format != Format::Fastq {
            bail!("{} is not FASTQ", path.display());
        }
        check.extend(reads.len()..reads.len() + batch.len().min(per_file));
        reads.extend(batch);
    }
    Ok((reads, check, Format::Fastq))
}

enum Output {
    Plain(BufWriter<File>),
    Gzip(GzEncoder<BufWriter<File>>),
    Stdout(BufWriter<io::Stdout>),
}

impl Write for Output {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        match self {
            Self::Plain(w) => w.write(data),
            Self::Gzip(w) => w.write(data),
            Self::Stdout(w) => w.write(data),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Plain(w) => w.flush(),
            Self::Gzip(w) => w.flush(),
            Self::Stdout(w) => w.flush(),
        }
    }
}

impl Output {
    fn finish(self) -> Result<()> {
        match self {
            Self::Plain(mut w) => w.flush()?,
            Self::Gzip(w) => w.finish()?.flush()?,
            Self::Stdout(mut w) => w.flush()?,
        }
        Ok(())
    }
}

fn writer(path: &Path, gzip: bool) -> Result<Output> {
    let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    if gzip {
        Ok(Output::Gzip(GzEncoder::new(
            BufWriter::new(file),
            Compression::default(),
        )))
    } else {
        Ok(Output::Plain(BufWriter::new(file)))
    }
}

pub struct OutputSettings<'a> {
    pub output: Option<&'a Path>,
    pub barcode_dir: Option<&'a Path>,
    pub format: Format,
    pub gzip: bool,
    pub min_size: usize,
    pub discard_middle: bool,
    pub untrimmed: bool,
    pub discard_unassigned: bool,
}

pub fn output(
    reads: &[ReadRecord],
    settings: OutputSettings<'_>,
) -> Result<BTreeMap<String, (usize, usize)>> {
    let OutputSettings {
        output,
        barcode_dir,
        format,
        gzip,
        min_size,
        discard_middle,
        untrimmed,
        discard_unassigned,
    } = settings;
    let mut counts = BTreeMap::new();
    if let Some(dir) = barcode_dir {
        fs::create_dir_all(dir)?;
        let mut files: BTreeMap<String, Output> = BTreeMap::new();
        for read in reads {
            if discard_unassigned && read.barcode == "none" {
                continue;
            }
            let result = read.render(format, min_size, discard_middle, untrimmed);
            if result.is_empty() {
                continue;
            }
            let count = counts.entry(read.barcode.clone()).or_insert((0, 0));
            count.0 += 1;
            count.1 += if untrimmed {
                read.seq.len()
            } else {
                read.trimmed().len()
            };
            if !files.contains_key(&read.barcode) {
                let ext = if format == Format::Fasta {
                    "fasta"
                } else {
                    "fastq"
                };
                let suffix = if gzip { ".gz" } else { "" };
                files.insert(
                    read.barcode.clone(),
                    writer(&dir.join(format!("{}.{ext}{suffix}", read.barcode)), gzip)?,
                );
            }
            files.get_mut(&read.barcode).unwrap().write_all(&result)?;
        }
        for (_, output) in files {
            output.finish()?;
        }
    } else {
        let mut out = if let Some(path) = output {
            writer(path, gzip)?
        } else {
            Output::Stdout(BufWriter::new(io::stdout()))
        };
        for read in reads {
            out.write_all(&read.render(format, min_size, discard_middle, false))?;
        }
        out.finish()?;
    }
    Ok(counts)
}
