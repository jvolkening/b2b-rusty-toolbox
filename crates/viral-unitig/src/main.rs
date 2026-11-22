use std::fs::File;
use std::io::{self, BufReader, BufWriter};
use std::process::ExitCode;
use viral_unitig::{Assembler, Options};

const HELP: &str = "\
viral-unitig - experimental, coverage-aware metagenomic unitigger

USAGE:
    viral-unitig [OPTIONS] [READS.fa|READS.fq ...] > unitigs.fa

Reads uncompressed FASTA/FASTQ files, or stdin when no files are given.
Use '-' once to include stdin among input files.

OPTIONS:
    -k, --kmer <ODD>         k-mer size, 15 to 31 [default: 21]
    --min-count <N>          Minimum k-mer occurrences [default: 3]
    --min-fraction <F>       Minimum local branch support, 0 <= F < 1 [default: 0.1]
    --clean-bubbles          Experimentally remove weak, similar bubble paths
    --min-length <N>         Minimum emitted unitig length [default: k]
    --max-kmers <N>          Maximum distinct input k-mers [default: 10000000]
    -h, --help               Print help

Research prototype: emits unphased unitigs, not complete viral genomes.
It is not validated for clinical or diagnostic use.
";

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let mut options = Options::default();
    let mut inputs = Vec::new();
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-h" | "--help" => {
                print!("{HELP}");
                return Ok(());
            }
            "--clean-bubbles" => options.clean_bubbles = true,
            "-k" | "--kmer" | "--min-count" | "--min-fraction" | "--min-length" | "--max-kmers" => {
                let value = args
                    .next()
                    .ok_or_else(|| format!("{arg} requires a value"))?;
                match arg.as_str() {
                    "-k" | "--kmer" => options.k = value.parse()?,
                    "--min-count" => options.min_count = value.parse()?,
                    "--min-fraction" => options.min_fraction = value.parse()?,
                    "--min-length" => options.min_length = value.parse()?,
                    "--max-kmers" => options.max_kmers = value.parse()?,
                    _ => unreachable!(),
                }
            }
            "-" => inputs.push(arg),
            _ if arg.starts_with('-') => return Err(format!("unknown option: {arg}").into()),
            _ => inputs.push(arg),
        }
    }
    if inputs.is_empty() {
        inputs.push("-".to_owned());
    }
    if inputs.iter().filter(|name| name.as_str() == "-").count() > 1 {
        return Err("stdin may only be specified once".into());
    }
    let mut assembler = Assembler::new(options)?;
    for name in inputs {
        let reader: Box<dyn io::BufRead> = if name == "-" {
            Box::new(BufReader::with_capacity(256 * 1024, io::stdin()))
        } else {
            Box::new(BufReader::with_capacity(256 * 1024, File::open(&name)?))
        };
        assembler
            .read_fastx(reader)
            .map_err(|err| io::Error::new(err.kind(), format!("{name}: {err}")))?;
    }
    let stdout = io::stdout();
    let (reads, bases) = (assembler.reads, assembler.bases);
    let stats = assembler.write_unitigs(BufWriter::with_capacity(256 * 1024, stdout.lock()))?;
    eprintln!(
        "processed {} reads, {} bases; wrote {} unitigs; removed {} weak bubble paths ({} k-mers)",
        reads, bases, stats.unitigs, stats.removed_bubble_paths, stats.removed_bubble_kmers
    );
    if stats.unitigs == 0 {
        eprintln!("warning: no unitigs passed the filters");
    }
    Ok(())
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("error: {err}");
            ExitCode::FAILURE
        }
    }
}
