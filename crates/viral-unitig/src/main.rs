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
    --min-fraction <F>       Minimum local branch support, 0 <= F < 1 [default: 0.2]
    --clean-bubbles          Experimentally remove weak, similar bubble paths
    --no-clean-errors        Disable experimental near-neighbor error filtering
    --clean-errors           Enable near-neighbor error filtering [default]
    --uncertain <FILE>       Save short, possibly redundant unitigs separately
    --min-length <N>         Minimum emitted unitig length [default: k]
    --max-kmers <N>          Maximum distinct input k-mers [default: 10000000]
    -h, --help               Print help

Research prototype: emits unphased unitigs, not complete viral genomes.
It is not validated for clinical or diagnostic use.
";

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let mut options = Options::default();
    let mut inputs = Vec::new();
    let mut uncertain_path = None;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-h" | "--help" => {
                print!("{HELP}");
                return Ok(());
            }
            "--clean-bubbles" => options.clean_bubbles = true,
            "--clean-errors" => options.clean_errors = true,
            "--no-clean-errors" => options.clean_errors = false,
            "-k" | "--kmer" | "--min-count" | "--min-fraction" | "--min-length" | "--max-kmers"
            | "--uncertain" => {
                let value = args
                    .next()
                    .ok_or_else(|| format!("{arg} requires a value"))?;
                match arg.as_str() {
                    "-k" | "--kmer" => options.k = value.parse()?,
                    "--min-count" => options.min_count = value.parse()?,
                    "--min-fraction" => options.min_fraction = value.parse()?,
                    "--min-length" => options.min_length = value.parse()?,
                    "--max-kmers" => options.max_kmers = value.parse()?,
                    "--uncertain" => {
                        if uncertain_path.replace(value).is_some() {
                            return Err("--uncertain may only be specified once".into());
                        }
                    }
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
    if uncertain_path.is_some() && !options.clean_errors {
        return Err("--uncertain requires --clean-errors".into());
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
    let mut uncertain = uncertain_path
        .map(|path| {
            File::create_new(&path)
                .map(BufWriter::new)
                .map_err(|err| io::Error::new(err.kind(), format!("cannot create {path}: {err}")))
        })
        .transpose()?;
    let stats = assembler.write_unitigs_with_uncertain(
        BufWriter::with_capacity(256 * 1024, stdout.lock()),
        uncertain.as_mut().map(|file| file as &mut dyn io::Write),
    )?;
    eprintln!(
        "processed {} reads, {} bases; wrote {} unitigs; removed {} error k-mers and {} weak bubble paths ({} k-mers); separated {} uncertain unitigs",
        reads,
        bases,
        stats.unitigs,
        stats.removed_error_kmers,
        stats.removed_bubble_paths,
        stats.removed_bubble_kmers,
        stats.separated_uncertain_unitigs
    );
    if stats.redundancy_search_limit > 0 {
        eprintln!(
            "warning: retained {} unitigs whose redundancy search hit its limit",
            stats.redundancy_search_limit
        );
    }
    if stats.discarded_uncertain_unitigs > 0 {
        eprintln!(
            "warning: discarded {} short, possibly redundant unitigs; use --uncertain FILE to save them",
            stats.discarded_uncertain_unitigs
        );
    }
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
