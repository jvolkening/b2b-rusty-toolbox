use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use viral_overlap::{Options, assemble_with_stats, auto_min_overlap, parse_reads, write_fasta};

const USAGE: &str = "\
viral-overlap - experimental overlap-based viral assembler

USAGE:
    viral-overlap [OPTIONS] READS...

OPTIONS:
    -t, --threads <N>       Worker threads [default: available CPUs, max 8]
        --min-overlap <N>   Minimum accepted read overlap [default: N90 / 3, 100-500]
        --min-identity <F>  Fixed overlap identity [default: adaptive 0.80-0.85]
        --max-secondary <N> Maximum minimap2 overlap hits per read [default: 20]
        --no-polish         Skip read-backed consensus (debugging only)
    -h, --help              Print this help

Reads uncompressed or gzip-compressed FASTA/FASTQ. Outputs polished contigs
as FASTA. This experimental prototype uses minimap2 all-vs-all overlaps,
a greedy non-branching read layout, and read-backed majority consensus.
";

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let mut options = Options::default();
    let mut min_overlap = None;
    let mut paths = Vec::new();
    let mut arguments = std::env::args().skip(1);
    while let Some(argument) = arguments.next() {
        match argument.as_str() {
            "-h" | "--help" => {
                print!("{USAGE}");
                return Ok(());
            }
            "-t" | "--threads" => {
                options.threads = parse_value(arguments.next(), "--threads")?;
            }
            "--min-overlap" => {
                min_overlap = Some(parse_value(arguments.next(), "--min-overlap")?);
            }
            "--min-identity" => {
                options.min_identity = parse_value(arguments.next(), "--min-identity")?;
                options.adaptive_identity = false;
            }
            "--max-secondary" => {
                options.max_secondary = parse_value(arguments.next(), "--max-secondary")?;
            }
            "--no-polish" => options.polish = false,
            "-" => return Err("stdin input is not supported yet; provide a read file".into()),
            _ if argument.starts_with('-') => {
                return Err(format!("unknown option: {argument}").into());
            }
            _ => paths.push(PathBuf::from(argument)),
        }
    }
    if paths.is_empty() {
        return Err("provide at least one FASTA/FASTQ read file".into());
    }
    if options.threads == 0 || options.max_secondary < 1 {
        return Err("threads and max-secondary must be positive".into());
    }
    if !(0.0..=1.0).contains(&options.min_identity) || !options.min_identity.is_finite() {
        return Err("min-identity must be between 0 and 1".into());
    }
    let reads = parse_reads(&paths)?;
    if reads.is_empty() {
        return Err("input contains no reads".into());
    }
    options.min_overlap = min_overlap.unwrap_or_else(|| auto_min_overlap(&reads));
    let identity = if options.adaptive_identity {
        "adaptive".to_owned()
    } else {
        format!("{:.3}", options.min_identity)
    };
    eprintln!(
        "reads={} bases={} min_overlap={} identity={identity}",
        reads.len(),
        reads.iter().map(|read| read.sequence.len()).sum::<usize>(),
        options.min_overlap,
    );
    let assembly = assemble_with_stats(&reads, options)?;
    eprintln!(
        "normalization retained {}/{} reads",
        assembly.retained_reads, assembly.input_reads
    );
    let stdout = io::stdout();
    let mut writer = BufWriter::new(stdout.lock());
    write_fasta(&assembly.contigs, &mut writer)?;
    writer.flush()?;
    eprintln!("contigs={}", assembly.contigs.len());
    Ok(())
}

fn parse_value<T: std::str::FromStr>(
    value: Option<String>,
    name: &str,
) -> Result<T, Box<dyn std::error::Error>> {
    value
        .ok_or_else(|| format!("missing value for {name}"))?
        .parse()
        .map_err(|_| format!("invalid value for {name}").into())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("viral-overlap: {error}");
        std::process::exit(2);
    }
}
