use ir_chimera::{Counts, Options, Result, filter};
use std::fs::{File, OpenOptions};
use std::io::{self, BufReader, BufWriter, Write};
use std::process;

const USAGE: &str = "\
ir-chimera - filter nanopore reads with inverted-repeat self-alignments

USAGE:
    ir-chimera [OPTIONS] [FILE ...] > output.fastx

INPUT:
    FASTA or four-line FASTQ, plain or gzip-compressed. All files must have
    the same format. Reads stdin when no files are given; use - for stdin.

OPTIONS:
    -c, --chimeric          Output chimeric rather than non-chimeric reads
    -l, --min-length <N>    Minimum aligned bases in each arm [default: 500]
    -i, --min-identity <P>  Minimum percent identity (0-100) [default: 80]
    -t, --threads <N>       Worker threads [default: available CPUs, max 8]
        --log <FILE>        Write summary TSV to a new file (must not exist)
    -h, --help              Print this help

Uses minimap2's ava-ont preset in-process, with reverse-strand alignments only.
Reads without a non-overlapping inverted 15-mer pair are kept without alignment.
Identity is exact A/C/G/T matches divided by alignment columns including gaps.
The two aligned arms must not overlap. This is a heuristic, not a proof of
chimera formation; repeats and sequencing errors can affect classification.
";

fn run() -> Result<()> {
    let mut options = Options::default();
    let mut files = Vec::new();
    let mut log = None;
    let mut threads = std::thread::available_parallelism().map_or(1, |n| n.get().min(8));
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-h" | "--help" => {
                print!("{USAGE}");
                return Ok(());
            }
            "-c" | "--chimeric" => options.chimeric = true,
            "--log" => {
                if log.is_some() {
                    return Err("--log may only be specified once".into());
                }
                let path = args.next().ok_or("missing value for --log")?;
                if path.starts_with('-') {
                    return Err("invalid --log path; use ./ for names beginning with '-'".into());
                }
                log = Some(path);
            }
            "-l" | "--min-length" => {
                let value = args.next().ok_or("missing value for --min-length")?;
                options.min_length = value.parse().map_err(|_| "invalid --min-length")?;
                if options.min_length == 0 {
                    return Err("--min-length must be positive".into());
                }
            }
            "-i" | "--min-identity" => {
                let value = args.next().ok_or("missing value for --min-identity")?;
                let percent: f64 = value.parse().map_err(|_| "invalid --min-identity")?;
                if !percent.is_finite() || !(0.0..=100.0).contains(&percent) {
                    return Err("--min-identity must be between 0 and 100".into());
                }
                options.min_identity = percent / 100.0;
            }
            "-t" | "--threads" => {
                let value = args.next().ok_or("missing value for --threads")?;
                threads = value.parse().map_err(|_| "invalid --threads")?;
                if threads == 0 {
                    return Err("--threads must be positive".into());
                }
            }
            "-" => files.push(arg),
            _ if arg.starts_with('-') => return Err(format!("unknown option: {arg}").into()),
            _ => files.push(arg),
        }
    }
    if files.is_empty() {
        files.push("-".to_owned());
    }
    if files.iter().filter(|f| f.as_str() == "-").count() > 1 {
        return Err("stdin may only be specified once".into());
    }
    let stdout = io::stdout();
    let mut writer = BufWriter::with_capacity(256 * 1024, stdout.lock());
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()?;
    let mut format = None;
    let mut counts = Counts::default();
    for file in files {
        if file == "-" {
            filter(
                BufReader::new(io::stdin()),
                &mut writer,
                options,
                &pool,
                &mut format,
                &mut counts,
            )
            .map_err(|e| format!("stdin: {e}"))?;
        } else {
            let input = File::open(&file).map_err(|e| format!("{file}: {e}"))?;
            filter(
                BufReader::new(input),
                &mut writer,
                options,
                &pool,
                &mut format,
                &mut counts,
            )
            .map_err(|e| format!("{file}: {e}"))?;
        }
    }
    writer.flush()?;
    if let Some(path) = log {
        let fraction = if counts.n_total == 0 {
            0.0
        } else {
            counts.n_kept as f64 / counts.n_total as f64
        };
        let mut file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .map_err(|e| format!("cannot create log {path}: {e}"))?;
        write!(
            file,
            "n_total\tn_kept\tfraction\n{}\t{}\t{fraction:?}\n",
            counts.n_total, counts.n_kept
        )
        .map_err(|e| format!("cannot write log {path}: {e}"))?;
        file.flush()
            .map_err(|e| format!("cannot flush log {path}: {e}"))?;
    }
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("Error: {e}");
        process::exit(1);
    }
}
