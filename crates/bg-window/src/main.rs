use bg_window::{Stat, read_bedgraph, read_fasta_lengths, run};
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::process;

const USAGE: &str = "\
bg-window - shrink a bedgraph with per-sequence adaptive sliding windows

USAGE:
    bg-window -b <in.bedgraph> -f <ref.fa> --min-entries <N> --max-entries <M> [OPTIONS]

OPTIONS:
    -b, --bedgraph <FILE>    Input bedgraph (uncovered bases are treated as depth 0)
    -f, --fasta <FILE>       Matching FASTA reference (defines sequence lengths)
    -n, --min-entries <N>    Minimum entries per sequence in output
    -N, --max-entries <M>    Maximum entries per sequence in output
    -s, --stat <STAT>        min, max, median, mean, sum [default: mean]
        --no-merge           Do not merge adjacent windows with identical values
    -o, --output <FILE>      Output file [default: stdout]
    -h, --help               Print this help message

Each sequence gets its own window size, tuned so the entry count falls within
[min, max]. Sequences too short to reach the minimum are emitted at the finest
resolution possible (1 bp windows).
";

fn die(msg: &str) -> ! {
    eprintln!("Error: {msg}\n\n{USAGE}");
    process::exit(2);
}

fn open(p: &str) -> BufReader<File> {
    match File::open(p) {
        Ok(f) => BufReader::new(f),
        Err(e) => {
            eprintln!("Error: cannot open {p}: {e}");
            process::exit(1);
        }
    }
}

fn main() -> io::Result<()> {
    let mut args = std::env::args().skip(1);
    let (mut bg, mut fa, mut min, mut max, mut out) = (None, None, None, None, None);
    let mut stat = Stat::Mean;
    let mut merge = true;
    while let Some(a) = args.next() {
        if a == "-h" || a == "--help" {
            print!("{USAGE}");
            return Ok(());
        }
        if a == "--no-merge" {
            merge = false;
            continue;
        }
        let v = args.next().unwrap_or_else(|| die(&format!("missing value for {a}")));
        let num = |v: &str| v.parse::<usize>().ok().filter(|n| *n > 0)
            .unwrap_or_else(|| die(&format!("{a} must be a positive integer")));
        match a.as_str() {
            "-b" | "--bedgraph" => bg = Some(v),
            "-f" | "--fasta" => fa = Some(v),
            "-n" | "--min-entries" => min = Some(num(&v)),
            "-N" | "--max-entries" => max = Some(num(&v)),
            "-s" | "--stat" => {
                stat = Stat::parse(&v).unwrap_or_else(|| die("unknown --stat value"))
            }
            "-o" | "--output" => out = Some(v),
            _ => die(&format!("unknown option {a}")),
        }
    }
    let (Some(bg), Some(fa), Some(min), Some(max)) = (bg, fa, min, max) else {
        die("--bedgraph, --fasta, --min-entries and --max-entries are required");
    };
    if min > max {
        die("--min-entries must not exceed --max-entries");
    }
    let fasta = read_fasta_lengths(open(&fa))?;
    let bedgraph = read_bedgraph(open(&bg))?;
    let w: Box<dyn Write> = match out {
        Some(p) => Box::new(File::create(p)?),
        None => Box::new(io::stdout().lock()),
    };
    let mut w = BufWriter::with_capacity(256 * 1024, w);
    if let Err(e) = run(&fasta, &bedgraph, min, max, stat, merge, &mut w) {
        eprintln!("Error: {e}");
        process::exit(1);
    }
    w.flush()
}
