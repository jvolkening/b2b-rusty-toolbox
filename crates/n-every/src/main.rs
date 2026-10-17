use std::env;
use std::io::{self, BufWriter};
use std::process;

const USAGE: &str = "\
n-every - insert an 'N' after every K sequence characters in FASTA input

USAGE:
    n-every <K> < input.fa > output.fa

ARGS:
    <k>    Positive integer; number of bases between inserted N's

OPTIONS:
    -h, --help    Print this help message
";

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().skip(1).collect();

    if args.iter().any(|a| a == "-h" || a == "--help") {
        print!("{USAGE}");
        return Ok(());
    }

    if args.len() != 1 {
        eprintln!("Error: Expected exactly one argument <k>.\n\n{USAGE}");
        process::exit(2);
    }

    let k: usize = match args[0].parse::<usize>() {
        Ok(val) if val > 0 => val,
        _ => {
            eprintln!("Error: 'k' must be a positive integer greater than 0.");
            process::exit(2);
        }
    };

    let stdin = io::stdin();
    let stdout = io::stdout();

    let reader = stdin.lock();
    // Large 256 KB internal I/O cache to maximize pipe throughput
    let writer = BufWriter::with_capacity(256 * 1024, stdout.lock());
    n_every::insert_n(reader, writer, k)
}
