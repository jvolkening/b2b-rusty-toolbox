use krakenuniq_filter::{Denominator, find_taxdb, read_taxonomy};
use std::env;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::process;

const USAGE: &str = "\
krakenuniq-filter - filter KrakenUniq classifications by hit proportion

USAGE:
    krakenuniq-filter [--db TAXDB] [--threshold NUM] [--denominator MODE] [FILE ...]

OPTIONS:
    --threshold <NUM>       Threshold between 0 and 1 [default: 0]
    --db <TAXDB>            Direct path to the KrakenUniq taxDB file
    --denominator <MODE>    all, unambiguous, or assigned [default: unambiguous]
    -h, --help              Print this help message
    --version               Print version information
";

fn die(msg: &str, code: i32) -> ! {
    eprintln!("{msg}\n\n{USAGE}");
    process::exit(code);
}

fn main() -> io::Result<()> {
    let mut threshold = 0.0f64;
    let mut db = None;
    let mut denominator = Denominator::Unambiguous;
    let mut files = Vec::new();
    let mut args = env::args().skip(1);

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-h" | "--help" => {
                print!("{USAGE}");
                return Ok(());
            }
            "--version" => {
                println!("KrakenUniq version 0.1.0");
                return Ok(());
            }
            "--threshold" => {
                let value = args
                    .next()
                    .unwrap_or_else(|| die("missing value for --threshold", 64));
                threshold = value
                    .parse()
                    .unwrap_or_else(|_| die("threshold must be in the interval [0,1].", 64));
            }
            "--db" => {
                db = Some(
                    args.next()
                        .unwrap_or_else(|| die("missing value for --db", 64)),
                );
            }
            "--denominator" => {
                let value = args
                    .next()
                    .unwrap_or_else(|| die("missing value for --denominator", 64));
                denominator = Denominator::parse(&value).unwrap_or_else(|| {
                    die("Denominator must be one of: all unambiguous assigned.", 64)
                });
            }
            value if value.starts_with('-') => {
                die(&format!("unknown option {value}"), 64);
            }
            value => files.push(value.to_string()),
        }
    }

    if !(0.0..=1.0).contains(&threshold) {
        die("threshold must be in the interval [0,1].", 64);
    }

    let taxdb = match find_taxdb(db.as_deref()) {
        Ok(path) => path,
        Err(message) => {
            eprintln!("{message}");
            process::exit(64);
        }
    };
    eprint!("Reading taxonomy ...");
    let taxonomy = match read_taxonomy(&taxdb) {
        Ok(taxonomy) => taxonomy,
        Err(error) => {
            eprintln!(" {error}");
            process::exit(1);
        }
    };
    eprintln!(" done.");

    let stdout = io::stdout();
    let mut writer = BufWriter::with_capacity(256 * 1024, stdout.lock());

    if files.is_empty() {
        let stdin = io::stdin();
        krakenuniq_filter::filter(
            BufReader::new(stdin.lock()),
            &mut writer,
            &taxonomy,
            threshold,
            denominator,
        )?;
    } else {
        for file in files {
            if file == "-" {
                let stdin = io::stdin();
                krakenuniq_filter::filter(
                    BufReader::new(stdin.lock()),
                    &mut writer,
                    &taxonomy,
                    threshold,
                    denominator,
                )?;
            } else {
                let input = match File::open(&file) {
                    Ok(input) => input,
                    Err(error) => {
                        eprintln!("Error: cannot open {file}: {error}");
                        process::exit(1);
                    }
                };
                krakenuniq_filter::filter(
                    BufReader::new(input),
                    &mut writer,
                    &taxonomy,
                    threshold,
                    denominator,
                )?;
            }
        }
    }
    writer.flush()
}
