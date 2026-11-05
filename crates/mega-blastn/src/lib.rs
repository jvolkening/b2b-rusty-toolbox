pub mod asn1;
pub mod blastdb;
pub mod extend;
pub mod fasta;
pub mod index;
pub mod output;
pub mod search;
pub mod stats;

use std::io::{BufRead, BufReader, BufWriter, Write};

pub const USAGE: &str = "\
mega-blastn: megablast-style nucleotide search against BLAST v4/v5 databases

USAGE:
    mega-blastn -db <name> [-query <fasta>|-] [options]

OPTIONS:
    -db <name>             BLAST nucleotide database (alias/multi-volume OK; $BLASTDB is searched)
    -query <file>          query FASTA (default: stdin)
    -out <file>            output file (default: stdout)
    -outfmt <fmt>          \"6 [fields...]\" or \"7 [fields...]\" (default: 6 std)
    -evalue <x>            E-value threshold (default: 10)
    -word_size <n>         exact-match seed length, >= 11 (default: 28)
    -max_target_seqs <n>   subjects reported per query (default: 500)
    -max_hsps <n>          HSPs reported per subject (default: unlimited)
    -perc_identity <x>     minimum percent identity (default: 0)
    -dust <yes|no>         DUST-mask queries for seeding (default: yes)
    -strand <both|plus|minus>
                           query strand(s) to search (default: both)
    -reward <1> -penalty <-2>
                           accepted only with the megablast defaults
    -num_threads <n>       worker threads (default: all cores)
    -query_batch_size <n>  query residues per database pass (default: 10000000)
    -h, -help              show this help
";

struct Opts {
    db: String,
    query: Option<String>,
    out: Option<String>,
    outfmt: String,
    evalue: f64,
    word_size: usize,
    max_target_seqs: usize,
    max_hsps: Option<usize>,
    perc_identity: f64,
    dust: bool,
    strand: String,
    threads: usize,
    batch: usize,
}

fn parse_args(args: &[String]) -> Result<Option<Opts>, String> {
    let mut o = Opts {
        db: String::new(),
        query: None,
        out: None,
        outfmt: "6".into(),
        evalue: 10.0,
        word_size: 28,
        max_target_seqs: 500,
        max_hsps: None,
        perc_identity: 0.0,
        dust: true,
        strand: "both".into(),
        threads: std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1),
        batch: 10_000_000,
    };
    let mut i = 0;
    while i < args.len() {
        let a = args[i].as_str();
        if a == "-h" || a == "-help" || a == "--help" {
            return Ok(None);
        }
        let mut val = || -> Result<String, String> {
            i += 1;
            args.get(i).cloned().ok_or_else(|| format!("missing value for {a}"))
        };
        fn num<T: std::str::FromStr>(s: &str, a: &str) -> Result<T, String> {
            s.parse().map_err(|_| format!("invalid value '{s}' for {a}"))
        }
        match a {
            "-db" => o.db = val()?,
            "-query" => o.query = Some(val()?),
            "-out" => o.out = Some(val()?),
            "-outfmt" => o.outfmt = val()?,
            "-evalue" => o.evalue = num(&val()?, a)?,
            "-word_size" => o.word_size = num(&val()?, a)?,
            "-max_target_seqs" => o.max_target_seqs = num(&val()?, a)?,
            "-max_hsps" => o.max_hsps = Some(num(&val()?, a)?),
            "-perc_identity" => o.perc_identity = num(&val()?, a)?,
            "-num_threads" => o.threads = num(&val()?, a)?,
            "-query_batch_size" => o.batch = num(&val()?, a)?,
            "-strand" => o.strand = val()?,
            "-dust" => {
                o.dust = match val()?.as_str() {
                    "yes" | "true" | "1" => true,
                    "no" | "false" | "0" => false,
                    v => return Err(format!("invalid value '{v}' for -dust")),
                }
            }
            "-reward" => {
                if num::<i32>(&val()?, a)? != 1 {
                    return Err("only -reward 1 is supported".into());
                }
            }
            "-penalty" => {
                if num::<i32>(&val()?, a)? != -2 {
                    return Err("only -penalty -2 is supported".into());
                }
            }
            "-task" => {
                if val()? != "megablast" {
                    return Err("only -task megablast is supported".into());
                }
            }
            _ => return Err(format!("unknown option '{a}'")),
        }
        i += 1;
    }
    if o.db.is_empty() {
        return Err("-db is required".into());
    }
    if o.word_size < 11 {
        return Err("-word_size must be at least 11".into());
    }
    if !matches!(o.strand.as_str(), "both" | "plus" | "minus") {
        return Err("-strand must be both, plus or minus".into());
    }
    Ok(Some(o))
}

pub fn run(args: &[String]) -> Result<(), String> {
    let Some(o) = parse_args(args)? else {
        print!("{USAGE}");
        return Ok(());
    };
    let (fmt, fields) = output::parse_outfmt(&o.outfmt)?;
    let db = blastdb::Database::open(&o.db).map_err(|e| e.to_string())?;
    let input: Box<dyn BufRead> = match o.query.as_deref() {
        None | Some("-") => Box::new(BufReader::new(std::io::stdin())),
        Some(p) => Box::new(BufReader::new(
            std::fs::File::open(p).map_err(|e| format!("{p}: {e}"))?,
        )),
    };
    let sink: Box<dyn Write> = match o.out.as_deref() {
        None | Some("-") => Box::new(std::io::stdout().lock()),
        Some(p) => Box::new(std::fs::File::create(p).map_err(|e| format!("{p}: {e}"))?),
    };
    let mut out = BufWriter::new(sink);
    let params = search::Params {
        word_size: o.word_size,
        evalue: o.evalue,
        max_target_seqs: o.max_target_seqs,
        max_hsps: o.max_hsps,
        perc_identity: o.perc_identity,
        threads: o.threads,
        plus: o.strand != "minus",
        minus: o.strand != "plus",
    };
    let printer = output::Printer { fields, fmt7: fmt == 7, db_arg: &o.db, db: &db };
    let mut reader = fasta::QueryReader::new(input, o.dust);
    let mut nq = 0usize;
    let io = |e: std::io::Error| e.to_string();
    loop {
        let queries = reader.next_batch(o.batch).map_err(io)?;
        if queries.is_empty() {
            break;
        }
        let ctxs = search::make_contexts(&queries, &db);
        let nkmers = 2 * queries.iter().map(|q| q.len()).sum::<usize>();
        let (k, stride) = index::choose_seed(o.word_size, nkmers);
        let index = index::Index::build(&queries, k, stride);
        let hits = if index.is_empty() {
            vec![Vec::new(); queries.len()]
        } else {
            search::search(&db, &queries, &ctxs, &index, &params)
        };
        drop(index);
        for (q, h) in queries.iter().zip(&hits) {
            printer.print_query(&mut out, q, h).map_err(io)?;
        }
        nq += queries.len();
    }
    if fmt == 7 {
        writeln!(out, "# BLAST processed {nq} queries").map_err(io)?;
    }
    out.flush().map_err(io)
}
