mod adapters;
mod align;
mod process;
mod reads;
mod verbose;

use anyhow::{Context, Result, bail};
use clap::Parser;
use process::Settings;
use reads::Format;
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    version,
    about = "Find and trim Oxford Nanopore adapters and demultiplex barcodes"
)]
struct Args {
    /// FASTA/FASTQ input file or directory of FASTQ files
    #[arg(short, long, required_unless_present = "write_default_adapters")]
    input: Option<PathBuf>,
    /// Output file (defaults to stdout)
    #[arg(short, long, conflicts_with = "barcode_dir")]
    output: Option<PathBuf>,
    /// Output format; auto chooses from output filename or input format
    #[arg(long, default_value = "auto", value_parser = ["auto", "fasta", "fastq", "fasta.gz", "fastq.gz"])]
    format: String,
    /// Level of progress information (0-3)
    #[arg(short, long, default_value_t = 1)]
    verbosity: u8,
    /// Number of threads to use for adapter alignment
    #[arg(short, long, default_value_t = default_threads())]
    threads: usize,
    /// Directory for barcode-specific files
    #[arg(short = 'b', long = "barcode_dir")]
    barcode_dir: Option<PathBuf>,
    /// Minimum barcode percent identity
    #[arg(
        long = "barcode_threshold",
        default_value_t = 75.0,
        allow_negative_numbers = true
    )]
    barcode_threshold: f64,
    /// Minimum difference between best and second-best barcodes
    #[arg(
        long = "barcode_diff",
        default_value_t = 5.0,
        allow_negative_numbers = true
    )]
    barcode_diff: f64,
    /// Require matching barcodes at both ends
    #[arg(long = "require_two_barcodes")]
    require_two_barcodes: bool,
    /// Bin reads without trimming them
    #[arg(long, requires = "barcode_dir")]
    untrimmed: bool,
    /// Discard reads not assigned to a barcode
    #[arg(long = "discard_unassigned")]
    discard_unassigned: bool,
    /// JSON file containing adapter sets (replaces defaults)
    #[arg(long = "adapter_file", conflicts_with = "write_default_adapters")]
    adapter_file: Option<PathBuf>,
    /// Write default adapter sets to a new JSON file and exit
    #[arg(long = "write_default_adapters", conflicts_with = "input")]
    write_default_adapters: Option<PathBuf>,
    /// Minimum identity for an adapter set to be considered present
    #[arg(
        long = "adapter_threshold",
        default_value_t = 90.0,
        allow_negative_numbers = true
    )]
    adapter_threshold: f64,
    /// Number of reads examined for adapter-set selection
    #[arg(long = "check_reads", default_value_t = 10000)]
    check_reads: usize,
    /// Alignment scores: match,mismatch,gap open,gap extend
    #[arg(
        long = "scoring_scheme",
        default_value = "3,-6,-5,-2",
        allow_hyphen_values = true
    )]
    scoring_scheme: String,
    /// Number of bases at each end to search
    #[arg(long = "end_size", default_value_t = 150)]
    end_size: usize,
    /// Minimum end adapter alignment length
    #[arg(long = "min_trim_size", default_value_t = 4)]
    min_trim_size: usize,
    /// Extra bases to trim after an end adapter
    #[arg(long = "extra_end_trim", default_value_t = 2)]
    extra_end_trim: usize,
    /// Minimum partial identity for an end adapter
    #[arg(
        long = "end_threshold",
        default_value_t = 75.0,
        allow_negative_numbers = true
    )]
    end_threshold: f64,
    /// Do not split reads at middle adapters
    #[arg(long = "no_split")]
    no_split: bool,
    /// Discard reads with middle adapters
    #[arg(long = "discard_middle")]
    discard_middle: bool,
    /// Minimum full identity for a middle adapter
    #[arg(
        long = "middle_threshold",
        default_value_t = 90.0,
        allow_negative_numbers = true
    )]
    middle_threshold: f64,
    /// Extra bases to remove on the good side of a middle adapter
    #[arg(long = "extra_middle_trim_good_side", default_value_t = 10)]
    extra_middle_trim_good_side: usize,
    /// Extra bases to remove on the bad side of a middle adapter
    #[arg(long = "extra_middle_trim_bad_side", default_value_t = 100)]
    extra_middle_trim_bad_side: usize,
    /// Minimum output length for split read pieces
    #[arg(long = "min_split_read_size", default_value_t = 1000)]
    min_split_read_size: usize,
}

fn default_threads() -> usize {
    std::thread::available_parallelism().map_or(1, |n| n.get().min(16))
}

fn message(args: &Args, text: &str) {
    if args.verbosity == 0 {
        return;
    }
    if args.output.is_none() && args.barcode_dir.is_none() {
        eprintln!("{text}");
    } else {
        println!("{text}");
    }
}

fn thousands(value: usize) -> String {
    let s = value.to_string();
    let mut result = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            result.push(',');
        }
        result.push(c);
    }
    result
}

fn run(args: Args) -> Result<()> {
    if let Some(path) = &args.write_default_adapters {
        return adapters::write_defaults(path);
    }
    if args.threads == 0 {
        bail!("at least one thread required");
    }
    if args.end_size == 0 {
        bail!("end_size must be greater than zero");
    }
    if [
        args.adapter_threshold,
        args.end_threshold,
        args.middle_threshold,
        args.barcode_threshold,
        args.barcode_diff,
    ]
    .iter()
    .any(|x| !x.is_finite())
    {
        bail!("thresholds must be finite");
    }
    let scores: Vec<i32> = args
        .scoring_scheme
        .split(',')
        .map(str::parse)
        .collect::<std::result::Result<_, _>>()
        .context("incorrectly formatted scoring scheme")?;
    let scores: [i32; 4] = scores
        .try_into()
        .map_err(|_| anyhow::anyhow!("incorrectly formatted scoring scheme"))?;
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;
    let s = Settings {
        verbosity: args.verbosity,
        end_size: args.end_size,
        end_threshold: args.end_threshold,
        min_trim_size: args.min_trim_size,
        extra_end_trim: args.extra_end_trim,
        middle_threshold: args.middle_threshold,
        extra_middle_trim_good_side: args.extra_middle_trim_good_side,
        extra_middle_trim_bad_side: args.extra_middle_trim_bad_side,
        barcode_threshold: args.barcode_threshold,
        barcode_diff: args.barcode_diff,
        require_two_barcodes: args.require_two_barcodes,
        scores,
    };
    let mut config = adapters::load(args.adapter_file.as_deref())?;
    let input = args.input.as_ref().unwrap();
    message(&args, &format!("\nLoading reads\n{}", input.display()));
    let (mut reads, check, input_format) = reads::load(input, args.check_reads)?;
    message(
        &args,
        &format!(
            "{} reads loaded\n\nLooking for known adapter sets",
            reads.len()
        ),
    );
    process::discover(&mut config.adapters, &reads, &check, &s);
    let mut matching: Vec<_> = config
        .adapters
        .into_iter()
        .filter(|a| {
            !a.name.contains("(full sequence)")
                && a.best_start.max(a.best_end) >= args.adapter_threshold
        })
        .collect();
    process::fix_up_1d2(&mut matching);
    let orientation = if args.barcode_dir.is_some() {
        let value = process::choose_orientation(&matching)?;
        message(
            &args,
            &format!("\nBarcodes determined to be in {value} orientation"),
        );
        Some(value)
    } else {
        None
    };
    if args.verbosity > 0 {
        for a in &matching {
            message(
                &args,
                &format!("  {}  {:.1}  {:.1}", a.name, a.best_start, a.best_end),
            );
        }
    }
    adapters::expand(&mut matching, &config.barcode_contexts);
    if !matching.is_empty() {
        message(&args, "\nTrimming adapters from read ends");
        if args.verbosity > 0 {
            for adapter in &matching {
                for seq in [&adapter.start, &adapter.end].into_iter().flatten() {
                    message(&args, &format!("  {}: {}", seq.name, seq.sequence));
                }
            }
        }
        process::trim_ends(&mut reads, &matching, &s, orientation);
        if args.verbosity >= 2 {
            for read in &reads {
                message(
                    &args,
                    &verbose::end_summary(
                        read,
                        args.end_size,
                        args.extra_end_trim,
                        orientation.is_some(),
                        args.verbosity,
                    ),
                );
            }
        }
        let start_count = reads.iter().filter(|r| r.start_trim > 0).count();
        let end_count = reads.iter().filter(|r| r.end_trim > 0).count();
        let start_total: usize = reads.iter().map(|r| r.start_trim).sum();
        let end_total: usize = reads.iter().map(|r| r.end_trim).sum();
        message(
            &args,
            &format!(
                "{start_count} / {} reads had adapters trimmed from their start ({start_total} bp removed)\n{end_count} / {} reads had adapters trimmed from their end ({end_total} bp removed)",
                reads.len(),
                reads.len()
            ),
        );
        if !args.no_split {
            let verb = if args.discard_middle || args.barcode_dir.is_some() {
                "Discarding"
            } else {
                "Splitting"
            };
            message(&args, &format!("\n{verb} reads containing middle adapters"));
            process::trim_middle(&mut reads, &matching, &s)?;
            if args.verbosity >= 2 {
                for read in &reads {
                    if !read.middle_hits.is_empty() {
                        message(&args, &verbose::middle_summary(read));
                    }
                }
            }
            let count = reads.iter().filter(|r| !r.middle.is_empty()).count();
            let result = if args.discard_middle || args.barcode_dir.is_some() {
                "discarded"
            } else {
                "split"
            };
            message(
                &args,
                &format!(
                    "{} / {} reads were {result} based on middle adapters",
                    thousands(count),
                    thousands(reads.len())
                ),
            );
        }
    } else {
        message(
            &args,
            "\nNo adapters found - output reads are unchanged from input reads",
        );
    }
    let name = args
        .output
        .as_ref()
        .map(|p| p.to_string_lossy().to_lowercase());
    let output_format = if args.format != "auto" {
        args.format.clone()
    } else if let Some(name) = name {
        if name.contains(".fasta.gz") {
            "fasta.gz".into()
        } else if name.contains(".fastq.gz") {
            "fastq.gz".into()
        } else if name.contains(".fasta") {
            "fasta".into()
        } else if name.contains(".fastq") {
            "fastq".into()
        } else if input_format == Format::Fasta {
            "fasta".into()
        } else {
            "fastq".into()
        }
    } else {
        let base = if input_format == Format::Fasta {
            "fasta"
        } else {
            "fastq"
        };
        if args.barcode_dir.is_some() && input.to_string_lossy().to_lowercase().ends_with(".gz") {
            format!("{base}.gz")
        } else {
            base.into()
        }
    };
    let format = if output_format == "fasta"
        || (output_format == "fasta.gz" && (args.output.is_some() || args.barcode_dir.is_some()))
    {
        Format::Fasta
    } else {
        Format::Fastq
    };
    let gzip =
        output_format.ends_with(".gz") && (args.output.is_some() || args.barcode_dir.is_some());
    let dest = if args.barcode_dir.is_some() {
        "barcode-specific files"
    } else if args.output.is_some() {
        "file"
    } else {
        "stdout"
    };
    let trimmed = if args.untrimmed {
        "untrimmed"
    } else {
        "trimmed"
    };
    let action = if args.output.is_none() && args.barcode_dir.is_none() {
        "Outputting"
    } else {
        "Saving"
    };
    message(&args, &format!("\n{action} {trimmed} reads to {dest}"));
    let counts = reads::output(
        &reads,
        reads::OutputSettings {
            output: args.output.as_deref(),
            barcode_dir: args.barcode_dir.as_deref(),
            format,
            gzip,
            min_size: args.min_split_read_size,
            discard_middle: args.discard_middle || args.barcode_dir.is_some(),
            untrimmed: args.untrimmed,
            discard_unassigned: args.discard_unassigned,
        },
    )?;
    if let Some(dir) = &args.barcode_dir {
        let name_width = counts.keys().map(String::len).max().unwrap_or(0).max(7);
        let reads_width = counts
            .values()
            .map(|(n, _)| thousands(*n).len())
            .max()
            .unwrap_or(0)
            .max(5);
        let bases_width = counts
            .values()
            .map(|(_, n)| thousands(*n).len())
            .max()
            .unwrap_or(0)
            .max(5);
        message(
            &args,
            &format!(
                "\n  {:<name_width$}  {:>reads_width$}  {:>bases_width$}  File",
                "Barcode", "Reads", "Bases"
            ),
        );
        let ext = if format == Format::Fasta {
            "fasta"
        } else {
            "fastq"
        };
        let suffix = if gzip { ".gz" } else { "" };
        for (name, (num_reads, bases)) in counts {
            let path = dir.join(format!("{name}.{ext}{suffix}"));
            message(
                &args,
                &format!(
                    "  {name:<name_width$}  {:>reads_width$}  {:>bases_width$}  {}",
                    thousands(num_reads),
                    thousands(bases),
                    path.display()
                ),
            );
        }
    }
    message(&args, "Done");
    Ok(())
}

fn main() {
    if let Err(error) = run(Args::parse()) {
        eprintln!("Error: {error:#}");
        std::process::exit(1);
    }
}
