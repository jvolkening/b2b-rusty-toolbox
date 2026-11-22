use flate2::{Compression, write::GzEncoder};
use std::io::Write;
use std::process::{Command, Output, Stdio};

fn dna(len: usize, mut state: u64) -> Vec<u8> {
    (0..len)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state & 3) as usize]
        })
        .collect()
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => panic!("unexpected base"),
        })
        .collect()
}

fn run(input: &[u8], args: &[&str]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_ir-chimera"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(input).unwrap();
    child.wait_with_output().unwrap()
}

fn hairpin(len: usize) -> Vec<u8> {
    let mut seq = dna(len, 0x1234567);
    seq.extend(revcomp(&seq));
    seq
}

fn temporary_path(label: &str) -> std::path::PathBuf {
    std::env::temp_dir().join(format!(
        "ir-chimera-{}-{}-{label}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ))
}

fn assert_log(path: &std::path::Path, total: u64, kept: u64, fraction: f64) {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    assert_eq!(lines.next(), Some("n_total\tn_kept\tfraction"));
    let fields: Vec<_> = lines.next().unwrap().split('\t').collect();
    assert_eq!(fields.len(), 3);
    assert_eq!(fields[0].parse::<u64>().unwrap(), total);
    assert_eq!(fields[1].parse::<u64>().unwrap(), kept);
    assert_eq!(fields[2].parse::<f64>().unwrap(), fraction);
    assert!(lines.next().is_none());
    std::fs::remove_file(path).unwrap();
}

#[test]
fn log_counts_all_batches_and_inputs_in_both_modes() {
    let fasta = temporary_path("reads.fa");
    let normal = b">normal\nACGT\n".repeat(12);
    std::fs::write(&fasta, &normal).unwrap();
    let chimera = format!(">chimera\n{}\n", String::from_utf8_lossy(&hairpin(700)));
    let log = temporary_path("summary.tsv");
    let args = [
        "--threads",
        "1",
        "--log",
        log.to_str().unwrap(),
        fasta.to_str().unwrap(),
        "-",
    ];
    let output = run(chimera.as_bytes(), &args);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, normal);
    assert_log(&log, 13, 12, 12.0 / 13.0);

    let args = [
        "--threads",
        "2",
        "--chimeric",
        "--log",
        log.to_str().unwrap(),
        fasta.to_str().unwrap(),
        "-",
    ];
    let output = run(chimera.as_bytes(), &args);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, chimera.as_bytes());
    assert_log(&log, 13, 1, 1.0 / 13.0);
    std::fs::remove_file(fasta).unwrap();
}

#[test]
fn empty_log_and_log_errors() {
    let log = temporary_path("empty.tsv");
    let output = run(b"", &["--log", log.to_str().unwrap()]);
    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert_log(&log, 0, 0, 0.0);

    std::fs::write(&log, b"existing").unwrap();
    let output = run(b">a\nACGT\n", &["--log", log.to_str().unwrap()]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("cannot create log"));
    assert_eq!(std::fs::read(&log).unwrap(), b"existing");
    std::fs::remove_file(&log).unwrap();

    let output = run(b">a\nACGT!\n", &["--log", log.to_str().unwrap()]);
    assert!(!output.status.success());
    assert!(!log.exists());
}

#[test]
fn selects_reads_and_preserves_fastq_records() {
    let normal = dna(1400, 0x87654321);
    let chimera = hairpin(700);
    let normal_record = format!(
        "@normal with spaces\n{}\n+normal\n{}\n",
        String::from_utf8_lossy(&normal),
        "I".repeat(normal.len())
    );
    let chimera_record = format!(
        "@chimera description\n{}\n+chimera\n{}\n",
        String::from_utf8_lossy(&chimera),
        "J".repeat(chimera.len())
    );
    let input = format!("{normal_record}{chimera_record}");
    let output = run(input.as_bytes(), &[]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, normal_record.as_bytes());
    let output = run(input.as_bytes(), &["--chimeric"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, chimera_record.as_bytes());
}

#[test]
fn alignment_report_is_independent_of_output_mode_and_keeps_stdout_fastq() {
    let report_path = temporary_path("alignments.txt");
    let seq = hairpin(700);
    let chimera = format!(
        "@chimera description\n{}\n+chimera\n{}\n",
        String::from_utf8_lossy(&seq),
        "I".repeat(seq.len())
    );
    let normal = "@normal\nACGT\n+\nIIII\n";
    let input = format!("{normal}{chimera}");
    let mut previous_report = None;
    for (mode, threads) in [(false, "1"), (true, "4")] {
        let mut args = vec![
            "--alignments",
            report_path.to_str().unwrap(),
            "--threads",
            threads,
        ];
        if mode {
            args.push("--chimeric");
        }
        let output = run(input.as_bytes(), &args);
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(
            output.stdout,
            if mode {
                chimera.as_bytes()
            } else {
                normal.as_bytes()
            }
        );
        let report = std::fs::read_to_string(&report_path).unwrap();
        assert!(report.contains("Read: chimera description\nLength: 1400\n"));
        assert!(!report.contains("Read: normal"));
        assert!(
            report.contains("Target (+): 1-700\nQuery (-): 1400-701\nIdentity: 700/700 (100.00%)")
        );
        assert!(report.contains(&format!(
            "Target    1  {}  60",
            String::from_utf8_lossy(&seq[..60])
        )));
        assert!(report.contains(&format!(
            "Query  1400  {}  1341",
            String::from_utf8_lossy(&seq[..60])
        )));
        if let Some(previous) = &previous_report {
            assert_eq!(previous, &report);
        }
        previous_report = Some(report);
        std::fs::remove_file(&report_path).unwrap();
    }
}

#[test]
fn alignment_report_handles_empty_inputs_and_path_errors() {
    let path = temporary_path("alignments-empty.txt");
    let output = run(
        b">normal\nACGT\n",
        &["--alignments", path.to_str().unwrap()],
    );
    assert!(output.status.success());
    assert_eq!(std::fs::read(&path).unwrap(), b"");

    let output = run(
        b">normal\nACGT\n",
        &["--alignments", path.to_str().unwrap()],
    );
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("cannot create alignment report"));
    assert_eq!(std::fs::read(&path).unwrap(), b"");
    std::fs::remove_file(&path).unwrap();
    let output = run(
        b"",
        &[
            "--log",
            path.to_str().unwrap(),
            "--alignments",
            path.to_str().unwrap(),
        ],
    );
    assert!(!output.status.success());
    assert!(!path.exists());
}

#[test]
fn alignment_report_respects_fractional_threshold_and_logs() {
    let path = temporary_path("fraction-alignments.txt");
    let log = temporary_path("fraction-log.tsv");
    let seq = hairpin(350);
    let input = format!(">small\n{}\n>normal\nACGT\n", String::from_utf8_lossy(&seq));
    let output = run(
        input.as_bytes(),
        &[
            "--alignments",
            path.to_str().unwrap(),
            "--log",
            log.to_str().unwrap(),
            "--min-fraction",
            "1",
        ],
    );
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, b">normal\nACGT\n");
    let report = std::fs::read_to_string(&path).unwrap();
    assert!(report.contains("Read: small\nLength: 700\nTarget (+): 1-350\nQuery (-): 700-351\nIdentity: 350/350 (100.00%)"));
    assert_log(&log, 2, 1, 0.5);
    std::fs::remove_file(path).unwrap();
}

#[test]
fn fasta_length_and_identity_thresholds() {
    let seq = hairpin(700);
    let mut imperfect = seq.clone();
    // Introduce substitutions and a deletion in one arm.
    let mask = dna(690, 0x467819);
    for (pos, window) in mask.windows(3).enumerate() {
        if window == b"AAA" {
            let base = &mut imperfect[710 + pos];
            *base = if *base == b'A' { b'C' } else { b'A' };
        }
    }
    imperfect.remove(1000);
    let input = format!(">imperfect\n{}\n", String::from_utf8_lossy(&imperfect));
    let output = run(input.as_bytes(), &["-c", "-i", "90"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, input.as_bytes());
    let output = run(input.as_bytes(), &["-c", "-i", "99"]);
    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    let output = run(input.as_bytes(), &["-c", "-l", "800"]);
    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    let output = run(input.as_bytes(), &["-c", "-l", "550", "-i", "90"]);
    assert!(output.status.success());
    assert_eq!(output.stdout, input.as_bytes());
}

#[test]
fn fractional_threshold_is_optional_and_either_length_rule_suffices() {
    let short = format!(">short\n{}\n", String::from_utf8_lossy(&hairpin(350)));
    let input = short.as_bytes();
    let output = run(input, &["--chimeric"]);
    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    let output = run(input, &["--chimeric", "--min-fraction", "1"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, input);
    let output = run(input, &["--min-fraction", "1"]);
    assert!(output.status.success());
    assert!(output.stdout.is_empty());

    let arm = dna(700, 0x4321);
    let seq = [
        dna(500, 0x9753),
        arm.clone(),
        revcomp(&arm),
        dna(500, 0x2468),
    ]
    .concat();
    let long = format!(">long\n{}\n", String::from_utf8_lossy(&seq));
    let output = run(long.as_bytes(), &["--chimeric", "--min-fraction", "1"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, long.as_bytes());
    let output = run(
        long.as_bytes(),
        &[
            "--chimeric",
            "--min-length",
            "1000",
            "--min-fraction",
            "0.5",
        ],
    );
    assert!(output.status.success());
    assert_eq!(output.stdout, long.as_bytes());
    let output = run(
        long.as_bytes(),
        &[
            "--chimeric",
            "--min-length",
            "1000",
            "--min-fraction",
            "0.8",
        ],
    );
    assert!(output.status.success());
    assert!(output.stdout.is_empty());
}

#[test]
fn reads_wrapped_crlf_fasta_and_gzip() {
    let seq = hairpin(600);
    let input = format!(
        ">id description\r\n{}\r\n{}\r\n",
        String::from_utf8_lossy(&seq[..600]),
        String::from_utf8_lossy(&seq[600..])
    );
    let mut gzip = GzEncoder::new(Vec::new(), Compression::default());
    gzip.write_all(input.as_bytes()).unwrap();
    let output = run(&gzip.finish().unwrap(), &["-c"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, input.as_bytes());
}

#[test]
fn empty_and_concatenated_gzip_are_read_and_truncation_is_rejected() {
    let gzip = |input: &[u8]| {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(input).unwrap();
        encoder.finish().unwrap()
    };
    let log = temporary_path("empty-gzip.tsv");
    let output = run(&gzip(b""), &["--log", log.to_str().unwrap()]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stdout.is_empty());
    assert_log(&log, 0, 0, 0.0);

    let mut concatenated = gzip(b">first\nACGT\n");
    concatenated.extend(gzip(b">second\nTGCA\n"));
    let mut with_empty_member = gzip(b"");
    with_empty_member.extend(&concatenated);
    for data in [&concatenated, &with_empty_member] {
        let output = run(data, &[]);
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(output.stdout, b">first\nACGT\n>second\nTGCA\n");
    }

    for mut corrupt in [gzip(b""), gzip(b">read\nACGT\n")] {
        corrupt.truncate(corrupt.len() - 3);
        let output = run(&corrupt, &[]);
        assert!(!output.status.success(), "truncated gzip was accepted");
        assert!(!output.stderr.is_empty());
    }
}

#[test]
fn synthetic_reads_have_expected_calls_across_modes_and_threads() {
    let path = temporary_path("synthetic.fa");
    let mut input = Vec::new();
    let mut expected = [Vec::new(), Vec::new()];
    for i in 1..=16 {
        let arm = dna(800, 0x5246 + i);
        let mut mutated = revcomp(&arm);
        for base in mutated.iter_mut().step_by(40) {
            *base = if *base == b'A' { b'C' } else { b'A' };
        }
        mutated.remove(420);
        let reads = [
            ("random", dna(2000, 0xa420 + i), false),
            (
                "direct",
                [arm.clone(), dna(200, 0x4490 + i), arm.clone()].concat(),
                false,
            ),
            ("clean", [arm.clone(), revcomp(&arm)].concat(), true),
            (
                "noisy",
                [
                    dna(200, 0x1234 + i),
                    arm.clone(),
                    dna(200, 0x9876 + i),
                    mutated,
                ]
                .concat(),
                true,
            ),
        ];
        for (kind, seq, chimeric) in reads {
            let record = format!(">{kind}-{i}\n{}\n", String::from_utf8_lossy(&seq));
            input.extend_from_slice(record.as_bytes());
            expected[usize::from(chimeric)].extend_from_slice(record.as_bytes());
        }
    }
    std::fs::write(&path, input).unwrap();
    for threads in ["1", "4"] {
        for chimeric in [false, true] {
            let mut command = Command::new(env!("CARGO_BIN_EXE_ir-chimera"));
            command.args(["--threads", threads]);
            if chimeric {
                command.arg("--chimeric");
            }
            let output = command.arg(&path).output().unwrap();
            assert!(
                output.status.success(),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
            assert_eq!(output.stdout, expected[usize::from(chimeric)]);
        }
    }
    std::fs::remove_file(path).unwrap();
}

#[test]
fn rejects_malformed_input_and_invalid_options() {
    for input in [
        b"@id\nACGT\n+\n!!!\n".as_slice(),
        b">id\nACGT!\n".as_slice(),
    ] {
        let output = run(input, &[]);
        assert!(!output.status.success());
        assert!(output.stdout.is_empty());
        assert!(!output.stderr.is_empty());
    }
    for args in [
        &["--min-length", "0"][..],
        &["--min-identity", "NaN"][..],
        &["--min-identity", "101"][..],
        &["--min-fraction", "0"][..],
        &["--min-fraction", "-0.1"][..],
        &["--min-fraction", "1.01"][..],
        &["--min-fraction", "NaN"][..],
        &["--min-fraction", "inf"][..],
        &["--min-fraction"][..],
        &["--log", "--chimeric"][..],
        &["--alignments"][..],
        &["--alignments", "--chimeric"][..],
        &["--invalid"][..],
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_ir-chimera"))
            .args(args)
            .output()
            .unwrap();
        assert!(!output.status.success());
    }
}

#[test]
fn short_and_empty_reads_are_nonchimeric() {
    let input = b">empty\n>short\nACGT\n";
    let output = run(input, &[]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, input);
}

#[test]
fn finds_separated_inverted_arms_with_multiple_workers() {
    let arm = dna(800, 0x1341);
    let mut chimera = dna(350, 0x113);
    chimera.extend(&arm);
    chimera.extend(dna(200, 0x87));
    chimera.extend(revcomp(&arm));
    chimera.extend(dna(350, 0x65));
    let normal = dna(chimera.len(), 0x567);
    let input = format!(
        ">normal\n{}\n>chimera\n{}\n",
        String::from_utf8_lossy(&normal),
        String::from_utf8_lossy(&chimera),
    );
    let output = run(input.as_bytes(), &["--threads", "4", "--chimeric"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        format!(">chimera\n{}\n", String::from_utf8_lossy(&chimera)).as_bytes()
    );
}

#[test]
fn supports_multiple_files_and_rejects_mixed_formats() {
    let path = std::env::temp_dir().join(format!(
        "ir-chimera-{}-{}.fa",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&path)
        .unwrap();
    file.write_all(b">first\nACGT\n").unwrap();
    drop(file);
    let path_str = path.to_str().unwrap();
    let output = run(b">second\nTGCA\n", &[path_str, "-"]);
    assert!(output.status.success());
    assert_eq!(output.stdout, b">first\nACGT\n>second\nTGCA\n");
    let output = run(b"@second\nTGCA\n+\nIIII\n", &[path_str, "-"]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("cannot mix"));
    std::fs::remove_file(path).unwrap();
}
