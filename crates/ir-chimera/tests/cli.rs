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
