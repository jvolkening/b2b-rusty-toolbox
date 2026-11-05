use std::path::PathBuf;
use std::process::Command;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(name)
}

fn run(extra: &[&str]) -> String {
    let out = Command::new(env!("CARGO_BIN_EXE_mega-blastn"))
        .arg("-db")
        .arg(data("tiny"))
        .arg("-query")
        .arg(data("q.fa"))
        .args(extra)
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout).unwrap()
}

#[test]
fn tabular_matches_blastn() {
    let expected = std::fs::read_to_string(data("expected.6")).unwrap();
    assert_eq!(run(&[]), expected);
    assert_eq!(run(&["-outfmt", "6", "-num_threads", "1"]), expected);
}

#[test]
fn commented_tabular_matches_blastn() {
    let expected = std::fs::read_to_string(data("expected.7")).unwrap();
    let got = run(&["-outfmt", "7"]);
    // the database line echoes the -db argument, which differs by location
    let norm = |s: &str| {
        s.lines()
            .filter(|l| !l.starts_with("# Database:"))
            .collect::<Vec<_>>()
            .join("\n")
    };
    assert_eq!(norm(&got), norm(&expected));
}

#[test]
fn custom_fields_and_strand() {
    let out = run(&["-outfmt", "6 qseqid sseqid sstrand qcovs"]);
    let lines: Vec<&str> = out.lines().collect();
    assert_eq!(lines.len(), 3);
    assert_eq!(lines[1], "q1\tref|NC_000001.1|\tminus\t100");
    let plus = run(&["-strand", "plus"]);
    assert!(!plus.contains("q1\t"));
}

#[test]
fn rejects_bad_options() {
    let out = Command::new(env!("CARGO_BIN_EXE_mega-blastn"))
        .args(["-db", "nonexistent", "-outfmt", "5"])
        .output()
        .unwrap();
    assert!(!out.status.success());
}
