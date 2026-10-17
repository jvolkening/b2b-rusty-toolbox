use std::fs;
use std::io::Write;
use std::process::{Command, Stdio};

fn run_bin(input: &[u8], k: &str) -> std::process::Output {
     let mut child = Command::new(env!("CARGO_BIN_EXE_n-every"))
         .arg(k)
         .stdin(Stdio::piped())
         .stdout(Stdio::piped())
         .stderr(Stdio::piped())
         .spawn()
         .unwrap();
     // The process may exit without reading stdin (e.g. invalid K), which
     // makes this write fail with BrokenPipe. That's expected, so ignore it;
     // tests assert on the exit status and output instead.
     let _ = child.stdin.take().unwrap().write_all(input);
     child.wait_with_output().unwrap()
 }

#[test]
fn fasta_fixtures_match_expected() {
    // Convention: tests/data/<name>.fa + <name>.k<K>.expected.fa
    let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data");
    let mut checked = 0;
    for entry in fs::read_dir(dir).unwrap() {
        let path = entry.unwrap().path();
        let fname = path.file_name().unwrap().to_string_lossy().to_string();
        if let Some((stem, rest)) = fname.split_once(".k") {
            if let Some(k) = rest.strip_suffix(".expected.fa") {
                let input = fs::read(format!("{dir}/{stem}.fa")).unwrap();
                let out = run_bin(&input, k);
                assert!(out.status.success(), "{fname}");
                assert_eq!(out.stdout, fs::read(&path).unwrap(), "{fname}");
                checked += 1;
            }
        }
    }
    assert!(checked > 0, "no fixture files found in tests/data");
}

#[test]
fn help_prints_usage() {
    let out = Command::new(env!("CARGO_BIN_EXE_n-every")).arg("--help").output().unwrap();
    assert!(out.status.success());
    assert!(String::from_utf8_lossy(&out.stdout).contains("USAGE"));
}

#[test]
fn rejects_zero_and_non_numeric_k() {
    for bad in ["0", "abc", "-1"] {
        assert!(!run_bin(b">a\nACGT\n", bad).status.success(), "{bad}");
    }
}
