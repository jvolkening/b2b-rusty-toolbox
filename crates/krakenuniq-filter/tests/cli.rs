use std::fs;
use std::process::{Command, Stdio};

#[test]
fn filters_stdin_with_database_fixture() {
    let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/tmp-db");
    fs::create_dir_all(dir).unwrap();
    let taxdb = format!("{dir}/taxDB");
    fs::write(&taxdb, "1\t1\n2\t1\n3\t1\n4\t2\n").unwrap();

    let mut child = Command::new(env!("CARGO_BIN_EXE_krakenuniq-filter"))
        .args(["--db", &taxdb, "--threshold", "0.5"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    use std::io::Write;
    child
        .stdin
        .take()
        .unwrap()
        .write_all(b"C\tread1\t4\t100\t4:5 A:5\n")
        .unwrap();
    let out = child.wait_with_output().unwrap();
    assert!(out.status.success(), "{:?}", out.stderr);
    assert_eq!(out.stdout, b"C\tread1\t4\t100\tP=1.000\t4:5 A:5\n");

    let _ = fs::remove_dir_all(dir);
}

#[test]
fn help_and_argument_errors() {
    let help = Command::new(env!("CARGO_BIN_EXE_krakenuniq-filter"))
        .arg("--help")
        .output()
        .unwrap();
    assert!(help.status.success());
    assert!(String::from_utf8_lossy(&help.stdout).contains("USAGE"));

    let bad = Command::new(env!("CARGO_BIN_EXE_krakenuniq-filter"))
        .args(["--threshold", "1.1"])
        .output()
        .unwrap();
    assert_eq!(bad.status.code(), Some(64));
}
