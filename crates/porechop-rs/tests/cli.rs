use flate2::read::MultiGzDecoder;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

struct TestDir(PathBuf);

impl TestDir {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "porechop-rs-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
}

impl Drop for TestDir {
    fn drop(&mut self) {
        for name in [
            "result.fastq.gz",
            "adapters.json",
            "reads.fasta",
            "bad.fastq",
            "BC01.fastq.gz",
            "BC02.fastq.gz",
            "BC03.fastq.gz",
            "none.fastq.gz",
        ] {
            let _ = fs::remove_file(self.0.join(name));
        }
        let _ = fs::remove_dir(&self.0);
    }
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn ungzip(path: &Path) -> Vec<u8> {
    let mut result = Vec::new();
    MultiGzDecoder::new(fs::File::open(path).unwrap())
        .read_to_end(&mut result)
        .unwrap();
    result
}

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_porechop-rs"))
        .args(args)
        .output()
        .unwrap()
}

fn assert_success(result: &Output) {
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
}

#[test]
fn trimming_matches_original_with_threads_and_gzip() {
    let input = fixture("trimming-input.fastq.gz");
    let expected = ungzip(&fixture("trimming-expected.fastq.gz"));
    let stdout = run(&["-i", input.to_str().unwrap(), "-t", "1", "-v", "0"]);
    assert_success(&stdout);
    assert_eq!(stdout.stdout, expected);

    let dir = TestDir::new();
    let output = dir.0.join("result.fastq.gz");
    let file = run(&[
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-t",
        "2",
        "-v",
        "0",
    ]);
    assert_success(&file);
    assert!(file.stdout.is_empty());
    assert_eq!(ungzip(&output), expected);
}

#[test]
fn barcode_bins_match_original() {
    let input = fixture("barcodes-input.fastq.gz");
    let dir = TestDir::new();
    let result = run(&[
        "-i",
        input.to_str().unwrap(),
        "-b",
        dir.0.to_str().unwrap(),
        "-t",
        "2",
        "-v",
        "0",
    ]);
    assert_success(&result);
    let mut actual = fs::read_dir(&dir.0)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect::<Vec<_>>();
    actual.sort();
    let expected = [
        "BC01.fastq.gz",
        "BC02.fastq.gz",
        "BC03.fastq.gz",
        "none.fastq.gz",
    ];
    assert_eq!(actual, expected);
    for name in expected {
        assert_eq!(
            ungzip(&dir.0.join(name)),
            ungzip(&fixture(&format!("barcodes/{name}")))
        );
    }
}

#[test]
fn exports_defaults_and_uses_a_custom_adapter() {
    let dir = TestDir::new();
    let config = dir.0.join("adapters.json");
    let exported = run(&["--write_default_adapters", config.to_str().unwrap()]);
    assert_success(&exported);
    assert!(
        !run(&["--write_default_adapters", config.to_str().unwrap()])
            .status
            .success()
    );

    let mut settings: serde_json::Value =
        serde_json::from_slice(&fs::read(&config).unwrap()).unwrap();
    let adapter = "AATGTACTTCGTTCAGTTACGTATTGCT";
    settings["adapters"] = serde_json::json!([{
        "name": "custom",
        "start": {"name": "custom_start", "sequence": adapter}
    }]);
    fs::write(&config, serde_json::to_vec(&settings).unwrap()).unwrap();
    let input = dir.0.join("reads.fasta");
    let payload = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";
    fs::write(&input, format!(">read\n{adapter}{payload}\n")).unwrap();
    let result = run(&[
        "-i",
        input.to_str().unwrap(),
        "--adapter_file",
        config.to_str().unwrap(),
        "-t",
        "1",
        "-v",
        "0",
        "--no_split",
        "--extra_end_trim",
        "0",
    ]);
    assert_success(&result);
    assert_eq!(result.stdout, format!(">read\n{payload}\n").as_bytes());
}

#[test]
fn rejects_truncated_fastq() {
    let dir = TestDir::new();
    let input = dir.0.join("bad.fastq");
    fs::write(&input, b"@read\nACGT\n+\n").unwrap();
    let result = run(&["-i", input.to_str().unwrap(), "-v", "0"]);
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("truncated FASTQ"));
}
