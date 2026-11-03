use bg_window::*;
use std::io::Cursor;

fn go(fa: &str, bg: &str, min: usize, max: usize, stat: Stat) -> String {
    let f = read_fasta_lengths(Cursor::new(fa)).unwrap();
    let b = read_bedgraph(Cursor::new(bg)).unwrap();
    let mut out = Vec::new();
    run(&f, &b, min, max, stat, &mut out).unwrap();
    String::from_utf8(out).unwrap()
}

#[test]
fn entries_within_range_per_sequence() {
    let fa = ">a\nACGTACGTAC\nGTACGTACGT\n>b\nACGTACGTACGTACGTACGTACGTACGTAC\n";
    let mut bg = String::new();
    for i in 0..20 { bg += &format!("a\t{i}\t{}\t{}\n", i + 1, i); }
    for i in 0..30 { bg += &format!("b\t{i}\t{}\t{}\n", i + 1, i * 2); }
    let out = go(&fa, &bg, 4, 6, Stat::Max);
    for s in ["a", "b"] {
        let n = out.lines().filter(|l| l.starts_with(&format!("{s}\t"))).count();
        assert!((4..=6).contains(&n), "{s}: {n}");
    }
}

#[test]
fn median_and_gaps() {
    let out = go(">a\nAAAAAAAA\n", "a\t0\t2\t10\na\t4\t6\t10\n", 1, 1, Stat::Median);
    // identical adjacent windows are merged
    assert_eq!(out, "a\t0\t8\t5\n");
}

#[test]
fn adaptive_sizing_escapes_width_cycles() {
    let out = go(
        ">a\nAAAAA\n",
        "a\t0\t1\t0\na\t1\t2\t1\na\t2\t3\t2\na\t3\t4\t0\na\t4\t5\t2\n",
        2,
        2,
        Stat::Mean,
    );
    assert_eq!(out.lines().count(), 2);
}
