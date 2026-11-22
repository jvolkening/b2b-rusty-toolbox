use std::io::Write;
use std::process::{Command, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};
use viral_unitig::{Assembler, AssemblyStats, Options};

fn reference(len: usize) -> Vec<u8> {
    let mut state = 0x92b7_a1e5_3f6c_8471u64;
    (0..len)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state & 3) as usize]
        })
        .collect()
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => panic!("unexpected base"),
        })
        .collect()
}

fn sequences(data: &[u8]) -> Vec<Vec<u8>> {
    let mut result: Vec<Vec<u8>> = Vec::new();
    for line in data.split(|&b| b == b'\n') {
        if line.starts_with(b">") {
            result.push(Vec::new());
        } else if !line.is_empty() {
            result.last_mut().unwrap().extend_from_slice(line);
        }
    }
    result
}

fn assemble(reads: &[Vec<u8>], options: Options) -> Vec<Vec<u8>> {
    assemble_with_stats(reads, options).0
}

fn assemble_with_stats(reads: &[Vec<u8>], options: Options) -> (Vec<Vec<u8>>, AssemblyStats) {
    let mut assembler = Assembler::new(options).unwrap();
    for read in reads {
        assembler.add_read(read).unwrap();
    }
    let mut output = Vec::new();
    let stats = assembler.write_unitigs(&mut output).unwrap();
    (sequences(&output), stats)
}

fn contains_either_strand(unitigs: &[Vec<u8>], sequence: &[u8]) -> bool {
    unitigs.iter().any(|seq| {
        seq.windows(sequence.len()).any(|window| window == sequence)
            || reverse_complement(seq)
                .windows(sequence.len())
                .any(|window| window == sequence)
    })
}

#[test]
fn assembles_short_overlapping_reads_from_both_strands() {
    let genome = reference(1400);
    let mut reads = Vec::new();
    for start in (0..=1200).step_by(25) {
        for i in 0..4 {
            let read = genome[start..start + 200].to_vec();
            reads.push(if i % 2 == 0 {
                read
            } else {
                reverse_complement(&read)
            });
        }
    }
    let output = assemble(&reads, Options::default());
    assert_eq!(output.len(), 1);
    assert!(output[0] == genome || output[0] == reverse_complement(&genome));
}

#[test]
fn retains_branches_for_two_supported_alleles() {
    let genome = reference(400);
    let mut variant = genome.clone();
    variant[200] = match variant[200] {
        b'A' => b'C',
        _ => b'A',
    };
    let mut reads = vec![genome.clone(); 30];
    reads.extend(vec![variant.clone(); 10]);
    let output = assemble(
        &reads,
        Options {
            clean_bubbles: true,
            clean_errors: false,
            ..Options::default()
        },
    );
    // Branches are emitted separately, not arbitrarily resolved into one strain.
    assert!(output.iter().any(|seq| {
        seq.windows(21).any(|w| w == &genome[190..211])
            || reverse_complement(seq)
                .windows(21)
                .any(|w| w == &genome[190..211])
    }));
    assert!(output.iter().any(|seq| {
        seq.windows(21).any(|w| w == &variant[190..211])
            || reverse_complement(seq)
                .windows(21)
                .any(|w| w == &variant[190..211])
    }));
    assert!(!output.iter().any(|seq| seq.len() == genome.len()));
}

#[test]
fn clean_mode_preserves_distinct_low_coverage_genome() {
    let genome = reference(500);
    let rare = reference(1200)[1000..1120].to_vec();
    let mut reads = vec![genome.clone(); 100];
    reads.extend(vec![rare.clone(); 3]);
    let (output, _) = assemble_with_stats(&reads, Options::default());
    assert!(contains_either_strand(&output, &genome));
    assert!(contains_either_strand(&output, &rare));
}

#[test]
fn clean_mode_keeps_well_supported_alternative_allele() {
    let genome = reference(400);
    let mut alternate = genome.clone();
    alternate[200] = if alternate[200] == b'A' { b'C' } else { b'A' };
    let mut reads = vec![genome.clone(); 100];
    reads.extend(vec![alternate.clone(); 30]);
    let (output, _) = assemble_with_stats(&reads, Options::default());
    assert!(contains_either_strand(&output, &genome[190..211]));
    assert!(contains_either_strand(&output, &alternate[190..211]));
}

#[test]
fn saves_ambiguous_similar_fragment_when_requested() {
    let genome = reference(500);
    let mut short = genome[200..225].to_vec();
    for pos in [8, 14] {
        short[pos] = if short[pos] == b'A' { b'C' } else { b'A' };
    }
    let mut assembler = Assembler::new(Options::default()).unwrap();
    for _ in 0..100 {
        assembler.add_read(&genome).unwrap();
    }
    for _ in 0..4 {
        assembler.add_read(&short).unwrap();
    }
    let mut output = Vec::new();
    let mut uncertain = Vec::new();
    let stats = assembler
        .write_unitigs_with_uncertain(&mut output, Some(&mut uncertain))
        .unwrap();
    assert!(contains_either_strand(&sequences(&output), &genome));
    assert!(contains_either_strand(&sequences(&uncertain), &short));
    assert!(stats.separated_uncertain_unitigs > 0);
    assert_eq!(stats.discarded_uncertain_unitigs, 0);
}

#[test]
fn cleans_weak_complete_bubble_but_preserves_unrelated_rare_sequence() {
    let genome = reference(400);
    let rare = reference(800)[600..720].to_vec();
    let mut error = genome.clone();
    error[200] = if error[200] == b'A' { b'C' } else { b'A' };
    let mut reads = vec![genome.clone(); 100];
    reads.extend(vec![error.clone(); 4]);
    reads.extend(vec![rare.clone(); 3]);

    let (cleaned, stats) = assemble_with_stats(
        &reads,
        Options {
            clean_bubbles: true,
            clean_errors: false,
            ..Options::default()
        },
    );
    assert!(stats.removed_bubble_paths > 0);
    assert!(stats.removed_bubble_kmers > 0);
    assert!(contains_either_strand(&cleaned, &genome));
    assert!(contains_either_strand(&cleaned, &rare));
    assert!(!contains_either_strand(&cleaned, &error[190..211]));

    let (uncleaned, stats) = assemble_with_stats(
        &reads,
        Options {
            clean_errors: false,
            ..Options::default()
        },
    );
    assert_eq!(stats.removed_bubble_paths, 0);
    assert!(contains_either_strand(&uncleaned, &error[190..211]));
    assert!(cleaned.len() < uncleaned.len());

    let (unfiltered, stats) = assemble_with_stats(
        &reads,
        Options {
            clean_bubbles: true,
            clean_errors: false,
            min_fraction: 0.0,
            ..Options::default()
        },
    );
    assert_eq!(stats.removed_bubble_paths, 0);
    assert!(contains_either_strand(&unfiltered, &error[190..211]));
}

#[test]
fn cleans_weak_single_insertion_bubble() {
    let genome = reference(400);
    let mut error = genome.clone();
    error.insert(200, if genome[200] == b'A' { b'C' } else { b'A' });
    let mut reads = vec![genome.clone(); 100];
    reads.extend(vec![error; 4]);
    let (output, stats) = assemble_with_stats(
        &reads,
        Options {
            clean_bubbles: true,
            clean_errors: false,
            ..Options::default()
        },
    );
    assert!(stats.removed_bubble_paths > 0);
    assert!(contains_either_strand(&output, &genome));
}

#[test]
fn does_not_clean_a_dissimilar_weak_bubble() {
    let genome = reference(400);
    let mut alternate = genome.clone();
    for base in &mut alternate[200..203] {
        *base = match *base {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
    }
    let mut reads = vec![genome; 100];
    reads.extend(vec![alternate.clone(); 4]);
    let (output, stats) = assemble_with_stats(
        &reads,
        Options {
            clean_bubbles: true,
            clean_errors: false,
            ..Options::default()
        },
    );
    assert_eq!(stats.removed_bubble_paths, 0);
    assert!(contains_either_strand(&output, &alternate[195..216]));
}

#[test]
fn circular_unitig_is_marked_and_contains_its_overlap() {
    let genome = reference(200);
    let mut wrapped = genome.clone();
    wrapped.extend_from_slice(&genome[..20]);
    let mut assembler = Assembler::new(Options::default()).unwrap();
    for _ in 0..3 {
        assembler.add_read(&wrapped).unwrap();
    }
    let mut output = Vec::new();
    assembler.write_unitigs(&mut output).unwrap();
    assert!(String::from_utf8_lossy(&output).contains("circular=true"));
    let seqs = sequences(&output);
    assert_eq!(seqs.len(), 1);
    assert_eq!(seqs[0].len(), genome.len() + 20);
    assert_eq!(&seqs[0][..20], &seqs[0][genome.len()..]);
    let mut doubled = genome.clone();
    doubled.extend_from_slice(&genome);
    assert!(
        doubled
            .windows(genome.len())
            .any(|w| w == &seqs[0][..genome.len()])
            || reverse_complement(&doubled)
                .windows(genome.len())
                .any(|w| w == &seqs[0][..genome.len()])
    );
}

#[test]
fn filters_substitution_errors_at_high_coverage() {
    let genome = reference(600);
    let mut reads = Vec::new();
    for i in 0..500 {
        let mut read = genome.clone();
        let position = (i * 43) % genome.len();
        read[position] = match read[position] {
            b'A' => b'C',
            _ => b'A',
        };
        reads.push(read);
    }
    let output = assemble(&reads, Options::default());
    assert!(
        output
            .iter()
            .any(|seq| *seq == genome || *seq == reverse_complement(&genome))
    );
}

#[test]
fn assembles_noisy_short_reads_at_high_depth() {
    let genome = reference(1200);
    let mut state = 0xc018_75bd_4aa6_3517u64;
    let mut random = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut reads = Vec::new();
    for _ in 0..6000 {
        let start = (random() as usize) % (genome.len() - 200 + 1);
        let mut read = Vec::new();
        for &original_base in &genome[start..start + 200] {
            let event = random() % 200;
            if event == 0 {
                continue;
            }
            if event == 1 {
                read.push(b"ACGT"[(random() % 4) as usize]);
            }
            let mut base = original_base;
            if random() % 100 < 1 {
                let offset = 1 + (random() % 3) as usize;
                let original = b"ACGT".iter().position(|b| *b == base).unwrap();
                base = b"ACGT"[(original + offset) % 4];
            }
            read.push(base);
        }
        reads.push(read);
    }
    let output = assemble(&reads, Options::default());
    assert_eq!(output.len(), 1, "unitigs: {}", output.len());
    let longest_exact = output
        .iter()
        .filter(|seq| {
            let reverse = reverse_complement(seq);
            genome.windows(seq.len()).any(|w| w == *seq)
                || genome.windows(seq.len()).any(|w| w == reverse)
        })
        .map(Vec::len)
        .max()
        .unwrap_or(0);
    assert!(
        longest_exact >= 1000,
        "longest exact unitig: {longest_exact}"
    );
}

#[test]
fn ambiguous_bases_break_paths_and_short_reads_are_accepted() {
    let mut sequence = reference(120);
    sequence[60] = b'N';
    let mut reads = vec![sequence; 3];
    reads.push(b"ACGT".to_vec());
    let output = assemble(&reads, Options::default());
    assert_eq!(output.len(), 2);
    assert_eq!(output.iter().map(Vec::len).sum::<usize>(), 119);
}

#[test]
fn validates_wrapped_fastq_and_fails_on_truncation() {
    let mut assembler = Assembler::new(Options::default()).unwrap();
    assembler
        .read_fastx(&b"@r\nACGTACGTAC\nGTACGTACGT\n+\nIIIIIIIIII\nIIIIIIIIII\n"[..])
        .unwrap();
    assert_eq!(assembler.reads, 1);
    assert!(assembler.read_fastx(&b"@r\nACGT\n+\nIII\n"[..]).is_err());
    assert!(assembler.read_fastx(&b"@r\nACGT\n+\nIIIII\n"[..]).is_err());
    assert!(assembler.read_fastx(&b">r\nACGT?\n"[..]).is_err());
}

#[test]
fn distinct_kmer_limit_fails_instead_of_dropping_data() {
    let mut assembler = Assembler::new(Options {
        max_kmers: 1,
        ..Options::default()
    })
    .unwrap();
    assert!(
        assembler
            .add_read(&reference(100))
            .unwrap_err()
            .to_string()
            .contains("limit")
    );
}

#[test]
fn cli_reports_bad_input_without_emitting_fasta() {
    let mut child = Command::new(env!("CARGO_BIN_EXE_viral-unitig"))
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(b"@read\nACGT\n+\nII\n")
        .unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("truncated FASTQ quality"));
}

#[test]
fn cli_refuses_to_overwrite_uncertain_output() {
    let existing = concat!(env!("CARGO_MANIFEST_DIR"), "/Cargo.toml");
    let mut child = Command::new(env!("CARGO_BIN_EXE_viral-unitig"))
        .args(["--uncertain", existing])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(b">read\nACGTACGTACGTACGTACGTACGT\n")
        .unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
}

#[test]
fn cli_writes_uncertain_fragments_separately() {
    let genome = reference(500);
    let mut short = genome[200..225].to_vec();
    for pos in [8, 14] {
        short[pos] = if short[pos] == b'A' { b'C' } else { b'A' };
    }
    let unique = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path =
        std::env::temp_dir().join(format!("viral-unitig-{}-{unique}.fa", std::process::id()));
    let mut child = Command::new(env!("CARGO_BIN_EXE_viral-unitig"))
        .arg("--uncertain")
        .arg(&path)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut input = Vec::new();
    for sequence in std::iter::repeat_n(&genome, 100).chain(std::iter::repeat_n(&short, 4)) {
        input.extend_from_slice(b">read\n");
        input.extend_from_slice(sequence);
        input.push(b'\n');
    }
    child.stdin.take().unwrap().write_all(&input).unwrap();
    let result = child.wait_with_output().unwrap();
    let uncertain = std::fs::read(&path).unwrap();
    std::fs::remove_file(&path).unwrap();
    assert!(result.status.success());
    assert!(contains_either_strand(&sequences(&result.stdout), &genome));
    assert!(contains_either_strand(&sequences(&uncertain), &short));
}
