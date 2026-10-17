use std::io::{self, BufRead, Write};

/// Insert an `N` after every `k` sequence characters in each FASTA record.
pub fn insert_n<R: BufRead, W: Write>(mut reader: R, mut writer: W, k: usize) -> io::Result<()> {
    assert!(k > 0, "k must be > 0");
    let mut line = Vec::new();
    let mut current_pos: usize = 0; 
    let mut inside_sequence = false;

    while reader.read_until(b'\n', &mut line)? > 0 {
        if line.starts_with(&[b'>']) {
            if inside_sequence {
                writer.write_all(b"\n")?;
                current_pos = 0;
            }
            writer.write_all(&line)?;
            inside_sequence = false;
        } else {
            inside_sequence = true;
            // Handle cross-platform line endings
            let mut len = line.len();
            while len > 0 && (line[len - 1] == b'\n' || line[len - 1] == b'\r') {
                len -= 1;
            }
            
            let seq_data = &line[..len];
            let mut idx = 0;
            let data_len = seq_data.len();

            while idx < data_len {
                if current_pos == k {
                    writer.write_all(b"N")?;
                    current_pos = 0;
                }
                let chunk_size = std::cmp::min(k - current_pos, data_len - idx);
                writer.write_all(&seq_data[idx..idx + chunk_size])?;
                idx += chunk_size;
                current_pos += chunk_size;
            }
        }
        line.clear();
    }

    if inside_sequence {
        writer.write_all(b"\n")?;
    }

    writer.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(input: &str, k: usize) -> String {
        let mut out = Vec::new();
        insert_n(input.as_bytes(), &mut out, k).unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn wraps_sequence_across_lines() {
        assert_eq!(run(">a\nACGT\nACGT\n", 3), ">a\nACGNTACNGT\n");
    }

    #[test]
    fn counter_resets_per_record() {
        assert_eq!(run(">a\nAC\n>b\nGT\n", 3), ">a\nAC\n>b\nGT\n");
    }
    #[test]
    fn uneven_line_lengths_count_across_lines() {
        // ACGTACG split as AC / GTAC / G; k=3 -> ACG N TAC N G
        assert_eq!(run(">a\nAC\nGTAC\nG\n", 3), ">a\nACGNTACNG\n");
    }
}
