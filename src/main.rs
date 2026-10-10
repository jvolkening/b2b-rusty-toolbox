use std::env;
use std::io::{self, BufRead, BufWriter, Write};
use std::process;

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();

    let k: usize = if args.len() < 2 {
        3
    } else {
        match args[1].parse::<usize>() {
            Ok(val) if val > 0 => val,
            _ => {
                eprintln!("Error: 'k' must be a positive integer greater than 0.");
                process::exit(1);
            }
        }
    };

    let stdin = io::stdin();
    let stdout = io::stdout();

    let mut reader = stdin.lock();
    // Large 256 KB internal I/O cache to maximize pipe throughput
    let mut writer = BufWriter::with_capacity(256 * 1024, stdout.lock());

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
            
            // Fast safe line-ending trimming 
            let mut len = line.len();
            while len > 0 && (line[len - 1] == b'\n' || line[len - 1] == b'\r') {
                len -= 1;
            }
            
            let seq_data = &line[..len];
            let mut idx = 0;
            let data_len = seq_data.len();

            while idx < data_len {
                // FIX: Check if we are sitting exactly on a k-boundary interval 
                // from a previous chunk iteration or a previous line block
                if current_pos == k {
                    writer.write_all(b"N")?;
                    current_pos = 0;
                }

                let space_left = k - current_pos;
                let chunk_size = std::cmp::min(space_left, data_len - idx);

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

    writer.flush()?;
    Ok(())
}
