use crate::reads::ReadRecord;
use std::fmt::Write;

const RED: &str = "\x1b[31m";
const YELLOW: &str = "\x1b[33m";
const RESET: &str = "\x1b[0m";

fn colored(seq: &[u8], start: usize, end: usize, style: impl Fn(usize) -> u8) -> String {
    let mut result = String::new();
    let mut previous = 0;
    for (i, &base) in seq.iter().enumerate().take(end).skip(start) {
        let current = style(i);
        if current != previous {
            if previous != 0 {
                result.push_str(RESET);
            }
            if current == 1 {
                result.push_str(RED);
            }
            if current == 2 {
                result.push_str(YELLOW);
            }
            previous = current;
        }
        result.push(base as char);
    }
    if previous != 0 {
        result.push_str(RESET);
    }
    result
}

fn ends(read: &ReadRecord, start: usize, end: usize, extra: usize) -> String {
    let len = read.seq.len();
    let left_core = read.start_trim.saturating_sub(extra);
    let right_core = read.end_trim.saturating_sub(extra);
    colored(&read.seq, start, end, |pos| {
        if (read.start_trim > 0 && pos < left_core)
            || (read.end_trim > 0 && pos >= len.saturating_sub(right_core))
        {
            1
        } else if (read.start_trim > 0 && pos < read.start_trim)
            || (read.end_trim > 0 && pos >= len.saturating_sub(read.end_trim))
        {
            2
        } else {
            0
        }
    })
}

fn best(scores: &[(String, f64)]) -> (&str, f64) {
    scores.iter().fold(("none", 0.0), |best, (name, score)| {
        if *score > best.1 {
            (name.as_str(), *score)
        } else {
            best
        }
    })
}

pub fn end_summary(
    read: &ReadRecord,
    end_size: usize,
    extra: usize,
    barcoded: bool,
    verbosity: u8,
) -> String {
    let len = read.seq.len();
    let start_end = len.min(end_size);
    let end_start = len.saturating_sub(end_size);
    let mut text = String::new();
    if verbosity == 2 {
        if barcoded {
            let (start_name, start_score) = best(&read.start_scores);
            let (end_name, end_score) = best(&read.end_scores);
            let _ = write!(
                text,
                "start: {start_name} ({start_score:.1}%), end: {end_name} ({end_score:.1}%), barcode call: {}   ",
                read.barcode
            );
        }
        if len <= end_size.saturating_mul(2) {
            text.push_str(&ends(read, 0, len, extra));
        } else {
            text.push_str(&ends(read, 0, start_end, extra));
            text.push_str("...");
            text.push_str(&ends(read, end_start, len, extra));
        }
    } else {
        let _ = writeln!(
            text,
            "{}\n  start: {}...",
            read.name,
            ends(read, 0, start_end, extra)
        );
        if !read.start_hits.is_empty() {
            text.push_str("    start alignments:\n");
            for (name, hit) in &read.start_hits {
                let _ = writeln!(
                    text,
                    "      {name}, full score={}, partial score={}, read position: {}-{}",
                    hit.full, hit.partial, hit.start, hit.end
                );
            }
        }
        let _ = writeln!(text, "  end:   ...{}", ends(read, end_start, len, extra));
        if !read.end_hits.is_empty() {
            text.push_str("    end alignments:\n");
            for (name, hit) in &read.end_hits {
                let _ = writeln!(
                    text,
                    "      {name}, full score={}, partial score={}, read position: {}-{}",
                    hit.full, hit.partial, hit.start, hit.end
                );
            }
        }
        if barcoded {
            let (start_name, start_score) = best(&read.start_scores);
            let (end_name, end_score) = best(&read.end_scores);
            let format_scores = |scores: &[(String, f64)]| {
                scores
                    .iter()
                    .map(|(name, score)| format!("{name} ({score:.1}%)"))
                    .collect::<Vec<_>>()
                    .join(", ")
            };
            let _ = writeln!(
                text,
                "  Barcodes:\n    start barcodes:        {}\n    end barcodes:          {}\n    best start barcode:    {start_name} ({start_score:.1}%)\n    best end barcode:      {end_name} ({end_score:.1}%)",
                format_scores(&read.start_scores),
                format_scores(&read.end_scores)
            );
            if let Some(call) = &read.albacore {
                let _ = writeln!(text, "    albacore barcode call: {call}");
            }
            let _ = writeln!(text, "    final barcode call:    {}", read.barcode);
        }
    }
    text
}

pub fn middle_summary(read: &ReadRecord) -> String {
    let mut text = format!("{}\n", read.name);
    for (name, hit) in &read.middle_hits {
        let _ = writeln!(
            text,
            "  {name} (read coords: {}-{}, identity: {:.1}%)",
            hit.start, hit.end, hit.full
        );
    }
    let seq = read.trimmed();
    if let (Some(first), Some(last)) = (
        read.middle.iter().map(|x| x.0).min(),
        read.middle.iter().map(|x| x.1).max(),
    ) {
        let start = first.saturating_sub(100).min(seq.len());
        let end = last.saturating_add(100).min(seq.len());
        if start > 0 {
            let _ = write!(text, "({start} bp)...");
        }
        text.push_str(&colored(seq, start, end, |pos| {
            if read
                .middle_hits
                .iter()
                .any(|(_, hit)| pos >= hit.start as usize && pos < hit.end as usize)
            {
                1
            } else if read.middle.iter().any(|&(a, b)| pos >= a && pos < b) {
                2
            } else {
                0
            }
        }));
        if end < seq.len() {
            let _ = write!(text, "...({} bp)", seq.len() - end);
        }
    }
    text
}
