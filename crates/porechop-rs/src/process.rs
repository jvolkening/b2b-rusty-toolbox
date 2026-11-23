use crate::adapters::Adapter;
use crate::align::align;
use crate::reads::ReadRecord;
use anyhow::{Result, bail};
use rayon::prelude::*;
use std::collections::HashSet;

pub struct Settings {
    pub verbosity: u8,
    pub end_size: usize,
    pub end_threshold: f64,
    pub min_trim_size: usize,
    pub extra_end_trim: usize,
    pub middle_threshold: f64,
    pub extra_middle_trim_good_side: usize,
    pub extra_middle_trim_bad_side: usize,
    pub barcode_threshold: f64,
    pub barcode_diff: f64,
    pub require_two_barcodes: bool,
    pub scores: [i32; 4],
}

pub fn discover(adapters: &mut [Adapter], reads: &[ReadRecord], indices: &[usize], s: &Settings) {
    if adapters.len() < rayon::current_num_threads() {
        for adapter in adapters {
            if adapter.name.contains("(full sequence)") {
                continue;
            }
            let (start, end) = indices
                .par_iter()
                .map(|&idx| {
                    let seq = &reads[idx].seq;
                    let start = adapter.start.as_ref().map_or(0.0, |a| {
                        align(
                            &seq[..seq.len().min(s.end_size)],
                            a.sequence.as_bytes(),
                            &s.scores,
                        )
                        .full
                    });
                    let end = adapter.end.as_ref().map_or(0.0, |a| {
                        align(
                            &seq[seq.len().saturating_sub(s.end_size)..],
                            a.sequence.as_bytes(),
                            &s.scores,
                        )
                        .full
                    });
                    (start, end)
                })
                .reduce(|| (0.0, 0.0), |a, b| (a.0.max(b.0), a.1.max(b.1)));
            adapter.best_start = start;
            adapter.best_end = end;
        }
        return;
    }
    adapters.par_iter_mut().for_each(|adapter| {
        if adapter.name.contains("(full sequence)") {
            return;
        }
        for &idx in indices {
            let seq = &reads[idx].seq;
            if let Some(start) = &adapter.start {
                let score = align(
                    &seq[..seq.len().min(s.end_size)],
                    start.sequence.as_bytes(),
                    &s.scores,
                )
                .full;
                adapter.best_start = adapter.best_start.max(score);
            }
            if let Some(end) = &adapter.end {
                let score = align(
                    &seq[seq.len().saturating_sub(s.end_size)..],
                    end.sequence.as_bytes(),
                    &s.scores,
                )
                .full;
                adapter.best_end = adapter.best_end.max(score);
            }
        }
    });
}

pub fn fix_up_1d2(adapters: &mut Vec<Adapter>) {
    let score = |name: &str| {
        adapters
            .iter()
            .find(|a| a.name == name)
            .map(|a| a.best_start.max(a.best_end))
    };
    if let (Some(first), Some(second), Some(short)) = (
        score("1D^2 part 1"),
        score("1D^2 part 2"),
        score("SQK-MAP006 Short"),
    ) && first >= short
        && second >= short
    {
        adapters.retain(|a| a.name != "SQK-MAP006 Short");
    }
}

pub fn choose_orientation(adapters: &[Adapter]) -> Result<&'static str> {
    let (mut f_or, mut r_or, mut f_and, mut r_and) = (0.0, 0.0, 0.0, 0.0);
    for a in adapters {
        let name = a.name.to_lowercase();
        if name.contains("barcode") && name.contains("(forward)") {
            f_or += a.best_start.max(a.best_end);
            f_and += a.best_start + a.best_end;
        } else if name.contains("barcode") && name.contains("(reverse)") {
            r_or += a.best_start.max(a.best_end);
            r_and += a.best_start + a.best_end;
        }
    }
    if f_or == 0.0 && r_or == 0.0 {
        bail!("no barcodes were found, so Porechop cannot perform barcode demultiplexing");
    }
    if f_or > r_or || (f_or == r_or && f_and > r_and) {
        Ok("forward")
    } else if r_or > f_or || (f_or == r_or && r_and > f_and) {
        Ok("reverse")
    } else {
        bail!("Porechop could not determine barcode orientation");
    }
}

fn update(scores: &mut Vec<(String, f64)>, name: String, score: f64) {
    if let Some(item) = scores.iter_mut().find(|item| item.0 == name) {
        item.1 = score;
    } else {
        scores.push((name, score));
    }
}

fn call_barcode(
    read: &mut ReadRecord,
    mut start: Vec<(String, f64)>,
    mut end: Vec<(String, f64)>,
    s: &Settings,
) {
    start.sort_by(|a, b| b.1.total_cmp(&a.1));
    end.sort_by(|a, b| b.1.total_cmp(&a.1));
    let none = ("none".to_owned(), 0.0);
    let best_start = start.first().unwrap_or(&none);
    let best_end = end.first().unwrap_or(&none);
    let second_start = start.get(1).unwrap_or(&none);
    let second_end = end.get(1).unwrap_or(&none);
    if s.require_two_barcodes {
        if best_start.1 >= s.barcode_threshold
            && best_end.1 >= s.barcode_threshold
            && best_start.1 >= second_start.1 + s.barcode_diff
            && best_end.1 >= second_end.1 + s.barcode_diff
            && best_start.0 == best_end.0
        {
            read.barcode = best_start.0.clone();
        }
    } else {
        start.extend(end);
        start.sort_by(|a, b| b.1.total_cmp(&a.1));
        let mut seen = HashSet::new();
        start.retain(|(name, _)| seen.insert(name.clone()));
        let best = start.first().unwrap_or(&none);
        let second = start.get(1).unwrap_or(&none);
        if best.1 >= s.barcode_threshold && best.1 >= second.1 + s.barcode_diff {
            read.barcode = best.0.clone();
        }
    }
    if read
        .albacore
        .as_ref()
        .is_some_and(|call| *call != read.barcode)
    {
        read.barcode = "none".into();
    }
}

pub fn trim_ends(
    reads: &mut [ReadRecord],
    adapters: &[Adapter],
    s: &Settings,
    orientation: Option<&str>,
) {
    reads.par_iter_mut().for_each(|read| {
        let start_seq = &read.seq[..read.seq.len().min(s.end_size)];
        let end_seq = &read.seq[read.seq.len().saturating_sub(s.end_size)..];
        let mut start_scores = Vec::new();
        let mut end_scores = Vec::new();
        for adapter in adapters {
            if let Some(start) = &adapter.start {
                let hit = align(start_seq, start.sequence.as_bytes(), &s.scores);
                if hit.partial > s.end_threshold
                    && hit.end != s.end_size as i32
                    && hit.end - hit.start >= s.min_trim_size as i32
                {
                    read.start_trim = read
                        .start_trim
                        .max((hit.end as usize).saturating_add(s.extra_end_trim));
                    if s.verbosity >= 3 {
                        read.start_hits.push((adapter.name.clone(), hit));
                    }
                }
                if orientation
                    .is_some_and(|o| adapter.barcode() && adapter.barcode_direction() == o)
                {
                    update(&mut start_scores, adapter.barcode_name(), hit.full);
                }
            }
            if let Some(end) = &adapter.end {
                let hit = align(end_seq, end.sequence.as_bytes(), &s.scores);
                if hit.partial > s.end_threshold
                    && hit.start > 0
                    && hit.end - hit.start >= s.min_trim_size as i32
                {
                    read.end_trim = read.end_trim.max(
                        s.end_size
                            .saturating_sub(hit.start as usize)
                            .saturating_add(s.extra_end_trim),
                    );
                    if s.verbosity >= 3 {
                        read.end_hits.push((adapter.name.clone(), hit));
                    }
                }
                if orientation
                    .is_some_and(|o| adapter.barcode() && adapter.barcode_direction() == o)
                {
                    update(&mut end_scores, adapter.barcode_name(), hit.full);
                }
            }
        }
        if s.verbosity >= 2 && orientation.is_some() {
            read.start_scores = start_scores.clone();
            read.end_scores = end_scores.clone();
        }
        if orientation.is_some() {
            call_barcode(read, start_scores, end_scores, s);
        }
    });
}

pub fn trim_middle(reads: &mut [ReadRecord], adapters: &[Adapter], s: &Settings) {
    let mut sequences = Vec::new();
    let mut start_names = HashSet::new();
    let mut end_names = HashSet::new();
    for adapter in adapters {
        if let Some(start) = &adapter.start {
            sequences.push(start);
            start_names.insert(start.name.as_str());
        }
        if let Some(end) = &adapter.end {
            if adapter
                .start
                .as_ref()
                .is_none_or(|start| start.sequence != end.sequence)
            {
                sequences.push(end);
            }
            end_names.insert(end.name.as_str());
        }
    }
    reads.par_iter_mut().for_each(|read| {
        let mut masked = read.trimmed().to_vec();
        for sequence in &sequences {
            loop {
                let hit = align(&masked, sequence.sequence.as_bytes(), &s.scores);
                if hit.full < s.middle_threshold || hit.start < 0 || hit.end <= hit.start {
                    break;
                }
                let start = hit.start as usize;
                let end = (hit.end as usize).min(masked.len());
                if start >= end || masked[start..end].iter().all(|&b| b == b'-') {
                    break;
                }
                masked[start..end].fill(b'-');
                if s.verbosity >= 2 {
                    read.middle_hits.push((sequence.name.clone(), hit));
                }
                let left = if start_names.contains(sequence.name.as_str()) {
                    s.extra_middle_trim_bad_side
                } else {
                    s.extra_middle_trim_good_side
                };
                let right = if end_names.contains(sequence.name.as_str()) {
                    s.extra_middle_trim_bad_side
                } else {
                    s.extra_middle_trim_good_side
                };
                read.middle
                    .push((start.saturating_sub(left), end.saturating_add(right)));
            }
        }
    });
}
