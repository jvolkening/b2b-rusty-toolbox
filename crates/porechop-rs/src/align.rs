#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct Alignment {
    pub start: i32,
    pub end: i32,
    pub partial: f64,
    pub full: f64,
}

unsafe extern "C" {
    fn porechop_align(
        read: *const u8,
        read_len: usize,
        adapter: *const u8,
        adapter_len: usize,
        scores: *const i32,
    ) -> Alignment;
    fn porechop_middle_align(
        read: *const u8,
        read_len: usize,
        adapter: *const u8,
        adapter_len: usize,
        backend: i32,
        result: *mut Alignment,
    ) -> i32;
}

pub fn align_middle(
    read: &[u8],
    adapter: &[u8],
    scores: &[i32; 4],
    threshold: f64,
) -> Result<Alignment> {
    if scores != &[3, -6, -5, -2] || threshold < 90.0 || adapter.len() < 9 {
        return Ok(align(read, adapter, scores));
    }
    let mut result = Alignment {
        start: -1,
        end: 0,
        partial: 0.0,
        full: 0.0,
    };
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    let backend = if std::is_x86_feature_detected!("avx2") {
        2
    } else if std::is_x86_feature_detected!("sse2") {
        1
    } else {
        0
    };
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    let backend = 0;
    let success = unsafe {
        porechop_middle_align(
            read.as_ptr(),
            read.len(),
            adapter.as_ptr(),
            adapter.len(),
            backend,
            &mut result,
        )
    };
    if success == 0 {
        bail!("Parasail middle-adapter alignment failed");
    }
    result.partial = (result.partial * 1_000_000.0).round() / 1_000_000.0;
    result.full = (result.full * 1_000_000.0).round() / 1_000_000.0;
    Ok(result)
}

pub fn align(read: &[u8], adapter: &[u8], scores: &[i32; 4]) -> Alignment {
    // Both input slices remain alive for the duration of this synchronous call.
    let mut result = unsafe {
        porechop_align(
            read.as_ptr(),
            read.len(),
            adapter.as_ptr(),
            adapter.len(),
            scores.as_ptr(),
        )
    };
    // The original C++/Python interface serializes identities to six decimal places.
    result.partial = (result.partial * 1_000_000.0).round() / 1_000_000.0;
    result.full = (result.full * 1_000_000.0).round() / 1_000_000.0;
    result
}
use anyhow::{Result, bail};

#[cfg(test)]
mod tests {
    use super::*;

    const ADAPTER: &[u8] = b"AATGTACTTCGTTCAGTTACGTATTGCT";
    const SCORES: [i32; 4] = [3, -6, -5, -2];

    #[test]
    fn middle_search_finds_and_masks_multiple_hits() {
        let mut read =
            b"CCCCCaAtGtAcTtCgTtCaGtTaCgTaTtGcTGGGGGAATGTACTTCGTTCAGTTACGTATTGCTCCCCC".to_vec();
        let first = align_middle(&read, ADAPTER, &SCORES, 90.0).unwrap();
        assert_eq!((first.start, first.end, first.full), (5, 33, 100.0));
        read[first.start as usize..first.end as usize].fill(b'-');
        let second = align_middle(&read, ADAPTER, &SCORES, 90.0).unwrap();
        assert_eq!((second.start, second.end, second.full), (38, 66, 100.0));
        read[second.start as usize..second.end as usize].fill(b'-');
        assert_eq!(
            align_middle(&read, ADAPTER, &SCORES, 90.0).unwrap().start,
            -1
        );
    }

    #[test]
    fn lower_threshold_and_custom_scores_use_unseeded_alignment() {
        let mut disrupted = ADAPTER.to_vec();
        for i in [1, 8, 15, 22] {
            disrupted[i] = if disrupted[i] == b'C' { b'A' } else { b'C' };
        }
        let mut read = b"GGGGG".to_vec();
        read.extend_from_slice(&disrupted);
        read.extend_from_slice(b"GGGGG");
        assert_eq!(
            align_middle(&read, ADAPTER, &SCORES, 90.0).unwrap().start,
            -1
        );
        let expected = align(&read, ADAPTER, &SCORES);
        let fallback = align_middle(&read, ADAPTER, &SCORES, 85.0).unwrap();
        assert!(fallback.full >= 85.0);
        assert_eq!(
            (fallback.start, fallback.end, fallback.full),
            (expected.start, expected.end, expected.full)
        );
        let scores = [4, -6, -5, -2];
        let expected = align(&read, ADAPTER, &scores);
        let fallback = align_middle(&read, ADAPTER, &scores, 90.0).unwrap();
        assert_eq!(
            (fallback.start, fallback.end, fallback.full),
            (expected.start, expected.end, expected.full)
        );
    }

    #[test]
    fn short_adapters_and_empty_reads_use_safe_paths() {
        let short = b"ACGTACGT";
        let expected = align(b"GGACGTACGTCC", short, &SCORES);
        let fallback = align_middle(b"GGACGTACGTCC", short, &SCORES, 90.0).unwrap();
        assert_eq!(
            (fallback.start, fallback.end, fallback.full),
            (expected.start, expected.end, expected.full)
        );
        assert_eq!(align_middle(b"", ADAPTER, &SCORES, 90.0).unwrap().start, -1);
    }

    #[test]
    fn available_native_backends_agree_on_exact_alignment() {
        let read = b"GGGGGAATGTACTTCGTTCAGTTACGTATTGCTCCCCC";
        #[allow(unused_mut)]
        let mut backends = vec![0];
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if std::is_x86_feature_detected!("sse2") {
                backends.push(1);
            }
            if std::is_x86_feature_detected!("avx2") {
                backends.push(2);
            }
        }
        for backend in backends {
            let mut result = Alignment {
                start: -1,
                end: 0,
                partial: 0.0,
                full: 0.0,
            };
            let status = unsafe {
                porechop_middle_align(
                    read.as_ptr(),
                    read.len(),
                    ADAPTER.as_ptr(),
                    ADAPTER.len(),
                    backend,
                    &mut result,
                )
            };
            assert_eq!(status, 1);
            assert_eq!((result.start, result.end, result.full), (5, 33, 100.0));
        }
    }
}
