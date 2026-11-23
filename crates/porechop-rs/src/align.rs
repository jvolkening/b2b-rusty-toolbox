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
