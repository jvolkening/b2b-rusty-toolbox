//! Karlin-Altschul statistics for megablast scoring (reward 1, penalty -2,
//! linear gap costs), plus BLAST-compatible number formatting.

pub const LAMBDA: f64 = 1.28;
pub const K: f64 = 0.46;
const ALPHA: f64 = 1.5;
const BETA: f64 = -2.0;

/// Port of `BLAST_ComputeLengthAdjustment`.
pub fn length_adjustment(query_len: usize, db_len: u64, nseq: u64) -> i64 {
    let m = query_len as f64;
    let n = db_len as f64;
    let nn = nseq as f64;
    let alpha_d_lambda = ALPHA / LAMBDA;
    let log_k = K.ln();
    let a = nn;
    let mb = m * nn + n;
    let c = n * m - m.max(n) / K;
    if c < 0.0 {
        return 0;
    }
    let mut ell_max = 2.0 * c / (mb + (mb * mb - 4.0 * a * c).sqrt());
    let mut ell_min = 0.0f64;
    let mut ell_next = 0.0f64;
    let mut converged = false;
    for i in 1..=20 {
        let ell = ell_next;
        let ss = (m - ell) * (n - nn * ell);
        let ell_bar = alpha_d_lambda * (log_k + ss.ln()) + BETA;
        if ell_bar >= ell {
            ell_min = ell;
            if ell_bar - ell_min <= 1.0 {
                converged = true;
                break;
            }
            if ell_min == ell_max {
                break;
            }
        } else {
            ell_max = ell;
        }
        if ell_min <= ell_bar && ell_bar <= ell_max {
            ell_next = ell_bar;
        } else {
            ell_next = if i == 1 { ell_max } else { (ell_min + ell_max) / 2.0 };
        }
    }
    let mut adj = ell_min as i64;
    if converged {
        let ell = ell_min.ceil();
        if ell <= ell_max {
            let ss = (m - ell) * (n - nn * ell);
            if alpha_d_lambda * (log_k + ss.ln()) + BETA >= ell {
                adj = ell as i64;
            }
        }
    }
    adj
}

/// Per-query search-space data.
#[derive(Clone, Copy)]
pub struct SearchSpace {
    pub searchsp: f64,
}

impl SearchSpace {
    pub fn new(query_len: usize, db_len: u64, nseq: u64) -> Self {
        let adj = length_adjustment(query_len, db_len, nseq);
        let eff_q = (query_len as i64 - adj).max(1) as f64;
        let eff_db = (db_len as i64 - nseq as i64 * adj).max(1) as f64;
        SearchSpace { searchsp: eff_q * eff_db }
    }

    /// `score` is in raw (integer) BLAST units.
    pub fn evalue(&self, score: i32) -> f64 {
        self.searchsp * K * (-LAMBDA * score as f64).exp()
    }
}

pub fn bit_score(score: i32) -> f64 {
    (LAMBDA * score as f64 - K.ln()) / std::f64::consts::LN_2
}

/// `%.{prec}e` with a two-digit minimum exponent, like C printf.
fn sci(v: f64, prec: usize) -> String {
    let s = format!("{:.*e}", prec, v);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{m}e{}{:02}", if e < 0 { '-' } else { '+' }, e.abs())
}

/// E-value formatting as done by BLAST tabular output.
pub fn fmt_evalue(e: f64) -> String {
    if e < 1.0e-180 {
        "0.0".to_string()
    } else if e < 0.0009 {
        sci(e, 2)
    } else if e < 0.1 {
        format!("{:4.3}", e)
    } else if e < 1.0 {
        format!("{:3.2}", e)
    } else if e < 10.0 {
        format!("{:2.1}", e)
    } else {
        format!("{:2.0}", e)
    }
}

pub fn fmt_bits(b: f64) -> String {
    if b > 99999.0 {
        sci(b, 3)
    } else if b > 99.9 {
        format!("{:3}", b as i64)
    } else {
        format!("{:4.1}", b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formats() {
        assert_eq!(fmt_evalue(0.0), "0.0");
        assert_eq!(fmt_evalue(1.5e-30), "1.50e-30");
        assert_eq!(fmt_evalue(0.05), "0.050");
        assert_eq!(fmt_evalue(3.7), "3.7");
        assert_eq!(fmt_bits(1234.9), "1234");
        assert_eq!(fmt_bits(50.04), "50.0");
    }
}
