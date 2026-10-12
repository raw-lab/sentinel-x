//! Composition bias correction ("Hauser 2016"), DIAMOND 2.1's default
//! `--comp-based-stats 1` (`stats/hauser_correction.cpp`).
//!
//! For each query position i with a standard residue r, the correction is
//! the expected BLOSUM62 score of r against the background minus its mean
//! score against the other residues in a ±20 window (window = 40):
//! `c_i = bg_score(r) − (Σ_window s(r, x_j) − s(r, r)) / (n − 1)`.
//! DIAMOND adds `round(c_i)` (int8) to every substitution score at query
//! position i whose target residue is a standard amino acid.

use super::alphabet::{score, BACKGROUND, TRUE_AA};

/// DIAMOND's `--cbs-window` default.
const WINDOW: usize = 40;

fn background_scores() -> [f64; TRUE_AA] {
    let mut out = [0f64; TRUE_AA];
    for (i, o) in out.iter_mut().enumerate() {
        for (j, &p) in BACKGROUND.iter().enumerate() {
            *o += p * score(i as u8, j as u8) as f64;
        }
    }
    out
}

/// Per-position corrections (float) for a (masked) query.
pub fn hauser_float(seq: &[u8]) -> Vec<f32> {
    let l = seq.len();
    let mut out = vec![0f32; l];
    if l == 0 {
        return out;
    }
    let bg = background_scores();
    let mut scores = [0i32; TRUE_AA];
    let add = |scores: &mut [i32; TRUE_AA], c: u8| {
        for (i, s) in scores.iter_mut().enumerate() {
            *s += score(c, i as u8);
        }
    };
    let sub = |scores: &mut [i32; TRUE_AA], c: u8| {
        for (i, s) in scores.iter_mut().enumerate() {
            *s -= score(c, i as u8);
        }
    };
    let window_half = (WINDOW / 2).min(l - 1);
    let (mut n, mut h, mut m, mut t) = (0usize, 0usize, 0usize, 0usize);
    let set = |out: &mut [f32], scores: &[i32; TRUE_AA], m: usize, n: usize| {
        let r = seq[m] as usize;
        if r < TRUE_AA {
            out[m] = bg[r] as f32 - (scores[r] - score(seq[m], seq[m])) as f32 / (n as f32 - 1.0);
        }
    };
    while n < window_half && h < l {
        n += 1;
        add(&mut scores, seq[h]);
        h += 1;
    }
    while n < WINDOW + 1 && h < l {
        n += 1;
        add(&mut scores, seq[h]);
        set(&mut out, &scores, m, n);
        h += 1;
        m += 1;
    }
    while h < l {
        add(&mut scores, seq[h]);
        sub(&mut scores, seq[t]);
        set(&mut out, &scores, m, n);
        h += 1;
        t += 1;
        m += 1;
    }
    while m < l && n > window_half + 1 {
        n -= 1;
        sub(&mut scores, seq[t]);
        set(&mut out, &scores, m, n);
        t += 1;
        m += 1;
    }
    while m < l {
        set(&mut out, &scores, m, n);
        m += 1;
    }
    out
}

/// Rounded (int8) corrections as DIAMOND's SIMD profiles use them.
pub fn hauser_int8(seq: &[u8]) -> Vec<i8> {
    hauser_float(seq).into_iter().map(|f| (if f < 0.0 { f - 0.5 } else { f + 0.5 }) as i8).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::protsearch::alphabet::encode;

    #[test]
    fn biased_composition_gets_negative_correction() {
        // a K/E-rich region: K scores well against K/E/R, so its local mean
        // exceeds the background expectation and the correction is negative
        let s = encode(b"MKTAYIAKQRQISFVKSHFSRQKEKEKKEEKKEKEKKEKEEKKEKEKEKEKLEERLGLIEVQAPILSRVGDGTQDNL");
        let c = hauser_float(&s);
        let mid = 40;
        assert!(c[mid] < -0.5, "{}", c[mid]);
        let i8s = hauser_int8(&s);
        assert_eq!(i8s.len(), s.len());
        assert!(i8s[mid] <= -1);
    }

    #[test]
    fn short_sequences_do_not_panic() {
        for n in 1..5 {
            let s = encode(&b"MKTAY"[..n]);
            let c = hauser_float(&s);
            assert_eq!(c.len(), n);
        }
    }
}
