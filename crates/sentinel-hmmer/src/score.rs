//! Numerically-stable log-space arithmetic for the Plan7 dynamic programming
//! recurrences. All profile scores in Sentinel's HMMER engine are kept in
//! **bits** (log base 2) of log-odds relative to a background null model, which
//! is the same convention HMMER3 reports.

/// Negative infinity sentinel used for impossible states / `*` transitions.
pub const NEG_INF: f32 = f32::NEG_INFINITY;

/// `log2(exp(a) + exp(b))` computed without overflow.
///
/// This is the inner kernel of the Forward algorithm; it is called billions of
/// times on large profiles, so it is kept branch-light.
#[inline(always)]
pub fn logsum2(a: f32, b: f32) -> f32 {
    if a == NEG_INF {
        return b;
    }
    if b == NEG_INF {
        return a;
    }
    let (hi, lo) = if a > b { (a, b) } else { (b, a) };
    // log2(2^hi + 2^lo) = hi + log2(1 + 2^(lo-hi))
    hi + (1.0 + exp2(lo - hi)).log2()
}

/// Fold a slice of log-space values into a single `logsum`.
#[inline]
pub fn logsum_slice(xs: &[f32]) -> f32 {
    let mut acc = NEG_INF;
    for &x in xs {
        acc = logsum2(acc, x);
    }
    acc
}

/// `2^x`. Pulled out so a faster approximation can be dropped in later without
/// touching call sites.
#[inline(always)]
pub fn exp2(x: f32) -> f32 {
    x.exp2()
}

/// Convert a probability to a log2 score, mapping `0` to [`NEG_INF`].
#[inline]
pub fn prob_to_log2(p: f64) -> f32 {
    if p <= 0.0 {
        NEG_INF
    } else {
        p.log2() as f32
    }
}

/// HMMER3 stores transition/emission parameters as `-ln(prob)`, with the token
/// `*` meaning probability zero. This converts that on-disk value to a log2
/// probability for the DP.
#[inline]
pub fn neg_ln_to_log2(neg_ln: f32) -> f32 {
    if neg_ln.is_infinite() {
        NEG_INF
    } else {
        // p = exp(-neg_ln);  log2(p) = -neg_ln / ln(2)
        -neg_ln / std::f32::consts::LN_2
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn logsum_matches_naive() {
        let a = -3.0f32;
        let b = -1.0f32;
        let naive = (exp2(a) + exp2(b)).log2();
        assert!((logsum2(a, b) - naive).abs() < 1e-5);
    }

    #[test]
    fn neg_inf_absorbs() {
        assert_eq!(logsum2(NEG_INF, -2.0), -2.0);
        assert_eq!(logsum2(-2.0, NEG_INF), -2.0);
        assert_eq!(logsum2(NEG_INF, NEG_INF), NEG_INF);
    }
}
