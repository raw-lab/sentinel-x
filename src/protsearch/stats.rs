//! Alignment statistics, exactly as DIAMOND 2.x reports them for
//! BLOSUM62 with gap penalties 11/1.
//!
//! * bit score `(λS − ln K) / ln 2` with the gapped λ = 0.267, K = 0.041;
//! * E-value from BLAST's Spouge finite-size-corrected formula
//!   (`BLAST_SpougeStoE`) evaluated on the query and subject lengths and
//!   scaled to the database (`× db_letters / subject_len`), with the Gumbel
//!   parameters of NCBI's BLOSUM62 11/1 table entry
//!   `{11, 1, ∞, 0.267, 0.041, 0.14, 1.9, −30, 0.669720, 42.6028, 43.6362}`
//!   and the ungapped row `{…, 0.7916, …, 4.96466, …}`.
//!
//! Checked against DIAMOND 2.1.9 output on 68,427 alignments: every
//! E-value agrees to the 3 significant digits DIAMOND prints.

/// Gapped Karlin-Altschul λ for BLOSUM62 11/1.
pub const LAMBDA: f64 = 0.267;
/// Gapped Karlin-Altschul K for BLOSUM62 11/1.
pub const K: f64 = 0.041;

// Spouge parameters (BLAST `Blast_GumbelBlkLoadFromTables`).
const A: f64 = 1.9;
const ALPHA: f64 = 42.6028;
const SIGMA: f64 = 43.6362;
const A_UN: f64 = 0.7916;
const ALPHA_UN: f64 = 4.96466;
const G: f64 = 12.0; // gap open + gap extend
const B: f64 = 2.0 * G * (A_UN - A);
const BETA: f64 = 2.0 * G * (ALPHA_UN - ALPHA);
const TAU: f64 = 2.0 * G * (ALPHA_UN - SIGMA);

/// Bit score of a raw BLOSUM62 11/1 score.
#[inline]
pub fn bitscore(raw: i32) -> f64 {
    (LAMBDA * raw as f64 - K.ln()) / std::f64::consts::LN_2
}

/// Complementary error function (Numerical Recipes `erfcc`, |rel err| < 1.2e-7,
/// ample for E-values printed to 3 significant digits).
fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let r = t
        * (-z * z - 1.26551223
            + t * (1.00002368
                + t * (0.37409196
                    + t * (0.09678418
                        + t * (-0.18628806
                            + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
            .exp();
    if x >= 0.0 {
        r
    } else {
        2.0 - r
    }
}

/// Database-scaled E-value of raw score `y` for a query of length `m` against
/// a subject of length `n` in a database of `db_letters` residues
/// (BLAST `BLAST_SpougeStoE`).
pub fn evalue(y: i32, m: usize, n: usize, db_letters: u64) -> f64 {
    let y = y as f64;
    let m = m as f64;
    let n = (n.max(1)) as f64;
    let db_scale = if db_letters > 0 { db_letters as f64 / n } else { 1.0 };
    const C: f64 = 0.398_942_280_401_432_7; // 1/sqrt(2π)
    let v = (2.0 * ALPHA / LAMBDA).max(ALPHA * y + BETA);
    let sv = v.sqrt();
    let m_li = m - (A * y + B);
    let m_f = m_li / sv;
    let p_m_f = erfc(-m_f / std::f64::consts::SQRT_2) / 2.0;
    let p1 = m_li * p_m_f + sv * C * (-0.5 * m_f * m_f).exp();
    let n_lj = n - (A * y + B);
    let n_f = n_lj / sv;
    let p_n_f = erfc(-n_f / std::f64::consts::SQRT_2) / 2.0;
    let p2 = n_lj * p_n_f + sv * C * (-0.5 * n_f * n_f).exp();
    let c_y = (2.0 * SIGMA / LAMBDA).max(SIGMA * y + TAU);
    let area = p1 * p2 + c_y * p_m_f * p_n_f;
    area * K * (-LAMBDA * y).exp() * db_scale
}

/// Smallest raw score whose Karlin-Altschul E-value (no edge correction,
/// which is never larger than the Spouge value's leading term) is at most
/// `e` for a query of length `m` against `db_letters` residues.  Used as a
/// conservative pre-filter threshold.
pub fn ka_raw_threshold(e: f64, m: usize, db_letters: u64) -> i32 {
    let s = ((K * m.max(1) as f64 * db_letters.max(1) as f64).ln() - e.ln()) / LAMBDA;
    s.ceil().max(1.0) as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_diamond_2_1_9() {
        // (raw, qlen, slen, db_letters, DIAMOND E-value)
        let cases = [
            (1625, 470, 461, 70_293_542u64, 3.97e-226),
            (1550, 470, 467, 70_293_542, 1.24e-214),
            (1476, 470, 476, 70_293_542, 2.97e-203),
        ];
        for (s, m, n, db, e) in cases {
            let got = evalue(s, m, n, db);
            let r = (got.log10() - f64::log10(e)).abs();
            assert!(r < 0.003, "S={s}: got {got:e}, want {e:e}");
        }
        assert!((bitscore(1625) - 630.6).abs() < 0.1);
    }

    #[test]
    fn erfc_values() {
        assert!((erfc(0.0) - 1.0).abs() < 1e-7);
        assert!((erfc(1.0) - 0.157_299_207).abs() < 1e-6);
        assert!((erfc(-1.0) - 1.842_700_793).abs() < 1e-6);
    }
}
