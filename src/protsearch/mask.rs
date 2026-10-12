//! tantan repeat masking, as DIAMOND 2.1 applies it to queries and targets
//! (`masking/tantan.cpp`; Frith, NAR 2011 39:e23).
//!
//! Two-state-family HMM: a background state and repeat states for offsets
//! 1..50, repeat start probability 0.005, end probability 0.05, offset
//! weights growing ×(1/0.9) towards short periods, likelihood ratios
//! `exp(λ·BLOSUM62)` with λ = 0.3240321673 (the value tantan's
//! `LambdaCalculator` finds for BLOSUM62). A residue is replaced by `X` when
//! its posterior repeat probability is ≥ 0.9 (`--tantan-minMaskProb`).
//! Single-precision arithmetic and the 16-residue rescaling cadence follow
//! DIAMOND so masking decisions agree with it.

use super::alphabet::{score, NCODES, X};

const WINDOW: usize = 50;
const P_REPEAT: f32 = 0.005;
const P_REPEAT_END: f32 = 0.05;
const REPEAT_GROWTH: f32 = 1.0 / 0.9;
/// DIAMOND's `--tantan-minMaskProb` default.
pub const MIN_MASK_PROB: f32 = 0.9;
const LAMBDA: f64 = 0.324_032_167_348_044;

/// Likelihood-ratio table `exp(λ s(a, b))` (rows = a).
fn lr_table() -> &'static [[f32; NCODES]; NCODES] {
    use std::sync::OnceLock;
    static T: OnceLock<[[f32; NCODES]; NCODES]> = OnceLock::new();
    T.get_or_init(|| {
        let mut t = [[0f32; NCODES]; NCODES];
        for (a, row) in t.iter_mut().enumerate() {
            for (b, v) in row.iter_mut().enumerate() {
                *v = (LAMBDA * score(a as u8, b as u8) as f64).exp() as f32;
            }
        }
        t
    })
}

/// Mask repeats in `seq` (codes) in place; returns the number of residues
/// masked.
pub fn tantan_mask(seq: &mut [u8]) -> usize {
    let len = seq.len();
    if len == 0 {
        return 0;
    }
    let lr = lr_table();
    let b2b = 1.0f32 - P_REPEAT;
    let f2f = 1.0f32 - P_REPEAT_END;
    let b2f0 = P_REPEAT * (1.0 - REPEAT_GROWTH) / (1.0 - REPEAT_GROWTH.powf(WINDOW as f32));
    // d[k]: background -> repeat-with-offset-(k+1) transition
    let mut d = [0f32; WINDOW];
    d[WINDOW - 1] = b2f0;
    for i in (0..WINDOW - 1).rev() {
        d[i] = d[i + 1] * REPEAT_GROWTH;
    }
    // e(i, k) = LR(seq[i], seq[i-1-k]) or 0 before the start
    let emit = |i: usize, k: usize| -> f32 {
        if k < i {
            lr[(seq[i] as usize).min(NCODES - 1)][(seq[i - 1 - k] as usize).min(NCODES - 1)]
        } else {
            0.0
        }
    };

    let mut f = [0f32; WINDOW];
    let mut b = 1.0f32;
    let mut pb = vec![0f32; len];
    let mut scale = vec![0f32; (len - 1) / 16 + 1];
    for i in 0..len {
        let s: f32 = f.iter().sum();
        for k in 0..WINDOW {
            f[k] = (f[k] * f2f + b * d[k]) * emit(i, k);
        }
        b = b * b2b + s * P_REPEAT_END;
        if (i & 15) == 15 {
            let sc = 1.0 / b;
            scale[i / 16] = sc;
            b *= sc;
            for v in f.iter_mut() {
                *v *= sc;
            }
        }
        pb[i] = b;
    }
    let z = b * b2b + f.iter().sum::<f32>() * P_REPEAT_END;

    let mut b = b2b;
    let mut f = [P_REPEAT_END; WINDOW];
    let mut masked = Vec::new();
    for i in (0..len).rev() {
        let pf = 1.0 - (pb[i] * b / z);
        if (i & 15) == 15 {
            let sc = scale[i / 16];
            b *= sc;
            for v in f.iter_mut() {
                *v *= sc;
            }
        }
        for (k, v) in f.iter_mut().enumerate() {
            *v *= emit(i, k);
        }
        if pf >= MIN_MASK_PROB {
            masked.push(i);
        }
        let mut t = 0f32;
        for k in 0..WINDOW {
            t += f[k] * d[k];
        }
        for v in f.iter_mut() {
            *v = *v * f2f + P_REPEAT_END * b;
        }
        b = b2b * b + t;
    }
    for &i in &masked {
        seq[i] = X;
    }
    masked.len()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::protsearch::alphabet::{decode, encode};

    #[test]
    fn masks_tandem_repeats_only() {
        let unique = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQ";
        let mut s = unique.to_vec();
        for _ in 0..8 {
            s.extend_from_slice(b"PEPTIDES");
        }
        s.extend_from_slice(b"WHTYLRDNCAGFKEMI");
        let mut e = encode(&s);
        let n = tantan_mask(&mut e);
        let d = decode(&e);
        assert!(n >= 40, "masked {n}: {}", String::from_utf8_lossy(&d));
        // the unique prefix stays unmasked
        assert_eq!(&d[..unique.len() - 5], &unique[..unique.len() - 5]);
        // low complexity / poly-residue runs are masked too
        let mut p = encode(b"MKTAYIAKQRQISFVKSHFSRQQQQQQQQQQQQQQQQQQQQQQQQQQLEERLGLIEVQ");
        assert!(tantan_mask(&mut p) > 15);
    }

    #[test]
    fn leaves_ordinary_proteins_alone() {
        let mut e = encode(b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDW");
        assert_eq!(tantan_mask(&mut e), 0);
    }
}

/// The 1,000 abundant 8-mer motifs active in DIAMOND 2.1.9 (`masking/motifs.cpp`;
/// the rest of that file is commented out).
const MOTIFS: &str = include_str!("motifs.txt");

fn motif_set() -> &'static std::collections::HashSet<[u8; 8]> {
    use std::sync::OnceLock;
    static S: OnceLock<std::collections::HashSet<[u8; 8]>> = OnceLock::new();
    S.get_or_init(|| {
        MOTIFS
            .lines()
            .filter(|l| l.len() == 8)
            .map(|l| {
                let e = super::alphabet::encode(l.as_bytes());
                let mut k = [0u8; 8];
                k.copy_from_slice(&e);
                k
            })
            .collect()
    })
}

/// Motif soft masking for seeding (DIAMOND default sensitivity,
/// `mask_motifs`): occurrences of the motif 8-mers are merged into ranges;
/// unless they cover half the sequence or more, ranges of ≤ 30 residues are
/// replaced by `X`. Returns a copy; DIAMOND restores the letters after seed
/// enumeration, so this only removes seeds.
pub fn motif_soft_mask(seq: &[u8]) -> Vec<u8> {
    let mut out = seq.to_vec();
    if seq.len() < 8 {
        return out;
    }
    let set = motif_set();
    let mut ranges: Vec<(usize, usize)> = Vec::new();
    let mut k = [0u8; 8];
    for i in 0..=seq.len() - 8 {
        k.copy_from_slice(&seq[i..i + 8]);
        if set.contains(&k) {
            match ranges.last_mut() {
                Some(r) if i <= r.1 => r.1 = i + 8,
                _ => ranges.push((i, i + 8)),
            }
        }
    }
    let n: usize = ranges.iter().map(|r| r.1 - r.0).sum();
    if n * 2 >= seq.len() {
        return out;
    }
    for (a, b) in ranges {
        if b - a <= 30 {
            for c in &mut out[a..b] {
                *c = X;
            }
        }
    }
    out
}

#[cfg(test)]
mod motif_tests {
    use super::*;
    use crate::protsearch::alphabet::{decode, encode};

    #[test]
    fn walker_a_motif_is_soft_masked() {
        // 'GAGKSTLL' (P-loop) is in DIAMOND's motif table
        let s = encode(b"MKTAYIAKQRQISFVKSHFSRQLEERGAGKSTLLRLGLIEVQAPILSRVGDGTQDNLSGAEKAVQ");
        let m = motif_soft_mask(&s);
        let d = String::from_utf8(decode(&m)).unwrap();
        assert!(d.contains("XXXXXXXX"), "{d}");
        assert_eq!(m.len(), s.len());
        assert_eq!(motif_set().len(), 1000);
    }
}
