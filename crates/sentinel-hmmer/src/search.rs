//! High-level search: score a set of (configured) profiles against a set of
//! protein sequences, compute E-values from the calibrated `STATS` parameters,
//! and apply bit-score or gathering-cutoff thresholds. This is the surface the
//! CheckV / CheckM2 / BUSCO / viralVerify modules call.

use rayon::prelude::*;

use crate::dp::score_seq;
use crate::model::Profile;

/// A single profile→sequence hit above threshold.
#[derive(Clone, Debug)]
pub struct Hit {
    /// Profile name (e.g. Pfam `PF00069.27` or BUSCO `OG…`).
    pub profile: String,
    /// Target sequence id.
    pub target: String,
    /// Forward bit score.
    pub bits: f32,
    /// Viterbi bit score.
    pub viterbi_bits: f32,
    /// E-value (per the calibrated tail), or `f64::NAN` if uncalibrated.
    pub evalue: f64,
    /// Best-alignment envelope on the target (1-based inclusive).
    pub env: Option<(usize, usize)>,
}

/// Thresholding policy for accepting hits.
#[derive(Clone, Copy, Debug)]
pub enum Threshold {
    /// Keep hits with E-value ≤ this (requires calibrated profiles).
    Evalue(f64),
    /// Keep hits with Forward bit score ≥ this.
    Bits(f32),
    /// Keep hits passing the per-profile `GA` gathering cutoff (Pfam style).
    GatheringCutoff,
}

/// Number of comparisons used to scale E-values (database size in residues is
/// folded into the calibration; HMMER multiplies the per-comparison P-value by
/// the number of target sequences).
fn evalue_forward(bits: f32, stats: &crate::model::Stats, n_targets: usize) -> f64 {
    if let Some((tau, lambda)) = stats.forward {
        // Exponential right tail: P(s > x) = exp(-lambda (x - tau)).
        let p = (-(lambda) * (bits as f64 - tau)).exp();
        (n_targets as f64) * p.min(1.0)
    } else {
        f64::NAN
    }
}

fn evalue_viterbi(bits: f32, stats: &crate::model::Stats, n_targets: usize) -> f64 {
    if let Some((mu, lambda)) = stats.viterbi {
        // Gumbel right tail: P(s > x) = 1 - exp(-exp(-lambda (x - mu))).
        let p = 1.0 - (-((-lambda * (bits as f64 - mu)).exp())).exp();
        (n_targets as f64) * p.min(1.0)
    } else {
        f64::NAN
    }
}

/// Search one profile against many sequences. `seqs` is `(id, residues)`.
pub fn search_profile(
    profile: &Profile,
    seqs: &[(String, Vec<u8>)],
    threshold: Threshold,
) -> Vec<Hit> {
    let n = seqs.len();
    seqs.par_iter()
        .filter_map(|(id, seq)| {
            let s = score_seq(profile, seq);
            let bits = s.forward_bits;
            let vbits = s.viterbi_bits;
            let ev = evalue_forward(bits, &profile.stats, n);
            let keep = match threshold {
                Threshold::Evalue(t) => ev.is_finite() && ev <= t,
                Threshold::Bits(t) => bits >= t,
                Threshold::GatheringCutoff => match profile.ga {
                    // GA holds (seq, domain) cutoffs; use the sequence cutoff.
                    Some((seq_ga, _dom_ga)) => bits as f64 >= seq_ga,
                    None => false,
                },
            };
            if keep {
                Some(Hit {
                    profile: profile.name.clone(),
                    target: id.clone(),
                    bits,
                    viterbi_bits: vbits,
                    evalue: ev,
                    env: s.env,
                })
            } else {
                None
            }
        })
        .collect()
}

/// Search many profiles against many sequences. Profiles are iterated serially
/// (each already parallelises over sequences) to bound peak memory on the
/// large Pfam/VOG databases.
pub fn search_all(
    profiles: &[Profile],
    seqs: &[(String, Vec<u8>)],
    threshold: Threshold,
) -> Vec<Hit> {
    let mut hits = Vec::new();
    for p in profiles {
        hits.extend(search_profile(p, seqs, threshold));
    }
    hits
}

/// Convenience: expose the Viterbi-based E-value too (used by callers that want
/// the MSV/Viterbi filtering semantics rather than Forward).
pub fn viterbi_evalue(bits: f32, profile: &Profile, n_targets: usize) -> f64 {
    evalue_viterbi(bits, &profile.stats, n_targets)
}
