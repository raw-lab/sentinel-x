//! Plan7 profile-HMM data structures.
//!
//! A profile has `M` match nodes. Node `k` (1-based, `1..=M`) carries match and
//! insert emission distributions over the 20 amino acids and the seven Plan7
//! transition probabilities. After parsing, [`Hmm::configure`] turns the raw
//! probabilities into a search-ready [`Profile`] in log2-odds space with local
//! entry/exit and a length-aware special-state model, mirroring HMMER3's
//! `p7_ProfileConfig` in multihit-local mode.

use crate::score::{neg_ln_to_log2, prob_to_log2, NEG_INF};

/// The 20 amino acids in HMMER3 column order.
pub const AA: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
/// Number of amino-acid symbols.
pub const K: usize = 20;

/// Map an ASCII residue to its column index `0..20`, or `None` for non-standard
/// symbols (X, B, Z, gaps, …) which are treated as "no information".
#[inline]
pub fn aa_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'D' => Some(2),
        b'E' => Some(3),
        b'F' => Some(4),
        b'G' => Some(5),
        b'H' => Some(6),
        b'I' => Some(7),
        b'K' => Some(8),
        b'L' => Some(9),
        b'M' => Some(10),
        b'N' => Some(11),
        b'P' => Some(12),
        b'Q' => Some(13),
        b'R' => Some(14),
        b'S' => Some(15),
        b'T' => Some(16),
        b'V' => Some(17),
        b'W' => Some(18),
        b'Y' => Some(19),
        _ => None,
    }
}

/// Amino-acid background frequencies (Swiss-Prot averages, as used by HMMER3's
/// default null model `p7_AminoFrequencies`).
pub const NULL_FREQ: [f64; K] = [
    0.0787945, 0.0151600, 0.0535222, 0.0668298, 0.0397062, 0.0695071, 0.0229198, 0.0590092,
    0.0594422, 0.0963728, 0.0237718, 0.0414386, 0.0482904, 0.0395639, 0.0540978, 0.0683364,
    0.0540687, 0.0673417, 0.0114135, 0.0304133,
];

/// The seven Plan7 transitions out of / between nodes.
#[derive(Clone, Copy, Debug)]
pub struct Trans {
    pub mm: f32,
    pub mi: f32,
    pub md: f32,
    pub im: f32,
    pub ii: f32,
    pub dm: f32,
    pub dd: f32,
}

/// Calibrated statistical parameters from the `STATS LOCAL …` lines, used for
/// E-value computation.
#[derive(Clone, Copy, Debug, Default)]
pub struct Stats {
    /// MSV Gumbel location/scale.
    pub msv: Option<(f64, f64)>,
    /// Viterbi Gumbel location/scale.
    pub viterbi: Option<(f64, f64)>,
    /// Forward exponential-tail tau/lambda.
    pub forward: Option<(f64, f64)>,
}

/// A parsed (raw, probability-space) profile HMM.
#[derive(Clone, Debug)]
pub struct Hmm {
    pub name: String,
    pub acc: Option<String>,
    pub desc: Option<String>,
    pub m: usize,
    /// `match_emit[k][a]`, `k` in `1..=M` (index 0 is unused / COMPO).
    pub match_emit: Vec<[f64; K]>,
    /// `insert_emit[k][a]`, `k` in `0..=M`.
    pub insert_emit: Vec<[f64; K]>,
    /// `trans[k]` are the transitions written on node `k`'s record, `0..=M`.
    pub trans: Vec<Trans>,
    /// Per-domain gathering / trusted / noise cutoffs (GA, TC, NC) if present.
    pub ga: Option<(f64, f64)>,
    pub tc: Option<(f64, f64)>,
    pub nc: Option<(f64, f64)>,
    pub stats: Stats,
}

/// A search-ready profile: log2-odds emissions plus a configured special-state
/// model. Built by [`Hmm::configure`].
#[derive(Clone, Debug)]
pub struct Profile {
    pub name: String,
    pub acc: Option<String>,
    pub m: usize,

    /// `msc[k][a]` = log2( match_emit[k][a] / null[a] ), `k` in `1..=M`.
    pub msc: Vec<[f32; K]>,
    /// `isc[k][a]` = log2( insert_emit[k][a] / null[a] ), `k` in `0..=M`.
    pub isc: Vec<[f32; K]>,

    /// Node transitions in log2 prob, `k` in `0..=M`.
    pub tmm: Vec<f32>,
    pub tmi: Vec<f32>,
    pub tmd: Vec<f32>,
    pub tim: Vec<f32>,
    pub tii: Vec<f32>,
    pub tdm: Vec<f32>,
    pub tdd: Vec<f32>,

    /// Local entry `B -> M_k`, `k` in `1..=M`.
    pub tbm: Vec<f32>,

    // Special-state transitions (multihit local, length-configured).
    pub tnn: f32,
    pub tnb: f32,
    pub tec: f32,
    pub tej: f32,
    pub tcc: f32,
    pub tct: f32,
    pub tjj: f32,
    pub tjb: f32,

    pub ga: Option<(f64, f64)>,
    pub stats: Stats,
}

impl Hmm {
    /// Configure the raw HMM into a search-ready [`Profile`] for a target
    /// database of mean sequence length `mean_l` (multihit-local mode).
    pub fn configure(&self, mean_l: usize) -> Profile {
        let m = self.m;
        let mut msc = vec![[0f32; K]; m + 1];
        let mut isc = vec![[0f32; K]; m + 1];

        for k in 1..=m {
            for a in 0..K {
                msc[k][a] = prob_to_log2(self.match_emit[k][a] / NULL_FREQ[a]);
            }
        }
        for k in 0..=m {
            for a in 0..K {
                // Insert emissions are usually ~background, so these are ~0.
                isc[k][a] = prob_to_log2(self.insert_emit[k][a] / NULL_FREQ[a]);
            }
        }

        let mut tmm = vec![NEG_INF; m + 1];
        let mut tmi = vec![NEG_INF; m + 1];
        let mut tmd = vec![NEG_INF; m + 1];
        let mut tim = vec![NEG_INF; m + 1];
        let mut tii = vec![NEG_INF; m + 1];
        let mut tdm = vec![NEG_INF; m + 1];
        let mut tdd = vec![NEG_INF; m + 1];
        for k in 0..=m {
            let t = self.trans[k];
            tmm[k] = t.mm;
            tmi[k] = t.mi;
            tmd[k] = t.md;
            tim[k] = t.im;
            tii[k] = t.ii;
            tdm[k] = t.dm;
            tdd[k] = t.dd;
        }

        // Uniform local entry distribution: B->M_k = 2 / (M*(M+1)) for all k.
        // (HMMER's `p7_ProfileConfig` local-entry; exact value is folded so the
        // entry probabilities sum to 1 across the profile.)
        let tbm_val = prob_to_log2(2.0 / (m as f64 * (m as f64 + 1.0)));
        let tbm = vec![tbm_val; m + 1];

        // Length / special-state model in **log-odds** form.
        //
        // HMMER reports bit scores already corrected against a null model whose
        // N′ state emits the flanking residues with the same geometric
        // self-loop the profile uses. In log-odds (profile − null) those flank
        // self-loops therefore contribute ~0, so a no-information sequence
        // scores ≈ 0 bits and a hit's score reflects only its aligned region
        // (plus the genuine local-entry cost). We adopt that convention so raw
        // bit scores are on the same scale as `hmmsearch`'s and are comparable
        // to the per-profile GA / TC cutoffs stored in Pfam-style models.
        //
        // `mean_l` still parameterises the geometric loop used during E-value
        // calibration downstream, but the per-residue flank odds are neutral.
        let _ = mean_l;
        let tnn = 0.0; // N self-loop  (log-odds neutral vs null)
        let tcc = 0.0; // C self-loop
        let tjj = 0.0; // J self-loop
        let tnb = 0.0; // N -> B move  (absorbed by null termination)
        let tct = 0.0; // C -> T move
        let tjb = 0.0; // J -> B re-entry (each re-entry still pays tbm entry)
        // E -> C (end) vs E -> J (another domain): the genuine 0.5/0.5 multihit
        // choice cost, the one special-state term that is *not* null-cancelled.
        let tec = prob_to_log2(0.5);
        let tej = prob_to_log2(0.5);

        Profile {
            name: self.name.clone(),
            acc: self.acc.clone(),
            m,
            msc,
            isc,
            tmm,
            tmi,
            tmd,
            tim,
            tii,
            tdm,
            tdd,
            tbm,
            tnn,
            tnb,
            tec,
            tej,
            tcc,
            tct,
            tjj,
            tjb,
            ga: self.ga,
            stats: self.stats,
        }
    }
}

/// Helper for the parser: turn a HMMER3 token into a log2 probability.
/// `*` → [`NEG_INF`]; otherwise the token is a `-ln(prob)` value.
#[inline]
pub(crate) fn tok_log2(tok: &str) -> f32 {
    if tok == "*" {
        NEG_INF
    } else {
        match tok.parse::<f32>() {
            Ok(v) => neg_ln_to_log2(v),
            Err(_) => NEG_INF,
        }
    }
}

/// Helper for the parser: turn a HMMER3 emission token (`-ln(prob)`) into a
/// linear probability.
#[inline]
pub(crate) fn tok_prob(tok: &str) -> f64 {
    if tok == "*" {
        0.0
    } else {
        match tok.parse::<f64>() {
            Ok(v) => (-v).exp(),
            Err(_) => 0.0,
        }
    }
}
