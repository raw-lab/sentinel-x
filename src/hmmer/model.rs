//! Plan7 profile-HMM data structures and HMMER3 profile configuration.
//!
//! [`Hmm`] is the parsed, probability-space model. [`Hmm::configure`] turns it
//! into a search-ready [`Profile`] exactly the way HMMER3's
//! `p7_ProfileConfig()` does in its default local, multihit mode:
//!
//! * match scores `ln(e_k(x) / f(x))`, with degenerate residues scored by the
//!   background-weighted expected score (`esl_abc_FExpectScVec`);
//! * insert scores hard-wired to 0 (background), with `I_M` impossible;
//! * occupancy-weighted local entry `B→M_k = occ[k] / Σ_i occ[i]·(M−i+1)`;
//! * free local exit `M_k→E = D_k→E = 0`;
//! * node transitions for `k = 1..M−1` only (no node 0, no node M exits).
//!
//! The length-dependent special-state transitions (N/C/J loops, E→C/E→J) are
//! not stored in the profile; they depend on each target's length and are
//! built per target with [`Xsc::new`], mirroring `p7_ReconfigLength()`.
//! All scores are natural-log (nats), as inside HMMER; conversion to bits
//! happens when final scores are reported.

use crate::hmmer::alphabet::{degenerate_members, BG, CODE_STAR, K, NCODE};

/// Negative infinity in score space.
pub const NEG_INF: f32 = f32::NEG_INFINITY;

/// Calibrated statistical parameters from the `STATS LOCAL …` header lines.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Stats {
    /// MSV Gumbel (mu, lambda).
    pub msv: Option<(f64, f64)>,
    /// Viterbi Gumbel (mu, lambda).
    pub viterbi: Option<(f64, f64)>,
    /// Forward exponential tail (tau, lambda).
    pub forward: Option<(f64, f64)>,
}

impl Stats {
    pub fn is_calibrated(&self) -> bool {
        self.msv.is_some() && self.viterbi.is_some() && self.forward.is_some()
    }
}

/// Index of each Plan7 transition in a node's 7-tuple (HMMER `p7H_*` order,
/// which is also the order on disk).
pub const T_MM: usize = 0;
pub const T_MI: usize = 1;
pub const T_MD: usize = 2;
pub const T_IM: usize = 3;
pub const T_II: usize = 4;
pub const T_DM: usize = 5;
pub const T_DD: usize = 6;

/// A parsed, probability-space profile HMM (one record of a `.hmm` file).
#[derive(Clone, Debug)]
pub struct Hmm {
    pub name: String,
    pub acc: Option<String>,
    pub desc: Option<String>,
    /// Number of match nodes.
    pub m: usize,
    /// Match emissions `mat[k][x]`, `k` in `1..=M` (row 0 unused).
    pub mat: Vec<[f32; K]>,
    /// Transitions `t[k][T_*]` for nodes `0..=M`.
    pub t: Vec<[f32; 7]>,
    /// Model composition (the `COMPO` line), if present.
    pub compo: Option<[f32; K]>,
    /// Gathering / trusted / noise cutoffs (sequence, domain) in bits.
    pub ga: Option<(f32, f32)>,
    pub tc: Option<(f32, f32)>,
    pub nc: Option<(f32, f32)>,
    pub stats: Stats,
}

/// A search-ready profile (HMMER's generic `P7_PROFILE`, local mode).
#[derive(Clone, Debug)]
pub struct Profile {
    pub name: String,
    pub acc: Option<String>,
    pub desc: Option<String>,
    pub m: usize,
    /// Match scores laid out residue-major: `msc[x * (M+1) + k]`, `k` in
    /// `1..=M` (column 0 is `-inf`), `x` in `0..NCODE`.
    pub msc: Vec<f32>,
    /// Node transitions (nats), indexed by node `k` in `0..=M`. `t_mm[k]` is
    /// `M_k → M_{k+1}` etc. Entries for `k = 0` and `k = M` are `-inf`.
    pub t_mm: Vec<f32>,
    pub t_mi: Vec<f32>,
    pub t_md: Vec<f32>,
    pub t_im: Vec<f32>,
    pub t_ii: Vec<f32>,
    pub t_dm: Vec<f32>,
    pub t_dd: Vec<f32>,
    /// Local entry `B → M_k`, indexed by `k` in `1..=M`.
    pub t_bm: Vec<f32>,
    /// MSV filter's uniform entry `ln(2 / (M(M+1)))`.
    pub msv_tbm: f32,
    /// Probability-space copies for the scaled Forward/Backward:
    /// match odds `exp(msc)` (same layout as `msc`) and transition probs.
    pub mo: Vec<f32>,
    pub p_mm: Vec<f32>,
    pub p_mi: Vec<f32>,
    pub p_md: Vec<f32>,
    pub p_im: Vec<f32>,
    pub p_ii: Vec<f32>,
    pub p_dm: Vec<f32>,
    pub p_dd: Vec<f32>,
    pub p_bm: Vec<f32>,
    /// Composition used by the bias filter (falls back to the background).
    pub compo: [f32; K],
    pub ga: Option<(f32, f32)>,
    pub tc: Option<(f32, f32)>,
    pub nc: Option<(f32, f32)>,
    pub stats: Stats,
}

impl Profile {
    /// Match-score row for residue code `x` (length `M+1`).
    #[inline(always)]
    pub fn msc_row(&self, x: u8) -> &[f32] {
        let w = self.m + 1;
        let x = x as usize;
        &self.msc[x * w..(x + 1) * w]
    }

    /// Match-odds row for residue code `x` (length `M+1`).
    #[inline(always)]
    pub fn mo_row(&self, x: u8) -> &[f32] {
        let w = self.m + 1;
        let x = x as usize;
        &self.mo[x * w..(x + 1) * w]
    }

    /// Insert score for node `k` and residue code `x` (HMMER hard-wires these
    /// to 0; `I_M` and non-residues are impossible).
    #[inline(always)]
    pub fn isc(&self, k: usize, x: u8) -> f32 {
        if k >= self.m || x == CODE_STAR {
            NEG_INF
        } else {
            0.0
        }
    }
}

impl Hmm {
    /// Configure into a local-mode search [`Profile`] (`p7_ProfileConfig`).
    pub fn configure(&self) -> Profile {
        let m = self.m;
        let w = m + 1;

        // ---- match scores, including degenerate codes ----
        let mut msc = vec![NEG_INF; NCODE * w];
        for k in 1..=m {
            let mut sc = [NEG_INF; NCODE];
            for x in 0..K {
                let p = self.mat[k][x] as f64;
                sc[x] = if p > 0.0 { (p / BG[x] as f64).ln() as f32 } else { NEG_INF };
            }
            // esl_abc_FExpectScVec: background-weighted expected score.
            for code in [20u8, 21, 22, 23] {
                let (mut num, mut den) = (0f32, 0f32);
                for &y in degenerate_members(code) {
                    num += sc[y] * BG[y];
                    den += BG[y];
                }
                sc[code as usize] = if den > 0.0 { num / den } else { 0.0 };
            }
            sc[CODE_STAR as usize] = NEG_INF;
            for x in 0..NCODE {
                msc[x * w + k] = sc[x];
            }
        }

        // ---- transitions: nodes 1..M-1 only ----
        let lnp = |p: f32| if p > 0.0 { (p as f64).ln() as f32 } else { NEG_INF };
        let mut t_mm = vec![NEG_INF; w];
        let mut t_mi = vec![NEG_INF; w];
        let mut t_md = vec![NEG_INF; w];
        let mut t_im = vec![NEG_INF; w];
        let mut t_ii = vec![NEG_INF; w];
        let mut t_dm = vec![NEG_INF; w];
        let mut t_dd = vec![NEG_INF; w];
        for k in 1..m {
            let t = &self.t[k];
            t_mm[k] = lnp(t[T_MM]);
            t_mi[k] = lnp(t[T_MI]);
            t_md[k] = lnp(t[T_MD]);
            t_im[k] = lnp(t[T_IM]);
            t_ii[k] = lnp(t[T_II]);
            t_dm[k] = lnp(t[T_DM]);
            t_dd[k] = lnp(t[T_DD]);
        }
        // D_1 does not exist (no D_0, no B->D_1 in local mode).
        if m >= 1 {
            t_dm[1] = NEG_INF;
            t_dd[1] = NEG_INF;
        }

        // ---- occupancy-weighted local entry ----
        let occ = self.occupancy();
        let mut z = 0f32;
        for (k, &o) in occ.iter().enumerate().skip(1) {
            z += o * (m - k + 1) as f32;
        }
        let mut t_bm = vec![NEG_INF; w];
        for k in 1..=m {
            t_bm[k] = if z > 0.0 && occ[k] > 0.0 { ((occ[k] / z) as f64).ln() as f32 } else { NEG_INF };
        }

        let msv_tbm = (2.0f32 / (m as f32 * (m as f32 + 1.0))).ln();
        let ex = |v: &Vec<f32>| v.iter().map(|&x| x.exp()).collect::<Vec<f32>>();
        let mo = ex(&msc);
        let (p_mm, p_mi, p_md, p_im, p_ii, p_dm, p_dd, p_bm) =
            (ex(&t_mm), ex(&t_mi), ex(&t_md), ex(&t_im), ex(&t_ii), ex(&t_dm), ex(&t_dd), ex(&t_bm));

        Profile {
            name: self.name.clone(),
            acc: self.acc.clone(),
            desc: self.desc.clone(),
            m,
            msc,
            t_mm,
            t_mi,
            t_md,
            t_im,
            t_ii,
            t_dm,
            t_dd,
            t_bm,
            msv_tbm,
            mo,
            p_mm,
            p_mi,
            p_md,
            p_im,
            p_ii,
            p_dm,
            p_dd,
            p_bm,
            compo: self.compo.unwrap_or(BG),
            ga: self.ga,
            tc: self.tc,
            nc: self.nc,
            stats: self.stats,
        }
    }

    /// Match-state occupancy (`p7_hmm_CalculateOccupancy`).
    pub fn occupancy(&self) -> Vec<f32> {
        let m = self.m;
        let mut occ = vec![0f32; m + 1];
        if m == 0 {
            return occ;
        }
        occ[1] = self.t[0][T_MI] + self.t[0][T_MM];
        for k in 2..=m {
            let t = &self.t[k - 1];
            occ[k] = occ[k - 1] * (t[T_MM] + t[T_MI]) + (1.0 - occ[k - 1]) * t[T_DM];
        }
        occ
    }
}

/// Length-dependent special-state transition scores for one target
/// (`p7_ReconfigLength` + multihit/unihit E transitions).
#[derive(Clone, Copy, Debug)]
pub struct Xsc {
    pub n_loop: f32,
    pub n_move: f32,
    pub j_loop: f32,
    pub j_move: f32,
    pub c_loop: f32,
    pub c_move: f32,
    pub e_move: f32,
    pub e_loop: f32,
    /// The same transitions as probabilities, computed directly the way
    /// `p7_oprofile_ReconfigLength` does (not via exp(ln)).
    p: XscP,
}

impl Xsc {
    /// Special transitions for target length `l`; `multihit` selects HMMER's
    /// default multihit mode (E→C = E→J = ½, NJC move `3/(L+3)`) versus
    /// unihit (E→C = 1, no J, move `2/(L+2)`).
    pub fn new(l: usize, multihit: bool) -> Xsc {
        let nj = if multihit { 1.0f32 } else { 0.0 };
        let pmove = (2.0f32 + nj) / (l as f32 + 2.0 + nj);
        let ploop = 1.0f32 - pmove;
        let lloop = (ploop as f64).ln() as f32;
        let lmove = (pmove as f64).ln() as f32;
        let (e_move, e_loop, pe_move, pe_loop) =
            if multihit { (-std::f32::consts::LN_2, -std::f32::consts::LN_2, 0.5f32, 0.5f32) } else { (0.0, NEG_INF, 1.0, 0.0) };
        Xsc {
            n_loop: lloop,
            n_move: lmove,
            j_loop: lloop,
            j_move: lmove,
            c_loop: lloop,
            c_move: lmove,
            e_move,
            e_loop,
            p: XscP {
                n_loop: ploop,
                n_move: pmove,
                j_loop: ploop,
                j_move: pmove,
                c_loop: ploop,
                c_move: pmove,
                e_move: pe_move,
                e_loop: pe_loop,
            },
        }
    }

    /// The same transitions as linear probabilities.
    pub fn probs(&self) -> XscP {
        self.p
    }
}

/// [`Xsc`] in probability space.
#[derive(Clone, Copy, Debug)]
pub struct XscP {
    pub n_loop: f32,
    pub n_move: f32,
    pub j_loop: f32,
    pub j_move: f32,
    pub c_loop: f32,
    pub c_move: f32,
    pub e_move: f32,
    pub e_loop: f32,
}

/// HMMER's null model (`p7_bg_NullOne`) for a target of length `l`, in nats.
pub fn null_score(l: usize) -> f32 {
    let p1 = l as f32 / (l as f32 + 1.0);
    ((l as f64) * (p1 as f64).ln() + (1.0 - p1 as f64).ln()) as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flat_hmm(m: usize) -> Hmm {
        let mut mat = vec![[0f32; K]; m + 1];
        for row in mat.iter_mut().skip(1) {
            *row = BG;
        }
        let tr = [0.9f32, 0.05, 0.05, 0.5, 0.5, 0.5, 0.5];
        Hmm {
            name: "flat".into(),
            acc: None,
            desc: None,
            m,
            mat,
            t: vec![tr; m + 1],
            compo: None,
            ga: None,
            tc: None,
            nc: None,
            stats: Stats::default(),
        }
    }

    #[test]
    fn local_entry_sums_to_one_over_paths() {
        // occ-weighted entry: Σ_k occ[k]/Z · (M-k+1)... the entry distribution
        // times the number of exits reachable is normalised.
        let h = flat_hmm(10);
        let p = h.configure();
        let occ = h.occupancy();
        let z: f32 = (1..=10).map(|k| occ[k] * (10 - k + 1) as f32).sum();
        let s: f32 = (1..=10).map(|k| p.t_bm[k].exp() * (10 - k + 1) as f32).sum();
        assert!((s - 1.0).abs() < 1e-4, "sum {s} (z={z})");
    }

    #[test]
    fn background_match_scores_zero_and_x_zero() {
        let p = flat_hmm(5).configure();
        for x in 0..24u8 {
            for &v in &p.msc_row(x)[1..] {
                assert!(v.abs() < 1e-5, "code {x} score {v}");
            }
        }
        assert!(p.msc_row(CODE_STAR)[1].is_infinite());
    }

    #[test]
    fn null_model_matches_hmmer_formula() {
        let l = 300usize;
        let p1 = 300.0f64 / 301.0;
        let want = 300.0 * p1.ln() + (1.0 - p1).ln();
        assert!((null_score(l) as f64 - want).abs() < 1e-3);
    }
}
