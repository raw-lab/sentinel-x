//! The Plan7 dynamic-programming core: log-odds **Viterbi** (optimal alignment
//! score + domain envelope) and **Forward** (full summed score), in
//! multihit-local mode. These mirror HMMER3's generic `p7_GViterbi` /
//! `p7_GForward` recurrences, kept in bits throughout.

use crate::model::{aa_index, Profile};
use crate::score::{logsum2, NEG_INF};

/// Outcome of scoring one sequence against one profile.
#[derive(Clone, Copy, Debug)]
pub struct DpScore {
    /// Forward bit score (summed over all alignments).
    pub forward_bits: f32,
    /// Viterbi bit score (best single alignment).
    pub viterbi_bits: f32,
    /// Best-alignment envelope on the target (1-based, inclusive); `None` if no
    /// positive-scoring alignment was found.
    pub env: Option<(usize, usize)>,
}

/// Encode a residue sequence to alphabet indices once (non-standard symbols map
/// to `usize::MAX`, which scores as "no emission information").
fn encode(seq: &[u8]) -> Vec<usize> {
    seq.iter()
        .map(|&b| aa_index(b).unwrap_or(usize::MAX))
        .collect()
}

#[inline]
fn msc(p: &Profile, k: usize, a: usize) -> f32 {
    if a == usize::MAX {
        0.0 // unknown residue: log-odds 0 (neutral)
    } else {
        p.msc[k][a]
    }
}
#[inline]
fn isc(p: &Profile, k: usize, a: usize) -> f32 {
    if a == usize::MAX {
        0.0
    } else {
        p.isc[k][a]
    }
}

/// Run both Viterbi and Forward in a single matrix sweep.
/// Score a target sequence against a configured profile, returning the Viterbi
/// (optimal-path) and Forward (full-sum) bit scores plus the best envelope.
///
/// The row-0 special-state initialisations include the E/J/C states for
/// documentation (they are `-inf` before any residue is emitted); the E state
/// happens to be recomputed every row before it is read, so its row-0 value is
/// not consumed — hence the targeted allow below.
#[allow(unused_assignments)]
pub fn score_seq(profile: &Profile, seq: &[u8]) -> DpScore {
    let m = profile.m;
    let xs = encode(seq);
    let l = xs.len();
    if l == 0 || m == 0 {
        return DpScore {
            forward_bits: NEG_INF,
            viterbi_bits: NEG_INF,
            env: None,
        };
    }

    // Two-row DP for the main states (current + previous), full width over k.
    let width = m + 1;
    let mut vm = vec![NEG_INF; width]; // Viterbi M
    let mut vi = vec![NEG_INF; width]; // Viterbi I
    let mut vd = vec![NEG_INF; width]; // Viterbi D
    let mut fm = vec![NEG_INF; width]; // Forward M
    let mut fi = vec![NEG_INF; width];
    let mut fd = vec![NEG_INF; width];
    let mut vm_prev = vec![NEG_INF; width];
    let mut vi_prev = vec![NEG_INF; width];
    let mut fm_prev = vec![NEG_INF; width];
    let mut fi_prev = vec![NEG_INF; width];

    // Special states (scalars per row).
    let (mut vn, mut vb, mut ve, mut vj, mut vc);
    let (mut fn_, mut fb, mut fe, mut fj, mut fc);

    // Row 0.
    vn = 0.0f32;
    fn_ = 0.0f32;
    vb = vn + profile.tnb;
    fb = fn_ + profile.tnb;
    ve = NEG_INF;
    fe = NEG_INF;
    vj = NEG_INF;
    fj = NEG_INF;
    vc = NEG_INF;
    fc = NEG_INF;

    // Track best Viterbi envelope end position.
    let mut best_v = NEG_INF;
    let mut best_end = 0usize;

    for i in 1..=l {
        let a = xs[i - 1];

        std::mem::swap(&mut vm, &mut vm_prev);
        std::mem::swap(&mut vi, &mut vi_prev);
        std::mem::swap(&mut fm, &mut fm_prev);
        std::mem::swap(&mut fi, &mut fi_prev);

        let vb_prev = vb;
        let fb_prev = fb;

        vm[0] = NEG_INF;
        vi[0] = NEG_INF;
        vd[0] = NEG_INF;
        fm[0] = NEG_INF;
        fi[0] = NEG_INF;
        fd[0] = NEG_INF;

        let mut ve_acc = NEG_INF;
        let mut fe_acc = NEG_INF;

        for k in 1..=m {
            // ---- Match ----
            let me = msc(profile, k, a);
            // Viterbi
            let v_from = {
                let mut best = vm_prev[k - 1] + profile.tmm[k - 1];
                best = best.max(vi_prev[k - 1] + profile.tim[k - 1]);
                best = best.max(vd[k - 1] + profile.tdm[k - 1]);
                best = best.max(vb_prev + profile.tbm[k]);
                best
            };
            vm[k] = if v_from == NEG_INF { NEG_INF } else { v_from + me };
            // Forward
            let f_from = logsum2(
                logsum2(
                    fm_prev[k - 1] + profile.tmm[k - 1],
                    fi_prev[k - 1] + profile.tim[k - 1],
                ),
                logsum2(fd[k - 1] + profile.tdm[k - 1], fb_prev + profile.tbm[k]),
            );
            fm[k] = if f_from == NEG_INF { NEG_INF } else { f_from + me };

            // ---- Insert ----
            let ie = isc(profile, k, a);
            let vi_from = (vm_prev[k] + profile.tmi[k]).max(vi_prev[k] + profile.tii[k]);
            vi[k] = if vi_from == NEG_INF {
                NEG_INF
            } else {
                vi_from + ie
            };
            let fi_from = logsum2(fm_prev[k] + profile.tmi[k], fi_prev[k] + profile.tii[k]);
            fi[k] = if fi_from == NEG_INF {
                NEG_INF
            } else {
                fi_from + ie
            };

            // ---- Delete (depends on current row k-1) ----
            vd[k] = (vm[k - 1] + profile.tmd[k - 1]).max(vd[k - 1] + profile.tdd[k - 1]);
            fd[k] = logsum2(
                fm[k - 1] + profile.tmd[k - 1],
                fd[k - 1] + profile.tdd[k - 1],
            );

            // ---- Local exit to E (uniform M_k->E = 0 in local mode) ----
            ve_acc = ve_acc.max(vm[k]);
            fe_acc = logsum2(fe_acc, fm[k]);
        }

        ve = ve_acc;
        fe = fe_acc;

        // J (multihit loop) — Viterbi & Forward.
        vj = (vj + profile.tjj).max(ve + profile.tej);
        fj = logsum2(fj + profile.tjj, fe + profile.tej);

        // C (to terminal).
        vc = (vc + profile.tcc).max(ve + profile.tec);
        fc = logsum2(fc + profile.tcc, fe + profile.tec);

        // N (leading flank).
        vn += profile.tnn;
        fn_ += profile.tnn;

        // B (begin), reachable from N and J.
        vb = (vn + profile.tnb).max(vj + profile.tjb);
        fb = logsum2(fn_ + profile.tnb, fj + profile.tjb);

        if ve > best_v {
            best_v = ve;
            best_end = i;
        }
    }

    let viterbi_bits = vc + profile.tct;
    let forward_bits = fc + profile.tct;

    let env = if best_v > 0.0 {
        Some((best_end.saturating_sub(1).max(1), best_end))
    } else {
        None
    };

    DpScore {
        forward_bits,
        viterbi_bits,
        env,
    }
}
