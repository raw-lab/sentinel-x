//! The per-target search pipeline of `hmmsearch` (`p7_Pipeline`): acceleration
//! filters, Forward/Backward, posterior domain definition, envelope rescoring
//! and null2 bias correction, producing HMMER-scale bit scores.
//!
//! Filter thresholds and domain-definition parameters are HMMER 3.4's
//! defaults (F1 = 0.02, F2 = 1e-3, F3 = 1e-5; rt1 = 0.25, rt2 = 0.10,
//! rt3 = 0.20, omega = 1/256).
//!
//! Multi-domain regions are resolved as HMMER does: 200 stochastic
//! tracebacks (Easel's LCG, seed 42, reseeded per region) are clustered by
//! single linkage and the cluster envelopes rescored; null2 for those regions
//! comes from the sampled traces. Alignment coordinates are taken from each
//! envelope's optimal-accuracy (maximum expected accuracy) alignment, as in
//! HMMER, including its tie-breaking order.

use crate::hmmer::alphabet::{K, NCODE};
use crate::hmmer::dp::{
    backward, backward_matrix, decode, forward, forward_matrix, logsum, msv, null2_from_posterior, optimal_accuracy_coords,
    stochastic_trace, viterbi_score, viterbi_trace, DpBuf, EslRng, Matrix, Posterior, Specials,
};
use crate::hmmer::model::{null_score, Profile, Xsc};

const LN2: f32 = std::f32::consts::LN_2;

/// Pipeline switches (HMMER defaults via [`Default`]).
#[derive(Clone, Copy, Debug)]
pub struct PipelineConfig {
    /// Run the MSV / bias / Viterbi / Forward acceleration filters
    /// (`false` is hmmsearch `--max`).
    pub filters: bool,
    /// Apply the null2 composition-bias correction (`false` = `--nonull2`).
    pub null2: bool,
    pub f1: f64,
    pub f2: f64,
    pub f3: f64,
}

impl Default for PipelineConfig {
    fn default() -> Self {
        PipelineConfig { filters: true, null2: true, f1: 0.02, f2: 1e-3, f3: 1e-5 }
    }
}

/// One domain envelope with HMMER-scale scores.
#[derive(Clone, Debug)]
pub struct RawDomain {
    pub ienv: usize,
    pub jenv: usize,
    pub iali: usize,
    pub jali: usize,
    pub kali: usize,
    pub mali: usize,
    /// Domain bit score (null2-corrected).
    pub bits: f32,
    /// Domain null2 bias (bits).
    pub bias: f32,
    /// ln P-value of the domain score.
    pub lnp: f64,
}

/// A target that survived the pipeline.
#[derive(Clone, Debug)]
pub struct RawHit {
    /// Per-sequence bit score (null2-corrected; HMMER's `score`).
    pub bits: f32,
    /// Bias correction in bits (`pre_score - score`).
    pub bias: f32,
    /// ln P-value of the sequence score (NaN for uncalibrated models).
    pub lnp: f64,
    pub domains: Vec<RawDomain>,
    pub nexpected: f32,
    pub nregions: usize,
}

/// Per-thread scratch space.
#[derive(Default)]
pub struct Workspace {
    buf: DpBuf,
    fsp: Specials,
    bsp: Specials,
    fx: Matrix,
    bx: Matrix,
    pp: Posterior,
}

/// Gumbel survival function (`esl_gumbel_surv`).
pub fn gumbel_surv(x: f64, mu: f64, lambda: f64) -> f64 {
    let y = lambda * (x - mu);
    let ey = -(-y).exp();
    if ey.abs() < 1e-5 {
        -ey
    } else {
        1.0 - ey.exp()
    }
}

/// Exponential-tail ln survival (`esl_exp_logsurv`).
pub fn exp_logsurv(x: f64, tau: f64, lambda: f64) -> f64 {
    if x < tau {
        0.0
    } else {
        -lambda * (x - tau)
    }
}

/// Score of the two-state composition-bias null (`p7_bg_FilterScore`), nats.
fn bias_filter_score(p: &Profile, dsq: &[u8]) -> f32 {
    use crate::hmmer::alphabet::{degenerate_members, BG, CODE_STAR};
    let l = dsq.len();
    let l0 = 400.0f64;
    let l1 = p.m as f64 / 8.0;
    let t = [[l0 / (l0 + 1.0), 1.0 / (l0 + 1.0)], [1.0 / (l1 + 1.0), l1 / (l1 + 1.0)]];
    let pi = [0.999f64, 0.001];
    // Emission odds per code for state 0 (background) and state 1 (model
    // composition), esl_hmm_Configure semantics.
    let mut eo = [[1f64; 2]; NCODE];
    for x in 0..K {
        eo[x][0] = 1.0;
        eo[x][1] = p.compo[x] as f64 / BG[x] as f64;
    }
    for code in [20u8, 21, 22, 23] {
        let mem = degenerate_members(code);
        let (mut e0, mut e1, mut z) = (0f64, 0f64, 0f64);
        for &y in mem {
            e0 += BG[y] as f64;
            e1 += p.compo[y] as f64;
            z += BG[y] as f64;
        }
        eo[code as usize] = [e0 / z, e1 / z];
    }
    eo[CODE_STAR as usize] = [1.0, 1.0];

    let mut logsc = 0f64;
    let x1 = dsq[0] as usize;
    let mut dp = [eo[x1][0] * pi[0], eo[x1][1] * pi[1]];
    let mx = dp[0].max(dp[1]).max(f64::MIN_POSITIVE);
    dp = [dp[0] / mx, dp[1] / mx];
    logsc += mx.ln();
    for &x in &dsq[1..] {
        let x = x as usize;
        let n0 = (dp[0] * t[0][0] + dp[1] * t[1][0]) * eo[x][0];
        let n1 = (dp[0] * t[0][1] + dp[1] * t[1][1]) * eo[x][1];
        let mx = n0.max(n1).max(f64::MIN_POSITIVE);
        dp = [n0 / mx, n1 / mx];
        logsc += mx.ln();
    }
    logsc += (dp[0] + dp[1]).ln();
    let p1 = l as f64 / (l as f64 + 1.0);
    (logsc + l as f64 * p1.ln() + (1.0 - p1).ln()) as f32
}

/// Run the full hmmsearch pipeline for one profile/target pair. Returns
/// `None` when the target is rejected by a filter or has no domains.
pub fn score_target(p: &Profile, dsq: &[u8], cfg: &PipelineConfig, ws: &mut Workspace) -> Option<RawHit> {
    let l = dsq.len();
    if l == 0 || p.m == 0 {
        return None;
    }
    let nullsc = null_score(l);
    let xs = Xsc::new(l, true);
    let calibrated = p.stats.is_calibrated();

    let fwdsc;
    if cfg.filters && calibrated {
        let (mmu, mlam) = p.stats.msv.unwrap();
        let usc = msv(p, dsq, &mut ws.buf);
        if gumbel_surv(((usc - nullsc) / LN2) as f64, mmu, mlam) > cfg.f1 {
            return None;
        }
        let filtersc = bias_filter_score(p, dsq);
        if gumbel_surv(((usc - filtersc) / LN2) as f64, mmu, mlam) > cfg.f1 {
            return None;
        }
        let (vmu, vlam) = p.stats.viterbi.unwrap();
        let vsc = viterbi_score(p, dsq, &xs, &mut ws.buf);
        if gumbel_surv(((vsc - filtersc) / LN2) as f64, vmu, vlam) > cfg.f2 {
            return None;
        }
        let (ftau, flam) = p.stats.forward.unwrap();
        fwdsc = forward(p, dsq, &xs, &mut ws.buf, Some(&mut ws.fsp));
        if exp_logsurv(((fwdsc - filtersc) / LN2) as f64, ftau, flam).exp() > cfg.f3 {
            return None;
        }
    } else {
        fwdsc = forward(p, dsq, &xs, &mut ws.buf, Some(&mut ws.fsp));
    }
    if !fwdsc.is_finite() {
        return None;
    }
    let bcksc = backward(p, dsq, &xs, &mut ws.buf, &mut ws.bsp);
    let overall = if bcksc.is_finite() { bcksc } else { fwdsc };

    // ---- posterior domain decoding (p7_GDomainDecoding) ----
    let f = &ws.fsp;
    let b = &ws.bsp;
    let mut btot = vec![0f32; l + 1];
    let mut etot = vec![0f32; l + 1];
    let mut mocc = vec![0f32; l + 1];
    for i in 1..=l {
        btot[i] = btot[i - 1] + (f.b[i - 1] + b.b[i - 1] - overall).exp();
        etot[i] = etot[i - 1] + (f.e[i] + b.e[i] - overall).exp();
        let njc = (f.n[i - 1] + b.n[i] + xs.n_loop - overall).exp()
            + (f.j[i - 1] + b.j[i] + xs.j_loop - overall).exp()
            + (f.c[i - 1] + b.c[i] + xs.c_loop - overall).exp();
        mocc[i] = 1.0 - njc;
    }
    let nexpected = btot[l];

    // ---- region identification (p7_domaindef_ByPosteriorHeuristics) ----
    let (rt1, rt2, rt3) = (0.25f32, 0.10f32, 0.20f32);
    let mut n2sc = vec![0f32; l + 1];
    let mut envs: Vec<Envelope> = Vec::new();
    let mut nregions = 0usize;
    let mut i: isize = -1;
    let mut triggered = false;
    for j in 1..=l {
        if !triggered {
            if mocc[j] - (btot[j] - btot[j - 1]) < rt2 || i == -1 {
                i = j as isize;
            }
            if mocc[j] >= rt1 {
                triggered = true;
            }
        } else if mocc[j] - (etot[j] - etot[j - 1]) < rt2 {
            let ri = i.max(1) as usize;
            nregions += 1;
            if is_multidomain(&btot, &etot, ri, j, rt3) {
                // HMMER: resolve by clustering 200 stochastic tracebacks
                // (seed 42, reseeded per region); null2 comes from the traces.
                let envs_r = region_trace_ensemble(p, dsq, ri, j, l, ws, &mut n2sc);
                for (ei, ej) in envs_r {
                    if let Some(d) = rescore(p, dsq, ei, ej, l, cfg, ws, &mut n2sc, true) {
                        envs.push(d);
                    }
                }
            } else if let Some(d) = rescore(p, dsq, ri, j, l, cfg, ws, &mut n2sc, false) {
                envs.push(d);
            }
            i = -1;
            triggered = false;
        }
    }
    if nregions == 0 || envs.is_empty() {
        return None;
    }

    // ---- sequence and domain scores (p7_Pipeline) ----
    let omega_ln = (1.0f32 / 256.0).ln();
    let seqbias = if cfg.null2 {
        let s: f32 = n2sc[1..].iter().sum();
        logsum(0.0, omega_ln + s)
    } else {
        0.0
    };
    let mut pre_score = (fwdsc - nullsc) / LN2;
    let mut seq_score = (fwdsc - (nullsc + seqbias)) / LN2;

    let ln_multi_loop = (l as f32 / (l as f32 + 3.0)).ln();
    let mut sum = 0f32;
    let mut sum_bias = 0f32;
    let mut ld_tot = 0usize;
    for &(ie, je, envsc, domcorr, _) in &envs {
        let keep = if cfg.null2 { envsc - domcorr > 0.0 } else { envsc > 0.0 };
        if keep {
            sum += envsc;
            ld_tot += je - ie + 1;
            sum_bias += domcorr;
        }
    }
    let sum_seqbias = if cfg.null2 { logsum(0.0, omega_ln + sum_bias) } else { 0.0 };
    sum += (l - ld_tot.min(l)) as f32 * ln_multi_loop;
    let pre2 = (sum - nullsc) / LN2;
    let sum_score = (sum - (nullsc + sum_seqbias)) / LN2;
    if ld_tot > 0 && sum_score > seq_score {
        seq_score = sum_score;
        pre_score = pre2;
    }

    let (ftau, flam) = p.stats.forward.unwrap_or((f64::NAN, f64::NAN));
    let lnp = if calibrated || p.stats.forward.is_some() { exp_logsurv(seq_score as f64, ftau, flam) } else { f64::NAN };

    let domains = envs
        .into_iter()
        .map(|(ie, je, envsc, domcorr, (ia, ja, ka, ma))| {
            let ld = je - ie + 1;
            let dombias = if cfg.null2 { logsum(0.0, omega_ln + domcorr) } else { 0.0 };
            let bits = (envsc + (l - ld) as f32 * ln_multi_loop - (nullsc + dombias)) / LN2;
            let lnp = if p.stats.forward.is_some() { exp_logsurv(bits as f64, ftau, flam) } else { f64::NAN };
            RawDomain { ienv: ie, jenv: je, iali: ia, jali: ja, kali: ka, mali: ma, bits, bias: dombias / LN2, lnp }
        })
        .collect();

    Some(RawHit { bits: seq_score, bias: pre_score - seq_score, lnp, domains, nexpected, nregions })
}

/// Sampled domain coordinates for clustering (HMMER `p7_spcoord_s`).
#[derive(Clone, Copy, Debug)]
struct SpCoord {
    idx: usize,
    i: usize,
    j: usize,
    k: usize,
    m: usize,
}

/// `link_spsamples` with HMMER's defaults (min_overlap 0.8 of the smaller,
/// max_diagdiff 4).
fn sp_link(h1: &SpCoord, h2: &SpCoord) -> bool {
    let nov = h1.j.min(h2.j) as i64 - h1.i.max(h2.i) as i64 + 1;
    let n = ((h1.j - h1.i + 1).min(h2.j - h2.i + 1)) as f32;
    if (nov as f32) / n < 0.8 {
        return false;
    }
    let nov = h1.m.min(h2.m) as i64 - h1.k.max(h2.k) as i64;
    let n = ((h1.m - h1.k + 1).min(h2.m - h2.k + 1)) as f32;
    if (nov as f32) / n < 0.8 {
        return false;
    }
    let d1 = h1.i as i64 - h1.k as i64;
    let d2 = h2.i as i64 - h2.k as i64;
    if (d1 - d2).abs() <= 4 {
        return true;
    }
    let d1 = h1.j as i64 - h1.m as i64;
    let d2 = h2.j as i64 - h2.m as i64;
    (d1 - d2).abs() <= 4
}

/// `esl_cluster_SingleLinkage` (same traversal order, so cluster ids match).
fn single_linkage(v: &[SpCoord]) -> (Vec<usize>, usize) {
    let n = v.len();
    let mut a: Vec<usize> = (0..n).map(|x| n - x - 1).collect();
    let mut b: Vec<usize> = Vec::with_capacity(n);
    let mut c = vec![0usize; n];
    let mut nc = 0usize;
    while let Some(v0) = a.pop() {
        b.push(v0);
        while let Some(x) = b.pop() {
            c[x] = nc;
            let mut ii = a.len();
            while ii > 0 {
                ii -= 1;
                if sp_link(&v[x], &v[a[ii]]) {
                    let w = a[ii];
                    let last = a.len() - 1;
                    a[ii] = a[last];
                    a.pop();
                    b.push(w);
                }
            }
        }
        nc += 1;
    }
    (c, nc)
}

/// HMMER `region_trace_ensemble` + `p7_spensemble_Cluster`: sample 200
/// traces of region `ireg..=jreg`, set position-specific null2 scores, and
/// return the cluster envelopes (sorted by start, dominated ones removed).
fn region_trace_ensemble(
    p: &Profile,
    dsq: &[u8],
    ireg: usize,
    jreg: usize,
    l: usize,
    ws: &mut Workspace,
    n2sc: &mut [f32],
) -> Vec<(usize, usize)> {
    const NSAMPLES: usize = 200;
    let sub = &dsq[ireg - 1..jreg];
    let lr = sub.len();
    let xs = Xsc::new(l, true);
    forward_matrix(p, sub, &xs, &mut ws.fx);
    let mut rng = EslRng::new(42);
    for v in &mut n2sc[ireg..=jreg] {
        *v = 0.0;
    }
    let mut sp: Vec<SpCoord> = Vec::new();
    let w = p.m + 1;
    let mut null2 = [0f32; NCODE];
    for t in 0..NSAMPLES {
        let doms = stochastic_trace(p, &ws.fx, &xs, &mut rng);
        let mut pos = 1usize;
        for d in &doms {
            sp.push(SpCoord { idx: t, i: d.sqfrom + ireg - 1, j: d.sqto + ireg - 1, k: d.hmmfrom, m: d.hmmto });
            // p7_Null2_ByTrace over this domain's emissions.
            let ld = d.nodes.len().max(1) as f32;
            for (x, slot) in null2.iter_mut().enumerate().take(K) {
                let row = &p.mo[x * w..(x + 1) * w];
                let mut v = 0f32;
                for &kk in &d.nodes {
                    v += row[kk];
                }
                *slot = v / ld;
            }
            for code in [20u8, 21, 22, 23] {
                let mem = crate::hmmer::alphabet::degenerate_members(code);
                let s: f32 = mem.iter().map(|&y| null2[y]).sum();
                null2[code as usize] = s / mem.len() as f32;
            }
            null2[crate::hmmer::alphabet::CODE_STAR as usize] = 1.0;
            while pos <= d.sqfrom {
                n2sc[ireg + pos - 1] += 1.0;
                pos += 1;
            }
            while pos <= d.sqto {
                n2sc[ireg + pos - 1] += null2[dsq[ireg + pos - 2] as usize];
                pos += 1;
            }
        }
        while pos <= lr {
            n2sc[ireg + pos - 1] += 1.0;
            pos += 1;
        }
    }
    for v in &mut n2sc[ireg..=jreg] {
        *v = (*v / NSAMPLES as f32).ln();
    }

    // ---- cluster ----
    let (assign, nc) = single_linkage(&sp);
    struct Sig {
        i: usize,
        j: usize,
        prob: f32,
    }
    let mut sigc: Vec<Sig> = Vec::new();
    for c in 0..nc {
        let mut ninc = 0usize;
        let mut last: Option<usize> = None;
        for (h, s) in sp.iter().enumerate() {
            if assign[h] == c {
                if last != Some(s.idx) {
                    ninc += 1;
                }
                last = Some(s.idx);
            }
        }
        if (ninc as f32) / (NSAMPLES as f32) < 0.25 {
            continue;
        }
        let members: Vec<&SpCoord> = sp.iter().enumerate().filter(|(h, _)| assign[*h] == c).map(|(_, s)| s).collect();
        let (imin, imax) = minmax(members.iter().map(|s| s.i));
        let (jmin, jmax) = minmax(members.iter().map(|s| s.j));
        let (kmin, kmax) = minmax(members.iter().map(|s| s.k));
        let (mmin, mmax) = minmax(members.iter().map(|s| s.m));
        let thr = ((ninc as f32) * 0.02).ceil() as usize;
        let best_i = leftmost(members.iter().map(|s| s.i), imin, imax, thr);
        let best_k = leftmost(members.iter().map(|s| s.k), kmin, kmax, thr);
        let best_j = rightmost(members.iter().map(|s| s.j), jmin, jmax, thr);
        let best_m = rightmost(members.iter().map(|s| s.m), mmin, mmax, thr);
        if best_i > best_j || best_k > best_m {
            continue;
        }
        sigc.push(Sig { i: best_i, j: best_j, prob: ninc as f32 / NSAMPLES as f32 });
    }
    sigc.sort_by_key(|s| s.i);
    // Remove dominated envelopes (>= 80% overlap of the smaller).
    let mut dominated = vec![false; sigc.len()];
    for d in 0..sigc.len() {
        for d2 in d + 1..sigc.len() {
            let nov = sigc[d].j.min(sigc[d2].j) as i64 - sigc[d].i.max(sigc[d2].i) as i64 + 1;
            if nov == 0 {
                break;
            }
            let n = (sigc[d].j - sigc[d].i + 1).min(sigc[d2].j - sigc[d2].i + 1) as f32;
            if nov as f32 / n >= 0.8 {
                if sigc[d].prob > sigc[d2].prob {
                    dominated[d2] = true;
                } else {
                    dominated[d] = true;
                }
            }
        }
    }
    sigc.iter().zip(dominated).filter(|(_, dom)| !dom).map(|(s, _)| (s.i, s.j)).collect()
}

fn minmax<I: Iterator<Item = usize>>(it: I) -> (usize, usize) {
    let mut lo = usize::MAX;
    let mut hi = 0usize;
    for v in it {
        lo = lo.min(v);
        hi = hi.max(v);
    }
    (lo, hi)
}

fn endpoint_counts<I: Iterator<Item = usize>>(it: I, lo: usize, hi: usize) -> Vec<usize> {
    let mut epc = vec![0usize; hi - lo + 1];
    for v in it {
        epc[v - lo] += 1;
    }
    epc
}

fn argmax(v: &[usize]) -> usize {
    let mut best = 0;
    for (i, &x) in v.iter().enumerate() {
        if x > v[best] {
            best = i;
        }
    }
    best
}

fn leftmost<I: Iterator<Item = usize>>(it: I, lo: usize, hi: usize, thr: usize) -> usize {
    let epc = endpoint_counts(it, lo, hi);
    for (o, &c) in epc.iter().enumerate() {
        if c >= thr {
            return lo + o;
        }
    }
    lo + argmax(&epc)
}

fn rightmost<I: Iterator<Item = usize>>(it: I, lo: usize, hi: usize, thr: usize) -> usize {
    let epc = endpoint_counts(it, lo, hi);
    for o in (0..epc.len()).rev() {
        if epc[o] >= thr {
            return lo + o;
        }
    }
    lo + argmax(&epc)
}

fn is_multidomain(btot: &[f32], etot: &[f32], i: usize, j: usize, rt3: f32) -> bool {
    let mut max = -1.0f32;
    for jj in i..=j {
        let v = (etot[jj] - etot[i - 1]).min(btot[j] - btot[jj - 1]);
        if v > max {
            max = v;
        }
    }
    max >= rt3
}

/// A rescored envelope: `(ienv, jenv, envsc, domcorrection, (iali, jali, kali, mali))`.
type Envelope = (usize, usize, f32, f32, (usize, usize, usize, usize));

/// Rescore an isolated envelope `i..=j` (1-based, inclusive) in unihit mode
/// with the full-length config (`rescore_isolated_domain`). Returns
/// `(ienv, jenv, envsc, domcorrection, (iali, jali, kali, mali))`.
#[allow(clippy::too_many_arguments)]
fn rescore(
    p: &Profile,
    dsq: &[u8],
    i: usize,
    j: usize,
    l: usize,
    cfg: &PipelineConfig,
    ws: &mut Workspace,
    n2sc: &mut [f32],
    null2_done: bool,
) -> Option<Envelope> {
    if i == 0 || j < i || j > l {
        return None;
    }
    let sub = &dsq[i - 1..j];
    let xs = Xsc::new(l, false);
    let envsc = forward_matrix(p, sub, &xs, &mut ws.fx);
    if !envsc.is_finite() {
        return None;
    }
    let bsc = backward_matrix(p, sub, &xs, &mut ws.bx);
    if !bsc.is_finite() {
        return None;
    }
    decode(p, &xs, &ws.fx, &ws.bx, &mut ws.pp);
    let mut domcorr = 0f32;
    if cfg.null2 {
        if !null2_done {
            let n2 = null2_from_posterior(p, &ws.pp);
            for pos in i..=j {
                let v = n2[dsq[pos - 1] as usize];
                n2sc[pos] = if v.is_finite() { v } else { 0.0 };
            }
        }
        for v in &n2sc[i..=j] {
            domcorr += *v;
        }
    }
    let ali = match optimal_accuracy_coords(p, &xs, &ws.pp) {
        Some(a) => (i - 1 + a.iali, i - 1 + a.jali, a.kali, a.mali),
        None => match viterbi_trace(p, sub, &xs).first() {
            Some(a) => (i - 1 + a.iali, i - 1 + a.jali, a.kali, a.mali),
            None => (i, j, 1, p.m),
        },
    };
    Some((i, j, envsc, domcorr, ali))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hmmer::alphabet::encode;
    use crate::hmmer::parse::parse_hmm_str;

    #[test]
    fn gumbel_and_exp_tails() {
        assert!((gumbel_surv(0.0, 0.0, 1.0) - (1.0 - (-1.0f64).exp())).abs() < 1e-12);
        assert_eq!(exp_logsurv(-5.0, 0.0, 0.7), 0.0);
        assert!((exp_logsurv(10.0, 0.0, 0.7) + 7.0).abs() < 1e-12);
    }

    #[test]
    fn finds_both_copies_of_a_motif() {
        let p = parse_hmm_str(include_str!("testdata/tiny.hmm")).unwrap()[0].configure();
        let seq = encode(b"GGGGGMKAGGGGGGGGGGGGGGGGGGGGMKAGGGGG");
        let cfg = PipelineConfig { filters: false, ..Default::default() };
        let mut ws = Workspace::default();
        let h = score_target(&p, &seq, &cfg, &mut ws).expect("hit");
        assert_eq!(h.domains.len(), 2, "{:?}", h.domains);
        assert!(h.bits > 0.0);
        assert_eq!((h.domains[0].iali, h.domains[0].jali), (6, 8));
    }
}
