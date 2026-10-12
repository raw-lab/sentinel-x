//! Plan7 dynamic programming, following HMMER3's generic implementations
//! (`p7_GMSV`, `p7_GViterbi`, `p7_GForward`, `p7_GBackward`, `p7_GDecoding`,
//! `p7_GNull2_ByExpectation`). All scores are natural-log (nats).
//!
//! * [`msv`] — ungapped multi-segment Viterbi, the first acceleration filter.
//! * [`viterbi_score`] — optimal-path score (second filter).
//! * [`forward`] / [`backward`] — full-sum scores over the whole target in
//!   linear memory, optionally keeping the special-state rows needed for
//!   posterior domain decoding.
//! * [`Matrix`], [`forward_matrix`], [`backward_matrix`], [`null2_by_expectation`]
//!   — full-matrix versions used to rescore individual domain envelopes.
//! * [`viterbi_trace`] — optimal alignment with traceback, giving domain
//!   coordinates on the model and the target.

use crate::hmmer::alphabet::{degenerate_members, CODE_STAR, K, NCODE};
use crate::hmmer::model::{Profile, Xsc, NEG_INF};

const LN2: f32 = std::f32::consts::LN_2;

#[inline(always)]
fn fmax(a: f32, b: f32) -> f32 {
    if a > b {
        a
    } else {
        b
    }
}

/// `ln(e^a + e^b)` without overflow; exact to f32 precision.
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)] // NaN-aware on purpose
pub fn logsum(a: f32, b: f32) -> f32 {
    let (hi, lo) = if a > b { (a, b) } else { (b, a) };
    let d = hi - lo;
    if !(d < 30.0) {
        // lo is -inf (d = inf), both -inf (d = NaN), or lo negligible.
        return hi;
    }
    hi + (-d).exp().ln_1p()
}

/// Reusable scratch rows for the linear-memory DP routines.
#[derive(Default)]
pub struct DpBuf {
    a: Vec<f32>,
    b: Vec<f32>,
    c: Vec<f32>,
    d: Vec<f32>,
    e: Vec<f32>,
    f: Vec<f32>,
}

impl DpBuf {
    fn rows(&mut self, w: usize) -> [&mut Vec<f32>; 6] {
        for v in [&mut self.a, &mut self.b, &mut self.c, &mut self.d, &mut self.e, &mut self.f] {
            v.clear();
            v.resize(w, NEG_INF);
        }
        [&mut self.a, &mut self.b, &mut self.c, &mut self.d, &mut self.e, &mut self.f]
    }
}

// ---------------------------------------------------------------------------
// MSV filter
// ---------------------------------------------------------------------------

/// MSV score (nats) of `dsq` against `p` (HMMER `p7_GMSV`, nu = 2): ungapped
/// local segments, uniform entry, multihit, length-configured to `dsq.len()`.
pub fn msv(p: &Profile, dsq: &[u8], buf: &mut DpBuf) -> f32 {
    #[cfg(target_arch = "x86_64")]
    {
        if std::is_x86_feature_detected!("avx2") {
            // SAFETY: the CPU supports AVX2 (checked above); the function only
            // enables the target feature for auto-vectorisation.
            return unsafe { msv_avx2(p, dsq, buf) };
        }
    }
    msv_generic(p, dsq, buf)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn msv_avx2(p: &Profile, dsq: &[u8], buf: &mut DpBuf) -> f32 {
    msv_kernel(p, dsq, buf)
}

fn msv_generic(p: &Profile, dsq: &[u8], buf: &mut DpBuf) -> f32 {
    msv_kernel(p, dsq, buf)
}

#[inline(always)]
fn msv_kernel(p: &Profile, dsq: &[u8], buf: &mut DpBuf) -> f32 {
    let m = p.m;
    let l = dsq.len();
    if l == 0 || m == 0 {
        return NEG_INF;
    }
    let lf = l as f32;
    let tloop = (lf / (lf + 3.0)).ln();
    let tmove = (3.0f32 / (lf + 3.0)).ln();
    let tbmk = p.msv_tbm;
    let (tec, tej) = (-LN2, -LN2);

    let [prev, cur, ..] = buf.rows(m + 1);
    let mut xn = 0.0f32;
    let mut xb = tmove;
    let mut xj = NEG_INF;
    let mut xc = NEG_INF;

    for &x in dsq {
        let rs = &p.msc_row(x)[1..=m];
        let be = xb + tbmk;
        let mut e8 = [NEG_INF; 8];
        {
            let pv = &prev[0..m];
            let cu = &mut cur[1..=m];
            let mut cch = cu.chunks_exact_mut(8);
            let mut rch = rs.chunks_exact(8);
            let mut pch = pv.chunks_exact(8);
            for ((c, r), pp) in (&mut cch).zip(&mut rch).zip(&mut pch) {
                for j in 0..8 {
                    let v = r[j] + fmax(pp[j], be);
                    c[j] = v;
                    e8[j] = fmax(e8[j], v);
                }
            }
            for ((c, r), pp) in cch.into_remainder().iter_mut().zip(rch.remainder()).zip(pch.remainder()) {
                let v = *r + fmax(*pp, be);
                *c = v;
                e8[0] = fmax(e8[0], v);
            }
        }
        cur[0] = NEG_INF;
        let mut xe = NEG_INF;
        for v in e8 {
            xe = fmax(xe, v);
        }
        xj = fmax(xj + tloop, xe + tej);
        xc = fmax(xc + tloop, xe + tec);
        xn += tloop;
        xb = fmax(xn + tmove, xj + tmove);
        std::mem::swap(prev, cur);
    }
    xc + tmove
}

// ---------------------------------------------------------------------------
// Viterbi score (filter)
// ---------------------------------------------------------------------------

/// Optimal single-path score (nats) under the length-configured profile.
pub fn viterbi_score(p: &Profile, dsq: &[u8], xs: &Xsc, buf: &mut DpBuf) -> f32 {
    let m = p.m;
    if dsq.is_empty() || m == 0 {
        return NEG_INF;
    }
    let w = m + 1;
    let [mp, ip, dp, mc, ic, dc] = buf.rows(w);
    let mut xn = 0.0f32;
    let mut xb = xs.n_move;
    let mut xj = NEG_INF;
    let mut xc = NEG_INF;
    for &x in dsq {
        let rs = p.msc_row(x);
        let ins_ok = x != CODE_STAR;
        mc[0] = NEG_INF;
        ic[0] = NEG_INF;
        dc[0] = NEG_INF;
        let mut xe = NEG_INF;
        for k in 1..=m {
            let sc =
                fmax(fmax(mp[k - 1] + p.t_mm[k - 1], ip[k - 1] + p.t_im[k - 1]), fmax(xb + p.t_bm[k], dp[k - 1] + p.t_dm[k - 1]));
            let mv = sc + rs[k];
            mc[k] = mv;
            ic[k] = if k < m && ins_ok { fmax(mp[k] + p.t_mi[k], ip[k] + p.t_ii[k]) } else { NEG_INF };
            let dv = fmax(mc[k - 1] + p.t_md[k - 1], dc[k - 1] + p.t_dd[k - 1]);
            dc[k] = dv;
            xe = fmax(xe, fmax(mv, dv));
        }
        xj = fmax(xj + xs.j_loop, xe + xs.e_loop);
        xc = fmax(xc + xs.c_loop, xe + xs.e_move);
        xn += xs.n_loop;
        xb = fmax(xn + xs.n_move, xj + xs.j_move);
        std::mem::swap(mp, mc);
        std::mem::swap(ip, ic);
        std::mem::swap(dp, dc);
    }
    xc + xs.c_move
}

// ---------------------------------------------------------------------------
// Forward / Backward (probability space with sparse rescaling, as in HMMER's
// optimised implementation), linear memory, with special-state rows
// ---------------------------------------------------------------------------

/// Rescaling threshold (HMMER rescales when xE exceeds 1e4).
const RESCALE: f32 = 1.0e4;

#[inline(always)]
fn lnf(x: f32) -> f32 {
    if x > 0.0 {
        x.ln()
    } else {
        NEG_INF
    }
}

/// Special-state rows `N, B, E, J, C` for rows `0..=L`, as true natural-log
/// values (scaling already folded in).
#[derive(Default, Clone)]
pub struct Specials {
    pub n: Vec<f32>,
    pub b: Vec<f32>,
    pub e: Vec<f32>,
    pub j: Vec<f32>,
    pub c: Vec<f32>,
}

impl Specials {
    fn reset(&mut self, l: usize) {
        for v in [&mut self.n, &mut self.b, &mut self.e, &mut self.j, &mut self.c] {
            v.clear();
            v.resize(l + 1, NEG_INF);
        }
    }
}

/// Forward score (nats). If `sp` is given, it receives the special rows.
pub fn forward(p: &Profile, dsq: &[u8], xs: &Xsc, buf: &mut DpBuf, sp: Option<&mut Specials>) -> f32 {
    let m = p.m;
    let l = dsq.len();
    if l == 0 || m == 0 {
        return NEG_INF;
    }
    let mut sp = sp;
    if let Some(s) = sp.as_mut() {
        s.reset(l);
    }
    let xp = xs.probs();
    let w = m + 1;
    let [mp, ip, dp, mc, ic, dc] = buf.rows(w);
    for v in [&mut *mp, &mut *ip, &mut *dp, &mut *mc, &mut *ic, &mut *dc] {
        v.iter_mut().for_each(|x| *x = 0.0);
    }
    let mut xn = 1.0f32;
    let mut xb = xp.n_move;
    let mut xj = 0.0f32;
    let mut xc = 0.0f32;
    let mut lnscale = 0.0f64;
    if let Some(s) = sp.as_mut() {
        s.n[0] = 0.0;
        s.b[0] = lnf(xb);
    }
    for (i0, &x) in dsq.iter().enumerate() {
        let eo = p.mo_row(x);
        let ins = if x != CODE_STAR { 1.0f32 } else { 0.0 };
        let mut xe = 0.0f32;
        for k in 1..=m {
            let mv = (xb * p.p_bm[k] + mp[k - 1] * p.p_mm[k - 1] + ip[k - 1] * p.p_im[k - 1] + dp[k - 1] * p.p_dm[k - 1]) * eo[k];
            mc[k] = mv;
            ic[k] = if k < m { (mp[k] * p.p_mi[k] + ip[k] * p.p_ii[k]) * ins } else { 0.0 };
            let dv = mc[k - 1] * p.p_md[k - 1] + dc[k - 1] * p.p_dd[k - 1];
            dc[k] = dv;
            xe += mv + dv;
        }
        xj = xj * xp.j_loop + xe * xp.e_loop;
        xc = xc * xp.c_loop + xe * xp.e_move;
        xn *= xp.n_loop;
        xb = xn * xp.n_move + xj * xp.j_move;
        if xe > RESCALE {
            let inv = 1.0 / xe;
            for k in 1..=m {
                mc[k] *= inv;
                ic[k] *= inv;
                dc[k] *= inv;
            }
            xn *= inv;
            xb *= inv;
            xj *= inv;
            xc *= inv;
            lnscale += (xe as f64).ln();
            xe = 1.0;
        }
        if let Some(s) = sp.as_mut() {
            let i = i0 + 1;
            let ls = lnscale as f32;
            s.n[i] = lnf(xn) + ls;
            s.b[i] = lnf(xb) + ls;
            s.e[i] = lnf(xe) + ls;
            s.j[i] = lnf(xj) + ls;
            s.c[i] = lnf(xc) + ls;
        }
        std::mem::swap(mp, mc);
        std::mem::swap(ip, ic);
        std::mem::swap(dp, dc);
    }
    (lnf(xc * xp.c_move) as f64 + lnscale) as f32
}

/// Backward score (nats); `sp` receives the special rows (`p7_GBackward`).
pub fn backward(p: &Profile, dsq: &[u8], xs: &Xsc, buf: &mut DpBuf, sp: &mut Specials) -> f32 {
    let m = p.m;
    let l = dsq.len();
    sp.reset(l);
    if l == 0 || m == 0 {
        return NEG_INF;
    }
    let xp = xs.probs();
    let w = m + 1;
    let [mn, inx, dn, mc, ic, dc] = buf.rows(w);
    for v in [&mut *mn, &mut *inx, &mut *dn, &mut *mc, &mut *ic, &mut *dc] {
        v.iter_mut().for_each(|x| *x = 0.0);
    }
    // Row L.
    let mut xc = xp.c_move;
    let mut xe = xc * xp.e_move;
    let mut xj = 0.0f32;
    let mut xn = 0.0f32;
    let mut xb;
    mn[m] = xe;
    dn[m] = xe;
    for k in (1..m).rev() {
        mn[k] = xe + dn[k + 1] * p.p_md[k];
        dn[k] = xe + dn[k + 1] * p.p_dd[k];
    }
    let mut lnscale = 0.0f64; // cumulative scale of rows >= i
    sp.c[l] = lnf(xc);
    sp.e[l] = lnf(xe);

    for i in (1..l).rev() {
        let xnext = dsq[i]; // residue i+1
        let eo = p.mo_row(xnext);
        let ins = if xnext != CODE_STAR { 1.0f32 } else { 0.0 };
        let mut b = 0.0f32;
        for k in 1..=m {
            b += mn[k] * p.p_bm[k] * eo[k];
        }
        xb = b;
        xj = xj * xp.j_loop + xb * xp.j_move;
        xc *= xp.c_loop;
        xe = xj * xp.e_loop + xc * xp.e_move;
        xn = xn * xp.n_loop + xb * xp.n_move;
        mc[m] = xe;
        dc[m] = xe;
        ic[m] = 0.0;
        let mut rowmax = xb.max(xe).max(xn).max(xj).max(xc);
        for k in (1..m).rev() {
            let mnext = mn[k + 1] * eo[k + 1];
            let iv = inx[k] * ins;
            let mv = mnext * p.p_mm[k] + iv * p.p_mi[k] + xe + dc[k + 1] * p.p_md[k];
            mc[k] = mv;
            ic[k] = mnext * p.p_im[k] + iv * p.p_ii[k];
            dc[k] = mnext * p.p_dm[k] + dc[k + 1] * p.p_dd[k] + xe;
            if mv > rowmax {
                rowmax = mv;
            }
        }
        if rowmax > RESCALE {
            let inv = 1.0 / rowmax;
            for k in 1..=m {
                mc[k] *= inv;
                ic[k] *= inv;
                dc[k] *= inv;
            }
            xb *= inv;
            xe *= inv;
            xn *= inv;
            xj *= inv;
            xc *= inv;
            lnscale += (rowmax as f64).ln();
        }
        let ls = lnscale as f32;
        sp.n[i] = lnf(xn) + ls;
        sp.b[i] = lnf(xb) + ls;
        sp.e[i] = lnf(xe) + ls;
        sp.j[i] = lnf(xj) + ls;
        sp.c[i] = lnf(xc) + ls;
        std::mem::swap(mn, mc);
        std::mem::swap(inx, ic);
        std::mem::swap(dn, dc);
    }
    // Row 0: only N and B.
    let eo = p.mo_row(dsq[0]);
    let mut b = 0.0f32;
    for k in 1..=m {
        b += mn[k] * p.p_bm[k] * eo[k];
    }
    let ls = lnscale as f32;
    sp.b[0] = lnf(b) + ls;
    sp.n[0] = logsum(sp.n[1] + xs.n_loop, sp.b[0] + xs.n_move);
    sp.n[0]
}

// ---------------------------------------------------------------------------
// Full-matrix Forward/Backward for envelope rescoring
// ---------------------------------------------------------------------------

/// Full DP matrix in probability space (rows `0..=L`, columns `0..=M`), with
/// the per-row cumulative ln scale (`true = stored · e^{lnscale[i]}`).
#[derive(Default)]
pub struct Matrix {
    pub m: usize,
    pub l: usize,
    pub mm: Vec<f32>,
    pub im: Vec<f32>,
    pub dm: Vec<f32>,
    pub lnscale: Vec<f64>,
    pub sp: Specials,
    /// Scaled (stored) special-state values and the per-row rescaling
    /// factor (1.0 when the row was not rescaled), HMMER's `xmx[i][SCALE]`.
    /// Filled by [`forward_matrix`] only; used by [`stochastic_trace`].
    pub xn: Vec<f32>,
    pub xb: Vec<f32>,
    pub xe: Vec<f32>,
    pub xj: Vec<f32>,
    pub xc: Vec<f32>,
    pub scale: Vec<f32>,
}

impl Matrix {
    fn reset(&mut self, m: usize, l: usize) {
        self.m = m;
        self.l = l;
        let n = (l + 1) * (m + 1);
        for v in [&mut self.mm, &mut self.im, &mut self.dm] {
            v.clear();
            v.resize(n, 0.0);
        }
        self.lnscale.clear();
        self.lnscale.resize(l + 2, 0.0);
        self.sp.reset(l);
        for v in [&mut self.xn, &mut self.xb, &mut self.xe, &mut self.xj, &mut self.xc] {
            v.clear();
            v.resize(l + 1, 0.0);
        }
        self.scale.clear();
        self.scale.resize(l + 1, 1.0);
    }
}

/// Forward over `dsq` keeping the full matrix.
pub fn forward_matrix(p: &Profile, dsq: &[u8], xs: &Xsc, fx: &mut Matrix) -> f32 {
    let m = p.m;
    let l = dsq.len();
    fx.reset(m, l);
    if l == 0 || m == 0 {
        return NEG_INF;
    }
    let xp = xs.probs();
    let w = m + 1;
    let mut xn = 1.0f32;
    let mut xb = xp.n_move;
    let mut xj = 0.0f32;
    let mut xc = 0.0f32;
    let mut lnscale = 0.0f64;
    fx.sp.n[0] = 0.0;
    fx.sp.b[0] = lnf(xb);
    fx.xn[0] = xn;
    fx.xb[0] = xb;
    for i in 1..=l {
        let x = dsq[i - 1];
        let eo = p.mo_row(x);
        let ins = if x != CODE_STAR { 1.0f32 } else { 0.0 };
        let (pb, cb) = ((i - 1) * w, i * w);
        let mut xe = 0.0f32;
        for k in 1..=m {
            // HMMER's summation order: B, M, I, D.
            let mv = (xb * p.p_bm[k]
                + fx.mm[pb + k - 1] * p.p_mm[k - 1]
                + fx.im[pb + k - 1] * p.p_im[k - 1]
                + fx.dm[pb + k - 1] * p.p_dm[k - 1])
                * eo[k];
            fx.mm[cb + k] = mv;
            fx.im[cb + k] = if k < m { (fx.mm[pb + k] * p.p_mi[k] + fx.im[pb + k] * p.p_ii[k]) * ins } else { 0.0 };
            let dv = fx.mm[cb + k - 1] * p.p_md[k - 1] + fx.dm[cb + k - 1] * p.p_dd[k - 1];
            fx.dm[cb + k] = dv;
            xe += mv + dv;
        }
        xj = xj * xp.j_loop + xe * xp.e_loop;
        xc = xc * xp.c_loop + xe * xp.e_move;
        xn *= xp.n_loop;
        xb = xn * xp.n_move + xj * xp.j_move;
        if xe > RESCALE {
            let inv = 1.0 / xe;
            for k in 1..=m {
                fx.mm[cb + k] *= inv;
                fx.im[cb + k] *= inv;
                fx.dm[cb + k] *= inv;
            }
            xn *= inv;
            xb *= inv;
            xj *= inv;
            xc *= inv;
            lnscale += (xe as f64).ln();
            fx.scale[i] = xe;
            xe = 1.0;
        }
        fx.xn[i] = xn;
        fx.xb[i] = xb;
        fx.xe[i] = xe;
        fx.xj[i] = xj;
        fx.xc[i] = xc;
        fx.lnscale[i] = lnscale;
        let ls = lnscale as f32;
        let s = &mut fx.sp;
        s.n[i] = lnf(xn) + ls;
        s.b[i] = lnf(xb) + ls;
        s.e[i] = lnf(xe) + ls;
        s.j[i] = lnf(xj) + ls;
        s.c[i] = lnf(xc) + ls;
    }
    (lnf(xc * xp.c_move) as f64 + lnscale) as f32
}

/// Backward over `dsq` keeping the full matrix.
pub fn backward_matrix(p: &Profile, dsq: &[u8], xs: &Xsc, bx: &mut Matrix) -> f32 {
    let m = p.m;
    let l = dsq.len();
    bx.reset(m, l);
    if l == 0 || m == 0 {
        return NEG_INF;
    }
    let xp = xs.probs();
    let w = m + 1;
    let mut xc = xp.c_move;
    let mut xe = xc * xp.e_move;
    let mut xj = 0.0f32;
    let mut xn = 0.0f32;
    let base = l * w;
    bx.mm[base + m] = xe;
    bx.dm[base + m] = xe;
    for k in (1..m).rev() {
        bx.mm[base + k] = xe + bx.dm[base + k + 1] * p.p_md[k];
        bx.dm[base + k] = xe + bx.dm[base + k + 1] * p.p_dd[k];
    }
    let mut lnscale = 0.0f64;
    bx.lnscale[l] = 0.0;
    bx.sp.c[l] = lnf(xc);
    bx.sp.e[l] = lnf(xe);
    for i in (1..l).rev() {
        let xnext = dsq[i];
        let eo = p.mo_row(xnext);
        let ins = if xnext != CODE_STAR { 1.0f32 } else { 0.0 };
        let (nb, cb) = ((i + 1) * w, i * w);
        let mut xb = 0.0f32;
        for k in 1..=m {
            xb += bx.mm[nb + k] * p.p_bm[k] * eo[k];
        }
        xj = xj * xp.j_loop + xb * xp.j_move;
        xc *= xp.c_loop;
        xe = xj * xp.e_loop + xc * xp.e_move;
        xn = xn * xp.n_loop + xb * xp.n_move;
        bx.mm[cb + m] = xe;
        bx.dm[cb + m] = xe;
        bx.im[cb + m] = 0.0;
        let mut rowmax = xb.max(xe).max(xn).max(xj).max(xc);
        for k in (1..m).rev() {
            let mnext = bx.mm[nb + k + 1] * eo[k + 1];
            let iv = bx.im[nb + k] * ins;
            let mv = mnext * p.p_mm[k] + iv * p.p_mi[k] + xe + bx.dm[cb + k + 1] * p.p_md[k];
            bx.mm[cb + k] = mv;
            bx.im[cb + k] = mnext * p.p_im[k] + iv * p.p_ii[k];
            bx.dm[cb + k] = mnext * p.p_dm[k] + bx.dm[cb + k + 1] * p.p_dd[k] + xe;
            if mv > rowmax {
                rowmax = mv;
            }
        }
        if rowmax > RESCALE {
            let inv = 1.0 / rowmax;
            for k in 1..=m {
                bx.mm[cb + k] *= inv;
                bx.im[cb + k] *= inv;
                bx.dm[cb + k] *= inv;
            }
            xb *= inv;
            xe *= inv;
            xn *= inv;
            xj *= inv;
            xc *= inv;
            lnscale += (rowmax as f64).ln();
        }
        bx.lnscale[i] = lnscale;
        let ls = lnscale as f32;
        let s = &mut bx.sp;
        s.n[i] = lnf(xn) + ls;
        s.b[i] = lnf(xb) + ls;
        s.e[i] = lnf(xe) + ls;
        s.j[i] = lnf(xj) + ls;
        s.c[i] = lnf(xc) + ls;
    }
    let eo = p.mo_row(dsq[0]);
    let mut b = 0.0f32;
    for k in 1..=m {
        b += bx.mm[w + k] * p.p_bm[k] * eo[k];
    }
    bx.lnscale[0] = lnscale;
    let s = &mut bx.sp;
    s.b[0] = lnf(b) + lnscale as f32;
    s.n[0] = logsum(s.n[1] + xs.n_loop, s.b[0] + xs.n_move);
    s.n[0]
}

/// Null2 odds ratios by posterior expectation (`p7_GNull2_ByExpectation`),
/// returned as `ln(null2[x])` for every digital code.
pub fn null2_by_expectation(p: &Profile, dsq: &[u8], xs: &Xsc, fx: &Matrix, bx: &Matrix) -> [f32; NCODE] {
    let m = p.m;
    let l = dsq.len();
    let w = m + 1;
    let overall = bx.sp.n[0] as f64;
    let mut pm = vec![0f64; w];
    let mut pi = vec![0f64; w];
    let mut pncj = 0f64;
    let mut rm = vec![0f64; w];
    let mut ri = vec![0f64; w];
    for i in 1..=l {
        let row = i * w;
        let fac = (fx.lnscale[i] + bx.lnscale[i] - overall).exp();
        let mut denom = 0f64;
        for k in 1..=m {
            let vm = fx.mm[row + k] as f64 * bx.mm[row + k] as f64 * fac;
            rm[k] = vm;
            denom += vm;
            let vi = if k < m { fx.im[row + k] as f64 * bx.im[row + k] as f64 * fac } else { 0.0 };
            ri[k] = vi;
            denom += vi;
        }
        let n = ((fx.sp.n[i - 1] + bx.sp.n[i] + xs.n_loop) as f64 - overall).exp();
        let j = ((fx.sp.j[i - 1] + bx.sp.j[i] + xs.j_loop) as f64 - overall).exp();
        let c = ((fx.sp.c[i - 1] + bx.sp.c[i] + xs.c_loop) as f64 - overall).exp();
        denom += n + j + c;
        if denom > 0.0 && denom.is_finite() {
            let inv = 1.0 / denom;
            for k in 1..=m {
                pm[k] += rm[k] * inv;
                pi[k] += ri[k] * inv;
            }
            pncj += (n + j + c) * inv;
        }
    }
    let ld = l.max(1) as f64;
    let mut null2 = [0f64; NCODE];
    let isum: f64 = pi.iter().sum::<f64>() / ld;
    let xfactor = pncj / ld;
    for (x, slot) in null2.iter_mut().enumerate().take(K) {
        let row = p.mo_row(x as u8);
        let mut v = 0f64;
        for k in 1..=m {
            v += (pm[k] / ld) * row[k] as f64;
        }
        *slot = v + isum + xfactor;
    }
    // Degenerate codes: plain average of member odds (esl_abc_FAvgScVec).
    for code in [20u8, 21, 22, 23] {
        let mem = degenerate_members(code);
        let s: f64 = mem.iter().map(|&y| null2[y]).sum();
        null2[code as usize] = s / mem.len() as f64;
    }
    null2[CODE_STAR as usize] = 1.0;
    let mut out = [0f32; NCODE];
    for x in 0..NCODE {
        out[x] = if null2[x] > 0.0 { null2[x].ln() as f32 } else { NEG_INF };
    }
    out
}
// ---------------------------------------------------------------------------
// Viterbi with traceback
// ---------------------------------------------------------------------------

/// One aligned domain from a traceback (1-based coordinates on the scored
/// subsequence and on the model).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AliCoords {
    pub iali: usize,
    pub jali: usize,
    pub kali: usize,
    pub mali: usize,
}

const P_MM: u8 = 0;
const P_IM: u8 = 1;
const P_DM: u8 = 2;
const P_BM: u8 = 3;
const P_MI: u8 = 0;
const P_II: u8 = 1;
const P_MD: u8 = 0;
const P_DD: u8 = 1;

/// Optimal alignment(s) of `dsq` to `p` under the special-state config `xs`;
/// returns one [`AliCoords`] per B…E segment of the Viterbi path.
#[allow(clippy::neg_cmp_op_on_partial_ord)] // NaN-aware on purpose
pub fn viterbi_trace(p: &Profile, dsq: &[u8], xs: &Xsc) -> Vec<AliCoords> {
    let m = p.m;
    let l = dsq.len();
    if l == 0 || m == 0 {
        return Vec::new();
    }
    let w = m + 1;
    let n = (l + 1) * w;
    let mut vm = vec![NEG_INF; n];
    let mut vi = vec![NEG_INF; n];
    let mut vd = vec![NEG_INF; n];
    let mut pm = vec![0u8; n];
    let mut pi = vec![0u8; n];
    let mut pd = vec![0u8; n];
    let mut xn = vec![NEG_INF; l + 1];
    let mut xb = vec![NEG_INF; l + 1];
    let mut xe = vec![NEG_INF; l + 1];
    let mut xj = vec![NEG_INF; l + 1];
    let mut xc = vec![NEG_INF; l + 1];
    // E source per row: (k, is_delete); B from J?; J from E?; C from E?
    let mut e_src = vec![(0usize, false); l + 1];
    let mut b_from_j = vec![false; l + 1];
    let mut j_from_e = vec![false; l + 1];
    let mut c_from_e = vec![false; l + 1];

    xn[0] = 0.0;
    xb[0] = xs.n_move;
    for i in 1..=l {
        let x = dsq[i - 1];
        let rs = p.msc_row(x);
        let ins_ok = x != CODE_STAR;
        let mut best_e = NEG_INF;
        let mut best_src = (0usize, false);
        for k in 1..=m {
            let pc = (i - 1) * w + k;
            let cc = i * w + k;
            let c0 = vm[pc - 1] + p.t_mm[k - 1];
            let c1 = vi[pc - 1] + p.t_im[k - 1];
            let c2 = vd[pc - 1] + p.t_dm[k - 1];
            let c3 = xb[i - 1] + p.t_bm[k];
            let (mut sc, mut ptr) = (c0, P_MM);
            if c1 > sc {
                sc = c1;
                ptr = P_IM;
            }
            if c2 > sc {
                sc = c2;
                ptr = P_DM;
            }
            if c3 > sc {
                sc = c3;
                ptr = P_BM;
            }
            vm[cc] = sc + rs[k];
            pm[cc] = ptr;
            if k < m && ins_ok {
                let a = vm[pc] + p.t_mi[k];
                let b = vi[pc] + p.t_ii[k];
                if a >= b {
                    vi[cc] = a;
                    pi[cc] = P_MI;
                } else {
                    vi[cc] = b;
                    pi[cc] = P_II;
                }
            }
            let a = vm[cc - 1] + p.t_md[k - 1];
            let b = vd[cc - 1] + p.t_dd[k - 1];
            if a >= b {
                vd[cc] = a;
                pd[cc] = P_MD;
            } else {
                vd[cc] = b;
                pd[cc] = P_DD;
            }
            if vm[cc] > best_e {
                best_e = vm[cc];
                best_src = (k, false);
            }
            if vd[cc] > best_e {
                best_e = vd[cc];
                best_src = (k, true);
            }
        }
        xe[i] = best_e;
        e_src[i] = best_src;
        let jl = xj[i - 1] + xs.j_loop;
        let je = best_e + xs.e_loop;
        if je > jl {
            xj[i] = je;
            j_from_e[i] = true;
        } else {
            xj[i] = jl;
        }
        let cl = xc[i - 1] + xs.c_loop;
        let ce = best_e + xs.e_move;
        if ce > cl {
            xc[i] = ce;
            c_from_e[i] = true;
        } else {
            xc[i] = cl;
        }
        xn[i] = xn[i - 1] + xs.n_loop;
        let bn = xn[i] + xs.n_move;
        let bj = xj[i] + xs.j_move;
        if bj > bn {
            xb[i] = bj;
            b_from_j[i] = true;
        } else {
            xb[i] = bn;
        }
    }
    if !(xc[l] > NEG_INF) {
        return Vec::new();
    }

    // Traceback.
    #[derive(Clone, Copy, PartialEq)]
    enum St {
        C,
        E,
        M(usize),
        I(usize),
        D(usize),
        B,
        J,
        N,
    }
    let mut out = Vec::new();
    let mut st = St::C;
    let mut i = l;
    let mut cur_end: Option<(usize, usize)> = None; // (jali, mali)
    let mut guard = 0usize;
    loop {
        guard += 1;
        if guard > 4 * (l + 2) * (m + 2) {
            break; // defensive: never loop forever on malformed models
        }
        match st {
            St::C => {
                if i == 0 {
                    break;
                }
                if c_from_e[i] {
                    st = St::E;
                } else {
                    i -= 1;
                }
            }
            St::E => {
                let (k, del) = e_src[i];
                if k == 0 {
                    break;
                }
                st = if del { St::D(k) } else { St::M(k) };
                cur_end = Some((i, k));
            }
            St::M(k) => {
                let ptr = pm[i * w + k];
                match ptr {
                    P_MM => {
                        st = St::M(k - 1);
                        i -= 1;
                    }
                    P_IM => {
                        st = St::I(k - 1);
                        i -= 1;
                    }
                    P_DM => {
                        st = St::D(k - 1);
                        i -= 1;
                    }
                    _ => {
                        if let Some((jali, mali)) = cur_end.take() {
                            out.push(AliCoords { iali: i, jali, kali: k, mali });
                        }
                        st = St::B;
                        i -= 1;
                    }
                }
            }
            St::I(k) => {
                if pi[i * w + k] == P_MI {
                    st = St::M(k);
                } else {
                    st = St::I(k);
                }
                i -= 1;
            }
            St::D(k) => {
                if pd[i * w + k] == P_MD {
                    st = St::M(k - 1);
                } else {
                    st = St::D(k - 1);
                }
            }
            St::B => {
                st = if b_from_j[i] { St::J } else { St::N };
            }
            St::J => {
                if j_from_e[i] {
                    st = St::E;
                } else {
                    i -= 1;
                }
            }
            St::N => {
                break;
            }
        }
    }
    out.reverse();
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hmmer::parse::parse_hmm_str;

    fn tiny() -> Profile {
        let txt = include_str!("testdata/tiny.hmm");
        parse_hmm_str(txt).unwrap()[0].configure()
    }

    #[test]
    fn forward_equals_backward() {
        let p = tiny();
        let seq = crate::hmmer::alphabet::encode(b"GGMKAGGGMKAGG");
        let xs = Xsc::new(seq.len(), true);
        let mut buf = DpBuf::default();
        let mut sf = Specials::default();
        let mut sb = Specials::default();
        let f = forward(&p, &seq, &xs, &mut buf, Some(&mut sf));
        let b = backward(&p, &seq, &xs, &mut buf, &mut sb);
        assert!((f - b).abs() < 1e-3, "fwd {f} bwd {b}");
        let mut fx = Matrix::default();
        let mut bx = Matrix::default();
        let f2 = forward_matrix(&p, &seq, &xs, &mut fx);
        let b2 = backward_matrix(&p, &seq, &xs, &mut bx);
        assert!((f2 - f).abs() < 1e-4 && (b2 - b).abs() < 1e-3);
    }

    #[test]
    fn viterbi_le_forward_and_trace_finds_two_domains() {
        let p = tiny();
        let seq = crate::hmmer::alphabet::encode(b"GGMKAGGGGGGMKAGG");
        let xs = Xsc::new(seq.len(), true);
        let mut buf = DpBuf::default();
        let v = viterbi_score(&p, &seq, &xs, &mut buf);
        let f = forward(&p, &seq, &xs, &mut buf, None);
        assert!(v <= f + 1e-4);
        let tr = viterbi_trace(&p, &seq, &xs);
        assert_eq!(tr.len(), 2, "{tr:?}");
        assert_eq!((tr[0].iali, tr[0].jali), (3, 5));
        assert_eq!((tr[1].iali, tr[1].jali), (12, 14));
        assert_eq!((tr[0].kali, tr[0].mali), (1, 3));
    }

    #[test]
    fn logsum_basics() {
        assert!((logsum(0.0, 0.0) - std::f32::consts::LN_2).abs() < 1e-6);
        assert_eq!(logsum(NEG_INF, -2.0), -2.0);
        assert_eq!(logsum(NEG_INF, NEG_INF), NEG_INF);
    }
}

// ---------------------------------------------------------------------------
// Stochastic traceback (HMMER `p7_StochasticTrace`) with Easel's RNG
// ---------------------------------------------------------------------------

/// Easel's "fast" RNG (`esl_randomness_CreateFast`): Knuth's LCG
/// `x = 69069·x + 1 (mod 2^32)`, seeded through `esl_mix3`.
#[derive(Clone, Debug)]
pub struct EslRng {
    seed: u32,
    x: u32,
}

fn esl_mix3(mut a: u32, mut b: u32, mut c: u32) -> u32 {
    a = a.wrapping_sub(b);
    a = a.wrapping_sub(c);
    a ^= c >> 13;
    b = b.wrapping_sub(c);
    b = b.wrapping_sub(a);
    b ^= a << 8;
    c = c.wrapping_sub(a);
    c = c.wrapping_sub(b);
    c ^= b >> 13;
    a = a.wrapping_sub(b);
    a = a.wrapping_sub(c);
    a ^= c >> 12;
    b = b.wrapping_sub(c);
    b = b.wrapping_sub(a);
    b ^= a << 16;
    c = c.wrapping_sub(a);
    c = c.wrapping_sub(b);
    c ^= b >> 5;
    a = a.wrapping_sub(b);
    a = a.wrapping_sub(c);
    a ^= c >> 3;
    b = b.wrapping_sub(c);
    b = b.wrapping_sub(a);
    b ^= a << 10;
    c = c.wrapping_sub(a);
    c = c.wrapping_sub(b);
    c ^= b >> 15;
    c
}

impl EslRng {
    /// `esl_randomness_CreateFast(seed)`; hmmsearch uses seed 42.
    pub fn new(seed: u32) -> EslRng {
        let mut r = EslRng { seed, x: 0 };
        r.reseed();
        r
    }

    /// Reset to the original seed (`esl_randomness_Init(r, GetSeed(r))`).
    pub fn reseed(&mut self) {
        self.x = esl_mix3(self.seed, 87654321, 12345678);
        if self.x == 0 {
            self.x = 42;
        }
    }

    /// `esl_random()`: uniform in [0, 1).
    #[inline]
    pub fn next(&mut self) -> f64 {
        self.x = self.x.wrapping_mul(69069).wrapping_add(1);
        self.x as f64 / 4294967296.0
    }
}

/// `esl_vec_FNorm` (Kahan-summed normalisation).
fn esl_fnorm(v: &mut [f32]) {
    let mut sum = 0f32;
    let mut c = 0f32;
    for &x in v.iter() {
        let y = x - c;
        let t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    if sum != 0.0 {
        for x in v.iter_mut() {
            *x /= sum;
        }
    } else {
        let n = v.len() as f32;
        for x in v.iter_mut() {
            *x = 1.0 / n;
        }
    }
}

/// `esl_rnd_FChoose`.
fn esl_fchoose(rng: &mut EslRng, p: &[f32]) -> usize {
    let norm: f64 = p.iter().map(|&x| x as f64).sum();
    let roll = rng.next();
    let mut sum = 0f64;
    for (i, &x) in p.iter().enumerate() {
        sum += x as f64;
        if roll < sum / norm {
            return i;
        }
    }
    p.len() - 1
}

/// One domain of a sampled trace (1-based coordinates on the region).
#[derive(Clone, Debug, Default)]
pub struct SampledDomain {
    pub sqfrom: usize,
    pub sqto: usize,
    pub hmmfrom: usize,
    pub hmmto: usize,
    /// Model node of every M or I emission in the domain (HMMER's
    /// `p7_Null2_ByTrace` counts I emissions at node k as M_k).
    pub nodes: Vec<usize>,
}

/// Sample one trace from a Forward matrix (`p7_StochasticTrace`, SSE
/// version's candidate order), returning its domains in sequence order.
pub fn stochastic_trace(p: &Profile, fx: &Matrix, xs: &Xsc, rng: &mut EslRng) -> Vec<SampledDomain> {
    #[derive(Clone, Copy, PartialEq)]
    enum St {
        C,
        E,
        M,
        I,
        D,
        B,
        J,
        N,
    }
    let xp = xs.probs();
    let m = p.m;
    let w = m + 1;
    let q_stripes = std::cmp::max(2, (m.saturating_sub(1)) / 4 + 1);
    let mut doms: Vec<SampledDomain> = Vec::new();
    let mut cur: Option<SampledDomain> = None;
    let mut s0 = St::C;
    let mut i = fx.l;
    let mut k = 0usize;
    let mut guard = 0usize;
    let limit = 8 * (fx.l + 2) * (m + 2) + 64;
    loop {
        guard += 1;
        if guard > limit {
            break;
        }
        let s1 = match s0 {
            St::C => {
                if i == 0 {
                    break;
                }
                let mut path = [fx.xc[i - 1] * xp.c_loop, fx.xe[i] * xp.e_move * fx.scale[i]];
                esl_fnorm(&mut path);
                if esl_fchoose(rng, &path) == 0 {
                    St::C
                } else {
                    St::E
                }
            }
            St::E => {
                let roll = rng.next();
                let norm = (1.0f64 / fx.xe[i] as f64) as f32;
                let row = i * w;
                let mut sum = 0f64;
                let mut pick: Option<(usize, St)> = None;
                let mut passes = 0;
                'outer: while pick.is_none() && passes < 1000 {
                    passes += 1;
                    for q in 0..q_stripes {
                        for r in 0..4 {
                            let kk = r * q_stripes + q + 1;
                            if kk <= m {
                                sum += (fx.mm[row + kk] * norm) as f64;
                                if roll < sum {
                                    pick = Some((kk, St::M));
                                    break 'outer;
                                }
                            }
                        }
                        for r in 0..4 {
                            let kk = r * q_stripes + q + 1;
                            if kk <= m {
                                sum += (fx.dm[row + kk] * norm) as f64;
                                if roll < sum {
                                    pick = Some((kk, St::D));
                                    break 'outer;
                                }
                            }
                        }
                    }
                    if sum <= 0.0 {
                        break;
                    }
                }
                match pick {
                    Some((kk, st)) => {
                        k = kk;
                        cur = Some(SampledDomain::default());
                        st
                    }
                    None => break,
                }
            }
            St::M => {
                if let Some(d) = cur.as_mut() {
                    if d.sqto == 0 {
                        d.sqto = i;
                        d.hmmto = k;
                    }
                    d.sqfrom = i;
                    d.hmmfrom = k;
                    d.nodes.push(k);
                }
                let prow = (i - 1) * w;
                let (mpv, ipv, dpv, tmm, tim, tdm) = if k >= 2 {
                    (fx.mm[prow + k - 1], fx.im[prow + k - 1], fx.dm[prow + k - 1], p.p_mm[k - 1], p.p_im[k - 1], p.p_dm[k - 1])
                } else {
                    (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
                };
                let mut path = [fx.xb[i - 1] * p.p_bm[k], mpv * tmm, ipv * tim, dpv * tdm];
                esl_fnorm(&mut path);
                let c = esl_fchoose(rng, &path);
                k -= 1;
                i -= 1;
                match c {
                    0 => St::B,
                    1 => St::M,
                    2 => St::I,
                    _ => St::D,
                }
            }
            St::D => {
                let row = i * w;
                let (mpv, dpv, tmd, tdd) = if k >= 2 {
                    (fx.mm[row + k - 1], fx.dm[row + k - 1], p.p_md[k - 1], p.p_dd[k - 1])
                } else {
                    (0.0, 0.0, 0.0, 0.0)
                };
                let mut path = [mpv * tmd, dpv * tdd];
                esl_fnorm(&mut path);
                let c = esl_fchoose(rng, &path);
                k -= 1;
                if c == 0 {
                    St::M
                } else {
                    St::D
                }
            }
            St::I => {
                if let Some(d) = cur.as_mut() {
                    d.nodes.push(k);
                }
                let prow = (i - 1) * w;
                let mut path = [fx.mm[prow + k] * p.p_mi[k], fx.im[prow + k] * p.p_ii[k]];
                esl_fnorm(&mut path);
                let c = esl_fchoose(rng, &path);
                i -= 1;
                if c == 0 {
                    St::M
                } else {
                    St::I
                }
            }
            St::B => {
                if let Some(d) = cur.take() {
                    doms.push(d);
                }
                let mut path = [fx.xn[i] * xp.n_move, fx.xj[i] * xp.j_move];
                esl_fnorm(&mut path);
                if esl_fchoose(rng, &path) == 0 {
                    St::N
                } else {
                    St::J
                }
            }
            St::J => {
                if i == 0 {
                    break;
                }
                let mut path = [fx.xj[i - 1] * xp.j_loop, fx.xe[i] * xp.e_loop * fx.scale[i]];
                esl_fnorm(&mut path);
                if esl_fchoose(rng, &path) == 0 {
                    St::J
                } else {
                    St::E
                }
            }
            St::N => {
                if i == 0 {
                    break;
                }
                St::N
            }
        };
        if matches!(s1, St::N | St::J | St::C) && s1 == s0 {
            i -= 1;
        }
        s0 = s1;
    }
    doms.reverse();
    doms
}

// ---------------------------------------------------------------------------
// Posterior decoding + optimal-accuracy alignment (HMMER `p7_Decoding`,
// `p7_OptimalAccuracy`, `p7_OATrace`, SSE semantics incl. tie-breaking)
// ---------------------------------------------------------------------------

/// Posterior probabilities of each emitting state, per residue.
#[derive(Default)]
pub struct Posterior {
    pub m: usize,
    pub l: usize,
    pub pm: Vec<f32>,
    pub pi: Vec<f32>,
    pub pn: Vec<f32>,
    pub pj: Vec<f32>,
    pub pc: Vec<f32>,
}

/// SSE `p7_Decoding` semantics (no per-row renormalisation).
pub fn decode(p: &Profile, xs: &Xsc, fx: &Matrix, bx: &Matrix, pp: &mut Posterior) {
    let m = p.m;
    let l = fx.l;
    let w = m + 1;
    pp.m = m;
    pp.l = l;
    for v in [&mut pp.pm, &mut pp.pi] {
        v.clear();
        v.resize((l + 1) * w, 0.0);
    }
    for v in [&mut pp.pn, &mut pp.pj, &mut pp.pc] {
        v.clear();
        v.resize(l + 1, 0.0);
    }
    let overall = bx.sp.n[0] as f64;
    for i in 1..=l {
        let row = i * w;
        let fac = (fx.lnscale[i] + bx.lnscale[i] - overall).exp();
        for k in 1..=m {
            pp.pm[row + k] = (fx.mm[row + k] as f64 * bx.mm[row + k] as f64 * fac) as f32;
            if k < m {
                pp.pi[row + k] = (fx.im[row + k] as f64 * bx.im[row + k] as f64 * fac) as f32;
            }
        }
        pp.pn[i] = ((fx.sp.n[i - 1] + bx.sp.n[i] + xs.n_loop) as f64 - overall).exp() as f32;
        pp.pj[i] = ((fx.sp.j[i - 1] + bx.sp.j[i] + xs.j_loop) as f64 - overall).exp() as f32;
        pp.pc[i] = ((fx.sp.c[i - 1] + bx.sp.c[i] + xs.c_loop) as f64 - overall).exp() as f32;
    }
}

/// Null2 by expectation from a posterior matrix (SSE `p7_Null2_ByExpectation`),
/// as `ln(null2[x])` per digital code.
pub fn null2_from_posterior(p: &Profile, pp: &Posterior) -> [f32; NCODE] {
    let m = p.m;
    let l = pp.l;
    let w = m + 1;
    let mut sm = vec![0f32; w];
    let mut si = 0f32;
    let (mut sn, mut sj, mut sc) = (0f32, 0f32, 0f32);
    for i in 1..=l {
        let row = i * w;
        for k in 1..=m {
            sm[k] += pp.pm[row + k];
            si += pp.pi[row + k];
        }
        sn += pp.pn[i];
        sj += pp.pj[i];
        sc += pp.pc[i];
    }
    let norm = 1.0 / l.max(1) as f32;
    let xfactor = sn * norm + sc * norm + sj * norm;
    let isum = si * norm;
    let mut null2 = [0f32; NCODE];
    for (x, slot) in null2.iter_mut().enumerate().take(K) {
        let row = p.mo_row(x as u8);
        let mut v = 0f32;
        for k in 1..=m {
            v += sm[k] * norm * row[k];
        }
        *slot = v + isum + xfactor;
    }
    for code in [20u8, 21, 22, 23] {
        let mem = degenerate_members(code);
        let s: f32 = mem.iter().map(|&y| null2[y]).sum();
        null2[code as usize] = s / mem.len() as f32;
    }
    null2[CODE_STAR as usize] = 1.0;
    let mut out = [0f32; NCODE];
    for x in 0..NCODE {
        out[x] = if null2[x] > 0.0 { null2[x].ln() } else { NEG_INF };
    }
    out
}

/// Optimal-accuracy alignment of one envelope (unihit), returning the
/// coordinates of the first-to-last match state (HMMER's alidisplay
/// `sqfrom..sqto`, `hmmfrom..hmmto`), 1-based on the envelope.
pub fn optimal_accuracy_coords(p: &Profile, xs: &Xsc, pp: &Posterior) -> Option<AliCoords> {
    let m = p.m;
    let l = pp.l;
    if l == 0 || m == 0 {
        return None;
    }
    let w = m + 1;
    let xp = xs.probs();
    let n = (l + 1) * w;
    let mut om = vec![NEG_INF; n];
    let mut oi = vec![NEG_INF; n];
    let mut od = vec![NEG_INF; n];
    let mut xe = vec![NEG_INF; l + 1];
    let mut xn = vec![NEG_INF; l + 1];
    let mut xj = vec![NEG_INF; l + 1];
    let mut xb = vec![NEG_INF; l + 1];
    let mut xc = vec![NEG_INF; l + 1];
    xn[0] = 0.0;
    xb[0] = 0.0;
    // masked "AND" semantics: a disallowed predecessor contributes 0.0.
    let gate = |t: f32, v: f32| if t > 0.0 { v } else { 0.0 };
    for i in 1..=l {
        let (pr, cr) = ((i - 1) * w, i * w);
        let mut emax = NEG_INF;
        for k in 1..=m {
            let mut sv = gate(p.p_bm[k], xb[i - 1]);
            if k >= 2 {
                sv = sv.max(gate(p.p_mm[k - 1], om[pr + k - 1]));
                sv = sv.max(gate(p.p_im[k - 1], oi[pr + k - 1]));
                sv = sv.max(gate(p.p_dm[k - 1], od[pr + k - 1]));
            } else {
                sv = sv.max(0.0);
            }
            let mv = sv + pp.pm[cr + k];
            om[cr + k] = mv;
            emax = emax.max(mv);
            let iv = gate(p.p_mi[k], om[pr + k]).max(gate(p.p_ii[k], oi[pr + k])) + pp.pi[cr + k];
            oi[cr + k] = iv;
        }
        // D: from M(i,k-1) and D(i,k-1); D(i,1) impossible.
        od[cr + 1] = NEG_INF;
        for k in 2..=m {
            let a = gate(p.p_md[k - 1], om[cr + k - 1]);
            let b = gate(p.p_dd[k - 1], od[cr + k - 1]);
            od[cr + k] = a.max(b);
            emax = emax.max(od[cr + k]);
        }
        xe[i] = emax;
        let t1 = if xp.j_loop == 0.0 { 0.0 } else { xj[i - 1] + pp.pj[i] };
        let t2 = if xp.e_loop == 0.0 { 0.0 } else { xe[i] };
        xj[i] = t1.max(t2);
        let t1 = if xp.c_loop == 0.0 { 0.0 } else { xc[i - 1] + pp.pc[i] };
        let t2 = if xp.e_move == 0.0 { 0.0 } else { xe[i] };
        xc[i] = t1.max(t2);
        xn[i] = if xp.n_loop == 0.0 { 0.0 } else { xn[i - 1] + pp.pn[i] };
        let t1 = if xp.n_move == 0.0 { 0.0 } else { xn[i] };
        let t2 = if xp.j_move == 0.0 { 0.0 } else { xj[i] };
        xb[i] = t1.max(t2);
    }

    // ---- traceback ----
    #[derive(Clone, Copy, PartialEq)]
    enum St {
        C,
        E,
        M,
        I,
        D,
        B,
        J,
        N,
    }
    let q_stripes = std::cmp::max(2, (m - 1) / 4 + 1);
    let mut s0 = St::C;
    let mut i = l;
    let mut k = 0usize;
    let mut first: Option<(usize, usize)> = None; // (i,k) of first M in forward order
    let mut last: Option<(usize, usize)> = None;
    let mut guard = 0usize;
    let limit = 8 * (l + 2) * (m + 2) + 64;
    let dis = |t: f32, v: f32| if t == 0.0 { NEG_INF } else { v };
    loop {
        guard += 1;
        if guard > limit {
            break;
        }
        let s1 = match s0 {
            St::C => {
                if i == 0 {
                    break;
                }
                let a = if xp.c_loop == 0.0 { NEG_INF } else { xc[i - 1] + pp.pc[i] };
                let b = if xp.e_move == 0.0 { NEG_INF } else { xe[i] };
                if a > b {
                    St::C
                } else {
                    St::E
                }
            }
            St::E => {
                let row = i * w;
                let mut max = NEG_INF;
                let mut pick: Option<(usize, St)> = None;
                for q in 0..q_stripes {
                    for r in 0..4 {
                        let kk = r * q_stripes + q + 1;
                        if kk <= m && om[row + kk] >= max {
                            max = om[row + kk];
                            pick = Some((kk, St::M));
                        }
                    }
                    for r in 0..4 {
                        let kk = r * q_stripes + q + 1;
                        if kk <= m && od[row + kk] > max {
                            max = od[row + kk];
                            pick = Some((kk, St::D));
                        }
                    }
                }
                match pick {
                    Some((kk, st)) => {
                        k = kk;
                        st
                    }
                    None => break,
                }
            }
            St::M => {
                if last.is_none() {
                    last = Some((i, k));
                }
                first = Some((i, k));
                let prow = (i - 1) * w;
                let (mv, iv, dv) = if k >= 2 {
                    (
                        dis(p.p_mm[k - 1], om[prow + k - 1]),
                        dis(p.p_im[k - 1], oi[prow + k - 1]),
                        dis(p.p_dm[k - 1], od[prow + k - 1]),
                    )
                } else {
                    (NEG_INF, NEG_INF, NEG_INF)
                };
                let bv = dis(p.p_bm[k], xb[i - 1]);
                let path = [mv, iv, dv, bv];
                let mut best = 0;
                for t in 1..4 {
                    if path[t] > path[best] {
                        best = t;
                    }
                }
                k -= 1;
                i -= 1;
                match best {
                    0 => St::M,
                    1 => St::I,
                    2 => St::D,
                    _ => St::B,
                }
            }
            St::D => {
                let row = i * w;
                let (a, b) = if k >= 2 {
                    (dis(p.p_md[k - 1], om[row + k - 1]), dis(p.p_dd[k - 1], od[row + k - 1]))
                } else {
                    (NEG_INF, NEG_INF)
                };
                k -= 1;
                if a >= b {
                    St::M
                } else {
                    St::D
                }
            }
            St::I => {
                let prow = (i - 1) * w;
                let a = dis(p.p_mi[k], om[prow + k]);
                let b = dis(p.p_ii[k], oi[prow + k]);
                i -= 1;
                if a >= b {
                    St::M
                } else {
                    St::I
                }
            }
            St::B => {
                let a = if xp.n_move == 0.0 { NEG_INF } else { xn[i] };
                let b = if xp.j_move == 0.0 { NEG_INF } else { xj[i] };
                if a > b {
                    St::N
                } else {
                    // A domain begins here; unihit traces have one domain.
                    break;
                }
            }
            St::J => break,
            St::N => break,
        };
        if matches!(s1, St::N | St::J | St::C) && s1 == s0 {
            i -= 1;
        }
        if s1 == St::N {
            break; // reached the leading flank: the domain is complete
        }
        s0 = s1;
    }
    match (first, last) {
        (Some((i1, k1)), Some((i2, k2))) => Some(AliCoords { iali: i1, jali: i2, kali: k1, mali: k2 }),
        _ => None,
    }
}
