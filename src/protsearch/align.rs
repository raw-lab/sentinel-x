//! Seed verification, ungapped extension and banded Smith-Waterman, following
//! DIAMOND 2.1's default-sensitivity code paths:
//!
//! * [`fingerprint_ok`] — stage-1 filter: ≥ 11 identical residues in the 48
//!   positions `[−16, +32)` around the seed (`search/finger_print.h`);
//! * [`window_ungapped`] — stage-2 filter: best ungapped local segment in a
//!   ±48 window (`dp/ungapped_align.cpp:ungapped_window`);
//! * [`xdrop_bias`] — X-drop ungapped extension with the composition bias
//!   (`xdrop_ungapped`, X = 12.3 bits = 20 raw);
//! * [`banded_score`] / [`banded_align`] — local alignment restricted to a
//!   diagonal band (`BandedSwipe`): BLOSUM62 + per-query-position bias, gap
//!   open 11 / extend 1, target-major iteration, ties resolved as DIAMOND's
//!   cell update does (diagonal, then target gap, then query gap; extension
//!   preferred over opening; first maximum wins).

use super::alphabet::{score, NCODES};

/// Cost of the first residue of a gap (BLAST open 11 + extend 1).
pub const GAP_FIRST: i32 = 12;
/// Cost of each further gap residue.
pub const GAP_EXT: i32 = 1;
/// DIAMOND `--xdrop 12.3` bits as a raw score.
pub const XDROP: i32 = 20;
/// DIAMOND `hamming_filter_id` for default sensitivity.
pub const MIN_FINGERPRINT_ID: u32 = 11;
/// DIAMOND `--ungapped-window`.
pub const UNGAPPED_WINDOW: usize = 48;

const NEG: i32 = i32::MIN / 4;

/// DIAMOND's sequence delimiter letter (one between consecutive sequences
/// of a set, 256 around the set).
pub const DELIM: u8 = 31;

/// Letter at relative position `p` of sequence `k` in a DIAMOND sequence
/// set laid out as `[delims] s0 D s1 D s2 … [delims]`: positions before or
/// after a sequence read the delimiter and then the neighbouring sequences,
/// exactly what DIAMOND's fixed-width SIMD loads see.
#[inline]
pub fn set_letter<'a, G: Fn(usize) -> &'a [u8]>(get: G, n: usize, mut k: usize, mut p: isize) -> u8 {
    loop {
        let s = get(k);
        let l = s.len() as isize;
        if p >= 0 && p < l {
            return s[p as usize];
        }
        if p == -1 || p == l {
            return DELIM;
        }
        if p < -1 {
            if k == 0 {
                return DELIM;
            }
            k -= 1;
            p += get(k).len() as isize + 1;
        } else {
            if k + 1 >= n {
                return DELIM;
            }
            p -= l + 1;
            k += 1;
        }
    }
}

/// Stage-1 fingerprint (`Byte_finger_print_48`): identical letters among the
/// 48 positions `[−16, +32)` around the seed starts `qp`/`sp`; positions
/// outside the sequences are read through `qx`/`sx` (delimiters and
/// neighbouring sequences, see [`set_letter`]).
pub fn fingerprint_ok<FQ: Fn(isize) -> u8, FS: Fn(isize) -> u8>(
    q: &[u8],
    qp: usize,
    s: &[u8],
    sp: usize,
    qx: FQ,
    sx: FS,
) -> bool {
    let mut n = 0u32;
    if qp >= 16 && qp + 32 <= q.len() && sp >= 16 && sp + 32 <= s.len() {
        let a = &q[qp - 16..qp + 32];
        let b = &s[sp - 16..sp + 32];
        for t in 0..48 {
            n += (a[t] == b[t]) as u32;
        }
    } else {
        let (ql, sl) = (q.len() as isize, s.len() as isize);
        for t in -16isize..32 {
            let a = qp as isize + t;
            let b = sp as isize + t;
            let ca = if a >= 0 && a < ql { q[a as usize] } else { qx(a) };
            let cb = if b >= 0 && b < sl { s[b as usize] } else { sx(b) };
            n += (ca == cb) as u32;
        }
    }
    n >= MIN_FINGERPRINT_ID
}

/// Stage-2 window score (`ungapped_window`): Kadane maximum of ungapped
/// scores along the seed's diagonal over the query window `[qp−48, qp+48)`
/// clipped to the query; target positions past its ends read the delimiter
/// (score −128) and the neighbouring sequences through `sx`.
pub fn window_ungapped<FS: Fn(isize) -> u8>(q: &[u8], qp: usize, s: &[u8], sp: usize, sx: FS) -> i32 {
    let w = UNGAPPED_WINDOW;
    let q0 = qp.saturating_sub(w);
    let q1 = (qp + w).min(q.len());
    let mut best = 0i32;
    let mut st = 0i32;
    let sl = s.len() as isize;
    for i in q0..q1 {
        let j = sp as isize + (i as isize - qp as isize);
        let c = if j >= 0 && j < sl { s[j as usize] } else { sx(j) };
        let v = if c == DELIM { -128 } else { score(q[i], c) };
        st = (st + v).max(0);
        if st > best {
            best = st;
        }
    }
    best
}

/// An ungapped diagonal segment (0-based query start `i`, subject start `j`).
#[derive(Clone, Copy, Debug)]
pub struct Segment {
    pub i: usize,
    pub j: usize,
    pub len: usize,
    pub score: i32,
}

impl Segment {
    /// Diagonal `i − j` (DIAMOND's convention).
    pub fn diag(&self) -> isize {
        self.i as isize - self.j as isize
    }
    pub fn subject_end(&self) -> usize {
        self.j + self.len
    }
}

/// X-drop ungapped extension around the anchor `(qa, sa)` with the int8
/// composition bias (left part ends at `qa−1`, right part starts at `qa`).
pub fn xdrop_bias(q: &[u8], bias: &[i8], s: &[u8], qa: usize, sa: usize) -> Segment {
    let xdrop = XDROP;
    let (mut score_best, mut st) = (0i32, 0i32);
    let mut delta = 0usize;
    let mut n = 1usize;
    let (mut i, mut j) = (qa as isize - 1, sa as isize - 1);
    while score_best - st < xdrop && i >= 0 && j >= 0 {
        st += score(q[i as usize], s[j as usize]) + bias[i as usize] as i32;
        if st > score_best {
            score_best = st;
            delta = n;
        }
        i -= 1;
        j -= 1;
        n += 1;
    }
    let (mut i, mut j) = (qa, sa);
    st = score_best;
    n = 1;
    let mut len = 0usize;
    while score_best - st < xdrop && i < q.len() && j < s.len() {
        st += score(q[i], s[j]) + bias[i] as i32;
        if st > score_best {
            score_best = st;
            len = n;
        }
        i += 1;
        j += 1;
        n += 1;
    }
    Segment { i: qa - delta, j: sa - delta, len: len + delta, score: score_best }
}

/// Query profile with DIAMOND's composition bias as its banded SWIPE kernel
/// applies it: `p[c * m + i] = BLOSUM62(q_i, c) + bias_i` for every target
/// code `c` (the bias vector is added per query row, whatever the target).
pub struct Profile {
    pub m: usize,
    pub p: Vec<i32>,
}

impl Profile {
    pub fn new(q: &[u8], bias: Option<&[i8]>) -> Profile {
        let m = q.len();
        let mut p = vec![0i32; NCODES * m];
        for c in 0..NCODES {
            let row = &mut p[c * m..(c + 1) * m];
            for (i, &a) in q.iter().enumerate() {
                let mut v = score(a, c as u8);
                if let Some(b) = bias {
                    v += b[i] as i32;
                }
                row[i] = v;
            }
        }
        Profile { m, p }
    }
    #[inline(always)]
    fn get(&self, c: u8, i: usize) -> i32 {
        self.p[(c as usize).min(NCODES - 1) * self.m + i]
    }
}

/// Diagonal band `[d0, d1)` with `d = i − j`, clamped to the matrix.
#[derive(Clone, Copy, Debug)]
pub struct Band {
    pub d0: isize,
    pub d1: isize,
}

impl Band {
    pub fn full(m: usize, n: usize) -> Band {
        Band { d0: -(n as isize - 1), d1: m as isize }
    }
    pub fn width(&self) -> usize {
        (self.d1 - self.d0).max(0) as usize
    }
}

/// Banded local alignment score.
pub fn banded_score(prof: &Profile, s: &[u8], band: Band) -> i32 {
    let m = prof.m as isize;
    let w = band.width();
    if w == 0 || s.is_empty() || m == 0 {
        return 0;
    }
    let mut hp = vec![0i32; w + 1];
    let mut ep = vec![NEG; w + 1];
    let mut hc = vec![0i32; w + 1];
    let mut ec = vec![NEG; w + 1];
    let mut best = 0i32;
    for (j, &c) in s.iter().enumerate() {
        let ji = j as isize;
        let mut f = NEG;
        let k0 = (-(ji + band.d0)).max(0) as usize; // first k with i >= 0
        let k1 = ((m - ji - band.d0).min(w as isize)).max(0) as usize; // i < m
        for k in 0..k0.min(w) {
            hc[k] = 0;
            ec[k] = NEG;
        }
        for k in k0..k1 {
            let i = (ji + band.d0 + k as isize) as usize;
            let diag = hp[k];
            let e = (ep[k + 1] - GAP_EXT).max(hp[k + 1] - GAP_FIRST);
            let mut h = diag + prof.get(c, i);
            if e > h {
                h = e;
            }
            if f > h {
                h = f;
            }
            if h < 0 {
                h = 0;
            }
            if h > best {
                best = h;
            }
            hc[k] = h;
            ec[k] = e;
            f = (f - GAP_EXT).max(h - GAP_FIRST);
        }
        for k in k1.max(k0)..w {
            hc[k] = 0;
            ec[k] = NEG;
        }
        std::mem::swap(&mut hp, &mut hc);
        std::mem::swap(&mut ep, &mut ec);
    }
    best
}

/// A traced local alignment (1-based inclusive coordinates).
#[derive(Clone, Debug, Default)]
pub struct Alignment {
    pub score: i32,
    pub qstart: usize,
    pub qend: usize,
    pub sstart: usize,
    pub send: usize,
    /// Alignment columns (matches + mismatches + gap columns).
    pub length: usize,
    pub identities: usize,
    pub mismatches: usize,
    pub gap_opens: usize,
    pub gaps: usize,
}

// Per-cell trace bits, DIAMOND `TracebackVectorMatrix` semantics.
const VGAP: u8 = 1; // H(i,j) == F_in(i,j)  (gap consuming query residues)
const HGAP: u8 = 2; // H(i,j) == E_in(i,j)  (gap consuming target residues)
const VOPEN: u8 = 4; // F_in(i,j) was opened from H(i-1,j)
const HOPEN: u8 = 8; // E_in(i,j) was opened from H(i,j-1)

/// Banded local alignment with traceback. Same recurrence as
/// [`banded_score`]; the reported path follows DIAMOND's traceback: the end
/// cell is the first target column reaching the maximum (last such query
/// row in that column); walking back, a cell equal to its incoming
/// query-gap value continues that gap, else one equal to its target-gap
/// value continues that one, else the step is diagonal; gaps end at the
/// first cell they could have been opened from.
pub fn banded_align(q: &[u8], prof: &Profile, s: &[u8], band: Band) -> Alignment {
    let m = q.len() as isize;
    let w = band.width();
    let n = s.len();
    if w == 0 || n == 0 || m == 0 {
        return Alignment::default();
    }
    let mut tr = vec![0u8; n * w];
    let mut hp = vec![0i32; w + 1];
    let mut ep = vec![NEG; w + 1];
    let mut hc = vec![0i32; w + 1];
    let mut ec = vec![NEG; w + 1];
    let (mut best, mut bi, mut bj) = (0i32, 0usize, 0usize);
    for (j, &c) in s.iter().enumerate() {
        let ji = j as isize;
        let mut f = NEG;
        let mut f_open = false;
        let k0 = (-(ji + band.d0)).max(0) as usize;
        let k1 = ((m - ji - band.d0).min(w as isize)).max(0) as usize;
        for k in 0..k0.min(w) {
            hc[k] = 0;
            ec[k] = NEG;
        }
        let row = j * w;
        let (mut col_best, mut col_row) = (0i32, 0usize);
        for k in k0..k1 {
            let i = (ji + band.d0 + k as isize) as usize;
            let e_ext = ep[k + 1] - GAP_EXT;
            let e_opn = hp[k + 1] - GAP_FIRST;
            let e = e_ext.max(e_opn);
            let e_open = e_opn >= e_ext;
            let mut h = hp[k] + prof.get(c, i);
            if e > h {
                h = e;
            }
            if f > h {
                h = f;
            }
            if h < 0 {
                h = 0;
            }
            let mut t = 0u8;
            if h == f {
                t |= VGAP;
            }
            if h == e {
                t |= HGAP;
            }
            if f_open {
                t |= VOPEN;
            }
            if e_open {
                t |= HOPEN;
            }
            tr[row + k] = t;
            if h >= col_best && h > 0 {
                col_best = h;
                col_row = i;
            }
            hc[k] = h;
            ec[k] = e;
            let f_ext_v = f - GAP_EXT;
            let f_opn_v = h - GAP_FIRST;
            f_open = f_opn_v >= f_ext_v;
            f = f_ext_v.max(f_opn_v);
        }
        if col_best > best {
            best = col_best;
            bi = col_row;
            bj = j;
        }
        for k in k1.max(k0)..w {
            hc[k] = 0;
            ec[k] = NEG;
        }
        std::mem::swap(&mut hp, &mut hc);
        std::mem::swap(&mut ep, &mut ec);
    }
    if best == 0 {
        return Alignment::default();
    }
    let at = |i: isize, j: isize| -> usize { j as usize * w + (i - j - band.d0) as usize };
    let mut al = Alignment { score: best, qend: bi + 1, send: bj + 1, ..Default::default() };
    let (mut i, mut j) = (bi as isize, bj as isize);
    let mut acc = 0i32;
    while i >= 0 && j >= 0 && acc < best {
        let t = tr[at(i, j)];
        if t & (VGAP | HGAP) == 0 {
            if q[i as usize] == s[j as usize] {
                al.identities += 1;
            } else {
                al.mismatches += 1;
            }
            acc += prof.get(s[j as usize], i as usize);
            al.length += 1;
            al.qstart = i as usize + 1;
            al.sstart = j as usize + 1;
            i -= 1;
            j -= 1;
        } else if t & VGAP != 0 {
            let mut l = 0i32;
            loop {
                let opened = tr[at(i, j)] & VOPEN != 0;
                l += 1;
                i -= 1;
                if opened || i <= 0 {
                    break;
                }
            }
            al.gap_opens += 1;
            al.gaps += l as usize;
            al.length += l as usize;
            acc -= GAP_FIRST - GAP_EXT + l * GAP_EXT;
        } else {
            let mut l = 0i32;
            loop {
                let opened = tr[at(i, j)] & HOPEN != 0;
                l += 1;
                j -= 1;
                if opened || j <= 0 {
                    break;
                }
            }
            al.gap_opens += 1;
            al.gaps += l as usize;
            al.length += l as usize;
            acc -= GAP_FIRST - GAP_EXT + l * GAP_EXT;
        }
    }
    debug_assert_eq!(acc, best, "traceback score mismatch");
    al
}

/// Full (unbanded) local alignment score — reference implementation.
pub fn sw_score(prof: &Profile, s: &[u8]) -> i32 {
    banded_score(prof, s, Band::full(prof.m, s.len()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::protsearch::alphabet::encode;

    fn rnd(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                ((x >> 33) % 20) as u8
            })
            .collect()
    }

    /// Plain O(mn) Gotoh reference (no band, no bias).
    fn gotoh(a: &[u8], b: &[u8]) -> i32 {
        let (m, n) = (a.len(), b.len());
        let mut h = vec![vec![0i32; n + 1]; m + 1];
        let mut e = vec![vec![NEG; n + 1]; m + 1];
        let mut f = vec![vec![NEG; n + 1]; m + 1];
        let mut best = 0;
        for i in 1..=m {
            for j in 1..=n {
                e[i][j] = (e[i][j - 1] - GAP_EXT).max(h[i][j - 1] - GAP_FIRST);
                f[i][j] = (f[i - 1][j] - GAP_EXT).max(h[i - 1][j] - GAP_FIRST);
                h[i][j] = (h[i - 1][j - 1] + score(a[i - 1], b[j - 1])).max(e[i][j]).max(f[i][j]).max(0);
                best = best.max(h[i][j]);
            }
        }
        best
    }

    #[test]
    fn full_band_matches_gotoh_and_traceback_is_consistent() {
        for t in 0..60u64 {
            let a = rnd(30 + (t as usize * 7) % 150, t + 1);
            let mut b = a.clone();
            for k in (0..b.len()).step_by(3 + (t as usize % 4)) {
                b[k] = ((b[k] as u64 + t) % 20) as u8;
            }
            if b.len() > 30 {
                b.drain(10..13 + (t as usize % 5));
            }
            if t % 3 == 0 && b.len() > 40 {
                let ins = rnd(4, t + 99);
                b.splice(25..25, ins);
            }
            let prof = Profile::new(&a, None);
            let full = Band::full(a.len(), b.len());
            let s1 = banded_score(&prof, &b, full);
            let al = banded_align(&a, &prof, &b, full);
            let g = gotoh(&a, &b);
            assert_eq!(s1, g, "case {t}");
            assert_eq!(al.score, g, "case {t}");
            assert_eq!(al.length, al.identities + al.mismatches + al.gaps);
            // recompute the score of the traced path from its statistics is not
            // possible without the path, but coordinates must be consistent
            assert!(al.qstart >= 1 && al.qstart <= al.qend && al.qend <= a.len());
            assert!(al.sstart >= 1 && al.sstart <= al.send && al.send <= b.len());
        }
    }

    #[test]
    fn narrow_band_is_lower_bound() {
        let a = rnd(200, 5);
        let mut b = a.clone();
        b.drain(100..130); // 30-residue deletion shifts the diagonal by 30
        let prof = Profile::new(&a, None);
        let full = banded_score(&prof, &b, Band::full(a.len(), b.len()));
        let narrow = banded_score(&prof, &b, Band { d0: -5, d1: 6 });
        assert!(narrow < full);
        let wide = banded_score(&prof, &b, Band { d0: -5, d1: 40 });
        assert_eq!(wide, full);
    }

    #[test]
    fn identical_and_gapped() {
        let a = encode(b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK");
        let prof = Profile::new(&a, None);
        let al = banded_align(&a, &prof, &a, Band::full(a.len(), a.len()));
        assert_eq!((al.qstart, al.qend, al.sstart, al.send), (1, a.len(), 1, a.len()));
        assert_eq!(al.identities, a.len());
        let mut b = a.clone();
        b.drain(20..24);
        let al = banded_align(&a, &prof, &b, Band::full(a.len(), b.len()));
        assert_eq!(al.gap_opens, 1);
        assert_eq!(al.gaps, 4);
        assert_eq!(al.identities, b.len());
    }

    #[test]
    fn xdrop_and_window() {
        let a = encode(b"MKTAYIAKQRQISFVKSHFSRQ");
        let zero = vec![0i8; a.len()];
        let seg = xdrop_bias(&a, &zero, &a, 5, 5);
        assert_eq!((seg.i, seg.j, seg.len), (0, 0, a.len()));
        let self_score: i32 = a.iter().map(|&c| score(c, c)).sum();
        assert_eq!(seg.score, self_score);
        let none = |_: isize| DELIM;
        assert_eq!(window_ungapped(&a, 5, &a, 5, none), self_score);
        assert!(fingerprint_ok(&a, 3, &a, 3, none, none));
        // far from the ends a mismatching window fails; near both starts the
        // aligned perimeter delimiters count as matches, as in DIAMOND
        let long = encode(b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ");
        let r = encode(b"WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
        assert!(!fingerprint_ok(&long, 20, &r, 20, none, none));
        assert!(fingerprint_ok(&long, 2, &r, 2, none, none));
    }

    #[test]
    fn sequence_set_layout() {
        let set = [encode(b"MKT"), encode(b"AY"), encode(b"GGW")];
        let get = |k: usize| set[k].as_slice();
        let at = |k: usize, p: isize| set_letter(get, set.len(), k, p);
        assert_eq!(at(1, 0), encode(b"A")[0]);
        assert_eq!(at(1, -1), DELIM);
        assert_eq!(at(1, -2), encode(b"T")[0]); // last letter of the previous sequence
        assert_eq!(at(1, -4), encode(b"M")[0]);
        assert_eq!(at(1, -5), DELIM); // perimeter before the first sequence
        assert_eq!(at(1, 2), DELIM);
        assert_eq!(at(1, 3), encode(b"G")[0]);
        assert_eq!(at(2, 3), DELIM);
        assert_eq!(at(2, 9), DELIM);
    }
}
