//! Spaced seeds over DIAMOND's reduced alphabet and the query seed index.
//!
//! DIAMOND 2.1's default sensitivity uses two weight-10 shapes,
//! `111101110111` and `111011010010111` (from `diamond blastp --log`), over
//! the 10-letter alphabet in [`super::alphabet::REDUCED`]. The query set is
//! small (one genome's proteins), so it is indexed in memory and the
//! database is streamed past it.

use super::alphabet::{REDUCED, RED_INVALID, TRUE_AA};

/// DIAMOND 2.1 default-mode seed shapes.
pub const DEFAULT_SHAPES: [&str; 2] = ["111101110111", "111011010010111"];

/// A compiled spaced-seed shape (span ≤ 16).
#[derive(Clone, Debug)]
pub struct Shape {
    pub span: usize,
    pub weight: usize,
    /// (shift, nibble-count) of each run of consecutive `1`s inside the
    /// 4-bit-per-letter rolling window (the newest letter is nibble 0).
    runs: Vec<(u32, u32)>,
    /// Bitmask over the window (bit j = letter j positions back) of the
    /// positions the shape reads; used for validity.
    used: u32,
}

impl Shape {
    pub fn parse(code: &str) -> Shape {
        let b = code.as_bytes();
        let span = b.len();
        assert!((1..=16).contains(&span), "seed span must be 1..=16");
        let mut runs = Vec::new();
        let mut used = 0u32;
        let mut o = 0usize;
        let mut weight = 0;
        while o < span {
            if b[o] == b'1' {
                let start = o;
                while o < span && b[o] == b'1' {
                    used |= 1 << (span - 1 - o);
                    o += 1;
                    weight += 1;
                }
                // offsets start..o-1 occupy nibbles (span-o) ..= (span-1-start)
                let low_nibble = (span - o) as u32;
                runs.push((4 * low_nibble, (o - start) as u32));
            } else {
                o += 1;
            }
        }
        Shape { span, weight, runs, used }
    }

    /// Key of the window `win` (4 bits per letter, newest in the low nibble).
    #[inline(always)]
    pub fn key(&self, win: u64) -> u64 {
        let mut k = 0u64;
        for &(shift, n) in &self.runs {
            let mask = (1u64 << (4 * n)) - 1;
            k = (k << (4 * n)) | ((win >> shift) & mask);
        }
        k
    }

    /// True when every letter the shape reads is valid (`inv` holds one bit
    /// per window position, set for invalid letters).
    #[inline(always)]
    pub fn valid(&self, inv: u32) -> bool {
        inv & self.used == 0
    }
}

#[inline(always)]
fn mix(key: u64) -> u64 {
    // splitmix64 finaliser: good avalanche for the packed 40-bit keys.
    let mut z = key.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Reduced letter of a code on the target side: like DIAMOND, B/J/Z (codes
/// 20–22) reduce to group 0, X and `*` break seeds.
#[inline(always)]
fn reduce_target(c: u8) -> u8 {
    let c = (c as usize).min(24);
    if c < TRUE_AA {
        REDUCED[c]
    } else if c <= 22 {
        0
    } else {
        RED_INVALID
    }
}

/// Iterate the seeds of one encoded sequence for `shape`, calling
/// `f(position, key)` for every valid window (position = first letter).
/// `query` selects the query-side rules (only the 20 standard residues may
/// seed; DIAMOND drops query seeds containing B/J/Z as non-complex).
#[inline]
pub fn for_each_seed<F: FnMut(usize, u64)>(seq: &[u8], shape: &Shape, query: bool, mut f: F) {
    if seq.len() < shape.span {
        return;
    }
    let mut win = 0u64;
    let mut inv = 0u32;
    let span = shape.span;
    let keep = if span >= 32 { u32::MAX } else { (1u32 << span) - 1 };
    for (e, &c) in seq.iter().enumerate() {
        let r = if query { REDUCED[(c as usize).min(24)] } else { reduce_target(c) };
        win = (win << 4) | (if r == RED_INVALID { 0 } else { r as u64 });
        inv = ((inv << 1) | (r == RED_INVALID) as u32) & keep;
        if e + 1 >= span && shape.valid(inv) {
            f(e + 1 - span, shape.key(win));
        }
    }
}

/// DIAMOND's seed complexity test (`search/seed_complexity.cpp`): the
/// log multinomial coefficient of the reduced-letter counts must reach
/// `seed_cut · ln 2 · weight` with seed_cut = 0.8 (default sensitivity).
#[inline]
fn complex_enough(key: u64, weight: usize) -> bool {
    let mut counts = [0usize; 16];
    let mut k = key;
    for _ in 0..weight {
        counts[(k & 0xF) as usize] += 1;
        k >>= 4;
    }
    let lnfact = |n: usize| -> f64 { (2..=n).map(|x| (x as f64).ln()).sum() };
    let mut entropy = lnfact(weight);
    for &c in &counts {
        entropy -= lnfact(c);
    }
    entropy >= 0.8 * std::f64::consts::LN_2 * weight as f64
}

/// In-memory index of the query seeds for one shape.
pub struct ShapeIndex {
    pub shape: Shape,
    filter_bits: u32,
    filter: Vec<u64>,
    table_bits: u32,
    /// Open-addressing slots holding `unique index + 1` (0 = empty).
    slots: Vec<u32>,
    keys: Vec<u64>,
    offs: Vec<u32>,
    /// `(query, position)` grouped by key.
    entries: Vec<(u32, u32)>,
}

impl ShapeIndex {
    pub fn build(shape: Shape, queries: &[Vec<u8>]) -> ShapeIndex {
        let mut tmp: Vec<(u64, u32, u32)> = Vec::new();
        for (qi, q) in queries.iter().enumerate() {
            for_each_seed(q, &shape, true, |p, k| {
                if complex_enough(k, shape.weight) {
                    tmp.push((k, qi as u32, p as u32));
                }
            });
        }
        tmp.sort_unstable();
        let mut keys = Vec::new();
        let mut offs = Vec::new();
        let mut entries = Vec::with_capacity(tmp.len());
        for &(k, q, p) in &tmp {
            if keys.last() != Some(&k) {
                keys.push(k);
                offs.push(entries.len() as u32);
            }
            entries.push((q, p));
        }
        offs.push(entries.len() as u32);
        drop(tmp);

        let n = keys.len().max(1);
        // presence filter: ~32 bits per key, 2^16 ..= 2^31 bits
        let filter_bits = (64 - ((n as u64) * 32).leading_zeros()).clamp(16, 31);
        let mut filter = vec![0u64; (1usize << filter_bits) / 64];
        let table_bits = (64 - ((n as u64) * 2).leading_zeros()).clamp(4, 31);
        let tmask = (1usize << table_bits) - 1;
        let mut slots = vec![0u32; 1usize << table_bits];
        for (ui, &k) in keys.iter().enumerate() {
            let h = mix(k);
            let fb = (h >> (64 - filter_bits)) as usize;
            filter[fb >> 6] |= 1u64 << (fb & 63);
            let mut s = (h as usize) & tmask;
            while slots[s] != 0 {
                s = (s + 1) & tmask;
            }
            slots[s] = ui as u32 + 1;
        }
        ShapeIndex { shape, filter_bits, filter, table_bits, slots, keys, offs, entries }
    }

    /// Query occurrences of `key`, or an empty slice.
    #[inline(always)]
    pub fn lookup(&self, key: u64) -> &[(u32, u32)] {
        let h = mix(key);
        let fb = (h >> (64 - self.filter_bits)) as usize;
        if self.filter[fb >> 6] & (1u64 << (fb & 63)) == 0 {
            return &[];
        }
        let tmask = (1usize << self.table_bits) - 1;
        let mut s = (h as usize) & tmask;
        loop {
            let v = self.slots[s];
            if v == 0 {
                return &[];
            }
            let ui = (v - 1) as usize;
            if self.keys[ui] == key {
                return &self.entries[self.offs[ui] as usize..self.offs[ui + 1] as usize];
            }
            s = (s + 1) & tmask;
        }
    }

    pub fn n_seeds(&self) -> usize {
        self.entries.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::protsearch::alphabet::encode;

    #[test]
    fn shape_parsing() {
        let s = Shape::parse("111101110111");
        assert_eq!(s.span, 12);
        assert_eq!(s.weight, 10);
        let s2 = Shape::parse("111011010010111");
        assert_eq!(s2.span, 15);
        assert_eq!(s2.weight, 10);
    }

    #[test]
    fn same_reduced_letters_same_key() {
        // I/L/V/M share a reduced letter; K/R one; E/D/N/Q one; S/T one.
        let a = encode(b"ILVMKREDQNSTAGH");
        let b = encode(b"LVMIRKDEQNTSAGH");
        let s = Shape::parse("111101110111");
        let mut ka = Vec::new();
        let mut kb = Vec::new();
        for_each_seed(&a, &s, true, |p, k| ka.push((p, k)));
        for_each_seed(&b, &s, true, |p, k| kb.push((p, k)));
        assert_eq!(ka, kb);
        assert_eq!(ka.len(), a.len() - 12 + 1);
    }

    #[test]
    fn spaced_positions_are_ignored() {
        // Position 4 (0-based) is a '0' in the first shape: changing it to
        // another valid letter must not change the key.
        let s = Shape::parse("111101110111");
        let a = encode(b"ACDEFGHIKLMN");
        let mut b = a.clone();
        b[4] = encode(b"W")[0];
        let mut ka = 0;
        let mut kb = 0;
        for_each_seed(&a, &s, true, |_, k| ka = k);
        for_each_seed(&b, &s, true, |_, k| kb = k);
        assert_eq!(ka, kb);
        // ... but changing a '1' position does.
        b[0] = encode(b"W")[0];
        for_each_seed(&b, &s, true, |_, k| kb = k);
        assert_ne!(ka, kb);
    }

    #[test]
    fn invalid_letters_break_seeds() {
        let s = Shape::parse("111101110111");
        let a = encode(b"ACDEFGHIKLMNXACDEFGHIKLMN");
        let mut n = 0;
        for_each_seed(&a, &s, true, |_, _| n += 1);
        // windows containing the X at index 12 in a '1' position are skipped
        assert!(n < a.len() - 12 + 1);
    }

    #[test]
    fn complexity_cut_matches_diamond() {
        let key = |counts: &[usize]| -> u64 {
            let mut k = 0u64;
            for (letter, &c) in counts.iter().enumerate() {
                for _ in 0..c {
                    k = (k << 4) | letter as u64;
                }
            }
            k
        };
        assert!(!complex_enough(key(&[10]), 10));
        assert!(!complex_enough(key(&[5, 5]), 10)); // ln 252 = 5.53 < 5.545
        assert!(!complex_enough(key(&[6, 4]), 10));
        assert!(complex_enough(key(&[7, 2, 1]), 10)); // ln 360 = 5.89
        assert!(complex_enough(key(&[1; 10]), 10));
    }

    #[test]
    fn index_lookup() {
        let qs = [encode(b"MKTAYIAKQRQISFVKSHFSRQ"), encode(b"GGGGGGGGGGGGGGGGGGGG")];
        let idx = ShapeIndex::build(Shape::parse("111101110111"), &qs);
        // poly-G is low complexity -> not indexed (DIAMOND's entropy cut)
        assert!(idx.entries.iter().all(|&(q, _)| q == 0));
        let mut found = 0;
        for_each_seed(&qs[0], &idx.shape, true, |p, k| {
            let hits = idx.lookup(k);
            assert!(hits.iter().any(|&(q, pp)| q == 0 && pp as usize == p));
            found += 1;
        });
        assert!(found > 0);
        let _ = idx.lookup(0x01_2345_6789);
    }
}
