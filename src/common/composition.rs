//! Sequence-composition statistics shared across validators:
//! GC content and GC skew, tetranucleotide frequency vectors (the basis of the
//! CHIMERA-SCAN coherence check and a CheckM2 metadata feature), and k-mer
//! spectra. All are reference-free and database-free.

/// Map a nucleotide to a 2-bit code, or `None` for ambiguous bases.
#[inline]
pub fn base2(b: u8) -> Option<u8> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// GC fraction over A/C/G/T (ambiguous bases excluded from the denominator).
pub fn gc_content(seq: &[u8]) -> f64 {
    let mut gc = 0u64;
    let mut at = 0u64;
    for &b in seq {
        match b {
            b'G' | b'g' | b'C' | b'c' => gc += 1,
            b'A' | b'a' | b'T' | b't' => at += 1,
            _ => {}
        }
    }
    let tot = gc + at;
    if tot == 0 {
        0.0
    } else {
        gc as f64 / tot as f64
    }
}

/// Number of canonical tetranucleotides (256 folded with their reverse
/// complements: 120 pairs + 16 palindromes).
pub const TETRA_DIM: usize = 136;

/// Index of each 4-mer's canonical class (`min(kmer, revcomp(kmer))`).
fn canonical_tetra_index() -> [u8; 256] {
    let mut map = [0u8; 256];
    let mut next = 0u8;
    let mut seen = [u8::MAX; 256];
    for k in 0..256usize {
        let mut rc = 0usize;
        for i in 0..4 {
            let b = (k >> (2 * i)) & 3;
            rc = (rc << 2) | (3 - b);
        }
        let c = k.min(rc);
        if seen[c] == u8::MAX {
            seen[c] = next;
            next += 1;
        }
        map[k] = seen[c];
    }
    map
}

/// Canonical tetranucleotide frequency vector (136-dim, L1-normalised), so
/// the fingerprint does not depend on which strand genes lie on. 4-mers
/// containing an ambiguous base are skipped.
pub fn tetranucleotide_freqs(seq: &[u8]) -> [f64; TETRA_DIM] {
    let map = canonical_tetra_index();
    let mut counts = [0u64; TETRA_DIM];
    let mut total = 0u64;
    let mut idx = 0usize;
    let mut valid = 0usize;
    for &b in seq {
        match base2(b) {
            Some(code) => {
                idx = ((idx << 2) | code as usize) & 0xFF;
                valid += 1;
                if valid >= 4 {
                    counts[map[idx] as usize] += 1;
                    total += 1;
                }
            }
            None => {
                valid = 0;
                idx = 0;
            }
        }
    }
    let mut out = [0f64; TETRA_DIM];
    if total > 0 {
        let t = total as f64;
        for i in 0..TETRA_DIM {
            out[i] = counts[i] as f64 / t;
        }
    }
    out
}

/// Cosine similarity between two equal-length frequency vectors, in `[-1, 1]`.
pub fn cosine(a: &[f64], b: &[f64]) -> f64 {
    let mut dot = 0.0;
    let mut na = 0.0;
    let mut nb = 0.0;
    for (x, y) in a.iter().zip(b.iter()) {
        dot += x * y;
        na += x * x;
        nb += y * y;
    }
    if na == 0.0 || nb == 0.0 {
        0.0
    } else {
        dot / (na.sqrt() * nb.sqrt())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gc_is_correct() {
        assert!((gc_content(b"GGCC") - 1.0).abs() < 1e-9);
        assert!((gc_content(b"ATAT") - 0.0).abs() < 1e-9);
        assert!((gc_content(b"ACGT") - 0.5).abs() < 1e-9);
        assert!((gc_content(b"ACGTNNNN") - 0.5).abs() < 1e-9); // N excluded
    }

    #[test]
    fn tetra_sums_to_one() {
        let f = tetranucleotide_freqs(b"ACGTACGTACGTACGT");
        let s: f64 = f.iter().sum();
        assert!((s - 1.0).abs() < 1e-9);
    }

    #[test]
    fn tetra_is_strand_independent() {
        let s = b"ACGGTCAGTTTGACCAGTAGGCATCGATCGGGAT";
        let rc: Vec<u8> = s
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect();
        let (a, b) = (tetranucleotide_freqs(s), tetranucleotide_freqs(&rc));
        for i in 0..TETRA_DIM {
            assert!((a[i] - b[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn cosine_self_is_one() {
        let f = tetranucleotide_freqs(b"ACGTACGTACGTTTAA");
        assert!((cosine(&f, &f) - 1.0).abs() < 1e-9);
    }
}
