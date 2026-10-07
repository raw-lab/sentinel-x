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

/// Cumulative GC skew `(G - C) / (G + C)` over a sliding window. Returns one
/// value per non-overlapping `window`-sized bin. The cumulative form is used by
/// CHIMERA-SCAN to locate replication origins / assembly artifacts.
pub fn gc_skew(seq: &[u8], window: usize) -> Vec<f64> {
    let window = window.max(1);
    let mut out = Vec::new();
    let mut g = 0i64;
    let mut c = 0i64;
    let mut i = 0usize;
    for &b in seq {
        match b {
            b'G' | b'g' => g += 1,
            b'C' | b'c' => c += 1,
            _ => {}
        }
        i += 1;
        if i % window == 0 {
            let denom = (g + c) as f64;
            out.push(if denom > 0.0 { (g - c) as f64 / denom } else { 0.0 });
            g = 0;
            c = 0;
        }
    }
    if i % window != 0 {
        let denom = (g + c) as f64;
        out.push(if denom > 0.0 { (g - c) as f64 / denom } else { 0.0 });
    }
    out
}

/// Cumulative GC-skew curve (running sum of per-base ±1), used to find the
/// global minimum/maximum (origin/terminus signatures).
pub fn cumulative_gc_skew(seq: &[u8]) -> Vec<i64> {
    let mut acc = 0i64;
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq {
        match b {
            b'G' | b'g' => acc += 1,
            b'C' | b'c' => acc -= 1,
            _ => {}
        }
        out.push(acc);
    }
    out
}

/// Number of distinct canonical tetranucleotides (256 raw; 136 after reverse-
/// complement folding). We use the full 256-dim vector for simplicity.
pub const TETRA_DIM: usize = 256;

/// Tetranucleotide frequency vector (256-dim, L1-normalised). Windows
/// containing an ambiguous base are skipped. This is the composition fingerprint
/// CHIMERA-SCAN compares between sub-regions of a contig.
pub fn tetranucleotide_freqs(seq: &[u8]) -> [f64; TETRA_DIM] {
    let mut counts = [0u64; TETRA_DIM];
    let mut total = 0u64;
    if seq.len() >= 4 {
        let mut k = 0usize;
        let mut idx = 0usize;
        let mut valid = 0usize;
        for &b in seq {
            match base2(b) {
                Some(code) => {
                    idx = ((idx << 2) | code as usize) & 0xFF;
                    valid += 1;
                    if valid >= 4 {
                        counts[idx] += 1;
                        total += 1;
                    }
                }
                None => {
                    valid = 0;
                    idx = 0;
                }
            }
            k += 1;
            let _ = k;
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

/// Euclidean distance between two equal-length vectors.
pub fn euclidean(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        .sqrt()
}

/// Count canonical k-mers (k ≤ 12) into a hashable spectrum. Returns
/// `(distinct_kmers, total_kmers)`; used for redundancy / duplication signals
/// (e.g. CheckV's `kmer_freq` warning for large duplications).
pub fn kmer_freq(seq: &[u8], k: usize) -> f64 {
    use std::collections::HashMap;
    if k == 0 || seq.len() < k {
        return 1.0;
    }
    let mut counts: HashMap<u64, u32> = HashMap::new();
    let mask = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let mut idx = 0u64;
    let mut valid = 0usize;
    let mut total = 0u64;
    for &b in seq {
        match base2(b) {
            Some(code) => {
                idx = ((idx << 2) | code as u64) & mask;
                valid += 1;
                if valid >= k {
                    *counts.entry(idx).or_insert(0) += 1;
                    total += 1;
                }
            }
            None => {
                valid = 0;
                idx = 0;
            }
        }
    }
    if counts.is_empty() {
        1.0
    } else {
        // mean multiplicity per distinct k-mer; >1.5 hints at large duplication
        total as f64 / counts.len() as f64
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
    fn cosine_self_is_one() {
        let f = tetranucleotide_freqs(b"ACGTACGTACGTTTAA");
        assert!((cosine(&f, &f) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn duplication_raises_kmer_freq() {
        let unique = b"ACGTACAGTCAGTCAGGGCTACGATCGATCGTAGCTAGCATCGAT";
        let dup: Vec<u8> = unique.iter().chain(unique.iter()).copied().collect();
        assert!(kmer_freq(&dup, 12) > kmer_freq(unique, 12));
    }
}
