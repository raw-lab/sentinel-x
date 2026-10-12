//! Terminal repeats and genome copy number (CheckV `complete_genomes.py`,
//! `utility.get_average_kmer_freq`).

use std::collections::HashMap;

/// Kind of terminal repeat.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TrKind {
    Dtr,
    Itr,
}

impl TrKind {
    pub fn label(&self) -> &'static str {
        match self {
            TrKind::Dtr => "DTR",
            TrKind::Itr => "ITR",
        }
    }
}

/// A terminal repeat with CheckV's filtering statistics.
#[derive(Clone, Debug)]
pub struct TerminalRepeat {
    pub kind: TrKind,
    pub seq: Vec<u8>,
    /// Frequency of the most common base in the repeat.
    pub mode_freq: f64,
    /// Fraction of `N` in the repeat.
    pub n_freq: f64,
    /// Non-overlapping occurrences of the repeat in the contig.
    pub count: usize,
}

impl TerminalRepeat {
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        self.seq.len()
    }
}

/// `fetch_dtr`: the longest suffix starting at an occurrence of the first
/// `min_len` bases in the second half of the contig (`pos ≥ L/2`, float
/// division) that equals the contig prefix of the same length. Matches are
/// overlapping and tried left to right.
pub fn fetch_dtr(seq: &[u8], min_len: usize) -> &[u8] {
    let l = seq.len();
    let k = min_len.min(l);
    let start = &seq[..k];
    if k == 0 {
        return &[];
    }
    let half = l as f64 / 2.0;
    let first = half.ceil() as usize;
    let mut pos = first;
    while pos + k <= l {
        if &seq[pos..pos + k] == start {
            let end = &seq[pos..];
            if &seq[..end.len()] == end {
                return end;
            }
        }
        pos += 1;
    }
    &[]
}

/// Reverse complement translating only `ACTG` (CheckV's `str.maketrans`).
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'T' => b'A',
            b'G' => b'C',
            o => o,
        })
        .collect()
}

/// `fetch_itr`: the prefix shared with the reverse complement, if at least
/// `min_len` long, extended up to `max_len`.
pub fn fetch_itr(seq: &[u8], min_len: usize, max_len: usize) -> &[u8] {
    let rev = reverse_complement(seq);
    let pre = |i: usize| -> (&[u8], &[u8]) { (&seq[..i.min(seq.len())], &rev[..i.min(rev.len())]) };
    let (a, b) = pre(min_len);
    if a != b {
        return &[];
    }
    let mut i = min_len + 1;
    loop {
        let (a, b) = pre(i);
        if a == b && i <= max_len {
            i += 1;
        } else {
            break;
        }
    }
    &seq[..(i - 1).min(seq.len())]
}

/// Detect and characterise a terminal repeat (`complete_genomes` steps 3–4).
pub fn detect(seq: &[u8], min_len: usize) -> Option<TerminalRepeat> {
    let dtr = fetch_dtr(seq, min_len);
    let itr = fetch_itr(seq, min_len, 1000);
    let (kind, tr) = if dtr.len() < min_len && itr.len() < min_len {
        return None;
    } else if dtr.len() >= itr.len() {
        (TrKind::Dtr, dtr)
    } else {
        (TrKind::Itr, itr)
    };
    let mut counts: HashMap<u8, usize> = HashMap::new();
    for &b in tr {
        *counts.entry(b).or_insert(0) += 1;
    }
    let mode = counts.values().copied().max().unwrap_or(0);
    let n = tr.iter().filter(|&&b| b == b'N').count();
    Some(TerminalRepeat {
        kind,
        seq: tr.to_vec(),
        mode_freq: mode as f64 / tr.len() as f64,
        n_freq: n as f64 / tr.len() as f64,
        count: count_non_overlapping(seq, tr),
    })
}

/// Python `str.count`: non-overlapping occurrences, left to right.
pub fn count_non_overlapping(hay: &[u8], needle: &[u8]) -> usize {
    if needle.is_empty() || needle.len() > hay.len() {
        return 0;
    }
    let mut n = 0;
    let mut i = 0;
    while i + needle.len() <= hay.len() {
        if &hay[i..i + needle.len()] == needle {
            n += 1;
            i += needle.len();
        } else {
            i += 1;
        }
    }
    n
}

/// CheckV's TR filters, in order (first match wins).
pub fn flag_reason(tr: &TerminalRepeat, max_ambig: f64, max_count: usize, max_basefreq: f64) -> Option<&'static str> {
    if tr.n_freq > max_ambig {
        Some("Too many ambiguous bases in TR")
    } else if tr.count > max_count {
        Some("Repetetive TR sequence")
    } else if tr.mode_freq > max_basefreq {
        Some("Low complexity TR")
    } else {
        None
    }
}

/// Mean count of the distinct canonical 21-mers (kcounter semantics:
/// k-mers with any non-ACGT base are skipped); NaN when there are none.
pub fn average_kmer_freq(seq: &[u8], k: usize) -> f64 {
    let enc = |b: u8| -> Option<u64> {
        match b {
            b'A' | b'a' => Some(0),
            b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(2),
            b'T' | b't' => Some(3),
            _ => None,
        }
    };
    if k == 0 || k > 32 || seq.len() < k {
        return f64::NAN;
    }
    let mask = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let shift = 2 * (k as u64 - 1);
    let mut counts: HashMap<u64, u32> = HashMap::new();
    let (mut fwd, mut rev, mut valid) = (0u64, 0u64, 0usize);
    for &b in seq {
        match enc(b) {
            Some(c) => {
                fwd = ((fwd << 2) | c) & mask;
                rev = (rev >> 2) | ((3 - c) << shift);
                valid += 1;
                if valid >= k {
                    *counts.entry(fwd.min(rev)).or_insert(0) += 1;
                }
            }
            None => {
                valid = 0;
                fwd = 0;
                rev = 0;
            }
        }
    }
    if counts.is_empty() {
        return f64::NAN;
    }
    // numpy.mean: pairwise summation over the counts (exact for integers)
    let total: u64 = counts.values().map(|&c| c as u64).sum();
    total as f64 / counts.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dtr_and_itr() {
        let rep = b"ACGTACGTACGTACGTACGTACGTAAGGCC";
        let mut seq = Vec::new();
        seq.extend_from_slice(rep);
        seq.extend_from_slice(b"TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGG");
        seq.extend_from_slice(rep);
        let tr = detect(&seq, 21).unwrap();
        assert_eq!(tr.kind, TrKind::Dtr);
        assert_eq!(tr.len(), 30);
        assert_eq!(tr.count, 2);

        let arm = b"ACGTACGTACGTACGTACGTAC";
        let mut s2 = Vec::new();
        s2.extend_from_slice(arm);
        s2.extend_from_slice(b"TTTTAGAGAGAGATTTTCCCCGGGGAAAA");
        s2.extend_from_slice(&reverse_complement(arm));
        assert_eq!(detect(&s2, 21).unwrap().kind, TrKind::Itr);
    }

    #[test]
    fn odd_length_dtr_uses_float_half() {
        // 101 bp periodic contig: CheckV reports a 26 bp DTR (matches start
        // at pos >= 50.5, i.e. 51 onwards)
        let unit = b"ACGTTGCAACGTTGCAACGATCGATCGATCG"; // arbitrary
        let mut s: Vec<u8> = Vec::new();
        while s.len() < 101 {
            s.extend_from_slice(&unit[..25]);
        }
        s.truncate(101);
        let d = fetch_dtr(&s, 21);
        assert_eq!(d.len(), 26);
    }

    #[test]
    fn itr_prefers_longer_and_short_palindromes() {
        // 30 bp perfect palindrome: whole sequence is an ITR
        let half = b"ACGTTGCAAGCTTAG";
        let mut s = half.to_vec();
        s.extend_from_slice(&reverse_complement(half));
        let t = detect(&s, 21).unwrap();
        assert_eq!(t.kind, TrKind::Itr);
        assert_eq!(t.len(), 30);
    }

    #[test]
    fn n_runs_are_flagged() {
        let mut s = vec![b'N'; 40];
        s.extend_from_slice(b"ACGTACGATCGATCGTAGCTAGCTAGCATCGATCGATGCAT");
        s.extend(vec![b'N'; 40]);
        let tr = detect(&s, 21).unwrap();
        assert_eq!(flag_reason(&tr, 0.2, 8, 0.70), Some("Too many ambiguous bases in TR"));
    }

    #[test]
    fn canonical_kmers() {
        // kcounter: k-mers with N skipped, reverse complements merged
        let f = average_kmer_freq(b"ACGTNACGTACGTAAAACGT", 4);
        // {TAAA:1, CGTA:3, AACG:1, AAAC:1, AAAA:1, GTAC:1, ACGT:4, GTAA:1} -> 13/8
        assert!((f - 13.0 / 8.0).abs() < 1e-12, "{f}");
        assert!(average_kmer_freq(b"NNNN", 21).is_nan());
        assert_eq!(count_non_overlapping(b"AAAA", b"AA"), 2);
    }
}
