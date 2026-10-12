//! ## TRACE — Terminal-Repeat & Architecture Completeness Estimator
//!
//! A reference-free completeness signal for viruses/phages that fuses three
//! database-independent lines of evidence:
//!
//! 1. **Terminal repeats** — a direct (DTR) or inverted (ITR) terminal
//!    repeat that passes CheckV's filters (≤ 20 % N, ≤ 8 copies, most common
//!    base ≤ 70 %) and no sign of multiple genome copies (mean canonical
//!    21-mer count ≤ 1.5);
//! 2. **Circular permutation** — a head/tail k-mer overlap beyond the
//!    terminal repeat;
//! 3. **Coding architecture** — complete phage genomes are gene-dense with
//!    few long non-coding gaps (coding density = union of gene spans /
//!    contig length, so overlapping genes are not double-counted).
//!
//! TRACE reports a confidence label; it is not a completeness estimate.

use crate::checkv::repeats::{average_kmer_freq, detect, flag_reason, TerminalRepeat};
use crate::common::pyfmt::np_round;
use crate::common::Gene;

/// Confidence label for a TRACE completeness call.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TraceConfidence {
    Complete,
    HighlyComplete,
    Partial,
    Fragment,
}

impl TraceConfidence {
    pub fn label(&self) -> &'static str {
        match self {
            TraceConfidence::Complete => "complete",
            TraceConfidence::HighlyComplete => "highly-complete",
            TraceConfidence::Partial => "partial",
            TraceConfidence::Fragment => "fragment",
        }
    }
}

/// TRACE output.
#[derive(Clone, Debug)]
pub struct TraceResult {
    pub terminal_repeat: Option<TerminalRepeat>,
    /// Why the repeat was rejected (CheckV's TR filters / copy number).
    pub repeat_flag: Option<&'static str>,
    /// Length of the detected head/tail k-mer overlap (circular signal), 0 if none.
    pub circular_overlap: usize,
    /// Fraction of the contig covered by predicted genes.
    pub coding_density: f64,
    /// Largest non-coding stretch as a fraction of contig length.
    pub max_gap_frac: f64,
    pub confidence: TraceConfidence,
}

/// Detect a circular-permutation overlap: the longest prefix (≥ `min_k`) that
/// re-occurs at the contig's tail (excluding an already-known terminal repeat).
fn circular_overlap(seq: &[u8], skip: usize, min_k: usize, max_k: usize) -> usize {
    let n = seq.len();
    if n < 2 * (min_k + skip) {
        return 0;
    }
    let mut best = 0;
    let upper = max_k.min((n - 2 * skip) / 2);
    for k in min_k..=upper {
        if seq[skip..skip + k] == seq[n - skip - k..n - skip] {
            best = k;
        }
    }
    best
}

/// Coding density (union of gene spans) and the largest uncovered stretch.
pub fn coding_architecture(len: usize, genes: &[Gene]) -> (f64, f64) {
    let n = len.max(1);
    let mut spans: Vec<(usize, usize)> = genes.iter().map(|g| (g.begin.min(g.end), g.begin.max(g.end).min(n))).collect();
    spans.sort_unstable();
    let (mut covered, mut last_end, mut max_gap) = (0usize, 0usize, 0usize);
    for (b, e) in spans {
        if b > last_end + 1 {
            max_gap = max_gap.max(b - last_end - 1);
        }
        if e > last_end {
            covered += e - (b.max(last_end + 1)) + 1;
            last_end = e;
        }
    }
    max_gap = max_gap.max(n.saturating_sub(last_end));
    (covered as f64 / n as f64, max_gap as f64 / n as f64)
}

/// Run TRACE on one contig's sequence and its predicted genes.
pub fn trace(seq: &[u8], genes: &[Gene]) -> TraceResult {
    let tr = detect(seq, 21);
    let mut flag = tr.as_ref().and_then(|t| flag_reason(t, 0.2, 8, 0.70));
    if tr.is_some() && flag.is_none() && np_round(average_kmer_freq(seq, 21), 2) > 1.5 {
        flag = Some("Multiple genome copies detected");
    }
    let skip = tr.as_ref().map(|t| t.len()).unwrap_or(0);
    let overlap = circular_overlap(seq, skip, 15, 500);
    let (coding_density, max_gap_frac) = coding_architecture(seq.len(), genes);

    let good_tr = tr.is_some() && flag.is_none();
    let confidence = if good_tr && coding_density > 0.6 {
        TraceConfidence::Complete
    } else if (overlap >= 15 && coding_density > 0.6) || (coding_density > 0.7 && max_gap_frac < 0.15) {
        TraceConfidence::HighlyComplete
    } else if coding_density > 0.4 {
        TraceConfidence::Partial
    } else {
        TraceConfidence::Fragment
    };
    TraceResult { terminal_repeat: tr, repeat_flag: flag, circular_overlap: overlap, coding_density, max_gap_frac, confidence }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::checkv::repeats::TrKind;

    fn gene(b: usize, e: usize) -> Gene {
        Gene {
            contig: "c".into(),
            contig_index: 0,
            index: 1,
            begin: b,
            end: e,
            strand: 1,
            aa: vec![b'M'; (e - b) / 3],
            partial_left: false,
            partial_right: false,
            stop: true,
        }
    }

    fn random_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x % 4) as usize]
            })
            .collect()
    }

    #[test]
    fn dtr_plus_dense_is_complete() {
        let repeat = b"ACGTTGCATGCAAGCTTAGCCGATAAGG";
        let mut seq = Vec::new();
        seq.extend_from_slice(repeat);
        seq.extend(random_seq(3000, 7));
        seq.extend_from_slice(repeat);
        let genes = vec![gene(1, (seq.len() as f64 * 0.8) as usize)];
        let r = trace(&seq, &genes);
        assert_eq!(r.terminal_repeat.as_ref().unwrap().kind, TrKind::Dtr);
        assert_eq!(r.confidence, TraceConfidence::Complete);
    }

    #[test]
    fn n_run_repeat_is_not_complete() {
        let mut seq = vec![b'N'; 60];
        seq.extend(random_seq(3000, 3));
        seq.extend(vec![b'N'; 60]);
        let genes = vec![gene(61, 2900)];
        let r = trace(&seq, &genes);
        assert_eq!(r.repeat_flag, Some("Too many ambiguous bases in TR"));
        assert_ne!(r.confidence, TraceConfidence::Complete);
    }

    #[test]
    fn overlapping_genes_are_not_double_counted() {
        let (d, gap) = coding_architecture(1000, &[gene(1, 600), gene(301, 900)]);
        assert!((d - 0.9).abs() < 1e-12, "{d}");
        assert!((gap - 0.1).abs() < 1e-12, "{gap}");
    }

    #[test]
    fn sparse_coding_is_fragment() {
        let seq = random_seq(5000, 11);
        let r = trace(&seq, &[gene(1, 300)]);
        assert_eq!(r.confidence, TraceConfidence::Fragment);
    }
}
