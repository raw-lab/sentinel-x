//! ## TRACE — Terminal-Repeat & Architecture Completeness Estimator
//!
//! A reference-free completeness signal for viruses/phages. It fuses three
//! database-independent lines of evidence:
//!
//! 1. **Terminal repeats** — a direct (DTR) or inverted (ITR) terminal repeat is
//!    strong evidence of a complete genome (re-using [`sentinel_checkv`]).
//! 2. **Circular permutation** — a k-mer overlap between the contig's head and
//!    tail beyond the terminal repeat suggests a circular / headful-packaged
//!    molecule, i.e. a complete genome.
//! 3. **Coding architecture** — complete phage genomes are gene-dense with few
//!    long non-coding gaps; very low coding density or a large terminal gap
//!    flags an incomplete or fragmented assembly.

use sentinel_checkv::{detect_terminal_repeat, TerminalRepeat, TrKind};
use sentinel_core::composition::base2;
use sentinel_core::Gene;

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
    pub terminal_repeat: TerminalRepeat,
    /// Length of the detected head/tail k-mer overlap (circular signal), 0 if none.
    pub circular_overlap: usize,
    /// Fraction of the contig covered by predicted genes.
    pub coding_density: f64,
    /// Largest inter-gene gap as a fraction of contig length.
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
    // Compare the prefix starting after `skip` to the suffix ending before `skip`.
    let pstart = skip;
    let mut best = 0;
    let upper = max_k.min((n - skip) / 2);
    let mut k = min_k;
    while k <= upper {
        let prefix = &seq[pstart..pstart + k];
        let suffix = &seq[n - skip - k..n - skip];
        if prefix == suffix {
            best = k;
        }
        k += 1;
    }
    best
}

/// Run TRACE on one contig's sequence and its predicted genes.
pub fn trace(seq: &[u8], genes: &[Gene]) -> TraceResult {
    let tr = detect_terminal_repeat(seq);

    // Circular overlap beyond any terminal repeat.
    let skip = tr.length;
    let overlap = circular_overlap(seq, skip, 15, 500);

    // Coding density & largest gap.
    let n = seq.len().max(1);
    let mut covered = 0usize;
    let mut spans: Vec<(usize, usize)> = genes
        .iter()
        .map(|g| (g.begin.min(g.end), g.begin.max(g.end)))
        .collect();
    spans.sort_unstable();
    let mut last_end = 0usize;
    let mut max_gap = 0usize;
    for (b, e) in &spans {
        let b = *b;
        let e = (*e).min(n);
        covered += e.saturating_sub(b);
        if b > last_end {
            max_gap = max_gap.max(b - last_end);
        }
        last_end = last_end.max(e);
    }
    max_gap = max_gap.max(n.saturating_sub(last_end));
    let coding_density = covered as f64 / n as f64;
    let max_gap_frac = max_gap as f64 / n as f64;

    // Decision.
    let has_tr = tr.kind != TrKind::None;
    let confidence = if has_tr && coding_density > 0.6 {
        TraceConfidence::Complete
    } else if overlap >= 15 && coding_density > 0.6 {
        TraceConfidence::HighlyComplete
    } else if coding_density > 0.7 && max_gap_frac < 0.15 {
        TraceConfidence::HighlyComplete
    } else if coding_density > 0.4 {
        TraceConfidence::Partial
    } else {
        TraceConfidence::Fragment
    };

    TraceResult {
        terminal_repeat: tr,
        circular_overlap: overlap,
        coding_density,
        max_gap_frac,
        confidence,
    }
}

/// GC-independent helper retained for callers that pre-screen non-ACGT contigs.
pub fn acgt_fraction(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let n = seq.iter().filter(|&&b| base2(b).is_some()).count();
    n as f64 / seq.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gene(b: usize, e: usize) -> Gene {
        Gene {
            contig: "c".into(),
            index: 1,
            begin: b,
            end: e,
            strand: 1,
            aa: vec![b'M'; (e - b) / 3],
            partial: false,
        }
    }

    #[test]
    fn dtr_plus_dense_is_complete() {
        let repeat = b"ACGTACGTACGTACGTACGTACGTAAGG";
        let mut seq = Vec::new();
        seq.extend_from_slice(repeat);
        seq.extend(std::iter::repeat(b'A').take(2000));
        seq.extend_from_slice(repeat);
        // genes covering ~80% of the contig
        let genes = vec![gene(1, (seq.len() as f64 * 0.8) as usize)];
        let r = trace(&seq, &genes);
        assert_eq!(r.terminal_repeat.kind, TrKind::Dtr);
        assert_eq!(r.confidence, TraceConfidence::Complete);
    }

    #[test]
    fn sparse_coding_is_fragment() {
        let seq = vec![b'A'; 5000];
        let genes = vec![gene(1, 300)]; // ~6% coding
        let r = trace(&seq, &genes);
        assert_eq!(r.confidence, TraceConfidence::Fragment);
    }
}
