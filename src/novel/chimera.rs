//! ## CHIMERA-SCAN — compositional junction detector
//!
//! A reference-free screen for contigs that join DNA of two different
//! compositions (mis-assembled inter-species chimeras). The contig is tiled
//! into windows (5 kb in `sentinel run`), each window gets its canonical
//! (reverse-complement folded) tetranucleotide profile, and every split of
//! the window series into a left and a right segment (≥ 2 windows each) is
//! scored by the between-segment sum of squares. At the best split the
//! junction strength is
//!
//! ```text
//! D = 1 − cos(mean profile of the left segment, mean profile of the right segment)
//! ```
//!
//! and the contig is flagged when `D` reaches the threshold (0.12 in
//! `sentinel run`). Using segment means rather than single windows is what
//! keeps genomic islands, rRNA operons and prophages inside an otherwise
//! homogeneous contig from being called: one atypical block barely moves
//! either segment mean, whereas the two halves of a chimera differ along
//! their whole length.
//!
//! Calibration (window 5 kb; three complete bacterial chromosomes of
//! 0.7–4.9 Mb and inter-species junctions between them; see AUDIT.md):
//! whole chromosomes give D = 0.001–0.035 (never flagged); window-aligned
//! pieces of 20 kb–1 Mb of the same chromosomes are flagged at a rate of
//! ≈ 1 % (pieces that end inside an island); junctions with ≥ 10 kb of each
//! genome are flagged 91–100 % of the time. Phage–host junctions are
//! usually NOT visible compositionally (phages ameliorate to their host);
//! CheckV's gene-based provirus detection is the tool for those.
//!
//! The per-window z-scores (most atypical window, fraction of windows beyond
//! 2σ) are reported for information only: in real genomes they mark islands,
//! not chimeras.

use crate::common::composition::{cosine, tetranucleotide_freqs, TETRA_DIM};

/// Minimum windows on each side of a split.
const MIN_SIDE: usize = 2;

/// CHIMERA-SCAN output.
#[derive(Clone, Debug)]
pub struct ChimeraResult {
    /// Whether the contig had enough windows (2 × [`MIN_SIDE`]) to be judged.
    pub assessed: bool,
    pub n_windows: usize,
    /// Junction strength `D` at the best split (0 when not assessed).
    pub junction_distance: f64,
    /// Fraction of the window-profile variance explained by the best split.
    pub junction_r2: f64,
    /// Base position of the best split (start of the right segment).
    pub split_position: Option<usize>,
    /// Largest per-window z-score of the distance to the contig centroid.
    pub max_zscore: f64,
    /// Fraction of windows beyond 2σ.
    pub outlier_fraction: f64,
    /// Breakpoint (base position) when flagged.
    pub breakpoint: Option<usize>,
    /// True when the junction strength reaches the threshold.
    pub flagged: bool,
}

/// Run CHIMERA-SCAN with a window size (bases, at least 500) and a junction
/// strength threshold `d_threshold` (see the module docs).
pub fn scan(seq: &[u8], window: usize, d_threshold: f64) -> ChimeraResult {
    let window = window.max(500);
    let m = seq.len() / window;
    let mut r = ChimeraResult {
        assessed: m >= 2 * MIN_SIDE,
        n_windows: m,
        junction_distance: 0.0,
        junction_r2: 0.0,
        split_position: None,
        max_zscore: 0.0,
        outlier_fraction: 0.0,
        breakpoint: None,
        flagged: false,
    };
    if !r.assessed {
        return r;
    }
    let windows: Vec<[f64; TETRA_DIM]> = (0..m).map(|j| tetranucleotide_freqs(&seq[j * window..(j + 1) * window])).collect();

    // centroid, total sum of squares, per-window z-scores (informational)
    let mut centroid = [0f64; TETRA_DIM];
    for w in &windows {
        for k in 0..TETRA_DIM {
            centroid[k] += w[k];
        }
    }
    for c in centroid.iter_mut() {
        *c /= m as f64;
    }
    let tss: f64 = windows.iter().map(|w| sq_dist(w, &centroid)).sum();
    let dists: Vec<f64> = windows.iter().map(|w| 1.0 - cosine(w, &centroid)).collect();
    let mean = dists.iter().sum::<f64>() / m as f64;
    let sd = (dists.iter().map(|d| (d - mean) * (d - mean)).sum::<f64>() / m as f64).sqrt();
    if sd > 1e-12 {
        r.max_zscore = dists.iter().map(|d| (d - mean) / sd).fold(0.0f64, f64::max);
        r.outlier_fraction = dists.iter().filter(|&&d| (d - mean) / sd > 2.0).count() as f64 / m as f64;
    }

    // best two-segment split by between-segment sum of squares
    let mut prefix = vec![[0f64; TETRA_DIM]; m + 1];
    for j in 0..m {
        for k in 0..TETRA_DIM {
            prefix[j + 1][k] = prefix[j][k] + windows[j][k];
        }
    }
    let mut best: Option<(f64, usize)> = None;
    for split in MIN_SIDE..=m - MIN_SIDE {
        let (nl, nr) = (split as f64, (m - split) as f64);
        let mut bss = 0.0;
        for k in 0..TETRA_DIM {
            let l = prefix[split][k] / nl;
            let rr = (prefix[m][k] - prefix[split][k]) / nr;
            bss += (l - rr) * (l - rr);
        }
        bss *= nl * nr / m as f64;
        if best.map_or(true, |(b, _)| bss > b) {
            best = Some((bss, split));
        }
    }
    let (bss, split) = best.expect("at least one split when assessed");
    let mut left = [0f64; TETRA_DIM];
    let mut right = [0f64; TETRA_DIM];
    for k in 0..TETRA_DIM {
        left[k] = prefix[split][k] / split as f64;
        right[k] = (prefix[m][k] - prefix[split][k]) / (m - split) as f64;
    }
    let both_informative = left.iter().any(|&x| x > 0.0) && right.iter().any(|&x| x > 0.0);
    r.junction_distance = if both_informative { 1.0 - cosine(&left, &right) } else { 0.0 };
    r.junction_r2 = if tss > 0.0 { bss / tss } else { 0.0 };
    r.split_position = Some(split * window);
    r.flagged = r.junction_distance >= d_threshold;
    r.breakpoint = if r.flagged { r.split_position } else { None };
    r
}

fn sq_dist(a: &[f64; TETRA_DIM], b: &[f64; TETRA_DIM]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| (x - y) * (x - y)).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Random sequence with a given G+C probability (xorshift, deterministic).
    fn random_seq(n: usize, gc: f64, seed: u64) -> Vec<u8> {
        let mut x = seed | 1;
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                let u = (x >> 11) as f64 / (1u64 << 53) as f64;
                let strong = u < gc;
                let pick = (x & 1) == 1;
                match (strong, pick) {
                    (true, true) => b'G',
                    (true, false) => b'C',
                    (false, true) => b'A',
                    (false, false) => b'T',
                }
            })
            .collect()
    }

    #[test]
    fn homogeneous_contig_is_not_flagged() {
        let r = scan(&random_seq(60_000, 0.5, 3), 5000, 0.12);
        assert!(r.assessed && !r.flagged, "D={}", r.junction_distance);
        assert!(r.junction_distance < 0.02);
    }

    #[test]
    fn short_contigs_are_not_assessed() {
        let r = scan(&random_seq(19_999, 0.5, 5), 5000, 0.12);
        assert!(!r.assessed && !r.flagged && r.split_position.is_none());
    }

    #[test]
    fn chimera_is_flagged_at_the_junction() {
        let mut seq = random_seq(30_000, 0.30, 7);
        seq.extend(random_seq(25_000, 0.68, 9));
        let r = scan(&seq, 5000, 0.12);
        assert!(r.flagged, "D={}", r.junction_distance);
        assert_eq!(r.breakpoint, Some(30_000));
        assert!(r.junction_r2 > 0.5);
    }

    #[test]
    fn internal_island_is_not_a_junction() {
        // 10 kb island of very different composition inside 100 kb
        let mut seq = random_seq(45_000, 0.5, 11);
        seq.extend(random_seq(10_000, 0.25, 13));
        seq.extend(random_seq(45_000, 0.5, 17));
        let r = scan(&seq, 5000, 0.12);
        assert!(r.max_zscore > 2.0, "the island windows are atypical");
        assert!(!r.flagged, "D={}", r.junction_distance);
    }

    #[test]
    fn all_n_contig_is_safe() {
        let r = scan(&vec![b'N'; 40_000], 5000, 0.12);
        assert!(r.assessed && !r.flagged);
        assert_eq!(r.junction_distance, 0.0);
    }
}
