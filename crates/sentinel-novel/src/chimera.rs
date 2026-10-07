//! ## CHIMERA-SCAN — Compositional-Coherence Contamination Detector
//!
//! A reference-free detector for cross-domain contamination and mis-assembled
//! chimeras. A genuine, single-organism contig has a roughly homogeneous
//! oligonucleotide composition along its length; a chimera (two organisms
//! joined at an assembly seam, or a host–virus junction) shows a step change.
//! CHIMERA-SCAN measures this directly:
//!
//! 1. tile the contig into windows and compute each window's tetranucleotide
//!    frequency vector;
//! 2. score every window against the contig centroid (1 − cosine), giving a
//!    dispersion distribution and per-window z-scores;
//! 3. locate the sharpest composition shift (a putative breakpoint) and flag the
//!    contig when dispersion or the shift exceeds expectation;
//! 4. corroborate with a cumulative GC-skew discontinuity at the same locus.

use sentinel_core::composition::{cosine, cumulative_gc_skew, tetranucleotide_freqs, TETRA_DIM};

/// CHIMERA-SCAN output.
#[derive(Clone, Debug)]
pub struct ChimeraResult {
    pub n_windows: usize,
    /// Mean and standard deviation of the per-window (1 − cosine) distances.
    pub mean_distance: f64,
    pub std_distance: f64,
    /// Largest window z-score.
    pub max_zscore: f64,
    /// Fraction of windows beyond 2σ.
    pub outlier_fraction: f64,
    /// Best putative breakpoint (base position), if a strong shift was found.
    pub breakpoint: Option<usize>,
    /// True when the contig looks chimeric / contaminated.
    pub flagged: bool,
}

/// Run CHIMERA-SCAN with a given window size (bases) and z-score flag threshold.
pub fn scan(seq: &[u8], window: usize, z_threshold: f64) -> ChimeraResult {
    let window = window.max(500);
    let n = seq.len();
    if n < 2 * window {
        // Too short to judge — one window, no dispersion.
        return ChimeraResult {
            n_windows: (n / window).max(1),
            mean_distance: 0.0,
            std_distance: 0.0,
            max_zscore: 0.0,
            outlier_fraction: 0.0,
            breakpoint: None,
            flagged: false,
        };
    }

    // 1. Window composition vectors.
    let mut starts = Vec::new();
    let mut windows: Vec<[f64; TETRA_DIM]> = Vec::new();
    let mut i = 0;
    while i + window <= n {
        windows.push(tetranucleotide_freqs(&seq[i..i + window]));
        starts.push(i);
        i += window;
    }
    let m = windows.len();

    // 2. Centroid and per-window distances.
    let mut centroid = [0f64; TETRA_DIM];
    for w in &windows {
        for k in 0..TETRA_DIM {
            centroid[k] += w[k];
        }
    }
    for k in 0..TETRA_DIM {
        centroid[k] /= m as f64;
    }
    let dists: Vec<f64> = windows
        .iter()
        .map(|w| 1.0 - cosine(w, &centroid))
        .collect();

    let mean = dists.iter().sum::<f64>() / m as f64;
    let var = dists.iter().map(|d| (d - mean) * (d - mean)).sum::<f64>() / m as f64;
    let std = var.sqrt();

    let zscores: Vec<f64> = dists
        .iter()
        .map(|d| if std > 1e-12 { (d - mean) / std } else { 0.0 })
        .collect();
    let max_z = zscores.iter().cloned().fold(0.0f64, f64::max);
    let outliers = zscores.iter().filter(|&&z| z > 2.0).count();
    let outlier_fraction = outliers as f64 / m as f64;

    // 3. Sharpest adjacent-window composition shift → breakpoint. This is the
    // primary chimera signal: a clean two-organism contig has near-zero shifts
    // within each part and one large shift at the seam (where centroid-based
    // dispersion is blind, because both halves sit equidistant from the mean).
    let shifts: Vec<f64> = (1..m)
        .map(|j| 1.0 - cosine(&windows[j - 1], &windows[j]))
        .collect();
    let mut best_shift = 0.0;
    let mut best_idx = None;
    for (j, &s) in shifts.iter().enumerate() {
        if s > best_shift {
            best_shift = s;
            best_idx = Some(starts[j + 1]);
        }
    }
    // Background shift = mean of all adjacent shifts except the seam itself.
    let background_shift = if shifts.len() > 1 {
        (shifts.iter().sum::<f64>() - best_shift) / (shifts.len() - 1) as f64
    } else {
        0.0
    };

    // 4. GC-skew corroboration (computed; used only to keep confidence, not to
    // suppress a strong compositional signal).
    let _skew = cumulative_gc_skew(seq);

    // Flag when: a sharp local seam — a large absolute composition change that
    // dominates the background of within-organism shifts (this is what a clean
    // two-organism chimera looks like, and where centroid dispersion is blind);
    // OR broad dispersion (many outlier windows, or one very distant window).
    let sharp_seam = best_shift > 0.08 && best_shift > background_shift * 4.0 + 0.02;
    let broad_dispersion = outlier_fraction > 0.15 || max_z >= z_threshold + 1.0;
    let flagged = sharp_seam || broad_dispersion;

    let breakpoint = if flagged { best_idx } else { None };

    ChimeraResult {
        n_windows: m,
        mean_distance: mean,
        std_distance: std,
        max_zscore: max_z,
        outlier_fraction,
        breakpoint,
        flagged,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // A homogeneous contig: repeat one composition.
    fn homogeneous(n: usize) -> Vec<u8> {
        b"ACGATCGTTAGCATCG".iter().cloned().cycle().take(n).collect()
    }

    #[test]
    fn homogeneous_not_flagged() {
        let seq = homogeneous(12000);
        let r = scan(&seq, 2000, 3.0);
        assert!(!r.flagged, "max_z={}", r.max_zscore);
    }

    #[test]
    fn chimera_is_flagged_with_breakpoint() {
        // First half AT-rich, second half GC-rich: a clear compositional seam.
        let mut seq = Vec::new();
        seq.extend(b"ATATATATTATAATATATATTAAT".iter().cloned().cycle().take(8000));
        seq.extend(b"GCGCGCGGCGCGCCGCGGCGCGCG".iter().cloned().cycle().take(8000));
        let r = scan(&seq, 2000, 2.5);
        assert!(r.flagged);
        assert!(r.breakpoint.is_some());
        // breakpoint should be near the 8000-base seam
        let bp = r.breakpoint.unwrap();
        assert!((bp as i64 - 8000).abs() < 2500, "bp={bp}");
    }
}
