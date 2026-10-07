//! ## SAGE — Single-copy Anchor-Gene & Codon-Coherence Estimator
//!
//! A reference-free contamination and self-consistency estimator. Rather than
//! matching genes to a marker database, SAGE asks whether a genome's predicted
//! proteins *look like they came from one organism*:
//!
//! * **Codon / amino-acid coherence** — proteins from a single genome share a
//!   characteristic amino-acid composition. SAGE builds each protein's 20-dim
//!   composition vector and measures how tightly they cluster around the
//!   proteome centroid; loose clustering betrays mixed organisms.
//! * **Compositional outliers** — individual proteins lying far from the
//!   centroid (high z-score) are candidate contaminants; their fraction is a
//!   database-free contamination proxy.
//! * **Single-copy violations** — near-identical proteins occurring more than
//!   once (very high composition cosine at similar length) indicate redundancy
//!   or a duplicated/contaminating lineage.
//!
//! When a domain-typical gene count is supplied, SAGE also reports a coarse
//! completeness estimate from observed-vs-expected gene density.

use sentinel_core::Gene;

/// The 20 standard amino acids (composition-vector axes).
const AA: [u8; 20] = *b"ACDEFGHIKLMNPQRSTVWY";

/// SAGE output.
#[derive(Clone, Debug)]
pub struct SageResult {
    pub n_genes: usize,
    /// Genes per kilobase of analysed sequence.
    pub gene_density: f64,
    /// Mean cosine of per-protein composition to the centroid (1 = perfectly
    /// homogeneous).
    pub codon_coherence: f64,
    /// Number / fraction of compositional-outlier proteins (z > 2).
    pub n_outliers: usize,
    pub outlier_fraction: f64,
    /// Near-identical protein pairs (single-copy redundancy).
    pub redundant_pairs: usize,
    /// Database-free contamination proxy (%), from outliers + redundancy.
    pub est_contamination: f64,
    /// Completeness estimate (%) when an expected gene count was supplied.
    pub est_completeness: Option<f64>,
}

/// Amino-acid composition (L1-normalised 20-vector) of a protein.
fn aa_composition(protein: &[u8]) -> [f64; 20] {
    let mut counts = [0f64; 20];
    let mut total = 0f64;
    for &b in protein {
        if let Some(p) = AA.iter().position(|&a| a == b) {
            counts[p] += 1.0;
            total += 1.0;
        }
    }
    if total > 0.0 {
        for c in counts.iter_mut() {
            *c /= total;
        }
    }
    counts
}

fn cosine20(a: &[f64; 20], b: &[f64; 20]) -> f64 {
    let mut dot = 0.0;
    let mut na = 0.0;
    let mut nb = 0.0;
    for k in 0..20 {
        dot += a[k] * b[k];
        na += a[k] * a[k];
        nb += b[k] * b[k];
    }
    if na == 0.0 || nb == 0.0 {
        0.0
    } else {
        dot / (na.sqrt() * nb.sqrt())
    }
}

/// Run SAGE. `analysed_bp` is the total sequence length the genes were called
/// on; `expected_genes` is an optional domain-typical gene count for the
/// completeness proxy.
pub fn sage(genes: &[Gene], analysed_bp: usize, expected_genes: Option<f64>) -> SageResult {
    let n = genes.len();
    if n == 0 {
        return SageResult {
            n_genes: 0,
            gene_density: 0.0,
            codon_coherence: 0.0,
            n_outliers: 0,
            outlier_fraction: 0.0,
            redundant_pairs: 0,
            est_contamination: 0.0,
            est_completeness: expected_genes.map(|_| 0.0),
        };
    }

    // Per-protein composition vectors (skip very short proteins).
    let comps: Vec<[f64; 20]> = genes
        .iter()
        .map(|g| aa_composition(&g.aa))
        .collect();

    // Centroid.
    let mut centroid = [0f64; 20];
    for c in &comps {
        for k in 0..20 {
            centroid[k] += c[k];
        }
    }
    for k in 0..20 {
        centroid[k] /= n as f64;
    }

    // Distances (1 - cosine) to centroid → coherence & outliers.
    let dists: Vec<f64> = comps.iter().map(|c| 1.0 - cosine20(c, &centroid)).collect();
    let coherence = comps.iter().map(|c| cosine20(c, &centroid)).sum::<f64>() / n as f64;

    let mean = dists.iter().sum::<f64>() / n as f64;
    let var = dists.iter().map(|d| (d - mean) * (d - mean)).sum::<f64>() / n as f64;
    let std = var.sqrt();
    let n_outliers = dists
        .iter()
        .filter(|&&d| std > 1e-12 && (d - mean) / std > 2.0)
        .count();
    let outlier_fraction = n_outliers as f64 / n as f64;

    // Single-copy redundancy: near-identical composition at similar length.
    // O(n^2) but n (genes per contig/genome) is modest; cap work for safety.
    let mut redundant_pairs = 0usize;
    let cap = 4000usize;
    if n <= cap {
        for i in 0..n {
            for j in (i + 1)..n {
                let li = genes[i].aa.len().max(1);
                let lj = genes[j].aa.len().max(1);
                let lr = li.min(lj) as f64 / li.max(lj) as f64;
                if lr > 0.9 && cosine20(&comps[i], &comps[j]) > 0.995 && li > 60 {
                    redundant_pairs += 1;
                }
            }
        }
    }

    let gene_density = if analysed_bp > 0 {
        n as f64 / (analysed_bp as f64 / 1000.0)
    } else {
        0.0
    };

    // Contamination proxy: outliers dominate; redundancy adds a smaller term.
    let redun_frac = redundant_pairs as f64 / n as f64;
    let est_contamination = ((outlier_fraction + 0.5 * redun_frac) * 100.0).min(100.0);

    let est_completeness = expected_genes.map(|e| {
        if e <= 0.0 {
            0.0
        } else {
            (n as f64 / e * 100.0).min(100.0)
        }
    });

    SageResult {
        n_genes: n,
        gene_density,
        codon_coherence: coherence,
        n_outliers,
        outlier_fraction,
        redundant_pairs,
        est_contamination,
        est_completeness,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gene_from(aa: &[u8]) -> Gene {
        Gene {
            contig: "c".into(),
            index: 1,
            begin: 1,
            end: aa.len() * 3,
            strand: 1,
            aa: aa.to_vec(),
            partial: false,
        }
    }

    #[test]
    fn coherent_proteome_high_coherence_low_contam() {
        // Proteins all drawn from a similar composition.
        let base = b"MAAALLLGGGVVVKKKDDDEEERRRSSSTTTNNNQQQ";
        let genes: Vec<Gene> = (0..30).map(|_| gene_from(base)).collect();
        let r = sage(&genes, 30 * base.len() * 3, None);
        assert!(r.codon_coherence > 0.98, "coh={}", r.codon_coherence);
        assert!(r.est_contamination < 5.0, "contam={}", r.est_contamination);
    }

    #[test]
    fn outlier_protein_raises_contamination() {
        let base = b"MAAALLLGGGVVVKKKDDDEEERRRSSSTTT";
        let mut genes: Vec<Gene> = (0..30).map(|_| gene_from(base)).collect();
        // add a wildly different composition protein (all tryptophan/cysteine)
        genes.push(gene_from(b"WWWWWWWWWWWWWWWWWWWWCCCCCCCCCCCCCCCC"));
        let r = sage(&genes, 1_000_000, None);
        assert!(r.n_outliers >= 1, "outliers={}", r.n_outliers);
    }

    #[test]
    fn completeness_from_expected() {
        let base = b"MAAALLLGGG";
        let genes: Vec<Gene> = (0..50).map(|_| gene_from(base)).collect();
        let r = sage(&genes, 100_000, Some(100.0));
        assert!((r.est_completeness.unwrap() - 50.0).abs() < 1e-9);
    }
}
