//! ## SAGE — Self-consistency And Genome-composition Estimator
//!
//! A reference-free contamination signal for genome bins. Contamination in a
//! metagenome-assembled genome arrives as whole contigs from other
//! organisms, and organisms differ in their tetranucleotide composition, so
//! SAGE asks which contigs share the genome's core composition:
//!
//! 1. every contig of at least 5 kb gets its canonical tetranucleotide
//!    profile;
//! 2. the core is seeded at the length-weighted medoid (the contig closest,
//!    by total length-weighted distance, to all others) and refined
//!    iteratively: a contig belongs to the core when `1 − cos(profile, core
//!    centroid)` is below a length-dependent threshold
//!    `τ(L) = 0.06 + 0.025·√(5 kb / L)` (short contigs are noisier), and the
//!    centroid is recomputed from the members (length-weighted) until the
//!    membership stops changing;
//! 3. **discordant bp %** = sequence outside the core / assessed sequence.
//!
//! Calibration (three complete bacterial genomes fragmented into contigs of
//! mean length 8–40 kb, with 5–45 % of contigs from another of the genomes
//! mixed in): estimates within −1.2 … +2.0 points of the true contaminant
//! fraction; single genomes score ≤ 2.7 % (mean ≤ 0.7 %). Thresholds come
//! from contig-to-genome-centroid distances: within-genome 99th percentile
//! 0.080 (5–10 kb contigs) … 0.007 (100 kb); between-genome 1st percentile
//! 0.089 … 0.119 (see AUDIT.md).
//!
//! Limits: contamination by a close relative (similar composition) is not
//! seen; plasmids, prophages and islands with atypical composition count as
//! discordant; contigs < 5 kb are not assessed; a single-contig genome has
//! nothing to compare (use CHIMERA-SCAN for within-contig junctions).
//!
//! Gene-level summaries are reported alongside: gene density, amino-acid
//! composition coherence (mean cosine of each protein's composition to the
//! proteome centroid) and near-identical protein pairs.

use crate::common::composition::{cosine, tetranucleotide_freqs, TETRA_DIM};
use crate::common::{Contig, Gene};
use rayon::prelude::*;

/// The 20 standard amino acids (composition-vector axes).
const AA: [u8; 20] = *b"ACDEFGHIKLMNPQRSTVWY";
/// Contigs shorter than this are not assessed (too noisy).
pub const MIN_CONTIG: usize = 5000;
/// At most this many contigs take part in the medoid search (longest first).
const MEDOID_CAP: usize = 2000;

/// SAGE output.
#[derive(Clone, Debug)]
pub struct SageResult {
    pub n_genes: usize,
    /// Genes per kilobase of analysed sequence.
    pub gene_density: f64,
    /// Mean cosine of per-protein amino-acid composition to the proteome
    /// centroid (proteins ≥ 60 aa; 1 = perfectly homogeneous).
    pub aa_coherence: f64,
    /// Near-identical protein pairs (composition cosine > 0.995, length
    /// ratio > 0.9), counted on a subsample of ≤ 4,000 proteins and scaled.
    pub redundant_pairs: usize,
    /// Contigs ≥ [`MIN_CONTIG`] and their total length.
    pub contigs_assessed: usize,
    pub assessed_bp: usize,
    /// Sequence outside the core composition, and those contigs (indices
    /// into the input slice, input order).
    pub discordant_bp: usize,
    pub discordant: Vec<usize>,
}

impl SageResult {
    /// Discordant bp as a percentage of the assessed sequence; `None` when
    /// fewer than two contigs could be assessed.
    pub fn discordant_pct(&self) -> Option<f64> {
        if self.contigs_assessed < 2 || self.assessed_bp == 0 {
            None
        } else {
            Some(100.0 * self.discordant_bp as f64 / self.assessed_bp as f64)
        }
    }
}

/// Length-dependent core-membership threshold.
fn tau(len: usize) -> f64 {
    0.06 + 0.025 * (MIN_CONTIG as f64 / len.max(1) as f64).sqrt()
}

/// Core membership for contig profiles (see the module docs).
fn core_members(profiles: &[[f64; TETRA_DIM]], lens: &[usize]) -> Vec<bool> {
    let n = profiles.len();
    if n == 0 {
        return Vec::new();
    }
    // medoid over the longest MEDOID_CAP contigs
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| (std::cmp::Reverse(lens[i]), i));
    order.truncate(MEDOID_CAP);
    let medoid = order
        .par_iter()
        .map(|&i| {
            let cost: f64 = order.iter().map(|&j| lens[j] as f64 * (1.0 - cosine(&profiles[i], &profiles[j]))).sum();
            (cost, i)
        })
        .min_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)))
        .map(|(_, i)| i)
        .unwrap_or(0);

    let mut centroid = profiles[medoid];
    let mut members: Vec<bool> = Vec::new();
    for _ in 0..50 {
        let next: Vec<bool> = (0..n).map(|i| 1.0 - cosine(&profiles[i], &centroid) < tau(lens[i])).collect();
        if next == members || !next.iter().any(|&m| m) {
            if members.is_empty() {
                members = next;
            }
            break;
        }
        members = next;
        let mut c = [0f64; TETRA_DIM];
        let mut w = 0f64;
        for i in (0..n).filter(|&i| members[i]) {
            for k in 0..TETRA_DIM {
                c[k] += lens[i] as f64 * profiles[i][k];
            }
            w += lens[i] as f64;
        }
        for x in c.iter_mut() {
            *x /= w;
        }
        centroid = c;
    }
    members
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

/// Run SAGE on a genome's contigs (may be empty, e.g. protein input) and its
/// predicted genes. `analysed_bp` is the sequence length the genes were
/// called on (0 when unknown).
pub fn sage(contigs: &[Contig], genes: &[Gene], analysed_bp: usize) -> SageResult {
    // ---- contig composition
    let assessed: Vec<usize> = (0..contigs.len()).filter(|&i| contigs[i].len() >= MIN_CONTIG).collect();
    let profiles: Vec<[f64; TETRA_DIM]> = assessed.par_iter().map(|&i| tetranucleotide_freqs(&contigs[i].seq)).collect();
    let lens: Vec<usize> = assessed.iter().map(|&i| contigs[i].len()).collect();
    let assessed_bp: usize = lens.iter().sum();
    let (discordant, discordant_bp) = if assessed.len() >= 2 {
        let members = core_members(&profiles, &lens);
        let d: Vec<usize> = assessed.iter().zip(&members).filter(|(_, &m)| !m).map(|(&i, _)| i).collect();
        let bp = d.iter().map(|&i| contigs[i].len()).sum();
        (d, bp)
    } else {
        (Vec::new(), 0)
    };

    // ---- gene-level summaries
    let used: Vec<&Gene> = genes.iter().filter(|g| g.aa.len() >= 60).collect();
    let comps: Vec<[f64; 20]> = used.iter().map(|g| aa_composition(&g.aa)).collect();
    let mut aa_coherence = 0.0;
    let mut redundant_pairs = 0usize;
    if !comps.is_empty() {
        let mut centroid = [0f64; 20];
        for c in &comps {
            for k in 0..20 {
                centroid[k] += c[k];
            }
        }
        for c in centroid.iter_mut() {
            *c /= comps.len() as f64;
        }
        aa_coherence = comps.iter().map(|c| cosine(c, &centroid)).sum::<f64>() / comps.len() as f64;
        let cap = 4000usize;
        let step = comps.len().div_ceil(cap).max(1);
        let sub: Vec<usize> = (0..comps.len()).step_by(step).collect();
        let pairs: usize = sub
            .par_iter()
            .enumerate()
            .map(|(a, &i)| {
                sub[a + 1..]
                    .iter()
                    .filter(|&&j| {
                        let (li, lj) = (used[i].aa.len(), used[j].aa.len());
                        li.min(lj) as f64 / li.max(lj) as f64 > 0.9 && cosine(&comps[i], &comps[j]) > 0.995
                    })
                    .count()
            })
            .sum();
        let frac = sub.len() as f64 / comps.len() as f64;
        redundant_pairs = (pairs as f64 / (frac * frac)).round() as usize;
    }

    SageResult {
        n_genes: genes.len(),
        gene_density: if analysed_bp > 0 { genes.len() as f64 / (analysed_bp as f64 / 1000.0) } else { 0.0 },
        aa_coherence,
        redundant_pairs,
        contigs_assessed: assessed.len(),
        assessed_bp,
        discordant_bp,
        discordant,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Random sequence with a given G+C probability (xorshift, deterministic).
    fn random_seq(n: usize, gc: f64, seed: u64) -> Vec<u8> {
        let mut x = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                let u = (x >> 11) as f64 / (1u64 << 53) as f64;
                match (u < gc, x & 1 == 1) {
                    (true, true) => b'G',
                    (true, false) => b'C',
                    (false, true) => b'A',
                    (false, false) => b'T',
                }
            })
            .collect()
    }

    fn contig(id: &str, seq: Vec<u8>) -> Contig {
        Contig { id: id.into(), header: id.into(), seq }
    }

    fn gene(aa: Vec<u8>) -> Gene {
        Gene {
            contig: "c".into(),
            contig_index: 0,
            index: 1,
            begin: 1,
            end: aa.len() * 3,
            strand: 1,
            aa,
            partial_left: false,
            partial_right: false,
            stop: true,
        }
    }

    #[test]
    fn single_composition_bin_is_clean() {
        let cs: Vec<Contig> =
            (0..12).map(|i| contig(&format!("c{i}"), random_seq(8_000 + 3_000 * i, 0.5, i as u64 + 1))).collect();
        let r = sage(&cs, &[], 0);
        assert_eq!(r.contigs_assessed, 12);
        assert_eq!(r.discordant_pct(), Some(0.0));
    }

    #[test]
    fn contaminant_contigs_are_found() {
        let mut cs: Vec<Contig> = (0..10).map(|i| contig(&format!("a{i}"), random_seq(20_000, 0.5, i as u64 + 1))).collect();
        cs.push(contig("b0", random_seq(30_000, 0.30, 100)));
        cs.push(contig("b1", random_seq(20_000, 0.30, 101)));
        cs.push(contig("short", random_seq(3_000, 0.30, 102))); // < 5 kb: not assessed
        let r = sage(&cs, &[], 0);
        assert_eq!(r.contigs_assessed, 12);
        assert_eq!(r.discordant, vec![10, 11]);
        let pct = r.discordant_pct().unwrap();
        assert!((pct - 100.0 * 50.0 / 250.0).abs() < 1e-9, "{pct}");
    }

    #[test]
    fn majority_defines_the_core_even_when_listed_last() {
        let mut cs: Vec<Contig> = (0..3).map(|i| contig(&format!("b{i}"), random_seq(10_000, 0.3, 50 + i as u64))).collect();
        cs.extend((0..8).map(|i| contig(&format!("a{i}"), random_seq(25_000, 0.55, i as u64 + 1))));
        let r = sage(&cs, &[], 0);
        assert_eq!(r.discordant, vec![0, 1, 2]);
    }

    #[test]
    fn nothing_to_compare_for_one_contig_or_proteins_only() {
        let r = sage(&[contig("c", random_seq(50_000, 0.5, 1))], &[], 50_000);
        assert_eq!(r.discordant_pct(), None);
        let genes: Vec<Gene> = (0..5).map(|i| gene(vec![b"ACDEFGHIKLMNPQRSTVWY"[i]; 80])).collect();
        let r = sage(&[], &genes, 0);
        assert_eq!((r.contigs_assessed, r.discordant_pct()), (0, None));
        assert_eq!(r.n_genes, 5);
    }

    #[test]
    fn gene_level_summaries() {
        let mut x = 7u64;
        let alphabet = b"MAAALLLGGGVVVKKKDDDEEERRRSSSTTTNNNQQQ";
        let genes: Vec<Gene> = (0..30)
            .map(|i| {
                gene(
                    (0..150 + 7 * i)
                        .map(|_| {
                            x ^= x << 13;
                            x ^= x >> 7;
                            x ^= x << 17;
                            alphabet[(x % alphabet.len() as u64) as usize]
                        })
                        .collect(),
                )
            })
            .collect();
        let r = sage(&[], &genes, 30_000);
        assert!(r.aa_coherence > 0.95, "coh={}", r.aa_coherence);
        assert!((r.gene_density - 1.0).abs() < 1e-12);
        // identical proteins are redundant pairs
        let same: Vec<Gene> = (0..4).map(|_| genes[0].clone()).collect();
        assert_eq!(sage(&[], &same, 0).redundant_pairs, 6);
    }
}
