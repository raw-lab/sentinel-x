//! # sentinel-busco
//!
//! Single-copy ortholog completeness in the spirit of **BUSCO**. Given a
//! lineage's marker orthologs (each an HMM with a score cutoff and an expected
//! protein length) and the HMM hits of a genome's predicted proteins, every
//! ortholog is classified as:
//!
//! * **Complete – Single**: exactly one protein passes the score cutoff at an
//!   acceptable length;
//! * **Complete – Duplicated**: two or more proteins pass;
//! * **Fragmented**: a protein passes the score cutoff but is too short;
//! * **Missing**: nothing passes.
//!
//! The lineage tables (`scores_cutoff`, `lengths_cutoff`) are small and may be
//! supplied at runtime; the categorisation logic itself is database-free.

use std::collections::{HashMap, HashSet};

use sentinel_core::{Domain, QcRecord, Quality};
use sentinel_hmmer::Hit;

/// Per-ortholog cutoffs from a BUSCO lineage dataset.
#[derive(Clone, Copy, Debug)]
pub struct OrthologCutoff {
    /// Minimum HMM bit score for a hit to count (`scores_cutoff`).
    pub score: f64,
    /// Expected protein length mean (`lengths_cutoff` zeta/sigma model).
    pub length_mean: f64,
    /// Length standard deviation; a hit shorter than `mean - 2*sigma` is
    /// Fragmented.
    pub length_sigma: f64,
}

/// A lineage: the set of marker orthologs and their cutoffs.
#[derive(Clone, Debug, Default)]
pub struct Lineage {
    pub name: String,
    pub cutoffs: HashMap<String, OrthologCutoff>,
}

impl Lineage {
    /// Build a lineage from parsed `scores_cutoff` and `lengths_cutoff` maps.
    /// `lengths` maps ortholog → (sigma, mean) as BUSCO stores them.
    pub fn from_tables(
        name: impl Into<String>,
        scores: &HashMap<String, f64>,
        lengths: &HashMap<String, (f64, f64)>,
    ) -> Self {
        let mut cutoffs = HashMap::new();
        for (og, &score) in scores {
            let (sigma, mean) = lengths.get(og).copied().unwrap_or((0.0, 0.0));
            cutoffs.insert(
                og.clone(),
                OrthologCutoff {
                    score,
                    length_mean: mean,
                    length_sigma: sigma,
                },
            );
        }
        Lineage {
            name: name.into(),
            cutoffs,
        }
    }

    /// Total number of markers in the lineage.
    pub fn total(&self) -> usize {
        self.cutoffs.len()
    }
}

/// Status of a single ortholog in a genome.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OrthologStatus {
    Single,
    Duplicated,
    Fragmented,
    Missing,
}

/// A BUSCO-style completeness summary.
#[derive(Clone, Debug, Default)]
pub struct BuscoSummary {
    pub total: usize,
    pub single: usize,
    pub duplicated: usize,
    pub fragmented: usize,
    pub missing: usize,
}

impl BuscoSummary {
    /// Complete = single + duplicated.
    pub fn complete(&self) -> usize {
        self.single + self.duplicated
    }
    /// Completeness percentage (Complete / total).
    pub fn completeness_pct(&self) -> f64 {
        if self.total == 0 {
            0.0
        } else {
            100.0 * self.complete() as f64 / self.total as f64
        }
    }
    /// Duplication percentage (Duplicated / total).
    pub fn duplication_pct(&self) -> f64 {
        if self.total == 0 {
            0.0
        } else {
            100.0 * self.duplicated as f64 / self.total as f64
        }
    }
    /// The compact BUSCO string, e.g. `C:96.4%[S:95.1%,D:1.3%],F:1.2%,M:2.4%`.
    pub fn busco_string(&self) -> String {
        let p = |n: usize| {
            if self.total == 0 {
                0.0
            } else {
                100.0 * n as f64 / self.total as f64
            }
        };
        format!(
            "C:{:.1}%[S:{:.1}%,D:{:.1}%],F:{:.1}%,M:{:.1}%,n:{}",
            self.completeness_pct(),
            p(self.single),
            p(self.duplicated),
            p(self.fragmented),
            p(self.missing),
            self.total
        )
    }
}

/// Map a hit's protein length: we approximate from the hit envelope when
/// present, else fall back to a supplied protein-length map.
fn hit_length(h: &Hit, protein_len: &HashMap<String, usize>) -> usize {
    if let Some((s, e)) = h.env {
        (e.saturating_sub(s) + 1).max(1)
    } else {
        protein_len.get(&h.target).copied().unwrap_or(0)
    }
}

/// Classify every ortholog in the lineage given HMM hits (profile = ortholog id,
/// target = gene id) and a gene→protein-length map. `min_len_frac` is the
/// fraction of `mean - 2*sigma` below which a passing hit is deemed Fragmented.
pub fn categorize(
    lineage: &Lineage,
    hits: &[Hit],
    protein_len: &HashMap<String, usize>,
) -> (BuscoSummary, HashMap<String, OrthologStatus>) {
    // For each ortholog, collect passing (gene, complete?) decisions.
    let mut complete_hits: HashMap<&str, HashSet<String>> = HashMap::new();
    let mut fragmented_hits: HashMap<&str, HashSet<String>> = HashMap::new();

    for h in hits {
        let og = h.profile.as_str();
        let Some(cut) = lineage.cutoffs.get(og) else {
            continue;
        };
        if (h.bits as f64) < cut.score {
            continue;
        }
        let len = hit_length(h, protein_len) as f64;
        let min_len = (cut.length_mean - 2.0 * cut.length_sigma).max(0.0);
        if cut.length_mean > 0.0 && len < min_len {
            fragmented_hits.entry(og).or_default().insert(h.target.clone());
        } else {
            complete_hits.entry(og).or_default().insert(h.target.clone());
        }
    }

    let mut statuses = HashMap::new();
    let mut summary = BuscoSummary {
        total: lineage.total(),
        ..Default::default()
    };

    for og in lineage.cutoffs.keys() {
        let n_complete = complete_hits.get(og.as_str()).map(|s| s.len()).unwrap_or(0);
        let n_frag = fragmented_hits.get(og.as_str()).map(|s| s.len()).unwrap_or(0);
        let status = if n_complete == 1 {
            summary.single += 1;
            OrthologStatus::Single
        } else if n_complete >= 2 {
            summary.duplicated += 1;
            OrthologStatus::Duplicated
        } else if n_frag >= 1 {
            summary.fragmented += 1;
            OrthologStatus::Fragmented
        } else {
            summary.missing += 1;
            OrthologStatus::Missing
        };
        statuses.insert(og.clone(), status);
    }

    (summary, statuses)
}

/// Build a single genome-level QC record from a BUSCO summary.
pub fn to_record(genome_id: &str, total_len: usize, domain: Domain, summary: &BuscoSummary) -> QcRecord {
    let mut rec = QcRecord::new(genome_id, total_len, domain, "busco");
    rec.completeness = Some(summary.completeness_pct());
    rec.contamination = Some(summary.duplication_pct()); // duplication ~ redundancy
    rec.method = "single-copy-orthologs".into();
    rec.quality = Quality::from_completeness(Some(summary.completeness_pct()));
    if summary.duplication_pct() > 10.0 {
        rec.warnings.push(format!(
            "high duplication ({:.1}%) — possible contamination or assembly redundancy",
            summary.duplication_pct()
        ));
    }
    rec.with_extra("busco", summary.busco_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lineage3() -> Lineage {
        let mut scores = HashMap::new();
        let mut lengths = HashMap::new();
        for og in ["OG1", "OG2", "OG3", "OG4"] {
            scores.insert(og.to_string(), 20.0);
            lengths.insert(og.to_string(), (10.0, 100.0)); // sigma=10, mean=100
        }
        Lineage::from_tables("test_odb", &scores, &lengths)
    }

    fn hit(profile: &str, target: &str, bits: f32, len: usize) -> Hit {
        Hit {
            profile: profile.into(),
            target: target.into(),
            bits,
            viterbi_bits: bits,
            evalue: 0.0,
            env: Some((1, len)),
        }
    }

    #[test]
    fn classifies_single_dup_frag_missing() {
        let lin = lineage3();
        let plen = HashMap::new();
        let hits = vec![
            hit("OG1", "g1", 50.0, 100), // single complete
            hit("OG2", "g2", 50.0, 100), // duplicated:
            hit("OG2", "g3", 50.0, 100),
            hit("OG3", "g4", 50.0, 40), // fragmented (len 40 < 100-20=80)
                                        // OG4 -> missing
        ];
        let (sum, st) = categorize(&lin, &hits, &plen);
        assert_eq!(sum.total, 4);
        assert_eq!(sum.single, 1);
        assert_eq!(sum.duplicated, 1);
        assert_eq!(sum.fragmented, 1);
        assert_eq!(sum.missing, 1);
        assert_eq!(st["OG1"], OrthologStatus::Single);
        assert_eq!(st["OG2"], OrthologStatus::Duplicated);
        assert_eq!(st["OG3"], OrthologStatus::Fragmented);
        assert_eq!(st["OG4"], OrthologStatus::Missing);
        assert!((sum.completeness_pct() - 50.0).abs() < 1e-9);
    }

    #[test]
    fn below_score_cutoff_is_missing() {
        let lin = lineage3();
        let plen = HashMap::new();
        let hits = vec![hit("OG1", "g1", 5.0, 100)]; // below cutoff 20
        let (sum, st) = categorize(&lin, &hits, &plen);
        assert_eq!(st["OG1"], OrthologStatus::Missing);
        assert_eq!(sum.complete(), 0);
    }
}
