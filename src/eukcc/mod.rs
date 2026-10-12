//! EukCC 2 completeness / contamination scoring.
//!
//! EukCC searches the predicted proteins against a marker set with
//! `hmmsearch --cut_ga --tblout` and counts **every** reported row
//! (`eukcc.py: compute_quality`):
//!
//! ```text
//! found          = Counter(query HMM of every tblout row)
//! completeness   = round(|found ∩ expected| / |expected| · 100, 2)
//! contamination  = round(|{m : found[m] ≠ 1} ∩ expected| / |expected| · 100, 2)
//! ```
//!
//! A protein that passes GA for two markers counts for both, exactly as in
//! EukCC. Every marker must carry a GA line (hmmsearch aborts otherwise).
//!
//! The marker set is supplied explicitly (`--markers <set.hmm>`, every model
//! in the file is expected). EukCC 2's automatic choice of the set by
//! phylogenetic placement on its reference tree (`epa-ng`) is not
//! implemented; pass the set EukCC would use for the lineage.

use std::collections::{HashMap, HashSet};
use std::path::Path;

use anyhow::{bail, Result};

use crate::common::pyfmt::py_round;
use crate::hmmer::{search_hmm_file_models, Cutoff, SearchConfig, Targets};

/// EukCC result for one genome.
#[derive(Clone, Debug)]
pub struct EukccSummary {
    pub marker_set: String,
    pub completeness: f64,
    pub contamination: f64,
    pub markers_expected: usize,
    /// Expected markers with at least one row.
    pub markers_found: usize,
    /// Expected markers with ≠ 1 rows (i.e. ≥ 2).
    pub markers_multicopy: usize,
    pub markers_missing: usize,
    /// Rows in the tblout (protein × marker pairs passing GA).
    pub rows: usize,
    pub orf_count: usize,
    /// Rows per found marker (marker file order).
    pub copies: Vec<(String, usize)>,
}

impl EukccSummary {
    pub fn mean_copies(&self) -> f64 {
        if self.copies.is_empty() {
            0.0
        } else {
            self.copies.iter().map(|c| c.1 as f64).sum::<f64>() / self.copies.len() as f64
        }
    }
}

/// EukCC's `compute_quality` on the query names of the tblout rows.
pub fn compute_quality(rows: &[&str], expected: &[String]) -> (f64, f64, usize, usize) {
    let exp: HashSet<&str> = expected.iter().map(|s| s.as_str()).collect();
    let mut found: HashMap<&str, usize> = HashMap::new();
    for q in rows {
        *found.entry(q).or_insert(0) += 1;
    }
    let found_marker = found.keys().filter(|k| exp.contains(*k)).count();
    let multiton = found.iter().filter(|(k, &v)| v != 1 && exp.contains(*k)).count();
    if exp.is_empty() {
        return (0.0, 0.0, 0, 0);
    }
    let n = exp.len() as f64;
    (py_round(found_marker as f64 / n * 100.0, 2), py_round(multiton as f64 / n * 100.0, 2), found_marker, multiton)
}

/// Search `proteins` against the marker set with GA cutoffs and score it.
pub fn run(proteins: &Targets, markers: &Path) -> Result<EukccSummary> {
    let cfg = SearchConfig { cutoff: Cutoff::Gathering, ..Default::default() };
    let (hits, models) = match search_hmm_file_models(markers, proteins, &cfg) {
        Ok(x) => x,
        Err(e) => bail!("EukCC marker search: {e:#} (every marker needs a GA line, as for hmmsearch --cut_ga)"),
    };
    if models.is_empty() {
        bail!("no HMMs in marker set {}", markers.display());
    }
    let expected: Vec<String> = models.iter().map(|m| m.name.clone()).collect();
    let rows: Vec<&str> = hits.iter().map(|h| h.profile.as_str()).collect();
    let (completeness, contamination, found, multi) = compute_quality(&rows, &expected);
    let mut counts: HashMap<&str, usize> = HashMap::new();
    for r in &rows {
        *counts.entry(r).or_insert(0) += 1;
    }
    let copies: Vec<(String, usize)> = expected.iter().filter_map(|m| counts.get(m.as_str()).map(|&c| (m.clone(), c))).collect();
    let distinct: HashSet<&str> = expected.iter().map(|s| s.as_str()).collect();
    Ok(EukccSummary {
        marker_set: markers.file_stem().and_then(|s| s.to_str()).unwrap_or("markers").to_string(),
        completeness,
        contamination,
        markers_expected: distinct.len(),
        markers_found: found,
        markers_multicopy: multi,
        markers_missing: distinct.len() - found,
        rows: rows.len(),
        orf_count: proteins.len(),
        copies,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_every_row_like_eukcc() {
        let expected: Vec<String> = ["m1", "m2", "m3", "m4"].iter().map(|s| s.to_string()).collect();
        // protein p1 passes GA for m1 and m2; m3 hit twice; extra non-expected row
        let rows = ["m1", "m2", "m3", "m3", "zz"];
        let (comp, cont, found, multi) = compute_quality(&rows, &expected);
        assert_eq!((comp, cont, found, multi), (75.0, 25.0, 3, 1));
        // 1/3 rounding as Python round(x, 2)
        let e3: Vec<String> = ["a", "b", "c"].iter().map(|s| s.to_string()).collect();
        assert_eq!(compute_quality(&["a"], &e3).0, 33.33);
    }

    #[test]
    fn missing_ga_is_an_error() {
        let dir = std::env::temp_dir().join(format!("sentinel_x_eukcc_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("m.hmm");
        let with_ga = include_str!("../hmmer/testdata/tiny.hmm");
        let t = Targets::new([("p1", &b"MKAGG"[..])]);
        std::fs::write(&p, with_ga).unwrap();
        let s = run(&t, &p).unwrap();
        assert_eq!(s.markers_expected, 1);
        let no_ga: String = with_ga.lines().filter(|l| !l.starts_with("GA ")).map(|l| format!("{l}\n")).collect();
        std::fs::write(&p, no_ga).unwrap();
        assert!(run(&t, &p).is_err());
        let _ = std::fs::remove_dir_all(&dir);
    }
}
