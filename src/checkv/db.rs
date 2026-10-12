//! The CheckV reference database (`checkv-db-v1.x`), read in its released
//! layout:
//!
//! ```text
//! genome_db/checkv_reps.faa      reference proteins (<checkv_id>_<n>)
//! genome_db/checkv_reps.dmnd     optional DIAMOND copy (read in place)
//! genome_db/checkv_reps.tsv      checkv_id, type (circular|genbank), length, …
//! genome_db/checkv_error.tsv     length, aai, cov, count, error
//! hmm_db/checkv_hmms/<n>.hmm     80 HMM files
//! hmm_db/checkv_hmms.tsv         hmm, category (viral|microbial), score_cutoff, database
//! hmm_db/genome_lengths.tsv      hmm, genomes, cv, lengths (comma list)
//! ```

use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

/// HMM annotation info (`checkv_hmms.tsv`).
#[derive(Clone, Debug)]
pub struct HmmInfo {
    pub viral: bool,
    /// `max(score_cutoff, 25)` bits.
    pub cutoff: f64,
    pub database: String,
}

/// Genome-length distribution of an HMM's database targets
/// (`genome_lengths.tsv`).
#[derive(Clone, Debug)]
pub struct HmmLengths {
    pub genomes: usize,
    pub cv: f64,
    pub lengths: Vec<f64>,
}

/// A reference genome (`checkv_reps.tsv`).
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct RefGenome {
    pub length: usize,
    pub kind: String,
    /// 1 for `circular`, 2 for `genbank`.
    pub weight: Option<f64>,
}

/// AAI error-rate grid (`checkv_error.tsv`).
#[derive(Clone, Debug, Default)]
pub struct ErrorTable {
    pub rates: HashMap<(i64, i64, i64), Option<f64>>,
    pub lengths: Vec<i64>,
    pub aais: Vec<i64>,
    pub covs: Vec<i64>,
}

/// `take_closest`: nearest grid value, ties to the smaller.
fn take_closest(list: &[i64], x: f64) -> i64 {
    let pos = list.partition_point(|&v| (v as f64) < x);
    if pos == 0 {
        return list[0];
    }
    if pos == list.len() {
        return list[list.len() - 1];
    }
    let before = list[pos - 1] as f64;
    let after = list[pos] as f64;
    if (after - x) < (x - before) {
        list[pos]
    } else {
        list[pos - 1]
    }
}

impl ErrorTable {
    /// Error rate (`None` = NA) and confidence for a contig length, AAI and
    /// aligned fraction (CheckV `lookup_error_rate`).
    pub fn lookup(&self, length: usize, aai: f64, af: f64) -> (Option<f64>, &'static str) {
        if self.lengths.is_empty() {
            return (None, "low");
        }
        let k = (take_closest(&self.lengths, length as f64), take_closest(&self.aais, aai), take_closest(&self.covs, af));
        // a key absent from the grid is a KeyError in CheckV; report NA
        let err = self.rates.get(&k).copied().flatten();
        let conf = match err {
            None => "low",
            Some(e) if e <= 5.0 => "high",
            Some(e) if e <= 10.0 => "medium",
            _ => "low",
        };
        (err, conf)
    }
}

/// A loaded CheckV database.
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct CheckvDb {
    pub dir: PathBuf,
    pub hmm_files: Vec<PathBuf>,
    pub hmm_info: HashMap<String, HmmInfo>,
    pub genome_lengths: HashMap<String, HmmLengths>,
    pub reps: HashMap<String, RefGenome>,
    pub error: ErrorTable,
    /// Reference proteins for the AAI search (FASTA, else `.dmnd`).
    pub reps_proteins: PathBuf,
}

/// Read a TSV with a header into maps (CheckV reads these with
/// `csv.DictReader`).
fn read_tsv(path: &Path) -> Result<Vec<BTreeMap<String, String>>> {
    let txt = std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let mut lines = txt.lines();
    let header: Vec<String> = match lines.next() {
        Some(h) => h.split('\t').map(|s| s.trim().to_string()).collect(),
        None => bail!("{} is empty", path.display()),
    };
    let mut out = Vec::new();
    for (n, l) in lines.enumerate() {
        if l.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() < header.len() {
            bail!("{}: line {} has {} fields, header has {}", path.display(), n + 2, f.len(), header.len());
        }
        out.push(header.iter().cloned().zip(f.iter().map(|s| s.to_string())).collect());
    }
    Ok(out)
}

fn get<'a>(r: &'a BTreeMap<String, String>, k: &str, path: &Path) -> Result<&'a str> {
    r.get(k).map(|s| s.as_str()).with_context(|| format!("{}: missing column '{k}'", path.display()))
}

/// Does `dir` look like a CheckV database root?
pub fn is_checkv_dir(dir: &Path) -> bool {
    dir.join("genome_db/checkv_reps.tsv").is_file()
        && (dir.join("hmm_db/checkv_hmms").is_dir() || dir.join("hmm_db/checkv_hmms.hmm").is_file())
}

/// Files a usable CheckV release must contain that are missing from `dir`.
pub fn missing_files(dir: &Path) -> Vec<String> {
    let mut m = Vec::new();
    for f in ["genome_db/checkv_reps.tsv", "genome_db/checkv_error.tsv", "hmm_db/checkv_hmms.tsv", "hmm_db/genome_lengths.tsv"] {
        if !dir.join(f).is_file() {
            m.push(f.to_string());
        }
    }
    if !(dir.join("genome_db/checkv_reps.faa").is_file() || dir.join("genome_db/checkv_reps.dmnd").is_file()) {
        m.push("genome_db/checkv_reps.faa".into());
    }
    let hdir = dir.join("hmm_db/checkv_hmms");
    let has_hmms = std::fs::read_dir(&hdir)
        .map(|rd| {
            rd.flatten().any(|e| {
                let n = e.file_name().to_string_lossy().into_owned();
                n.ends_with(".hmm") || n.ends_with(".hmm.gz")
            })
        })
        .unwrap_or(false);
    if !(has_hmms || dir.join("hmm_db/checkv_hmms.hmm").is_file()) {
        m.push("hmm_db/checkv_hmms/*.hmm".into());
    }
    m
}

/// Find CheckV database roots under `root` (no symlink following); the
/// newest release (by version-aware name order) comes first.
pub fn find_checkv_dirs(root: &Path) -> Vec<PathBuf> {
    let mut found = Vec::new();
    let mut stack = vec![(root.to_path_buf(), 0usize)];
    while let Some((d, depth)) = stack.pop() {
        if is_checkv_dir(&d) {
            found.push(d);
            continue;
        }
        if depth >= 4 {
            continue;
        }
        let Ok(rd) = std::fs::read_dir(&d) else { continue };
        for e in rd.flatten() {
            let Ok(ft) = e.file_type() else { continue };
            // a link directly inside the root may point at an existing
            // installation; links deeper down are not followed
            if ft.is_dir() || (depth == 0 && ft.is_symlink() && e.path().is_dir()) {
                stack.push((e.path(), depth + 1));
            }
        }
    }
    found.sort_by_key(|p| std::cmp::Reverse(version_key(p)));
    found
}

/// Natural-order key of a path's file name (`checkv-db-v1.10` > `v1.5`).
fn version_key(p: &Path) -> Vec<(u64, String)> {
    let name = p.file_name().and_then(|s| s.to_str()).unwrap_or("");
    let mut out = Vec::new();
    let mut num = String::new();
    let mut txt = String::new();
    for c in name.chars() {
        if c.is_ascii_digit() {
            if !txt.is_empty() {
                out.push((0, std::mem::take(&mut txt)));
            }
            num.push(c);
        } else {
            if !num.is_empty() {
                out.push((num.parse().unwrap_or(0), String::new()));
                num.clear();
            }
            txt.push(c);
        }
    }
    if !num.is_empty() {
        out.push((num.parse().unwrap_or(0), String::new()));
    }
    if !txt.is_empty() {
        out.push((0, txt));
    }
    out
}

impl CheckvDb {
    /// Load a database root.
    pub fn load(dir: &Path) -> Result<CheckvDb> {
        // HMM files: hmm_db/checkv_hmms/*.hmm (numeric order), or a single
        // combined hmm_db/checkv_hmms.hmm
        let hdir = dir.join("hmm_db/checkv_hmms");
        let mut hmm_files: Vec<PathBuf> = Vec::new();
        if hdir.is_dir() {
            for e in std::fs::read_dir(&hdir)?.flatten() {
                let p = e.path();
                let n = p.file_name().and_then(|s| s.to_str()).unwrap_or("");
                if n.ends_with(".hmm") || n.ends_with(".hmm.gz") {
                    hmm_files.push(p);
                }
            }
            hmm_files.sort_by_key(|p| version_key(p));
        } else if dir.join("hmm_db/checkv_hmms.hmm").is_file() {
            hmm_files.push(dir.join("hmm_db/checkv_hmms.hmm"));
        }
        if hmm_files.is_empty() {
            bail!("no HMM files under {}", hdir.display());
        }

        let p = dir.join("hmm_db/checkv_hmms.tsv");
        let mut hmm_info = HashMap::new();
        for r in read_tsv(&p)? {
            let cutoff: f64 = get(&r, "score_cutoff", &p)?.parse().unwrap_or(0.0);
            let cat = get(&r, "category", &p)?;
            if cat != "viral" && cat != "microbial" {
                bail!("{}: unknown category '{cat}'", p.display());
            }
            hmm_info.insert(
                get(&r, "hmm", &p)?.to_string(),
                HmmInfo {
                    viral: cat == "viral",
                    cutoff: cutoff.max(25.0),
                    database: r.get("database").cloned().unwrap_or_else(|| "NA".into()),
                },
            );
        }

        let p = dir.join("hmm_db/genome_lengths.tsv");
        let mut genome_lengths = HashMap::new();
        if p.is_file() {
            for r in read_tsv(&p)? {
                let lengths: Vec<f64> = get(&r, "lengths", &p)?
                    .split(',')
                    .filter(|s| !s.trim().is_empty())
                    .map(|s| s.trim().parse::<i64>().map(|v| v as f64))
                    .collect::<std::result::Result<_, _>>()
                    .with_context(|| format!("{}: bad lengths", p.display()))?;
                genome_lengths.insert(
                    get(&r, "hmm", &p)?.to_string(),
                    HmmLengths {
                        genomes: get(&r, "genomes", &p)?.parse().unwrap_or(0),
                        cv: get(&r, "cv", &p)?.parse().unwrap_or(0.0),
                        lengths,
                    },
                );
            }
        }

        let p = dir.join("genome_db/checkv_reps.tsv");
        let mut reps = HashMap::new();
        for r in read_tsv(&p)? {
            let kind = get(&r, "type", &p)?.to_string();
            let weight = match kind.as_str() {
                "circular" => Some(1.0),
                "genbank" => Some(2.0),
                _ => None,
            };
            reps.insert(
                get(&r, "checkv_id", &p)?.to_string(),
                RefGenome { length: get(&r, "length", &p)?.parse().unwrap_or(0), kind, weight },
            );
        }

        let p = dir.join("genome_db/checkv_error.tsv");
        let mut error = ErrorTable::default();
        if p.is_file() {
            for r in read_tsv(&p)? {
                let key = (
                    get(&r, "length", &p)?.parse::<i64>()?,
                    get(&r, "aai", &p)?.parse::<i64>()?,
                    get(&r, "cov", &p)?.parse::<i64>()?,
                );
                let count: i64 = get(&r, "count", &p)?.parse().unwrap_or(0);
                let e: Option<f64> = if count > 100 { get(&r, "error", &p)?.parse().ok() } else { None };
                error.rates.insert(key, e);
            }
            let mut l: Vec<i64> = error.rates.keys().map(|k| k.0).collect();
            let mut a: Vec<i64> = error.rates.keys().map(|k| k.1).collect();
            let mut c: Vec<i64> = error.rates.keys().map(|k| k.2).collect();
            for v in [&mut l, &mut a, &mut c] {
                v.sort_unstable();
                v.dedup();
            }
            error.lengths = l;
            error.aais = a;
            error.covs = c;
        }

        // the FASTA ships with every release; a locally built .dmnd may be
        // from any DIAMOND version, so it is only the fallback
        let dmnd = dir.join("genome_db/checkv_reps.dmnd");
        let faa = dir.join("genome_db/checkv_reps.faa");
        let reps_proteins = if faa.is_file() {
            faa
        } else if dmnd.is_file() {
            dmnd
        } else {
            bail!("neither genome_db/checkv_reps.dmnd nor checkv_reps.faa found in {}", dir.display());
        };

        Ok(CheckvDb { dir: dir.to_path_buf(), hmm_files, hmm_info, genome_lengths, reps, error, reps_proteins })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn closest_ties_to_smaller() {
        let v = [10, 20, 30];
        assert_eq!(take_closest(&v, 15.0), 10);
        assert_eq!(take_closest(&v, 16.0), 20);
        assert_eq!(take_closest(&v, 1.0), 10);
        assert_eq!(take_closest(&v, 99.0), 30);
        assert_eq!(take_closest(&v, 20.0), 20);
    }

    #[test]
    fn version_order() {
        let a = PathBuf::from("/x/checkv-db-v1.10");
        let b = PathBuf::from("/x/checkv-db-v1.5");
        assert!(version_key(&a) > version_key(&b));
    }
}
