//! `hmmsearch`-equivalent search layer: profiles × target proteins, with
//! HMMER's reporting thresholds (E-value, bit score, or the per-model GA / NC /
//! TC cutoffs), sequence and domain E-values, and streaming over multi-GB HMM
//! libraries.
//!
//! Semantics follow HMMER 3.4:
//! * a target is reported when its sequence E-value `≤ E` (or score `≥ T`, or
//!   score `≥ cutoff[0]` under `--cut_*`);
//! * a reported target's domain is reported when its conditional E-value
//!   (P × number of reported targets) `≤ domE` (or domain score `≥ domT`, or
//!   `≥ cutoff[1]`);
//! * `Z` (database size for E-values) is the number of target sequences.

use std::path::Path;

use anyhow::{bail, Result};
use rayon::prelude::*;

use crate::hmmer::alphabet::encode;
use crate::hmmer::model::Profile;
use crate::hmmer::parse::open_hmm;
use crate::hmmer::pipeline::{score_target, PipelineConfig, RawHit, Workspace};

/// Reporting threshold (hmmsearch `-E/--domE`, `-T/--domT`, `--cut_ga|nc|tc`).
#[derive(Clone, Copy, Debug)]
pub enum Cutoff {
    Evalue { seq: f64, dom: f64 },
    Bits { seq: f32, dom: f32 },
    Gathering,
    Noise,
    Trusted,
}

impl Default for Cutoff {
    fn default() -> Self {
        Cutoff::Evalue { seq: 10.0, dom: 10.0 }
    }
}

/// Search configuration.
#[derive(Clone, Copy, Debug, Default)]
pub struct SearchConfig {
    pub cutoff: Cutoff,
    pub pipeline: PipelineConfig,
    /// Override the database size used for E-values (hmmsearch `-Z`).
    pub z: Option<f64>,
}

/// Pre-digitised target sequences.
pub struct Targets {
    pub ids: Vec<String>,
    pub dsq: Vec<Vec<u8>>,
}

impl Targets {
    pub fn new<'a, I>(pairs: I) -> Targets
    where
        I: IntoIterator<Item = (&'a str, &'a [u8])>,
    {
        let mut ids = Vec::new();
        let mut dsq = Vec::new();
        for (id, s) in pairs {
            ids.push(id.to_string());
            dsq.push(encode(s));
        }
        Targets { ids, dsq }
    }

    pub fn from_pairs(pairs: &[(String, Vec<u8>)]) -> Targets {
        Targets::new(pairs.iter().map(|(a, b)| (a.as_str(), b.as_slice())))
    }

    pub fn len(&self) -> usize {
        self.ids.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }
}

/// One domain of a reported hit (hmmsearch `--domtblout` row).
#[derive(Clone, Debug)]
pub struct DomainHit {
    pub score: f32,
    pub bias: f32,
    pub c_evalue: f64,
    pub i_evalue: f64,
    pub hmm_from: usize,
    pub hmm_to: usize,
    pub ali_from: usize,
    pub ali_to: usize,
    pub env_from: usize,
    pub env_to: usize,
    /// Passes the domain reporting threshold.
    pub reported: bool,
}

/// A reported target (hmmsearch `--tblout` row).
#[derive(Clone, Debug)]
pub struct Hit {
    pub profile: String,
    pub acc: Option<String>,
    /// Model length M.
    pub qlen: usize,
    pub target: String,
    pub target_index: usize,
    pub tlen: usize,
    /// Full-sequence bit score.
    pub score: f32,
    /// Full-sequence null2 bias (bits).
    pub bias: f32,
    pub evalue: f64,
    pub nexpected: f32,
    pub domains: Vec<DomainHit>,
}

impl Hit {
    /// Domains that pass the domain reporting threshold.
    pub fn reported_domains(&self) -> impl Iterator<Item = &DomainHit> {
        self.domains.iter().filter(|d| d.reported)
    }

    /// Best domain score (any domain).
    pub fn best_domain_score(&self) -> f32 {
        self.domains.iter().map(|d| d.score).fold(f32::NEG_INFINITY, f32::max)
    }
}

/// Resolve the (sequence, domain) thresholds for one profile.
fn thresholds(p: &Profile, c: Cutoff) -> Result<Threshold> {
    Ok(match c {
        Cutoff::Evalue { seq, dom } => Threshold::E(seq, dom),
        Cutoff::Bits { seq, dom } => Threshold::T(seq, dom),
        Cutoff::Gathering => match p.ga {
            Some((a, b)) => Threshold::T(a, b),
            None => bail!("GA bit thresholds unavailable on model {}", p.name),
        },
        Cutoff::Noise => match p.nc {
            Some((a, b)) => Threshold::T(a, b),
            None => bail!("NC bit thresholds unavailable on model {}", p.name),
        },
        Cutoff::Trusted => match p.tc {
            Some((a, b)) => Threshold::T(a, b),
            None => bail!("TC bit thresholds unavailable on model {}", p.name),
        },
    })
}

#[derive(Clone, Copy)]
enum Threshold {
    E(f64, f64),
    T(f32, f32),
}

/// Turn the raw pipeline survivors of one profile into reported hits.
fn report(p: &Profile, raws: Vec<(usize, RawHit)>, targets: &Targets, th: Threshold, z: f64) -> Vec<Hit> {
    let mut hits: Vec<Hit> = Vec::new();
    for (ti, r) in raws {
        let evalue = r.lnp.exp() * z;
        let keep = match th {
            Threshold::E(e, _) => evalue <= e,
            Threshold::T(t, _) => r.bits >= t,
        };
        if !keep {
            continue;
        }
        hits.push(Hit {
            profile: p.name.clone(),
            acc: p.acc.clone(),
            qlen: p.m,
            target: targets.ids[ti].clone(),
            target_index: ti,
            tlen: targets.dsq[ti].len(),
            score: r.bits,
            bias: r.bias,
            evalue,
            nexpected: r.nexpected,
            domains: r
                .domains
                .iter()
                .map(|d| DomainHit {
                    score: d.bits,
                    bias: d.bias,
                    c_evalue: f64::NAN,
                    i_evalue: d.lnp.exp() * z,
                    hmm_from: d.kali,
                    hmm_to: d.mali,
                    ali_from: d.iali,
                    ali_to: d.jali,
                    env_from: d.ienv,
                    env_to: d.jenv,
                    reported: false,
                })
                .collect(),
        });
    }
    let dom_z = hits.len().max(1) as f64;
    for h in &mut hits {
        for d in &mut h.domains {
            let p_dom = d.i_evalue / z;
            d.c_evalue = p_dom * dom_z;
            d.reported = match th {
                Threshold::E(_, de) => d.c_evalue <= de,
                Threshold::T(_, dt) => d.score >= dt,
            };
        }
    }
    hits
}

/// Search in-memory profiles against the targets.
pub fn search_profiles(profiles: &[Profile], targets: &Targets, cfg: &SearchConfig) -> Result<Vec<Hit>> {
    let z = cfg.z.unwrap_or(targets.len() as f64).max(1.0);
    let ths: Vec<Threshold> = profiles.iter().map(|p| thresholds(p, cfg.cutoff)).collect::<Result<_>>()?;
    let nt = targets.len();
    let pairs: Vec<(usize, usize, RawHit)> = (0..profiles.len() * nt)
        .into_par_iter()
        .map_init(Workspace::default, |ws, idx| {
            let (pi, ti) = (idx / nt, idx % nt);
            score_target(&profiles[pi], &targets.dsq[ti], &cfg.pipeline, ws).map(|r| (pi, ti, r))
        })
        .flatten()
        .collect();
    let mut per: Vec<Vec<(usize, RawHit)>> = (0..profiles.len()).map(|_| Vec::new()).collect();
    for (pi, ti, r) in pairs {
        per[pi].push((ti, r));
    }
    let mut out = Vec::new();
    for (pi, raws) in per.into_iter().enumerate() {
        out.extend(report(&profiles[pi], raws, targets, ths[pi], z));
    }
    Ok(out)
}

/// Stream every profile of an HMM file (optionally gzip-compressed) and search
/// it against the targets. Returns `(hits, number_of_models)`.
pub fn search_hmm_file(path: &Path, targets: &Targets, cfg: &SearchConfig) -> Result<(Vec<Hit>, usize)> {
    let (hits, models) = search_hmm_file_models(path, targets, cfg)?;
    Ok((hits, models.len()))
}

/// Identity of a searched model (name, accession, length).
#[derive(Clone, Debug)]
pub struct ModelInfo {
    pub name: String,
    pub acc: Option<String>,
    pub m: usize,
}

/// [`search_hmm_file`], also returning every model's identity in file order.
pub fn search_hmm_file_models(path: &Path, targets: &Targets, cfg: &SearchConfig) -> Result<(Vec<Hit>, Vec<ModelInfo>)> {
    let mut reader = open_hmm(path)?;
    let mut hits = Vec::new();
    let mut models = Vec::new();
    loop {
        let mut batch = Vec::with_capacity(256);
        for _ in 0..256 {
            match reader.next() {
                Some(h) => batch.push(h?.configure()),
                None => break,
            }
        }
        if batch.is_empty() {
            break;
        }
        models.extend(batch.iter().map(|p| ModelInfo { name: p.name.clone(), acc: p.acc.clone(), m: p.m }));
        hits.extend(search_profiles(&batch, targets, cfg)?);
    }
    Ok((hits, models))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hmmer::parse::parse_hmm_str;

    fn tiny() -> Profile {
        parse_hmm_str(include_str!("testdata/tiny.hmm")).unwrap()[0].configure()
    }

    #[test]
    fn cut_ga_requires_cutoffs() {
        let mut p = tiny();
        let t = Targets::new([("s", &b"MKA"[..])]);
        let cfg = SearchConfig { cutoff: Cutoff::Trusted, ..Default::default() };
        assert!(search_profiles(std::slice::from_ref(&p), &t, &cfg).is_err());
        p.tc = Some((1.0, 1.0));
        assert!(search_profiles(std::slice::from_ref(&p), &t, &cfg).is_ok());
    }

    #[test]
    fn bits_threshold_reports_motif_carrier() {
        let p = tiny();
        let t = Targets::new([("hit", &b"GGGGMKAGGGG"[..]), ("miss", &b"WWWWWWWWWWW"[..])]);
        let cfg = SearchConfig {
            cutoff: Cutoff::Bits { seq: 0.0, dom: 0.0 },
            pipeline: PipelineConfig { filters: false, ..Default::default() },
            z: None,
        };
        let hits = search_profiles(&[p], &t, &cfg).unwrap();
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].target, "hit");
        assert!(hits[0].domains[0].reported);
    }
}

#[cfg(test)]
mod hmmer34_parity {
    //! Regression against HMMER 3.4 `hmmsearch --max --domtblout` output for
    //! two simulated protein families (multi-domain, fragment, remote,
    //! degenerate-residue, low-complexity and random targets).
    use super::*;
    use crate::hmmer::parse::parse_hmm_str;

    #[test]
    fn matches_hmmsearch_scores_and_coordinates() {
        let profiles: Vec<Profile> =
            parse_hmm_str(include_str!("testdata/hmmer34_pair.hmm")).unwrap().iter().map(|h| h.configure()).collect();
        let fa = include_str!("testdata/hmmer34_targets.faa");
        let mut pairs: Vec<(String, Vec<u8>)> = Vec::new();
        for line in fa.lines() {
            if let Some(h) = line.strip_prefix('>') {
                pairs.push((h.trim().to_string(), Vec::new()));
            } else if let Some(last) = pairs.last_mut() {
                last.1.extend(line.trim().bytes());
            }
        }
        let t = Targets::from_pairs(&pairs);
        let cfg = SearchConfig {
            cutoff: Cutoff::Evalue { seq: 10.0, dom: 10.0 },
            pipeline: PipelineConfig { filters: false, ..Default::default() },
            z: None,
        };
        let hits = search_profiles(&profiles, &t, &cfg).unwrap();
        let expected = include_str!("testdata/hmmer34_expected.tsv");
        let mut n = 0;
        for row in expected.lines() {
            let f: Vec<&str> = row.split('\t').collect();
            let (tgt, q) = (f[0], f[1]);
            let num = |i: usize| f[i].parse::<f32>().unwrap();
            let h = hits.iter().find(|h| h.target == tgt && h.profile == q).unwrap_or_else(|| panic!("missing hit {tgt} {q}"));
            assert!((h.score - num(2)).abs() <= 0.1, "{tgt} {q}: seq {} vs {}", h.score, num(2));
            assert!((h.bias - num(3)).abs() <= 0.1, "{tgt} {q}: bias {} vs {}", h.bias, num(3));
            let d: Vec<&DomainHit> = h.reported_domains().collect();
            let di = f[4].parse::<usize>().unwrap() - 1;
            let dd = d[di];
            assert!((dd.score - num(5)).abs() <= 0.1, "{tgt} {q} dom{}: {} vs {}", di + 1, dd.score, num(5));
            let coords = [dd.hmm_from, dd.hmm_to, dd.ali_from, dd.ali_to, dd.env_from, dd.env_to];
            let want: Vec<usize> = (7..13).map(|i| f[i].parse().unwrap()).collect();
            assert_eq!(coords.to_vec(), want, "{tgt} {q} dom{} coords", di + 1);
            n += 1;
        }
        assert_eq!(n, 20);
    }
}
