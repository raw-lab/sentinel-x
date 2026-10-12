//! # BUSCO (protein-level), ported from BUSCO 6.1.0
//!
//! Each BUSCO's profile (`hmms/<id>.hmm`) is searched against the predicted
//! proteins (hmmsearch defaults, `--domtblout`), and the hits are processed
//! as `busco_tools/hmmer.py` does:
//!
//! * targets ordered by HMM-matched length (union of the domains' hmm
//!   coordinates, BUSCO's `sum_hmm_len`) then score, both descending;
//! * in genome mode a gene overlapping an already accepted gene on the same
//!   contig is skipped; hits below `scores_cutoff` (full-sequence score) are
//!   dropped;
//! * length class per dataset version — odb10: ζ = (mean − matched)/σ,
//!   complete for |ζ| ≤ 2, very large for ζ < −2; odb12.2+ (`lengths_cutoff`
//!   = median, scaled MAD): protein length within [min(0.9·med, med −
//!   min(25, 2·sMAD), 0.8·M), 1.5·med]; odb12 (no `lengths_cutoff`):
//!   complete when matched length ≥ 0.8·M;
//! * a gene used by a higher class (complete > very large > fragment) is
//!   removed from lower classes, a gene matching several BUSCOs in one class
//!   keeps only its best BUSCO, matches below 0.85× a BUSCO's best are
//!   dropped, one gene = single copy, more = duplicated, and only the best
//!   fragment is kept;
//! * `C:{}%[S:{}%,D:{}%],F:{}%,M:{}%,n:{}` with Python `round(x, 1)`.
//!
//! Gene prediction differs from BUSCO's (DYNAMITE instead of Prodigal for
//! prokaryotes, MetaEuk/miniprot for eukaryotes), so counts can differ where
//! the gene calls do.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

use crate::common::pyfmt::{py_float, py_round, reparse};
use crate::common::{Domain, QcRecord, Quality};
use crate::hmmer::{search_hmm_file_models, Cutoff, SearchConfig, Targets};

/// Per-BUSCO cutoffs (`scores_cutoff`, `lengths_cutoff`).
#[derive(Clone, Copy, Debug, Default)]
pub struct BuscoCutoff {
    pub score: Option<f64>,
    /// odb10: mean length; odb12.2+: median length.
    pub length: Option<f64>,
    /// odb10: σ (0 → 1); odb12.2+: scaled MAD.
    pub sigma: Option<f64>,
}

/// A BUSCO lineage dataset directory.
#[derive(Clone, Debug)]
pub struct Lineage {
    pub dir: PathBuf,
    pub name: String,
    /// `odb10`, `odb12`, `odb12.2`, …
    pub version: String,
    /// `prokaryota`, `eukaryota` or `viruses`.
    pub domain: String,
    pub creation_date: Option<String>,
    pub has_lengths_file: bool,
    /// BUSCO id → cutoffs (BUSCO's `cutoff_dict`).
    pub cutoffs: BTreeMap<String, BuscoCutoff>,
}

impl Lineage {
    /// Load `dataset.cfg`, `lengths_cutoff` and `scores_cutoff`.
    pub fn load(dir: &Path) -> Result<Lineage> {
        let cfg_path = dir.join("dataset.cfg");
        let cfg = std::fs::read_to_string(&cfg_path).with_context(|| format!("reading {}", cfg_path.display()))?;
        let mut kv: HashMap<String, String> = HashMap::new();
        for l in cfg.lines() {
            if let Some((k, v)) = l.split_once('=') {
                kv.insert(k.trim().to_string(), v.trim().to_string());
            }
        }
        let name =
            kv.get("name").cloned().unwrap_or_else(|| dir.file_name().and_then(|s| s.to_str()).unwrap_or("lineage").to_string());
        let version = match name.rsplit_once("_odb") {
            Some((_, v)) => format!("odb{v}"),
            None => bail!("{}: cannot read the dataset version from name '{name}'", cfg_path.display()),
        };
        let domain = kv.get("domain").cloned().unwrap_or_default();
        let mut cutoffs: BTreeMap<String, BuscoCutoff> = BTreeMap::new();
        let lpath = dir.join("lengths_cutoff");
        let has_lengths_file = lpath.is_file();
        if has_lengths_file {
            let txt = std::fs::read_to_string(&lpath)?;
            for l in txt.lines() {
                let f: Vec<&str> = l.split_whitespace().collect();
                if f.is_empty() {
                    continue;
                }
                let num = |i: usize| -> Result<f64> {
                    f.get(i)
                        .with_context(|| format!("{}: short line '{l}'", lpath.display()))?
                        .parse::<f64>()
                        .with_context(|| format!("{}: bad number in '{l}'", lpath.display()))
                };
                let e = cutoffs.entry(f[0].to_string()).or_default();
                if version == "odb10" {
                    let sd = num(2)?;
                    e.sigma = Some(if sd == 0.0 { 1.0 } else { sd });
                    e.length = Some(num(3)?);
                } else {
                    e.length = Some(num(1)?);
                    e.sigma = Some(num(2)?);
                }
            }
        }
        let spath = dir.join("scores_cutoff");
        let txt = std::fs::read_to_string(&spath).with_context(|| format!("reading {}", spath.display()))?;
        for l in txt.lines() {
            let f: Vec<&str> = l.split_whitespace().collect();
            if f.is_empty() {
                continue;
            }
            let s: f64 = f
                .get(1)
                .with_context(|| format!("{}: short line '{l}'", spath.display()))?
                .parse()
                .with_context(|| format!("{}: bad score in '{l}'", spath.display()))?;
            cutoffs.entry(f[0].to_string()).or_default().score = Some(s);
        }
        Ok(Lineage {
            dir: dir.to_path_buf(),
            name,
            version,
            domain,
            creation_date: kv.get("creation_date").cloned(),
            has_lengths_file,
            cutoffs,
        })
    }

    /// Number of BUSCOs (`n`).
    pub fn total(&self) -> usize {
        self.cutoffs.len()
    }

    /// Resolve `--lineage`: a dataset directory, or a name under
    /// `<root>/busco/` (searched up to `lineages/<name>`).
    pub fn resolve(spec: &str, busco_root: Option<&Path>) -> Result<PathBuf> {
        let p = PathBuf::from(spec);
        if p.join("dataset.cfg").is_file() {
            return Ok(p);
        }
        if let Some(root) = busco_root {
            for cand in [root.join(spec), root.join("lineages").join(spec)] {
                if cand.join("dataset.cfg").is_file() {
                    return Ok(cand);
                }
            }
        }
        bail!("BUSCO lineage '{spec}' not found (give a dataset directory or run `sentinel download busco --lineage {spec}`)")
    }
}

/// A predicted gene as BUSCO sees it in genome mode.
#[derive(Clone, Debug)]
pub struct BuscoGene {
    /// `contig:start-end|strand` (BUSCO's protein id).
    pub target: String,
    pub contig: String,
    pub start: usize,
    pub end: usize,
    /// Protein, with Prodigal's terminal `*`.
    pub protein: Vec<u8>,
}

impl BuscoGene {
    pub fn new(contig: &str, start: usize, end: usize, strand: i8, protein: Vec<u8>) -> BuscoGene {
        BuscoGene {
            target: format!("{}:{}-{}|{}", contig, start, end, if strand >= 0 { '+' } else { '-' }),
            contig: contig.to_string(),
            start,
            end,
            protein,
        }
    }
    fn gene_id(&self) -> String {
        format!("{}:{}-{}", self.contig, self.start, self.end)
    }
}

#[derive(Clone, Debug)]
struct Match {
    bitscore: f64,
    length: f64,
}

type ClassDict = BTreeMap<String, Vec<(String, Vec<Match>)>>; // busco -> [(gene, matches)] in insertion order

/// BUSCO's result.
#[derive(Clone, Debug, Default)]
pub struct BuscoSummary {
    pub lineage: String,
    pub single: usize,
    pub duplicated: usize,
    pub fragmented: usize,
    pub missing: usize,
    pub total: usize,
    /// (BUSCO id, status, gene, score, length) rows of `full_table.tsv`.
    pub table: Vec<(String, &'static str, String, f64, f64)>,
}

impl BuscoSummary {
    fn pct(&self, n: usize) -> f64 {
        if self.total == 0 {
            0.0
        } else {
            py_round(n as f64 / self.total as f64 * 100.0, 1).abs()
        }
    }
    pub fn complete_pct(&self) -> f64 {
        self.pct(self.single + self.duplicated)
    }
    pub fn single_pct(&self) -> f64 {
        self.pct(self.single)
    }
    pub fn duplicated_pct(&self) -> f64 {
        self.pct(self.duplicated)
    }
    pub fn fragmented_pct(&self) -> f64 {
        self.pct(self.fragmented)
    }
    pub fn missing_pct(&self) -> f64 {
        self.pct(self.missing)
    }
    /// `C:…%[S:…%,D:…%],F:…%,M:…%,n:…`
    pub fn one_line(&self) -> String {
        format!(
            "C:{}%[S:{}%,D:{}%],F:{}%,M:{}%,n:{}",
            py_float(self.complete_pct()),
            py_float(self.single_pct()),
            py_float(self.duplicated_pct()),
            py_float(self.fragmented_pct()),
            py_float(self.missing_pct()),
            self.total
        )
    }
}

/// `sum_hmm_len`: BUSCO's merge of hmm coordinate ranges.
fn sum_hmm_len(coords: &[(usize, usize)]) -> usize {
    let mut v = coords.to_vec();
    v.sort_by_key(|c| c.0);
    let mut used: Vec<[usize; 2]> = Vec::new();
    for r in v {
        let mut merged = false;
        for u in used.iter_mut() {
            if u[0] <= r.0 && r.0 <= u[1] {
                if r.1 > u[1] {
                    u[1] = r.1;
                }
                merged = true;
                break;
            }
        }
        if !merged {
            used.push([r.0, r.1]);
        }
    }
    used.iter().map(|u| u[1] - u[0] + 1).sum()
}

/// `_check_overlap` (same contig only).
fn overlaps(accepted: &[(String, usize, usize)], contig: &str, s2: usize, e2: usize) -> bool {
    for (c1, s1, e1) in accepted {
        if c1 != contig {
            continue;
        }
        let (mut a1, mut b1, mut a2, mut b2) = (*s1 as i64, *e1 as i64, s2 as i64, e2 as i64);
        if b2 - a2 > b1 - a1 {
            std::mem::swap(&mut a1, &mut a2);
            std::mem::swap(&mut b1, &mut b2);
        }
        if (a1 <= a2 && a2 <= b1) || (a1 <= b2 && b2 <= b1) {
            return true;
        }
    }
    false
}

/// Run BUSCO's protein search and classification.
pub fn run(lineage: &Lineage, genes: &[BuscoGene], genome_mode: bool) -> Result<BuscoSummary> {
    let targets = Targets::new(genes.iter().map(|g| (g.target.as_str(), g.protein.as_slice())));
    let cfg = SearchConfig { cutoff: Cutoff::Evalue { seq: 10.0, dom: 10.0 }, ..Default::default() };
    let by_target: HashMap<&str, &BuscoGene> = genes.iter().map(|g| (g.target.as_str(), g)).collect();

    let mut complete: ClassDict = BTreeMap::new();
    let mut vlarge: ClassDict = BTreeMap::new();
    let mut fragment: ClassDict = BTreeMap::new();
    // gene -> buscos, per class (insertion order)
    let mut mg: [Vec<(String, Vec<String>)>; 3] = [Vec::new(), Vec::new(), Vec::new()];
    let push_mg = |mg: &mut Vec<(String, Vec<String>)>, gene: &str, busco: &str| match mg.iter_mut().find(|x| x.0 == gene) {
        Some(x) => x.1.push(busco.to_string()),
        None => mg.push((gene.to_string(), vec![busco.to_string()])),
    };

    // BUSCO processes its result files in sorted order
    for (busco, cut) in &lineage.cutoffs {
        let hmm = lineage.dir.join("hmms").join(format!("{busco}.hmm"));
        if !hmm.is_file() {
            continue;
        }
        let Some(score_cut) = cut.score else { continue };
        let (hits, models) = search_hmm_file_models(&hmm, &targets, &cfg)?;
        let profile_len = models.first().map(|m| m.m).unwrap_or(0) as f64;
        // per target: printed score, tlen, matched length (domtblout rows)
        let mut rows: Vec<(String, f64, usize, usize)> = Vec::new();
        for h in &hits {
            let coords: Vec<(usize, usize)> = h.reported_domains().map(|d| (d.hmm_from, d.hmm_to)).collect();
            if coords.is_empty() {
                continue; // not in --domtblout
            }
            rows.push((h.target.clone(), reparse(&format!("{:.1}", h.score)), h.tlen, sum_hmm_len(&coords)));
        }
        rows.sort_by(|a, b| a.0.cmp(&b.0)); // np.unique
        rows.sort_by(|a, b| b.3.cmp(&a.3).then(b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal)));
        let mut accepted: Vec<(String, usize, usize)> = Vec::new();
        let mut records: Vec<(String, f64, usize, usize)> = Vec::new(); // gene_id, score, tlen, matched
        for (target, score, tlen, matched) in rows {
            let g = by_target[target.as_str()];
            let gene_id = if genome_mode { g.gene_id() } else { target.clone() };
            if genome_mode && overlaps(&accepted, &g.contig, g.start, g.end) {
                continue;
            }
            if score < score_cut {
                continue;
            }
            records.push((gene_id, score, tlen, matched));
            accepted.push((g.contig.clone(), g.start, g.end));
        }
        // _sort_matches
        for (gene_id, score, tlen, matched) in records {
            let (class, size) = if lineage.version == "odb10" {
                let (Some(len), Some(sigma)) = (cut.length, cut.sigma) else { continue };
                let size = matched as f64;
                let zeta = (len - size) / sigma;
                (
                    if (-2.0..=2.0).contains(&zeta) {
                        0
                    } else if zeta < -2.0 {
                        1
                    } else {
                        2
                    },
                    size,
                )
            } else if lineage.has_lengths_file {
                let (Some(med), Some(smad)) = (cut.length, cut.sigma) else { continue };
                let margin = 25.0f64.min(2.0 * smad);
                let size = tlen as f64;
                let lower = (0.9 * med).min(med - margin).min(0.8 * profile_len);
                let upper = 1.5 * med;
                (
                    if lower <= size && size <= upper {
                        0
                    } else if size > upper {
                        1
                    } else {
                        2
                    },
                    size,
                )
            } else {
                let size = matched as f64;
                (if size >= 0.8 * profile_len { 0 } else { 2 }, size)
            };
            let d = match class {
                0 => &mut complete,
                1 => &mut vlarge,
                _ => &mut fragment,
            };
            let e = d.entry(busco.clone()).or_default();
            match e.iter_mut().find(|x| x.0 == gene_id) {
                Some(x) => x.1.push(Match { bitscore: score, length: size }),
                None => e.push((gene_id.clone(), vec![Match { bitscore: score, length: size }])),
            }
            push_mg(&mut mg[class], &gene_id, busco);
        }
    }

    // _remove_duplicates
    let mut used: HashSet<String> = complete.values().flat_map(|v| v.iter().map(|x| x.0.clone())).collect();
    remove_lower_ranked(&mut vlarge, &mut mg[1], &complete.keys().cloned().collect(), &used);
    used.extend(vlarge.values().flat_map(|v| v.iter().map(|x| x.0.clone())));
    let higher: HashSet<String> = complete.keys().chain(vlarge.keys()).cloned().collect();
    remove_lower_ranked(&mut fragment, &mut mg[2], &higher, &used);
    let [mg0, mg1, mg2] = &mut mg;
    remove_remaining_duplicates(&mut complete, mg0);
    remove_remaining_duplicates(&mut vlarge, mg1);
    remove_remaining_duplicates(&mut fragment, mg2);
    // filter: 85 % of the best per BUSCO
    for d in [&mut complete, &mut vlarge, &mut fragment] {
        remove_low_scoring(d);
    }

    // consolidate
    let mut s = BuscoSummary { lineage: lineage.name.clone(), total: lineage.total(), ..Default::default() };
    let mut single_ids: BTreeMap<String, Vec<(String, Vec<Match>)>> = BTreeMap::new();
    let mut multi_ids: BTreeMap<String, Vec<(String, Vec<Match>)>> = BTreeMap::new();
    for d in [&complete, &vlarge] {
        for (b, gm) in d {
            if gm.len() == 1 {
                single_ids.insert(b.clone(), gm.clone());
            } else {
                multi_ids.insert(b.clone(), gm.clone());
            }
        }
    }
    let mut frag_ids: BTreeMap<String, Vec<(String, Vec<Match>)>> = BTreeMap::new();
    for (b, gm) in &fragment {
        if gm.len() > 1 {
            let best = best_match(gm).map(|x| x.0).unwrap_or(0);
            frag_ids.insert(b.clone(), vec![gm[best].clone()]);
        } else {
            frag_ids.insert(b.clone(), gm.clone());
        }
    }
    s.single = single_ids.len();
    s.duplicated = multi_ids.len();
    s.fragmented = frag_ids.len();
    s.missing = s.total.saturating_sub(s.single + s.duplicated + s.fragmented);
    for (label, d) in [("Complete", &single_ids), ("Duplicated", &multi_ids), ("Fragmented", &frag_ids)] {
        for (b, gm) in d {
            for (g, ms) in gm {
                for m in ms {
                    s.table.push((b.clone(), label, g.clone(), m.bitscore, m.length));
                }
            }
        }
    }
    for b in lineage.cutoffs.keys() {
        if !single_ids.contains_key(b) && !multi_ids.contains_key(b) && !frag_ids.contains_key(b) {
            s.table.push((b.clone(), "Missing", String::new(), f64::NAN, f64::NAN));
        }
    }
    s.table.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(s)
}

fn remove_from_mg(mg: &mut Vec<(String, Vec<String>)>, gene: &str, busco: &str) {
    if let Some(pos) = mg.iter().position(|x| x.0 == gene) {
        mg[pos].1.retain(|b| b != busco);
        if mg[pos].1.is_empty() {
            mg.remove(pos);
        }
    }
}

/// `_remove_lower_ranked_duplicates`.
fn remove_lower_ranked(d: &mut ClassDict, mg: &mut Vec<(String, Vec<String>)>, higher: &HashSet<String>, used: &HashSet<String>) {
    let keys: Vec<String> = d.keys().cloned().collect();
    for b in keys {
        if higher.contains(&b) {
            if let Some(matches) = d.remove(&b) {
                for (g, _) in matches {
                    remove_from_mg(mg, &g, &b);
                }
            }
            continue;
        }
        let genes: Vec<String> = d[&b].iter().map(|x| x.0.clone()).collect();
        for g in genes {
            if used.contains(&g) {
                if let Some(v) = d.get_mut(&b) {
                    v.retain(|x| x.0 != g);
                    if v.is_empty() {
                        d.remove(&b);
                    }
                }
                remove_from_mg(mg, &g, &b);
            }
        }
    }
}

/// `_remove_remaining_duplicate_matches`: a gene matched by several BUSCOs
/// within one class keeps only its best-scoring BUSCO.
fn remove_remaining_duplicates(d: &mut ClassDict, mg: &mut Vec<(String, Vec<String>)>) {
    let mut to_remove: Vec<(String, String)> = Vec::new();
    for (gene, buscos) in mg.iter() {
        if buscos.len() <= 1 {
            continue;
        }
        let mut scores: Vec<f64> = Vec::new();
        let mut order: Vec<String> = Vec::new();
        for b in buscos {
            if order.contains(b) {
                continue;
            }
            if let Some(gm) = d.get(b).and_then(|v| v.iter().find(|x| &x.0 == gene)) {
                for m in &gm.1 {
                    scores.push(m.bitscore);
                }
            }
            order.push(b.clone());
        }
        let distinct: HashSet<&String> = buscos.iter().collect();
        if distinct.len() == 1 || scores.is_empty() {
            continue;
        }
        // first maximum, as max(range, key=...) returns
        let mut best = 0;
        for (i, &s) in scores.iter().enumerate() {
            if s > scores[best] {
                best = i;
            }
        }
        let keep = order.get(best).cloned().unwrap_or_default();
        let mut dups: Vec<&String> = buscos.iter().filter(|b| **b != keep).collect();
        dups.sort();
        dups.dedup();
        for dup in dups {
            if let Some(v) = d.get_mut(dup) {
                v.retain(|x| &x.0 != gene);
                if v.is_empty() {
                    d.remove(dup);
                }
            }
            to_remove.push((gene.clone(), dup.clone()));
        }
    }
    for (g, b) in to_remove {
        if let Some(pos) = mg.iter().position(|x| x.0 == g) {
            if let Some(i) = mg[pos].1.iter().position(|x| *x == b) {
                mg[pos].1.remove(i);
            }
            if mg[pos].1.is_empty() {
                mg.remove(pos);
            }
        }
    }
}

fn best_match(gm: &[(String, Vec<Match>)]) -> Option<(usize, f64)> {
    let mut best: Option<(usize, f64)> = None;
    for (i, (_, ms)) in gm.iter().enumerate() {
        for m in ms {
            if best.map(|b| m.bitscore > b.1).unwrap_or(true) {
                best = Some((i, m.bitscore));
            }
        }
    }
    best
}

/// `_remove_low_scoring_matches`: drop matches below 0.85× the BUSCO's best.
fn remove_low_scoring(d: &mut ClassDict) {
    for gm in d.values_mut() {
        if gm.len() <= 1 {
            continue;
        }
        let Some((_, max)) = best_match(gm) else { continue };
        for (_, ms) in gm.iter_mut() {
            ms.retain(|m| m.bitscore >= 0.85 * max);
        }
        gm.retain(|x| !x.1.is_empty());
    }
}

/// Unified-report record. BUSCO makes no contamination claim; the
/// duplicated fraction is reported in the extra column.
pub fn to_record(genome_id: &str, total_len: usize, domain: Domain, n_genes: usize, s: &BuscoSummary) -> QcRecord {
    let mut r = QcRecord::new(genome_id, total_len, domain, "busco");
    r.gene_count = n_genes;
    r.completeness = Some(s.complete_pct());
    r.contamination = None;
    r.method = format!("BUSCO {}", s.lineage);
    r.quality = Quality::from_completeness(Some(s.complete_pct()));
    r = r
        .with_extra("busco", s.one_line())
        .with_extra("single", s.single.to_string())
        .with_extra("duplicated", s.duplicated.to_string())
        .with_extra("duplicated_pct", py_float(s.duplicated_pct()))
        .with_extra("fragmented", s.fragmented.to_string())
        .with_extra("missing", s.missing.to_string())
        .with_extra("n", s.total.to_string());
    if s.duplicated_pct() >= 5.0 {
        r.warnings
            .push(format!("BUSCO duplicated {}% (possible contamination or recent duplication)", py_float(s.duplicated_pct())));
    }
    r
}

/// BUSCO-style `full_table.tsv` (Busco id, Status, Sequence, Score, Length).
pub fn full_table_tsv(s: &BuscoSummary) -> String {
    let mut out = String::from("# Busco id\tStatus\tSequence\tScore\tLength\n");
    for (b, st, g, sc, len) in &s.table {
        if *st == "Missing" {
            out.push_str(&format!("{b}\tMissing\n"));
        } else {
            out.push_str(&format!("{b}\t{st}\t{g}\t{}\t{}\n", py_float(*sc), *len as i64));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hmm_length_merge_quirk() {
        assert_eq!(sum_hmm_len(&[(1, 10), (5, 20), (30, 40)]), 20 + 11);
        // a region starting inside the first used region extends it
        assert_eq!(sum_hmm_len(&[(1, 10), (2, 5)]), 10);
    }

    #[test]
    fn overlap_rule() {
        let acc = vec![("c".to_string(), 100, 400)];
        assert!(overlaps(&acc, "c", 350, 900));
        assert!(!overlaps(&acc, "c", 401, 900));
        assert!(!overlaps(&acc, "d", 100, 400));
    }

    #[test]
    fn one_line_summary() {
        let s = BuscoSummary { single: 120, duplicated: 1, fragmented: 1, missing: 2, total: 124, ..Default::default() };
        assert_eq!(s.one_line(), "C:97.6%[S:96.8%,D:0.8%],F:0.8%,M:1.6%,n:124");
    }

    #[test]
    fn duplicate_resolution_keeps_best_busco() {
        let mut d: ClassDict = BTreeMap::new();
        d.insert("b1".into(), vec![("g".into(), vec![Match { bitscore: 50.0, length: 1.0 }])]);
        d.insert("b2".into(), vec![("g".into(), vec![Match { bitscore: 80.0, length: 1.0 }])]);
        let mut mg = vec![("g".to_string(), vec!["b1".to_string(), "b2".to_string()])];
        remove_remaining_duplicates(&mut d, &mut mg);
        assert!(!d.contains_key("b1"));
        assert!(d.contains_key("b2"));
        assert_eq!(mg[0].1, vec!["b2".to_string()]);
    }
}
