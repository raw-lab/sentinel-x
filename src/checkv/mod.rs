//! # CheckV 1.1.1, ported
//!
//! Viral genome quality following CheckV's `end_to_end` pipeline step by
//! step, with Sentinel's own engines in place of the external tools:
//!
//! 1. **Gene annotation** ([`annotate`]) — every protein is searched against
//!    each `hmm_db/checkv_hmms/*.hmm` file with Sentinel's HMMER 3.4 engine
//!    (`hmmsearch -E 10`); a gene's best hit at `max(score_cutoff, 25)` bits
//!    labels it viral (+1) or microbial (−1).
//! 2. **Contamination** ([`regions`]) — host/virus boundaries from CheckV's
//!    recursive breakpoint scan; contamination = host bp / contig bp.
//! 3. **Completeness** — AAI-based (DIAMOND-compatible search of the
//!    reference proteins with Sentinel's `protsearch`, CheckV's AAI and
//!    50%-of-top-hit weighting, error-table confidence) and HMM-based
//!    (genome-length quantiles of the hit HMMs).
//! 4. **Complete genomes** ([`repeats`]) — DTR/ITR detection with CheckV's
//!    N/count/low-complexity filters and the 21-mer copy-number check;
//!    complete proviruses.
//! 5. **Quality summary** — CheckV's tiers, MIUViG quality and warnings.
//!
//! [`report`] writes CheckV's output files byte-compatibly.

pub mod db;
pub mod regions;
pub mod repeats;
pub mod report;

use std::collections::{HashMap, HashSet};

use anyhow::Result;
use rayon::prelude::*;

pub use db::CheckvDb;
use regions::{GeneView, Region, RegionParams};
use repeats::TerminalRepeat;

use crate::common::pyfmt::{c_fmt_g, np_quantile, np_round, py_float, py_fsum, py_round, py_sum, reparse, PyNum};
use crate::common::{Contig, Domain, GeneSet, QcRecord, Quality};

/// A predicted protein-coding gene on a contig.
#[derive(Clone, Debug)]
pub struct ViralGene {
    /// `<contig>_<n>` (1-based n in coordinate order).
    pub id: String,
    pub start: usize,
    pub end: usize,
    pub strand: i8,
    /// Protein (ASCII) including the terminal `*` of genes that end in a
    /// stop codon, as pyrodigal-gv writes it for CheckV.
    pub protein: Vec<u8>,
    /// Prodigal `partial=` flags: gene runs off the left / right contig edge.
    pub partial: [bool; 2],
}

/// One input contig with its genes.
#[derive(Clone, Debug)]
pub struct ContigInput {
    pub id: String,
    /// Uppercase sequence.
    pub seq: Vec<u8>,
    pub genes: Vec<ViralGene>,
}

/// CheckV options (defaults as `checkv end_to_end`).
#[derive(Clone, Debug)]
pub struct Params {
    pub regions: RegionParams,
    pub percent_of_top_hit: f64,
    pub max_aai: Option<f64>,
    pub exclude_identical: bool,
    pub tr_min_len: usize,
    pub tr_max_count: usize,
    pub tr_max_ambig: f64,
    /// 0.70 under `end_to_end` (0.75 for `checkv complete_genomes`).
    pub tr_max_basefreq: f64,
    pub kmer_max_freq: f64,
    /// Group DIAMOND hits the way CheckV 1.1.1 actually does: consecutive
    /// lines of DIAMOND's output per contig. DIAMOND 2.x prints queries
    /// sorted by length when query and subject coverage are equal (as
    /// CheckV sets them), so a contig's genes are split into many groups
    /// and AAI is computed per fragment. Off by default (all of a contig's
    /// hits form one group, as CheckV intends); on for byte-level parity.
    pub strict_diamond_order: bool,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            regions: RegionParams::default(),
            percent_of_top_hit: 50.0,
            max_aai: None,
            exclude_identical: false,
            tr_min_len: 21,
            tr_max_count: 8,
            tr_max_ambig: 0.2,
            tr_max_basefreq: 0.70,
            kmer_max_freq: 1.5,
            strict_diamond_order: false,
        }
    }
}

/// A gene's best HMM hit.
#[derive(Clone, Debug)]
pub struct HmmHit {
    pub hmm: String,
    pub database: String,
    pub evalue: f64,
    pub score: f32,
}

impl HmmHit {
    /// Full-sequence E-value as `--tblout` prints it (`%9.2g`) and CheckV
    /// reads it back.
    pub fn printed_evalue(&self) -> f64 {
        reparse(&c_fmt_g(self.evalue, 2))
    }
    /// Full-sequence score as `--tblout` prints it (`%6.1f`).
    pub fn printed_score(&self) -> f64 {
        reparse(&format!("{:.1}", self.score))
    }
}

/// Per-gene annotation.
#[derive(Clone, Debug)]
pub struct GeneAnnotation {
    pub gc: f64,
    pub cat: i8,
    pub hit: Option<HmmHit>,
}

/// One AAI record (`tmp/aai.tsv`).
#[derive(Clone, Debug)]
pub struct AaiRecord {
    pub target: String,
    pub aligned_genes: usize,
    pub percent_genes: f64,
    pub aligned_length: usize,
    pub percent_length: f64,
    pub identity: f64,
    pub score: f64,
    /// Position in CheckV's `tmp/aai.tsv` (rows are written group by group).
    pub order: usize,
}

/// Completeness estimates for a contig.
#[derive(Clone, Debug, Default)]
pub struct Completeness {
    pub viral_length: Option<usize>,
    pub aai_hits: Vec<AaiRecord>,
    pub aai_expected_length: Option<f64>,
    pub aai_completeness: Option<f64>,
    pub aai_confidence: Option<&'static str>,
    /// `Some(None)` = error rate NA.
    pub aai_error: Option<Option<f64>>,
    pub hmm_hits: usize,
    pub hmm_completeness: Option<(f64, f64)>,
    pub kmer_freq: f64,
}

/// Complete-genome call (`complete_genomes.tsv` row).
#[derive(Clone, Debug)]
pub struct CompleteGenome {
    pub prediction: String,
    pub confidence: &'static str,
    pub reason: String,
    pub tr: Option<TerminalRepeat>,
}

/// All CheckV results for one contig.
#[derive(Clone, Debug)]
pub struct ContigResult {
    pub id: String,
    pub length: usize,
    pub seq: Vec<u8>,
    pub genes: Vec<ViralGene>,
    pub annotations: Vec<GeneAnnotation>,
    /// `None` without a database (contamination not computed).
    pub viral_genes: Option<usize>,
    pub host_genes: Option<usize>,
    pub regions: Vec<Region>,
    pub completeness: Completeness,
    pub complete: Option<CompleteGenome>,
    pub quality: &'static str,
    pub miuvig: &'static str,
    pub completeness_value: Option<f64>,
    pub completeness_method: String,
    pub contamination: Option<f64>,
    pub warnings: Vec<String>,
    pub aai_table: Vec<AaiRecord>,
}

impl ContigResult {
    pub fn provirus(&self) -> bool {
        !self.regions.is_empty()
    }
    pub fn proviral_length(&self) -> Option<usize> {
        if self.regions.is_empty() {
            None
        } else {
            Some(self.regions.iter().filter(|r| r.viral).map(|r| r.length).sum())
        }
    }
    pub fn host_length(&self) -> Option<usize> {
        if self.regions.is_empty() {
            None
        } else {
            Some(self.regions.iter().filter(|r| !r.viral).map(|r| r.length).sum())
        }
    }
    pub fn region_types(&self) -> String {
        self.regions.iter().map(|r| r.type_label()).collect::<Vec<_>>().join(",")
    }
}

/// CheckV inputs from Sentinel's gene calls (pyrodigal-gv conventions:
/// `<contig>_<n>` in coordinate order, `*`-terminated proteins).
pub fn inputs_from_genes(contigs: &[Contig], gs: &GeneSet) -> Vec<ContigInput> {
    let per = gs.per_contig(contigs.len());
    contigs
        .iter()
        .zip(per)
        .map(|(c, idx)| ContigInput {
            id: c.id.clone(),
            seq: c.seq.clone(),
            genes: idx
                .into_iter()
                .map(|i| {
                    let g = &gs.genes[i];
                    ViralGene {
                        id: g.id(),
                        start: g.begin,
                        end: g.end,
                        strand: g.strand,
                        protein: g.protein_with_stop(),
                        partial: [g.partial_left, g.partial_right],
                    }
                })
                .collect(),
        })
        .collect()
}

/// CheckV inputs from a Prodigal-format protein FASTA (`>id # start # end #
/// strand # …`, ids `<contig>_<n>`), read the way CheckV reads
/// `tmp/proteins.faa` (genes kept in file order).
pub fn inputs_from_faa(contigs: &[Contig], faa: &std::path::Path) -> Result<Vec<ContigInput>> {
    let prots = crate::common::read_fasta(faa)?;
    let index: HashMap<&str, usize> = contigs.iter().enumerate().map(|(i, c)| (c.id.as_str(), i)).collect();
    let mut out: Vec<ContigInput> =
        contigs.iter().map(|c| ContigInput { id: c.id.clone(), seq: c.seq.clone(), genes: Vec::new() }).collect();
    for p in prots {
        let f: Vec<&str> = p.header.split_whitespace().collect();
        let num = |i: usize| -> Result<i64> {
            f.get(i)
                .and_then(|x| x.parse::<i64>().ok())
                .ok_or_else(|| anyhow::anyhow!("{}: header '{}' is not in Prodigal format", faa.display(), p.header))
        };
        let (start, end, strand) = (num(2)?, num(4)?, num(6)?);
        let genome = p.id.rsplit_once('_').map(|x| x.0).unwrap_or(&p.id);
        let Some(&ci) = index.get(genome) else {
            anyhow::bail!("{}: protein '{}' refers to unknown contig '{genome}'", faa.display(), p.id);
        };
        let partial = f
            .iter()
            .find_map(|x| x.strip_prefix("partial="))
            .map(|x| {
                let b = x.as_bytes();
                [b.first() == Some(&b'1'), b.get(1) == Some(&b'1')]
            })
            .unwrap_or([false, false]);
        if start < 1 || end < start {
            anyhow::bail!("{}: bad coordinates in '{}'", faa.display(), p.header);
        }
        out[ci].genes.push(ViralGene {
            id: p.id.clone(),
            start: start as usize,
            end: end as usize,
            strand: if strand < 0 { -1 } else { 1 },
            protein: p.seq,
            partial,
        });
    }
    Ok(out)
}

/// Unified-report records (one per contig).
pub fn to_records(results: &[ContigResult], domain: Domain, have_db: bool) -> Vec<QcRecord> {
    results
        .iter()
        .map(|r| {
            let mut q = QcRecord::new(&r.id, r.length, domain, "checkv");
            q.gene_count = r.genes.len();
            q.completeness = r.completeness_value;
            q.contamination = r.contamination;
            q.quality = match r.quality {
                "Complete" => Quality::Complete,
                "High-quality" => Quality::High,
                "Medium-quality" => Quality::Medium,
                "Low-quality" => Quality::Low,
                _ => Quality::NotDetermined,
            };
            q.method = r.completeness_method.clone();
            q.warnings = r.warnings.clone();
            let c = &r.completeness;
            let opt = |v: Option<f64>| v.map(py_float).unwrap_or_else(|| "NA".into());
            if have_db {
                q = q
                    .with_extra("provirus", if r.provirus() { "Yes" } else { "No" })
                    .with_extra("proviral_length", r.proviral_length().map(|x| x.to_string()).unwrap_or_else(|| "NA".into()))
                    .with_extra("viral_genes", r.viral_genes.unwrap_or(0).to_string())
                    .with_extra("host_genes", r.host_genes.unwrap_or(0).to_string());
            }
            q = q.with_extra("miuvig_quality", r.miuvig).with_extra("kmer_freq", report::kmer_freq_str(c.kmer_freq));
            if let Some(top) = c.aai_hits.first() {
                q = q
                    .with_extra("aai_completeness", opt(c.aai_completeness))
                    .with_extra("aai_confidence", c.aai_confidence.unwrap_or("NA"))
                    .with_extra("aai_expected_length", opt(c.aai_expected_length))
                    .with_extra("aai_top_hit", top.target.clone())
                    .with_extra("aai_id", py_float(top.identity))
                    .with_extra("aai_af", py_float(top.percent_length));
            }
            if let Some((lo, hi)) = c.hmm_completeness {
                q = q
                    .with_extra("hmm_completeness_lower", py_float(lo))
                    .with_extra("hmm_completeness_upper", py_float(hi))
                    .with_extra("hmm_num_hits", c.hmm_hits.to_string());
            }
            if let Some(cg) = &r.complete {
                q = q
                    .with_extra("complete_prediction", cg.prediction.clone())
                    .with_extra("complete_confidence", cg.confidence)
                    .with_extra("complete_reason", cg.reason.clone());
                if let Some(t) = &cg.tr {
                    q = q.with_extra("repeat_length", t.len().to_string()).with_extra("repeat_count", t.count.to_string());
                }
            }
            q
        })
        .collect()
}

/// `compute_gc`: percent G+C with every base (including N) counted.
fn gene_gc(seq: &[u8], start: usize, end: usize) -> f64 {
    let s = &seq[start.saturating_sub(1).min(seq.len())..end.min(seq.len())];
    if s.is_empty() {
        return 0.0;
    }
    let gc = s.iter().filter(|&&b| b == b'G' || b == b'C').count();
    100.0 * gc as f64 / s.len() as f64
}

/// Step 1: HMM annotation of every gene.
pub fn annotate(contigs: &[ContigInput], db: &CheckvDb) -> Result<Vec<Vec<Option<HmmHit>>>> {
    use crate::hmmer::{search_hmm_file, Cutoff, SearchConfig, Targets};
    let mut ids: Vec<&str> = Vec::new();
    let mut seqs: Vec<&[u8]> = Vec::new();
    let mut index: HashMap<&str, (usize, usize)> = HashMap::new();
    for (ci, c) in contigs.iter().enumerate() {
        for (gi, g) in c.genes.iter().enumerate() {
            ids.push(&g.id);
            seqs.push(&g.protein);
            index.insert(&g.id, (ci, gi));
        }
    }
    let targets = Targets::new(ids.iter().copied().zip(seqs.iter().copied()));
    let cfg = SearchConfig { cutoff: Cutoff::Evalue { seq: 10.0, dom: 10.0 }, ..Default::default() };
    let mut best: Vec<Vec<Option<HmmHit>>> = contigs.iter().map(|c| vec![None; c.genes.len()]).collect();
    for f in &db.hmm_files {
        let (hits, _) = search_hmm_file(f, &targets, &cfg)?;
        for h in hits {
            // name if known, else accession (CheckV contamination.py:46)
            let key = if db.hmm_info.contains_key(&h.profile) {
                Some(h.profile.clone())
            } else {
                h.acc.clone().filter(|a| db.hmm_info.contains_key(a))
            };
            let Some(key) = key else { continue };
            let info = &db.hmm_info[&key];
            // CheckV compares the score as `--tblout` prints it (`%6.1f`)
            let printed = reparse(&format!("{:.1}", h.score));
            if printed < info.cutoff {
                continue;
            }
            let Some(&(ci, gi)) = index.get(h.target.as_str()) else { continue };
            let slot = &mut best[ci][gi];
            if slot.as_ref().map(|b| printed > b.printed_score()).unwrap_or(true) {
                *slot = Some(HmmHit { hmm: key, database: info.database.clone(), evalue: h.evalue, score: h.score });
            }
        }
    }
    Ok(best)
}

/// Run the full CheckV pipeline. Without a database only the
/// database-free parts (terminal repeats, k-mer copy number) are computed
/// and every quality call is `Not-determined`.
pub fn run(contigs: &[ContigInput], db: Option<&CheckvDb>, p: &Params) -> Result<Vec<ContigResult>> {
    // ---- 1. annotation
    let hits = match db {
        Some(d) => Some(annotate(contigs, d)?),
        None => None,
    };
    let mut results: Vec<ContigResult> = contigs
        .par_iter()
        .enumerate()
        .map(|(ci, c)| {
            let annotations: Vec<GeneAnnotation> = c
                .genes
                .iter()
                .enumerate()
                .map(|(gi, g)| {
                    let hit = hits.as_ref().and_then(|h| h[ci][gi].clone());
                    let cat = match (&hit, db) {
                        (Some(hh), Some(d)) => {
                            if d.hmm_info[&hh.hmm].viral {
                                1
                            } else {
                                -1
                            }
                        }
                        _ => 0,
                    };
                    GeneAnnotation { gc: gene_gc(&c.seq, g.start, g.end), cat, hit }
                })
                .collect();
            // ---- 2. host regions
            let regions = if db.is_some() {
                let views: Vec<GeneView> = c
                    .genes
                    .iter()
                    .zip(&annotations)
                    .map(|(g, a)| GeneView { start: g.start, end: g.end, gc: a.gc, cat: a.cat })
                    .collect();
                regions::recursive_regions(&views, c.seq.len(), &p.regions)
            } else {
                Vec::new()
            };
            let (viral_genes, host_genes) = if db.is_some() {
                (
                    Some(annotations.iter().filter(|a| a.cat == 1).count()),
                    Some(annotations.iter().filter(|a| a.cat == -1).count()),
                )
            } else {
                (None, None)
            };
            ContigResult {
                id: c.id.clone(),
                length: c.seq.len(),
                seq: c.seq.clone(),
                genes: c.genes.clone(),
                annotations,
                viral_genes,
                host_genes,
                regions,
                completeness: Completeness { kmer_freq: repeats::average_kmer_freq(&c.seq, 21), ..Default::default() },
                complete: None,
                quality: "Not-determined",
                miuvig: "Genome-fragment",
                completeness_value: None,
                completeness_method: "NA".into(),
                contamination: None,
                warnings: Vec::new(),
                aai_table: Vec::new(),
            }
        })
        .collect();

    // ---- 3. completeness
    if let Some(d) = db {
        for r in results.iter_mut() {
            // viral_length: proviral length when the region list has a viral region
            r.completeness.viral_length = if r.regions.iter().any(|x| x.viral) { r.proviral_length() } else { None };
        }
        aai_completeness(&mut results, d, p)?;
        hmm_completeness(&mut results, d);
    }

    // ---- 4. complete genomes, 5. quality summary
    for r in results.iter_mut() {
        complete_genome(r, p);
        quality_summary(r, db.is_some());
    }
    Ok(results)
}

/// CheckV `compute_aai` for one group of DIAMOND lines of a contig: per
/// reference genome, the best line per gene (strictly higher printed bit
/// score replaces), AAI weighted by protein length, rows sorted by score
/// (stable).
fn aai_rows(r: &ContigResult, hs: &[&crate::protsearch::Hit], owner: &[(usize, usize)]) -> Vec<AaiRecord> {
    let protlen_total: usize = r.genes.iter().map(|g| g.protein.len()).sum();
    let ngenes = r.genes.len();
    let mut targets: Vec<String> = Vec::new();
    let mut best: HashMap<String, Vec<(usize, &crate::protsearch::Hit)>> = HashMap::new();
    for h in hs {
        let gi = owner[h.query].1;
        let tgt = match h.subject.rsplit_once('_') {
            Some((a, _)) => a.to_string(),
            None => h.subject.clone(),
        };
        let entry = best.entry(tgt.clone()).or_insert_with(|| {
            targets.push(tgt.clone());
            Vec::new()
        });
        // the first line is kept unless a later one has a strictly higher
        // bit score as DIAMOND prints it
        match entry.iter_mut().find(|(g, _)| *g == gi) {
            None => entry.push((gi, h)),
            Some(slot) => {
                if h.printed_bitscore() > slot.1.printed_bitscore() {
                    *slot = (gi, h);
                }
            }
        }
    }
    let mut aai: Vec<AaiRecord> = Vec::new();
    for t in &targets {
        let alns = &best[t];
        let lens: Vec<usize> = alns.iter().map(|(g, _)| r.genes[*g].protein.len()).collect();
        let aligned_length: usize = lens.iter().sum();
        if aligned_length == 0 {
            continue;
        }
        // pident as parsed from DIAMOND's table (one decimal, `llround`);
        // Python's sum() of floats is compensated (CPython >= 3.12)
        let num: f64 = py_fsum(alns.iter().zip(&lens).map(|((_, h), &l)| h.printed_pident() * l as f64));
        let identity = py_round(num / aligned_length as f64, 2);
        aai.push(AaiRecord {
            target: t.clone(),
            aligned_genes: alns.len(),
            percent_genes: py_round(100.0 * alns.len() as f64 / ngenes.max(1) as f64, 2),
            aligned_length,
            percent_length: py_round(100.0 * aligned_length as f64 / protlen_total.max(1) as f64, 2),
            identity,
            score: py_round(identity * aligned_length as f64 / 100.0, 2),
            order: 0,
        });
    }
    aai.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal));
    aai
}

/// AAI-based completeness (CheckV `compute_aai`, `store_aai`,
/// `aai_based_completeness`).
fn aai_completeness(results: &mut [ContigResult], db: &CheckvDb, p: &Params) -> Result<()> {
    let mut queries: Vec<(String, Vec<u8>)> = Vec::new();
    let mut owner: Vec<(usize, usize)> = Vec::new();
    for (ci, r) in results.iter().enumerate() {
        for (gi, g) in r.genes.iter().enumerate() {
            queries.push((g.id.clone(), g.protein.clone()));
            owner.push((ci, gi));
        }
    }
    if queries.is_empty() {
        return Ok(());
    }
    // hits come back in DIAMOND's output order
    let (hits, _) = crate::protsearch::search_path(&queries, &db.reps_proteins, &crate::protsearch::Params::checkv())?;
    let mut groups: Vec<(usize, Vec<&crate::protsearch::Hit>)> = Vec::new();
    if p.strict_diamond_order {
        // yield_query_alns: consecutive lines of the same contig
        for h in &hits {
            let ci = owner[h.query].0;
            match groups.last_mut() {
                Some((c, v)) if *c == ci => v.push(h),
                _ => groups.push((ci, vec![h])),
            }
        }
    } else {
        let mut per: Vec<Vec<&crate::protsearch::Hit>> = results.iter().map(|_| Vec::new()).collect();
        for h in &hits {
            per[owner[h.query].0].push(h);
        }
        groups = per.into_iter().enumerate().filter(|(_, v)| !v.is_empty()).collect();
    }
    let mut order = 0usize;
    for (ci, hs) in &groups {
        let mut rows = aai_rows(&results[*ci], hs, &owner);
        for r in rows.iter_mut() {
            r.order = order;
            order += 1;
        }
        results[*ci].aai_table.extend(rows);
    }

    for r in results.iter_mut() {
        if r.aai_table.is_empty() {
            continue;
        }
        // store_aai: rows in file order
        let mut kept: Vec<AaiRecord> = Vec::new();
        for rec in &r.aai_table {
            if p.max_aai.map(|m| rec.identity > m).unwrap_or(false) {
                continue;
            }
            if !db.reps.contains_key(&rec.target) {
                continue;
            }
            if p.exclude_identical && rec.identity == 100.0 && rec.percent_length == 100.0 {
                continue;
            }
            match kept.first() {
                None => kept.push(rec.clone()),
                Some(top) => {
                    let ts = top.score;
                    if ts != 0.0 && 100.0 * (ts - rec.score) / ts <= p.percent_of_top_hit {
                        kept.push(rec.clone());
                    }
                }
            }
        }
        if kept.is_empty() {
            continue;
        }
        let (err, conf) = db.error.lookup(r.length, kept[0].identity, kept[0].percent_length);
        let mut nums: Vec<f64> = Vec::new();
        let mut dens: Vec<f64> = Vec::new();
        let mut ok = true;
        for k in &kept {
            let rg = &db.reps[&k.target];
            match rg.weight {
                Some(w) => {
                    nums.push(rg.length as f64 * (w * k.score));
                    dens.push(w * k.score);
                }
                None => ok = false,
            }
        }
        let den = py_fsum(dens);
        if !ok || den == 0.0 {
            continue;
        }
        let expected = py_fsum(nums) / den;
        let qlen = r.completeness.viral_length.unwrap_or(r.length) as f64;
        r.completeness.aai_expected_length = Some(expected);
        r.completeness.aai_completeness = Some((100.0 * qlen / expected).min(100.0));
        r.completeness.aai_confidence = Some(conf);
        r.completeness.aai_error = Some(err);
        r.completeness.aai_hits = kept;
    }
    Ok(())
}

/// HMM-based completeness (CheckV `hmm_based_completeness`).
fn hmm_completeness(results: &mut [ContigResult], db: &CheckvDb) {
    for r in results.iter_mut() {
        let mut hmms: Vec<&str> = Vec::new();
        let mut seen: HashSet<&str> = HashSet::new();
        for a in &r.annotations {
            if let Some(h) = &a.hit {
                if let Some(gl) = db.genome_lengths.get(&h.hmm) {
                    if gl.genomes >= 10 && seen.insert(&h.hmm) {
                        hmms.push(&h.hmm);
                    }
                }
            }
        }
        r.completeness.hmm_hits = hmms.len();
        if hmms.is_empty() {
            continue;
        }
        let qlen = r.completeness.viral_length.unwrap_or(r.length) as f64;
        let mut comps: Vec<(f64, f64, f64)> = Vec::new();
        let mut weights: Vec<PyNum> = Vec::new();
        for h in &hmms {
            let gl = &db.genome_lengths[*h];
            let q1 = np_quantile(&gl.lengths, 0.05);
            let q2 = np_quantile(&gl.lengths, 0.95);
            // `round()` of a numpy.float64 is NumPy's rint(x·100)/100
            let c1 = np_round((100.0 * qlen / q2).min(100.0), 2);
            let c2 = np_round((100.0 * qlen / q1).min(100.0), 2);
            // min([1 / cv, 50]) keeps 1/cv (a float) unless the int 50 is smaller
            let w = if gl.cv != 0.0 && 1.0 / gl.cv <= 50.0 {
                weights.push(PyNum::F(1.0 / gl.cv));
                1.0 / gl.cv
            } else {
                weights.push(PyNum::I(50));
                50.0
            };
            comps.push((w, c1, c2));
        }
        // CheckV iterates a Python set (hash-randomised order), so its own
        // results can differ from run to run in the last bit; the weight
        // total is a compensated Python sum, the products (NumPy scalars)
        // are summed plainly
        let wt: f64 = py_sum(&weights);
        let lo: f64 = comps.iter().map(|c| c.0 * c.1).sum::<f64>() / wt;
        let hi: f64 = comps.iter().map(|c| c.0 * c.2).sum::<f64>() / wt;
        r.completeness.hmm_completeness = Some((lo, hi));
    }
}

/// `complete_genomes`: terminal repeats, flags and confidence; complete
/// proviruses.
fn complete_genome(r: &mut ContigResult, p: &Params) {
    let tr = repeats::detect(&r.seq, p.tr_min_len);
    let mut flagged = false;
    let mut reason = String::new();
    if let Some(t) = &tr {
        if let Some(why) = repeats::flag_reason(t, p.tr_max_ambig, p.tr_max_count, p.tr_max_basefreq) {
            flagged = true;
            reason = why.to_string();
        }
    }
    // completeness.tsv carries kmer_freq rounded to 2 dp
    let kmer = np_round(r.completeness.kmer_freq, 2);
    if kmer > p.kmer_max_freq {
        flagged = true;
        reason = "Multiple genome copies detected".into();
    }
    let complete_provirus = r.region_types() == "host,viral,host";
    let prediction = match &tr {
        Some(t) => t.kind.label().to_string(),
        None if complete_provirus => "Provirus".to_string(),
        None => return,
    };
    let aai = match r.completeness.aai_confidence {
        Some("medium") | Some("high") => r.completeness.aai_completeness,
        _ => None,
    };
    // HMM completeness only when hmm_num_hits > 0
    let hmm = if r.completeness.hmm_hits > 0 { r.completeness.hmm_completeness } else { None };
    let (confidence, why): (&'static str, String) = if tr.is_some() && flagged {
        ("low", reason)
    } else if aai.map(|a| a != 0.0 && a > 90.0).unwrap_or(false) {
        ("high", "AAI-based completeness > 90%".into())
    } else if aai.map(|a| a != 0.0 && a > 80.0).unwrap_or(false) {
        ("medium", "AAI-based completeness > 80%".into())
    } else if aai.map(|a| a != 0.0 && a < 80.0).unwrap_or(false) {
        ("low", "AAI-based completeness < 80%".into())
    } else if hmm.map(|h| h.0 > 90.0).unwrap_or(false) {
        ("high", "HMM-based completeness (lower-bound) > 90%".into())
    } else if hmm.map(|h| h.0 > 80.0).unwrap_or(false) {
        ("medium", "HMM-based completeness (lower-bound) > 80%".into())
    } else if hmm.map(|h| h.1 < 80.0).unwrap_or(false) {
        ("low", "HMM-based completeness (upper-bound) < 80%".into())
    } else {
        ("not-determined", "NA".into())
    };
    r.complete = Some(CompleteGenome { prediction, confidence, reason: why, tr });
}

/// `quality_summary`: completeness choice, tiers and warnings.
fn quality_summary(r: &mut ContigResult, have_db: bool) {
    let mut warnings = Vec::new();
    if have_db {
        if r.viral_genes == Some(0) {
            warnings.push("no viral genes detected".to_string());
        }
        if r.regions.iter().filter(|x| x.viral).count() > 1 {
            warnings.push(">1 viral region detected".to_string());
        }
        r.contamination = Some(match r.host_length() {
            Some(h) => py_round(100.0 * h as f64 / r.length as f64, 2),
            None => 0.0,
        });
    }
    let c = &r.completeness;
    let kmer = np_round(c.kmer_freq, 2);
    if kmer >= 1.5 {
        warnings.push("high kmer_freq may indicate large duplication".to_string());
    }
    let aai_ok = matches!(c.aai_confidence, Some("high") | Some("medium"));
    if aai_ok && (r.length as f64) > 1.5 * c.aai_expected_length.unwrap_or(f64::INFINITY) {
        warnings.push("contig >1.5x longer than expected genome length".to_string());
    }
    let (mut value, mut method) = if aai_ok {
        (c.aai_completeness.map(|v| py_round(v, 2)), format!("AAI-based ({}-confidence)", c.aai_confidence.unwrap_or("NA")))
    } else if let Some((lo, _)) = c.hmm_completeness {
        (Some(py_round(lo, 2)), "HMM-based (lower-bound)".to_string())
    } else {
        (None, "NA".to_string())
    };
    let mut complete = false;
    if let Some(cg) = &r.complete {
        if cg.confidence == "medium" || cg.confidence == "high" {
            complete = true;
            value = Some(100.0);
            method = format!("{} ({}-confidence)", cg.prediction, cg.confidence);
        } else {
            warnings.push(format!("low-confidence {}", cg.prediction));
        }
    }
    let (q, m) = if complete {
        ("Complete", "High-quality")
    } else {
        match value {
            None => ("Not-determined", "Genome-fragment"),
            Some(v) if v >= 90.0 => ("High-quality", "High-quality"),
            Some(v) if v >= 50.0 => ("Medium-quality", "Genome-fragment"),
            Some(_) => ("Low-quality", "Genome-fragment"),
        }
    };
    r.quality = q;
    r.miuvig = m;
    r.completeness_value = value;
    r.completeness_method = method;
    r.warnings = warnings;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gc_counts_n_in_denominator() {
        assert!((gene_gc(b"GGNN", 1, 4) - 50.0).abs() < 1e-12);
        assert!((gene_gc(b"ATGCATGC", 3, 6) - 50.0).abs() < 1e-12);
    }

    fn contig(seq: &[u8]) -> ContigInput {
        ContigInput { id: "c1".into(), seq: seq.to_vec(), genes: Vec::new() }
    }

    fn random_seq(n: usize, mut x: u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x % 4) as usize]
            })
            .collect()
    }

    #[test]
    fn n_run_terminal_repeat_is_not_complete() {
        // CheckV's five test contigs with all-N "DTRs" are low-confidence
        let mut s = vec![b'N'; 60];
        s.extend(random_seq(3000, 99));
        s.extend(vec![b'N'; 60]);
        let r = run(&[contig(&s)], None, &Params::default()).unwrap();
        let cg = r[0].complete.as_ref().unwrap();
        assert_eq!(cg.confidence, "low");
        assert_eq!(cg.reason, "Too many ambiguous bases in TR");
        assert_eq!(r[0].quality, "Not-determined");
        assert!(r[0].warnings.iter().any(|w| w == "low-confidence DTR"));
    }

    #[test]
    fn no_database_means_not_determined() {
        let s = random_seq(5000, 7);
        let r = run(&[contig(&s)], None, &Params::default()).unwrap();
        assert_eq!(r[0].quality, "Not-determined");
        assert!(r[0].contamination.is_none());
        assert!(r[0].complete.is_none());
    }
}
