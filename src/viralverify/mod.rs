//! # viralVerify, ported
//!
//! Naive-Bayes classification of each contig as **Virus**, **Plasmid**,
//! **Chromosome** or **Uncertain** from the set of Pfam domains on its genes
//! (Antipov et al. 2020), following `bin/viralverify` (v1.1) step by step:
//!
//! 1. **Circularity** — an identical 51–200-mer at both ends of a contig of
//!    at least 500 bp marks it circular (`+`), and the gene caller sees the
//!    contig extended by `seq[k:5000]` so genes can span the junction; genes
//!    starting beyond the original end are dropped.
//! 2. **Genes** — metagenomic Prodigal with closed ends (`prodigal -p meta
//!    -c`), here DYNAMITE's meta preset with `closed = true`.
//! 3. **Domains** — `hmmsearch --cut_nc` against Pfam (sequence ≥ NC1 and
//!    domain ≥ NC2). Per gene the domain rows are sorted by (gene, domain
//!    score, alignment start *as text*) descending; the first is kept, and
//!    a later one is kept only if it does not overlap that first domain.
//! 4. **Score** — each contig's *set* of HMM names: viral score
//!    `Σ ln v − Σ ln pc`, plasmid score `Σ ln p − Σ ln c` over names in the
//!    classifier table; threshold 7 by default.
//! 5. **Rescue** — "too short" contigs (and contigs without hits) become
//!    "Uncertain - viral or bacterial" when longer than 3 kb or circular.

use std::collections::{HashMap, HashSet};
use std::path::Path;

use anyhow::Result;

use crate::common::genes::{call_genes_with, CallOptions};
use crate::common::pyfmt::reparse;
use crate::common::{Contig, Domain, QcRecord, Quality};
use crate::hmmer::{search_hmm_file, Cutoff, SearchConfig, Targets};

/// The embedded classifier table (`db/classifier_table_uniq.txt`).
const TABLE_TSV: &str = include_str!("../../assets/viralverify_pfam_table.tsv");

/// Frequencies (p, c, v, p+c) per HMM name.
#[derive(Clone, Debug)]
pub struct ClassifierTable {
    rows: HashMap<String, [f64; 4]>,
}

impl Default for ClassifierTable {
    fn default() -> Self {
        Self::embedded()
    }
}

impl ClassifierTable {
    pub fn embedded() -> Self {
        Self::parse(TABLE_TSV)
    }

    /// `name p_count c_count v_count p_freq c_freq v_freq pc_freq`, split on
    /// any whitespace as viralVerify does.
    pub fn parse(text: &str) -> Self {
        let mut rows = HashMap::new();
        for line in text.lines() {
            let f: Vec<&str> = line.split_whitespace().collect();
            if f.len() < 8 {
                continue;
            }
            let v = |i: usize| f[i].parse::<f64>().unwrap_or(f64::NAN);
            rows.insert(f[0].to_string(), [v(4), v(5), v(6), v(7)]);
        }
        ClassifierTable { rows }
    }

    #[cfg(test)]
    pub fn len(&self) -> usize {
        self.rows.len()
    }
}

/// viralVerify's verdicts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Class {
    Virus,
    Plasmid,
    Chromosome,
    UncertainViralOrBacterial,
    UncertainPlasmidOrChromosome,
    UncertainTooShort,
}

impl Class {
    pub fn label(&self) -> &'static str {
        match self {
            Class::Virus => "Virus",
            Class::Plasmid => "Plasmid",
            Class::Chromosome => "Chromosome",
            Class::UncertainViralOrBacterial => "Uncertain - viral or bacterial",
            Class::UncertainPlasmidOrChromosome => "Uncertain - plasmid or chromosomal",
            Class::UncertainTooShort => "Uncertain - too short",
        }
    }
}

/// One contig's result (a `_result_table.csv` row).
#[derive(Clone, Debug)]
pub struct Verdict {
    pub contig_id: String,
    pub length: usize,
    pub circular: bool,
    pub class: Class,
    /// The printed score (`"{:.2f}"`), `None` for contigs without hits.
    pub score: Option<String>,
    /// Distinct HMM names on the contig (sorted).
    pub hmms: Vec<String>,
}

/// `naive_bayes` for one contig's HMM-name set.
pub fn naive_bayes(names: &[String], table: &ClassifierTable, thr: f64) -> (Class, String) {
    let (mut p, mut c, mut v, mut pc) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    for n in names {
        if let Some(r) = table.rows.get(n) {
            p += r[0].ln();
            c += r[1].ln();
            v += r[2].ln();
            pc += r[3].ln();
        }
    }
    let vs = v - pc;
    let ps = p - c;
    if vs > thr {
        (Class::Virus, format!("{vs:.2}"))
    } else if vs > -thr {
        if names.len() > 2 {
            (Class::UncertainViralOrBacterial, format!("{vs:.2}"))
        } else {
            (Class::UncertainTooShort, format!("{vs:.2}"))
        }
    } else if ps > thr {
        (Class::Plasmid, format!("{ps:.2}"))
    } else if c - p > thr {
        (Class::Chromosome, format!("{ps:.2}"))
    } else {
        (Class::UncertainPlasmidOrChromosome, format!("{ps:.2}"))
    }
}

/// `check_circular`: the k (200 down to 51) of an identical k-mer at both
/// ends of a contig of ≥ 500 bp.
pub fn check_circular(seq: &[u8]) -> Option<usize> {
    if seq.len() < 500 {
        return None;
    }
    (51..=200).rev().find(|&k| k < seq.len() && seq[..k] == seq[seq.len() - k..])
}

/// A domain row (`--domtblout` fields viralVerify reads).
#[derive(Clone, Debug)]
struct DomRow {
    gene: String,
    hmm: String,
    score: f64,
    ali_from: usize,
    ali_to: usize,
}

/// Run viralVerify on the contigs with a Pfam HMM library.
pub fn run(contigs: &[Contig], pfam: &Path, threads: usize, thr: f64, table: &ClassifierTable) -> Result<Vec<Verdict>> {
    // 1. circularity and extended sequences
    let circ: Vec<Option<usize>> = contigs.iter().map(|c| check_circular(&c.seq)).collect();
    let extended: Vec<Contig> = contigs
        .iter()
        .zip(&circ)
        .map(|(c, k)| {
            let mut seq = c.seq.clone();
            if let Some(k) = *k {
                let end = c.seq.len().min(5000);
                if k < end {
                    seq.extend_from_slice(&c.seq[k..end]);
                }
            }
            Contig { id: c.id.clone(), header: c.header.clone(), seq }
        })
        .collect();

    // 2. genes: prodigal -p meta -c; drop genes starting past the original end
    let opt = CallOptions { mode: Some(dynamite::model::Mode::Meta), closed: true, code: None };
    let gs = call_genes_with(&extended, Domain::Mixed, threads, opt)?;
    let mut ids: Vec<String> = Vec::new();
    let mut prots: Vec<Vec<u8>> = Vec::new();
    for g in &gs.genes {
        if g.begin < contigs[g.contig_index].seq.len() {
            ids.push(g.id());
            prots.push(g.protein_with_stop());
        }
    }

    // validation aid: the proteins as searched
    if let Ok(path) = std::env::var("SENTINEL_VIRALVERIFY_PROTEINS") {
        let mut t = String::new();
        for (id, p) in ids.iter().zip(&prots) {
            t.push('>');
            t.push_str(id);
            t.push('\n');
            t.push_str(&String::from_utf8_lossy(p));
            t.push('\n');
        }
        std::fs::write(path, t)?;
    }

    // 3. Pfam --cut_nc, per-gene domain selection
    let mut rows: Vec<DomRow> = Vec::new();
    if !ids.is_empty() {
        let targets = Targets::new(ids.iter().map(|s| s.as_str()).zip(prots.iter().map(|p| p.as_slice())));
        let cfg = SearchConfig { cutoff: Cutoff::Noise, ..Default::default() };
        let (hits, _) = search_hmm_file(pfam, &targets, &cfg)?;
        for h in &hits {
            for d in h.reported_domains() {
                rows.push(DomRow {
                    gene: h.target.clone(),
                    hmm: h.profile.clone(),
                    score: reparse(&format!("{:.1}", d.score)),
                    ali_from: d.ali_from,
                    ali_to: d.ali_to,
                });
            }
        }
    }
    // sort by (gene, score, ali_from as text), descending, stable
    rows.sort_by(|a, b| {
        b.gene
            .cmp(&a.gene)
            .then(b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal))
            .then(b.ali_from.to_string().cmp(&a.ali_from.to_string()))
    });
    let mut kept: Vec<(String, Vec<String>)> = Vec::new(); // gene -> hmm names
    let mut first: HashMap<String, (usize, usize)> = HashMap::new();
    for r in &rows {
        match first.get(&r.gene) {
            None => {
                first.insert(r.gene.clone(), (r.ali_from, r.ali_to));
                kept.push((r.gene.clone(), vec![r.hmm.clone()]));
            }
            Some(&(s, e)) => {
                // only the first kept domain is checked (both branches break)
                if (r.ali_to as f64) <= s as f64 || (r.ali_from as f64) >= e as f64 {
                    if let Some(k) = kept.iter_mut().rev().find(|k| k.0 == r.gene) {
                        k.1.push(r.hmm.clone());
                    }
                }
            }
        }
    }
    // per-contig sets
    let mut sets: HashMap<String, HashSet<String>> = HashMap::new();
    for (gene, hmms) in kept {
        let contig = gene.rsplit_once('_').map(|x| x.0.to_string()).unwrap_or(gene);
        sets.entry(contig).or_default().extend(hmms);
    }

    // 4–5. classify every input contig
    let mut out = Vec::with_capacity(contigs.len());
    for (c, k) in contigs.iter().zip(&circ) {
        let len = c.seq.len();
        let circular = k.is_some();
        let rescue = len > 3000 || circular;
        let v = match sets.get(&c.id) {
            Some(set) => {
                let mut names: Vec<String> = set.iter().cloned().collect();
                names.sort();
                let (mut class, score) = naive_bayes(&names, table, thr);
                if class == Class::UncertainTooShort && rescue {
                    class = Class::UncertainViralOrBacterial;
                }
                Verdict { contig_id: c.id.clone(), length: len, circular, class, score: Some(score), hmms: names }
            }
            None => Verdict {
                contig_id: c.id.clone(),
                length: len,
                circular,
                class: if rescue { Class::UncertainViralOrBacterial } else { Class::UncertainTooShort },
                score: None,
                hmms: Vec::new(),
            },
        };
        out.push(v);
    }
    Ok(out)
}

/// viralVerify's `<name>_result_table.csv`.
pub fn result_table_csv(verdicts: &[Verdict]) -> String {
    let mut s = String::from("Contig name,Prediction,Length,Circular,Score,Pfam hits\n");
    for v in verdicts {
        let field = |x: &str| -> String {
            if x.contains(',') || x.contains('"') {
                format!("\"{}\"", x.replace('"', "\"\""))
            } else {
                x.to_string()
            }
        };
        s.push_str(&format!(
            "{},{},{},{},{}",
            field(&v.contig_id),
            v.class.label(),
            v.length,
            if v.circular { "+" } else { "-" },
            v.score.as_deref().unwrap_or("-")
        ));
        if v.score.is_some() {
            s.push(',');
            s.push_str(&field(&v.hmms.join(" ")));
        }
        s.push('\n');
    }
    s
}

/// Unified-report records.
pub fn to_records(verdicts: &[Verdict], domain: Domain) -> Vec<QcRecord> {
    verdicts
        .iter()
        .map(|v| {
            let mut r = QcRecord::new(&v.contig_id, v.length, domain, "viralverify");
            r.method = "pfam-naive-bayes".into();
            r.quality = Quality::NotDetermined;
            r = r
                .with_extra("prediction", v.class.label())
                .with_extra("circular", if v.circular { "+" } else { "-" })
                .with_extra("score", v.score.clone().unwrap_or_else(|| "-".into()))
                .with_extra("pfam_hits", v.hmms.len().to_string());
            r
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_and_bayes() {
        let t = ClassifierTable::embedded();
        assert!(t.len() > 8000);
        // no names in the table -> score 0.00, too short with <= 2 names
        let (c, s) = naive_bayes(&["NotAPfam".into()], &t, 7.0);
        assert_eq!((c, s.as_str()), (Class::UncertainTooShort, "0.00"));
        let three: Vec<String> = vec!["x".into(), "y".into(), "z".into()];
        assert_eq!(naive_bayes(&three, &t, 7.0).0, Class::UncertainViralOrBacterial);
    }

    #[test]
    fn circularity_window() {
        let mut x = 12345u64;
        let mut s: Vec<u8> = (0..1000)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x % 4) as usize]
            })
            .collect();
        assert_eq!(check_circular(&s), None);
        let head = s[..60].to_vec();
        s.extend_from_slice(&head);
        assert_eq!(check_circular(&s), Some(60));
        assert_eq!(check_circular(&s[..400]), None);
    }
}
