//! # sentinel-viralverify
//!
//! A naive-Bayes classifier that labels each contig **Virus**, **Plasmid**,
//! **Chromosome**, or **Uncertain** from the Pfam domains its genes carry — the
//! viralVerify method. For every Pfam HMM, the shipped frequency table holds
//! how often that domain occurs in viral, plasmid, chromosomal, and
//! (plasmid+chromosomal) reference sequences; the classifier sums the
//! log-frequencies of all domains hitting a contig and compares the resulting
//! log-likelihoods.
//!
//! The frequency table (`classifier_table_uniq.txt` from viralVerify) is
//! embedded, so the only external requirement for a full run is a Pfam HMM
//! database to generate the per-gene hits.

use std::collections::HashMap;

use sentinel_core::{Domain, GeneSet, QcRecord, Quality};
use sentinel_hmmer::Hit;

/// The embedded Pfam frequency table.
const TABLE_TSV: &str = include_str!("../../../assets/viralverify_pfam_table.tsv");

/// Per-HMM log-frequencies for the four reference classes.
#[derive(Clone, Copy, Debug)]
struct Row {
    plasmid: f64,
    chromosome: f64,
    virus: f64,
    /// plasmid+chromosome combined ("not viral").
    combined: f64,
}

/// The naive-Bayes frequency table keyed by Pfam HMM name.
#[derive(Clone, Debug)]
pub struct ClassifierTable {
    rows: HashMap<String, Row>,
}

impl Default for ClassifierTable {
    fn default() -> Self {
        Self::embedded()
    }
}

impl ClassifierTable {
    /// Load the table embedded in the binary.
    pub fn embedded() -> Self {
        Self::parse(TABLE_TSV)
    }

    /// Parse a viralVerify classifier table:
    /// `name  p_count  c_count  v_count  p_freq  c_freq  v_freq  pc_freq`.
    pub fn parse(text: &str) -> Self {
        let mut rows = HashMap::new();
        for line in text.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 8 {
                continue;
            }
            let name = strip_version(f[0]).to_string();
            let p = f[4].parse::<f64>().unwrap_or(0.0);
            let c = f[5].parse::<f64>().unwrap_or(0.0);
            let v = f[6].parse::<f64>().unwrap_or(0.0);
            let pc = f[7].parse::<f64>().unwrap_or(0.0);
            const FLOOR: f64 = 1e-9;
            rows.insert(
                name,
                Row {
                    plasmid: p.max(FLOOR).ln(),
                    chromosome: c.max(FLOOR).ln(),
                    virus: v.max(FLOOR).ln(),
                    combined: pc.max(FLOOR).ln(),
                },
            );
        }
        ClassifierTable { rows }
    }

    /// Number of HMMs in the table.
    pub fn len(&self) -> usize {
        self.rows.len()
    }
    /// True if the table is empty.
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }
}

/// Strip a trailing `.NN` Pfam version so `PF00069.27` matches `PF00069`.
fn strip_version(s: &str) -> &str {
    match s.split_once('.') {
        Some((base, ver)) if ver.chars().all(|c| c.is_ascii_digit()) => base,
        _ => s,
    }
}

/// The four viralVerify verdicts plus the uncertainty bands.
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

/// Outcome for one contig.
#[derive(Clone, Debug)]
pub struct Verdict {
    pub contig_id: String,
    pub class: Class,
    pub viral_score: f64,
    pub plasmid_score: f64,
    pub n_domains: usize,
}

/// Classify a single set of Pfam HMM names (one contig's domains).
pub fn classify_domains(
    domains: &[String],
    table: &ClassifierTable,
    threshold: f64,
) -> (Class, f64, f64, usize) {
    let mut vir = 0.0;
    let mut comb = 0.0;
    let mut plas = 0.0;
    let mut chrom = 0.0;
    let mut n = 0usize;
    for d in domains {
        if let Some(r) = table.rows.get(strip_version(d)) {
            vir += r.virus;
            comb += r.combined;
            plas += r.plasmid;
            chrom += r.chromosome;
            n += 1;
        }
    }
    let viral_score = vir - comb;
    let plasmid_score = plas - chrom;

    let class = if n == 0 {
        Class::UncertainTooShort
    } else if viral_score > threshold {
        Class::Virus
    } else if viral_score > -threshold {
        Class::UncertainViralOrBacterial
    } else if plasmid_score > threshold {
        Class::Plasmid
    } else if -plasmid_score > threshold {
        Class::Chromosome
    } else {
        Class::UncertainPlasmidOrChromosome
    };
    (class, viral_score, plasmid_score, n)
}

/// Run viralVerify over all contigs given the per-gene Pfam hits.
pub fn classify(
    geneset: &GeneSet,
    hits: &[Hit],
    table: &ClassifierTable,
    threshold: f64,
) -> Vec<Verdict> {
    let mut by_contig: HashMap<String, Vec<String>> = HashMap::new();
    let gene_to_contig: HashMap<String, String> = geneset
        .genes
        .iter()
        .map(|g| (g.id(), g.contig.clone()))
        .collect();

    for h in hits {
        if let Some(contig) = gene_to_contig.get(&h.target) {
            by_contig
                .entry(contig.clone())
                .or_default()
                .push(h.profile.clone());
        }
    }

    let mut seen: Vec<&str> = Vec::new();
    for g in &geneset.genes {
        if !seen.contains(&g.contig.as_str()) {
            seen.push(g.contig.as_str());
        }
    }
    let mut out = Vec::new();
    for contig in seen {
        let domains = by_contig.get(contig).cloned().unwrap_or_default();
        let (class, vs, ps, n) = classify_domains(&domains, table, threshold);
        out.push(Verdict {
            contig_id: contig.to_string(),
            class,
            viral_score: vs,
            plasmid_score: ps,
            n_domains: n,
        });
    }
    out
}

/// Convert verdicts to the unified report schema.
pub fn to_records(verdicts: &[Verdict], lengths: &HashMap<String, usize>) -> Vec<QcRecord> {
    verdicts
        .iter()
        .map(|v| {
            let len = lengths.get(&v.contig_id).copied().unwrap_or(0);
            let mut r = QcRecord::new(v.contig_id.clone(), len, Domain::Virus, "viralverify");
            r.method = "naive-bayes".into();
            r.quality = match v.class {
                Class::Virus => Quality::High,
                _ => Quality::NotDetermined,
            };
            r = r
                .with_extra("prediction", v.class.label())
                .with_extra("viral_score", format!("{:.2}", v.viral_score))
                .with_extra("plasmid_score", format!("{:.2}", v.plasmid_score))
                .with_extra("pfam_domains", v.n_domains.to_string());
            r
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_loads() {
        let t = ClassifierTable::embedded();
        assert!(t.len() > 1000, "expected a populated Pfam table, got {}", t.len());
    }

    #[test]
    fn no_domains_is_too_short() {
        let t = ClassifierTable::embedded();
        let (c, _, _, n) = classify_domains(&[], &t, 7.0);
        assert_eq!(c, Class::UncertainTooShort);
        assert_eq!(n, 0);
    }

    #[test]
    fn strongly_viral_domains_classify_virus() {
        let t = ClassifierTable::embedded();
        let mut scored: Vec<(String, f64)> = t
            .rows
            .iter()
            .map(|(k, r)| (k.clone(), r.virus - r.combined))
            .collect();
        scored.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        let top: Vec<String> = scored.iter().take(5).map(|(k, _)| k.clone()).collect();
        let (c, vs, _, _) = classify_domains(&top, &t, 7.0);
        assert!(vs > 0.0);
        assert_eq!(c, Class::Virus, "viral_score={vs}");
    }
}
