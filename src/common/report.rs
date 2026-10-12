//! Shared result types. Every validator emits a [`QcRecord`]; the CLI collects
//! them into a [`QcReport`] and renders TSV / JSON. Keeping one record schema
//! lets viruses, bacteria, archaea and eukaryotes share an output table.

use serde::{Deserialize, Serialize};

use crate::common::seqio::Domain;

/// A coarse quality tier, harmonised across the validators (CheckV's
/// Complete/High/Medium/Low; CheckM2's pass/fail; BUSCO's C/F/M summaries all
/// fold onto this).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Quality {
    Complete,
    High,
    Medium,
    Low,
    NotDetermined,
}

impl Quality {
    pub fn label(&self) -> &'static str {
        match self {
            Quality::Complete => "Complete",
            Quality::High => "High-quality",
            Quality::Medium => "Medium-quality",
            Quality::Low => "Low-quality",
            Quality::NotDetermined => "Not-determined",
        }
    }

    /// Standard tiering from a completeness percentage (CheckV's tiers,
    /// which make no contamination claim).
    pub fn from_completeness(comp: Option<f64>) -> Quality {
        match comp {
            None => Quality::NotDetermined,
            Some(c) if c >= 90.0 => Quality::High,
            Some(c) if c >= 50.0 => Quality::Medium,
            Some(_) => Quality::Low,
        }
    }

    /// MIMAG tier from completeness *and* contamination.
    pub fn from_mimag(comp: Option<f64>, cont: Option<f64>) -> Quality {
        match (comp, cont) {
            (Some(c), Some(x)) => match mimag_tier(c, x) {
                "High-quality" => Quality::High,
                "Medium-quality" => Quality::Medium,
                _ => Quality::Low,
            },
            _ => Quality::NotDetermined,
        }
    }
}

/// MIMAG draft-genome tier (Bowers et al. 2017) from completeness and
/// contamination: high > 90 % / < 5 %, medium ≥ 50 % / < 10 %, otherwise
/// low. The rRNA/tRNA requirements of the high tier are not assessed.
pub fn mimag_tier(comp: f64, cont: f64) -> &'static str {
    if comp > 90.0 && cont < 5.0 {
        "High-quality"
    } else if comp >= 50.0 && cont < 10.0 {
        "Medium-quality"
    } else {
        "Low-quality"
    }
}

/// Make a value safe for one TSV cell / `key=value;` list: tabs and line
/// breaks become spaces; `;` and `=` are percent-encoded (`%3B`, `%3D`).
/// A `%` is encoded (`%25`) only when two hex digits follow it, so text like
/// BUSCO's `C:98.4%[S:97.1%,…]` stays readable while decoding every `%XX`
/// still restores the original exactly.
pub fn escape_cell(v: &str) -> String {
    let chars: Vec<char> = v.chars().collect();
    let mut out = String::with_capacity(v.len());
    for (i, &ch) in chars.iter().enumerate() {
        match ch {
            '\t' | '\n' | '\r' => out.push(' '),
            '%' if chars.get(i + 1).is_some_and(|c| c.is_ascii_hexdigit())
                && chars.get(i + 2).is_some_and(|c| c.is_ascii_hexdigit()) =>
            {
                out.push_str("%25")
            }
            ';' => out.push_str("%3B"),
            '=' => out.push_str("%3D"),
            c => out.push(c),
        }
    }
    out
}

fn clean(v: &str) -> String {
    v.replace(['\t', '\n', '\r'], " ")
}

/// One validated contig.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct QcRecord {
    pub contig_id: String,
    pub length: usize,
    pub domain: Domain,
    /// Which validator produced this verdict (e.g. "checkv", "checkm2").
    pub validator: String,
    pub gene_count: usize,
    pub completeness: Option<f64>,
    pub contamination: Option<f64>,
    pub quality: Quality,
    /// Free-form method note (e.g. "AAI-based", "HMM", "ML:specific").
    pub method: String,
    /// Warnings / flags raised during validation.
    pub warnings: Vec<String>,
    /// Validator-specific extra fields, rendered as `key=value;…`.
    pub extra: Vec<(String, String)>,
}

impl QcRecord {
    /// Construct a minimal record.
    pub fn new(contig_id: impl Into<String>, length: usize, domain: Domain, validator: impl Into<String>) -> Self {
        QcRecord {
            contig_id: contig_id.into(),
            length,
            domain,
            validator: validator.into(),
            gene_count: 0,
            completeness: None,
            contamination: None,
            quality: Quality::NotDetermined,
            method: String::new(),
            warnings: Vec::new(),
            extra: Vec::new(),
        }
    }

    /// Add a key/value to the extra column.
    pub fn with_extra(mut self, k: impl Into<String>, v: impl Into<String>) -> Self {
        self.extra.push((k.into(), v.into()));
        self
    }
}

/// A collection of per-contig records.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct QcReport {
    pub records: Vec<QcRecord>,
}

impl QcReport {
    pub fn new() -> Self {
        QcReport::default()
    }
    pub fn push(&mut self, r: QcRecord) {
        self.records.push(r);
    }
    pub fn extend(&mut self, rs: impl IntoIterator<Item = QcRecord>) {
        self.records.extend(rs);
    }

    /// Render a tab-separated table with a header row.
    pub fn to_tsv(&self) -> String {
        let mut out = String::new();
        out.push_str(
            "contig_id\tlength\tdomain\tvalidator\tgene_count\tcompleteness\tcontamination\tquality\tmethod\twarnings\textra\n",
        );
        for r in &self.records {
            let comp = r.completeness.map(|c| format!("{c:.2}")).unwrap_or_else(|| "NA".into());
            let cont = r.contamination.map(|c| format!("{c:.2}")).unwrap_or_else(|| "NA".into());
            let warns = clean(&r.warnings.join("; "));
            let extra =
                r.extra.iter().map(|(k, v)| format!("{}={}", escape_cell(k), escape_cell(v))).collect::<Vec<_>>().join(";");
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                clean(&r.contig_id),
                r.length,
                r.domain.label(),
                r.validator,
                r.gene_count,
                comp,
                cont,
                r.quality.label(),
                clean(&r.method),
                warns,
                extra
            ));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tiers_and_escaping() {
        assert_eq!(mimag_tier(95.0, 40.0), "Low-quality");
        assert_eq!(mimag_tier(95.0, 4.9), "High-quality");
        assert_eq!(mimag_tier(60.0, 9.0), "Medium-quality");
        assert_eq!(escape_cell("a;b=c%\td"), "a%3Bb%3Dc% d");
        assert_eq!(escape_cell("C:98.4%[S:97.1%,D:1.3%]"), "C:98.4%[S:97.1%,D:1.3%]");
        assert_eq!(escape_cell("literal %3B and 100%25"), "literal %253B and 100%2525");
        let mut r = QcRecord::new("c\t1", 10, Domain::Virus, "x");
        r = r.with_extra("k", "v;w");
        let t = QcReport { records: vec![r] }.to_tsv();
        assert!(t.lines().nth(1).unwrap().starts_with("c 1\t10\t"));
        assert!(t.contains("k=v%3Bw"));
    }
}
