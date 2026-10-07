//! Shared result types. Every validator emits a [`QcRecord`]; the CLI collects
//! them into a [`QcReport`] and renders TSV / JSON. Keeping one record schema
//! lets viruses, bacteria, archaea and eukaryotes share an output table.

use serde::{Deserialize, Serialize};

use crate::seqio::Domain;

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

    /// Standard tiering from a completeness percentage.
    pub fn from_completeness(comp: Option<f64>) -> Quality {
        match comp {
            None => Quality::NotDetermined,
            Some(c) if c >= 90.0 => Quality::High,
            Some(c) if c >= 50.0 => Quality::Medium,
            Some(_) => Quality::Low,
        }
    }
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
            let comp = r
                .completeness
                .map(|c| format!("{c:.2}"))
                .unwrap_or_else(|| "NA".into());
            let cont = r
                .contamination
                .map(|c| format!("{c:.2}"))
                .unwrap_or_else(|| "NA".into());
            let warns = r.warnings.join("; ");
            let extra = r
                .extra
                .iter()
                .map(|(k, v)| format!("{k}={v}"))
                .collect::<Vec<_>>()
                .join(";");
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                r.contig_id,
                r.length,
                r.domain.label(),
                r.validator,
                r.gene_count,
                comp,
                cont,
                r.quality.label(),
                r.method,
                warns,
                extra
            ));
        }
        out
    }
}
