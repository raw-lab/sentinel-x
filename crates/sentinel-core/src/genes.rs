//! ORF / gene calling, delegated to the **DYNAMITE** engine (vendored). Sentinel
//! calls genes exactly once per run; every protein-level validator (CheckV,
//! CheckM2, BUSCO, viralVerify) then annotates the same predicted proteins via
//! the shared pure-Rust HMMER engine. This module is the bridge from Sentinel's
//! [`Contig`](crate::Contig) type to DYNAMITE and back.

use anyhow::Result;

use crate::seqio::{Contig, Domain};

/// A predicted protein-coding gene.
#[derive(Clone, Debug)]
pub struct Gene {
    /// Contig the gene was called on.
    pub contig: String,
    /// Per-contig gene index (1-based), used to build `<contig>_<n>` ids.
    pub index: usize,
    /// 1-based inclusive start on the contig (forward coordinates).
    pub begin: usize,
    /// 1-based inclusive end.
    pub end: usize,
    /// `+1` or `-1`.
    pub strand: i8,
    /// Translated amino-acid sequence (stop codon stripped).
    pub aa: Vec<u8>,
    /// Whether the gene runs off a contig end (partial).
    pub partial: bool,
}

impl Gene {
    /// Stable gene id of the form `<contig>_<index>` (the convention CheckV /
    /// CheckM2 use to map proteins back to their parent contig).
    pub fn id(&self) -> String {
        format!("{}_{}", self.contig, self.index)
    }
    /// Protein length in residues.
    pub fn protein_len(&self) -> usize {
        self.aa.len()
    }
}

/// All genes called in a run, plus the resolved translation table.
#[derive(Clone, Debug, Default)]
pub struct GeneSet {
    pub genes: Vec<Gene>,
    pub trans_table: u32,
}

impl GeneSet {
    /// Group gene indices by contig id, preserving call order.
    pub fn by_contig(&self) -> std::collections::HashMap<&str, Vec<usize>> {
        let mut map: std::collections::HashMap<&str, Vec<usize>> = std::collections::HashMap::new();
        for (i, g) in self.genes.iter().enumerate() {
            map.entry(g.contig.as_str()).or_default().push(i);
        }
        map
    }

    /// Return `(id, residues)` pairs for HMMER search.
    pub fn protein_records(&self) -> Vec<(String, Vec<u8>)> {
        self.genes
            .iter()
            .map(|g| (g.id(), g.aa.clone()))
            .collect()
    }
}

/// Call genes on a set of contigs using the DYNAMITE preset for `domain`.
pub fn call_genes(contigs: &[Contig], domain: Domain, threads: usize) -> Result<GeneSet> {
    let records: Vec<dynamite_core::fasta::Record> = contigs
        .iter()
        .map(|c| dynamite_core::fasta::Record {
            id: c.id.clone(),
            desc: c.header.clone(),
            seq: c.seq.clone(),
        })
        .collect();

    let cfg = dynamite_core::RunConfig {
        mode: domain.dynamite_mode(),
        engine: None,
        code: dynamite_core::model::CodeChoice::Default,
        min_len: 0,
        closed: false,
        threads,
    };

    let out = dynamite_core::run_records(&records, &cfg)?;

    // DYNAMITE returns genes in call order; assign per-contig 1-based indices.
    let mut per_contig: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let mut genes = Vec::with_capacity(out.genes.len());
    for g in &out.genes {
        let n = per_contig.entry(g.contig.clone()).or_insert(0);
        *n += 1;
        genes.push(Gene {
            contig: g.contig.clone(),
            index: *n,
            begin: g.begin,
            end: g.end,
            strand: g.strand,
            aa: g.aa.as_bytes().to_vec(),
            partial: g.partial_left || g.partial_right,
        });
    }

    Ok(GeneSet {
        genes,
        trans_table: out.trans_table,
    })
}
