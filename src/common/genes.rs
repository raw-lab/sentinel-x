//! ORF / gene calling, delegated to the **DYNAMITE** engine (`dynamite-orfs`
//! on crates.io). Sentinel calls genes once per run; the protein-level
//! validators (CheckV, CheckM2, BUSCO, EukCC) then annotate the same
//! predicted proteins. This module is the bridge from [`Contig`] to DYNAMITE
//! and back, normalising DYNAMITE's output to Prodigal conventions:
//!
//! * genes are numbered `<contig>_<n>` in coordinate order;
//! * `partial` flags are Prodigal's left/right *contig-edge* flags;
//! * the protein never carries a stop; [`Gene::stop`] records whether the
//!   gene ends in one, and [`Gene::protein_with_stop`] gives Prodigal's
//!   `*`-terminated translation;
//! * proteins are re-translated from the contig as Prodigal does (`X` for
//!   codons with ambiguous bases; a 5′-partial gene's first codon is
//!   translated as itself, not `M`).

use anyhow::Result;

use crate::common::seqio::{Contig, Domain};

/// A predicted protein-coding gene.
#[derive(Clone, Debug)]
pub struct Gene {
    /// Contig the gene was called on.
    pub contig: String,
    /// Index of the contig in the input order.
    pub contig_index: usize,
    /// Per-contig gene index (1-based, in coordinate order).
    pub index: usize,
    /// 1-based inclusive start on the contig (forward coordinates).
    pub begin: usize,
    /// 1-based inclusive end (includes the stop codon of complete genes).
    pub end: usize,
    /// `+1` or `-1`.
    pub strand: i8,
    /// Translated amino-acid sequence, without a stop.
    pub aa: Vec<u8>,
    /// Gene runs off the left contig edge.
    pub partial_left: bool,
    /// Gene runs off the right contig edge.
    pub partial_right: bool,
    /// The gene ends in a stop codon (3′ complete).
    pub stop: bool,
}

impl Gene {
    /// Stable gene id of the form `<contig>_<index>` (the convention CheckV /
    /// CheckM2 use to map proteins back to their parent contig).
    pub fn id(&self) -> String {
        format!("{}_{}", self.contig, self.index)
    }
    /// Prodigal's translation: the protein followed by `*` when the gene
    /// ends in a stop codon.
    pub fn protein_with_stop(&self) -> Vec<u8> {
        let mut p = self.aa.clone();
        if self.stop {
            p.push(b'*');
        }
        p
    }
}

/// All genes called in a run, plus the resolved translation table.
#[derive(Clone, Debug, Default)]
pub struct GeneSet {
    pub genes: Vec<Gene>,
    pub trans_table: u32,
}

impl GeneSet {
    /// Gene indices per contig (by contig input index), in coordinate order.
    pub fn per_contig(&self, n_contigs: usize) -> Vec<Vec<usize>> {
        let mut v: Vec<Vec<usize>> = vec![Vec::new(); n_contigs];
        for (i, g) in self.genes.iter().enumerate() {
            if g.contig_index < n_contigs {
                v[g.contig_index].push(i);
            }
        }
        v
    }

    /// Return `(id, residues)` pairs (no stops).
    pub fn protein_records(&self) -> Vec<(String, Vec<u8>)> {
        self.genes.iter().map(|g| (g.id(), g.aa.clone())).collect()
    }

    /// Return `(id, residues)` pairs with Prodigal's terminal `*`.
    pub fn protein_records_with_stop(&self) -> Vec<(String, Vec<u8>)> {
        self.genes.iter().map(|g| (g.id(), g.protein_with_stop())).collect()
    }
}

/// Gene-calling options beyond the domain preset.
#[derive(Clone, Copy, Debug, Default)]
pub struct CallOptions {
    /// Override the preset (`None` = the domain's preset).
    pub mode: Option<dynamite::model::Mode>,
    /// Disallow genes running off contig ends (Prodigal `-c`).
    pub closed: bool,
    /// Force a translation table (Prodigal `-g`).
    pub code: Option<u32>,
}

/// Call genes on a set of contigs using the DYNAMITE preset for `domain`.
pub fn call_genes(contigs: &[Contig], domain: Domain, threads: usize) -> Result<GeneSet> {
    call_genes_with(contigs, domain, threads, CallOptions::default())
}

/// [`call_genes`] with explicit options.
pub fn call_genes_with(contigs: &[Contig], domain: Domain, threads: usize, opt: CallOptions) -> Result<GeneSet> {
    let records: Vec<dynamite::fasta::Record> = contigs
        .iter()
        .map(|c| dynamite::fasta::Record { id: c.id.clone(), desc: c.header.clone(), seq: c.seq.clone() })
        .collect();

    let cfg = dynamite::RunConfig {
        mode: opt.mode.unwrap_or_else(|| domain.dynamite_mode()),
        engine: None,
        code: match opt.code {
            Some(t) => dynamite::model::CodeChoice::Fixed(t),
            None => dynamite::model::CodeChoice::Default,
        },
        min_len: 0,
        closed: opt.closed,
        threads,
    };

    let out = dynamite::run_records(&records, &cfg)?;
    let index: std::collections::HashMap<&str, usize> = contigs.iter().enumerate().map(|(i, c)| (c.id.as_str(), i)).collect();

    let mut genes: Vec<Gene> = Vec::with_capacity(out.genes.len());
    for g in &out.genes {
        let Some(&ci) = index.get(g.contig.as_str()) else { continue };
        genes.push(normalise(g, ci, &contigs[ci].seq));
    }
    // coordinate order within each contig, contigs in input order
    genes.sort_by_key(|g| (g.contig_index, g.begin, g.end, g.strand));
    let mut cur = usize::MAX;
    let mut n = 0usize;
    for g in genes.iter_mut() {
        if g.contig_index != cur {
            cur = g.contig_index;
            n = 0;
        }
        n += 1;
        g.index = n;
    }
    Ok(GeneSet { genes, trans_table: out.trans_table })
}

/// Convert one DYNAMITE gene to Prodigal conventions and re-translate it
/// from the contig the way Prodigal does. DYNAMITE's engines set
/// `partial_left/right` from the contig edges but its translation treats
/// them as 5′/3′ flags (a minus-strand gene running off the right edge keeps
/// its stop and may get `M` for a non-ATG edge codon), and it translates
/// codons containing `N` as if the `N` were `A` instead of `X`.
fn normalise(g: &dynamite::gene::Gene, contig_index: usize, contig: &[u8]) -> Gene {
    let (five_partial, three_partial) =
        if g.strand >= 0 { (g.partial_left, g.partial_right) } else { (g.partial_right, g.partial_left) };
    let fallback = || {
        let mut aa: Vec<u8> = g.aa.as_bytes().to_vec();
        while aa.last() == Some(&b'*') {
            aa.pop();
        }
        aa
    };
    let (aa, stop) = match dynamite::gencode::GeneticCode::new(g.trans_table) {
        Some(code) if g.begin >= 1 && g.end <= contig.len() && g.begin <= g.end => {
            translate(contig, g.begin, g.end, g.strand, &code, !five_partial, !three_partial)
        }
        _ => (fallback(), !three_partial),
    };
    Gene {
        contig: g.contig.clone(),
        contig_index,
        index: 0,
        begin: g.begin,
        end: g.end,
        strand: if g.strand >= 0 { 1 } else { -1 },
        aa,
        partial_left: g.partial_left,
        partial_right: g.partial_right,
        stop,
    }
}

/// Prodigal's translation of `contig[begin-1..end]` on `strand`: `X` for
/// codons with an ambiguous base, `M` for the start codon of a 5′-complete
/// gene, and the terminal stop removed (reported in the second value).
pub fn translate(
    contig: &[u8],
    begin: usize,
    end: usize,
    strand: i8,
    code: &dynamite::gencode::GeneticCode,
    five_complete: bool,
    three_complete: bool,
) -> (Vec<u8>, bool) {
    let span = &contig[begin - 1..end];
    let coding: Vec<u8> = if strand >= 0 {
        span.to_vec()
    } else {
        span.iter()
            .rev()
            .map(|&b| match b.to_ascii_uppercase() {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' | b'U' => b'A',
                _ => b'N',
            })
            .collect()
    };
    let n = coding.len() / 3;
    let mut aa = Vec::with_capacity(n);
    for k in 0..n {
        let c = &coding[3 * k..3 * k + 3];
        let (b1, a1) = dynamite::seq::encode_base(c[0]);
        let (b2, a2) = dynamite::seq::encode_base(c[1]);
        let (b3, a3) = dynamite::seq::encode_base(c[2]);
        let r = if a1 || a2 || a3 {
            b'X'
        } else if k == 0 && five_complete && code.is_start(b1, b2, b3) {
            b'M'
        } else {
            code.translate(b1, b2, b3)
        };
        aa.push(r);
    }
    let stop = three_complete && aa.last() == Some(&b'*');
    if stop {
        aa.pop();
    }
    (aa, stop)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn genes_are_numbered_in_coordinate_order_with_stops() {
        // a simple ORF-rich synthetic contig
        let unit = b"ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA";
        let mut seq = Vec::new();
        for _ in 0..40 {
            seq.extend_from_slice(unit);
            seq.extend_from_slice(b"TTTTTTAGGAGGTTTTTT");
        }
        let c = Contig { id: "ctg".into(), header: "ctg".into(), seq };
        let gs = call_genes(std::slice::from_ref(&c), Domain::Virus, 1).unwrap();
        let mut last = 0;
        for (k, g) in gs.genes.iter().enumerate() {
            assert_eq!(g.index, k + 1);
            assert!(g.begin >= last);
            last = g.begin;
            assert!(!g.aa.contains(&b'*'), "internal or trailing stop in {}", g.id());
            assert_eq!(g.protein_with_stop().last() == Some(&b'*'), g.stop);
        }
    }
}
