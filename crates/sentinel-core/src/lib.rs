//! # sentinel-core
//!
//! The shared foundation for Sentinel: sequence IO, the biological [`Domain`]
//! model, composition statistics, ORF/gene calling via DYNAMITE, the unified
//! [`QcReport`] result schema, and the SENTINEL [`banner`].
//!
//! Every validator crate (`sentinel-checkv`, `sentinel-checkm2`,
//! `sentinel-busco`, `sentinel-viralverify`, `sentinel-novel`) builds on these
//! types so a single run produces one coherent, multi-domain report.

pub mod banner;
pub mod composition;
pub mod genes;
pub mod report;
pub mod seqio;

pub use banner::{banner, ART, VERSION};
pub use genes::{call_genes, Gene, GeneSet};
pub use report::{QcRecord, QcReport, Quality};
pub use seqio::{read_fasta, read_fasta_reader, Contig, Domain};

/// N50 of a set of contig lengths (a standard assembly-quality summary).
pub fn n50(lengths: &[usize]) -> usize {
    if lengths.is_empty() {
        return 0;
    }
    let mut ls: Vec<usize> = lengths.to_vec();
    ls.sort_unstable_by(|a, b| b.cmp(a));
    let total: usize = ls.iter().sum();
    let mut acc = 0usize;
    for l in ls {
        acc += l;
        if acc * 2 >= total {
            return l;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn n50_basic() {
        // lengths 2,3,4,5,6 -> total 20; sorted desc 6,5,..; 6 then 11 (>=10) -> 5
        assert_eq!(n50(&[2, 3, 4, 5, 6]), 5);
        assert_eq!(n50(&[]), 0);
    }
}
