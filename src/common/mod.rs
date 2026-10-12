//! Shared foundation: sequence IO, the biological [`Domain`], composition
//! statistics, gene calling via DYNAMITE, the unified [`QcReport`] schema,
//! Python/NumPy-compatible number formatting, and the banner.

pub mod banner;
pub mod composition;
pub mod genes;
pub mod pyfmt;
pub mod report;
pub mod seqio;
#[cfg(feature = "rustyomestats")]
pub mod xstats;

pub use banner::{banner, VERSION};
pub use genes::{call_genes, Gene, GeneSet};
pub use report::{QcRecord, QcReport, Quality};
pub use seqio::{read_fasta, Contig, Domain};

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
