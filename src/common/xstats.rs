//! Optional cross-check of Sentinel-X's assembly statistics against the
//! crates.io [`rustyomestats`] crate (`--features rustyomestats`).
//!
//! rustyomestats 0.2.x needs rustc >= 1.88 and pulls in polars, plotters and
//! rust-bio, so it is not part of the default build. With the feature on,
//! `sentinel-x run` recomputes sequence count, total length, base composition
//! and N50/L50 for every input file with rustyomestats and reports whether
//! the two implementations agree.
//!
//! Known difference: rustyomestats' `compute_nl` takes the 50 % threshold as
//! `total * 50 / 100` (integer division). For an odd total it can stop one
//! contig early, giving a larger N50 / smaller L50 than the standard
//! definition (smallest set of longest contigs covering at least half of the
//! assembly), e.g. lengths 4, 3, 2 -> rustyomestats N50 4 / L50 1, standard
//! N50 3 / L50 2.

use crate::common::Contig;

/// Base counts and N50/L50 for a set of sequences.
#[derive(Debug, PartialEq, Eq, Clone, Copy, Default)]
pub struct BasicCounts {
    pub seqs: usize,
    pub total: usize,
    pub gc: usize,
    pub at: usize,
    pub n: usize,
    pub n50: usize,
    pub l50: usize,
}

/// Sentinel-X's own numbers. Bases are counted as rustyomestats counts them
/// (G/C; A/T/U; N; anything else ignored); N50/L50 use the standard
/// definition (cumulative length * 2 >= total).
pub fn sentinel_x_counts(contigs: &[Contig]) -> BasicCounts {
    let mut c = BasicCounts { seqs: contigs.len(), ..Default::default() };
    let mut lengths: Vec<usize> = Vec::with_capacity(contigs.len());
    for ctg in contigs {
        c.total += ctg.len();
        lengths.push(ctg.len());
        for &b in &ctg.seq {
            match b.to_ascii_uppercase() {
                b'G' | b'C' => c.gc += 1,
                b'A' | b'T' | b'U' => c.at += 1,
                b'N' => c.n += 1,
                _ => {}
            }
        }
    }
    lengths.sort_unstable_by(|a, b| b.cmp(a));
    let mut acc = 0usize;
    for (i, &l) in lengths.iter().enumerate() {
        acc += l;
        if acc * 2 >= c.total {
            c.n50 = l;
            c.l50 = i + 1;
            break;
        }
    }
    c
}

/// The same numbers from rustyomestats (`stats::compute_basic` +
/// `stats::compute_nl`). Records are converted in batches so the copy made
/// for rust-bio stays bounded.
pub fn rustyomestats_counts(contigs: &[Contig]) -> BasicCounts {
    const BATCH_BASES: usize = 64 << 20;
    let mut c = BasicCounts::default();
    let mut lengths: Vec<usize> = Vec::with_capacity(contigs.len());
    let mut start = 0usize;
    while start < contigs.len() {
        let mut end = start;
        let mut bases = 0usize;
        while end < contigs.len() && (end == start || bases + contigs[end].len() <= BATCH_BASES) {
            bases += contigs[end].len();
            end += 1;
        }
        let recs: Vec<bio::io::fasta::Record> =
            contigs[start..end].iter().map(|x| bio::io::fasta::Record::with_attrs(&x.id, None, &x.seq)).collect();
        let b = rustyomestats::stats::compute_basic(&recs);
        c.seqs += b.num_seq;
        c.total += b.total_bp;
        c.gc += b.gc_count;
        c.at += b.at_count;
        c.n += b.n_count;
        lengths.extend(b.lengths);
        start = end;
    }
    lengths.sort_unstable_by(|a, b| b.cmp(a));
    let (nx, lx) = rustyomestats::stats::compute_nl(&lengths);
    c.n50 = nx[1];
    c.l50 = lx[1];
    c
}

/// Human-readable differences between the two (empty when they agree).
pub fn differences(s: &BasicCounts, r: &BasicCounts) -> Vec<String> {
    let mut d = Vec::new();
    let mut cmp = |what: &str, a: usize, b: usize| {
        if a != b {
            d.push(format!("{what}: sentinel-x {a} vs rustyomestats {b}"));
        }
    };
    cmp("sequences", s.seqs, r.seqs);
    cmp("total bp", s.total, r.total);
    cmp("G+C", s.gc, r.gc);
    cmp("A+T", s.at, r.at);
    cmp("N", s.n, r.n);
    cmp("N50", s.n50, r.n50);
    cmp("L50", s.l50, r.l50);
    d
}

/// Run the cross-check for one input file and log the outcome.
pub fn crosscheck(label: &str, contigs: &[Contig]) {
    let s = sentinel_x_counts(contigs);
    let r = rustyomestats_counts(contigs);
    let d = differences(&s, &r);
    let gc = if s.gc + s.at > 0 { 100.0 * s.gc as f64 / (s.gc + s.at) as f64 } else { 0.0 };
    if d.is_empty() {
        eprintln!(
            "[sentinel-x] rustyomestats cross-check {label}: agrees ({} sequence{}, {} bp, GC {gc:.2}%, N50 {} / L50 {})",
            s.seqs,
            if s.seqs == 1 { "" } else { "s" },
            s.total,
            s.n50,
            s.l50
        );
    } else {
        let note = if (s.n50, s.l50) != (r.n50, r.l50) && s.total % 2 == 1 {
            " (rustyomestats floors the 50% threshold for an odd total length)"
        } else {
            ""
        };
        eprintln!("[sentinel-x] rustyomestats cross-check {label}: DIFFERS — {}{note}", d.join("; "));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn contigs(lens: &[usize]) -> Vec<Contig> {
        lens.iter()
            .enumerate()
            .map(|(i, &l)| Contig {
                id: format!("c{i}"),
                header: format!("c{i}"),
                seq: (0..l).map(|j| b"ACGTN"[(i + j) % 5]).collect(),
            })
            .collect()
    }

    #[test]
    fn agrees_on_even_totals() {
        let cs = contigs(&[1000, 800, 600, 400]);
        let (s, r) = (sentinel_x_counts(&cs), rustyomestats_counts(&cs));
        assert_eq!(s, r);
        assert!(differences(&s, &r).is_empty());
    }

    #[test]
    fn documents_the_odd_total_n50_difference() {
        let cs = contigs(&[4, 3, 2]);
        let (s, r) = (sentinel_x_counts(&cs), rustyomestats_counts(&cs));
        assert_eq!((s.n50, s.l50), (3, 2));
        assert_eq!((r.n50, r.l50), (4, 1));
        assert_eq!((s.seqs, s.total, s.gc, s.at, s.n), (r.seqs, r.total, r.gc, r.at, r.n));
        assert_eq!(differences(&s, &r).len(), 2);
    }
}
