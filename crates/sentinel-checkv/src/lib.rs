//! # sentinel-checkv
//!
//! Viral quality assessment in the spirit of **CheckV**, split into the parts
//! that need a reference database and the parts that do not:
//!
//! * **Terminal repeats** ([`find_dtr`], [`find_itr`]) — direct/inverted
//!   terminal repeats are strong, *database-free* evidence of a complete genome.
//! * **Provirus / host contamination** ([`find_host_boundary`]) — a sliding
//!   gene-window scan combining viral/host annotation with a GC shift to locate
//!   the host–virus junction in an integrated provirus.
//! * **Completeness** ([`hmm_completeness`]) — when an HMM → reference-length
//!   table is supplied, the expected genome length gives a completeness
//!   estimate (the CheckV HMM-based path). Without it, completeness is left
//!   undetermined and the verdict rests on terminal-repeat evidence.
//! * **Quality tiering** ([`quality_call`]) — folds the above into the standard
//!   Complete / High / Medium / Low tiers.

use std::collections::HashMap;

use sentinel_core::composition::{gc_content, kmer_freq};
use sentinel_core::{Contig, Domain, GeneSet, QcRecord, Quality};

/// Reverse complement of a nucleotide slice (`N` preserved).
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'G' | b'g' => b'C',
            b'C' | b'c' => b'G',
            other => other,
        })
        .collect()
}

/// Kind of terminal repeat found.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TrKind {
    /// Direct terminal repeat (same sequence at both ends) — typical of
    /// complete, circularly-permuted / headful phage genomes.
    Dtr,
    /// Inverted terminal repeat (end is reverse-complement of start).
    Itr,
    None,
}

/// A detected terminal repeat.
#[derive(Clone, Debug)]
pub struct TerminalRepeat {
    pub kind: TrKind,
    pub length: usize,
}

/// Find a **direct terminal repeat**: the contig begins and ends with the same
/// sequence (≥ `min_length`). Mirrors CheckV's `fetch_dtr`.
pub fn find_dtr(seq: &[u8], min_length: usize) -> TerminalRepeat {
    if seq.len() < 2 * min_length {
        return TerminalRepeat { kind: TrKind::None, length: 0 };
    }
    let startseq = &seq[..min_length];
    let half = seq.len() / 2;
    // Find every start position of `startseq` within the 2nd half of the contig.
    let mut i = half;
    while i + min_length <= seq.len() {
        if &seq[i..i + min_length] == startseq {
            // does the match extend to the contig end, repeating the prefix?
            let endseq = &seq[i..];
            if endseq.len() <= seq.len() && &seq[..endseq.len()] == endseq {
                return TerminalRepeat { kind: TrKind::Dtr, length: endseq.len() };
            }
        }
        i += 1;
    }
    TerminalRepeat { kind: TrKind::None, length: 0 }
}

/// Find an **inverted terminal repeat**: the contig's prefix equals the
/// reverse-complement's prefix. Mirrors CheckV's `fetch_itr`.
pub fn find_itr(seq: &[u8], min_len: usize, max_len: usize) -> TerminalRepeat {
    if seq.len() < 2 * min_len {
        return TerminalRepeat { kind: TrKind::None, length: 0 };
    }
    let rev = reverse_complement(seq);
    if seq[..min_len] == rev[..min_len] {
        let mut i = min_len + 1;
        while i <= max_len && i <= seq.len() && seq[..i] == rev[..i] {
            i += 1;
        }
        return TerminalRepeat { kind: TrKind::Itr, length: i - 1 };
    }
    TerminalRepeat { kind: TrKind::None, length: 0 }
}

/// Detect the best terminal repeat (prefers DTR, then ITR).
pub fn detect_terminal_repeat(seq: &[u8]) -> TerminalRepeat {
    let dtr = find_dtr(seq, 21);
    if dtr.kind == TrKind::Dtr {
        return dtr;
    }
    find_itr(seq, 21, 1000)
}

/// Per-gene annotation category used by the contamination scan.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GeneCat {
    /// Annotated by a viral marker (CheckV `+1`).
    Viral,
    /// Annotated by a microbial/host marker (CheckV `-1`).
    Host,
    /// No confident annotation (`0`).
    Unknown,
}

impl GeneCat {
    fn score(self) -> i32 {
        match self {
            GeneCat::Viral => 1,
            GeneCat::Host => -1,
            GeneCat::Unknown => 0,
        }
    }
}

/// A host/viral partition of a contig.
#[derive(Clone, Debug)]
pub struct Boundary {
    /// Gene index of the breakpoint (genes `[0, idx)` vs `[idx, n)`).
    pub gene_index: usize,
    /// Number of viral / host genes either side.
    pub viral_left: usize,
    pub host_left: usize,
    pub viral_right: usize,
    pub host_right: usize,
    /// Combined annotation+GC score of the breakpoint.
    pub score: f64,
}

/// Find the single best host–virus breakpoint by scanning gene windows and
/// combining the viral/host annotation delta with the GC-content shift, after
/// CheckV's contamination scoring. `gc_per_gene` is the GC fraction of each
/// gene's nucleotide span; `gc_weight` scales the GC contribution.
///
/// Returns `None` when no breakpoint separates a meaningfully host-enriched
/// region from a viral one (i.e. the contig looks purely viral or purely host).
pub fn find_host_boundary(
    cats: &[GeneCat],
    gc_per_gene: &[f64],
    window: usize,
    gc_weight: f64,
    min_host_genes: usize,
) -> Option<Boundary> {
    let n = cats.len();
    if n < 2 {
        return None;
    }
    let w = window.max(1);
    let mut best: Option<Boundary> = None;

    for idx in 1..n {
        let l0 = idx.saturating_sub(w);
        let l1 = idx; // left window [l0, idx)
        let r0 = idx;
        let r1 = (idx + w).min(n); // right window [idx, r1)

        let (mut vl, mut hl) = (0usize, 0usize);
        let mut ann_left = 0i32;
        let mut gc_left = 0.0;
        for i in l0..l1 {
            ann_left += cats[i].score();
            gc_left += gc_per_gene.get(i).copied().unwrap_or(0.5);
            match cats[i] {
                GeneCat::Viral => vl += 1,
                GeneCat::Host => hl += 1,
                GeneCat::Unknown => {}
            }
        }
        let (mut vr, mut hr) = (0usize, 0usize);
        let mut ann_right = 0i32;
        let mut gc_right = 0.0;
        for i in r0..r1 {
            ann_right += cats[i].score();
            gc_right += gc_per_gene.get(i).copied().unwrap_or(0.5);
            match cats[i] {
                GeneCat::Viral => vr += 1,
                GeneCat::Host => hr += 1,
                GeneCat::Unknown => {}
            }
        }
        let nl = (l1 - l0).max(1) as f64;
        let nr = (r1 - r0).max(1) as f64;
        let ann_delta = (ann_right as f64 / nr) - (ann_left as f64 / nl);
        let gc_delta = ((gc_right / nr) - (gc_left / nl)).abs();
        // A host→viral or viral→host transition shows a large |annotation delta|
        // reinforced by a GC shift.
        let score = ann_delta.abs() + gc_weight * gc_delta;

        // Require one side to be host-enriched.
        let host_side_ok = hl >= min_host_genes || hr >= min_host_genes;
        if host_side_ok {
            let cand = Boundary {
                gene_index: idx,
                viral_left: vl,
                host_left: hl,
                viral_right: vr,
                host_right: hr,
                score,
            };
            if best.as_ref().map(|b| score > b.score).unwrap_or(true) {
                best = Some(cand);
            }
        }
    }
    best
}

/// HMM-based completeness: given the set of CheckV HMMs a contig's genes hit and
/// a table of the expected (median) host-genome length for each HMM, the
/// consensus expected length yields a completeness percentage. Returns `None`
/// when no hit HMM is in the table.
pub fn hmm_completeness(
    hit_hmms: &[String],
    contig_len: usize,
    ref_len_table: &HashMap<String, f64>,
) -> Option<f64> {
    let lens: Vec<f64> = hit_hmms
        .iter()
        .filter_map(|h| ref_len_table.get(h).copied())
        .collect();
    if lens.is_empty() {
        return None;
    }
    // Median expected length (robust to a few off-target HMMs).
    let mut sorted = lens.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let expected = sorted[sorted.len() / 2];
    if expected <= 0.0 {
        return None;
    }
    Some((contig_len as f64 / expected * 100.0).min(100.0))
}

/// Final quality tier from the available evidence.
pub fn quality_call(
    completeness: Option<f64>,
    contamination: f64,
    tr: &TerminalRepeat,
    kmer_dup: f64,
) -> (Quality, Vec<String>) {
    let mut warnings = Vec::new();
    if kmer_dup > 1.5 {
        warnings.push(format!("high k-mer redundancy ({kmer_dup:.2}x) — possible duplication"));
    }
    if contamination > 5.0 {
        warnings.push(format!("host contamination {contamination:.1}%"));
    }

    // A clean terminal repeat with low contamination ⇒ Complete.
    if tr.kind != TrKind::None && contamination < 5.0 && kmer_dup < 1.5 {
        return (Quality::Complete, warnings);
    }
    match completeness {
        Some(c) if c >= 90.0 => (Quality::High, warnings),
        Some(c) if c >= 50.0 => (Quality::Medium, warnings),
        Some(_) => (Quality::Low, warnings),
        None => {
            if tr.kind != TrKind::None {
                (Quality::High, warnings)
            } else {
                (Quality::NotDetermined, warnings)
            }
        }
    }
}

/// Convenience: assess one contig end-to-end (database-free parts only) given
/// its genes, per-gene annotation categories, and optional HMM completeness.
#[allow(clippy::too_many_arguments)]
pub fn assess_contig(
    contig: &Contig,
    gene_indices: &[usize],
    geneset: &GeneSet,
    cats: &[GeneCat],
    hit_hmms: &[String],
    ref_len_table: &HashMap<String, f64>,
) -> QcRecord {
    let tr = detect_terminal_repeat(&contig.seq);
    let kmer_dup = kmer_freq(&contig.seq, 15);

    // Per-gene GC for the boundary scan.
    let gc_per_gene: Vec<f64> = gene_indices
        .iter()
        .map(|&gi| {
            let g = &geneset.genes[gi];
            let (b, e) = (g.begin.min(g.end), g.begin.max(g.end));
            let s = &contig.seq[b.saturating_sub(1)..e.min(contig.seq.len())];
            gc_content(s)
        })
        .collect();

    let boundary = find_host_boundary(cats, &gc_per_gene, 40, 8.0, 2);
    let (viral_genes, host_genes) = cats.iter().fold((0usize, 0usize), |(v, h), c| match c {
        GeneCat::Viral => (v + 1, h),
        GeneCat::Host => (v, h + 1),
        GeneCat::Unknown => (v, h),
    });

    // Contamination ≈ fraction of the contig assigned to a host region. If a
    // boundary was found, estimate host fraction from gene counts on the
    // host-enriched side.
    let contamination = match &boundary {
        Some(b) => {
            let host_total = host_genes as f64;
            let denom = (viral_genes + host_genes).max(1) as f64;
            let frac = host_total / denom;
            // only count as contamination if a real boundary separated them
            if b.host_left.max(b.host_right) >= 2 {
                frac * 100.0
            } else {
                0.0
            }
        }
        None => 0.0,
    };

    let completeness = hmm_completeness(hit_hmms, contig.len(), ref_len_table);
    let (quality, mut warnings) = quality_call(completeness, contamination, &tr, kmer_dup);

    let mut rec = QcRecord::new(contig.id.clone(), contig.len(), Domain::Virus, "checkv");
    rec.gene_count = gene_indices.len();
    rec.completeness = completeness;
    rec.contamination = Some(contamination);
    rec.quality = quality;
    rec.method = if completeness.is_some() {
        "HMM-completeness".into()
    } else if tr.kind != TrKind::None {
        "terminal-repeat".into()
    } else {
        "structural".into()
    };
    if tr.kind != TrKind::None {
        rec = rec.with_extra(
            "terminal_repeat",
            format!("{:?}:{}bp", tr.kind, tr.length),
        );
    }
    rec = rec
        .with_extra("viral_genes", viral_genes.to_string())
        .with_extra("host_genes", host_genes.to_string());
    if let Some(b) = boundary {
        rec = rec.with_extra("host_boundary_gene", b.gene_index.to_string());
    }
    rec.warnings.append(&mut warnings);
    rec
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detects_direct_terminal_repeat() {
        // Build a contig that begins and ends with the same 30-mer.
        let repeat = b"ACGTACGTACGTACGTACGTACGTAAGGCC";
        let mut seq = Vec::new();
        seq.extend_from_slice(repeat);
        seq.extend_from_slice(b"TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGG");
        seq.extend_from_slice(repeat);
        let tr = find_dtr(&seq, 21);
        assert_eq!(tr.kind, TrKind::Dtr);
        assert!(tr.length >= 30);
    }

    #[test]
    fn detects_inverted_terminal_repeat() {
        // start with X; end with revcomp(X) so the contig prefix == revcomp prefix
        let arm = b"ACGTACGTACGTACGTACGTAC";
        let mut seq = Vec::new();
        seq.extend_from_slice(arm);
        seq.extend_from_slice(b"TTTTAGAGAGAGATTTTCCCCGGGGAAAA");
        seq.extend_from_slice(&reverse_complement(arm));
        let tr = find_itr(&seq, 21, 1000);
        assert_eq!(tr.kind, TrKind::Itr);
    }

    #[test]
    fn no_tr_on_random() {
        let seq = b"ACGTTGCAACGTTGCAACGATCGATCGATCGTAGCTAGCATCGGCTAGCTAGCTAGCAT";
        assert_eq!(detect_terminal_repeat(seq).kind, TrKind::None);
    }

    #[test]
    fn finds_provirus_boundary() {
        // host genes (-1) then viral genes (+1): boundary should be at the switch.
        let mut cats = vec![GeneCat::Host; 6];
        cats.extend(vec![GeneCat::Viral; 6]);
        let gc = vec![0.6; 12];
        let b = find_host_boundary(&cats, &gc, 40, 8.0, 2).unwrap();
        assert!(b.gene_index >= 4 && b.gene_index <= 8, "idx={}", b.gene_index);
        assert!(b.host_left >= 2);
    }

    #[test]
    fn completeness_from_ref_table() {
        let mut t = HashMap::new();
        t.insert("HMMa".to_string(), 50000.0);
        t.insert("HMMb".to_string(), 50000.0);
        let c = hmm_completeness(&["HMMa".into(), "HMMb".into()], 25000, &t).unwrap();
        assert!((c - 50.0).abs() < 1e-6);
    }
}
