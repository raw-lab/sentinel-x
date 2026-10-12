//! `protsearch` — a DIAMOND-compatible `blastp` in Rust, used wherever the
//! reference tools call DIAMOND: CheckV's AAI search and CheckM2's KO
//! annotation. It reimplements DIAMOND 2.1.9's default-sensitivity code
//! path so the same hits, scores and statistics come out:
//!
//! * **Masking** — tantan repeat masking of queries and targets
//!   ([`mask::tantan_mask`]); motif soft-masking of seeds
//!   ([`mask::motif_soft_mask`]).
//! * **Seeds** — two weight-10 spaced shapes over Murphy's 10-letter
//!   alphabet with DIAMOND's seed-complexity cut ([`seeds`]); queries are
//!   indexed in memory, the database is streamed in blocks.
//! * **Seed verification** — the 48-residue identity fingerprint (≥ 11) and
//!   the ±48 ungapped window score against the 10⁴-E-value cutoff, reading
//!   across sequence boundaries as DIAMOND's sequence sets do ([`align`]).
//! * **Extension** — X-drop ungapped segments, then banded Smith-Waterman
//!   (BLOSUM62 as DIAMOND ships it, gaps 11/1, `banded-fast` band widths)
//!   with the Hauser composition-bias correction (`--comp-based-stats 1`,
//!   [`bias`]) and DIAMOND's traceback tie rules.
//! * **Statistics** — bit scores and Spouge finite-size E-values exactly as
//!   DIAMOND computes them ([`stats`]).
//! * **Reporting** — `--evalue`, `--id`, `--query-cover`, `--subject-cover`
//!   filter before `--max-target-seqs`; ranking by E-value, then score.
//!
//! Both DIAMOND seeding regimes are reproduced ([`Algo`]): query-indexed
//! (lazy target masking; DIAMOND's choice for a genome against a large
//! database) and double-indexed. On 4,869 simulated queries against 227,425
//! proteins, CheckM2's options give the same top hit and raw score as
//! DIAMOND 2.1.9 for 4,658/4,658 queries; with CheckV's options 99.99% of
//! DIAMOND's 72,638 hits are reproduced with identical scores.
//!
//! References are protein FASTA (plain/gzip) or DIAMOND `.dmnd` files read
//! in place ([`db`]); no DIAMOND binary is involved.

// Engine API: some accessors and record fields exist for completeness (BLAST tabular fields, alternative sources) and are used by tests only.
#![allow(dead_code)]

pub mod align;
pub mod alphabet;
pub mod bias;
pub mod db;
pub mod mask;
pub mod seeds;
pub mod stats;

use std::path::Path;

use anyhow::Result;
use rayon::prelude::*;

use align::{banded_align, banded_score, fingerprint_ok, set_letter, window_ungapped, xdrop_bias, Band, Profile};
use db::{Block, SeqSource};
use seeds::{for_each_seed, Shape, ShapeIndex, DEFAULT_SHAPES};

/// Search options (DIAMOND option names in comments).
#[derive(Clone, Debug)]
pub struct Params {
    /// `--evalue` (maximum reported E-value)
    pub evalue: f64,
    /// `--id` (minimum percent identity, 0..100)
    pub min_identity: f64,
    /// `--query-cover` (percent of the query covered by the alignment)
    pub query_cover: f64,
    /// `--subject-cover` (percent of the subject covered)
    pub subject_cover: f64,
    /// `--max-target-seqs` / `-k` (0 = unlimited)
    pub max_target_seqs: usize,
    /// `--masking tantan` (queries and targets)
    pub masking: bool,
    /// `--comp-based-stats 1` (Hauser composition bias correction)
    pub comp_based_stats: bool,
    /// DIAMOND's banded extension (`--ext banded-fast`); `false` = full SW.
    pub banded: bool,
    /// Seeding regime, as DIAMOND's `--algo` (auto = DIAMOND's own choice).
    pub algo: Algo,
    /// Reference residues per block (memory bound, cf. `--block-size`).
    pub block_letters: usize,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            evalue: 0.001,
            min_identity: 0.0,
            query_cover: 0.0,
            subject_cover: 0.0,
            max_target_seqs: 25,
            masking: true,
            comp_based_stats: true,
            banded: true,
            algo: Algo::Auto,
            block_letters: 400_000_000,
        }
    }
}

impl Params {
    /// CheckM2 1.1.0: `--max-target-seqs 1 --query-cover 80 --subject-cover 80
    /// --id 30 --evalue 1e-05`.
    pub fn checkm2() -> Params {
        Params {
            evalue: 1e-5,
            min_identity: 30.0,
            query_cover: 80.0,
            subject_cover: 80.0,
            max_target_seqs: 1,
            ..Default::default()
        }
    }
    /// CheckV 1.1.1: `--evalue 1.0e-5 --query-cover 50 --subject-cover 50 -k 10000`.
    pub fn checkv() -> Params {
        Params { evalue: 1e-5, query_cover: 50.0, subject_cover: 50.0, max_target_seqs: 10_000, ..Default::default() }
    }

    /// DIAMOND's `min_length_ratio` (run/config.cpp): with equal query and
    /// subject coverage thresholds, `max(cover/100 − 0.05, 0)`. When it is
    /// positive DIAMOND only seeds query/target pairs whose length ratio is
    /// at least this in both directions, and it processes (and prints) the
    /// queries of each block sorted by length.
    pub fn min_length_ratio(&self) -> f64 {
        if self.query_cover == self.subject_cover {
            (self.query_cover / 100.0 - 0.05).max(0.0)
        } else {
            0.0
        }
    }
}

/// DIAMOND's two seeding regimes, which differ in what the targets look like
/// while seeds are found and verified.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Algo {
    /// DIAMOND's automatic choice: query-indexed when the queries have at
    /// most 32 Mi residues and the database at least 256 Mi residues.
    Auto,
    /// Targets are seeded unmasked and tantan-masked lazily once they have
    /// a verified seed hit (DIAMOND `--algo query-indexed`).
    QueryIndexed,
    /// Targets are tantan-masked before seeding and motif-masked for seed
    /// enumeration (DIAMOND `--algo 0`).
    DoubleIndexed,
}

/// One reported alignment (BLAST tabular fields).
#[derive(Clone, Debug)]
pub struct Hit {
    /// Index of the query in the input slice.
    pub query: usize,
    pub subject: String,
    /// Index of the subject in the database (file order).
    pub subject_index: u64,
    pub qlen: usize,
    pub slen: usize,
    pub pident: f64,
    pub length: usize,
    pub mismatch: usize,
    pub gapopen: usize,
    pub qstart: usize,
    pub qend: usize,
    pub sstart: usize,
    pub send: usize,
    pub evalue: f64,
    pub bitscore: f64,
    pub raw: i32,
}

impl Hit {
    /// `--outfmt 6` line (qseqid sseqid pident length mismatch gapopen qstart
    /// qend sstart send evalue bitscore), formatted as DIAMOND prints it.
    pub fn outfmt6(&self, qid: &str) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            qid,
            self.subject,
            diamond_double(self.pident),
            self.length,
            self.mismatch,
            self.gapopen,
            self.qstart,
            self.qend,
            self.sstart,
            self.send,
            diamond_evalue(self.evalue),
            diamond_double(self.bitscore)
        )
    }

    /// `pident` as a parser of DIAMOND's table sees it (CheckV's AAI).
    pub fn printed_pident(&self) -> f64 {
        diamond_double(self.pident).parse().unwrap_or(self.pident)
    }

    /// `bitscore` as a parser of DIAMOND's table sees it.
    pub fn printed_bitscore(&self) -> f64 {
        diamond_double(self.bitscore).parse().unwrap_or(self.bitscore)
    }
}

/// DIAMOND's `Util::String::format_double`: values ≥ 100 truncated to an
/// integer ("for keeping output compatible with BLAST"), others printed
/// with one decimal after `llround(x·10)` (ties away from zero).
pub fn diamond_double(x: f64) -> String {
    if x >= 100.0 {
        return format!("{}", x.floor() as i64);
    }
    let i = (x * 10.0).round() as i64;
    let (q, r) = (i / 10, i % 10);
    if i < 0 && q == 0 {
        format!("-0.{}", r.abs())
    } else {
        format!("{}.{}", q, r.abs())
    }
}

/// DIAMOND's `print_e`: `"0.0"` for zero, else C `%.2e`.
pub fn diamond_evalue(x: f64) -> String {
    if x == 0.0 {
        "0.0".into()
    } else {
        crate::common::pyfmt::c_fmt_e(x, 2)
    }
}

/// Search statistics.
#[derive(Clone, Debug, Default)]
pub struct Stats {
    pub db_letters: u64,
    pub db_sequences: u64,
    pub blocks: usize,
    pub seed_hits: u64,
    pub verified_hits: u64,
    pub candidates: u64,
    pub extensions: u64,
    pub hits: u64,
    pub masked_query_letters: u64,
    /// Seeding regime actually used.
    pub algo: Option<Algo>,
}

struct Prepared {
    seqs: Vec<Vec<u8>>,
    bias: Vec<Vec<i8>>,
    profiles: Vec<Profile>,
    shapes: Vec<ShapeIndex>,
    ug_cutoff: Vec<i32>,
}

/// DIAMOND's stage-2 ungapped cutoff (`CutoffTable` at `ungapped_evalue`
/// 10,000 normalised to 10⁹ letters, query length rounded down to a power
/// of two; 25 bits for queries ≤ 60 residues). Hits need a window score
/// strictly above it.
pub fn ungapped_cutoff(qlen: usize) -> i32 {
    let raw = |bits: f64| -> i32 { ((bits * std::f64::consts::LN_2 + stats::K.ln()) / stats::LAMBDA).ceil() as i32 };
    if qlen <= 60 {
        return raw(25.0);
    }
    let l2 = 1u64 << (63 - (qlen as u64).leading_zeros());
    raw(-(10_000.0f64 / 1e9 / l2 as f64).log2())
}

/// DIAMOND `--ext banded-fast` band half-width for a query length.
pub fn band_for(qlen: usize) -> isize {
    match qlen {
        0..=49 => 12,
        50..=99 => 16,
        100..=249 => 30,
        250..=349 => 40,
        _ => 64,
    }
}

/// Search `queries` (ASCII protein sequences) against a reference database
/// file (`.dmnd` or FASTA).
pub fn search_path(queries: &[(String, Vec<u8>)], db: &Path, p: &Params) -> Result<(Vec<Hit>, Stats)> {
    let mut src = db::open(db)?;
    search(queries, src.as_mut(), p)
}

/// Search `queries` (ASCII protein sequences) against any [`SeqSource`].
/// Hits are returned grouped by query (input order), best first.
pub fn search(queries: &[(String, Vec<u8>)], src: &mut dyn SeqSource, p: &Params) -> Result<(Vec<Hit>, Stats)> {
    let db_letters = src.letters();
    let mut st = Stats { db_letters, db_sequences: src.sequences().unwrap_or(0), ..Default::default() };
    let masked: Vec<(Vec<u8>, usize)> = queries
        .par_iter()
        .map(|q| {
            let mut e = alphabet::encode(&q.1);
            let n = if p.masking { mask::tantan_mask(&mut e) } else { 0 };
            (e, n)
        })
        .collect();
    st.masked_query_letters = masked.iter().map(|m| m.1 as u64).sum();
    let seqs: Vec<Vec<u8>> = masked.into_iter().map(|m| m.0).collect();
    let bias: Vec<Vec<i8>> =
        seqs.par_iter().map(|q| if p.comp_based_stats { bias::hauser_int8(q) } else { vec![0i8; q.len()] }).collect();
    let profiles: Vec<Profile> = seqs
        .par_iter()
        .zip(bias.par_iter())
        .map(|(q, b)| Profile::new(q, if p.comp_based_stats { Some(b.as_slice()) } else { None }))
        .collect();
    // Query seeds are enumerated on the motif-soft-masked queries (DIAMOND
    // restores the motif letters after seed enumeration).
    let seed_queries: Vec<Vec<u8>> = seqs.par_iter().map(|q| mask::motif_soft_mask(q)).collect();
    let shapes: Vec<ShapeIndex> = DEFAULT_SHAPES.par_iter().map(|c| ShapeIndex::build(Shape::parse(c), &seed_queries)).collect();
    drop(seed_queries);
    let ug_cutoff: Vec<i32> = seqs.iter().map(|q| ungapped_cutoff(q.len())).collect();
    let query_letters: u64 = seqs.iter().map(|q| q.len() as u64).sum();
    let prep = Prepared { seqs, bias, profiles, shapes, ug_cutoff };
    let algo = match p.algo {
        Algo::Auto => {
            if query_letters <= 32 << 20 && db_letters >= 256 << 20 {
                Algo::QueryIndexed
            } else {
                Algo::DoubleIndexed
            }
        }
        a => a,
    };
    st.algo = Some(algo);

    let k = if p.max_target_seqs == 0 { usize::MAX } else { p.max_target_seqs };
    let mut results: Vec<Vec<Hit>> = (0..queries.len()).map(|_| Vec::new()).collect();
    let mut block = Block::default();
    while src.next_block(p.block_letters.max(1), &mut block)? {
        st.blocks += 1;
        // Double-indexed: targets are tantan-masked before seeding.
        let masked_block = if algo == Algo::DoubleIndexed && p.masking { Some(mask_block(&block)) } else { None };
        let vblock = masked_block.as_ref().unwrap_or(&block);
        let (hits, seed_hits) = scan_block(&prep, vblock, algo == Algo::DoubleIndexed, p.min_length_ratio());
        st.seed_hits += seed_hits;
        st.verified_hits += hits.len() as u64;
        // Query-indexed: lazily mask the targets that have verified seed hits.
        let mut subj: Vec<u32> = hits.iter().map(|h| h.subject).collect();
        subj.sort_unstable();
        subj.dedup();
        let lazy: Vec<Vec<u8>> = if masked_block.is_none() {
            subj.par_iter()
                .map(|&s| {
                    let mut t = block.get(s as usize).to_vec();
                    if p.masking {
                        mask::tantan_mask(&mut t);
                    }
                    t
                })
                .collect()
        } else {
            Vec::new()
        };
        let target_of = |s: u32| -> &[u8] {
            match &masked_block {
                Some(b) => b.get(s as usize),
                None => &lazy[subj.binary_search(&s).expect("masked target")],
            }
        };
        // group verified seed hits by query, then subject
        let mut by_q: Vec<Vec<SeedHit>> = (0..queries.len()).map(|_| Vec::new()).collect();
        for h in hits {
            by_q[h.query as usize].push(h);
        }
        let block_hits: Vec<(usize, Vec<Hit>, u64, u64)> = by_q
            .into_par_iter()
            .enumerate()
            .filter(|(_, v)| !v.is_empty())
            .map(|(qi, mut v)| {
                v.sort_unstable_by_key(|h| (h.subject, h.qpos as i64 - h.spos as i64, h.spos));
                let (h, cands, ext) = extend_query(&prep, qi, &v, &block, &target_of, db_letters, p, k);
                (qi, h, cands, ext)
            })
            .collect();
        for (qi, hs, cands, ext) in block_hits {
            st.candidates += cands;
            st.extensions += ext;
            let r = &mut results[qi];
            r.extend(hs);
            if r.len() > k {
                sort_hits(r);
                r.truncate(k);
            }
        }
    }
    for r in results.iter_mut() {
        sort_hits(r);
        if r.len() > k {
            r.truncate(k);
        }
    }
    let mut out = Vec::new();
    for qi in output_order(&prep.seqs, p) {
        out.append(&mut results[qi]);
    }
    st.hits = out.len() as u64;
    Ok((out, st))
}

/// The order in which DIAMOND prints queries: input order, except that with
/// a positive `min_length_ratio` each query block (input order, closed once
/// it reaches the block size) is sorted by length, longest first, ties to
/// the later query (`Block::length_sorted`, `greater<pair<Loc, BlockId>>`).
fn output_order(seqs: &[Vec<u8>], p: &Params) -> Vec<usize> {
    let n = seqs.len();
    if p.min_length_ratio() <= 0.0 {
        return (0..n).collect();
    }
    let mut order = Vec::with_capacity(n);
    let mut start = 0;
    while start < n {
        let mut end = start;
        let mut letters = 0usize;
        // DIAMOND loads query blocks of `--block-size` (default 2.0) billion letters
        while end < n && letters < 2_000_000_000 {
            letters += seqs[end].len();
            end += 1;
        }
        let mut blk: Vec<usize> = (start..end).collect();
        blk.sort_by(|&a, &b| (seqs[b].len(), b).cmp(&(seqs[a].len(), a)));
        order.extend(blk);
        start = end;
    }
    order
}

/// DIAMOND's ranking (`Match::cmp_evalue`): E-value ascending, then raw
/// score descending (equal raw scores therefore favour the shorter target);
/// complete ties, which DIAMOND leaves to an unstable sort, go to the later
/// database entry as DIAMOND 2.1.9 was observed to do.
fn rank(e1: f64, r1: i32, i1: u64, e2: f64, r2: i32, i2: u64) -> std::cmp::Ordering {
    e1.partial_cmp(&e2).unwrap_or(std::cmp::Ordering::Equal).then(r2.cmp(&r1)).then(i2.cmp(&i1))
}

fn sort_hits(v: &mut [Hit]) {
    v.sort_by(|a, b| rank(a.evalue, a.raw, a.subject_index, b.evalue, b.raw, b.subject_index));
}

#[derive(Clone, Copy, Debug)]
struct SeedHit {
    query: u32,
    subject: u32,
    qpos: u32,
    spos: u32,
}

/// Tantan-mask every sequence of a block (double-indexed regime).
fn mask_block(block: &Block) -> Block {
    let parts: Vec<Vec<u8>> = (0..block.len())
        .into_par_iter()
        .map(|i| {
            let mut t = block.get(i).to_vec();
            mask::tantan_mask(&mut t);
            t
        })
        .collect();
    let mut seq = Vec::with_capacity(block.seq.len());
    for p in &parts {
        seq.extend_from_slice(p);
    }
    Block { first: block.first, ids: Vec::new(), seq, off: block.off.clone() }
}

/// Seed stage over one block: seed matches, DIAMOND's fingerprint (stage 1)
/// and ungapped-window (stage 2) verification. `motif` soft-masks target
/// motifs for seed enumeration (double-indexed regime).
fn scan_block(prep: &Prepared, block: &Block, motif: bool, mlr: f64) -> (Vec<SeedHit>, u64) {
    let n = block.len();
    let chunk = 256usize;
    let parts: Vec<(Vec<SeedHit>, u64)> = (0..n.div_ceil(chunk))
        .into_par_iter()
        .map(|c| {
            let mut out = Vec::new();
            let mut nseed = 0u64;
            for si in c * chunk..((c + 1) * chunk).min(n) {
                let s = block.get(si);
                let nb = block.len();
                let sx = |p: isize| set_letter(|k| block.get(k), nb, si, p);
                let seed_view =
                    if motif { std::borrow::Cow::Owned(mask::motif_soft_mask(s)) } else { std::borrow::Cow::Borrowed(s) };
                for idx in &prep.shapes {
                    for_each_seed(&seed_view, &idx.shape, false, |sp, key| {
                        for &(q, qp) in idx.lookup(key) {
                            nseed += 1;
                            let qs = &prep.seqs[q as usize];
                            // stage1_mutual_cov: pairs whose lengths differ too much are never seeded
                            if mlr > 0.0 {
                                let (ql, sl) = (qs.len() as f64, s.len() as f64);
                                if ql / sl < mlr || sl / ql < mlr {
                                    continue;
                                }
                            }
                            let nq = prep.seqs.len();
                            let qx = |p: isize| set_letter(|k| prep.seqs[k].as_slice(), nq, q as usize, p);
                            if !fingerprint_ok(qs, qp as usize, s, sp, qx, sx) {
                                continue;
                            }
                            if window_ungapped(qs, qp as usize, s, sp, sx) <= prep.ug_cutoff[q as usize] {
                                continue;
                            }
                            out.push(SeedHit { query: q, subject: si as u32, qpos: qp, spos: sp as u32 });
                        }
                    });
                }
            }
            (out, nseed)
        })
        .collect();
    let mut all = Vec::new();
    let mut seeds = 0;
    for (v, s) in parts {
        seeds += s;
        all.extend(v);
    }
    (all, seeds)
}

/// Ungapped + gapped stages for one query against its candidate targets in a
/// block (hits sorted by subject, diagonal, subject position).
#[allow(clippy::too_many_arguments)]
fn extend_query<'a>(
    prep: &Prepared,
    qi: usize,
    hits: &[SeedHit],
    block: &Block,
    target_of: &(dyn Fn(u32) -> &'a [u8] + Sync),
    db_letters: u64,
    p: &Params,
    k: usize,
) -> (Vec<Hit>, u64, u64) {
    let q = &prep.seqs[qi];
    let bias = &prep.bias[qi];
    let prof = &prep.profiles[qi];
    let m = q.len();
    let bw = band_for(m);
    // (raw score, subject, band, E-value)
    let mut scored: Vec<(i32, u32, Band, f64)> = Vec::new();
    let mut ncand = 0u64;
    let mut i = 0;
    while i < hits.len() {
        let s_idx = hits[i].subject;
        let mut j = i;
        while j < hits.len() && hits[j].subject == s_idx {
            j += 1;
        }
        let t = target_of(s_idx);
        // ungapped stage (DIAMOND align/ungapped.cpp): one X-drop extension
        // per seed hit unless covered by the previous segment on its diagonal
        let mut segs: Vec<align::Segment> = Vec::new();
        for h in &hits[i..j] {
            let d = h.qpos as isize - h.spos as isize;
            if let Some(last) = segs.last() {
                if last.diag() == d && last.subject_end() >= h.spos as usize {
                    continue;
                }
            }
            let seg = xdrop_bias(q, bias, t, h.qpos as usize, h.spos as usize);
            if seg.score > 0 {
                segs.push(seg);
            }
        }
        i = j;
        if segs.is_empty() {
            continue;
        }
        ncand += 1;
        let slen = t.len() as isize;
        let band = if p.banded {
            let dmin = segs.iter().map(|s| s.diag()).min().unwrap();
            let dmax = segs.iter().map(|s| s.diag()).max().unwrap();
            Band { d0: (dmin - bw).max(-(slen - 1)), d1: (dmax + 1 + bw).min(m as isize) }
        } else {
            Band::full(m, t.len())
        };
        let raw = banded_score(prof, t, band);
        if raw <= 0 {
            continue;
        }
        let ev = stats::evalue(raw, m, t.len(), db_letters);
        if ev <= p.evalue {
            scored.push((raw, s_idx, band, ev));
        }
    }
    scored.sort_by(|a, b| rank(a.3, a.0, a.1 as u64, b.3, b.0, b.1 as u64));
    let mut out = Vec::new();
    let mut ext = 0u64;
    for &(raw, s, band, ev) in &scored {
        if out.len() >= k {
            break;
        }
        let t = target_of(s);
        let al = banded_align(q, prof, t, band);
        ext += 1;
        let len = al.length.max(1);
        let pident = 100.0 * al.identities as f64 / len as f64;
        if pident < p.min_identity {
            continue;
        }
        let qcov = 100.0 * (al.qend + 1 - al.qstart) as f64 / m.max(1) as f64;
        let scov = 100.0 * (al.send + 1 - al.sstart) as f64 / t.len().max(1) as f64;
        if qcov < p.query_cover || scov < p.subject_cover {
            continue;
        }
        out.push(Hit {
            query: qi,
            subject: block.ids[s as usize].clone(),
            subject_index: block.first + s as u64,
            qlen: m,
            slen: t.len(),
            pident,
            length: al.length,
            mismatch: al.mismatches,
            gapopen: al.gap_opens,
            qstart: al.qstart,
            qend: al.qend,
            sstart: al.sstart,
            send: al.send,
            evalue: ev,
            bitscore: stats::bitscore(raw),
            raw,
        });
    }
    (out, ncand, ext)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Substitute ~`rate` of the residues at pseudo-random positions (a
    /// periodic pattern would be masked by tantan as a repeat).
    fn mutate(s: &[u8], rate: f64, seed: u64) -> Vec<u8> {
        let mut x = seed;
        s.iter()
            .map(|&c| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                if ((x >> 33) % 10_000) as f64 / 10_000.0 < rate {
                    b"ACDEFGHIKLMNPQRSTVWY"[((x >> 20) % 20) as usize]
                } else {
                    c
                }
            })
            .collect()
    }

    const Q: &[u8] = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQCGVWPAACRESVPALL";

    #[test]
    fn finds_homologs_and_rejects_random() {
        let hom = mutate(Q, 0.2, 11);
        let far = mutate(Q, 0.4, 12);
        let mut x = 7u64;
        let random: Vec<u8> = (0..Q.len())
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACDEFGHIKLMNPQRSTVWY"[((x >> 33) % 20) as usize]
            })
            .collect();
        let db = vec![("rnd".to_string(), random), ("hom".to_string(), hom), ("far".to_string(), far)];
        let queries = vec![("q".to_string(), Q.to_vec())];
        let mut src = db::MemSource::new(db);
        let (hits, st) = search(&queries, &mut src, &Params { evalue: 1e-5, ..Default::default() }).unwrap();
        assert!(st.seed_hits > 0);
        let names: Vec<&str> = hits.iter().map(|h| h.subject.as_str()).collect();
        assert_eq!(names[0], "hom");
        assert!(names.contains(&"far"));
        assert!(!names.contains(&"rnd"));
        let h = &hits[0];
        assert!(h.pident > 78.0 && h.pident < 88.0, "{}", h.pident);
        assert!(h.qstart <= 3 && h.qend + 3 >= Q.len(), "{} {}", h.qstart, h.qend);
    }

    #[test]
    fn checkm2_filters_and_top1() {
        // A short local match covering < 80% of the query must be filtered,
        // leaving the full-length (lower-scoring) homolog as the top hit.
        let part = Q[..120].to_vec();
        let full = mutate(Q, 0.25, 13);
        let db = vec![("part".to_string(), part), ("full".to_string(), full)];
        let queries = vec![("q".to_string(), Q.to_vec())];
        let mut src = db::MemSource::new(db);
        let (hits, _) = search(&queries, &mut src, &Params::checkm2()).unwrap();
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].subject, "full");
    }

    /// DIAMOND 2.1.9 parity on a fixture of 16 simulated queries against
    /// 1,100 proteins (`testdata/`, outputs from `diamond blastp` with
    /// CheckM2's and CheckV's options).
    mod diamond_parity {
        use super::super::*;
        use std::collections::HashMap;

        fn read_fasta(txt: &str) -> Vec<(String, Vec<u8>)> {
            let mut v: Vec<(String, Vec<u8>)> = Vec::new();
            for l in txt.lines() {
                if let Some(h) = l.strip_prefix('>') {
                    v.push((h.to_string(), Vec::new()));
                } else if let Some(last) = v.last_mut() {
                    last.1.extend_from_slice(l.trim().as_bytes());
                }
            }
            v
        }

        fn run(params: &Params, expected: &str) -> (usize, usize, usize) {
            let queries = read_fasta(include_str!("testdata/queries.faa"));
            let db = read_fasta(include_str!("testdata/db.faa"));
            let mut src = db::MemSource::new(db);
            let (hits, _) = search(&queries, &mut src, params).unwrap();
            let mine: HashMap<(String, String), i32> =
                hits.iter().map(|h| ((queries[h.query].0.clone(), h.subject.clone()), h.raw)).collect();
            let (mut same, mut total) = (0, 0);
            for l in expected.lines() {
                let f: Vec<&str> = l.split('\t').collect();
                total += 1;
                if mine.get(&(f[0].to_string(), f[1].to_string())) == Some(&f[12].parse::<i32>().unwrap()) {
                    same += 1;
                }
            }
            (same, total, mine.len())
        }

        #[test]
        fn checkm2_options() {
            let (same, total, n) = run(&Params::checkm2(), include_str!("testdata/diamond_checkm2.tsv"));
            assert_eq!((same, n), (total, total), "top hits differ from DIAMOND");
        }

        #[test]
        fn checkv_options() {
            let (same, total, n) = run(&Params::checkv(), include_str!("testdata/diamond_checkv.tsv"));
            assert!(same * 100 >= total * 99, "{same}/{total} DIAMOND hits reproduced");
            assert!(n <= total + total / 50, "{n} hits vs DIAMOND's {total}");
        }
    }
}
