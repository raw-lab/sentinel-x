//! Host–virus boundary detection (CheckV 1.1.1 `contamination.py`:
//! `compute_delta`, `define_regions`, `define_recursive_regions`).
//!
//! Genes carry a category (+1 viral, −1 microbial, 0 unannotated) and a GC
//! percentage. Breakpoints are scored between a cumulative left window
//! (from the previous breakpoint) and a right window of
//! `min(max(15, round(0.3·n)), 50)` genes; the procedure is re-applied to
//! each viral island up to five times and adjacent same-type regions are
//! merged.

use crate::common::pyfmt::{np_mean, round_ties_even};

/// Minimal per-gene view used by the boundary scan.
#[derive(Clone, Copy, Debug)]
pub struct GeneView {
    pub start: usize,
    pub end: usize,
    pub gc: f64,
    pub cat: i8,
}

/// A host or viral region (1-based inclusive coordinates; gene indexes are
/// 0-based, `end_gene` exclusive).
#[derive(Clone, Debug, PartialEq)]
pub struct Region {
    pub viral: bool,
    pub start_pos: usize,
    pub end_pos: usize,
    pub start_gene: usize,
    pub end_gene: usize,
    pub size: usize,
    pub length: usize,
    pub host_genes: usize,
    pub viral_genes: usize,
}

impl Region {
    pub fn type_label(&self) -> &'static str {
        if self.viral {
            "viral"
        } else {
            "host"
        }
    }
}

/// Scan parameters (CheckV defaults in [`RegionParams::default`]).
#[derive(Clone, Copy, Debug)]
pub struct RegionParams {
    pub gc_weight: f64,
    pub delta_cutoff: f64,
    pub min_host_fract: f64,
    pub min_host_genes: usize,
    pub min_viral_genes: usize,
    pub max_iterations: usize,
}

impl Default for RegionParams {
    fn default() -> Self {
        RegionParams {
            gc_weight: 0.02,
            delta_cutoff: 1.285,
            min_host_fract: 0.30,
            min_host_genes: 2,
            min_viral_genes: 2,
            max_iterations: 5,
        }
    }
}

struct Delta {
    delta: f64,
    coords: [usize; 4],
    v1: Vec<i8>,
    v2: Vec<i8>,
    win1_fract_host: f64,
    win2_fract_host: f64,
}

fn mean_i8(v: &[i8]) -> f64 {
    v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64
}

fn mean_f64(v: &[f64]) -> f64 {
    np_mean(v)
}

fn sign(x: f64) -> f64 {
    if x > 0.0 {
        1.0
    } else if x < 0.0 {
        -1.0
    } else {
        0.0
    }
}

fn compute_delta(g: &[GeneView], mut s1: usize, e1: usize, s2: usize, mut e2: usize, gc_weight: f64) -> Delta {
    // extend windows to ensure at least 1 annotated gene
    if g[s1..e1].iter().all(|x| x.cat == 0) {
        for j in (0..s1).rev() {
            s1 = j;
            if g[s1].cat != 0 {
                break;
            }
        }
    }
    if g[s2..e2].iter().all(|x| x.cat == 0) {
        // CheckV starts the extension at e2 + 1 (and slices [s2, e2) after)
        for j in (e2 + 1)..g.len() {
            e2 = j;
            if g[e2].cat != 0 {
                break;
            }
        }
    }
    let win1 = &g[s1..e1];
    let win2 = &g[s2..e2.min(g.len())];
    let v1: Vec<i8> = win1.iter().filter(|x| x.cat != 0).map(|x| x.cat).collect();
    let v2: Vec<i8> = win2.iter().filter(|x| x.cat != 0).map(|x| x.cat).collect();
    let g1: Vec<f64> = win1.iter().map(|x| x.gc).collect();
    let g2: Vec<f64> = win2.iter().map(|x| x.gc).collect();
    let delta = if !v1.is_empty() && !v2.is_empty() {
        let dv = mean_i8(&v1) - mean_i8(&v2);
        let dg = (mean_f64(&g1) - mean_f64(&g2)).abs();
        (dv.abs() + dg * gc_weight) * sign(dv)
    } else {
        0.0
    };
    let host = |v: &[i8]| v.iter().filter(|&&c| c == -1).count() as f64;
    Delta {
        delta,
        coords: [s1, e1, s2, e2],
        win1_fract_host: if win1.is_empty() { 0.0 } else { host(&v1) / win1.len() as f64 },
        win2_fract_host: if win2.is_empty() { 0.0 } else { host(&v2) / win2.len() as f64 },
        v1,
        v2,
    }
}

fn count(v: &[i8], c: i8) -> usize {
    v.iter().filter(|&&x| x == c).count()
}

/// One pass of `define_regions` over a genome of `length` bp.
pub fn define_regions(g: &[GeneView], length: usize, p: &RegionParams) -> Vec<Region> {
    let n = g.len();
    if n == 0 {
        return Vec::new();
    }
    let win = (round_ties_even(0.30 * n as f64) as usize).clamp(15, 50);
    let mut breaks: Vec<Delta> = Vec::new();
    loop {
        let s1 = breaks.last().map(|b| b.coords[2]).unwrap_or(0);
        let mut deltas = Vec::new();
        for i in 1..n {
            let e1 = i;
            let (s2, e2) = (i, (i + win).min(n));
            let size1 = e1 as isize - s1 as isize;
            let size2 = e2 - s2;
            let edge1 = s1 == 0;
            let edge2 = e2 == n;
            if (size1 < win as isize && !edge1) || (size2 < win && !edge2) {
                continue;
            }
            // (e1 > s1 here: a left window not starting at the contig edge
            // must hold `win` genes)
            deltas.push(compute_delta(g, s1, e1, s2, e2, p.gc_weight));
        }
        let big = n > 10;
        let filtered: Vec<Delta> = deltas
            .into_iter()
            .filter(|d| {
                if d.delta.abs() < p.delta_cutoff {
                    return false;
                }
                if d.delta < 0.0 && (count(&d.v1, -1) == 0 || count(&d.v2, 1) == 0) {
                    return false;
                }
                if d.delta > 0.0 && (count(&d.v1, 1) == 0 || count(&d.v2, -1) == 0) {
                    return false;
                }
                if d.delta < 0.0 && count(&d.v1, -1) < p.min_host_genes && big {
                    return false;
                }
                if d.delta > 0.0 && count(&d.v2, -1) < p.min_host_genes && big {
                    return false;
                }
                if d.delta > 0.0 && count(&d.v1, 1) < p.min_viral_genes && big {
                    return false;
                }
                if d.delta < 0.0 && count(&d.v2, 1) < p.min_viral_genes && big {
                    return false;
                }
                if d.win1_fract_host < p.min_host_fract && d.win2_fract_host < p.min_host_fract {
                    return false;
                }
                true
            })
            .collect();
        let mut selected: Option<Delta> = None;
        for d in filtered {
            match &selected {
                None => selected = Some(d),
                Some(s) if sign(d.delta) == sign(s.delta) => {
                    if d.delta.abs() > s.delta.abs() {
                        selected = Some(d);
                    }
                }
                Some(_) => break,
            }
        }
        match selected {
            None => break,
            Some(d) => breaks.push(d),
        }
        if breaks.len() > n {
            break; // defensive: CheckV cannot loop forever either (s1 advances)
        }
    }
    if let Some(last) = breaks.last_mut() {
        last.coords[3] = n;
    }
    let mut regions: Vec<Region> = Vec::new();
    for b in &breaks {
        let [s1, e1, s2, e2] = b.coords;
        let d = compute_delta(g, s1, e1, s2, e2, p.gc_weight);
        let start_pos = regions.last().map(|r| r.end_pos + 1).unwrap_or(1);
        let end_pos = g[e1 - 1].end;
        regions.push(Region {
            viral: d.delta >= 0.0,
            start_pos,
            end_pos,
            start_gene: s1,
            end_gene: e1,
            size: e1 - s1,
            length: end_pos + 1 - start_pos,
            host_genes: count(&d.v1, -1),
            viral_genes: count(&d.v1, 1),
        });
    }
    if let Some(last) = regions.last().cloned() {
        let (s1, e1) = (last.start_gene, last.end_gene);
        let (s2, e2) = (e1, n);
        let d = compute_delta(g, s1, e1, s2, e2, p.gc_weight);
        let start_pos = last.end_pos + 1;
        regions.push(Region {
            viral: !last.viral,
            start_pos,
            end_pos: length,
            start_gene: s2,
            end_gene: e2,
            size: e2 - s2,
            length: (length + 1).saturating_sub(start_pos),
            host_genes: count(&d.v2, -1),
            viral_genes: count(&d.v2, 1),
        });
    }
    regions
}

/// Recursive host-region detection with merging (`define_recursive_regions`
/// + `finalize_regions`). Returns no regions when no host region is found.
pub fn recursive_regions(genes: &[GeneView], length: usize, p: &RegionParams) -> Vec<Region> {
    #[derive(Clone, Copy)]
    struct Span {
        viral: bool,
        start: usize,
        end: usize,
    }
    let mut finished: Vec<Span> = Vec::new();
    let mut active = vec![Span { viral: true, start: 1, end: length }];
    for _ in 0..p.max_iterations {
        if active.is_empty() {
            break;
        }
        let mut next = Vec::new();
        for isl in &active {
            let offset = isl.start - 1;
            let view: Vec<GeneView> = genes
                .iter()
                .filter(|g| g.start >= isl.start && g.end <= isl.end)
                .map(|g| GeneView { start: g.start - offset, end: g.end - offset, gc: g.gc, cat: g.cat })
                .collect();
            let regions = define_regions(&view, isl.end - isl.start + 1, p);
            if regions.is_empty() {
                finished.push(*isl);
                continue;
            }
            for r in regions {
                let span = Span { viral: r.viral, start: isl.start + r.start_pos - 1, end: isl.start + r.end_pos - 1 };
                if span.viral {
                    next.push(span);
                } else {
                    finished.push(span);
                }
            }
        }
        active = next;
    }
    finished.extend(active);
    // merge_regions: sort by (start, end, type) with "host" < "viral"
    finished.sort_by_key(|a| (a.start, a.end, a.viral));
    let mut merged: Vec<Span> = Vec::new();
    for r in finished {
        match merged.last_mut() {
            Some(m) if m.viral == r.viral && r.start <= m.end + 1 => {
                m.end = m.end.max(r.end);
            }
            _ => merged.push(r),
        }
    }
    if !merged.iter().any(|r| !r.viral) {
        return Vec::new();
    }
    merged
        .into_iter()
        .map(|r| {
            let idx: Vec<usize> =
                genes.iter().enumerate().filter(|(_, g)| g.start >= r.start && g.end <= r.end).map(|(i, _)| i).collect();
            Region {
                viral: r.viral,
                start_pos: r.start,
                end_pos: r.end,
                start_gene: idx.first().copied().unwrap_or(0),
                end_gene: idx.last().map(|&i| i + 1).unwrap_or(0),
                size: idx.len(),
                length: r.end + 1 - r.start,
                host_genes: idx.iter().filter(|&&i| genes[i].cat == -1).count(),
                viral_genes: idx.iter().filter(|&&i| genes[i].cat == 1).count(),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn genes(cats: &[i8], gc: f64) -> Vec<GeneView> {
        cats.iter().enumerate().map(|(i, &c)| GeneView { start: i * 1000 + 1, end: i * 1000 + 900, gc, cat: c }).collect()
    }

    #[test]
    fn host_then_virus() {
        // 20 host genes then 25 viral genes
        let mut cats = vec![-1i8; 20];
        cats.extend(vec![1i8; 25]);
        let g = genes(&cats, 50.0);
        let r = recursive_regions(&g, 45_000, &RegionParams::default());
        assert_eq!(r.len(), 2, "{r:?}");
        assert!(!r[0].viral && r[1].viral);
        assert_eq!(r[0].host_genes, 20);
        assert_eq!(r[1].viral_genes, 25);
        assert_eq!(r[1].end_pos, 45_000);
    }

    #[test]
    fn purely_viral_has_no_regions() {
        let g = genes(&[1i8; 30], 40.0);
        assert!(recursive_regions(&g, 30_000, &RegionParams::default()).is_empty());
    }

    #[test]
    fn integrated_provirus_host_viral_host() {
        let mut cats = vec![-1i8; 18];
        cats.extend(vec![1i8; 24]);
        cats.extend(vec![-1i8; 18]);
        let g = genes(&cats, 50.0);
        let r = recursive_regions(&g, 60_000, &RegionParams::default());
        let types: Vec<&str> = r.iter().map(|x| x.type_label()).collect();
        assert_eq!(types, vec!["host", "viral", "host"]);
    }

    #[test]
    fn window_size_rounds_half_to_even() {
        // n = 55 -> 0.3*55 = 16.5 -> 16 (Python round)
        assert_eq!(round_ties_even(0.30 * 55.0) as usize, 16);
        assert_eq!(round_ties_even(0.30 * 75.0) as usize, 22);
    }
}
