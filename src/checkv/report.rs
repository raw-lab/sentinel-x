//! CheckV's output files, byte-compatible with `checkv end_to_end` 1.1.1:
//!
//! ```text
//! quality_summary.tsv   completeness.tsv   contamination.tsv
//! complete_genomes.tsv  viruses.fna        proviruses.fna
//! tmp/gene_features.tsv tmp/aai.tsv        tmp/tr.fna   tmp/proteins.faa
//! ```
//!
//! Numbers are printed the way CheckV's Python prints them (`str(float)`,
//! NumPy rounding where CheckV rounds NumPy scalars, HMMER/DIAMOND table
//! precision where CheckV re-reads those tables).

use std::fmt::Write as _;
use std::path::Path;

use anyhow::{Context, Result};

use super::ContigResult;
use crate::common::pyfmt::{np_round, py_float, py_round};

fn na<T: ToString>(v: Option<T>) -> String {
    v.map(|x| x.to_string()).unwrap_or_else(|| "NA".into())
}

fn na_f(v: Option<f64>) -> String {
    v.map(py_float).unwrap_or_else(|| "NA".into())
}

/// `kmer_freq` column: `str(np.round(kmer_freq, 2))`.
pub fn kmer_freq_str(k: f64) -> String {
    py_float(np_round(k, 2))
}

/// `contamination.tsv`.
pub fn contamination_tsv(results: &[ContigResult]) -> String {
    let mut out = String::from(
        "contig_id\tcontig_length\ttotal_genes\tviral_genes\thost_genes\tprovirus\tproviral_length\thost_length\t\
         region_types\tregion_lengths\tregion_coords_bp\tregion_coords_genes\tregion_viral_genes\tregion_host_genes\n",
    );
    for r in results {
        let _ = write!(out, "{}\t{}\t{}\t{}\t{}", r.id, r.length, r.genes.len(), na(r.viral_genes), na(r.host_genes));
        if r.regions.is_empty() {
            out.push_str("\tNo");
            for _ in 0..8 {
                out.push_str("\tNA");
            }
        } else {
            let join = |f: &dyn Fn(&super::Region) -> String| -> String { r.regions.iter().map(f).collect::<Vec<_>>().join(",") };
            let _ = write!(
                out,
                "\tYes\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                r.proviral_length().unwrap_or(0),
                r.host_length().unwrap_or(0),
                r.region_types(),
                join(&|x| x.length.to_string()),
                join(&|x| format!("{}-{}", x.start_pos, x.end_pos)),
                join(&|x| format!("{}-{}", x.start_gene + 1, x.end_gene)),
                join(&|x| x.viral_genes.to_string()),
                join(&|x| x.host_genes.to_string()),
            );
        }
        out.push('\n');
    }
    out
}

/// `tmp/gene_features.tsv`.
pub fn gene_features_tsv(results: &[ContigResult]) -> String {
    let mut out = String::from("contig_id\tgene_num\tstart\tend\tstrand\tgc\thmm_cat\thmm_db\thmm_name\tevalue\tscore\n");
    for r in results {
        for (g, a) in r.genes.iter().zip(&r.annotations) {
            let num = g.id.rsplit('_').next().unwrap_or(&g.id);
            let _ = write!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                r.id,
                num,
                g.start,
                g.end,
                g.strand,
                py_float(py_round(a.gc, 1)),
                a.cat
            );
            match &a.hit {
                Some(h) => {
                    let _ = write!(
                        out,
                        "\t{}\t{}\t{}\t{}",
                        h.database,
                        h.hmm,
                        py_float(h.printed_evalue()),
                        py_float(h.printed_score())
                    );
                }
                None => out.push_str("\tNA\tNA\tNA\tNA"),
            }
            out.push('\n');
        }
    }
    out
}

/// `tmp/aai.tsv` (every target genome, before the 50%-of-top filter).
pub fn aai_tsv(results: &[ContigResult]) -> String {
    let mut out = String::from("query\ttarget\taligned_genes\tpercent_genes\taligned_length\tpercent_length\tidentity\tscore\n");
    let mut rows: Vec<(&str, &super::AaiRecord)> =
        results.iter().flat_map(|r| r.aai_table.iter().map(move |a| (r.id.as_str(), a))).collect();
    rows.sort_by_key(|x| x.1.order);
    for (id, a) in rows {
        {
            let _ = writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                id,
                a.target,
                a.aligned_genes,
                py_float(a.percent_genes),
                a.aligned_length,
                py_float(a.percent_length),
                py_float(a.identity),
                py_float(a.score)
            );
        }
    }
    out
}

/// `completeness.tsv`.
pub fn completeness_tsv(results: &[ContigResult]) -> String {
    let mut out = String::from(
        "contig_id\tcontig_length\tviral_length\taai_expected_length\taai_completeness\taai_confidence\taai_error\t\
         aai_num_hits\taai_top_hit\taai_id\taai_af\thmm_completeness_lower\thmm_completeness_upper\thmm_num_hits\tkmer_freq\n",
    );
    for r in results {
        let c = &r.completeness;
        let top = c.aai_hits.first();
        let _ = writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.id,
            r.length,
            na(c.viral_length),
            na_f(c.aai_expected_length),
            na_f(c.aai_completeness),
            na(c.aai_confidence),
            na_f(c.aai_error.flatten()),
            c.aai_hits.len(),
            na(top.map(|t| t.target.clone())),
            na_f(top.map(|t| t.identity)),
            na_f(top.map(|t| t.percent_length)),
            na_f(c.hmm_completeness.map(|h| h.0)),
            na_f(c.hmm_completeness.map(|h| h.1)),
            c.hmm_hits,
            kmer_freq_str(c.kmer_freq)
        );
    }
    out
}

/// `complete_genomes.tsv`.
pub fn complete_genomes_tsv(results: &[ContigResult]) -> String {
    let mut out = String::from(
        "contig_id\tcontig_length\tkmer_freq\tprediction_type\tconfidence_level\tconfidence_reason\t\
         repeat_length\trepeat_count\trepeat_n_freq\trepeat_mode_base_freq\trepeat_seq\n",
    );
    for r in results {
        let Some(cg) = &r.complete else { continue };
        let _ = write!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}",
            r.id,
            r.length,
            kmer_freq_str(r.completeness.kmer_freq),
            cg.prediction,
            cg.confidence,
            cg.reason
        );
        match &cg.tr {
            Some(t) => {
                let _ = writeln!(
                    out,
                    "\t{}\t{}\t{}\t{}\t{}",
                    t.len(),
                    t.count,
                    py_float(t.n_freq),
                    py_float(t.mode_freq),
                    String::from_utf8_lossy(&t.seq)
                );
            }
            None => out.push_str("\tNA\tNA\tNA\tNA\tNA\n"),
        }
    }
    out
}

/// `quality_summary.tsv`. Without a database (`have_db = false`) the
/// contamination-module columns are `NA`, as CheckV leaves them when
/// `contamination.tsv` is absent.
pub fn quality_summary_tsv(results: &[ContigResult], have_db: bool) -> String {
    let mut out = String::from(
        "contig_id\tcontig_length\tprovirus\tproviral_length\tgene_count\tviral_genes\thost_genes\tcheckv_quality\t\
         miuvig_quality\tcompleteness\tcompleteness_method\tcontamination\tkmer_freq\twarnings\n",
    );
    for r in results {
        let (prov, plen, ngenes) = if have_db {
            (if r.provirus() { "Yes" } else { "No" }.to_string(), na(r.proviral_length()), r.genes.len().to_string())
        } else {
            ("NA".into(), "NA".into(), "NA".into())
        };
        let _ = writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.id,
            r.length,
            prov,
            plen,
            ngenes,
            na(r.viral_genes),
            na(r.host_genes),
            r.quality,
            r.miuvig,
            na_f(r.completeness_value),
            r.completeness_method,
            na_f(r.contamination),
            kmer_freq_str(r.completeness.kmer_freq),
            r.warnings.join("; ")
        );
    }
    out
}

/// `viruses.fna`: contigs without host regions.
pub fn viruses_fna(results: &[ContigResult]) -> String {
    let mut out = String::new();
    for r in results.iter().filter(|r| r.regions.is_empty()) {
        out.push('>');
        out.push_str(&r.id);
        out.push('\n');
        out.push_str(&String::from_utf8_lossy(&r.seq));
        out.push('\n');
    }
    out
}

/// `proviruses.fna`: the viral regions of proviruses (every region when a
/// contig has none labelled viral), `>{id}_{i} {start}-{end}/{length}`.
pub fn proviruses_fna(results: &[ContigResult]) -> String {
    let mut out = String::new();
    for r in results.iter().filter(|r| !r.regions.is_empty()) {
        let viral: Vec<&super::Region> = r.regions.iter().filter(|x| x.viral).collect();
        let regs: Vec<&super::Region> = if viral.is_empty() { r.regions.iter().collect() } else { viral };
        for (i, x) in regs.iter().enumerate() {
            let s = x.start_pos.saturating_sub(1).min(r.seq.len());
            let e = x.end_pos.min(r.seq.len()).max(s);
            let _ = writeln!(out, ">{}_{} {}-{}/{}", r.id, i + 1, x.start_pos, x.end_pos, r.length);
            out.push_str(&String::from_utf8_lossy(&r.seq[s..e]));
            out.push('\n');
        }
    }
    out
}

/// `tmp/tr.fna`: detected terminal repeats.
pub fn tr_fna(results: &[ContigResult]) -> String {
    let mut out = String::new();
    for r in results {
        if let Some(t) = r.complete.as_ref().and_then(|c| c.tr.as_ref()) {
            let _ = writeln!(out, ">{}_{}", t.kind.label(), r.id);
            out.push_str(&String::from_utf8_lossy(&t.seq));
            out.push('\n');
        }
    }
    out
}

/// `tmp/proteins.faa` with Prodigal-style headers (the gene coordinates
/// CheckV parses from fields 3, 5 and 7), 60 residues per line.
pub fn proteins_faa(results: &[ContigResult]) -> String {
    let mut out = String::new();
    for r in results {
        for g in &r.genes {
            let _ = writeln!(
                out,
                ">{} # {} # {} # {} # partial={}{}",
                g.id, g.start, g.end, g.strand, g.partial[0] as u8, g.partial[1] as u8
            );
            for chunk in g.protein.chunks(60) {
                out.push_str(&String::from_utf8_lossy(chunk));
                out.push('\n');
            }
        }
    }
    out
}

/// Write every CheckV output file under `dir` (and `dir/tmp`). With
/// `append`, rows are added to existing files without repeating headers
/// (chunked runs).
pub fn write_all(dir: &Path, results: &[ContigResult], have_db: bool, append: bool) -> Result<()> {
    use std::io::Write as _;
    let tmp = dir.join("tmp");
    std::fs::create_dir_all(&tmp).with_context(|| format!("creating {}", tmp.display()))?;
    let w = |name: &Path, body: String, has_header: bool| -> Result<()> {
        if !append {
            return std::fs::write(name, body).with_context(|| format!("writing {}", name.display()));
        }
        let rest = if has_header { body.split_once('\n').map(|x| x.1.to_string()).unwrap_or_default() } else { body };
        let mut f = std::fs::OpenOptions::new().append(true).create(true).open(name)?;
        f.write_all(rest.as_bytes()).with_context(|| format!("appending to {}", name.display()))
    };
    w(&dir.join("quality_summary.tsv"), quality_summary_tsv(results, have_db), true)?;
    w(&dir.join("completeness.tsv"), completeness_tsv(results), true)?;
    w(&dir.join("complete_genomes.tsv"), complete_genomes_tsv(results), true)?;
    if have_db {
        w(&dir.join("contamination.tsv"), contamination_tsv(results), true)?;
        w(&dir.join("viruses.fna"), viruses_fna(results), false)?;
        w(&dir.join("proviruses.fna"), proviruses_fna(results), false)?;
        w(&tmp.join("gene_features.tsv"), gene_features_tsv(results), true)?;
        w(&tmp.join("aai.tsv"), aai_tsv(results), true)?;
    }
    w(&tmp.join("tr.fna"), tr_fna(results), false)?;
    w(&tmp.join("proteins.faa"), proteins_faa(results), false)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::*;
    use super::*;

    fn result_with_regions() -> ContigResult {
        let seq: Vec<u8> = (0..1000).map(|i| b"ACGT"[i % 4]).collect();
        ContigResult {
            id: "c1".into(),
            length: 1000,
            seq,
            genes: vec![ViralGene {
                id: "c1_1".into(),
                start: 1,
                end: 300,
                strand: 1,
                protein: b"MK*".to_vec(),
                partial: [false, false],
            }],
            annotations: vec![GeneAnnotation {
                gc: 51.66666,
                cat: -1,
                hit: Some(HmmHit { hmm: "PF001".into(), database: "Pfam".into(), evalue: 3.71e-05, score: 25.04 }),
            }],
            viral_genes: Some(0),
            host_genes: Some(1),
            regions: vec![
                Region {
                    viral: false,
                    start_pos: 1,
                    end_pos: 400,
                    start_gene: 0,
                    end_gene: 1,
                    size: 1,
                    length: 400,
                    host_genes: 1,
                    viral_genes: 0,
                },
                Region {
                    viral: true,
                    start_pos: 401,
                    end_pos: 1000,
                    start_gene: 1,
                    end_gene: 1,
                    size: 0,
                    length: 600,
                    host_genes: 0,
                    viral_genes: 0,
                },
            ],
            completeness: Completeness { kmer_freq: 1.0, ..Default::default() },
            complete: None,
            quality: "Not-determined",
            miuvig: "Genome-fragment",
            completeness_value: None,
            completeness_method: "NA".into(),
            contamination: Some(40.0),
            warnings: vec!["no viral genes detected".into()],
            aai_table: Vec::new(),
        }
    }

    #[test]
    fn contamination_rows() {
        let r = result_with_regions();
        let t = contamination_tsv(std::slice::from_ref(&r));
        let row = t.lines().nth(1).unwrap();
        assert_eq!(row, "c1\t1000\t1\t0\t1\tYes\t600\t400\thost,viral\t400,600\t1-400,401-1000\t1-1,2-1\t0,0\t1,0");
        let p = proviruses_fna(std::slice::from_ref(&r));
        assert!(p.starts_with(">c1_1 401-1000/1000\n"));
        assert_eq!(p.lines().nth(1).unwrap().len(), 600);
        assert!(viruses_fna(std::slice::from_ref(&r)).is_empty());
    }

    #[test]
    fn gene_feature_numbers_are_python_formatted() {
        let r = result_with_regions();
        let t = gene_features_tsv(std::slice::from_ref(&r));
        let row = t.lines().nth(1).unwrap();
        // gc round(51.66666, 1); tblout "%9.2g" e-value and "%6.1f" score
        assert_eq!(row, "c1\t1\t1\t300\t1\t51.7\t-1\tPfam\tPF001\t3.7e-05\t25.0");
    }

    #[test]
    fn quality_summary_row() {
        let r = result_with_regions();
        let t = quality_summary_tsv(std::slice::from_ref(&r), true);
        let row = t.lines().nth(1).unwrap();
        assert_eq!(
            row,
            "c1\t1000\tYes\t600\t1\t0\t1\tNot-determined\tGenome-fragment\tNA\tNA\t40.0\t1.0\tno viral genes detected"
        );
        let f = proteins_faa(std::slice::from_ref(&r));
        assert_eq!(f, ">c1_1 # 1 # 300 # 1 # partial=00\nMK*\n");
    }
}
