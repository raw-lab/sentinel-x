//! The `sentinel run` pipeline: read contigs, call genes once, route each
//! contig/genome to the right validators, and emit a unified report.

use std::collections::{HashMap, HashSet};
use std::io::Write;

use anyhow::{Context, Result};
use rayon::prelude::*;

use sentinel_checkm2::{FeatureSpace, Models};
use sentinel_checkv::{assess_contig, GeneCat};
use sentinel_core::composition::gc_content;
use sentinel_core::{call_genes, read_fasta, Contig, Domain, Gene, GeneSet, QcRecord, QcReport, Quality};
use sentinel_gpu::Backend;
use sentinel_novel::{chimera_scan, run_sage, run_trace};

use crate::Options;

/// Entry point for `sentinel run`.
pub fn run(mut opts: Options) -> Result<()> {
    // Banner + diagnostics go to stderr so stdout stays clean for the TSV.
    eprintln!("{}", sentinel_core::banner());

    // If a database root was given, report readiness and auto-discover the
    // CheckM2 models when they live under <root>/checkm2/.
    if let Some(root) = &opts.db_root {
        match sentinel_db::DbRoot::resolve(Some(root)) {
            Ok(dbroot) => {
                eprintln!("[sentinel] database root: {}", dbroot.root.display());
                for (db, st) in dbroot.status_all() {
                    eprintln!("[sentinel]   {:<8} {}", db.id(), st.symbol());
                }
                if opts.checkm2_comp.is_none() {
                    if let Some(p) = find_under(&dbroot.dir(sentinel_db::Database::CheckM2), "general_model_COMP.gbm") {
                        eprintln!("[sentinel] using CheckM2 completeness model {}", p);
                        opts.checkm2_comp = Some(p);
                    }
                }
                if opts.checkm2_cont.is_none() {
                    if let Some(p) = find_under(&dbroot.dir(sentinel_db::Database::CheckM2), "model_CONT.gbm") {
                        eprintln!("[sentinel] using CheckM2 contamination model {}", p);
                        opts.checkm2_cont = Some(p);
                    }
                }
            }
            Err(e) => eprintln!("[sentinel] (could not resolve --db: {e})"),
        }
    }

    let inputs = resolve_inputs(&opts)?;

    // Read all contigs.
    let mut contigs: Vec<Contig> = Vec::new();
    for p in &inputs {
        let cs = read_fasta(p).with_context(|| format!("reading FASTA '{p}'"))?;
        contigs.extend(cs);
    }
    let n_before = contigs.len();
    contigs.retain(|c| c.len() >= opts.min_len);
    eprintln!(
        "[sentinel] read {n_before} contigs; {} pass min-len {}",
        contigs.len(),
        opts.min_len
    );
    if contigs.is_empty() {
        anyhow::bail!("no contigs passed the length filter");
    }

    let backend = if opts.gpu { Backend::best() } else { Backend::Cpu };
    eprintln!("[sentinel] compute backend: {}", backend.describe());

    // ---- One gene call for everything ----
    eprintln!(
        "[sentinel] calling genes (dynamite, domain = {})...",
        opts.domain.label()
    );
    let geneset = call_genes(&contigs, opts.domain, opts.threads)?;
    let by_contig = geneset.by_contig();
    eprintln!(
        "[sentinel] {} genes called; translation table {}",
        geneset.genes.len(),
        geneset.trans_table
    );

    let mut report = QcReport::new();

    // ---- Per-contig analyses ----
    let viral_like = matches!(opts.domain, Domain::Virus | Domain::Mixed);
    eprintln!(
        "[sentinel] per-contig analysis ({}{})...",
        if viral_like { "CheckV-structural + " } else { "" },
        if opts.novel { "novel TRACE/CHIMERA" } else { "composition only" }
    );

    let per_contig: Vec<Vec<QcRecord>> = if opts.scale {
        contigs
            .par_iter()
            .map(|c| analyse_contig(c, &geneset, &by_contig, &opts, viral_like))
            .collect()
    } else {
        contigs
            .iter()
            .map(|c| analyse_contig(c, &geneset, &by_contig, &opts, viral_like))
            .collect()
    };
    for recs in per_contig {
        report.records.extend(recs);
    }

    // ---- Genome-level analyses (single-genome domains) ----
    if matches!(
        opts.domain,
        Domain::Bacteria | Domain::Archaea | Domain::Eukaryote
    ) {
        let total_len: usize = contigs.iter().map(|c| c.len()).sum();
        let genome_id = genome_label(&opts.input);

        // SAGE self-consistency (always, database-free).
        if opts.novel {
            let s = run_sage(&geneset.genes, total_len, None);
            let mut rec = QcRecord::new(&genome_id, total_len, opts.domain, "sage");
            rec.gene_count = s.n_genes;
            rec.contamination = Some(s.est_contamination);
            rec.completeness = s.est_completeness;
            rec.method = "codon-coherence".into();
            rec.quality = Quality::NotDetermined;
            rec = rec
                .with_extra("codon_coherence", format!("{:.4}", s.codon_coherence))
                .with_extra("gene_density_per_kb", format!("{:.3}", s.gene_density))
                .with_extra("compositional_outliers", s.n_outliers.to_string())
                .with_extra("redundant_pairs", s.redundant_pairs.to_string());
            if s.est_contamination > 10.0 {
                rec.warnings.push(format!(
                    "SAGE contamination proxy {:.1}% (compositional outliers / redundancy)",
                    s.est_contamination
                ));
            }
            report.push(rec);
        }

        // CheckM2 (bacteria/archaea).
        if matches!(opts.domain, Domain::Bacteria | Domain::Archaea) {
            match checkm2_genome(&genome_id, total_len, &geneset, &opts) {
                Ok(rec) => report.push(rec),
                Err(e) => eprintln!("[sentinel] CheckM2 skipped: {e:#}"),
            }
        }

        if matches!(opts.domain, Domain::Eukaryote) {
            eprintln!(
                "[sentinel] BUSCO requires a lineage dataset (not shipped); reporting SAGE/composition only for eukaryote input"
            );
        }
    }

    // ---- Output ----
    let tsv = report.to_tsv();
    match &opts.output {
        Some(path) => {
            std::fs::write(path, &tsv).with_context(|| format!("writing '{path}'"))?;
            eprintln!("[sentinel] wrote {} records to {path}", report.records.len());
        }
        None => {
            let mut out = std::io::stdout().lock();
            out.write_all(tsv.as_bytes())?;
        }
    }
    if let Some(jpath) = &opts.json {
        let json = report_to_json(&report);
        std::fs::write(jpath, serde_json::to_string_pretty(&json)?)
            .with_context(|| format!("writing JSON '{jpath}'"))?;
        eprintln!("[sentinel] wrote JSON report to {jpath}");
    }

    // Summary line.
    let complete = report
        .records
        .iter()
        .filter(|r| r.quality == Quality::Complete)
        .count();
    eprintln!(
        "[sentinel] done — {} records, {} flagged Complete",
        report.records.len(),
        complete
    );
    Ok(())
}

/// Analyse one contig: viral structural (CheckV) + novel TRACE + CHIMERA.
fn analyse_contig(
    c: &Contig,
    geneset: &GeneSet,
    by_contig: &HashMap<&str, Vec<usize>>,
    opts: &Options,
    viral_like: bool,
) -> Vec<QcRecord> {
    let idxs: Vec<usize> = by_contig.get(c.id.as_str()).cloned().unwrap_or_default();
    let genes: Vec<Gene> = idxs.iter().map(|&i| geneset.genes[i].clone()).collect();
    let mut recs = Vec::new();

    if viral_like {
        // CheckV database-free assessment (terminal repeats + structure).
        let cats = vec![GeneCat::Unknown; idxs.len()];
        let rec = assess_contig(c, &idxs, geneset, &cats, &[], &HashMap::new());
        recs.push(rec);

        if opts.novel {
            // TRACE: terminal-repeat & architecture completeness.
            let t = run_trace(&c.seq, &genes);
            let mut r = QcRecord::new(c.id.clone(), c.len(), opts.domain, "trace");
            r.gene_count = genes.len();
            r.method = "architecture".into();
            r.quality = match t.confidence {
                sentinel_novel::TraceConfidence::Complete => Quality::Complete,
                sentinel_novel::TraceConfidence::HighlyComplete => Quality::High,
                sentinel_novel::TraceConfidence::Partial => Quality::Medium,
                sentinel_novel::TraceConfidence::Fragment => Quality::Low,
            };
            r = r
                .with_extra("trace_confidence", t.confidence.label().to_string())
                .with_extra("coding_density", format!("{:.3}", t.coding_density))
                .with_extra("max_gap_frac", format!("{:.3}", t.max_gap_frac));
            if t.terminal_repeat.kind != sentinel_checkv::TrKind::None {
                r = r.with_extra(
                    "terminal_repeat",
                    format!("{:?}:{}bp", t.terminal_repeat.kind, t.terminal_repeat.length),
                );
            }
            if t.circular_overlap > 0 {
                r = r.with_extra("circular_overlap_bp", t.circular_overlap.to_string());
            }
            recs.push(r);
        }
    }

    if opts.novel {
        // CHIMERA-SCAN: compositional coherence (contamination / mis-assembly).
        let ch = chimera_scan(&c.seq, 5000, 3.0);
        let mut r = QcRecord::new(c.id.clone(), c.len(), opts.domain, "chimera-scan");
        r.gene_count = genes.len();
        r.method = "tetranucleotide-coherence".into();
        r.contamination = Some(ch.outlier_fraction * 100.0);
        r.quality = Quality::NotDetermined;
        r = r
            .with_extra("n_windows", ch.n_windows.to_string())
            .with_extra("max_zscore", format!("{:.2}", ch.max_zscore))
            .with_extra("outlier_fraction", format!("{:.3}", ch.outlier_fraction))
            .with_extra("gc", format!("{:.3}", gc_content(&c.seq)));
        if ch.flagged {
            let bp = ch
                .breakpoint
                .map(|b| b.to_string())
                .unwrap_or_else(|| "n/a".into());
            r.warnings
                .push(format!("possible chimera/contamination near base {bp}"));
            if let Some(b) = ch.breakpoint {
                r = r.with_extra("breakpoint", b.to_string());
            }
        }
        recs.push(r);
    }

    recs
}

/// Build the CheckM2 genome record (LightGBM when models are supplied).
fn checkm2_genome(
    genome_id: &str,
    total_len: usize,
    geneset: &GeneSet,
    opts: &Options,
) -> Result<QcRecord> {
    let fs = FeatureSpace::load().context("loading CheckM2 feature space")?;

    // KO presence set from --ko (one KO id per line), else empty.
    let mut ko_present: HashSet<String> = HashSet::new();
    if let Some(kof) = &opts.ko_file {
        let txt = std::fs::read_to_string(kof).with_context(|| format!("reading KO file '{kof}'"))?;
        for line in txt.lines() {
            let k = line.trim();
            if !k.is_empty() {
                ko_present.insert(k.to_string());
            }
        }
    }
    let ko_annotated = !ko_present.is_empty();
    let features = fs.build_vector(&geneset.genes, &ko_present);

    match (&opts.checkm2_comp, &opts.checkm2_cont) {
        (Some(comp), Some(cont)) => {
            let models = Models::load(comp, cont).context("loading CheckM2 models")?;
            Ok(sentinel_checkm2::assess(
                genome_id,
                total_len,
                opts.domain,
                geneset.genes.len(),
                &features,
                &models,
                ko_annotated,
            ))
        }
        _ => {
            // No models: emit a metadata-only record.
            let mut rec = QcRecord::new(genome_id, total_len, opts.domain, "checkm2");
            rec.gene_count = geneset.genes.len();
            rec.method = "metadata-only".into();
            rec.quality = Quality::NotDetermined;
            rec.warnings.push(
                "CheckM2 LightGBM models not provided (--checkm2-comp/--checkm2-cont); \
                 feature vector built but no prediction made"
                    .into(),
            );
            rec = rec.with_extra("feature_vector_len", fs.len().to_string());
            if !ko_annotated {
                rec.warnings
                    .push("no KO annotations (--ko) — KO feature blocks are zero".into());
            }
            Ok(rec)
        }
    }
}

/// Resolve the input path(s), optionally chunking a large file via cleaver.
fn resolve_inputs(opts: &Options) -> Result<Vec<String>> {
    #[cfg(feature = "reads")]
    {
        if let Some(chunk_size) = opts.chunk_size {
            use std::path::Path;
            let dest = std::env::temp_dir().join(format!("sentinel_chunks_{}", std::process::id()));
            let cfg = cleaver_core::Config {
                chunk_size,
                mode: cleaver_core::Mode::Delimiter {
                    delim: b">".to_vec(),
                    matching: cleaver_core::Match::StartsWith,
                },
            };
            let stats = cleaver_core::chunk_file(Path::new(&opts.input), &dest, &cfg)
                .context("chunking input with cleaver")?;
            eprintln!(
                "[sentinel] cleaver split input into {} chunk(s) in {}",
                stats.chunks.len(),
                dest.display()
            );
            return Ok(stats
                .chunks
                .iter()
                .map(|p| p.to_string_lossy().into_owned())
                .collect());
        }
    }
    #[cfg(not(feature = "reads"))]
    {
        if opts.chunk_size.is_some() {
            anyhow::bail!("--chunk-size requires the `reads` feature");
        }
    }
    Ok(vec![opts.input.clone()])
}

/// Recursively search for a file by name beneath `base`, returning its path.
fn find_under(base: &std::path::Path, name: &str) -> Option<String> {
    if !base.exists() {
        return None;
    }
    let mut stack = vec![base.to_path_buf()];
    while let Some(d) = stack.pop() {
        let rd = std::fs::read_dir(&d).ok()?;
        for entry in rd.flatten() {
            let p = entry.path();
            if p.is_dir() {
                stack.push(p);
            } else if p.file_name().and_then(|s| s.to_str()) == Some(name) {
                return Some(p.to_string_lossy().into_owned());
            }
        }
    }
    None
}

fn genome_label(input: &str) -> String {
    std::path::Path::new(input)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("genome")
        .to_string()
}

/// Serialise a report to a JSON value (no serde derive needed on the records).
fn report_to_json(report: &QcReport) -> serde_json::Value {
    let recs: Vec<serde_json::Value> = report
        .records
        .iter()
        .map(|r| {
            serde_json::json!({
                "contig_id": r.contig_id,
                "length": r.length,
                "domain": r.domain.label(),
                "validator": r.validator,
                "gene_count": r.gene_count,
                "completeness": r.completeness,
                "contamination": r.contamination,
                "quality": r.quality.label(),
                "method": r.method,
                "warnings": r.warnings,
                "extra": r.extra.iter().cloned().collect::<HashMap<String, String>>(),
            })
        })
        .collect();
    serde_json::json!({ "sentinel_version": sentinel_core::VERSION, "records": recs })
}
