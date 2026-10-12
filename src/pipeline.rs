//! The `sentinel-x run` pipeline: read contigs, call genes, route the input
//! to the validators for its domain, and write the unified report plus each
//! tool's own output tables.
//!
//! * **virus / mixed** — per contig: CheckV, TRACE, CHIMERA-SCAN, and
//!   viralVerify with `--viralverify`;
//! * **bacteria / archaea** — per genome (each input FASTA file): CheckM2,
//!   SAGE, BUSCO with `--lineage`; CHIMERA-SCAN per contig;
//! * **eukaryote** — per genome: EukCC with `--markers`, BUSCO with
//!   `--lineage`, SAGE; CHIMERA-SCAN per contig.

use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use rayon::prelude::*;

use crate::common::{call_genes, read_fasta, Contig, Domain, GeneSet, QcRecord, QcReport, Quality};
use crate::db::DbRoot;
use crate::novel::{chimera_scan, run_sage, run_trace};
use crate::Options;

/// Use `threads` worker threads for everything that runs on rayon.
fn init_threads(threads: usize) {
    let n = if threads == 0 { std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1) } else { threads };
    let _ = rayon::ThreadPoolBuilder::new().num_threads(n).build_global();
}

/// An input genome (a FASTA file) or, for per-contig domains, all contigs.
struct Genome {
    name: String,
    path: PathBuf,
    contigs: Vec<Contig>,
}

fn is_fasta_name(p: &Path) -> bool {
    let n = p.file_name().and_then(|s| s.to_str()).unwrap_or("").to_ascii_lowercase();
    let n = n.strip_suffix(".gz").unwrap_or(&n);
    [".fa", ".fna", ".fasta", ".fas", ".ffn", ".contigs", ".faa"].iter().any(|e| n.ends_with(e))
}

/// CheckM2's bin name: the file name without its (last) extension, after
/// removing `.gz`.
fn genome_name(p: &Path) -> String {
    let n = p.file_name().and_then(|s| s.to_str()).unwrap_or("genome");
    let n = n.strip_suffix(".gz").unwrap_or(n);
    match n.rsplit_once('.') {
        Some((stem, _)) if !stem.is_empty() => stem.to_string(),
        _ => n.to_string(),
    }
}

fn read_inputs(opts: &Options) -> Result<Vec<Genome>> {
    let input = Path::new(&opts.input);
    let mut files: Vec<PathBuf> = Vec::new();
    if input.is_dir() {
        for e in std::fs::read_dir(input).with_context(|| format!("listing {}", input.display()))?.flatten() {
            let p = e.path();
            if p.is_file() && is_fasta_name(&p) && std::fs::metadata(&p).map(|m| m.len() > 0).unwrap_or(false) {
                files.push(p);
            }
        }
        files.sort();
        if files.is_empty() {
            bail!("no FASTA files (.fa/.fna/.fasta[.gz]) in {}", input.display());
        }
    } else {
        files.push(input.to_path_buf());
    }
    let mut out = Vec::new();
    for f in files {
        let mut contigs = read_fasta(&f).with_context(|| format!("reading FASTA '{}'", f.display()))?;
        let before = contigs.len();
        if opts.min_len > 0 {
            contigs.retain(|c| c.len() >= opts.min_len);
        }
        if contigs.len() < before {
            eprintln!("[sentinel-x] {}: {} of {} contigs pass --min-len {}", f.display(), contigs.len(), before, opts.min_len);
        }
        out.push(Genome { name: genome_name(&f), path: f, contigs });
    }
    Ok(out)
}

/// Entry point for `sentinel-x run`.
pub fn run(opts: Options) -> Result<()> {
    eprintln!("{}", crate::common::banner());
    init_threads(opts.threads);
    let dbroot = DbRoot::resolve(opts.db_root.as_deref());
    if dbroot.root.is_dir() {
        eprintln!("[sentinel-x] database root: {}", dbroot.root.display());
        for (db, st) in dbroot.status_all() {
            eprintln!("[sentinel-x]   {:<8} {}", db.id(), st.symbol());
        }
    } else {
        eprintln!(
            "[sentinel-x] no database root at {} (set --db or $SENTINEL_X_DB); database-backed validators are skipped",
            dbroot.root.display()
        );
    }
    if opts.gpu {
        eprintln!("[sentinel-x] note: --gpu has no effect (no GPU kernels are used); running on CPU");
    }
    let genomes = read_inputs(&opts)?;
    let n_contigs: usize = genomes.iter().map(|g| g.contigs.len()).sum();
    if n_contigs == 0 {
        bail!("no contigs to analyse");
    }
    eprintln!("[sentinel-x] {} input file(s), {} contigs; domain {}", genomes.len(), n_contigs, opts.domain.label());

    // options that do not apply to the chosen domain
    let per_contig = matches!(opts.domain, Domain::Virus | Domain::Mixed);
    if opts.markers.is_some() && opts.domain != Domain::Eukaryote {
        eprintln!("[sentinel-x] warning: --markers (EukCC) is only used with -d eukaryote");
    }
    if opts.lineage.is_some() && per_contig {
        eprintln!("[sentinel-x] warning: --lineage (BUSCO) is only used for genome-level domains (bacteria, archaea, eukaryote)");
    }
    if opts.viralverify && !per_contig {
        eprintln!("[sentinel-x] warning: --viralverify is only used with -d virus or -d mixed");
    }
    if opts.min_len > 0 && !per_contig {
        eprintln!(
            "[sentinel-x] note: --min-len {} removes contigs before CheckM2/BUSCO/EukCC (CheckM2 itself uses every contig)",
            opts.min_len
        );
    }
    if let Some(d) = &opts.outdir {
        std::fs::create_dir_all(d).with_context(|| format!("creating {d}"))?;
    }

    #[cfg(feature = "rustyomestats")]
    {
        if per_contig && opts.chunk_size.is_some() {
            eprintln!("[sentinel-x] rustyomestats cross-check skipped with --chunk-size (needs whole files)");
        } else if !opts.genes {
            for g in &genomes {
                crate::common::xstats::crosscheck(&g.name, &g.contigs);
            }
        }
    }

    let mut report = QcReport::new();
    let mut stats = StatsAcc::default();
    let mut per_genome: Vec<crate::html_report::GenomeStats> = Vec::new();
    if per_contig && opts.chunk_size.is_some() {
        drop(genomes);
        run_chunked(&opts, &dbroot, &mut report, &mut stats)?;
    } else if per_contig {
        let mut contigs: Vec<Contig> = Vec::new();
        let mut seen: HashMap<String, PathBuf> = HashMap::new();
        for g in genomes {
            for c in g.contigs {
                if let Some(prev) = seen.insert(c.id.clone(), g.path.clone()) {
                    bail!("contig id '{}' occurs in both {} and {}", c.id, prev.display(), g.path.display());
                }
                contigs.push(c);
            }
        }
        let gs = run_per_contig(&contigs, &opts, &dbroot, &mut report, true)?;
        stats.add(&contigs, &gs);
    } else {
        let mut cm2_results = Vec::new();
        let mut cm2 = Checkm2Context::new(&opts)?;
        for g in &genomes {
            if g.contigs.is_empty() {
                eprintln!("[sentinel-x] {}: no contigs, skipped", g.name);
                continue;
            }
            let gs = run_genome(g, &opts, &dbroot, &mut report, &mut cm2, &mut cm2_results)?;
            if !opts.genes {
                let mut one = StatsAcc::default();
                one.add(&g.contigs, &gs);
                per_genome.push((g.name.clone(), one.table()));
                stats.merge(one);
            }
        }
        if let (Some(d), false) = (&opts.outdir, cm2_results.is_empty()) {
            let dir = Path::new(d).join("checkm2");
            std::fs::create_dir_all(&dir)?;
            std::fs::write(dir.join("quality_report.tsv"), crate::checkm2::quality_report_tsv(&cm2_results))?;
            eprintln!("[sentinel-x] wrote {}", dir.join("quality_report.tsv").display());
        }
    }

    write_outputs(&opts, &report, &stats, &per_genome)?;
    let complete = report.records.iter().filter(|r| r.quality == Quality::Complete).count();
    eprintln!("[sentinel-x] done — {} records, {} flagged Complete", report.records.len(), complete);
    Ok(())
}

/// Per-contig validators (virus / mixed): CheckV, TRACE, CHIMERA-SCAN and,
/// with `--viralverify`, viralVerify.
fn run_per_contig(contigs: &[Contig], opts: &Options, dbroot: &DbRoot, report: &mut QcReport, first: bool) -> Result<GeneSet> {
    // ---- genes (prodigal-gv preset, or supplied proteins)
    let gs = call_genes(contigs, opts.domain, opts.threads)?;
    eprintln!("[sentinel-x] {} genes called (DYNAMITE prodigal-gv preset, translation table {})", gs.genes.len(), gs.trans_table);

    // ---- CheckV
    let checkv_db = match dbroot.checkv_dir() {
        Some(d) => match crate::checkv::CheckvDb::load(&d) {
            Ok(db) => {
                eprintln!("[sentinel-x] CheckV database: {} ({} HMM files)", d.display(), db.hmm_files.len());
                Some(db)
            }
            Err(e) => {
                eprintln!("[sentinel-x] CheckV database at {} unusable: {e:#}", d.display());
                None
            }
        },
        None => {
            eprintln!("[sentinel-x] no CheckV database: terminal repeats and copy number only (`sentinel-x download checkv`)");
            None
        }
    };
    let inputs = match &opts.proteins {
        Some(faa) => {
            eprintln!("[sentinel-x] CheckV: using proteins from {faa}");
            crate::checkv::inputs_from_faa(contigs, Path::new(faa))?
        }
        None => crate::checkv::inputs_from_genes(contigs, &gs),
    };
    eprintln!("[sentinel-x] CheckV: annotating, scanning host regions and estimating completeness");
    let params = crate::checkv::Params { strict_diamond_order: opts.checkv_strict, ..Default::default() };
    let results = crate::checkv::run(&inputs, checkv_db.as_ref(), &params)?;
    report.extend(crate::checkv::to_records(&results, opts.domain, checkv_db.is_some()));
    if let Some(d) = &opts.outdir {
        let dir = Path::new(d).join("checkv");
        crate::checkv::report::write_all(&dir, &results, checkv_db.is_some(), !first)?;
        if first {
            eprintln!("[sentinel-x] writing CheckV tables to {}", dir.display());
        }
    }

    // ---- novel modules
    if opts.novel {
        let per = gs.per_contig(contigs.len());
        let recs: Vec<Vec<QcRecord>> = contigs
            .par_iter()
            .zip(per.par_iter())
            .map(|(c, idx)| {
                let genes: Vec<crate::common::Gene> = idx.iter().map(|&i| gs.genes[i].clone()).collect();
                let mut v = vec![trace_record(c, &genes, opts.domain)];
                v.push(chimera_record(c, genes.len(), opts.domain));
                v
            })
            .collect();
        for r in recs {
            report.extend(r);
        }
    }

    // ---- viralVerify (opt-in: a Pfam-A search is the heaviest step)
    if opts.viralverify {
        match dbroot.pfam_hmm() {
            Some(pfam) => {
                eprintln!("[sentinel-x] viralVerify: Pfam-A search ({})", pfam.display());
                let table = crate::viralverify::ClassifierTable::embedded();
                let verdicts = crate::viralverify::run(contigs, &pfam, opts.threads, 7.0, &table)?;
                report.extend(crate::viralverify::to_records(&verdicts, opts.domain));
                if let Some(d) = &opts.outdir {
                    let dir = Path::new(d).join("viralverify");
                    std::fs::create_dir_all(&dir)?;
                    let stem = genome_name(Path::new(&opts.input));
                    let csv = crate::viralverify::result_table_csv(&verdicts);
                    append_or_write(&dir.join(format!("{stem}_result_table.csv")), &csv, !first)?;
                }
            }
            None => eprintln!("[sentinel-x] viralVerify skipped: no Pfam-A library (`sentinel-x download pfam`)"),
        }
    }
    Ok(gs)
}

fn trace_record(c: &Contig, genes: &[crate::common::Gene], domain: Domain) -> QcRecord {
    let t = run_trace(&c.seq, genes);
    let mut r = QcRecord::new(&c.id, c.len(), domain, "trace");
    r.gene_count = genes.len();
    r.method = "architecture".into();
    r.quality = Quality::NotDetermined; // a confidence label, not a completeness estimate
    r = r
        .with_extra("trace_confidence", t.confidence.label())
        .with_extra("coding_density", format!("{:.3}", t.coding_density))
        .with_extra("max_gap_frac", format!("{:.3}", t.max_gap_frac));
    if let Some(tr) = &t.terminal_repeat {
        r = r.with_extra("terminal_repeat", format!("{}:{}bp", tr.kind.label(), tr.len()));
    }
    if let Some(f) = t.repeat_flag {
        r = r.with_extra("repeat_flag", f);
    }
    if t.circular_overlap > 0 {
        r = r.with_extra("circular_overlap_bp", t.circular_overlap.to_string());
    }
    r
}

/// The unified record for a SAGE result (`contigs` = the genome's contigs,
/// empty for protein input).
fn sage_record(name: &str, total_len: usize, domain: Domain, s: &crate::novel::sage::SageResult, contigs: &[Contig]) -> QcRecord {
    let mut rec = QcRecord::new(name, total_len, domain, "sage");
    rec.gene_count = s.n_genes;
    rec.method = "contig-composition (reference-free)".into();
    rec.quality = Quality::NotDetermined;
    rec.contamination = s.discordant_pct().map(|p| (p * 100.0).round() / 100.0);
    rec = rec
        .with_extra("contigs_assessed", s.contigs_assessed.to_string())
        .with_extra("assessed_bp", s.assessed_bp.to_string())
        .with_extra("discordant_bp", s.discordant_bp.to_string())
        .with_extra("discordant_contigs", s.discordant.len().to_string());
    if let Some(p) = s.discordant_pct() {
        rec = rec.with_extra("discordant_bp_pct", format!("{p:.2}"));
        if p >= 5.0 {
            let mut ids: Vec<&str> = s.discordant.iter().take(10).map(|&i| contigs[i].id.as_str()).collect();
            if s.discordant.len() > 10 {
                ids.push("…");
            }
            rec.warnings.push(format!(
                "SAGE: {p:.1}% of the assessed sequence departs from the core composition ({} contig{}: {})",
                s.discordant.len(),
                if s.discordant.len() == 1 { "" } else { "s" },
                ids.join(", ")
            ));
        }
    }
    rec.with_extra("aa_composition_coherence", format!("{:.4}", s.aa_coherence))
        .with_extra("gene_density_per_kb", format!("{:.3}", s.gene_density))
        .with_extra("near_identical_protein_pairs", s.redundant_pairs.to_string())
}

fn chimera_record(c: &Contig, n_genes: usize, domain: Domain) -> QcRecord {
    let ch = chimera_scan(&c.seq, 5000, 0.12);
    let mut r = QcRecord::new(&c.id, c.len(), domain, "chimera-scan");
    r.gene_count = n_genes;
    r.method = "tetranucleotide-junction".into();
    r.quality = Quality::NotDetermined;
    r = r
        .with_extra("assessed", if ch.assessed { "yes" } else { "no (< 4 windows of 5 kb)" })
        .with_extra("n_windows", ch.n_windows.to_string())
        .with_extra("gc_percent", format!("{:.2}", 100.0 * crate::common::composition::gc_content(&c.seq)));
    if ch.assessed {
        r = r
            .with_extra("junction_distance", format!("{:.4}", ch.junction_distance))
            .with_extra("junction_r2", format!("{:.3}", ch.junction_r2))
            .with_extra("max_window_zscore", format!("{:.2}", ch.max_zscore))
            .with_extra("outlier_window_fraction", format!("{:.3}", ch.outlier_fraction));
    }
    if ch.flagged {
        let bp = ch.breakpoint.map(|b| b.to_string()).unwrap_or_else(|| "n/a".into());
        r.warnings.push(format!("compositional junction near base {bp} (D = {:.3}): possible chimera", ch.junction_distance));
        if let Some(b) = ch.breakpoint {
            r = r.with_extra("breakpoint", b.to_string());
        }
    }
    r
}

/// Models and annotation sources shared across genomes.
struct Checkm2Context {
    fs: Option<crate::checkm2::FeatureSpace>,
    models: Option<crate::checkm2::Models>,
    nn: Option<crate::checkm2::nn::SpecificModel>,
    loaded: bool,
}

impl Checkm2Context {
    fn new(opts: &Options) -> Result<Self> {
        if opts.checkm2_comp.is_some() != opts.checkm2_cont.is_some() {
            bail!("--checkm2-comp and --checkm2-cont must be given together");
        }
        Ok(Checkm2Context { fs: None, models: None, nn: None, loaded: false })
    }

    fn ensure(&mut self, opts: &Options) -> Result<()> {
        if self.loaded {
            return Ok(());
        }
        self.loaded = true;
        self.fs = Some(crate::checkm2::FeatureSpace::load()?);
        self.models = match (&opts.checkm2_comp, &opts.checkm2_cont) {
            (Some(a), Some(b)) => Some(crate::checkm2::Models::load(a, b)?),
            _ => crate::checkm2::Models::load_embedded_opt().transpose()?,
        };
        #[cfg(feature = "embed-models")]
        {
            self.nn = Some(crate::checkm2::nn::SpecificModel::embedded()?);
        }
        Ok(())
    }
}

/// Genome-level validators for one input genome.
fn run_genome(
    g: &Genome,
    opts: &Options,
    dbroot: &DbRoot,
    report: &mut QcReport,
    cm2: &mut Checkm2Context,
    cm2_results: &mut Vec<crate::checkm2::Checkm2Result>,
) -> Result<GeneSet> {
    let total_len: usize = g.contigs.iter().map(|c| c.len()).sum();
    let prokaryote = matches!(opts.domain, Domain::Bacteria | Domain::Archaea);
    if opts.genes {
        return run_genome_proteins(g, opts, dbroot, report, cm2, cm2_results);
    }
    eprintln!(
        "[sentinel-x] genome {} ({} contig{}, {} bp)",
        g.name,
        g.contigs.len(),
        if g.contigs.len() == 1 { "" } else { "s" },
        total_len
    );

    // ---- genes
    let (gs, calls) = if prokaryote {
        let calls = crate::checkm2::call_genes(&g.contigs, opts.threads, opts.ttable)?;
        eprintln!(
            "[sentinel-x]   {} genes (Prodigal-style {} mode, translation table {}, coding density {:.3})",
            calls.genes.genes.len(),
            calls.procedure,
            calls.ttable,
            calls.coding_density
        );
        (calls.genes.clone(), Some(calls))
    } else {
        let gs = call_genes(&g.contigs, opts.domain, opts.threads)?;
        eprintln!("[sentinel-x]   {} genes", gs.genes.len());
        (gs, None)
    };

    // ---- CheckM2
    if let Some(calls) = &calls {
        let proteins = calls.genes.protein_records_with_stop();
        if let Some(d) = &opts.outdir {
            let dir = Path::new(d).join("checkm2").join("protein_files");
            std::fs::create_dir_all(&dir)?;
            std::fs::write(dir.join(format!("{}.faa", g.name)), prodigal_faa(&calls.genes))?;
        }
        let stats = crate::checkm2::genome_stats(&g.contigs, calls);
        match checkm2_genome(&g.name, proteins, Some(stats), opts, dbroot, cm2) {
            Ok(r) => {
                report.push(crate::checkm2::to_record(&r, total_len, opts.domain));
                cm2_results.push(r);
            }
            Err(e) => eprintln!("[sentinel-x]   CheckM2 failed: {e:#}"),
        }
    }

    // ---- SAGE
    if opts.novel {
        let s = run_sage(&g.contigs, &gs.genes, total_len);
        report.push(sage_record(&g.name, total_len, opts.domain, &s, &g.contigs));
        if let (Some(d), false) = (&opts.outdir, s.discordant.is_empty()) {
            let dir = Path::new(d).join("sage");
            std::fs::create_dir_all(&dir)?;
            let list: String = s.discordant.iter().map(|&i| format!("{}\t{}\n", g.contigs[i].id, g.contigs[i].len())).collect();
            std::fs::write(dir.join(format!("{}.discordant_contigs.tsv", g.name)), format!("contig_id\tlength\n{list}"))?;
        }
    }

    // ---- EukCC
    if opts.domain == Domain::Eukaryote {
        match &opts.markers {
            Some(m) => {
                let prots = gs.protein_records();
                let targets = crate::hmmer::Targets::from_pairs(&prots);
                let s = crate::eukcc::run(&targets, Path::new(m))?;
                eprintln!(
                    "[sentinel-x]   EukCC: completeness {} contamination {} ({}/{} markers, {} multi-copy)",
                    s.completeness, s.contamination, s.markers_found, s.markers_expected, s.markers_multicopy
                );
                let mut rec = QcRecord::new(&g.name, total_len, opts.domain, "eukcc");
                rec.gene_count = s.orf_count;
                rec.completeness = Some(s.completeness);
                rec.contamination = Some(s.contamination);
                rec.method = format!("marker set {}", s.marker_set);
                rec.quality = Quality::from_mimag(rec.completeness, rec.contamination);
                rec = rec
                    .with_extra("markers_expected", s.markers_expected.to_string())
                    .with_extra("markers_found", s.markers_found.to_string())
                    .with_extra("markers_multicopy", s.markers_multicopy.to_string())
                    .with_extra("markers_missing", s.markers_missing.to_string())
                    .with_extra("tblout_rows", s.rows.to_string())
                    .with_extra("mean_marker_copies", format!("{:.2}", s.mean_copies()));
                if s.contamination >= 10.0 {
                    rec.warnings
                        .push(format!("EukCC contamination {}% ({} multi-copy markers)", s.contamination, s.markers_multicopy));
                }
                report.push(rec);
            }
            None => eprintln!("[sentinel-x]   EukCC needs a marker set (--markers <set.hmm>)"),
        }
    }

    // ---- BUSCO
    if let Some(spec) = &opts.lineage {
        let busco_root = dbroot.dir(crate::db::Database::Busco);
        let dir = crate::busco::Lineage::resolve(spec, Some(&busco_root))?;
        let lineage = crate::busco::Lineage::load(&dir)?;
        lineage_domain_check(&lineage, opts.domain);
        let genes: Vec<crate::busco::BuscoGene> = gs
            .genes
            .iter()
            .map(|x| crate::busco::BuscoGene::new(&x.contig, x.begin, x.end, x.strand, x.protein_with_stop()))
            .collect();
        eprintln!("[sentinel-x]   BUSCO {} ({} BUSCOs)", lineage.name, lineage.total());
        let s = crate::busco::run(&lineage, &genes, true)?;
        eprintln!("[sentinel-x]   BUSCO: {}", s.one_line());
        report.push(crate::busco::to_record(&g.name, total_len, opts.domain, gs.genes.len(), &s));
        if let Some(d) = &opts.outdir {
            let dir = Path::new(d).join("busco").join(&g.name);
            std::fs::create_dir_all(&dir)?;
            std::fs::write(dir.join("full_table.tsv"), crate::busco::full_table_tsv(&s))?;
            std::fs::write(
                dir.join("short_summary.txt"),
                format!(
                    "# The lineage dataset is: {} (Creation date: {}, number of BUSCOs: {})\n\t{}\n",
                    lineage.name,
                    lineage.creation_date.as_deref().unwrap_or("NA"),
                    lineage.total(),
                    s.one_line()
                ),
            )?;
        }
    }

    // ---- CHIMERA-SCAN per contig
    if opts.novel {
        let per = gs.per_contig(g.contigs.len());
        let recs: Vec<QcRecord> = g
            .contigs
            .par_iter()
            .zip(per.par_iter())
            .map(|(c, idx)| chimera_record(c, idx.len(), opts.domain).with_extra("genome", &g.name))
            .collect();
        report.extend(recs);
    }
    Ok(gs)
}

/// Warn when a BUSCO lineage is for another domain of life.
fn lineage_domain_check(l: &crate::busco::Lineage, d: Domain) {
    let want = match d {
        Domain::Bacteria | Domain::Archaea => "prokaryota",
        Domain::Eukaryote => "eukaryota",
        _ => return,
    };
    if !l.domain.is_empty() && !l.domain.eq_ignore_ascii_case(want) {
        eprintln!("[sentinel-x]   warning: BUSCO lineage {} is for {}, the input is {}", l.name, l.domain, d.label());
    }
}

/// Proteins in Prodigal's `.faa` layout (`>id # start # end # strand # …`).
fn prodigal_faa(gs: &GeneSet) -> String {
    let mut out = String::new();
    for g in &gs.genes {
        out.push_str(&format!(
            ">{} # {} # {} # {} # partial={}{}\n",
            g.id(),
            g.begin,
            g.end,
            g.strand,
            g.partial_left as u8,
            g.partial_right as u8
        ));
        for chunk in g.protein_with_stop().chunks(60) {
            out.push_str(&String::from_utf8_lossy(chunk));
            out.push('\n');
        }
    }
    out
}

/// Genome-level validators on supplied proteins (`--genes`).
fn run_genome_proteins(
    g: &Genome,
    opts: &Options,
    dbroot: &DbRoot,
    report: &mut QcReport,
    cm2: &mut Checkm2Context,
    cm2_results: &mut Vec<crate::checkm2::Checkm2Result>,
) -> Result<GeneSet> {
    let proteins: Vec<(String, Vec<u8>)> = g.contigs.iter().map(|c| (c.id.clone(), c.seq.clone())).collect();
    let aa_total: usize = proteins.iter().map(|p| p.1.len()).sum();
    eprintln!("[sentinel-x] genome {} ({} proteins, {} residues)", g.name, proteins.len(), aa_total);
    let genes: Vec<crate::common::Gene> = proteins
        .iter()
        .enumerate()
        .map(|(i, (id, p))| {
            let stop = p.last() == Some(&b'*');
            let aa: Vec<u8> = if stop { p[..p.len() - 1].to_vec() } else { p.clone() };
            crate::common::Gene {
                contig: id.clone(),
                contig_index: i,
                index: 0,
                begin: 1,
                end: 3 * p.len(),
                strand: 1,
                aa,
                partial_left: false,
                partial_right: false,
                stop,
            }
        })
        .collect();
    let gs = GeneSet { genes, trans_table: 0 };
    if matches!(opts.domain, Domain::Bacteria | Domain::Archaea) {
        match checkm2_genome(&g.name, proteins.clone(), None, opts, dbroot, cm2) {
            Ok(r) => {
                report.push(crate::checkm2::to_record(&r, 0, opts.domain));
                cm2_results.push(r);
            }
            Err(e) => eprintln!("[sentinel-x]   CheckM2 failed: {e:#}"),
        }
    }
    if opts.novel {
        let s = run_sage(&[], &gs.genes, 0);
        report.push(sage_record(&g.name, 0, opts.domain, &s, &[]));
    }
    if opts.domain == Domain::Eukaryote {
        if let Some(m) = &opts.markers {
            let targets = crate::hmmer::Targets::from_pairs(&proteins);
            let s = crate::eukcc::run(&targets, Path::new(m))?;
            let mut rec = QcRecord::new(&g.name, 0, opts.domain, "eukcc");
            rec.gene_count = s.orf_count;
            rec.completeness = Some(s.completeness);
            rec.contamination = Some(s.contamination);
            rec.method = format!("marker set {}", s.marker_set);
            rec.quality = Quality::from_mimag(rec.completeness, rec.contamination);
            report.push(rec);
        }
    }
    if let Some(spec) = &opts.lineage {
        let busco_root = dbroot.dir(crate::db::Database::Busco);
        let lineage = crate::busco::Lineage::load(&crate::busco::Lineage::resolve(spec, Some(&busco_root))?)?;
        lineage_domain_check(&lineage, opts.domain);
        let genes: Vec<crate::busco::BuscoGene> = proteins
            .iter()
            .map(|(id, p)| crate::busco::BuscoGene {
                target: id.clone(),
                contig: id.clone(),
                start: 0,
                end: 0,
                protein: p.clone(),
            })
            .collect();
        let s = crate::busco::run(&lineage, &genes, false)?;
        eprintln!("[sentinel-x]   BUSCO (proteins): {}", s.one_line());
        report.push(crate::busco::to_record(&g.name, 0, opts.domain, proteins.len(), &s));
        if let Some(d) = &opts.outdir {
            let dir = Path::new(d).join("busco").join(&g.name);
            std::fs::create_dir_all(&dir)?;
            std::fs::write(dir.join("full_table.tsv"), crate::busco::full_table_tsv(&s))?;
        }
    }
    Ok(gs)
}

fn checkm2_genome(
    name: &str,
    proteins: Vec<(String, Vec<u8>)>,
    stats: Option<crate::checkm2::GenomeStats>,
    opts: &Options,
    dbroot: &DbRoot,
    cm2: &mut Checkm2Context,
) -> Result<crate::checkm2::Checkm2Result> {
    cm2.ensure(opts)?;
    let fs = cm2.fs.as_ref().unwrap();
    let models = cm2
        .models
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("no CheckM2 models (build with `embed-models` or pass --checkm2-comp/--checkm2-cont)"))?;
    let (kos, source): (HashMap<String, u32>, String) = if let Some(f) = &opts.ko_file {
        let t = std::fs::read_to_string(f).with_context(|| format!("reading KO file '{f}'"))?;
        (crate::checkm2::parse_ko_file(&t), format!("--ko {f}"))
    } else {
        match dbroot.checkm2_db() {
            Some(db) => {
                eprintln!("[sentinel-x]   CheckM2: annotating {} proteins against {}", proteins.len(), db.display());
                let (k, st) = crate::checkm2::annotate_kos(&proteins, &db)?;
                eprintln!("[sentinel-x]   CheckM2: {} proteins with a KO hit ({} reference letters)", st.hits, st.db_letters);
                (k, format!("DIAMOND-compatible search of {}", db.file_name().and_then(|s| s.to_str()).unwrap_or("db")))
            }
            None => (HashMap::new(), "CheckM2 database not installed (`sentinel-x download checkm2`)".into()),
        }
    };
    let ko = fs.ko_counts(kos.iter().map(|(k, &n)| (k.as_str(), n)));
    let kos_annotated = ko.iter().filter(|&&c| c > 0).count();
    let plain: Vec<Vec<u8>> = proteins.into_iter().map(|p| p.1).collect();
    let features = fs.build_vector(&fs.metadata(&plain), &ko);
    let prediction = crate::checkm2::predict(&features, models, cm2.nn.as_ref());
    // validation aid: dump the feature vector and raw predictions
    if let Ok(dir) = std::env::var("SENTINEL_X_CHECKM2_DEBUG") {
        let dir = Path::new(&dir);
        std::fs::create_dir_all(dir)?;
        let mut t = String::new();
        for v in &features {
            t.push_str(&crate::common::pyfmt::py_float(*v));
            t.push('\n');
        }
        std::fs::write(dir.join(format!("{name}.features.txt")), t)?;
        std::fs::write(
            dir.join(format!("{name}.prediction.txt")),
            format!(
                "general\t{}\nspecific\t{}\ncontamination\t{}\ncosine\t{}\n",
                crate::common::pyfmt::py_float(prediction.general),
                prediction.specific.map(|x| crate::common::pyfmt::py_float(x as f64)).unwrap_or_default(),
                crate::common::pyfmt::py_float(prediction.contamination),
                prediction.cosine.map(crate::common::pyfmt::py_float).unwrap_or_default()
            ),
        )?;
    }
    eprintln!(
        "[sentinel-x]   CheckM2: completeness {:.2} ({}), contamination {:.2}",
        prediction.completeness, prediction.model, prediction.contamination
    );
    Ok(crate::checkm2::Checkm2Result { name: name.to_string(), prediction, stats, kos_annotated, ko_source: source })
}

/// Write `body`, or append it without its header line.
fn append_or_write(path: &Path, body: &str, append: bool) -> Result<()> {
    if !append {
        return std::fs::write(path, body).with_context(|| format!("writing {}", path.display()));
    }
    let rest = body.split_once('\n').map(|x| x.1).unwrap_or("");
    let mut f = std::fs::OpenOptions::new().append(true).create(true).open(path)?;
    f.write_all(rest.as_bytes())?;
    Ok(())
}

/// `--chunk-size`: split a large FASTA with Cleaver-X and run the
/// per-contig validators chunk by chunk (bounded memory); chunks are
/// deleted afterwards.
#[cfg(feature = "reads")]
fn run_chunked(opts: &Options, dbroot: &DbRoot, report: &mut QcReport, stats: &mut StatsAcc) -> Result<()> {
    let size = opts.chunk_size.unwrap_or(0);
    if size == 0 {
        bail!("--chunk-size must be > 0");
    }
    let input = Path::new(&opts.input);
    if input.is_dir() {
        bail!("--chunk-size needs a single FASTA file");
    }
    let dest = std::env::temp_dir().join(format!("sentinel_x_chunks_{}", std::process::id()));
    let cfg = cleaver::Config {
        chunk_size: size,
        mode: cleaver::Mode::Delimiter { delim: b">".to_vec(), matching: cleaver::Match::StartsWith },
    };
    let res = (|| -> Result<()> {
        let split = cleaver::chunk_file(input, &dest, &cfg).context("chunking input with Cleaver-X")?;
        eprintln!("[sentinel-x] Cleaver-X split the input into {} chunk(s)", split.chunks.len());
        let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
        for (i, chunk) in split.chunks.iter().enumerate() {
            let mut contigs = read_fasta(chunk)?;
            if opts.min_len > 0 {
                contigs.retain(|c| c.len() >= opts.min_len);
            }
            for c in &contigs {
                if !seen.insert(c.id.clone()) {
                    bail!("duplicate contig id '{}' in the input", c.id);
                }
            }
            if contigs.is_empty() {
                continue;
            }
            eprintln!("[sentinel-x] chunk {}/{}: {} contigs", i + 1, split.chunks.len(), contigs.len());
            let gs = run_per_contig(&contigs, opts, dbroot, report, stats.contigs == 0)?;
            stats.add(&contigs, &gs);
        }
        Ok(())
    })();
    let _ = std::fs::remove_dir_all(&dest);
    res
}

#[cfg(not(feature = "reads"))]
fn run_chunked(_opts: &Options, _dbroot: &DbRoot, _report: &mut QcReport, _stats: &mut StatsAcc) -> Result<()> {
    bail!("--chunk-size requires the `reads` feature (Cleaver-X)")
}

/// Assembly statistics accumulated over all analysed contigs.
#[derive(Default)]
struct StatsAcc {
    contigs: usize,
    lengths: Vec<usize>,
    gc: u64,
    at: u64,
    covered: usize,
    genes: usize,
    aa: usize,
}

impl StatsAcc {
    fn add(&mut self, contigs: &[Contig], gs: &GeneSet) {
        let per = gs.per_contig(contigs.len());
        for (c, idx) in contigs.iter().zip(per) {
            self.contigs += 1;
            self.lengths.push(c.len());
            for &b in &c.seq {
                match b {
                    b'G' | b'C' => self.gc += 1,
                    b'A' | b'T' => self.at += 1,
                    _ => {}
                }
            }
            let genes: Vec<crate::common::Gene> = idx.iter().map(|&i| gs.genes[i].clone()).collect();
            let (d, _) = crate::novel::trace::coding_architecture(c.len(), &genes);
            self.covered += (d * c.len() as f64).round() as usize;
        }
        self.genes += gs.genes.len();
        self.aa += gs.genes.iter().map(|g| g.aa.len()).sum::<usize>();
    }

    fn merge(&mut self, o: StatsAcc) {
        self.contigs += o.contigs;
        self.lengths.extend(o.lengths);
        self.gc += o.gc;
        self.at += o.at;
        self.covered += o.covered;
        self.genes += o.genes;
        self.aa += o.aa;
    }

    /// Size, contigs, N50 (standard definition), GC (% of A/C/G/T), coding
    /// density (union of gene spans), genes, mean protein length.
    fn table(&self) -> Vec<(String, String)> {
        let total: usize = self.lengths.iter().sum();
        vec![
            ("genome_size".into(), total.to_string()),
            ("num_contigs".into(), self.contigs.to_string()),
            ("contig_n50".into(), crate::common::n50(&self.lengths).to_string()),
            ("max_contig_length".into(), self.lengths.iter().max().copied().unwrap_or(0).to_string()),
            (
                "gc_percent".into(),
                format!("{:.2}", if self.gc + self.at > 0 { 100.0 * self.gc as f64 / (self.gc + self.at) as f64 } else { 0.0 }),
            ),
            ("coding_density".into(), format!("{:.3}", if total > 0 { self.covered as f64 / total as f64 } else { 0.0 })),
            ("total_coding_sequences".into(), self.genes.to_string()),
            (
                "average_gene_length_aa".into(),
                format!("{:.1}", if self.genes > 0 { self.aa as f64 / self.genes as f64 } else { 0.0 }),
            ),
        ]
    }
}

/// Paths that must not be overwritten by an automatic HTML report.
fn is_special_path(p: &str) -> bool {
    p == "-" || p.starts_with("/dev/") || p.starts_with("/proc/")
}

fn write_outputs(
    opts: &Options,
    report: &QcReport,
    stats: &StatsAcc,
    per_genome: &[crate::html_report::GenomeStats],
) -> Result<()> {
    let tsv = report.to_tsv();
    match opts.output.as_deref() {
        Some(path) if path != "-" => {
            std::fs::write(path, &tsv).with_context(|| format!("writing '{path}'"))?;
            eprintln!("[sentinel-x] wrote {} records to {path}", report.records.len());
        }
        _ => {
            let mut out = std::io::stdout().lock();
            out.write_all(tsv.as_bytes())?;
        }
    }
    if let Some(jpath) = &opts.json {
        std::fs::write(jpath, serde_json::to_string_pretty(&report_to_json(report))?)
            .with_context(|| format!("writing JSON '{jpath}'"))?;
        eprintln!("[sentinel-x] wrote JSON report to {jpath}");
    }
    // HTML: --html, else next to -o (never over the table or the JSON, never
    // for special files)
    let html_path: Option<String> = opts.html.clone().or_else(|| {
        if opts.no_html {
            return None;
        }
        let o = opts.output.as_deref()?;
        if is_special_path(o) {
            return None;
        }
        let mut p = PathBuf::from(o);
        p.set_extension("html");
        let p = p.to_string_lossy().into_owned();
        if p == o || opts.json.as_deref() == Some(p.as_str()) {
            return None;
        }
        Some(p)
    });
    if let Some(hpath) = html_path {
        if opts.output.as_deref() == Some(hpath.as_str()) || opts.json.as_deref() == Some(hpath.as_str()) {
            bail!("--html {hpath} would overwrite another output");
        }
        let html = crate::html_report::render(report, &opts.input, opts.domain, &stats.table(), per_genome);
        std::fs::write(&hpath, html).with_context(|| format!("writing HTML '{hpath}'"))?;
        eprintln!("[sentinel-x] wrote HTML report to {hpath}");
    }
    Ok(())
}

/// The unified report as JSON (extra fields in their recorded order).
fn report_to_json(report: &QcReport) -> serde_json::Value {
    let recs: Vec<serde_json::Value> = report
        .records
        .iter()
        .map(|r| {
            let mut extra = serde_json::Map::new();
            for (k, v) in &r.extra {
                extra.insert(k.clone(), serde_json::Value::String(v.clone()));
            }
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
                "extra": serde_json::Value::Object(extra),
            })
        })
        .collect();
    serde_json::json!({ "sentinel_x_version": crate::common::VERSION, "records": recs })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_and_special_paths() {
        assert_eq!(genome_name(Path::new("/x/bin.1.fna.gz")), "bin.1");
        assert_eq!(genome_name(Path::new("TEST1.tst")), "TEST1");
        assert!(is_fasta_name(Path::new("a.fasta.gz")));
        assert!(!is_fasta_name(Path::new("a.txt")));
        assert!(is_special_path("/dev/null"));
        assert!(is_special_path("-"));
    }
}
