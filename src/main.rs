//! # Sentinel-X CLI
//!
//! Unified genome quality-control: one gene call, every validator. See
//! [`run`](crate) `--help` for usage.

// Index loops are the clearest form for the DP and table code.
#![allow(clippy::needless_range_loop)]

// ---- engines & validators (one crate, one binary) ----
mod busco;
mod checkm2;
mod checkv;
mod common;
mod db;
mod eukcc;
mod gbm;
mod hmmer;
mod novel;
mod protsearch;
mod viralverify;
// ---- CLI ----
mod doctor;
mod html_report;
mod pipeline;

use crate::common::{banner, Domain, VERSION};
use anyhow::{bail, Result};

const RUN_HELP: &str = r#"`sentinel-x run` OPTIONS
    -i, --input <FASTA|DIR>    Nucleotide FASTA (.fa/.fna[.gz]); a directory = one genome per file
    -d, --domain <D>           virus | bacteria | archaea | eukaryote | mixed      [default: mixed]
                               (`auto` is accepted as an alias of `mixed`)
        --db <DIR>             Database root                     [default: $SENTINEL_X_DB or ~/.sentinel-x/db]
    -o, --output <TSV>         Unified QC table ("-" = stdout)                   [default: stdout]
        --json <PATH>          Also write the table as JSON
        --html <PATH>          HTML report path        [default: next to -o; never for stdout or /dev/*]
        --no-html              Do not write the automatic HTML report
        --outdir <DIR>         Also write each tool's own tables (CheckV, CheckM2 quality_report.tsv,
                               BUSCO full_table.tsv, viralVerify result table) under DIR
    -t, --threads <N>          Worker threads for every stage (0 = all CPUs)      [default: 1]
        --min-len <N>          Ignore contigs shorter than N bases                [default: 0]
        --no-novel             Skip the novel TRACE / CHIMERA-SCAN / SAGE modules
      virus / mixed
        --proteins <FAA>       CheckV: use these Prodigal-format proteins instead of calling genes
        --checkv-strict        Reproduce CheckV 1.1.1 exactly, including its grouping of DIAMOND
                               hits by output order (AAI computed per fragment; see AUDIT.md)
        --viralverify          Run viralVerify (Pfam-A search; needs `sentinel-x download pfam`)
        --chunk-size <BYTES>   Process a huge FASTA in Cleaver-X chunks (bounded memory)
      bacteria / archaea
        --ttable <N>           CheckM2 translation table (default: choose 4 or 11 as CheckM2 does)
        --ko <FILE>            KO counts instead of the DIAMOND-compatible search ("K00001" per
                               protein, or "K00001 <count>")
        --checkm2-comp <GBM>   Replacement LightGBM completeness model (with --checkm2-cont)
        --checkm2-cont <GBM>   Replacement LightGBM contamination model
      genome-level (bacteria / archaea / eukaryote)
        --lineage <NAME|DIR>   BUSCO lineage (e.g. bacteria_odb12.2) under <db>/busco, or a dataset dir
        --genes                Inputs are protein FASTA (.faa) files, one per genome (CheckM2 --genes;
                               BUSCO protein mode)
      eukaryote
        --markers <HMM>        EukCC marker set (every model must carry GA cutoffs)
"#;

fn print_help() {
    println!("{}", banner());
    println!(
        r#"Sentinel-X {ver} — unified genome quality control
Completeness, contamination and structural checks for viruses, bacteria,
archaea and eukaryotes: CheckV, CheckM2, BUSCO, EukCC and viralVerify
re-implemented in one Rust binary, plus database-free TRACE, CHIMERA-SCAN
and SAGE.

USAGE
    sentinel-x <COMMAND> [OPTIONS]

ANALYSIS
    run                    Assess a FASTA of contigs/genomes (options below)

DATABASES
    doctor [--build] [--verify]   Check the environment and reference databases
    download <db> [--lineage N]   Download, verify and install a database
                                  (<db> = all | checkv | checkm2 | pfam | busco)
    convert <in.dmnd> <out.faa>   Extract a DIAMOND database to protein FASTA
    blastp -q <faa> -d <db>       DIAMOND-compatible protein search (BLAST tabular output)

INFORMATION
    info | version | banner | help

CLUSTER
    cluster-status <addr> [port]  Query a HydraMPP cluster (feature `scale`)

{run_help}
DATABASES (detail)
    Root: --db <DIR>  >  $SENTINEL_X_DB  >  ~/.sentinel-x/db
    CheckV   portal.nersc.gov/CheckV            (CURRENT_RELEASE.txt → release tarball)
    CheckM2  doi.org/10.5281/zenodo.14897628    (MD5 + CONTENTS.json + SHA-256 verified)
    Pfam-A   ftp.ebi.ac.uk current_release      (MD5 verified)
    BUSCO    busco-data.ezlab.org/v5/data       (file_versions.tsv + MD5; needs --lineage)

ROUTING BY DOMAIN
    virus / mixed      per contig: CheckV, TRACE, CHIMERA-SCAN (+ viralVerify)
    bacteria / archaea per genome: CheckM2, SAGE (+ BUSCO); CHIMERA-SCAN per contig
    eukaryote          per genome: EukCC (--markers), BUSCO (--lineage), SAGE; CHIMERA-SCAN per contig

EXAMPLES
    sentinel-x download checkv
    sentinel-x run -i phages.fna -d virus -o viral.tsv --outdir viral_out -t 8
    sentinel-x download checkm2
    sentinel-x run -i bins/ -d bacteria -o bins.tsv --outdir bins_out -t 16
    sentinel-x download busco --lineage bacteria_odb12.2
    sentinel-x run -i genome.fna -d bacteria --lineage bacteria_odb12.2 -o qc.tsv

Docs: github.com/raw-lab/sentinel-x   ·   License: CC-BY-NC-4.0
"#,
        ver = VERSION,
        run_help = RUN_HELP
    );
}

fn print_version() {
    println!("{}", banner());
    let feat = |on: bool| if on { "on" } else { "off" };
    println!(
        r#"sentinel-x {ver}

build
    edition           Rust 2021 (MSRV 1.75)
    profile           {profile}
    features          embed-models={em}  reads={reads}  scale={scale}  rustyomestats={ro}

components
    validators        CheckV 1.1.1 · CheckM2 1.1.0 · BUSCO 6.1.0 · EukCC 2 · viralVerify 1.1
    novel modules     TRACE · CHIMERA-SCAN · SAGE
    engines           dynamite-orfs · HMMER 3.4-compatible search · DIAMOND 2.1-compatible search
                      · LightGBM evaluator · Cleaver-X · HydraMPP · rayon

License: CC-BY-NC-4.0   ·   github.com/raw-lab/sentinel-x"#,
        ver = VERSION,
        profile = if cfg!(debug_assertions) { "debug" } else { "release" },
        em = feat(cfg!(feature = "embed-models")),
        reads = feat(cfg!(feature = "reads")),
        scale = feat(cfg!(feature = "scale")),
        ro = feat(cfg!(feature = "rustyomestats")),
    );
}

fn print_info() {
    println!("{}", banner());
    println!(
        r#"Sentinel-X {ver} — genome quality control for every domain of life

REFERENCE VALIDATORS (re-implemented in Rust)
    CheckV        viral completeness (AAI + HMM), host contamination, proviruses,
                  complete genomes (DTR/ITR), quality tiers — CheckV 1.1.1 tables
    CheckM2       completeness/contamination from 21,241 KEGG features: LightGBM
                  general models, neural-network specific model, model selection
    BUSCO         single-copy orthologs (Single / Duplicated / Fragmented / Missing)
    EukCC         marker-set completeness/contamination (EukCC 2 scoring)
    viralVerify   virus / plasmid / chromosome classification from Pfam domains

NOVEL MODULES (database-free)
    TRACE         terminal-repeat & coding-architecture confidence for phages
    CHIMERA-SCAN  compositional junctions inside a contig (inter-species chimeras)
    SAGE          contigs off a bin's core composition (reference-free contamination)

ENGINES
    dynamite-orfs   gene calling (Prodigal / prodigal-gv style presets)
    HMMER engine    hmmsearch 3.4-compatible scores, E-values, domains, GA/NC/TC
    protsearch      DIAMOND 2.1.9-compatible blastp (seeds, masking, X-drop, banded SW)
    gbm             LightGBM model evaluator (CheckM2's models)
    Cleaver-X       chunked reading of huge inputs            [feature: reads]
    HydraMPP        cluster status                            [feature: scale]
    rustyomestats   assembly-statistics cross-check           [feature: rustyomestats]

Run `sentinel-x help` for usage and `sentinel-x doctor` for database status."#,
        ver = VERSION
    );
}

fn main() {
    // Restore default SIGPIPE handling so that piping Sentinel-X's output into a
    // consumer that closes early (e.g. `sentinel-x help | head`) terminates the
    // process quietly, like every other Unix tool, instead of turning the
    // broken pipe into a panic/abort (Rust ignores SIGPIPE by default).
    #[cfg(unix)]
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }

    if let Err(e) = real_main() {
        eprintln!("sentinel-x: error: {e:#}");
        std::process::exit(1);
    }
}

fn real_main() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let cmd = args.first().map(String::as_str).unwrap_or("help");

    match cmd {
        "banner" => {
            println!("{}", banner());
            Ok(())
        }
        "version" | "--version" | "-V" => {
            print_version();
            Ok(())
        }
        "info" | "--info" => {
            print_info();
            Ok(())
        }
        "help" | "--help" | "-h" | "" => {
            print_help();
            Ok(())
        }
        "cluster-status" => cluster_status(&args[1..]),
        "doctor" => doctor::run(&args[1..]),
        "download" => doctor::download_cmd(&args[1..]),
        "convert" => convert_cmd(&args[1..]),
        "blastp" => blastp_cmd(&args[1..]),
        "run" => match Options::parse(&args[1..])? {
            Some(opts) => pipeline::run(opts),
            None => Ok(()),
        },
        other => {
            bail!("unknown command '{other}'. Try `sentinel-x help`.");
        }
    }
}

/// `sentinel-x convert <in.dmnd> <out.fasta>` — extract a DIAMOND database to a
/// plain protein FASTA (headers preserved), for use with other tools.
/// Sentinel-X itself searches `.dmnd` files in place; no DIAMOND binary is
/// invoked.
fn convert_cmd(args: &[String]) -> Result<()> {
    if args.iter().any(|a| a == "-h" || a == "--help") {
        println!("usage: sentinel-x convert <in.dmnd> <out.fasta>\n");
        println!("Extract a DIAMOND (.dmnd) database to protein FASTA (headers kept, e.g. CheckM2's `accession~KO`).");
        println!("Sentinel-X reads .dmnd files directly; conversion is only needed for other tools.");
        return Ok(());
    }
    if args.len() != 2 {
        bail!("usage: sentinel-x convert <in.dmnd> <out.fasta>");
    }
    let input = std::path::Path::new(&args[0]);
    let output = std::path::Path::new(&args[1]);
    eprintln!("Converting {} -> {}", input.display(), output.display());
    let conv = crate::db::dmnd::dmnd_to_fasta(input, output)?;
    eprintln!("Wrote {} sequences ({} residues) to {}", conv.sequences, conv.letters, output.display());
    Ok(())
}

/// `sentinel-x blastp`: the DIAMOND-compatible protein search on its own
/// (BLAST tabular output, DIAMOND's formatting and order).
fn blastp_cmd(args: &[String]) -> Result<()> {
    let usage = "usage: sentinel-x blastp -q <queries.faa> -d <db.dmnd|db.faa> [--preset checkv|checkm2] \
                 [--evalue E] [--id P] [--query-cover P] [--subject-cover P] [-k N] [-t N]";
    let mut q = None;
    let mut db = None;
    let mut p = crate::protsearch::Params::default();
    let mut threads = 1usize;
    let mut i = 0;
    let val = |i: &mut usize| -> Result<String> {
        *i += 1;
        args.get(*i).cloned().ok_or_else(|| anyhow::anyhow!("{usage}"))
    };
    while i < args.len() {
        match args[i].as_str() {
            "-q" | "--query" => q = Some(val(&mut i)?),
            "-d" | "--db" => db = Some(val(&mut i)?),
            "--preset" => {
                p = match val(&mut i)?.as_str() {
                    "checkv" => crate::protsearch::Params::checkv(),
                    "checkm2" => crate::protsearch::Params::checkm2(),
                    other => bail!("unknown preset '{other}'"),
                }
            }
            "--evalue" | "-e" => p.evalue = val(&mut i)?.parse()?,
            "--id" => p.min_identity = val(&mut i)?.parse()?,
            "--query-cover" => p.query_cover = val(&mut i)?.parse()?,
            "--subject-cover" => p.subject_cover = val(&mut i)?.parse()?,
            "-k" | "--max-target-seqs" => p.max_target_seqs = val(&mut i)?.parse()?,
            "-t" | "--threads" => threads = val(&mut i)?.parse()?,
            "-h" | "--help" => {
                println!("{usage}");
                return Ok(());
            }
            other => bail!("unknown option '{other}'\n{usage}"),
        }
        i += 1;
    }
    let (Some(q), Some(db)) = (q, db) else { bail!("{usage}") };
    let _ = rayon::ThreadPoolBuilder::new().num_threads(threads.max(1)).build_global();
    let queries: Vec<(String, Vec<u8>)> = crate::common::read_fasta(&q)?.into_iter().map(|c| (c.id, c.seq)).collect();
    let (hits, st) = crate::protsearch::search_path(&queries, std::path::Path::new(&db), &p)?;
    let mut out = std::io::BufWriter::new(std::io::stdout().lock());
    use std::io::Write;
    for h in &hits {
        writeln!(out, "{}", h.outfmt6(&queries[h.query].0))?;
    }
    out.flush()?;
    eprintln!("[sentinel-x blastp] {} queries, {} hits, {} reference letters", queries.len(), st.hits, st.db_letters);
    Ok(())
}

/// Parsed options for `sentinel-x run`.
#[derive(Debug, Clone)]
pub struct Options {
    pub input: String,
    pub domain: Domain,
    pub output: Option<String>,
    pub json: Option<String>,
    pub html: Option<String>,
    pub no_html: bool,
    pub outdir: Option<String>,
    pub min_len: usize,
    pub threads: usize,
    pub db_root: Option<String>,
    pub novel: bool,
    pub proteins: Option<String>,
    /// Inputs are protein FASTA files (CheckM2 `--genes`).
    pub genes: bool,
    pub checkv_strict: bool,
    pub viralverify: bool,
    pub chunk_size: Option<u64>,
    pub ttable: Option<u32>,
    pub ko_file: Option<String>,
    pub checkm2_comp: Option<String>,
    pub checkm2_cont: Option<String>,
    pub lineage: Option<String>,
    pub markers: Option<String>,
    /// Accepted for compatibility; has no effect.
    pub gpu: bool,
}

impl Options {
    fn parse(args: &[String]) -> Result<Option<Options>> {
        let mut o = Options {
            input: String::new(),
            domain: Domain::Mixed,
            output: None,
            json: None,
            html: None,
            no_html: false,
            outdir: None,
            min_len: 0,
            threads: 1,
            db_root: None,
            novel: true,
            proteins: None,
            genes: false,
            checkv_strict: false,
            viralverify: false,
            chunk_size: None,
            ttable: None,
            ko_file: None,
            checkm2_comp: None,
            checkm2_cont: None,
            lineage: None,
            markers: None,
            gpu: false,
        };
        let mut i = 0;
        let value = |i: &mut usize, flag: &str| -> Result<String> {
            *i += 1;
            match args.get(*i) {
                Some(v) if v == "-" || !v.starts_with('-') => Ok(v.clone()),
                _ => bail!("{flag} requires a value"),
            }
        };
        let number = |s: String, flag: &str| -> Result<u64> {
            s.parse::<u64>().map_err(|_| anyhow::anyhow!("{flag}: '{s}' is not a non-negative integer"))
        };
        while i < args.len() {
            let a = args[i].as_str();
            match a {
                "-h" | "--help" => {
                    println!("{}", RUN_HELP);
                    return Ok(None);
                }
                "-i" | "--input" => o.input = value(&mut i, a)?,
                "-d" | "--domain" => {
                    let v = value(&mut i, a)?;
                    o.domain = if v.eq_ignore_ascii_case("auto") {
                        Domain::Mixed
                    } else {
                        Domain::parse(&v).ok_or_else(|| {
                            anyhow::anyhow!("unknown domain '{v}' (virus | bacteria | archaea | eukaryote | mixed)")
                        })?
                    };
                }
                "-o" | "--output" => o.output = Some(value(&mut i, a)?),
                "--json" => o.json = Some(value(&mut i, a)?),
                "--html" => o.html = Some(value(&mut i, a)?),
                "--no-html" => o.no_html = true,
                "--outdir" => o.outdir = Some(value(&mut i, a)?),
                "--min-len" => o.min_len = number(value(&mut i, a)?, a)? as usize,
                "-t" | "--threads" => o.threads = number(value(&mut i, a)?, a)? as usize,
                "--db" => o.db_root = Some(value(&mut i, a)?),
                "--no-novel" => o.novel = false,
                "--proteins" => o.proteins = Some(value(&mut i, a)?),
                "--viralverify" => o.viralverify = true,
                "--checkv-strict" => o.checkv_strict = true,
                "--genes" => o.genes = true,
                "--chunk-size" => {
                    let n = number(value(&mut i, a)?, a)?;
                    if n == 0 {
                        bail!("--chunk-size must be > 0");
                    }
                    o.chunk_size = Some(n);
                }
                "--ttable" => {
                    let t = number(value(&mut i, a)?, a)? as u32;
                    if dynamite::gencode::GeneticCode::new(t).is_none() {
                        bail!("--ttable {t}: unsupported translation table");
                    }
                    o.ttable = Some(t);
                }
                "--ko" => o.ko_file = Some(value(&mut i, a)?),
                "--checkm2-comp" => o.checkm2_comp = Some(value(&mut i, a)?),
                "--checkm2-cont" => o.checkm2_cont = Some(value(&mut i, a)?),
                "--lineage" => o.lineage = Some(value(&mut i, a)?),
                "--markers" => o.markers = Some(value(&mut i, a)?),
                "--gpu" => o.gpu = true,
                "--scale" => {
                    eprintln!("[sentinel-x] note: --scale is no longer needed (contigs are always processed in parallel)")
                }
                other if !other.starts_with('-') && o.input.is_empty() => o.input = other.to_string(),
                other => bail!("unknown option '{other}' (see `sentinel-x run --help`)"),
            }
            i += 1;
        }
        if o.input.is_empty() {
            bail!("no input FASTA given (use -i/--input)");
        }
        if o.checkm2_comp.is_some() != o.checkm2_cont.is_some() {
            bail!("--checkm2-comp and --checkm2-cont must be given together");
        }
        if o.genes && matches!(o.domain, Domain::Virus | Domain::Mixed) {
            bail!("--genes applies to genome-level domains; for CheckV give contigs with --proteins <faa>");
        }
        Ok(Some(o))
    }
}

#[cfg(feature = "scale")]
fn cluster_status(args: &[String]) -> Result<()> {
    use std::time::Duration;
    let addr = args.first().cloned().unwrap_or_else(|| "127.0.0.1".to_string());
    let port: u16 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(8765);
    match hydra_mpp::query_status(&addr, port, Duration::from_secs(5)) {
        Ok(s) => {
            println!("hydra-mpp cluster @ {addr}:{port}");
            println!("{s:#?}");
            Ok(())
        }
        Err(e) => {
            bail!("could not reach hydra-mpp cluster at {addr}:{port}: {e}");
        }
    }
}

#[cfg(not(feature = "scale"))]
fn cluster_status(_args: &[String]) -> Result<()> {
    bail!("cluster-status requires the `scale` feature (rebuild with --features scale)");
}
