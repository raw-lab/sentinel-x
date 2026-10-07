//! # Sentinel CLI
//!
//! Unified genome quality-control: one gene call, every validator. See
//! [`run`](crate) `--help` for usage.

mod doctor;
mod pipeline;

use anyhow::{bail, Result};
use sentinel_core::{banner, Domain, VERSION};

fn print_help() {
    println!("{}", banner());
    println!(
        r#"Sentinel {ver} — unified genome quality control
Completeness, contamination & structural validity for viruses, bacteria,
archaea and eukaryotes — one gene call, every validator.

USAGE
    sentinel <COMMAND> [OPTIONS]

ANALYSIS
    run                    Assess a FASTA of contigs/genomes (options below)

DATABASES
    doctor [--build]       Check environment + reference databases;
                           --build creates the database folder + manifest
    download <db>          Download a reference database into the DB root
                           (use `download all` to fetch every database at once)

INFORMATION
    info                   Capabilities, validators, novel modules & credits
    version                Version and build details
    banner                 Print the Sentinel banner
    help                   Show this help

CLUSTER
    cluster-status <addr>  Query a hydra-mpp cluster (requires `scale` feature)

`sentinel run` OPTIONS
    -i, --input <FASTA>        Input nucleotide FASTA (.fa/.fna[.gz])           [required]
    -d, --domain <D>           virus | bacteria | archaea | eukaryote | auto    [default: auto]
        --db <DIR>             Reference-database root (auto-discovers models)
    -o, --output <TSV>         Write the QC table here                          [default: stdout]
        --json <PATH>          Also write a JSON report
        --min-len <N>          Ignore contigs shorter than N bases             [default: 0]
    -t, --threads <N>          Worker threads for gene calling                  [default: 1]
        --checkm2-comp <GBM>   CheckM2 completeness LightGBM model (.gbm)
        --checkm2-cont <GBM>   CheckM2 contamination LightGBM model (.gbm)
        --ko <FILE>            KOs present (one KO id per line) for CheckM2
        --gpu                  Use the GPU backend if available (else CPU)
        --no-novel             Skip the novel TRACE/CHIMERA/SAGE modules
        --chunk-size <BYTES>   Pre-split huge inputs (requires `reads` feature)
        --scale                Parallelise contigs across threads

DATABASES (detail)
    <db> is one of:  all | checkv | checkm2 | pfam | kofam | busco
    `download all` fetches every database (BUSCO included only with --lineage)
    DB root resolves in order:   --db <DIR>   >   $SENTINEL_DB   >   ~/.sentinel/db
    `download busco` also takes   --lineage <name>   (e.g. bacteria_odb10)
    Sources: CheckV (portal.nersc.gov), CheckM2 (Zenodo), Pfam-A (EBI),
             KOfam (genome.jp), BUSCO (busco-data.ezlab.org).

ROUTING BY DOMAIN
    virus / mixed      per-contig CheckV structural + TRACE + CHIMERA-SCAN
    bacteria / archaea genome-level CheckM2 (+ SAGE) + per-contig CHIMERA-SCAN
    eukaryote          genome-level SAGE/composition (BUSCO needs a lineage)

OUTPUT
    One row per (contig|genome) x validator with: length, domain, gene_count,
    completeness, contamination, quality tier, method, warnings, and
    validator-specific fields in the `extra` column. `--json` mirrors this.

EXAMPLES
    # Viral contigs — database-free structural + novel QC
    sentinel run -i phage.fna -d virus -o viral.tsv --json viral.json

    # Bacterial genome — build the DB folder, fetch CheckM2, auto-discover models
    sentinel doctor --build
    sentinel download checkm2
    sentinel run -i genome.fna -d bacteria --db ~/.sentinel/db -o qc.tsv

    # See what's installed and which validators are ready
    sentinel doctor

Reference validators use their databases when present; structural CheckV
signals and all novel modules always run database-free.
Docs: github.com/raw-lab   ·   License: GPL-3.0-or-later
"#,
        ver = VERSION
    );
}

fn print_version() {
    println!("{}", banner());
    let backend = sentinel_gpu::Backend::best();
    let feat = |on: bool| if on { "on" } else { "off" };
    println!(
        r#"sentinel {ver}

build
    edition           Rust 2021 (MSRV 1.75)
    profile           {profile}
    compute backend   {backend}
    features          reads={reads}  scale={scale}  gpu={gpu}  polars={polars}

components
    validators        CheckV · CheckM2 · BUSCO · ViralVerify
    novel modules     TRACE · CHIMERA-SCAN · SAGE
    engines           dynamite · sentinel-hmmer · sentinel-gbm · cleaver · hydra-mpp

License: GPL-3.0-or-later   ·   github.com/raw-lab
Run `sentinel info` for capabilities, `sentinel help` for usage."#,
        ver = VERSION,
        profile = if cfg!(debug_assertions) { "debug" } else { "release" },
        backend = backend.describe(),
        reads = feat(cfg!(feature = "reads")),
        scale = feat(cfg!(feature = "scale")),
        gpu = feat(cfg!(feature = "gpu")),
        polars = feat(cfg!(feature = "polars")),
    );
}

fn print_info() {
    println!("{}", banner());
    println!(
        r#"Sentinel {ver} — genome integrity sentinel for every domain of life

WHAT IT DOES
    Assesses completeness, contamination and structural validity of assembled
    contigs and genomes — viruses, bacteria, archaea and eukaryotes — behind a
    single CLI, a single gene call, and a single unified report.

REFERENCE VALIDATORS  (the wrapped tools, re-implemented in Rust)
    CheckV        viral completeness/contamination — terminal repeats (DTR/ITR),
                  provirus host-boundary, HMM-based completeness
    CheckM2       bacterial/archaeal prediction — a 21,241-feature KEGG vector
                  scored by LightGBM for completeness & contamination
    BUSCO         single-copy orthologs — Single/Duplicated/Fragmented/Missing
    ViralVerify   virus / plasmid / chromosome classification from Pfam content

NOVEL MODULES  (unique to Sentinel — always database-free)
    TRACE         terminal-repeat & architecture completeness for phages
    CHIMERA-SCAN  compositional-coherence contamination / mis-assembly detection
    SAGE          single-copy anchor-gene & codon-coherence self-consistency

ENGINES
    dynamite          ORF calling, shared across every validator
    sentinel-hmmer    pure-Rust Plan-7 profile-HMM search (Forward/Viterbi)
    sentinel-gbm      LightGBM evaluator (exact CheckM2 reproduction)
    cleaver           FASTA/FASTQ/BAM I/O + chunking          [feature: reads]
    hydra-mpp         distributed scaling across nodes         [feature: scale]
    wgpu              optional GPU compute kernels             [feature: gpu]

VALIDATED AGAINST THE ORIGINAL TOOLS
    CheckV terminal repeats   identical on all 46 CheckV test contigs
    CheckM2 metadata          22/22 amino-acid/length/CDS features identical
    CheckM2 KEGG aggregation  exact on all 208 pathway/module/category values
    CheckM2 LightGBM          completeness & contamination bit-identical (0.00e+00)

DATABASES  (managed under $SENTINEL_DB, default ~/.sentinel/db)
    sentinel doctor --build       create the folder skeleton + manifest
    sentinel download all         fetch every reference database at once
    sentinel download <db>        fetch one of checkv|checkm2|pfam|kofam|busco
    sentinel doctor               verify installs + per-validator readiness

Author: RAW Lab (github.com/raw-lab)    License: GPL-3.0-or-later
Run `sentinel help` for full usage and examples."#,
        ver = VERSION
    );
}

fn main() {
    // Restore default SIGPIPE handling so that piping Sentinel's output into a
    // consumer that closes early (e.g. `sentinel help | head`) terminates the
    // process quietly, like every other Unix tool, instead of turning the
    // broken pipe into a panic/abort (Rust ignores SIGPIPE by default).
    #[cfg(unix)]
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }

    if let Err(e) = real_main() {
        eprintln!("sentinel: error: {e:#}");
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
        "run" => {
            let opts = Options::parse(&args[1..])?;
            pipeline::run(opts)
        }
        other => {
            bail!("unknown command '{other}'. Try `sentinel help`.");
        }
    }
}

/// Parsed options for `sentinel run`.
#[derive(Debug, Clone)]
pub struct Options {
    pub input: String,
    pub domain: Domain,
    pub output: Option<String>,
    pub json: Option<String>,
    pub min_len: usize,
    pub threads: usize,
    pub db_root: Option<String>,
    pub checkm2_comp: Option<String>,
    pub checkm2_cont: Option<String>,
    pub ko_file: Option<String>,
    pub gpu: bool,
    pub novel: bool,
    pub chunk_size: Option<u64>,
    pub scale: bool,
}

impl Options {
    fn parse(args: &[String]) -> Result<Options> {
        let mut o = Options {
            input: String::new(),
            domain: Domain::Mixed,
            output: None,
            json: None,
            min_len: 0,
            threads: 1,
            checkm2_comp: None,
            checkm2_cont: None,
            ko_file: None,
            gpu: false,
            novel: true,
            chunk_size: None,
            db_root: None,
            scale: false,
        };
        let mut i = 0;
        let next = |i: &mut usize, flag: &str| -> Result<String> {
            *i += 1;
            args.get(*i)
                .cloned()
                .ok_or_else(|| anyhow::anyhow!("{flag} requires a value"))
        };
        while i < args.len() {
            let a = args[i].as_str();
            match a {
                "-i" | "--input" => o.input = next(&mut i, a)?,
                "-d" | "--domain" => {
                    let v = next(&mut i, a)?;
                    o.domain = Domain::parse(&v)
                        .ok_or_else(|| anyhow::anyhow!("unknown domain '{v}'"))?;
                }
                "-o" | "--output" => o.output = Some(next(&mut i, a)?),
                "--json" => o.json = Some(next(&mut i, a)?),
                "--min-len" => o.min_len = next(&mut i, a)?.parse()?,
                "-t" | "--threads" => o.threads = next(&mut i, a)?.parse()?,
                "--db" => o.db_root = Some(next(&mut i, a)?),
                "--checkm2-comp" => o.checkm2_comp = Some(next(&mut i, a)?),
                "--checkm2-cont" => o.checkm2_cont = Some(next(&mut i, a)?),
                "--ko" => o.ko_file = Some(next(&mut i, a)?),
                "--gpu" => o.gpu = true,
                "--no-novel" => o.novel = false,
                "--chunk-size" => o.chunk_size = Some(next(&mut i, a)?.parse()?),
                "--scale" => o.scale = true,
                other if !other.starts_with('-') && o.input.is_empty() => {
                    o.input = other.to_string();
                }
                other => bail!("unknown option '{other}' (see `sentinel run --help`)"),
            }
            i += 1;
        }
        if o.input.is_empty() {
            bail!("no input FASTA given (use -i/--input)");
        }
        Ok(o)
    }
}

#[cfg(feature = "scale")]
fn cluster_status(args: &[String]) -> Result<()> {
    use std::time::Duration;
    let addr = args
        .first()
        .cloned()
        .unwrap_or_else(|| "127.0.0.1".to_string());
    let port: u16 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(8765);
    match hydra_mpp_core::query_status(&addr, port, Duration::from_secs(5)) {
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
