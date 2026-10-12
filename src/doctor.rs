//! `sentinel-x doctor` (environment and database diagnostics) and
//! `sentinel-x download`.

use anyhow::{bail, Result};

use crate::db::{download, Database, DbRoot, Status};

fn take_value(args: &[String], i: &mut usize, flag: &str) -> Result<String> {
    *i += 1;
    match args.get(*i) {
        Some(v) if !v.starts_with('-') => Ok(v.clone()),
        _ => bail!("{flag} requires a value"),
    }
}

/// Run `sentinel-x doctor [--db DIR] [--build] [--verify]`.
pub fn run(args: &[String]) -> Result<()> {
    let mut db_path: Option<String> = None;
    let mut build = false;
    let mut verify = false;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--db" => db_path = Some(take_value(args, &mut i, "--db")?),
            "--build" | "--build-db-dir" => build = true,
            "--verify" => verify = true,
            "-h" | "--help" => {
                println!("usage: sentinel-x doctor [--db DIR] [--build] [--verify]\n");
                println!("  --build    create the database folder skeleton + manifest");
                println!("  --verify   re-hash installed databases against their recorded checksums");
                return Ok(());
            }
            other => bail!("unknown doctor option '{other}'"),
        }
        i += 1;
    }

    println!("{}", crate::common::banner());
    println!("Sentinel-X doctor — environment & database check\n");

    println!("[ runtime ]");
    println!("  sentinel-x version     {}", crate::common::VERSION);
    println!("  threads available      {}", std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1));
    match download::fetcher_available() {
        Some(t) => println!("  database downloader    {t} (found)"),
        None => println!("  database downloader    MISSING — install curl or wget to fetch databases"),
    }
    println!("  engines                built in (HMMER 3.4-compatible profile search, DIAMOND-compatible protein search)");

    println!("\n[ bundled data (in the binary) ]");
    println!("  viralVerify classifier table   embedded");
    println!("  CheckM2 KEGG definitions       embedded");
    println!(
        "  CheckM2 LightGBM models        {}",
        if cfg!(feature = "embed-models") { "embedded" } else { "not embedded (pass --checkm2-comp/--checkm2-cont)" }
    );

    let root = DbRoot::resolve(db_path.as_deref());
    println!("\n[ reference databases ]");
    println!("  root                   {}", root.root.display());
    if build {
        root.build_skeleton()?;
        println!("  (created folder skeleton + manifest.json)");
    }
    if !root.root.exists() {
        println!("  status                 root does not exist — run `sentinel-x doctor --build` or `sentinel-x download <db>`");
    }
    let mut ready = 0;
    for (db, st) in root.status_all() {
        if st.ok() {
            ready += 1;
        }
        let detail = match &st {
            Status::Installed { detail } => detail.clone(),
            Status::Unverified { detail } => format!("{detail} (no install record: placed by hand?)"),
            Status::Partial { missing } => format!("missing: {}", missing.join(", ")),
            Status::Empty => "folder present, nothing installed".into(),
            Status::NotInstalled => "not downloaded".into(),
        };
        println!("  {:<10} {:<8} used by {:<12} {}", db.id(), st.symbol(), db.used_by().join(","), detail);
    }
    if verify {
        println!("\n[ verification ]");
        if let Some(p) = root.checkm2_db() {
            if p.extension().map(|e| e == "dmnd").unwrap_or(false) {
                let sha = download::sha256_hex(&p)?;
                let ok = sha == crate::db::expected_checkm2_sha256();
                println!("  checkm2   {}  SHA-256 {}", if ok { "OK  " } else { "FAIL" }, sha);
            }
        }
        if let Some(d) = root.checkv_dir() {
            match crate::checkv::db::CheckvDb::load(&d) {
                Ok(db) => println!(
                    "  checkv    OK    {} HMM files, {} reference genomes, {} error-table cells",
                    db.hmm_files.len(),
                    db.reps.len(),
                    db.error.rates.len()
                ),
                Err(e) => println!("  checkv    FAIL  {e:#}"),
            }
        }
        for (name, dir) in root.busco_lineages() {
            match crate::busco::Lineage::load(&dir) {
                Ok(l) => println!("  busco     OK    {name}: {} BUSCOs ({})", l.total(), l.version),
                Err(e) => println!("  busco     FAIL  {name}: {e:#}"),
            }
        }
    }

    println!("\n[ validator readiness ]");
    let st = |db| root.status(db).ok();
    line("CheckV (viral completeness/contamination)", st(Database::CheckV), "terminal repeats and copy number run without it");
    line("CheckM2 (bacteria/archaea)", st(Database::CheckM2), "KO annotation needs the DIAMOND database");
    line("viralVerify (virus/plasmid/chromosome)", st(Database::Pfam), "needs Pfam-A; enable with --viralverify");
    line("BUSCO (single-copy orthologs)", st(Database::Busco), "needs a lineage; select with --lineage");
    println!("  {:<44} needs a marker set: --markers <set.hmm>", "EukCC (eukaryotes)");
    println!("  {:<44} READY    (database-free)", "Novel modules (TRACE / CHIMERA-SCAN / SAGE)");
    println!(
        "\nSummary: {ready}/{} reference databases usable. Fetch them with `sentinel-x download all` \
         (add --lineage <name> for BUSCO) or `sentinel-x download <{}>`.",
        Database::all().len(),
        Database::all().iter().map(|d| d.id()).collect::<Vec<_>>().join("|")
    );
    Ok(())
}

fn line(label: &str, ready: bool, note: &str) {
    println!("  {:<44} {}   ({})", label, if ready { "READY   " } else { "needs DB" }, note);
}

/// `sentinel-x download <db|all> [--db DIR] [--lineage NAME]`
pub fn download_cmd(args: &[String]) -> Result<()> {
    let ids: Vec<&str> = Database::all().iter().map(|d| d.id()).collect();
    if args.is_empty() || args.iter().any(|a| a == "-h" || a == "--help") {
        println!("usage: sentinel-x download <all|{}> [--db DIR] [--lineage NAME]", ids.join("|"));
        println!("\n  BUSCO lineages are named as in BUSCO's file_versions.tsv, e.g. bacteria_odb12.2");
        if args.is_empty() {
            bail!("no database given");
        }
        return Ok(());
    }
    let target = args[0].to_ascii_lowercase();
    let mut db_path: Option<String> = None;
    let mut lineages: Vec<String> = Vec::new();
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--db" => db_path = Some(take_value(args, &mut i, "--db")?),
            "--lineage" => lineages.push(take_value(args, &mut i, "--lineage")?),
            other => bail!("unknown download option '{other}'"),
        }
        i += 1;
    }
    let root = DbRoot::resolve(db_path.as_deref());
    println!("{}", crate::common::banner());
    if download::fetcher_available().is_none() {
        bail!("downloading needs curl or wget on PATH");
    }

    if target == "all" {
        return download_all(&root, &lineages);
    }
    if target == "kofam" {
        bail!("KOfam is no longer used: CheckM2 annotates KOs with its own DIAMOND database (`sentinel-x download checkm2`)");
    }
    let db = Database::parse(&target)
        .ok_or_else(|| anyhow::anyhow!("unknown database '{}' (expected: all|{})", target, ids.join("|")))?;
    if db == Database::Busco {
        if lineages.is_empty() {
            bail!("BUSCO needs --lineage <name> (e.g. bacteria_odb12.2, eukaryota_odb12)");
        }
        for l in &lineages {
            println!("Downloading BUSCO lineage {l} into {}\n", root.dir(db).display());
            root.install(db, Some(l))?;
        }
    } else {
        println!("Downloading {} into {}\n", db.display(), root.dir(db).display());
        root.install(db, None)?;
    }
    let st = root.status(db);
    println!("\nResult: {} is now {}", db.id(), st.symbol());
    Ok(())
}

/// Download every database, continuing past failures; BUSCO only with
/// `--lineage`. Databases already installed (verified) are skipped.
fn download_all(root: &DbRoot, lineages: &[String]) -> Result<()> {
    root.build_skeleton()?;
    println!("Downloading all reference databases into {}\n", root.root.display());
    let mut results: Vec<(String, String)> = Vec::new();
    let mut failed = false;
    for &db in Database::all() {
        if db == Database::Busco {
            if lineages.is_empty() {
                results.push(("busco".into(), "SKIPPED (needs --lineage)".into()));
                continue;
            }
            let have: Vec<String> = root.busco_lineages().into_iter().map(|x| x.0).collect();
            for l in lineages {
                if have.contains(l) {
                    results.push((format!("busco:{l}"), "OK (already present)".into()));
                    continue;
                }
                println!("==> BUSCO lineage {l}");
                match root.install(db, Some(l)) {
                    Ok(()) => results.push((format!("busco:{l}"), "OK".into())),
                    Err(e) => {
                        failed = true;
                        eprintln!("    failed: {e:#}\n");
                        results.push((format!("busco:{l}"), format!("FAILED: {e}")));
                    }
                }
            }
            continue;
        }
        if matches!(root.status(db), Status::Installed { .. }) {
            results.push((db.id().into(), "OK (already present)".into()));
            continue;
        }
        println!("==> {} ({})", db.display(), db.id());
        match root.install(db, None) {
            Ok(()) => results.push((db.id().into(), root.status(db).symbol().into())),
            Err(e) => {
                failed = true;
                eprintln!("    failed: {e:#}\n");
                results.push((db.id().into(), format!("FAILED: {e}")));
            }
        }
    }
    println!("\n==================== download all: summary ====================");
    for (db, status) in &results {
        println!("  {:<24} {}", db, status);
    }
    if failed {
        bail!("one or more databases did not install (see the summary above)");
    }
    Ok(())
}
