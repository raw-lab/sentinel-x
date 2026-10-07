//! `sentinel doctor` — environment and database diagnostics, plus the database
//! folder builder and downloader entry points.

use anyhow::{bail, Result};
use sentinel_db::{download, Database, DbRoot, Status};

/// External tools that improve (but are not required for) Sentinel runs.
const OPTIONAL_TOOLS: &[(&str, &str)] = &[
    ("diamond", "protein alignment for CheckV/CheckM2 reference search"),
    ("hmmsearch", "external HMMER (Sentinel has a built-in engine; this is a cross-check)"),
];

/// Run `sentinel doctor`.
///
/// Flags:
///   --db <DIR>     database root (else $SENTINEL_DB, else ~/.sentinel/db)
///   --build        create the database folder skeleton
pub fn run(args: &[String]) -> Result<()> {
    let mut db_path: Option<String> = None;
    let mut build = false;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--db" => {
                i += 1;
                db_path = Some(
                    args.get(i)
                        .cloned()
                        .ok_or_else(|| anyhow::anyhow!("--db requires a path"))?,
                );
            }
            "--build" | "--build-db-dir" => build = true,
            other => bail!("unknown doctor option '{other}'"),
        }
        i += 1;
    }

    println!("{}", sentinel_core::banner());
    println!("Sentinel doctor — environment & database check\n");

    // ---- 1. Toolchain / runtime ----
    println!("[ runtime ]");
    println!("  sentinel version       {}", sentinel_core::VERSION);
    let backend = sentinel_gpu::Backend::best();
    println!("  compute backend        {}", backend.describe());
    match download::fetcher_available() {
        Some(t) => println!("  database downloader    {t} (found)"),
        None => println!("  database downloader    MISSING — install curl or wget to fetch databases"),
    }

    // ---- 2. Optional external tools ----
    println!("\n[ optional external tools ]");
    for (tool, why) in OPTIONAL_TOOLS {
        let present = which(tool);
        println!(
            "  {:<22} {}   ({})",
            tool,
            if present { "found" } else { "not found" },
            why
        );
    }

    // ---- 3. Bundled data (always present) ----
    println!("\n[ bundled tables (shipped in the binary) ]");
    println!("  ViralVerify Pfam table  embedded");
    println!("  CheckM2 KEGG definitions embedded");

    // ---- 4. Database root ----
    let root = DbRoot::resolve(db_path.as_deref())?;
    println!("\n[ reference databases ]");
    println!("  root                   {}", root.root.display());
    if build {
        root.build_skeleton()?;
        println!("  (created folder skeleton + manifest.json)");
    }
    let exists = root.root.exists();
    if !exists {
        println!("  status                 root does not exist — run `sentinel doctor --build` or `sentinel download <db>`");
    }

    let mut ready = 0;
    let mut total = 0;
    for (db, st) in root.status_all() {
        total += 1;
        if st.ok() {
            ready += 1;
        }
        let detail = match &st {
            Status::Installed => String::new(),
            Status::Partial { missing } => format!("  (missing: {})", missing.join(", ")),
            Status::Empty => "  (folder present, no files)".to_string(),
            Status::NotInstalled => "  (not downloaded)".to_string(),
        };
        println!(
            "  {:<10} {:<8} used by {:<24}{}",
            db.id(),
            st.symbol(),
            db.used_by().join(","),
            detail
        );
    }

    // ---- 5. Per-validator readiness summary ----
    println!("\n[ validator readiness ]");
    let cv = root.status(Database::CheckV).ok();
    let cm = root.status(Database::CheckM2).ok();
    let pf = root.status(Database::Pfam).ok();
    let bu = root.status(Database::Busco).ok();
    line("CheckV (viral, reference completeness)", cv, "terminal-repeat + structural signals still run without it");
    line("CheckM2 (bacteria/archaea prediction)", cm, "needs the DIAMOND DB + models for full prediction");
    line("ViralVerify (virus/plasmid/chromosome)", pf, "needs Pfam-A to annotate contigs");
    line("BUSCO (eukaryote orthologs)", bu, "needs a lineage dataset");
    println!("  Novel modules (TRACE/CHIMERA/SAGE)       READY   (always database-free)");

    println!(
        "\nSummary: {ready}/{total} reference databases installed. \
         Download everything with `sentinel download all`, or one at a time \
         with `sentinel download <{}>`.",
        Database::all()
            .iter()
            .map(|d| d.id())
            .collect::<Vec<_>>()
            .join("|")
    );
    Ok(())
}

fn line(label: &str, ready: bool, note: &str) {
    println!(
        "  {:<40} {}   ({})",
        label,
        if ready { "READY" } else { "needs DB" },
        note
    );
}

/// `sentinel download <db|all> [--db DIR] [--lineage NAME]`
pub fn download_cmd(args: &[String]) -> Result<()> {
    let ids: Vec<&str> = Database::all().iter().map(|d| d.id()).collect();
    if args.is_empty() {
        bail!(
            "usage: sentinel download <all|{}> [--db DIR] [--lineage NAME]",
            ids.join("|")
        );
    }

    // Parse the target (first positional) and options.
    let target = args[0].to_ascii_lowercase();
    let mut db_path: Option<String> = None;
    let mut lineage: Option<String> = None;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--db" => {
                i += 1;
                db_path = args.get(i).cloned();
            }
            "--lineage" => {
                i += 1;
                lineage = args.get(i).cloned();
            }
            other => bail!("unknown download option '{other}'"),
        }
        i += 1;
    }

    let root = DbRoot::resolve(db_path.as_deref())?;
    println!("{}", sentinel_core::banner());

    // ---- download all ----
    if target == "all" {
        return download_all(&root, lineage.as_deref());
    }

    // ---- download one ----
    let db = Database::parse(&target)
        .ok_or_else(|| anyhow::anyhow!("unknown database '{}' (expected: all|{})", target, ids.join("|")))?;
    println!(
        "Downloading {} into {}\n",
        db.display(),
        root.dir(db).display()
    );
    root.install(db, lineage.as_deref())?;
    let st = root.status(db);
    println!("\nResult: {} is now {}", db.id(), st.symbol());
    if !st.ok() {
        println!("(some required files are still missing — see `sentinel doctor`)");
    }
    Ok(())
}

/// Download every managed database into `root`, continuing past failures and
/// printing a summary. BUSCO is included only when a lineage is supplied
/// (it cannot be fetched without one).
fn download_all(root: &DbRoot, lineage: Option<&str>) -> Result<()> {
    root.build_skeleton()?;
    println!("Downloading all reference databases into {}\n", root.root.display());

    let mut results: Vec<(Database, String)> = Vec::new();
    let mut had_error = false;

    for &db in Database::all() {
        // BUSCO needs a lineage; skip it in `all` unless one was given.
        if matches!(db, Database::Busco) && lineage.is_none() {
            println!(
                "[{}] skipped — BUSCO needs a lineage; run `sentinel download busco --lineage <name>`\n",
                db.id()
            );
            results.push((db, "SKIPPED (needs --lineage)".to_string()));
            continue;
        }

        // Skip ones already fully installed.
        if root.status(db).ok() {
            println!("[{}] already installed — skipping\n", db.id());
            results.push((db, "OK (already present)".to_string()));
            continue;
        }

        println!("==> {} ({})", db.display(), db.id());
        match root.install(db, lineage) {
            Ok(()) => {
                let st = root.status(db);
                if !st.ok() {
                    had_error = true;
                }
                results.push((db, st.symbol().to_string()));
                println!("    {} -> {}\n", db.id(), st.symbol());
            }
            Err(e) => {
                had_error = true;
                results.push((db, format!("FAILED: {e}")));
                eprintln!("    {} download failed: {e:#}\n", db.id());
            }
        }
    }

    // Summary.
    println!("\n==================== download all: summary ====================");
    for (db, status) in &results {
        println!("  {:<10} {}", db.id(), status);
    }
    let ok = results.iter().filter(|(db, _)| root.status(*db).ok()).count();
    println!("--------------------------------------------------------------");
    println!("{ok}/{} databases installed. Run `sentinel doctor` for details.", results.len());

    if had_error {
        bail!("one or more databases did not install completely (see summary above)");
    }
    Ok(())
}

/// Minimal `which`: is an executable on PATH?
fn which(tool: &str) -> bool {
    let Some(paths) = std::env::var_os("PATH") else {
        return false;
    };
    std::env::split_paths(&paths).any(|dir| {
        let p = dir.join(tool);
        p.is_file()
    })
}
