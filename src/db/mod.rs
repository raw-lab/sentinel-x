//! # Reference-database management
//!
//! The heavy reference databases are fetched from the same upstream sources
//! the original tools use, verified, and installed atomically under one root
//! (`--db`, else `$SENTINEL_DB`, else `~/.sentinel/db`):
//!
//! | database | source | verification | layout |
//! |---|---|---|---|
//! | CheckV  | `portal.nersc.gov/CheckV/CURRENT_RELEASE.txt` → `<release>.tar.gz` | none upstream; SHA-256 + size recorded | `checkv/<release>/{genome_db,hmm_db}` |
//! | CheckM2 | DOI `10.5281/zenodo.14897628` → Zenodo API | MD5 per file, `CONTENTS.json` MD5s, SHA-256 of the `.dmnd` | `checkm2/CheckM2_database/uniref100.KO.1.dmnd` |
//! | Pfam-A  | EBI `current_release/Pfam-A.hmm.gz` | `md5_checksums` | `pfam/Pfam-A.hmm` |
//! | BUSCO   | `busco-data.ezlab.org/v5/data/file_versions.tsv` → `lineages/<name>.<date>.tar.gz` | MD5 from the manifest | `busco/lineages/<name>/` |
//!
//! Downloads land in `<root>/.staging/`, are verified and unpacked there, and
//! are moved into place only when complete; a `.sentinel-complete.json`
//! marker (version, checksums) records a verified install. The DIAMOND
//! databases are read in place by Sentinel's `protsearch` (no conversion, no
//! DIAMOND binary).

pub mod dmnd;
pub mod download;

use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result};

/// CheckM2's DIAMOND database v3 (CheckM2 ≥ 1.1), from
/// `version_hashes_1.1.0.json`.
pub const CHECKM2_DOI: &str = "10.5281/zenodo.14897628";
pub const CHECKM2_DMND_SHA256: &str = "1b86ef3eac0813c1853f53182c17657045e3763d66f384ec95747261a63ae46f";
const MARKER: &str = ".sentinel-complete.json";

/// The CheckM2 database hash to require; `SENTINEL_CHECKM2_SHA256`
/// overrides it (mirrors of other database releases, tests).
pub fn expected_checkm2_sha256() -> String {
    match std::env::var("SENTINEL_CHECKM2_SHA256") {
        Ok(v) if !v.trim().is_empty() => v.trim().to_ascii_lowercase(),
        _ => CHECKM2_DMND_SHA256.to_string(),
    }
}

/// A managed reference database.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Database {
    /// CheckV genome + HMM database (viral completeness/contamination).
    CheckV,
    /// CheckM2 UniRef100/KO DIAMOND database (KO annotation).
    CheckM2,
    /// Pfam-A HMM library (viralVerify).
    Pfam,
    /// BUSCO lineage datasets.
    Busco,
}

impl Database {
    pub fn id(&self) -> &'static str {
        match self {
            Database::CheckV => "checkv",
            Database::CheckM2 => "checkm2",
            Database::Pfam => "pfam",
            Database::Busco => "busco",
        }
    }

    pub fn display(&self) -> &'static str {
        match self {
            Database::CheckV => "CheckV reference database",
            Database::CheckM2 => "CheckM2 UniRef100/KO DIAMOND database",
            Database::Pfam => "Pfam-A HMM library",
            Database::Busco => "BUSCO lineage datasets",
        }
    }

    pub fn used_by(&self) -> &'static [&'static str] {
        match self {
            Database::CheckV => &["checkv"],
            Database::CheckM2 => &["checkm2"],
            Database::Pfam => &["viralverify"],
            Database::Busco => &["busco"],
        }
    }

    pub fn subdir(&self) -> &'static str {
        self.id()
    }

    pub fn all() -> &'static [Database] {
        &[Database::CheckV, Database::CheckM2, Database::Pfam, Database::Busco]
    }

    pub fn parse(s: &str) -> Option<Database> {
        match s.to_ascii_lowercase().as_str() {
            "checkv" => Some(Database::CheckV),
            "checkm2" => Some(Database::CheckM2),
            "pfam" => Some(Database::Pfam),
            "busco" => Some(Database::Busco),
            _ => None,
        }
    }
}

/// Upstream base URLs (overridable for mirrors and tests).
pub struct Sources;

impl Sources {
    fn env_or(key: &str, default: &str) -> String {
        match std::env::var(key) {
            Ok(v) if !v.trim().is_empty() => {
                let v = v.trim().to_string();
                if v.ends_with('/') {
                    v
                } else {
                    format!("{v}/")
                }
            }
            _ => default.to_string(),
        }
    }
    pub fn checkv() -> String {
        Self::env_or("SENTINEL_CHECKV_URL", "https://portal.nersc.gov/CheckV/")
    }
    pub fn doi_resolver() -> String {
        Self::env_or("SENTINEL_DOI_RESOLVER", "https://doi.org/")
    }
    pub fn zenodo_api() -> String {
        Self::env_or("SENTINEL_ZENODO_API", "https://zenodo.org/api/records/")
    }
    pub fn pfam() -> String {
        Self::env_or("SENTINEL_PFAM_URL", "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/")
    }
    pub fn busco() -> Vec<String> {
        match std::env::var("SENTINEL_BUSCO_URL") {
            Ok(v) if !v.trim().is_empty() => vec![Self::env_or("SENTINEL_BUSCO_URL", "")],
            _ => vec!["https://busco-data.ezlab.org/v5/data/".into(), "https://busco-data2.ezlab.org/v5/data/".into()],
        }
    }
}

/// Installation status of one database.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Status {
    /// Present and verified by a Sentinel install marker.
    Installed { detail: String },
    /// All required files are present but no install marker (placed by hand).
    Unverified { detail: String },
    /// Some required files are missing.
    Partial { missing: Vec<String> },
    /// The folder exists but holds nothing usable.
    Empty,
    /// The folder does not exist.
    NotInstalled,
}

impl Status {
    /// Usable by the validators.
    pub fn ok(&self) -> bool {
        matches!(self, Status::Installed { .. } | Status::Unverified { .. })
    }
    pub fn symbol(&self) -> &'static str {
        match self {
            Status::Installed { .. } => "OK",
            Status::Unverified { .. } => "OK*",
            Status::Partial { .. } => "PARTIAL",
            Status::Empty => "EMPTY",
            Status::NotInstalled => "MISSING",
        }
    }
}

/// The Sentinel database root directory.
#[derive(Clone, Debug)]
pub struct DbRoot {
    pub root: PathBuf,
}

/// Subdirectories of `dir`. Symbolic links to directories are followed only
/// when `follow` is set (used for the direct children of a database folder,
/// where a user may link an existing installation); walks never follow
/// links deeper down, and archives never create links.
fn subdirs_of(dir: &Path, follow: bool) -> Vec<PathBuf> {
    let mut v = Vec::new();
    if let Ok(rd) = std::fs::read_dir(dir) {
        for e in rd.flatten() {
            let is_dir = match e.file_type() {
                Ok(t) if t.is_dir() => true,
                Ok(t) if t.is_symlink() && follow => e.path().is_dir(),
                _ => false,
            };
            if is_dir {
                v.push(e.path());
            }
        }
    }
    v.sort();
    v
}

fn subdirs(dir: &Path) -> Vec<PathBuf> {
    subdirs_of(dir, true)
}

/// Is `p` a regular file (not following a final symlink to elsewhere is not
/// required here: users may link their own copies)?
fn is_file(p: &Path) -> bool {
    p.is_file()
}

/// Find a file by name below `base` (depth-limited, symlinked directories
/// not followed).
pub fn find_file(base: &Path, name: &str, max_depth: usize) -> Option<PathBuf> {
    let mut stack = vec![(base.to_path_buf(), 0usize)];
    while let Some((d, depth)) = stack.pop() {
        let direct = d.join(name);
        if is_file(&direct) {
            return Some(direct);
        }
        if depth < max_depth {
            for s in subdirs_of(&d, depth == 0).into_iter().rev() {
                if !s.file_name().map(|n| n.to_string_lossy().starts_with('.')).unwrap_or(false) {
                    stack.push((s, depth + 1));
                }
            }
        }
    }
    None
}

fn read_marker(dir: &Path) -> Option<serde_json::Value> {
    let t = std::fs::read_to_string(dir.join(MARKER)).ok()?;
    serde_json::from_str(&t).ok()
}

fn write_marker(dir: &Path, v: &serde_json::Value) -> Result<()> {
    let p = dir.join(MARKER);
    std::fs::write(&p, serde_json::to_string_pretty(v)?).with_context(|| format!("writing {}", p.display()))
}

fn now() -> u64 {
    std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0)
}

/// Move `staged` to `dest`, replacing an existing `dest` only after the
/// staged copy is in place.
fn swap_into_place(staged: &Path, dest: &Path) -> Result<()> {
    if let Some(p) = dest.parent() {
        std::fs::create_dir_all(p)?;
    }
    let old = dest.with_extension(format!("old-{}", std::process::id()));
    let had_old = dest.exists();
    if had_old {
        std::fs::rename(dest, &old).with_context(|| format!("moving aside {}", dest.display()))?;
    }
    if let Err(e) = std::fs::rename(staged, dest) {
        if had_old {
            let _ = std::fs::rename(&old, dest);
        }
        return Err(e).with_context(|| format!("installing {}", dest.display()));
    }
    if had_old {
        if old.is_dir() {
            let _ = std::fs::remove_dir_all(&old);
        } else {
            let _ = std::fs::remove_file(&old);
        }
    }
    Ok(())
}

/// A release / dataset name safe to use as a path component.
fn safe_name(s: &str) -> Result<String> {
    let t = s.trim();
    if t.is_empty() || t == "." || t == ".." || !t.chars().all(|c| c.is_ascii_alphanumeric() || "._-".contains(c)) {
        bail!("unexpected release name '{t}'");
    }
    Ok(t.to_string())
}

impl DbRoot {
    /// `~/.sentinel/db`.
    pub fn default_root() -> PathBuf {
        match std::env::var_os("HOME") {
            Some(h) if !h.is_empty() => PathBuf::from(h).join(".sentinel").join("db"),
            _ => PathBuf::from(".sentinel").join("db"),
        }
    }

    /// Resolve the DB root: explicit `path`, else `$SENTINEL_DB`, else
    /// `~/.sentinel/db`.
    pub fn resolve(path: Option<&str>) -> DbRoot {
        if let Some(p) = path {
            return DbRoot { root: PathBuf::from(p) };
        }
        if let Ok(env) = std::env::var("SENTINEL_DB") {
            if !env.trim().is_empty() {
                return DbRoot { root: PathBuf::from(env) };
            }
        }
        DbRoot { root: Self::default_root() }
    }

    pub fn dir(&self, db: Database) -> PathBuf {
        self.root.join(db.subdir())
    }

    fn staging(&self, db: Database) -> PathBuf {
        self.root.join(".staging").join(db.id())
    }

    /// Create the folder skeleton and manifest.
    pub fn build_skeleton(&self) -> Result<()> {
        std::fs::create_dir_all(&self.root).with_context(|| format!("creating DB root {}", self.root.display()))?;
        for db in Database::all() {
            std::fs::create_dir_all(self.dir(*db)).with_context(|| format!("creating {} dir", db.id()))?;
        }
        self.write_manifest()
    }

    // ---------------------------------------------------------------- lookup

    /// The CheckV release directory to use: the release recorded by the
    /// install marker if valid, else the newest valid release found.
    pub fn checkv_dir(&self) -> Option<PathBuf> {
        let dir = self.dir(Database::CheckV);
        if let Some(m) = read_marker(&dir) {
            if let Some(rel) = m.get("release").and_then(|r| r.as_str()) {
                let p = dir.join(rel);
                if crate::checkv::db::missing_files(&p).is_empty() {
                    return Some(p);
                }
            }
        }
        if crate::checkv::db::missing_files(&dir).is_empty() {
            return Some(dir);
        }
        crate::checkv::db::find_checkv_dirs(&dir).into_iter().find(|p| crate::checkv::db::missing_files(p).is_empty())
    }

    /// CheckM2's DIAMOND database (`uniref100.KO.1.dmnd`, or a FASTA copy).
    pub fn checkm2_db(&self) -> Option<PathBuf> {
        let dir = self.dir(Database::CheckM2);
        for name in ["uniref100.KO.1.dmnd", "uniref100.KO.1.fasta", "uniref100.KO.1.faa"] {
            if let Some(p) = find_file(&dir, name, 3) {
                return Some(p);
            }
        }
        None
    }

    /// Pfam-A HMM library (plain or gzip).
    pub fn pfam_hmm(&self) -> Option<PathBuf> {
        let dir = self.dir(Database::Pfam);
        find_file(&dir, "Pfam-A.hmm", 2).or_else(|| find_file(&dir, "Pfam-A.hmm.gz", 2))
    }

    /// Installed BUSCO lineages (`busco/lineages/<name>/dataset.cfg`, or
    /// directly under `busco/`).
    pub fn busco_lineages(&self) -> Vec<(String, PathBuf)> {
        let base = self.dir(Database::Busco);
        let mut out = Vec::new();
        for d in [base.join("lineages"), base.clone()] {
            for s in subdirs(&d) {
                if s.join("dataset.cfg").is_file() {
                    let name = s.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
                    if !out.iter().any(|(n, _)| *n == name) {
                        out.push((name, s));
                    }
                }
            }
        }
        out
    }

    // ---------------------------------------------------------------- status

    pub fn status(&self, db: Database) -> Status {
        let dir = self.dir(db);
        if !dir.is_dir() {
            return Status::NotInstalled;
        }
        let marker = read_marker(&dir);
        match db {
            Database::CheckV => match self.checkv_dir() {
                Some(p) => {
                    let rel = p.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
                    let verified = marker
                        .as_ref()
                        .and_then(|m| m.get("release").and_then(|r| r.as_str()).map(|r| r == rel))
                        .unwrap_or(false);
                    if verified {
                        Status::Installed { detail: rel }
                    } else {
                        Status::Unverified { detail: rel }
                    }
                }
                None => {
                    let cands = crate::checkv::db::find_checkv_dirs(&dir);
                    match cands.first() {
                        Some(c) => Status::Partial { missing: crate::checkv::db::missing_files(c) },
                        None => {
                            if subdirs(&dir).is_empty() {
                                Status::Empty
                            } else {
                                Status::Partial { missing: crate::checkv::db::missing_files(&dir) }
                            }
                        }
                    }
                }
            },
            Database::CheckM2 => match self.checkm2_db() {
                Some(p) => {
                    let name = p.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
                    let verified = marker
                        .as_ref()
                        .and_then(|m| m.get("sha256").and_then(|s| s.as_str()).map(|s| s == expected_checkm2_sha256()))
                        .unwrap_or(false);
                    if verified {
                        Status::Installed { detail: format!("{name} (DIAMOND DB v3)") }
                    } else {
                        Status::Unverified { detail: name }
                    }
                }
                None => Status::Empty,
            },
            Database::Pfam => match self.pfam_hmm() {
                Some(p) => {
                    let name = p.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
                    match marker.as_ref().and_then(|m| m.get("release").and_then(|r| r.as_str()).map(String::from)) {
                        Some(rel) => Status::Installed { detail: format!("{name} (Pfam {rel})") },
                        None => Status::Unverified { detail: name },
                    }
                }
                None => Status::Empty,
            },
            Database::Busco => {
                let l = self.busco_lineages();
                if l.is_empty() {
                    return Status::Empty;
                }
                let names: Vec<String> = l.iter().map(|(n, _)| n.clone()).collect();
                if l.iter().all(|(_, p)| read_marker(p).is_some()) {
                    Status::Installed { detail: names.join(", ") }
                } else {
                    Status::Unverified { detail: names.join(", ") }
                }
            }
        }
    }

    pub fn status_all(&self) -> Vec<(Database, Status)> {
        Database::all().iter().map(|db| (*db, self.status(*db))).collect()
    }

    // ---------------------------------------------------------------- install

    /// Download, verify and install a database (BUSCO needs a lineage).
    pub fn install(&self, db: Database, lineage: Option<&str>) -> Result<()> {
        std::fs::create_dir_all(self.dir(db)).with_context(|| format!("creating {}", self.dir(db).display()))?;
        let staging = self.staging(db);
        std::fs::create_dir_all(&staging)?;
        let r = match db {
            Database::CheckV => self.install_checkv(&staging),
            Database::CheckM2 => self.install_checkm2(&staging),
            Database::Pfam => self.install_pfam(&staging),
            Database::Busco => {
                let l = lineage.ok_or_else(|| anyhow!("BUSCO needs a lineage name, e.g. --lineage bacteria_odb12.2"))?;
                self.install_busco(&staging, l)
            }
        };
        // Clean up after success. After a failure keep the (verified or
        // partial) downloads so a retry can resume, but drop what can be
        // regenerated from them: the extraction tree and a decompressed
        // Pfam-A.hmm (gigabytes each for CheckV, CheckM2 and Pfam).
        if r.is_ok() {
            let _ = std::fs::remove_dir_all(&staging);
            let _ = std::fs::remove_dir(self.root.join(".staging")); // only if empty
        } else {
            let _ = std::fs::remove_dir_all(staging.join("extract"));
            let _ = std::fs::remove_file(staging.join("Pfam-A.hmm"));
        }
        let _ = self.write_manifest();
        r
    }

    fn install_checkv(&self, staging: &Path) -> Result<()> {
        let base = Sources::checkv();
        let rel = safe_name(
            download::fetch_text(&format!("{base}CURRENT_RELEASE.txt"))
                .context("discovering the current CheckV release")?
                .lines()
                .next()
                .unwrap_or(""),
        )?;
        let archive = staging.join(format!("{rel}.tar.gz"));
        if !archive.is_file() {
            eprintln!("[sentinel-db] downloading CheckV {rel}");
            download::fetch(&format!("{base}{rel}.tar.gz"), &archive)?;
        }
        let size = std::fs::metadata(&archive)?.len();
        eprintln!("[sentinel-db] hashing {} ({} MB)", archive.display(), size >> 20);
        let sha = download::sha256_hex(&archive)?;
        let extract = staging.join("extract");
        let _ = std::fs::remove_dir_all(&extract);
        eprintln!("[sentinel-db] extracting");
        if let Err(e) = download::extract_tar_gz(&archive, &extract) {
            let _ = std::fs::remove_file(&archive); // corrupt: force a fresh download next time
            return Err(e);
        }
        let root = if crate::checkv::db::missing_files(&extract.join(&rel)).is_empty() {
            extract.join(&rel)
        } else {
            crate::checkv::db::find_checkv_dirs(&extract)
                .into_iter()
                .find(|p| crate::checkv::db::missing_files(p).is_empty())
                .ok_or_else(|| anyhow!("the CheckV archive does not contain a complete database"))?
        };
        let dest = self.dir(Database::CheckV).join(&rel);
        swap_into_place(&root, &dest)?;
        write_marker(
            &self.dir(Database::CheckV),
            &serde_json::json!({
                "database": "checkv", "release": rel, "archive": format!("{rel}.tar.gz"),
                "sha256": sha, "size": size, "source": base, "installed": now(),
            }),
        )?;
        eprintln!("[sentinel-db] CheckV {rel} installed in {}", dest.display());
        Ok(())
    }

    fn install_checkm2(&self, staging: &Path) -> Result<()> {
        let resolved = download::resolve_redirect(&format!("{}{}", Sources::doi_resolver(), CHECKM2_DOI))
            .with_context(|| format!("resolving DOI {CHECKM2_DOI}"))?;
        let id = download::zenodo_record_id(&resolved)?;
        let api = format!("{}{}", Sources::zenodo_api(), id);
        let rec = download::parse_zenodo_record(&id, &download::fetch_text(&api).with_context(|| format!("querying {api}"))?)?;
        if rec.files.is_empty() {
            bail!("Zenodo record {id} lists no files");
        }
        for f in &rec.files {
            let dest = staging.join(&f.name);
            for attempt in 0..2 {
                if !dest.is_file() {
                    eprintln!("[sentinel-db] downloading {} (Zenodo record {id})", f.name);
                    download::fetch(&f.url, &dest)?;
                }
                match download::verify_checksum(&dest, &f.checksum) {
                    Ok(()) => break,
                    Err(e) => {
                        let _ = std::fs::remove_file(&dest);
                        if attempt == 1 {
                            return Err(e);
                        }
                        eprintln!("[sentinel-db] {e}; downloading again");
                    }
                }
            }
        }
        let extract = staging.join("extract");
        let _ = std::fs::remove_dir_all(&extract);
        for f in rec.files.iter().filter(|f| f.name.ends_with(".tar.gz")) {
            eprintln!("[sentinel-db] extracting {}", f.name);
            download::extract_tar_gz(&staging.join(&f.name), &extract)?;
        }
        verify_backpack(&extract, &rec.version)?;
        let dmnd = find_file(&extract, "uniref100.KO.1.dmnd", 3)
            .ok_or_else(|| anyhow!("uniref100.KO.1.dmnd not found in the CheckM2 archive"))?;
        eprintln!("[sentinel-db] verifying SHA-256 of {}", dmnd.display());
        let sha = download::sha256_hex(&dmnd)?;
        let expected = expected_checkm2_sha256();
        if sha != expected {
            bail!("uniref100.KO.1.dmnd SHA-256 {sha} does not match CheckM2's database v3 ({expected})");
        }
        let payload = dmnd.parent().map(Path::to_path_buf).unwrap_or(extract.clone());
        let dest = self.dir(Database::CheckM2).join("CheckM2_database");
        swap_into_place(&payload, &dest)?;
        write_marker(
            &self.dir(Database::CheckM2),
            &serde_json::json!({
                "database": "checkm2", "doi": CHECKM2_DOI, "zenodo_record": id, "version": rec.version,
                "sha256": sha, "installed": now(),
            }),
        )?;
        eprintln!("[sentinel-db] CheckM2 database installed in {}", dest.display());
        Ok(())
    }

    fn install_pfam(&self, staging: &Path) -> Result<()> {
        let base = Sources::pfam();
        // release number and published MD5 (both optional: report if absent)
        let release = {
            let gz = staging.join("Pfam.version.gz");
            let _ = std::fs::remove_file(&gz);
            match download::fetch(&format!("{base}Pfam.version.gz"), &gz) {
                Ok(()) => {
                    let txt = String::from_utf8_lossy(&crate::common::seqio::read_all(&gz)?).into_owned();
                    txt.lines()
                        .find_map(|l| l.strip_prefix("Pfam release").map(|r| r.trim_start_matches([' ', ':']).trim().to_string()))
                        .unwrap_or_else(|| "unknown".into())
                }
                Err(e) => {
                    eprintln!("[sentinel-db] (Pfam.version.gz unavailable: {e})");
                    "unknown".into()
                }
            }
        };
        let md5 = download::fetch_text(&format!("{base}md5_checksums")).ok().and_then(|t| {
            t.lines().find_map(|l| {
                let f: Vec<&str> = l.split_whitespace().collect();
                (f.len() >= 2 && f[1] == "Pfam-A.hmm.gz").then(|| f[0].to_string())
            })
        });
        let gz = staging.join("Pfam-A.hmm.gz");
        for attempt in 0..2 {
            if !gz.is_file() {
                eprintln!("[sentinel-db] downloading Pfam-A.hmm.gz (release {release})");
                download::fetch(&format!("{base}Pfam-A.hmm.gz"), &gz)?;
            }
            match &md5 {
                Some(m) => match download::verify_checksum(&gz, &format!("md5:{m}")) {
                    Ok(()) => break,
                    Err(e) => {
                        // a resumed transfer can splice two releases together
                        let _ = std::fs::remove_file(&gz);
                        if attempt == 1 {
                            return Err(e);
                        }
                        eprintln!("[sentinel-db] {e}; downloading again");
                    }
                },
                None => {
                    eprintln!("[sentinel-db] (no published MD5 found; integrity checked by decompression only)");
                    break;
                }
            }
        }
        let hmm = staging.join("Pfam-A.hmm");
        eprintln!("[sentinel-db] decompressing");
        download::gunzip(&gz, &hmm)?;
        let mut head = [0u8; 6];
        use std::io::Read;
        std::fs::File::open(&hmm)?.read_exact(&mut head).context("reading Pfam-A.hmm")?;
        if &head != b"HMMER3" {
            bail!("Pfam-A.hmm does not start with an HMMER3 header");
        }
        swap_into_place(&hmm, &self.dir(Database::Pfam).join("Pfam-A.hmm"))?;
        write_marker(
            &self.dir(Database::Pfam),
            &serde_json::json!({"database": "pfam", "release": release, "md5": md5, "source": base, "installed": now()}),
        )?;
        eprintln!("[sentinel-db] Pfam-A {release} installed");
        Ok(())
    }

    fn install_busco(&self, staging: &Path, lineage: &str) -> Result<()> {
        let lineage = safe_name(lineage)?;
        let mut last_err = None;
        for base in Sources::busco() {
            match self.install_busco_from(staging, &lineage, &base) {
                Ok(()) => return Ok(()),
                Err(e) => {
                    eprintln!("[sentinel-db] BUSCO mirror {base} failed: {e:#}");
                    last_err = Some(e);
                }
            }
        }
        Err(last_err.unwrap_or_else(|| anyhow!("no BUSCO mirror configured")))
    }

    fn install_busco_from(&self, staging: &Path, lineage: &str, base: &str) -> Result<()> {
        let manifest = download::fetch_text(&format!("{base}file_versions.tsv"))?;
        let rows: Vec<Vec<&str>> = manifest.lines().map(|l| l.split('\t').collect::<Vec<_>>()).filter(|f| f.len() >= 5).collect();
        let row = rows.iter().find(|f| f[0] == lineage && f[4] == "lineages").ok_or_else(|| {
            let similar: Vec<&str> = rows
                .iter()
                .filter(|f| f[4] == "lineages" && f[0].starts_with(lineage.split("_odb").next().unwrap_or(lineage)))
                .map(|f| f[0])
                .collect();
            anyhow!(
                "BUSCO lineage '{lineage}' is not in file_versions.tsv{}",
                if similar.is_empty() { String::new() } else { format!(" (available: {})", similar.join(", ")) }
            )
        })?;
        let (date, md5) = (safe_name(row[1])?, row[2].to_string());
        let file = format!("{lineage}.{date}.tar.gz");
        let archive = staging.join(&file);
        for attempt in 0..3 {
            if !archive.is_file() {
                eprintln!("[sentinel-db] downloading BUSCO {file}");
                download::fetch(&format!("{base}lineages/{file}"), &archive)?;
            }
            match download::verify_checksum(&archive, &format!("md5:{md5}")) {
                Ok(()) => break,
                Err(e) => {
                    let _ = std::fs::remove_file(&archive);
                    if attempt == 2 {
                        return Err(e);
                    }
                    eprintln!("[sentinel-db] {e}; downloading again");
                }
            }
        }
        let extract = staging.join("extract");
        let _ = std::fs::remove_dir_all(&extract);
        download::extract_tar_gz(&archive, &extract)?;
        let src = if extract.join(lineage).join("dataset.cfg").is_file() {
            extract.join(lineage)
        } else {
            subdirs(&extract)
                .into_iter()
                .find(|d| d.join("dataset.cfg").is_file())
                .ok_or_else(|| anyhow!("{file} does not contain a BUSCO dataset (no dataset.cfg)"))?
        };
        write_marker(
            &src,
            &serde_json::json!({"database": "busco", "dataset": lineage, "date": date, "md5": md5, "source": base, "installed": now()}),
        )?;
        let dest = self.dir(Database::Busco).join("lineages").join(lineage);
        swap_into_place(&src, &dest)?;
        eprintln!("[sentinel-db] BUSCO {lineage} ({date}) installed in {}", dest.display());
        Ok(())
    }

    /// Write `manifest.json` capturing the current status of each database.
    pub fn write_manifest(&self) -> Result<()> {
        let entries: Vec<serde_json::Value> = Database::all()
            .iter()
            .map(|db| {
                let st = self.status(*db);
                let detail = match &st {
                    Status::Installed { detail } | Status::Unverified { detail } => detail.clone(),
                    Status::Partial { missing } => format!("missing: {}", missing.join(", ")),
                    _ => String::new(),
                };
                serde_json::json!({
                    "id": db.id(), "name": db.display(), "subdir": db.subdir(),
                    "used_by": db.used_by(), "status": st.symbol(), "detail": detail,
                })
            })
            .collect();
        let manifest = serde_json::json!({
            "sentinel_db_version": 2,
            "root": self.root.to_string_lossy(),
            "databases": entries,
        });
        std::fs::create_dir_all(&self.root)?;
        let path = self.root.join("manifest.json");
        std::fs::write(&path, serde_json::to_string_pretty(&manifest)?).with_context(|| format!("writing {}", path.display()))
    }
}

/// Verify an extracted Zenodo backpack (`CONTENTS.json`): backpack format 1,
/// the record's version, and the MD5 of every payload file.
pub fn verify_backpack(dir: &Path, record_version: &str) -> Result<()> {
    let p = dir.join("CONTENTS.json");
    let txt = std::fs::read_to_string(&p).context("the CheckM2 archive has no CONTENTS.json")?;
    let mut v: serde_json::Map<String, serde_json::Value> = serde_json::from_str(&txt).context("parsing CONTENTS.json")?;
    let version = v.remove("version").map(|x| x.as_str().map(String::from).unwrap_or_else(|| x.to_string())).unwrap_or_default();
    let bp = v.remove("zenodo_backpack_version").and_then(|x| x.as_i64()).unwrap_or(-1);
    if bp != 1 {
        bail!("CONTENTS.json: unsupported zenodo_backpack_version {bp}");
    }
    if !record_version.is_empty() && version.trim() != record_version.trim() {
        bail!("CONTENTS.json version '{version}' does not match the Zenodo record version '{record_version}'");
    }
    for (file, md5) in v {
        let rel = Path::new(&file);
        if rel.is_absolute() || rel.components().any(|c| matches!(c, std::path::Component::ParentDir)) {
            bail!("CONTENTS.json: unsafe path '{file}'");
        }
        let want = md5.as_str().ok_or_else(|| anyhow!("CONTENTS.json: bad checksum for {file}"))?;
        eprintln!("[sentinel-db] verifying {file}");
        download::verify_checksum(&dir.join(rel), &format!("md5:{want}"))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("sentdb_{tag}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn ids_and_parsing_roundtrip() {
        for db in Database::all() {
            assert_eq!(Database::parse(db.id()), Some(*db));
        }
        assert_eq!(Database::parse("CheckV"), Some(Database::CheckV));
        assert_eq!(Database::parse("kofam"), None);
    }

    #[test]
    fn resolve_order() {
        assert_eq!(DbRoot::resolve(Some("/x/db")).root, PathBuf::from("/x/db"));
        assert!(DbRoot::default_root().ends_with(".sentinel/db"));
    }

    fn fake_checkv(dir: &Path) {
        std::fs::create_dir_all(dir.join("genome_db")).unwrap();
        std::fs::create_dir_all(dir.join("hmm_db/checkv_hmms")).unwrap();
        for f in [
            "genome_db/checkv_reps.tsv",
            "genome_db/checkv_reps.faa",
            "genome_db/checkv_error.tsv",
            "hmm_db/checkv_hmms.tsv",
            "hmm_db/genome_lengths.tsv",
            "hmm_db/checkv_hmms/1.hmm",
        ] {
            std::fs::write(dir.join(f), b"x").unwrap();
        }
    }

    #[test]
    fn checkv_layout_status() {
        let t = tmp("cv");
        let root = DbRoot { root: t.clone() };
        root.build_skeleton().unwrap();
        assert_eq!(root.status(Database::CheckV), Status::Empty);
        let rel = root.dir(Database::CheckV).join("checkv-db-v1.5");
        fake_checkv(&rel);
        assert!(matches!(root.status(Database::CheckV), Status::Unverified { .. }));
        write_marker(&root.dir(Database::CheckV), &serde_json::json!({"release": "checkv-db-v1.5"})).unwrap();
        assert!(matches!(root.status(Database::CheckV), Status::Installed { .. }));
        assert_eq!(root.checkv_dir().unwrap(), rel);
        std::fs::remove_file(rel.join("hmm_db/checkv_hmms/1.hmm")).unwrap();
        match root.status(Database::CheckV) {
            Status::Partial { missing } => assert!(missing.iter().any(|m| m.contains("checkv_hmms")), "{missing:?}"),
            other => panic!("{other:?}"),
        }
        std::fs::remove_dir_all(&t).ok();
    }

    #[test]
    fn busco_and_backpack() {
        let t = tmp("bp");
        std::fs::write(t.join("a.txt"), b"abc").unwrap();
        std::fs::write(
            t.join("CONTENTS.json"),
            br#"{"version": "3", "zenodo_backpack_version": 1, "a.txt": "900150983cd24fb0d6963f7d28e17f72"}"#,
        )
        .unwrap();
        assert!(verify_backpack(&t, "3").is_ok());
        assert!(verify_backpack(&t, "2").is_err());
        let root = DbRoot { root: t.join("root") };
        let l = root.dir(Database::Busco).join("lineages").join("bacteria_odb12.2");
        std::fs::create_dir_all(&l).unwrap();
        std::fs::write(l.join("dataset.cfg"), b"name=bacteria_odb12.2\n").unwrap();
        assert_eq!(root.busco_lineages().len(), 1);
        assert!(root.status(Database::Busco).ok());
        std::fs::remove_dir_all(&t).ok();
    }

    #[test]
    fn symlinked_dirs_are_not_walked_below_the_top() {
        let t = tmp("sym");
        std::fs::create_dir_all(t.join("a")).unwrap();
        #[cfg(unix)]
        std::os::unix::fs::symlink("/etc", t.join("a/etc")).unwrap();
        assert!(find_file(&t, "passwd", 4).is_none());
        std::fs::remove_dir_all(&t).ok();
    }
}
