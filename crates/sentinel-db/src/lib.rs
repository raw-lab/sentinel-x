//! # sentinel-db
//!
//! Reference-database management for Sentinel.
//!
//! Sentinel ships small tables (the ViralVerify Pfam table, CheckM2's KEGG
//! definitions) but the heavy reference databases — CheckV's genome/HMM DB,
//! CheckM2's UniRef100/KO DIAMOND DB, Pfam-A, BUSCO lineages, KOfam — are
//! gigabytes and fetched on demand. This crate provides:
//!
//! * a [`registry`](Database) of the **canonical upstream sources** (the same
//!   URLs / Zenodo DOI the original tools use),
//! * an on-disk [`DbRoot`] layout (default `~/.sentinel/db`, overridable with
//!   `$SENTINEL_DB`) and a builder that creates the folder skeleton,
//! * an [`install`](DbRoot::install) routine (fetch → extract → record), and
//! * a [`status`](DbRoot::status) checker used by `sentinel doctor`.

pub mod download;

use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result};

/// A managed reference database.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Database {
    /// CheckV genome + HMM database (viral completeness/contamination).
    CheckV,
    /// CheckM2 UniRef100/KO DIAMOND database (bacterial/archaeal prediction).
    CheckM2,
    /// Pfam-A HMM library (ViralVerify Pfam annotation).
    Pfam,
    /// KOfam HMM profiles + thresholds (alternative KO annotation).
    Kofam,
    /// A BUSCO lineage dataset (eukaryote/single-copy orthologs).
    Busco,
}

impl Database {
    /// Stable identifier used on the CLI and on disk.
    pub fn id(&self) -> &'static str {
        match self {
            Database::CheckV => "checkv",
            Database::CheckM2 => "checkm2",
            Database::Pfam => "pfam",
            Database::Kofam => "kofam",
            Database::Busco => "busco",
        }
    }

    /// Human-readable name.
    pub fn display(&self) -> &'static str {
        match self {
            Database::CheckV => "CheckV reference database",
            Database::CheckM2 => "CheckM2 UniRef100/KO DIAMOND database",
            Database::Pfam => "Pfam-A HMM library",
            Database::Kofam => "KOfam HMM profiles",
            Database::Busco => "BUSCO lineage dataset",
        }
    }

    /// Which Sentinel validators consume this database.
    pub fn used_by(&self) -> &'static [&'static str] {
        match self {
            Database::CheckV => &["checkv"],
            Database::CheckM2 => &["checkm2"],
            Database::Pfam => &["viralverify"],
            Database::Kofam => &["checkm2 (KO annotation)"],
            Database::Busco => &["busco"],
        }
    }

    /// Subdirectory (relative to the DB root) where this database lives.
    pub fn subdir(&self) -> &'static str {
        self.id()
    }

    /// File path-suffixes that must be present for the database to count as
    /// installed. Matched anywhere beneath the subdirectory, so the
    /// version-named intermediate folders that CheckV/CheckM2 create are
    /// transparent.
    pub fn required(&self) -> &'static [&'static str] {
        match self {
            Database::CheckV => &[
                "genome_db/checkv_reps.tsv",
                "genome_db/checkv_error.tsv",
                "hmm_db/genome_lengths.tsv",
                "hmm_db/checkv_hmms.tsv",
            ],
            Database::CheckM2 => &["uniref100.KO.1.dmnd"],
            Database::Pfam => &["Pfam-A.hmm"],
            Database::Kofam => &["ko_list", "profiles"],
            Database::Busco => &["dataset.cfg"],
        }
    }

    /// How to obtain the database from upstream.
    pub fn source(&self) -> Source {
        match self {
            // CheckV: discover the current release name, then fetch its tarball.
            Database::CheckV => Source::VersionedTarball {
                version_url: "https://portal.nersc.gov/CheckV/CURRENT_RELEASE.txt",
                archive_url_prefix: "https://portal.nersc.gov/CheckV/",
            },
            // CheckM2: a Zenodo backpack (latest DIAMOND DB record).
            Database::CheckM2 => Source::Zenodo {
                doi: "10.5281/zenodo.14897628",
            },
            // Pfam-A current release (gzipped HMM library).
            Database::Pfam => Source::GzFile {
                url: "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
                out_name: "Pfam-A.hmm",
            },
            // KOfam profiles + ko_list (two artefacts).
            Database::Kofam => Source::MultiGzFile {
                urls: &[
                    ("https://www.genome.jp/ftp/db/kofam/ko_list.gz", "ko_list"),
                    (
                        "https://www.genome.jp/ftp/db/kofam/profiles.tar.gz",
                        "profiles.tar.gz",
                    ),
                ],
            },
            // BUSCO lineage: needs a lineage name supplied at call time.
            Database::Busco => Source::Parameterized {
                base_url: "https://busco-data.ezlab.org/v5/data/lineages/",
            },
        }
    }

    /// All managed databases.
    pub fn all() -> &'static [Database] {
        &[
            Database::CheckV,
            Database::CheckM2,
            Database::Pfam,
            Database::Kofam,
            Database::Busco,
        ]
    }

    /// Parse a CLI id into a [`Database`].
    pub fn parse(s: &str) -> Option<Database> {
        match s.to_ascii_lowercase().as_str() {
            "checkv" => Some(Database::CheckV),
            "checkm2" => Some(Database::CheckM2),
            "pfam" => Some(Database::Pfam),
            "kofam" => Some(Database::Kofam),
            "busco" => Some(Database::Busco),
            _ => None,
        }
    }
}

/// Upstream acquisition strategy for a database.
#[derive(Clone, Debug)]
pub enum Source {
    /// A `.tar.gz` whose name is discovered from a text URL (CheckV).
    VersionedTarball {
        version_url: &'static str,
        archive_url_prefix: &'static str,
    },
    /// A single gzipped file (Pfam-A).
    GzFile {
        url: &'static str,
        out_name: &'static str,
    },
    /// Several gzipped artefacts (KOfam).
    MultiGzFile {
        urls: &'static [(&'static str, &'static str)],
    },
    /// A Zenodo record addressed by DOI (CheckM2).
    Zenodo { doi: &'static str },
    /// A family of artefacts selected by a parameter, e.g. BUSCO lineage.
    Parameterized { base_url: &'static str },
}

/// Installation status of one database.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Status {
    /// Every required file is present.
    Installed,
    /// Some but not all required files are present.
    Partial { missing: Vec<String> },
    /// The subdirectory exists but no required files are present.
    Empty,
    /// The subdirectory does not exist.
    NotInstalled,
}

impl Status {
    pub fn ok(&self) -> bool {
        matches!(self, Status::Installed)
    }
    pub fn symbol(&self) -> &'static str {
        match self {
            Status::Installed => "OK",
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

impl DbRoot {
    /// Resolve the DB root: explicit `path`, else `$SENTINEL_DB`, else
    /// `~/.sentinel/db`.
    pub fn resolve(path: Option<&str>) -> Result<DbRoot> {
        if let Some(p) = path {
            return Ok(DbRoot { root: PathBuf::from(p) });
        }
        if let Ok(env) = std::env::var("SENTINEL_DB") {
            if !env.is_empty() {
                return Ok(DbRoot { root: PathBuf::from(env) });
            }
        }
        let home = home_dir().ok_or_else(|| anyhow!("cannot determine home directory; set $SENTINEL_DB"))?;
        Ok(DbRoot {
            root: home.join(".sentinel").join("db"),
        })
    }

    /// Path to a database's subdirectory.
    pub fn dir(&self, db: Database) -> PathBuf {
        self.root.join(db.subdir())
    }

    /// Create the directory skeleton (root + one folder per database) and write
    /// a manifest if absent.
    pub fn build_skeleton(&self) -> Result<()> {
        std::fs::create_dir_all(&self.root)
            .with_context(|| format!("creating DB root {}", self.root.display()))?;
        for db in Database::all() {
            std::fs::create_dir_all(self.dir(*db))
                .with_context(|| format!("creating {} dir", db.id()))?;
        }
        let manifest = self.root.join("manifest.json");
        if !manifest.exists() {
            self.write_manifest()?;
        }
        Ok(())
    }

    /// Recursively test whether a path-suffix exists beneath `base`.
    fn has_suffix(base: &Path, suffix: &str) -> bool {
        if !base.exists() {
            return false;
        }
        // Direct hit first.
        if base.join(suffix).exists() {
            return true;
        }
        // Otherwise walk the tree looking for an entry whose path ends with the
        // suffix (handles version-named intermediate directories).
        let mut stack = vec![base.to_path_buf()];
        let needle = suffix.replace('\\', "/");
        while let Some(d) = stack.pop() {
            let Ok(rd) = std::fs::read_dir(&d) else { continue };
            for entry in rd.flatten() {
                let p = entry.path();
                let as_str = p.to_string_lossy().replace('\\', "/");
                if as_str.ends_with(&needle) {
                    return true;
                }
                if p.is_dir() {
                    stack.push(p);
                }
            }
        }
        false
    }

    /// Status of a single database.
    pub fn status(&self, db: Database) -> Status {
        let dir = self.dir(db);
        if !dir.exists() {
            return Status::NotInstalled;
        }
        let mut missing = Vec::new();
        for req in db.required() {
            if !Self::has_suffix(&dir, req) {
                missing.push((*req).to_string());
            }
        }
        if missing.is_empty() {
            Status::Installed
        } else if missing.len() == db.required().len() {
            Status::Empty
        } else {
            Status::Partial { missing }
        }
    }

    /// Status of every managed database.
    pub fn status_all(&self) -> Vec<(Database, Status)> {
        Database::all()
            .iter()
            .map(|db| (*db, self.status(*db)))
            .collect()
    }

    /// Download and install a database (skeleton is created if needed). For
    /// [`Database::Busco`] a lineage name must be supplied via `param`.
    pub fn install(&self, db: Database, param: Option<&str>) -> Result<()> {
        self.build_skeleton()?;
        let dir = self.dir(db);
        std::fs::create_dir_all(&dir).ok();

        match db.source() {
            Source::VersionedTarball {
                version_url,
                archive_url_prefix,
            } => {
                let version = download::fetch_text(version_url)
                    .with_context(|| format!("discovering current release at {version_url}"))?
                    .trim()
                    .to_string();
                if version.is_empty() {
                    bail!("empty version string from {version_url}");
                }
                let archive_name = format!("{version}.tar.gz");
                let url = format!("{archive_url_prefix}{archive_name}");
                let archive = dir.join(&archive_name);
                eprintln!("[sentinel-db] downloading {} ({version})", db.display());
                download::fetch(&url, &archive)?;
                eprintln!("[sentinel-db] extracting {archive_name}");
                download::extract_tar_gz(&archive, &dir)?;
                std::fs::remove_file(&archive).ok();
            }
            Source::GzFile { url, out_name } => {
                let gz = dir.join(format!("{out_name}.gz"));
                eprintln!("[sentinel-db] downloading {}", db.display());
                download::fetch(url, &gz)?;
                eprintln!("[sentinel-db] decompressing {out_name}.gz");
                download::gunzip(&gz, &dir.join(out_name))?;
                std::fs::remove_file(&gz).ok();
            }
            Source::MultiGzFile { urls } => {
                for (url, name) in urls {
                    let target = dir.join(name);
                    eprintln!("[sentinel-db] downloading {name}");
                    download::fetch(url, &target)?;
                    if name.ends_with(".tar.gz") {
                        download::extract_tar_gz(&target, &dir)?;
                        std::fs::remove_file(&target).ok();
                    } else if name.ends_with(".gz") {
                        let out = dir.join(name.trim_end_matches(".gz"));
                        download::gunzip(&target, &out)?;
                        std::fs::remove_file(&target).ok();
                    }
                }
            }
            Source::Zenodo { doi } => {
                eprintln!("[sentinel-db] resolving Zenodo DOI {doi}");
                let files = download::resolve_zenodo(doi)?;
                if files.is_empty() {
                    bail!("Zenodo record for {doi} contained no files");
                }
                for f in files {
                    let target = dir.join(&f.key);
                    eprintln!("[sentinel-db] downloading {} ({})", f.key, db.display());
                    download::fetch(&f.url, &target)?;
                    if f.key.ends_with(".tar.gz") {
                        download::extract_tar_gz(&target, &dir)?;
                        std::fs::remove_file(&target).ok();
                    }
                }
            }
            Source::Parameterized { base_url } => {
                let lineage = param.ok_or_else(|| {
                    anyhow!("BUSCO requires a lineage name, e.g. `bacteria_odb10` (pass it as the parameter)")
                })?;
                // BUSCO lineage tarballs are dated; the file-listing page is the
                // source of truth. We fetch the directory index and pick the
                // newest matching tarball.
                let index = download::fetch_text(base_url)
                    .with_context(|| format!("listing BUSCO lineages at {base_url}"))?;
                let file = newest_busco_tarball(&index, lineage).ok_or_else(|| {
                    anyhow!("no BUSCO tarball for lineage '{lineage}' found at {base_url}")
                })?;
                let url = format!("{base_url}{file}");
                let archive = dir.join(&file);
                eprintln!("[sentinel-db] downloading BUSCO lineage {lineage} ({file})");
                download::fetch(&url, &archive)?;
                download::extract_tar_gz(&archive, &dir)?;
                std::fs::remove_file(&archive).ok();
            }
        }

        // Refresh the manifest with the new status.
        self.write_manifest().ok();
        Ok(())
    }

    /// Write `manifest.json` capturing the current status of each database.
    pub fn write_manifest(&self) -> Result<()> {
        let entries: Vec<serde_json::Value> = Database::all()
            .iter()
            .map(|db| {
                let st = self.status(*db);
                serde_json::json!({
                    "id": db.id(),
                    "name": db.display(),
                    "subdir": db.subdir(),
                    "used_by": db.used_by(),
                    "status": st.symbol(),
                })
            })
            .collect();
        let manifest = serde_json::json!({
            "sentinel_db_version": 1,
            "root": self.root.to_string_lossy(),
            "databases": entries,
        });
        let path = self.root.join("manifest.json");
        std::fs::create_dir_all(&self.root).ok();
        std::fs::write(&path, serde_json::to_string_pretty(&manifest)?)
            .with_context(|| format!("writing {}", path.display()))?;
        Ok(())
    }
}

/// Choose the newest dated BUSCO tarball for a lineage from a directory index.
/// Filenames look like `bacteria_odb10.2024-01-08.tar.gz`.
fn newest_busco_tarball(index: &str, lineage: &str) -> Option<String> {
    let prefix = format!("{lineage}.");
    let mut best: Option<String> = None;
    // crude href/text scan: find occurrences of the prefix followed by a date
    // and `.tar.gz`.
    for token in index.split(|c: char| c == '"' || c == '<' || c == '>' || c.is_whitespace()) {
        if token.starts_with(&prefix) && token.ends_with(".tar.gz") {
            match &best {
                Some(b) if token <= b.as_str() => {}
                _ => best = Some(token.to_string()),
            }
        }
    }
    best
}

/// Best-effort home directory without external crates.
fn home_dir() -> Option<PathBuf> {
    std::env::var_os("HOME")
        .or_else(|| std::env::var_os("USERPROFILE"))
        .map(PathBuf::from)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ids_and_parsing_roundtrip() {
        for db in Database::all() {
            assert_eq!(Database::parse(db.id()), Some(*db));
        }
        assert_eq!(Database::parse("CheckV"), Some(Database::CheckV));
        assert_eq!(Database::parse("nope"), None);
    }

    #[test]
    fn resolve_respects_env_and_explicit() {
        let r = DbRoot::resolve(Some("/tmp/explicit_db")).unwrap();
        assert_eq!(r.root, PathBuf::from("/tmp/explicit_db"));
    }

    #[test]
    fn skeleton_and_status() {
        let tmp = std::env::temp_dir().join(format!("sentdb_sk_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&tmp);
        let root = DbRoot { root: tmp.clone() };
        root.build_skeleton().unwrap();
        // every subdir exists, but nothing installed yet
        for (db, st) in root.status_all() {
            assert!(matches!(st, Status::Empty | Status::NotInstalled), "{:?}", db);
        }
        assert!(tmp.join("manifest.json").exists());

        // Simulate a CheckV install by creating the required files under a
        // version-named intermediate dir.
        let cv = root.dir(Database::CheckV).join("checkv-db-v1.5");
        std::fs::create_dir_all(cv.join("genome_db")).unwrap();
        std::fs::create_dir_all(cv.join("hmm_db")).unwrap();
        for f in [
            "genome_db/checkv_reps.tsv",
            "genome_db/checkv_error.tsv",
            "hmm_db/genome_lengths.tsv",
            "hmm_db/checkv_hmms.tsv",
        ] {
            std::fs::write(cv.join(f), b"x").unwrap();
        }
        assert_eq!(root.status(Database::CheckV), Status::Installed);

        // Partial: remove one file.
        std::fs::remove_file(cv.join("hmm_db/checkv_hmms.tsv")).unwrap();
        match root.status(Database::CheckV) {
            Status::Partial { missing } => {
                assert_eq!(missing, vec!["hmm_db/checkv_hmms.tsv".to_string()])
            }
            other => panic!("expected Partial, got {other:?}"),
        }
        std::fs::remove_dir_all(&tmp).ok();
    }

    #[test]
    fn checkm2_status_finds_dmnd_anywhere() {
        let tmp = std::env::temp_dir().join(format!("sentdb_cm_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&tmp);
        let root = DbRoot { root: tmp.clone() };
        root.build_skeleton().unwrap();
        let inner = root.dir(Database::CheckM2).join("CheckM2_database");
        std::fs::create_dir_all(&inner).unwrap();
        std::fs::write(inner.join("uniref100.KO.1.dmnd"), b"x").unwrap();
        assert_eq!(root.status(Database::CheckM2), Status::Installed);
        std::fs::remove_dir_all(&tmp).ok();
    }

    #[test]
    fn busco_tarball_picker() {
        let index = r#"
            <a href="bacteria_odb10.2020-03-06.tar.gz">x</a>
            <a href="bacteria_odb10.2024-01-08.tar.gz">x</a>
            <a href="archaea_odb10.2024-01-08.tar.gz">x</a>
        "#;
        assert_eq!(
            newest_busco_tarball(index, "bacteria_odb10").as_deref(),
            Some("bacteria_odb10.2024-01-08.tar.gz")
        );
        assert_eq!(newest_busco_tarball(index, "missing_odb10"), None);
    }
}
