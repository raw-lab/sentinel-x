//! Network fetch, archive extraction and checksum helpers.
//!
//! HTTP transfer is delegated to `curl` (falling back to `wget`): the
//! reference databases are multi-gigabyte files that need resume, redirects
//! and retries, and both tools honour the system proxy settings
//! (`HTTPS_PROXY`, CA bundles). Extraction, decompression and hashing are
//! pure Rust (`flate2`, `tar`, `sha2`, `md-5`).
//!
//! Every download goes to a `*.part` file that is renamed only once complete,
//! and archives are unpacked entry by entry: symbolic links, hard links,
//! absolute paths and `..` components are rejected, so an archive can never
//! write outside its destination or plant links that later directory walks
//! would follow.

use std::io::{Read, Write};
use std::path::{Component, Path, PathBuf};
use std::process::Command;

use anyhow::{anyhow, bail, Context, Result};

/// Is a usable HTTP fetch tool available on PATH?
pub fn fetcher_available() -> Option<&'static str> {
    ["curl", "wget"].into_iter().find(|tool| {
        Command::new(tool)
            .arg("--version")
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .status()
            .map(|s| s.success())
            .unwrap_or(false)
    })
}

fn tool() -> Result<&'static str> {
    fetcher_available().ok_or_else(|| anyhow!("neither `curl` nor `wget` found on PATH"))
}

/// Turn a fetch failure into a readable message (HTTP status, proxy denial).
fn fetch_error(tool: &str, url: &str, stderr: &[u8]) -> anyhow::Error {
    let msg = String::from_utf8_lossy(stderr).trim().to_string();
    let hint = if msg.contains("403") {
        " — access denied (HTTP 403); the host may be blocked by a proxy or firewall"
    } else if msg.contains("404") {
        " — not found (HTTP 404)"
    } else {
        ""
    };
    anyhow!("{tool} could not fetch {url}{hint}: {msg}")
}

/// Download `url` to `dest`, resuming a previous `dest.part` if present. The
/// final file appears only when the transfer completed.
pub fn fetch(url: &str, dest: &Path) -> Result<()> {
    if let Some(parent) = dest.parent() {
        std::fs::create_dir_all(parent).with_context(|| format!("creating {}", parent.display()))?;
    }
    let part = part_path(dest);
    let t = tool()?;
    let out = match t {
        "curl" => Command::new("curl")
            .args(["-fL", "--retry", "5", "--retry-delay", "3", "--connect-timeout", "30", "-C", "-", "-sS", "-o"])
            .arg(&part)
            .arg(url)
            .output()
            .context("spawning curl")?,
        _ => Command::new("wget")
            .args(["-q", "-c", "--tries=5", "--timeout=60", "-O"])
            .arg(&part)
            .arg(url)
            .output()
            .context("spawning wget")?,
    };
    if !out.status.success() {
        return Err(fetch_error(t, url, &out.stderr));
    }
    std::fs::rename(&part, dest).with_context(|| format!("finalising {}", dest.display()))?;
    Ok(())
}

/// `<file>.part` next to `dest`.
pub fn part_path(dest: &Path) -> PathBuf {
    let mut s = dest.as_os_str().to_owned();
    s.push(".part");
    PathBuf::from(s)
}

/// Fetch a small text resource (version files, manifests, JSON API
/// responses), with retries and a time limit.
pub fn fetch_text(url: &str) -> Result<String> {
    let t = tool()?;
    let out = match t {
        "curl" => Command::new("curl")
            .args(["-fsSL", "--retry", "3", "--retry-delay", "2", "--connect-timeout", "20", "--max-time", "300", url])
            .output()
            .context("spawning curl")?,
        _ => Command::new("wget").args(["-q", "--tries=3", "--timeout=60", "-O", "-", url]).output().context("spawning wget")?,
    };
    if !out.status.success() {
        return Err(fetch_error(t, url, &out.stderr));
    }
    Ok(String::from_utf8_lossy(&out.stdout).into_owned())
}

/// Follow redirects and return the final URL (used to turn a DOI into its
/// Zenodo record, as CheckM2 does).
pub fn resolve_redirect(url: &str) -> Result<String> {
    let t = tool()?;
    let out = match t {
        "curl" => Command::new("curl")
            .args([
                "-fsSL",
                "--retry",
                "3",
                "--connect-timeout",
                "20",
                "--max-time",
                "120",
                "-o",
                "/dev/null",
                "-w",
                "%{url_effective}",
                url,
            ])
            .output()
            .context("spawning curl")?,
        _ => bail!("resolving redirects needs curl"),
    };
    if !out.status.success() {
        return Err(fetch_error(t, url, &out.stderr));
    }
    Ok(String::from_utf8_lossy(&out.stdout).trim().to_string())
}

/// Reject archive entry paths that are absolute or climb out of the
/// destination.
fn safe_relative(p: &Path) -> Result<PathBuf> {
    let mut out = PathBuf::new();
    for c in p.components() {
        match c {
            Component::Normal(x) => out.push(x),
            Component::CurDir => {}
            _ => bail!("unsafe path in archive: {}", p.display()),
        }
    }
    Ok(out)
}

/// Extract a `.tar.gz` (or plain `.tar`) archive into `dest_dir`, refusing
/// links and unsafe paths. Returns the number of files written.
pub fn extract_tar_gz(archive: &Path, dest_dir: &Path) -> Result<usize> {
    std::fs::create_dir_all(dest_dir).with_context(|| format!("creating {}", dest_dir.display()))?;
    let reader = crate::common::seqio::open_maybe_gz(archive)?;
    let mut tar = tar::Archive::new(reader);
    tar.set_preserve_permissions(false);
    let mut n = 0usize;
    for entry in tar.entries().with_context(|| format!("reading {}", archive.display()))? {
        let mut entry = entry.with_context(|| format!("reading {}", archive.display()))?;
        let path = entry.path()?.into_owned();
        let rel = safe_relative(&path)?;
        let kind = entry.header().entry_type();
        if kind.is_symlink() || kind.is_hard_link() {
            bail!("{}: refusing link entry {}", archive.display(), path.display());
        }
        if !(kind.is_file() || kind.is_dir()) {
            continue; // devices, fifos, pax headers
        }
        if rel.as_os_str().is_empty() {
            continue;
        }
        let ok = entry.unpack_in(dest_dir).with_context(|| format!("unpacking {}", path.display()))?;
        if !ok {
            bail!("{}: entry {} escapes the destination", archive.display(), path.display());
        }
        if kind.is_file() {
            n += 1;
        }
    }
    Ok(n)
}

/// Decompress a gzip file to `out` by streaming (via `out.part`).
pub fn gunzip(archive: &Path, out: &Path) -> Result<u64> {
    let f = std::fs::File::open(archive).with_context(|| format!("opening {}", archive.display()))?;
    let mut gz = flate2::read::MultiGzDecoder::new(std::io::BufReader::with_capacity(1 << 20, f));
    let part = part_path(out);
    let mut w = std::io::BufWriter::with_capacity(1 << 20, std::fs::File::create(&part)?);
    let n = std::io::copy(&mut gz, &mut w).with_context(|| format!("decompressing {}", archive.display()))?;
    w.flush()?;
    drop(w);
    std::fs::rename(&part, out)?;
    Ok(n)
}

fn hash_file<D: sha2::Digest>(path: &Path, mut h: D) -> Result<String> {
    let mut f = std::fs::File::open(path).with_context(|| format!("opening {} for hashing", path.display()))?;
    let mut buf = vec![0u8; 1 << 20];
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        h.update(&buf[..n]);
    }
    Ok(hex(&h.finalize()))
}

/// SHA-256 of a file as lowercase hex.
pub fn sha256_hex(path: &Path) -> Result<String> {
    hash_file(path, sha2::Sha256::default())
}

/// MD5 of a file as lowercase hex (Zenodo, BUSCO and Pfam publish MD5).
pub fn md5_hex(path: &Path) -> Result<String> {
    hash_file(path, md5::Md5::default())
}

fn hex(bytes: &[u8]) -> String {
    let mut s = String::with_capacity(bytes.len() * 2);
    for b in bytes {
        s.push_str(&format!("{b:02x}"));
    }
    s
}

/// A file entry from a Zenodo record.
#[derive(Clone, Debug)]
pub struct ZenodoFile {
    /// Base name of the record's key (never a path).
    pub name: String,
    pub url: String,
    /// `algorithm:hex`, e.g. `md5:abcdef…`.
    pub checksum: String,
}

/// A Zenodo record: its `metadata.version` and files.
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct ZenodoRecord {
    pub id: String,
    pub version: String,
    pub files: Vec<ZenodoFile>,
}

/// Parse a Zenodo `api/records/<id>` response.
pub fn parse_zenodo_record(id: &str, body: &str) -> Result<ZenodoRecord> {
    let v: serde_json::Value = serde_json::from_str(body).context("parsing Zenodo JSON")?;
    let files = v.get("files").and_then(|f| f.as_array()).ok_or_else(|| anyhow!("Zenodo record {id} has no files array"))?;
    let mut out = Vec::new();
    for f in files {
        let key = f.get("key").and_then(|k| k.as_str()).unwrap_or("");
        let name = key.rsplit('/').next().unwrap_or("").to_string();
        if name.is_empty() || name == "." || name == ".." {
            bail!("Zenodo record {id}: unusable file key '{key}'");
        }
        let url = f.get("links").and_then(|l| l.get("self")).and_then(|s| s.as_str()).unwrap_or("").to_string();
        let checksum = f.get("checksum").and_then(|c| c.as_str()).unwrap_or("").to_string();
        out.push(ZenodoFile { name, url, checksum });
    }
    let version = v
        .get("metadata")
        .and_then(|m| m.get("version"))
        .map(|x| match x {
            serde_json::Value::String(s) => s.clone(),
            other => other.to_string(),
        })
        .unwrap_or_default();
    Ok(ZenodoRecord { id: id.to_string(), version, files: out })
}

/// The record id at the end of a resolved DOI URL
/// (`https://zenodo.org/records/14897628` → `14897628`).
pub fn zenodo_record_id(resolved: &str) -> Result<String> {
    let tail = resolved.split(['?', '#']).next().unwrap_or("").trim_end_matches('/');
    let id = tail.rsplit('/').next().unwrap_or("");
    if id.is_empty() || !id.chars().all(|c| c.is_ascii_digit()) {
        bail!("cannot read a Zenodo record id from '{resolved}'");
    }
    Ok(id.to_string())
}

/// Verify a file against an `algorithm:hex` checksum (md5 or sha256).
pub fn verify_checksum(path: &Path, checksum: &str) -> Result<()> {
    let (alg, want) = checksum.split_once(':').unwrap_or(("md5", checksum));
    let got = match alg.to_ascii_lowercase().as_str() {
        "md5" => md5_hex(path)?,
        "sha256" => sha256_hex(path)?,
        other => bail!("unsupported checksum algorithm '{other}'"),
    };
    if !got.eq_ignore_ascii_case(want.trim()) {
        bail!("{} checksum mismatch for {}: expected {}, got {}", alg, path.display(), want.trim(), got);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmpdir(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("sentdb_{tag}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn hashes() {
        let dir = tmpdir("hash");
        let p = dir.join("x.txt");
        std::fs::write(&p, b"abc").unwrap();
        assert_eq!(sha256_hex(&p).unwrap(), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
        assert_eq!(md5_hex(&p).unwrap(), "900150983cd24fb0d6963f7d28e17f72");
        assert!(verify_checksum(&p, "md5:900150983cd24fb0d6963f7d28e17f72").is_ok());
        assert!(verify_checksum(&p, "md5:00000000000000000000000000000000").is_err());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn zenodo_parsing() {
        assert_eq!(zenodo_record_id("https://zenodo.org/records/14897628").unwrap(), "14897628");
        assert_eq!(zenodo_record_id("https://zenodo.org/records/14897628/?x=1").unwrap(), "14897628");
        assert!(zenodo_record_id("https://example.org/doi/abc").is_err());
        let body =
            r#"{"metadata":{"version":"3"},"files":[{"key":"../../etc/x.tar.gz","links":{"self":"u"},"checksum":"md5:ab"}]}"#;
        let r = parse_zenodo_record("1", body).unwrap();
        assert_eq!(r.version, "3");
        assert_eq!(r.files[0].name, "x.tar.gz");
    }

    fn build_tar(path: &Path, entries: &[(&str, &[u8])], link: Option<(&str, &str)>) {
        let f = std::fs::File::create(path).unwrap();
        let enc = flate2::write::GzEncoder::new(f, flate2::Compression::fast());
        let mut b = tar::Builder::new(enc);
        for (name, data) in entries {
            let mut h = tar::Header::new_gnu();
            h.set_size(data.len() as u64);
            h.set_mode(0o644);
            h.set_cksum();
            b.append_data(&mut h, name, *data).unwrap();
        }
        if let Some((name, target)) = link {
            let mut h = tar::Header::new_gnu();
            h.set_entry_type(tar::EntryType::Symlink);
            h.set_size(0);
            b.append_link(&mut h, name, target).unwrap();
        }
        b.into_inner().unwrap().finish().unwrap();
    }

    #[test]
    fn safe_extraction() {
        let dir = tmpdir("tar");
        let a = dir.join("a.tar.gz");
        build_tar(&a, &[("top/payload.txt", b"hello")], None);
        let out = dir.join("out");
        assert_eq!(extract_tar_gz(&a, &out).unwrap(), 1);
        assert_eq!(std::fs::read(out.join("top/payload.txt")).unwrap(), b"hello");
        // symlink entries are refused
        let b = dir.join("b.tar.gz");
        build_tar(&b, &[("top/x.txt", b"x")], Some(("top/evil", "/")));
        assert!(extract_tar_gz(&b, &dir.join("out2")).is_err());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn streaming_gunzip() {
        let dir = tmpdir("gz");
        let gz = dir.join("f.gz");
        let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        e.write_all(&vec![b'A'; 100_000]).unwrap();
        std::fs::write(&gz, e.finish().unwrap()).unwrap();
        assert_eq!(gunzip(&gz, &dir.join("f")).unwrap(), 100_000);
        assert!(!part_path(&dir.join("f")).exists());
        std::fs::remove_dir_all(&dir).ok();
    }
}
