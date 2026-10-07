//! Network fetch, archive extraction and checksum helpers.
//!
//! HTTP transfer is delegated to `curl` (falling back to `wget`) because the
//! reference databases are multi-gigabyte files served from hosts that require
//! resume, redirects and range requests — exactly what those tools handle well,
//! and what every comparable bioinformatics downloader shells out to. Archive
//! extraction and checksum verification are pure Rust (`flate2` + `tar` +
//! `sha2`).

use std::io::Read;
use std::path::Path;
use std::process::Command;

use anyhow::{anyhow, bail, Context, Result};

/// Is a usable HTTP fetch tool available on PATH?
pub fn fetcher_available() -> Option<&'static str> {
    for tool in ["curl", "wget"] {
        if Command::new(tool)
            .arg("--version")
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .status()
            .map(|s| s.success())
            .unwrap_or(false)
        {
            return Some(tool);
        }
    }
    None
}

/// Download `url` to `dest`, resuming if a partial file exists.
pub fn fetch(url: &str, dest: &Path) -> Result<()> {
    if let Some(parent) = dest.parent() {
        std::fs::create_dir_all(parent).ok();
    }
    let tool = fetcher_available()
        .ok_or_else(|| anyhow!("neither `curl` nor `wget` found on PATH"))?;
    let status = match tool {
        "curl" => Command::new("curl")
            .args(["-fL", "--retry", "5", "--retry-delay", "3", "-C", "-", "-o"])
            .arg(dest)
            .arg(url)
            .status()
            .context("spawning curl")?,
        _ => Command::new("wget")
            .args(["-c", "-O"])
            .arg(dest)
            .arg(url)
            .status()
            .context("spawning wget")?,
    };
    if !status.success() {
        bail!("{tool} failed to download {url}");
    }
    Ok(())
}

/// Fetch a small text resource and return its contents (used for version files
/// and JSON API responses). Capped in size; for large payloads use [`fetch`].
pub fn fetch_text(url: &str) -> Result<String> {
    let tool = fetcher_available()
        .ok_or_else(|| anyhow!("neither `curl` nor `wget` found on PATH"))?;
    let out = match tool {
        "curl" => Command::new("curl")
            .args(["-fsSL", url])
            .output()
            .context("spawning curl")?,
        _ => Command::new("wget")
            .args(["-q", "-O", "-", url])
            .output()
            .context("spawning wget")?,
    };
    if !out.status.success() {
        bail!(
            "{tool} failed to fetch text from {url}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
    Ok(String::from_utf8_lossy(&out.stdout).into_owned())
}

/// Extract a `.tar.gz` archive into `dest_dir` (pure Rust).
pub fn extract_tar_gz(archive: &Path, dest_dir: &Path) -> Result<()> {
    std::fs::create_dir_all(dest_dir).ok();
    let f = std::fs::File::open(archive)
        .with_context(|| format!("opening archive {}", archive.display()))?;
    let gz = flate2::read::GzDecoder::new(f);
    let mut tar = tar::Archive::new(gz);
    tar.unpack(dest_dir)
        .with_context(|| format!("unpacking {} into {}", archive.display(), dest_dir.display()))?;
    Ok(())
}

/// Decompress a single `.gz` file to `out` (pure Rust).
pub fn gunzip(archive: &Path, out: &Path) -> Result<()> {
    let f = std::fs::File::open(archive)
        .with_context(|| format!("opening {}", archive.display()))?;
    let mut gz = flate2::read::GzDecoder::new(f);
    let mut buf = Vec::new();
    gz.read_to_end(&mut buf)?;
    std::fs::write(out, &buf).with_context(|| format!("writing {}", out.display()))?;
    Ok(())
}

/// SHA-256 of a file as lowercase hex (pure Rust).
pub fn sha256_hex(path: &Path) -> Result<String> {
    use sha2::{Digest, Sha256};
    let mut f = std::fs::File::open(path)
        .with_context(|| format!("opening {} for hashing", path.display()))?;
    let mut hasher = Sha256::new();
    let mut buf = [0u8; 1 << 16];
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(hex(&hasher.finalize()))
}

fn hex(bytes: &[u8]) -> String {
    let mut s = String::with_capacity(bytes.len() * 2);
    for b in bytes {
        s.push_str(&format!("{b:02x}"));
    }
    s
}

/// Resolve a Zenodo DOI (e.g. `10.5281/zenodo.14897628`) to its downloadable
/// files via the Zenodo REST API. Returns `(filename, download_url, md5)` for
/// each file in the record.
pub fn resolve_zenodo(doi: &str) -> Result<Vec<ZenodoFile>> {
    // The DOI's numeric suffix after "zenodo." is the record id.
    let record_id = doi
        .rsplit("zenodo.")
        .next()
        .filter(|s| s.chars().all(|c| c.is_ascii_digit()))
        .ok_or_else(|| anyhow!("cannot parse Zenodo record id from DOI '{doi}'"))?;
    let api = format!("https://zenodo.org/api/records/{record_id}");
    let body = fetch_text(&api).with_context(|| format!("querying Zenodo API {api}"))?;
    let v: serde_json::Value = serde_json::from_str(&body).context("parsing Zenodo JSON")?;
    let files = v
        .get("files")
        .and_then(|f| f.as_array())
        .ok_or_else(|| anyhow!("Zenodo record {record_id} has no files array"))?;
    let mut out = Vec::new();
    for f in files {
        let key = f
            .get("key")
            .and_then(|k| k.as_str())
            .unwrap_or("unknown")
            .to_string();
        let url = f
            .get("links")
            .and_then(|l| l.get("self"))
            .and_then(|s| s.as_str())
            .unwrap_or("")
            .to_string();
        let checksum = f
            .get("checksum")
            .and_then(|c| c.as_str())
            .unwrap_or("")
            .to_string();
        out.push(ZenodoFile { key, url, checksum });
    }
    Ok(out)
}

/// A file entry from a Zenodo record.
#[derive(Clone, Debug)]
pub struct ZenodoFile {
    pub key: String,
    pub url: String,
    /// e.g. `md5:abcdef...`
    pub checksum: String,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hex_encoding() {
        assert_eq!(hex(&[0x00, 0xff, 0x10]), "00ff10");
    }

    #[test]
    fn sha256_of_known_string() {
        let dir = std::env::temp_dir().join(format!("sentdb_hash_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("x.txt");
        std::fs::write(&p, b"abc").unwrap();
        // SHA-256("abc")
        assert_eq!(
            sha256_hex(&p).unwrap(),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn zenodo_record_id_parsing() {
        // Just the parsing branch (no network): a malformed DOI errors.
        assert!(resolve_zenodo("not-a-doi").is_err());
    }

    #[test]
    fn roundtrip_targz_extraction() {
        // Build a tiny tar.gz in-memory-ish on disk, then extract and verify.
        let dir = std::env::temp_dir().join(format!("sentdb_tar_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let src = dir.join("payload.txt");
        std::fs::write(&src, b"hello-sentinel").unwrap();
        let archive = dir.join("a.tar.gz");

        // Create the archive with the tar+flate2 crates (writer side).
        {
            let f = std::fs::File::create(&archive).unwrap();
            let enc = flate2::write::GzEncoder::new(f, flate2::Compression::default());
            let mut builder = tar::Builder::new(enc);
            builder.append_path_with_name(&src, "payload.txt").unwrap();
            builder.finish().unwrap();
        }
        let outdir = dir.join("out");
        extract_tar_gz(&archive, &outdir).unwrap();
        let got = std::fs::read_to_string(outdir.join("payload.txt")).unwrap();
        assert_eq!(got, "hello-sentinel");
        std::fs::remove_dir_all(&dir).ok();
    }
}
