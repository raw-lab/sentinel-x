//! Core sequence types: the biological [`Domain`], a [`Contig`], and FASTA
//! reading (plain or gzip-compressed, detected by magic bytes).

use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// The kingdom/domain a sequence is being validated as. Sentinel routes each
/// contig to the appropriate validators based on this (a `Virus` goes to
/// CheckV + viralVerify + TRACE; `Bacteria`/`Archaea` to CheckM2 + SAGE; a
/// `Eukaryote` to EukCC/BUSCO).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Domain {
    Virus,
    Bacteria,
    Archaea,
    Eukaryote,
    /// Mixed / unknown input (a metagenome or a set of viral contigs of
    /// unknown origin): the per-contig validators run on every contig.
    Mixed,
}

impl Domain {
    /// Parse from a CLI string. Single-letter abbreviations are not accepted
    /// (`-d a` used to mean archaea, which surprised users expecting "auto").
    pub fn parse(s: &str) -> Option<Domain> {
        match s.to_ascii_lowercase().as_str() {
            "virus" | "viral" | "viruses" => Some(Domain::Virus),
            "bacteria" | "bacterium" | "bact" => Some(Domain::Bacteria),
            "archaea" | "archaeon" | "arch" => Some(Domain::Archaea),
            "eukaryote" | "eukaryotic" | "euk" => Some(Domain::Eukaryote),
            "mixed" | "meta" | "metagenome" => Some(Domain::Mixed),
            _ => None,
        }
    }

    /// Human label.
    pub fn label(&self) -> &'static str {
        match self {
            Domain::Virus => "virus",
            Domain::Bacteria => "bacteria",
            Domain::Archaea => "archaea",
            Domain::Eukaryote => "eukaryote",
            Domain::Mixed => "mixed",
        }
    }

    /// Map to the DYNAMITE gene-calling preset most appropriate for the domain
    /// (prodigal-gv for viral and mixed input, as CheckV uses).
    pub fn dynamite_mode(&self) -> dynamite::model::Mode {
        use dynamite::model::Mode;
        match self {
            Domain::Virus | Domain::Mixed => Mode::Virus,
            Domain::Bacteria => Mode::Bacteria,
            Domain::Archaea => Mode::Archaea,
            Domain::Eukaryote => Mode::Eukaryote,
        }
    }
}

/// A single input sequence.
#[derive(Clone, Debug)]
pub struct Contig {
    /// Identifier (first whitespace-delimited token of the FASTA header).
    pub id: String,
    /// Full FASTA header line (without `>`).
    pub header: String,
    /// Uppercased nucleotide sequence (`ACGTN`...).
    pub seq: Vec<u8>,
}

impl Contig {
    /// Sequence length in bases.
    pub fn len(&self) -> usize {
        self.seq.len()
    }
}

/// Open a path as a buffered reader, transparently decompressing gzip
/// (recognised by its `1f 8b` magic, whatever the file is called).
pub fn open_maybe_gz(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut br = BufReader::with_capacity(1 << 16, file);
    let is_gz = {
        let buf = br.fill_buf().with_context(|| format!("reading {}", path.display()))?;
        buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b
    };
    if is_gz {
        Ok(Box::new(BufReader::with_capacity(1 << 16, flate2::read::MultiGzDecoder::new(br))))
    } else {
        Ok(Box::new(br))
    }
}

/// Read all records from a FASTA file (plain or gzip). See
/// [`read_fasta_reader`] for the validation rules.
pub fn read_fasta<P: AsRef<Path>>(path: P) -> Result<Vec<Contig>> {
    let p = path.as_ref();
    let reader = open_maybe_gz(p)?;
    read_fasta_reader(reader).with_context(|| format!("parsing FASTA {}", p.display()))
}

/// Read FASTA from any reader, as raw bytes (non-UTF-8 headers are kept
/// lossily). Sequence letters are uppercased and whitespace stripped
/// (`\r\n` line ends are fine).
///
/// Rejected, as CheckV does: text before the first `>`, records with an
/// empty identifier, and duplicate identifiers. Records with an empty
/// sequence are skipped with a warning.
pub fn read_fasta_reader<R: BufRead>(mut reader: R) -> Result<Vec<Contig>> {
    let mut contigs: Vec<Contig> = Vec::new();
    let mut header: Option<String> = None;
    let mut seq: Vec<u8> = Vec::new();
    let mut line: Vec<u8> = Vec::new();
    let mut lineno = 0usize;
    let mut seen: HashSet<String> = HashSet::new();
    let mut skipped_empty = 0usize;

    let mut finish = |h: String, s: Vec<u8>, contigs: &mut Vec<Contig>, ln: usize| -> Result<()> {
        let id = h.split_whitespace().next().unwrap_or("").to_string();
        if id.is_empty() {
            bail!("record ending at line {ln} has an empty identifier ('>' with no name)");
        }
        if !seen.insert(id.clone()) {
            bail!("duplicate sequence identifier '{id}' (identifiers must be unique)");
        }
        if s.is_empty() {
            skipped_empty += 1;
            return Ok(());
        }
        contigs.push(Contig { id, header: h, seq: s });
        Ok(())
    };

    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line).context("reading FASTA")?;
        if n == 0 {
            break;
        }
        lineno += 1;
        while matches!(line.last(), Some(b'\n') | Some(b'\r')) {
            line.pop();
        }
        if let Some(rest) = line.strip_prefix(b">") {
            if let Some(h) = header.take() {
                finish(h, std::mem::take(&mut seq), &mut contigs, lineno - 1)?;
            }
            header = Some(String::from_utf8_lossy(rest).into_owned());
        } else if header.is_none() {
            if line.iter().any(|b| !b.is_ascii_whitespace()) {
                bail!("line {lineno}: sequence data before the first '>' header (not a FASTA file?)");
            }
        } else {
            seq.extend(line.iter().filter(|b| !b.is_ascii_whitespace()).map(|b| b.to_ascii_uppercase()));
        }
    }
    if let Some(h) = header.take() {
        finish(h, seq, &mut contigs, lineno)?;
    }
    if skipped_empty > 0 {
        eprintln!("[sentinel] skipped {skipped_empty} FASTA record(s) with an empty sequence");
    }
    Ok(contigs)
}

/// Read a whole (possibly gzip-compressed) file into memory.
pub fn read_all(path: &Path) -> Result<Vec<u8>> {
    let mut r = open_maybe_gz(path)?;
    let mut v = Vec::new();
    r.read_to_end(&mut v).with_context(|| format!("reading {}", path.display()))?;
    Ok(v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_multi_record_fasta() {
        let fa = ">c1 desc here\r\nACGT\r\nacgt\n\n>c2\nNNNN\n";
        let cs = read_fasta_reader(fa.as_bytes()).unwrap();
        assert_eq!(cs.len(), 2);
        assert_eq!(cs[0].id, "c1");
        assert_eq!(cs[0].header, "c1 desc here");
        assert_eq!(cs[0].seq, b"ACGTACGT");
        assert_eq!(cs[1].id, "c2");
        assert_eq!(cs[1].len(), 4);
    }

    #[test]
    fn rejects_malformed_input() {
        assert!(read_fasta_reader(&b"ACGT\n>c1\nACGT\n"[..]).is_err());
        assert!(read_fasta_reader(&b">c1\nACGT\n>c1 again\nAC\n"[..]).is_err());
        assert!(read_fasta_reader(&b">\nACGT\n"[..]).is_err());
        // empty records are skipped, not fatal
        let cs = read_fasta_reader(&b">e\n>c1\nAC\n"[..]).unwrap();
        assert_eq!(cs.len(), 1);
        // non-UTF-8 header bytes are tolerated
        let cs = read_fasta_reader(&b">c\xff1\nAC\n"[..]).unwrap();
        assert_eq!(cs[0].seq, b"AC");
    }

    #[test]
    fn gzip_by_magic() {
        use std::io::Write;
        let dir = std::env::temp_dir().join(format!("sentinel_gz_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("contigs.fasta"); // no .gz extension
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        enc.write_all(b">g1\nACGTN\n").unwrap();
        std::fs::write(&p, enc.finish().unwrap()).unwrap();
        let cs = read_fasta(&p).unwrap();
        assert_eq!(cs[0].seq, b"ACGTN");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn domain_roundtrip() {
        assert_eq!(Domain::parse("viral"), Some(Domain::Virus));
        assert_eq!(Domain::parse("ARCH"), Some(Domain::Archaea));
        assert_eq!(Domain::parse("a"), None);
        assert_eq!(Domain::Bacteria.label(), "bacteria");
    }
}
