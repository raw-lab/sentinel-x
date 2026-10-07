//! Core sequence types: the biological [`Domain`], a [`Contig`], and FASTA
//! reading (plain or gzip-compressed).

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// The kingdom/domain a sequence is being validated as. Sentinel routes each
/// contig to the appropriate validators based on this (e.g. a `Virus` goes to
/// CheckV + viralVerify + TRACE; `Bacteria`/`Archaea` to CheckM2 + SAGE; a
/// `Eukaryote` to BUSCO).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Domain {
    Virus,
    Bacteria,
    Archaea,
    Eukaryote,
    /// Mixed / unknown input (a metagenome); run the full panel and let the
    /// per-domain validators speak for themselves.
    Mixed,
}

impl Domain {
    /// Parse from a CLI string.
    pub fn parse(s: &str) -> Option<Domain> {
        match s.to_ascii_lowercase().as_str() {
            "virus" | "viral" | "viruses" | "v" => Some(Domain::Virus),
            "bacteria" | "bacterium" | "bact" | "b" => Some(Domain::Bacteria),
            "archaea" | "archaeon" | "arch" | "a" => Some(Domain::Archaea),
            "eukaryote" | "eukaryotic" | "euk" | "e" => Some(Domain::Eukaryote),
            "mixed" | "meta" | "metagenome" | "m" | "auto" => Some(Domain::Mixed),
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

    /// Map to the DYNAMITE gene-calling preset most appropriate for the domain.
    pub fn dynamite_mode(&self) -> dynamite_core::model::Mode {
        use dynamite_core::model::Mode;
        match self {
            Domain::Virus => Mode::Virus,
            Domain::Bacteria => Mode::Bacteria,
            Domain::Archaea => Mode::Archaea,
            Domain::Eukaryote => Mode::Eukaryote,
            Domain::Mixed => Mode::Meta,
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
    /// True if the sequence is empty.
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
}

/// Open a path as a buffered reader, transparently decompressing `.gz`.
fn open_maybe_gz(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("opening {}", path.display()))?;
    let is_gz = path
        .extension()
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);
    if is_gz {
        let dec = flate2::read::MultiGzDecoder::new(file);
        Ok(Box::new(BufReader::new(dec)))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Read all records from a FASTA file (plain or gzip). Sequence letters are
/// uppercased; whitespace within sequence lines is stripped.
pub fn read_fasta<P: AsRef<Path>>(path: P) -> Result<Vec<Contig>> {
    let reader = open_maybe_gz(path.as_ref())?;
    read_fasta_reader(reader)
}

/// Read FASTA from any reader.
pub fn read_fasta_reader<R: BufRead>(reader: R) -> Result<Vec<Contig>> {
    let mut contigs = Vec::new();
    let mut header: Option<String> = None;
    let mut seq: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line.context("reading FASTA line")?;
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(h) = header.take() {
                contigs.push(finish(h, std::mem::take(&mut seq)));
            }
            header = Some(rest.to_string());
        } else {
            for &b in line.as_bytes() {
                if !b.is_ascii_whitespace() {
                    seq.push(b.to_ascii_uppercase());
                }
            }
        }
    }
    if let Some(h) = header.take() {
        contigs.push(finish(h, seq));
    }
    Ok(contigs)
}

fn finish(header: String, seq: Vec<u8>) -> Contig {
    let id = header
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_string();
    Contig { id, header, seq }
}

/// Convenience: read FASTA from raw bytes already in memory.
pub fn read_fasta_bytes(bytes: &[u8]) -> Result<Vec<Contig>> {
    read_fasta_reader(bytes)
}

// Allow reading from a generic `Read` (used in tests).
#[allow(dead_code)]
fn read_to_string<R: Read>(mut r: R) -> Result<String> {
    let mut s = String::new();
    r.read_to_string(&mut s)?;
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_multi_record_fasta() {
        let fa = ">c1 desc here\nACGT\nACGT\n>c2\nNNNN\n";
        let cs = read_fasta_reader(fa.as_bytes()).unwrap();
        assert_eq!(cs.len(), 2);
        assert_eq!(cs[0].id, "c1");
        assert_eq!(cs[0].seq, b"ACGTACGT");
        assert_eq!(cs[1].id, "c2");
        assert_eq!(cs[1].len(), 4);
    }

    #[test]
    fn domain_roundtrip() {
        assert_eq!(Domain::parse("viral"), Some(Domain::Virus));
        assert_eq!(Domain::parse("ARCH"), Some(Domain::Archaea));
        assert_eq!(Domain::Bacteria.label(), "bacteria");
    }
}
