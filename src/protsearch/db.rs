//! Reference sequence sources for [`super::search`], read in bounded blocks:
//! protein FASTA (plain or gzip, streamed), DIAMOND `.dmnd` (memory-mapped,
//! read in place — no conversion needed) and in-memory lists.

use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

use super::alphabet::{ENCODE, X};

/// A block of reference sequences (encoded, concatenated).
#[derive(Default)]
pub struct Block {
    /// Global index of the first sequence in the database.
    pub first: u64,
    pub ids: Vec<String>,
    pub seq: Vec<u8>,
    /// `off[i]..off[i+1]` is sequence `i` in `seq`.
    pub off: Vec<usize>,
}

impl Block {
    pub fn len(&self) -> usize {
        self.off.len().saturating_sub(1)
    }
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    #[inline]
    pub fn get(&self, i: usize) -> &[u8] {
        &self.seq[self.off[i]..self.off[i + 1]]
    }
    fn clear(&mut self) {
        self.ids.clear();
        self.seq.clear();
        self.off.clear();
        self.off.push(0);
    }
}

/// A database the search can stream.
pub trait SeqSource {
    /// Total residues (DIAMOND's database size for E-values).
    fn letters(&self) -> u64;
    /// Number of sequences, if known.
    fn sequences(&self) -> Option<u64>;
    /// Fill `out` with the next block of at most ~`max_letters` residues;
    /// returns `false` at the end of the database.
    fn next_block(&mut self, max_letters: usize, out: &mut Block) -> Result<bool>;
}

/// Open a reference database by path: `.dmnd` (by magic) or protein FASTA.
pub fn open(path: &Path) -> Result<Box<dyn SeqSource>> {
    let mut magic = [0u8; 8];
    {
        let mut f = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let n = f.read(&mut magic)?;
        if n == 8 && u64::from_le_bytes(magic) == DMND_MAGIC {
            return Ok(Box::new(DmndSource::open(path)?));
        }
    }
    Ok(Box::new(FastaSource::open(path)?))
}

// ---------------------------------------------------------------- FASTA

/// Streaming protein FASTA (plain or gzip by magic bytes).
pub struct FastaSource {
    path: PathBuf,
    reader: Option<Box<dyn BufRead + Send>>,
    letters: u64,
    nseq: u64,
    pending: Option<String>,
    next_index: u64,
    line: Vec<u8>,
}

fn open_text(path: &Path) -> Result<Box<dyn BufRead + Send>> {
    let mut f = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut m = [0u8; 2];
    let n = f.read(&mut m)?;
    drop(f);
    let f = std::fs::File::open(path)?;
    if n == 2 && m == [0x1f, 0x8b] {
        Ok(Box::new(BufReader::with_capacity(1 << 20, flate2::read::MultiGzDecoder::new(f))))
    } else {
        Ok(Box::new(BufReader::with_capacity(1 << 20, f)))
    }
}

impl FastaSource {
    /// Open and pre-count residues (needed up front for E-values).
    pub fn open(path: &Path) -> Result<FastaSource> {
        let mut r = open_text(path)?;
        let mut letters = 0u64;
        let mut nseq = 0u64;
        let mut line = Vec::new();
        let mut first = true;
        loop {
            line.clear();
            if r.read_until(b'\n', &mut line)? == 0 {
                break;
            }
            if line.first() == Some(&b'>') {
                nseq += 1;
                first = false;
            } else {
                if first && line.iter().any(|b| !b.is_ascii_whitespace()) {
                    bail!("{}: sequence data before the first '>' header", path.display());
                }
                letters += line.iter().filter(|b| !b.is_ascii_whitespace()).count() as u64;
            }
        }
        Ok(FastaSource {
            path: path.to_path_buf(),
            reader: Some(open_text(path)?),
            letters,
            nseq,
            pending: None,
            next_index: 0,
            line: Vec::new(),
        })
    }

    /// Restart from the beginning (multi-pass searches).
    pub fn rewind(&mut self) -> Result<()> {
        self.reader = Some(open_text(&self.path)?);
        self.pending = None;
        self.next_index = 0;
        Ok(())
    }
}

impl SeqSource for FastaSource {
    fn letters(&self) -> u64 {
        self.letters
    }
    fn sequences(&self) -> Option<u64> {
        Some(self.nseq)
    }
    fn next_block(&mut self, max_letters: usize, out: &mut Block) -> Result<bool> {
        out.clear();
        out.first = self.next_index;
        let Some(r) = self.reader.as_mut() else { return Ok(false) };
        if let Some(h) = self.pending.take() {
            out.ids.push(h);
        }
        loop {
            self.line.clear();
            if r.read_until(b'\n', &mut self.line)? == 0 {
                self.reader = None;
                break;
            }
            if self.line.first() == Some(&b'>') {
                let h = String::from_utf8_lossy(&self.line[1..]);
                let id = h.split_whitespace().next().unwrap_or("").to_string();
                if !out.ids.is_empty() {
                    out.off.push(out.seq.len());
                    if out.seq.len() >= max_letters {
                        self.pending = Some(id);
                        break;
                    }
                }
                out.ids.push(id);
            } else if !out.ids.is_empty() {
                out.seq.extend(self.line.iter().filter(|b| !b.is_ascii_whitespace()).map(|&b| ENCODE[b as usize]));
            }
        }
        if out.off.len() == out.ids.len() {
            out.off.push(out.seq.len());
        }
        self.next_index += out.ids.len() as u64;
        Ok(!out.ids.is_empty())
    }
}

// ---------------------------------------------------------------- DMND

const DMND_MAGIC: u64 = 0x024a_f8a4_15ee_186d;
/// DIAMOND's residue alphabet, index -> ASCII.
const DMND_ALPHABET: &[u8; 25] = b"ARNDCQEGHILKMFPSTWYVBJZX*";

/// DIAMOND v3 `.dmnd` database, memory-mapped. Layout (little endian):
/// magic u64 @0, db_version u32 @12 (= 3), sequences u64 @16, letters u64
/// @24, pos_array_offset u64 @32; records `0xff residues 0xff title 0x00`;
/// position array of N+1 `(u64 offset, u32 length, u32 pad)` entries.
pub struct DmndSource {
    map: memmap2::Mmap,
    nseq: u64,
    letters: u64,
    pos_off: usize,
    next: u64,
    remap: [u8; 128],
}

impl DmndSource {
    pub fn open(path: &Path) -> Result<DmndSource> {
        let f = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
        // SAFETY: read-only mapping of a file we do not modify.
        let map = unsafe { memmap2::Mmap::map(&f) }.with_context(|| format!("mapping {}", path.display()))?;
        if map.len() < 96 || u64::from_le_bytes(map[0..8].try_into().unwrap()) != DMND_MAGIC {
            bail!("{} is not a DIAMOND database", path.display());
        }
        let version = u32::from_le_bytes(map[12..16].try_into().unwrap());
        if version != 3 {
            bail!("{}: unsupported .dmnd db_version {version} (need 3, DIAMOND 2.x)", path.display());
        }
        let nseq = u64::from_le_bytes(map[16..24].try_into().unwrap());
        let letters = u64::from_le_bytes(map[24..32].try_into().unwrap());
        let pos_off = u64::from_le_bytes(map[32..40].try_into().unwrap()) as usize;
        let need = (nseq as usize).checked_add(1).and_then(|n| n.checked_mul(16)).and_then(|b| b.checked_add(pos_off));
        match need {
            Some(end) if end <= map.len() => {}
            _ => bail!("{}: truncated .dmnd (position array out of bounds)", path.display()),
        }
        let mut remap = [X; 128];
        for (i, &c) in DMND_ALPHABET.iter().enumerate() {
            remap[i] = ENCODE[c as usize];
        }
        Ok(DmndSource { map, nseq, letters, pos_off, next: 0, remap })
    }

    fn entry(&self, i: u64) -> (usize, usize) {
        let o = self.pos_off + 16 * i as usize;
        let off = u64::from_le_bytes(self.map[o..o + 8].try_into().unwrap()) as usize;
        let len = u32::from_le_bytes(self.map[o + 8..o + 12].try_into().unwrap()) as usize;
        (off, len)
    }
}

impl SeqSource for DmndSource {
    fn letters(&self) -> u64 {
        self.letters
    }
    fn sequences(&self) -> Option<u64> {
        Some(self.nseq)
    }
    fn next_block(&mut self, max_letters: usize, out: &mut Block) -> Result<bool> {
        out.clear();
        out.first = self.next;
        let data = &self.map[..];
        while self.next < self.nseq {
            let (off, len) = self.entry(self.next);
            let res = off.checked_add(1);
            let end = res.and_then(|r| r.checked_add(len));
            let (res, end) = match (res, end) {
                (Some(r), Some(e)) if e < data.len() && data[off] == 0xff && data[e] == 0xff => (r, e),
                _ => bail!("malformed .dmnd record {}", self.next),
            };
            out.seq.extend(data[res..end].iter().map(|&b| self.remap[(b & 0x7f) as usize]));
            out.off.push(out.seq.len());
            let t0 = end + 1;
            let mut t1 = t0;
            while t1 < data.len() && data[t1] != 0 {
                t1 += 1;
            }
            let title = String::from_utf8_lossy(&data[t0..t1]);
            out.ids.push(title.split_whitespace().next().unwrap_or("").to_string());
            self.next += 1;
            if out.seq.len() >= max_letters {
                break;
            }
        }
        Ok(!out.ids.is_empty())
    }
}

// ---------------------------------------------------------------- memory

/// In-memory reference list.
pub struct MemSource {
    seqs: Vec<(String, Vec<u8>)>,
    letters: u64,
    next: usize,
}

impl MemSource {
    /// From ASCII sequences.
    pub fn new(seqs: Vec<(String, Vec<u8>)>) -> MemSource {
        let seqs: Vec<(String, Vec<u8>)> = seqs.into_iter().map(|(id, s)| (id, super::alphabet::encode(&s))).collect();
        let letters = seqs.iter().map(|s| s.1.len() as u64).sum();
        MemSource { seqs, letters, next: 0 }
    }
}

impl SeqSource for MemSource {
    fn letters(&self) -> u64 {
        self.letters
    }
    fn sequences(&self) -> Option<u64> {
        Some(self.seqs.len() as u64)
    }
    fn next_block(&mut self, max_letters: usize, out: &mut Block) -> Result<bool> {
        out.clear();
        out.first = self.next as u64;
        while self.next < self.seqs.len() {
            let (id, s) = &self.seqs[self.next];
            out.ids.push(id.clone());
            out.seq.extend_from_slice(s);
            out.off.push(out.seq.len());
            self.next += 1;
            if out.seq.len() >= max_letters {
                break;
            }
        }
        Ok(!out.ids.is_empty())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fasta_blocks() {
        let dir = std::env::temp_dir().join(format!("ps_db_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("x.faa");
        std::fs::write(&p, ">a desc\nMKT\nAY\n>b\nGG\n>c\n\n>d\nWWWW\n").unwrap();
        let mut src = FastaSource::open(&p).unwrap();
        assert_eq!(src.letters(), 11);
        assert_eq!(src.sequences(), Some(4));
        let mut b = Block::default();
        assert!(src.next_block(3, &mut b).unwrap());
        assert_eq!(b.ids, vec!["a"]);
        assert_eq!(b.get(0).len(), 5);
        assert!(src.next_block(3, &mut b).unwrap());
        assert_eq!(b.ids, vec!["b", "c", "d"]);
        assert_eq!(b.get(1).len(), 0);
        assert_eq!(b.get(2).len(), 4);
        assert_eq!(b.first, 1);
        assert!(!src.next_block(3, &mut b).unwrap());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn rejects_text_before_header() {
        let dir = std::env::temp_dir().join(format!("ps_db2_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("x.faa");
        std::fs::write(&p, "MKT\n>a\nAA\n").unwrap();
        assert!(FastaSource::open(&p).is_err());
        std::fs::remove_dir_all(&dir).ok();
    }
}
