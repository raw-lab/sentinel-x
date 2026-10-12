//! Convert a DIAMOND `.dmnd` database (format version 3, DIAMOND 2.x) to a
//! protein FASTA, for tools that need one. Sentinel-X itself reads `.dmnd`
//! files in place.
//!
//! Layout (little endian): magic `u64` @0, `db_version` `u32` @12, sequence
//! count N `u64` @16, letters `u64` @24, position-array offset `u64` @32;
//! each record is `0xff <residues> 0xff <title> 0x00`; the position array
//! holds N+1 `(u64 offset, u32 length, u32 padding)` entries. Residues index
//! DIAMOND's alphabet `ARNDCQEGHILKMFPSTWYVBJZX*` with bit 7 marking a
//! soft-masked residue.

use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{bail, Context, Result};

const DMND_MAGIC: u64 = 0x024af8a415ee186d;
const AA_ALPHABET: &[u8; 25] = b"ARNDCQEGHILKMFPSTWYVBJZX*";

/// Summary of a converted database.
#[derive(Clone, Copy, Debug)]
pub struct DmndConversion {
    pub sequences: usize,
    pub letters: u64,
}

fn u64_at(b: &[u8], off: usize) -> Result<u64> {
    let end = off.checked_add(8).filter(|&e| e <= b.len()).context("truncated .dmnd")?;
    Ok(u64::from_le_bytes(b[off..end].try_into().unwrap()))
}

fn u32_at(b: &[u8], off: usize) -> Result<u32> {
    let end = off.checked_add(4).filter(|&e| e <= b.len()).context("truncated .dmnd")?;
    Ok(u32::from_le_bytes(b[off..end].try_into().unwrap()))
}

/// Convert `dmnd_path` to FASTA at `fasta_out` (written to a temporary file
/// and renamed when complete; 60 residues per line).
pub fn dmnd_to_fasta(dmnd_path: &Path, fasta_out: &Path) -> Result<DmndConversion> {
    let f = std::fs::File::open(dmnd_path).with_context(|| format!("opening {}", dmnd_path.display()))?;
    // SAFETY: read-only mapping of a file we do not modify.
    let map = unsafe { memmap2::Mmap::map(&f) }.with_context(|| format!("mapping {}", dmnd_path.display()))?;
    let data: &[u8] = &map;
    if data.len() < 96 || u64_at(data, 0)? != DMND_MAGIC {
        bail!("{} is not a DIAMOND database", dmnd_path.display());
    }
    let version = u32_at(data, 12)?;
    if version != 3 {
        bail!("{}: .dmnd db_version {version} is not supported (DIAMOND 2.x writes version 3)", dmnd_path.display());
    }
    let n = u64_at(data, 16)? as usize;
    let letters = u64_at(data, 24)?;
    let pos = u64_at(data, 32)? as usize;
    let end = n.checked_add(1).and_then(|k| k.checked_mul(16)).and_then(|b| b.checked_add(pos)).filter(|&e| e <= data.len());
    if end.is_none() {
        bail!("{}: position array out of bounds (truncated file?)", dmnd_path.display());
    }
    let part = crate::db::download::part_path(fasta_out);
    let mut w =
        BufWriter::with_capacity(1 << 20, std::fs::File::create(&part).with_context(|| format!("creating {}", part.display()))?);
    let mut total = 0u64;
    let mut line = Vec::with_capacity(61);
    for i in 0..n {
        let e = pos + 16 * i;
        let off = u64_at(data, e)? as usize;
        let len = u32_at(data, e + 8)? as usize;
        let res = off.checked_add(1).context("bad record offset")?;
        let stop = res.checked_add(len).filter(|&s| s < data.len()).with_context(|| format!("record {i} out of bounds"))?;
        if data[off] != 0xff || data[stop] != 0xff {
            bail!("{}: malformed record {i}", dmnd_path.display());
        }
        let t0 = stop + 1;
        let t1 = data[t0..].iter().position(|&b| b == 0).map(|p| t0 + p).context("unterminated title")?;
        w.write_all(b">")?;
        w.write_all(&data[t0..t1])?;
        w.write_all(b"\n")?;
        for chunk in data[res..stop].chunks(60) {
            line.clear();
            for &b in chunk {
                let idx = (b & 0x7f) as usize;
                line.push(*AA_ALPHABET.get(idx).unwrap_or(&b'X'));
            }
            line.push(b'\n');
            w.write_all(&line)?;
        }
        total += len as u64;
    }
    w.flush()?;
    drop(w);
    std::fs::rename(&part, fasta_out).with_context(|| format!("finalising {}", fasta_out.display()))?;
    if total != letters {
        eprintln!("[sentinel-x] note: {} residues written, header says {}", total, letters);
    }
    Ok(DmndConversion { sequences: n, letters: total })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a minimal v3 .dmnd with two sequences.
    fn tiny_dmnd() -> Vec<u8> {
        let seqs: [(&str, &[u8]); 2] = [("ref1~K00001 desc", b"MKV"), ("ref2~K00002", b"AW*")];
        let enc = |c: u8| AA_ALPHABET.iter().position(|&a| a == c).unwrap() as u8;
        let mut body = vec![0u8; 96];
        let mut pos = Vec::new();
        for (t, s) in seqs {
            pos.push((body.len() as u64, s.len() as u32));
            body.push(0xff);
            body.extend(s.iter().map(|&c| enc(c)));
            body.push(0xff);
            body.extend_from_slice(t.as_bytes());
            body.push(0);
        }
        let pos_off = body.len() as u64;
        for (o, l) in &pos {
            body.extend_from_slice(&o.to_le_bytes());
            body.extend_from_slice(&l.to_le_bytes());
            body.extend_from_slice(&0u32.to_le_bytes());
        }
        body.extend_from_slice(&(body.len() as u64).to_le_bytes());
        body.extend_from_slice(&[0u8; 8]);
        body[0..8].copy_from_slice(&DMND_MAGIC.to_le_bytes());
        body[12..16].copy_from_slice(&3u32.to_le_bytes());
        body[16..24].copy_from_slice(&2u64.to_le_bytes());
        body[24..32].copy_from_slice(&6u64.to_le_bytes());
        body[32..40].copy_from_slice(&pos_off.to_le_bytes());
        body
    }

    #[test]
    fn converts_and_rejects_garbage() {
        let dir = std::env::temp_dir().join(format!("sentdmnd_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let d = dir.join("x.dmnd");
        std::fs::write(&d, tiny_dmnd()).unwrap();
        let out = dir.join("x.faa");
        let c = dmnd_to_fasta(&d, &out).unwrap();
        assert_eq!(c.sequences, 2);
        assert_eq!(std::fs::read_to_string(&out).unwrap(), ">ref1~K00001 desc\nMKV\n>ref2~K00002\nAW*\n");
        let mut bad = tiny_dmnd();
        bad.truncate(120);
        std::fs::write(&d, bad).unwrap();
        assert!(dmnd_to_fasta(&d, &out).is_err());
        std::fs::remove_dir_all(&dir).ok();
    }
}
