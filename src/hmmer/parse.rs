//! Streaming parser for HMMER3 profile files (`HMMER3/b` … `HMMER3/f`), the
//! format written by `hmmbuild` and shipped by Pfam, CheckV, BUSCO, EukCC and
//! viralVerify.
//!
//! Profiles are parsed one record at a time ([`HmmReader`]), so multi-gigabyte
//! libraries such as Pfam-A can be searched without loading them into memory.
//! Gzip-compressed files are detected by their magic bytes and decompressed on
//! the fly.
//!
//! Record layout:
//! ```text
//! HMMER3/f [3.4 | Aug 2023]
//! NAME  PF00001.25
//! ACC   …            DESC  …          LENG  M         ALPH  amino
//! GA    25.00 25.00;  TC …;  NC …;    (Pfam cutoff lines end in ';')
//! STATS LOCAL MSV       -9.9  0.71
//! STATS LOCAL VITERBI  -10.8  0.71
//! STATS LOCAL FORWARD   -4.2  0.71
//! HMM   A C D E …
//!       m->m m->i m->d i->m i->i d->m d->d
//!   COMPO <20 -ln p>
//!         <20 -ln p>      (node-0 insert emissions)
//!         <7 -ln p>       (node-0 transitions)
//!   1     <20 -ln p> <map> <cons> <rf> <mm> <cs>
//!         <20 -ln p>
//!         <7 -ln p>
//!   …
//! //
//! ```
//! Probabilities are stored as `-ln(p)`; `*` means `p = 0`.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};

use crate::hmmer::alphabet::K;
use crate::hmmer::model::{Hmm, Stats};

/// Open a (possibly gzip-compressed) HMM file for streaming.
pub fn open_hmm(path: &Path) -> Result<HmmReader<Box<dyn BufRead + Send>>> {
    let mut f = File::open(path).with_context(|| format!("opening HMM file {}", path.display()))?;
    let mut magic = [0u8; 2];
    let n = f.read(&mut magic)?;
    drop(f);
    let f = File::open(path)?;
    let reader: Box<dyn BufRead + Send> = if n == 2 && magic == [0x1f, 0x8b] {
        Box::new(BufReader::with_capacity(1 << 20, flate2::read::MultiGzDecoder::new(BufReader::new(f))))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, f))
    };
    Ok(HmmReader::new(reader).with_source(path.display().to_string()))
}

/// Parse every profile in a file (convenience for small files).
pub fn parse_hmm_path(path: &Path) -> Result<Vec<Hmm>> {
    open_hmm(path)?.collect()
}

/// Parse every profile in an in-memory string.
pub fn parse_hmm_str(text: &str) -> Result<Vec<Hmm>> {
    HmmReader::new(text.as_bytes()).collect()
}

/// Iterator over the profiles of a HMMER3 stream.
pub struct HmmReader<R: BufRead> {
    inner: R,
    line: String,
    lineno: usize,
    source: String,
    peeked: Option<String>,
    done: bool,
}

impl<R: BufRead> HmmReader<R> {
    pub fn new(inner: R) -> Self {
        HmmReader { inner, line: String::new(), lineno: 0, source: "<stream>".into(), peeked: None, done: false }
    }

    pub fn with_source(mut self, s: String) -> Self {
        self.source = s;
        self
    }

    fn next_line(&mut self) -> Result<Option<String>> {
        if let Some(l) = self.peeked.take() {
            return Ok(Some(l));
        }
        self.line.clear();
        let n = self.inner.read_line(&mut self.line).with_context(|| format!("reading {}", self.source))?;
        if n == 0 {
            return Ok(None);
        }
        self.lineno += 1;
        Ok(Some(self.line.trim_end_matches(['\n', '\r']).to_string()))
    }

    fn need_line(&mut self, what: &str) -> Result<String> {
        self.next_line()?.ok_or_else(|| anyhow!("{}: unexpected end of file ({what})", self.source))
    }

    fn err(&self, msg: impl std::fmt::Display) -> anyhow::Error {
        anyhow!("{} line {}: {msg}", self.source, self.lineno)
    }

    fn parse_record(&mut self) -> Result<Option<Hmm>> {
        // Skip blank lines before the magic.
        let first = loop {
            match self.next_line()? {
                None => return Ok(None),
                Some(l) if l.trim().is_empty() => continue,
                Some(l) => break l,
            }
        };
        if !first.starts_with("HMMER3") {
            if first.starts_with("HMMER2") {
                return Err(self.err("HMMER2 format is not supported; convert with `hmmconvert`"));
            }
            return Err(self.err(format!("expected 'HMMER3/x' header, found {first:?}")));
        }

        let mut name = String::new();
        let mut acc = None;
        let mut desc = None;
        let mut m = 0usize;
        let mut ga = None;
        let mut tc = None;
        let mut nc = None;
        let mut stats = Stats::default();
        let mut alph_ok = false;

        loop {
            let line = self.need_line("header")?;
            let trimmed = line.trim_start();
            if trimmed.starts_with("HMM ") || trimmed == "HMM" || trimmed.starts_with("HMM\t") {
                break;
            }
            let mut parts = trimmed.splitn(2, char::is_whitespace);
            let tag = parts.next().unwrap_or("");
            let rest = parts.next().unwrap_or("").trim();
            match tag {
                "NAME" => name = rest.to_string(),
                "ACC" => acc = Some(rest.to_string()),
                "DESC" => desc = Some(rest.to_string()),
                "LENG" => m = rest.parse().map_err(|_| self.err(format!("bad LENG value {rest:?}")))?,
                "ALPH" => {
                    if !rest.eq_ignore_ascii_case("amino") {
                        return Err(self.err(format!(
                            "profile '{name}' has alphabet '{rest}'; only protein (amino) profiles are supported"
                        )));
                    }
                    alph_ok = true;
                }
                "GA" => ga = Some(cutoff_pair(rest).ok_or_else(|| self.err("bad GA line"))?),
                "TC" => tc = Some(cutoff_pair(rest).ok_or_else(|| self.err("bad TC line"))?),
                "NC" => nc = Some(cutoff_pair(rest).ok_or_else(|| self.err("bad NC line"))?),
                "STATS" => {
                    let f: Vec<&str> = rest.split_whitespace().collect();
                    if f.len() >= 4 {
                        let a: f64 = f[2].parse().map_err(|_| self.err("bad STATS value"))?;
                        let b: f64 = f[3].parse().map_err(|_| self.err("bad STATS value"))?;
                        match f[1] {
                            "MSV" => stats.msv = Some((a, b)),
                            "VITERBI" => stats.viterbi = Some((a, b)),
                            "FORWARD" => stats.forward = Some((a, b)),
                            _ => {}
                        }
                    }
                }
                _ => {}
            }
        }
        if m == 0 {
            return Err(self.err(format!("profile '{name}' has no LENG")));
        }
        if !alph_ok {
            return Err(self.err(format!("profile '{name}' has no ALPH line")));
        }

        // Transition-order line.
        self.need_line("transition header")?;

        let mut mat = vec![[0f32; K]; m + 1];
        let mut t = vec![[0f32; 7]; m + 1];
        let mut compo = None;

        // Optional COMPO line.
        let l = self.need_line("model body")?;
        let body_first = if l.trim_start().starts_with("COMPO") {
            let mut c = [0f32; K];
            let toks: Vec<&str> = l.split_whitespace().skip(1).collect();
            fill_probs(&toks, &mut c).map_err(|e| self.err(format!("COMPO: {e}")))?;
            compo = Some(c);
            self.need_line("node-0 insert line")?
        } else {
            l
        };
        // node-0 insert emissions (not used by HMMER scoring) -> validate count.
        if body_first.split_whitespace().count() < K {
            return Err(self.err("truncated node-0 insert line"));
        }
        let tr0 = self.need_line("node-0 transitions")?;
        fill_trans(&tr0, &mut t[0]).map_err(|e| self.err(format!("node 0 transitions: {e}")))?;

        for k in 1..=m {
            let ml = self.need_line("match line")?;
            let toks: Vec<&str> = ml.split_whitespace().collect();
            let idx: usize = toks.first().and_then(|s| s.parse().ok()).ok_or_else(|| self.err("bad match-line node index"))?;
            if idx != k {
                return Err(self.err(format!("node index {idx}, expected {k}")));
            }
            if toks.len() < 1 + K {
                return Err(self.err(format!("truncated match emissions at node {k}")));
            }
            fill_probs(&toks[1..1 + K], &mut mat[k]).map_err(|e| self.err(format!("node {k}: {e}")))?;
            let il = self.need_line("insert line")?;
            if il.split_whitespace().count() < K {
                return Err(self.err(format!("truncated insert emissions at node {k}")));
            }
            let tl = self.need_line("transition line")?;
            fill_trans(&tl, &mut t[k]).map_err(|e| self.err(format!("node {k} transitions: {e}")))?;
        }

        // Record terminator.
        loop {
            match self.next_line()? {
                None => break,
                Some(l) if l.trim().is_empty() => continue,
                Some(l) if l.trim() == "//" => break,
                Some(l) => return Err(self.err(format!("expected '//', found {l:?}"))),
            }
        }

        Ok(Some(Hmm { name, acc, desc, m, mat, t, compo, ga, tc, nc, stats }))
    }
}

impl<R: BufRead> Iterator for HmmReader<R> {
    type Item = Result<Hmm>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        match self.parse_record() {
            Ok(Some(h)) => Some(Ok(h)),
            Ok(None) => {
                self.done = true;
                None
            }
            Err(e) => {
                self.done = true;
                Some(Err(e))
            }
        }
    }
}

/// Parse a `GA`/`TC`/`NC` value pair. Pfam writes a trailing `;`.
pub fn cutoff_pair(rest: &str) -> Option<(f32, f32)> {
    let cleaned = rest.replace(';', " ");
    let mut it = cleaned.split_whitespace();
    let a = it.next()?.parse::<f32>().ok()?;
    let b = it.next().and_then(|s| s.parse::<f32>().ok()).unwrap_or(a);
    Some((a, b))
}

fn tok_prob(tok: &str) -> Result<f32> {
    if tok == "*" {
        return Ok(0.0);
    }
    let v: f64 = tok.parse().map_err(|_| anyhow!("non-numeric probability token {tok:?}"))?;
    if !v.is_finite() || v < 0.0 {
        bail!("invalid -ln(p) value {tok:?}");
    }
    Ok((-v).exp() as f32)
}

fn fill_probs(toks: &[&str], out: &mut [f32; K]) -> Result<()> {
    if toks.len() < K {
        bail!("expected {K} values, found {}", toks.len());
    }
    for (o, t) in out.iter_mut().zip(toks.iter()) {
        *o = tok_prob(t)?;
    }
    Ok(())
}

fn fill_trans(line: &str, out: &mut [f32; 7]) -> Result<()> {
    let toks: Vec<&str> = line.split_whitespace().collect();
    if toks.len() < 7 {
        bail!("expected 7 transition values, found {}", toks.len());
    }
    for (o, t) in out.iter_mut().zip(toks.iter()) {
        *o = tok_prob(t)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    pub(crate) const TINY: &str = "HMMER3/f [3.4 | Aug 2023]
NAME  tiny
ACC   PF99999.1
DESC  a three node test model
LENG  3
ALPH  amino
GA    20.00 18.50;
NC    10.0 9.0;
STATS LOCAL MSV       -8.5000  0.71000
STATS LOCAL VITERBI   -9.0000  0.71000
STATS LOCAL FORWARD   -3.5000  0.71000
HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y
            m->m     m->i     m->d     i->m     i->i     d->m     d->d
  COMPO   2.6 4.2 2.9 2.7 3.2 2.7 3.8 2.8 2.8 2.3 3.7 3.2 3.0 3.2 2.9 2.7 2.9 2.7 4.5 3.5
          2.6 4.2 2.9 2.7 3.2 2.7 3.8 2.8 2.8 2.3 3.7 3.2 3.0 3.2 2.9 2.7 2.9 2.7 4.5 3.5
          0.02 4.5 4.5 0.6 0.8 0.0 *
      1   5 5 5 5 5 5 5 5 5 5 0.1 5 5 5 5 5 5 5 5 5 1 m - - -
          2.6 4.2 2.9 2.7 3.2 2.7 3.8 2.8 2.8 2.3 3.7 3.2 3.0 3.2 2.9 2.7 2.9 2.7 4.5 3.5
          0.02 4.5 4.5 0.6 0.8 0.6 0.8
      2   5 5 5 5 5 5 5 5 0.1 5 5 5 5 5 5 5 5 5 5 5 2 k - - -
          2.6 4.2 2.9 2.7 3.2 2.7 3.8 2.8 2.8 2.3 3.7 3.2 3.0 3.2 2.9 2.7 2.9 2.7 4.5 3.5
          0.02 4.5 4.5 0.6 0.8 0.6 0.8
      3   0.1 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 3 a - - -
          2.6 4.2 2.9 2.7 3.2 2.7 3.8 2.8 2.8 2.3 3.7 3.2 3.0 3.2 2.9 2.7 2.9 2.7 4.5 3.5
          0.0 * * 0.0 * * *
//
";

    #[test]
    fn parses_header_cutoffs_and_stats() {
        let h = parse_hmm_str(TINY).unwrap();
        assert_eq!(h.len(), 1);
        let h = &h[0];
        assert_eq!(h.name, "tiny");
        assert_eq!(h.acc.as_deref(), Some("PF99999.1"));
        assert_eq!(h.m, 3);
        assert_eq!(h.ga, Some((20.0, 18.5)));
        assert_eq!(h.nc, Some((10.0, 9.0)));
        assert!(h.tc.is_none());
        assert!(h.stats.is_calibrated());
        assert!((h.mat[1][10] - (-0.1f32).exp()).abs() < 1e-6);
        assert!(h.compo.is_some());
        assert_eq!(h.t[3][1], 0.0); // '*'
    }

    #[test]
    fn ga_line_with_trailing_semicolon() {
        assert_eq!(cutoff_pair("25.00 25.00;"), Some((25.0, 25.0)));
        assert_eq!(cutoff_pair("  792.1 640.3 ;"), Some((792.1, 640.3)));
        assert_eq!(cutoff_pair("30.0 30.0"), Some((30.0, 30.0)));
    }

    #[test]
    fn rejects_nucleotide_models_and_garbage() {
        let dna = TINY.replace("ALPH  amino", "ALPH  DNA");
        assert!(parse_hmm_str(&dna).is_err());
        let bad = TINY.replace("0.02 4.5 4.5 0.6 0.8 0.6 0.8", "0.02 x 4.5 0.6 0.8 0.6 0.8");
        assert!(parse_hmm_str(&bad).is_err());
    }

    #[test]
    fn parses_two_concatenated_records() {
        let two = format!("{TINY}{}", TINY.replace("NAME  tiny", "NAME  tiny2"));
        let h = parse_hmm_str(&two).unwrap();
        assert_eq!(h.len(), 2);
        assert_eq!(h[1].name, "tiny2");
    }
}
