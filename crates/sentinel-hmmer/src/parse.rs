//! Parser for HMMER3/f profile-HMM files (the text format emitted by
//! `hmmbuild` / shipped in Pfam-A.hmm, dbCAN, VOG, etc.).
//!
//! The format, per profile:
//! ```text
//! HMMER3/f [version]
//! NAME  <id>
//! ACC   <accession>
//! DESC  <description>
//! LENG  <M>
//! ALPH  amino
//! ... (GA/TC/NC/STATS lines) ...
//! HMM   A C D E F ...
//!       m->m m->i m->d i->m i->i d->m d->d
//!   COMPO <20 emission logs>
//!         <20 insert-0 emission logs>
//!         <7 node-0 transition logs>
//!   1     <20 match emission logs> <map> <cons> <rf> <mm> <cs>
//!         <20 insert emission logs>
//!         <7 transition logs>
//!   ...
//! //
//! ```
//! Emission and transition fields are `-ln(prob)`, with `*` meaning zero.

use anyhow::{anyhow, bail, Context, Result};
use std::io::BufRead;

use crate::model::{tok_log2, tok_prob, Hmm, Stats, Trans, K};

/// Parse every profile in a HMMER3 `.hmm` stream.
pub fn parse_hmm_file<R: BufRead>(mut reader: R) -> Result<Vec<Hmm>> {
    let mut text = String::new();
    reader
        .read_to_string(&mut text)
        .context("reading HMM file")?;
    parse_hmm_str(&text)
}

/// Parse every profile from an in-memory HMMER3 string.
pub fn parse_hmm_str(text: &str) -> Result<Vec<Hmm>> {
    let mut hmms = Vec::new();
    let mut lines = text.lines().peekable();

    while lines.peek().is_some() {
        // Skip blank lines between records.
        while let Some(l) = lines.peek() {
            if l.trim().is_empty() {
                lines.next();
            } else {
                break;
            }
        }
        if lines.peek().is_none() {
            break;
        }
        match parse_one(&mut lines)? {
            Some(h) => hmms.push(h),
            None => break,
        }
    }
    Ok(hmms)
}

fn two_floats(rest: &str) -> Option<(f64, f64)> {
    let mut it = rest.split_whitespace();
    let a = it.next()?.parse::<f64>().ok()?;
    let b = it.next()?.parse::<f64>().ok()?;
    Some((a, b))
}

fn parse_one<'a, I>(lines: &mut std::iter::Peekable<I>) -> Result<Option<Hmm>>
where
    I: Iterator<Item = &'a str>,
{
    let mut name = String::new();
    let mut acc = None;
    let mut desc = None;
    let mut m = 0usize;
    let mut ga = None;
    let mut tc = None;
    let mut nc = None;
    let mut stats = Stats::default();

    // ---- Header block, terminated by the `HMM` line. ----
    let first = match lines.next() {
        Some(l) => l,
        None => return Ok(None),
    };
    if !first.starts_with("HMMER3") {
        bail!("expected 'HMMER3/x' magic, found: {first:?}");
    }

    loop {
        let line = lines
            .next()
            .ok_or_else(|| anyhow!("unexpected EOF in HMM header"))?;
        let trimmed = line.trim_start();
        if trimmed.starts_with("HMM ") || trimmed == "HMM" {
            break;
        }
        let mut parts = trimmed.splitn(2, char::is_whitespace);
        let tag = parts.next().unwrap_or("");
        let rest = parts.next().unwrap_or("").trim();
        match tag {
            "NAME" => name = rest.to_string(),
            "ACC" => acc = Some(rest.to_string()),
            "DESC" => desc = Some(rest.to_string()),
            "LENG" => m = rest.parse().context("LENG")?,
            "GA" => ga = two_floats(rest),
            "TC" => tc = two_floats(rest),
            "NC" => nc = two_floats(rest),
            "STATS" => {
                // `STATS LOCAL MSV  -9.9 0.70`
                let mut it = rest.split_whitespace();
                let _scope = it.next();
                let kind = it.next().unwrap_or("");
                let mu = it.next().and_then(|s| s.parse::<f64>().ok());
                let lambda = it.next().and_then(|s| s.parse::<f64>().ok());
                if let (Some(mu), Some(lambda)) = (mu, lambda) {
                    match kind {
                        "MSV" => stats.msv = Some((mu, lambda)),
                        "VITERBI" => stats.viterbi = Some((mu, lambda)),
                        "FORWARD" => stats.forward = Some((mu, lambda)),
                        _ => {}
                    }
                }
            }
            _ => { /* ALPH, RF, MM, CONS, CS, DATE, … ignored */ }
        }
    }

    if m == 0 {
        bail!("HMM '{name}' has no LENG / zero length");
    }

    // The transition-order line after `HMM`.
    let _torder = lines
        .next()
        .ok_or_else(|| anyhow!("missing transition-order line"))?;

    // Optional COMPO line, then the node-0 insert + transition lines.
    let mut match_emit = vec![[0f64; K]; m + 1];
    let mut insert_emit = vec![[0f64; K]; m + 1];
    let mut trans = vec![
        Trans {
            mm: crate::score::NEG_INF,
            mi: crate::score::NEG_INF,
            md: crate::score::NEG_INF,
            im: crate::score::NEG_INF,
            ii: crate::score::NEG_INF,
            dm: crate::score::NEG_INF,
            dd: crate::score::NEG_INF,
        };
        m + 1
    ];

    // Peek: COMPO is optional.
    let peeked = lines
        .peek()
        .ok_or_else(|| anyhow!("unexpected EOF before model body"))?;
    if peeked.trim_start().starts_with("COMPO") {
        lines.next(); // consume COMPO (background; not needed for scoring)
    }
    // node-0 insert emissions
    let ins0 = lines
        .next()
        .ok_or_else(|| anyhow!("missing node-0 insert line"))?;
    parse_emission_into(ins0, &mut insert_emit[0])?;
    // node-0 transitions
    let tr0 = lines
        .next()
        .ok_or_else(|| anyhow!("missing node-0 transition line"))?;
    trans[0] = parse_trans(tr0)?;

    // ---- Per-node records 1..=M ----
    for k in 1..=m {
        let mrec = lines
            .next()
            .ok_or_else(|| anyhow!("missing match line for node {k}"))?;
        parse_match_into(mrec, &mut match_emit[k], k)?;

        let irec = lines
            .next()
            .ok_or_else(|| anyhow!("missing insert line for node {k}"))?;
        parse_emission_into(irec, &mut insert_emit[k])?;

        let trec = lines
            .next()
            .ok_or_else(|| anyhow!("missing transition line for node {k}"))?;
        trans[k] = parse_trans(trec)?;
    }

    // Consume the `//` record terminator if present.
    while let Some(l) = lines.peek() {
        if l.trim() == "//" {
            lines.next();
            break;
        }
        if l.trim().is_empty() {
            lines.next();
        } else {
            break;
        }
    }

    Ok(Some(Hmm {
        name,
        acc,
        desc,
        m,
        match_emit,
        insert_emit,
        trans,
        ga,
        tc,
        nc,
        stats,
    }))
}

/// A match line: `  k  e1 e2 .. e20  MAP CONS RF MM CS`. We take the first 20
/// numeric fields after the node index.
fn parse_match_into(line: &str, out: &mut [f64; K], k: usize) -> Result<()> {
    let mut it = line.split_whitespace();
    let idx = it
        .next()
        .ok_or_else(|| anyhow!("empty match line"))?
        .parse::<usize>()
        .with_context(|| format!("match node index (line: {line:?})"))?;
    if idx != k {
        bail!("match node index {idx} != expected {k}");
    }
    for a in out.iter_mut().take(K) {
        let tok = it
            .next()
            .ok_or_else(|| anyhow!("truncated match emissions at node {k}"))?;
        *a = tok_prob(tok);
    }
    Ok(())
}

/// An insert/COMPO-style line: 20 emission tokens (no leading index).
fn parse_emission_into(line: &str, out: &mut [f64; K]) -> Result<()> {
    let mut it = line.split_whitespace();
    // Some COMPO lines begin with the literal "COMPO"; that is handled by the
    // caller. Here we just read 20 tokens.
    for a in out.iter_mut().take(K) {
        let tok = it
            .next()
            .ok_or_else(|| anyhow!("truncated emission line: {line:?}"))?;
        *a = tok_prob(tok);
    }
    Ok(())
}

/// A transition line: 7 tokens `m->m m->i m->d i->m i->i d->m d->d`.
fn parse_trans(line: &str) -> Result<Trans> {
    let mut it = line.split_whitespace();
    let mut next = || -> Result<f32> {
        let tok = it
            .next()
            .ok_or_else(|| anyhow!("truncated transition line: {line:?}"))?;
        Ok(tok_log2(tok))
    };
    Ok(Trans {
        mm: next()?,
        mi: next()?,
        md: next()?,
        im: next()?,
        ii: next()?,
        dm: next()?,
        dd: next()?,
    })
}
