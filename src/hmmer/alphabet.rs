//! Digital amino-acid alphabet, matching Easel's `eslAMINO` semantics.
//!
//! Codes `0..20` are the canonical residues in HMMER column order
//! (`ACDEFGHIKLMNPQRSTVWY`). Degenerate symbols get their own codes so they can
//! be scored the way HMMER scores them (background-weighted expected score for
//! match states, averaged odds for null2):
//!
//! | code | symbol | members |
//! |------|--------|---------|
//! | 20   | `B`    | D, N    |
//! | 21   | `J`    | I, L    |
//! | 22   | `Z`    | E, Q    |
//! | 23   | `X`    | all 20  |
//! | 24   | `*`    | nonresidue (stop / gap / missing) |
//!
//! Easel defines `U` (selenocysteine) and `O` (pyrrolysine) as one-member
//! degeneracies of `C` and `K`; scoring them is exactly equivalent to scoring
//! `C` / `K`, so they are mapped onto those codes at encoding time.

/// Number of canonical residues.
pub const K: usize = 20;
/// Number of digital codes (canonical + B, J, Z, X, `*`).
pub const NCODE: usize = 25;

pub const CODE_B: u8 = 20;
pub const CODE_J: u8 = 21;
pub const CODE_Z: u8 = 22;
pub const CODE_X: u8 = 23;
/// Nonresidue: `*` (stop codon), gap characters and anything unscorable.
pub const CODE_STAR: u8 = 24;

/// The canonical residues in HMMER column order.
pub const AA: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

/// HMMER3's default amino-acid background (`p7_AminoFrequencies`).
pub const BG: [f32; K] = [
    0.0787945, 0.0151600, 0.0535222, 0.0668298, 0.0397062, 0.0695071, 0.0229198, 0.0590092, 0.0594422, 0.0963728, 0.0237718,
    0.0414386, 0.0482904, 0.0395639, 0.0540978, 0.0683364, 0.0540687, 0.0673417, 0.0114135, 0.0304133,
];

/// Members of each degenerate code (canonical indices).
pub fn degenerate_members(code: u8) -> &'static [usize] {
    const B: [usize; 2] = [2, 11]; // D N
    const J: [usize; 2] = [7, 9]; // I L
    const Z: [usize; 2] = [3, 13]; // E Q
    const X: [usize; 20] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19];
    match code {
        CODE_B => &B,
        CODE_J => &J,
        CODE_Z => &Z,
        CODE_X => &X,
        _ => &[],
    }
}

/// Digitise one residue (case-insensitive).
#[inline]
pub fn encode_residue(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => 0,
        b'C' | b'U' => 1,
        b'D' => 2,
        b'E' => 3,
        b'F' => 4,
        b'G' => 5,
        b'H' => 6,
        b'I' => 7,
        b'K' | b'O' => 8,
        b'L' => 9,
        b'M' => 10,
        b'N' => 11,
        b'P' => 12,
        b'Q' => 13,
        b'R' => 14,
        b'S' => 15,
        b'T' => 16,
        b'V' => 17,
        b'W' => 18,
        b'Y' => 19,
        b'B' => CODE_B,
        b'J' => CODE_J,
        b'Z' => CODE_Z,
        b'X' => CODE_X,
        b'*' | b'-' | b'.' | b'_' | b'~' => CODE_STAR,
        c if c.is_ascii_alphabetic() => CODE_X,
        _ => CODE_STAR,
    }
}

/// Digitise a protein sequence. Whitespace is dropped.
pub fn encode(seq: &[u8]) -> Vec<u8> {
    seq.iter().filter(|b| !b.is_ascii_whitespace()).map(|&b| encode_residue(b)).collect()
}

/// Map a canonical column index back to its residue letter.
pub fn residue_letter(i: usize) -> u8 {
    AA.get(i).copied().unwrap_or(b'X')
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encodes_canonical_and_degenerate() {
        assert_eq!(encode(b"ACDEFGHIKLMNPQRSTVWY"), (0u8..20).collect::<Vec<_>>());
        assert_eq!(encode(b"uo"), vec![1, 8]);
        assert_eq!(encode(b"BJZX*"), vec![20, 21, 22, 23, 24]);
        assert_eq!(encode(b"m k\n"), vec![10, 8]);
    }

    #[test]
    fn background_sums_to_one() {
        let s: f32 = BG.iter().sum();
        assert!((s - 1.0).abs() < 1e-4);
    }
}
