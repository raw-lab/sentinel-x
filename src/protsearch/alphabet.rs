//! Protein alphabet, DIAMOND's BLOSUM62 and its reduced seed alphabet.
//!
//! Residues use DIAMOND's encoding `ARNDCQEGHILKMFPSTWYVBJZX*` (codes
//! 0..24). `U`, `O` and `-` become `X` as in DIAMOND; any other byte also
//! maps to `X` (DIAMOND would reject the sequence).

/// DIAMOND's amino-acid alphabet (without its internal delimiter).
pub const ALPHABET: &[u8; 25] = b"ARNDCQEGHILKMFPSTWYVBJZX*";
/// Code of `X` (also DIAMOND's mask letter).
pub const X: u8 = 23;
/// Code of the stop symbol `*`.
pub const STOP: u8 = 24;
/// Number of codes.
pub const NCODES: usize = 25;
/// The 20 standard amino acids occupy codes 0..20.
pub const TRUE_AA: usize = 20;

const fn build_encoder() -> [u8; 256] {
    let mut t = [X; 256];
    let mut i = 0;
    while i < ALPHABET.len() {
        let c = ALPHABET[i];
        t[c as usize] = i as u8;
        if c.is_ascii_uppercase() {
            t[(c + 32) as usize] = i as u8;
        }
        i += 1;
    }
    t
}

/// ASCII -> code table.
pub static ENCODE: [u8; 256] = build_encoder();

/// Encode an ASCII protein sequence (whitespace dropped).
pub fn encode(seq: &[u8]) -> Vec<u8> {
    seq.iter().filter(|b| !b.is_ascii_whitespace()).map(|&b| ENCODE[b as usize]).collect()
}

/// Decode codes back to ASCII.
pub fn decode(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&c| ALPHABET[(c as usize).min(NCODES - 1)]).collect()
}

/// DIAMOND 2.1's BLOSUM62 (`stats/matrices/blosum62.cpp`), which differs
/// from NCBI's table in the ambiguity codes (X scores −1 against
/// everything, B/Z against N/D/Q/E use the newer NCBI values, J added).
#[rustfmt::skip]
pub const BLOSUM62: [[i8; 25]; 25] = [
    //A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   J   Z   X   *
    [ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1, -1, -1, -4], // A
    [-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1, -2,  0, -1, -4], // R
    [-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  4, -3,  0, -1, -4], // N
    [-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4, -3,  1, -1, -4], // D
    [ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -1, -3, -1, -4], // C
    [-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0, -2,  4, -1, -4], // Q
    [-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1, -3,  4, -1, -4], // E
    [ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -4, -2, -1, -4], // G
    [-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0, -3,  0, -1, -4], // H
    [-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3,  3, -3, -1, -4], // I
    [-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4,  3, -3, -1, -4], // L
    [-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0, -3,  1, -1, -4], // K
    [-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3,  2, -1, -1, -4], // M
    [-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3,  0, -3, -1, -4], // F
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -3, -1, -1, -4], // P
    [ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0, -2,  0, -1, -4], // S
    [ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1, -1, -1, -4], // T
    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -2, -2, -1, -4], // W
    [-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -1, -2, -1, -4], // Y
    [ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3,  2, -2, -1, -4], // V
    [-2, -1,  4,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4, -3,  0, -1, -4], // B
    [-1, -2, -3, -3, -1, -2, -3, -4, -3,  3,  3, -3,  2,  0, -3, -2, -1, -2, -1,  2, -3,  3, -3, -1, -4], // J
    [-1,  0,  0,  1, -3,  4,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -2, -2, -2,  0, -3,  4, -1, -4], // Z
    [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -4], // X
    [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1], // *
];

/// Score of codes `a` (query) and `b` (subject).
#[inline(always)]
pub fn score(a: u8, b: u8) -> i32 {
    BLOSUM62[(a as usize).min(NCODES - 1)][(b as usize).min(NCODES - 1)] as i32
}

/// BLOSUM62 background frequencies (DIAMOND `blosum62.background_freqs`,
/// the row sums of the BLOSUM62 joint probabilities).
pub const BACKGROUND: [f64; 20] = [
    0.07421620506799341,
    0.05161448614128464,
    0.044645808512757915,
    0.05362600083855441,
    0.02468745716794485,
    0.03425965059141602,
    0.0543119256845875,
    0.074146941452645,
    0.026212984805266227,
    0.06791736761895376,
    0.09890786849715096,
    0.05815568230307968,
    0.02499019757964311,
    0.04741845974228475,
    0.038538003320306206,
    0.05722902947649442,
    0.05089136455028703,
    0.013029956129972148,
    0.03228151231375858,
    0.07291909820561925,
];

/// Invalid reduced letter (breaks a seed).
pub const RED_INVALID: u8 = 0xFF;

/// DIAMOND 2.1's default reduction, Murphy et al.'s 10-letter alphabet
/// `[A] [KR] [EDNQ] [C] [G] [H] [ILVM] [FYW] [P] [ST]` (`basic/value.cpp`),
/// by code. B, J, Z, X and `*` cannot seed.
#[rustfmt::skip]
pub const REDUCED: [u8; 25] = [
    0, // A
    1, // R  (KR)
    2, // N  (EDNQ)
    2, // D
    3, // C
    2, // Q
    2, // E
    4, // G
    5, // H
    6, // I  (ILVM)
    6, // L
    1, // K
    6, // M
    7, // F  (FYW)
    8, // P
    9, // S  (ST)
    9, // T
    7, // W
    7, // Y
    6, // V
    RED_INVALID, RED_INVALID, RED_INVALID, RED_INVALID, RED_INVALID, // B J Z X *
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encoding_round_trip() {
        let e = encode(b"ARND cqeg*JUO-");
        assert_eq!(e[..4], [0, 1, 2, 3]);
        assert_eq!(e[4..8], [4, 5, 6, 7]);
        assert_eq!(e[8], STOP);
        assert_eq!(e[9], 21); // J
        assert_eq!(e[10..], [X, X, X]);
        assert_eq!(decode(&encode(b"MKTAYX*")), b"MKTAYX*".to_vec());
    }

    #[test]
    fn matrix_is_symmetric_and_matches_ncbi_core() {
        // NCBI BLOSUM62, residues ARNDCQEGHILKMFPSTWYV (lower triangle)
        const NCBI: &str = "4 -1 5 -2 0 6 -2 -2 1 6 0 -3 -3 -3 9 -1 1 0 0 -3 5 -1 0 0 2 -4 2 5 \
            0 -2 0 -1 -3 -2 -2 6 -2 0 1 -1 -3 0 0 -2 8 -1 -3 -3 -3 -1 -3 -3 -4 -3 4 \
            -1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 \
            -1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 -2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 \
            -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 \
            0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 \
            -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 \
            -2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 \
            0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4";
        let v: Vec<i8> = NCBI.split_whitespace().map(|t| t.parse().unwrap()).collect();
        assert_eq!(v.len(), 210);
        let mut k = 0;
        for i in 0..20 {
            for j in 0..=i {
                assert_eq!(BLOSUM62[i][j], v[k], "{i},{j}");
                k += 1;
            }
        }
        for i in 0..25 {
            for j in 0..25 {
                assert_eq!(BLOSUM62[i][j], BLOSUM62[j][i], "{i},{j}");
            }
        }
        let s: f64 = BACKGROUND.iter().sum();
        assert!((s - 1.0).abs() < 1e-9);
    }
}
