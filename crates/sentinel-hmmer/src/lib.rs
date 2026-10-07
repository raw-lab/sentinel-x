//! # sentinel-hmmer
//!
//! A pure-Rust Plan7 profile-HMM engine — a from-scratch reimplementation of
//! the parts of HMMER3 that Sentinel needs, with **no C, no `hmmsearch`
//! subprocess, and no FFI**. It provides:
//!
//! * a parser for HMMER3/f `.hmm` files ([`parse`]),
//! * the Plan7 model and its log-odds configuration ([`model`]),
//! * log-space **Viterbi** and **Forward** dynamic programming ([`dp`]),
//! * a thresholded search layer with E-values ([`search`]).
//!
//! All four quality-control validators in Sentinel (CheckV, CheckM2, BUSCO and
//! viralVerify) annotate predicted proteins by calling into this one engine,
//! so there is a single, auditable scoring path.
//!
//! ```no_run
//! use sentinel_hmmer::{parse_hmm_str, search_profile, Threshold};
//! # fn demo(hmm_text: &str) -> anyhow::Result<()> {
//! let hmms = parse_hmm_str(hmm_text)?;            // parse a .hmm file
//! let profile = hmms[0].configure(350);           // mean target length
//! let seqs = vec![("orf1".to_string(), b"MKAILV".to_vec())];
//! let hits = search_profile(&profile, &seqs, Threshold::Bits(20.0));
//! # let _ = hits; Ok(()) }
//! ```

pub mod dp;
pub mod model;
pub mod parse;
pub mod score;
pub mod search;

pub use dp::{score_seq, DpScore};
pub use model::{Hmm, Profile, Stats};
pub use parse::{parse_hmm_file, parse_hmm_str};
pub use search::{search_all, search_profile, Hit, Threshold};

/// Configure a batch of parsed HMMs for a given mean target length.
pub fn configure_all(hmms: &[Hmm], mean_len: usize) -> Vec<Profile> {
    hmms.iter().map(|h| h.configure(mean_len)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A tiny hand-written 3-node profile that strongly prefers the residues
    /// `M`, `K`, `A` at its three match positions. Built directly (rather than
    /// parsed) so the test is self-contained.
    fn toy_profile() -> Profile {
        use model::{Trans, K, NULL_FREQ};
        let m = 3;
        let mut match_emit = vec![[0f64; K]; m + 1];
        let insert_emit = vec![NULL_FREQ; m + 1];
        // Favoured residue gets 0.8, rest share 0.2.
        let favoured = [10usize /*M*/, 8 /*K*/, 0 /*A*/];
        for (k, &fav) in favoured.iter().enumerate() {
            let k = k + 1;
            for a in 0..K {
                match_emit[k][a] = 0.2 / (K as f64 - 1.0);
            }
            match_emit[k][fav] = 0.8;
        }
        // Strong M->M chain, weak indels.
        let strong = Trans {
            mm: prob(0.98),
            mi: prob(0.01),
            md: prob(0.01),
            im: prob(0.5),
            ii: prob(0.5),
            dm: prob(0.5),
            dd: prob(0.5),
        };
        let trans = vec![strong; m + 1];
        let hmm = Hmm {
            name: "TOY".into(),
            acc: None,
            desc: None,
            m,
            match_emit,
            insert_emit,
            trans,
            ga: None,
            tc: None,
            nc: None,
            stats: Stats::default(),
        };
        hmm.configure(50)
    }

    fn prob(p: f64) -> f32 {
        (p.log2()) as f32
    }

    #[test]
    fn cognate_scores_higher_than_random() {
        let p = toy_profile();
        let cognate = score_seq(&p, b"MKA");
        let random = score_seq(&p, b"WWW");
        assert!(
            cognate.forward_bits > random.forward_bits,
            "cognate {} should beat random {}",
            cognate.forward_bits,
            random.forward_bits
        );
        // Cognate should be net positive (better than null).
        assert!(cognate.forward_bits > 0.0);
        // Viterbi never exceeds Forward.
        assert!(cognate.viterbi_bits <= cognate.forward_bits + 1e-3);
    }

    #[test]
    fn finds_motif_embedded_in_flanks() {
        let p = toy_profile();
        let embedded = score_seq(&p, b"GGGGGMKAGGGGG");
        assert!(embedded.forward_bits > 0.0);
        assert!(embedded.env.is_some());
    }
}
