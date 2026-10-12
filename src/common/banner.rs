//! The SENTINEL-X banner — the wordmark printed on `sentinel-x --help`, at the
//! start of a run, and in the `doctor` report. A watchtower stands guard beside
//! the SENTINEL-X wordmark, scanning a double helix, in the same spirit as the
//! DYNAMITE banner it sits alongside.

/// Crate version, surfaced in the banner.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Full banner: wordmark + tagline footer.
pub fn banner() -> String {
    format!(
        "{ART}\nSentinel-X v{VERSION}  \u{00b7}  genome integrity for every domain of life\nviruses \u{00b7} bacteria \u{00b7} archaea \u{00b7} eukaryotes   \u{2014}   completeness \u{00b7} contamination \u{00b7} validation\nvalidators  checkv \u{00b7} checkm2 \u{00b7} busco \u{00b7} eukcc \u{00b7} viralverify \u{00b7} trace \u{00b7} chimera-scan \u{00b7} sage\nbuilt on    dynamite-orfs \u{00b7} Cleaver-X \u{00b7} hydra-mpp \u{00b7} rayon\n",
        ART = ART,
        VERSION = VERSION
    )
}

/// Just the wordmark art (no footer): a watchtower beside the SENTINEL-X figlet
/// with a double-helix readout.
pub const ART: &str = r#"            .-.
            |█|            ____  _____ _   _ _____ ___ _   _ _____ _             __  __
           _|=|_          / ___|| ____| \ | |_   _|_ _| \ | | ____| |            \ \/ /
          |     |         \___ \|  _| |  \| | | |  | ||  \| |  _| | |     _____   \  /
          | (O) |          ___) | |___| |\  | | |  | || |\  | |___| |___ |_____|  /  \
          |     |         |____/|_____|_| \_| |_| |___|_| \_|_____|_____|        /_/\_\
          | (O) |        ──────────────────────────────────────────────────────────────
          |_____|         guarding the genome, base by base
         /|=====|\        A=T  G≡C   ⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇
        /_|=====|_\       T=A  C≡G   ⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇⌇
          |_____|"#;

#[cfg(test)]
mod tests {
    #[test]
    fn banner_contains_wordmark_and_tagline() {
        let b = super::banner();
        // The figlet's distinctive S-row appears verbatim.
        assert!(b.contains("/ ___||"));
        assert!(b.contains("guarding the genome"));
        assert!(b.contains("completeness"));
        assert!(b.contains("eukcc"));
    }
}
