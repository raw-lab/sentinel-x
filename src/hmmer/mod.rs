//! # Sentinel HMMER engine
//!
//! A pure-Rust reimplementation of the parts of HMMER3 (3.4) that Sentinel's
//! validators need — no C, no FFI, no `hmmsearch` subprocess:
//!
//! * [`parse`] — streaming HMMER3 `.hmm` parser (plain or gzip);
//! * [`model`] — Plan7 model and HMMER's local/multihit profile configuration,
//!   null model and per-target length model;
//! * [`dp`] — MSV, Viterbi (score and traceback), Forward, Backward, posterior
//!   decoding and null2-by-expectation;
//! * [`pipeline`] — the hmmsearch pipeline: MSV → bias → Viterbi → Forward
//!   filters, domain definition, envelope rescoring and null2 correction;
//! * [`search`] — reporting thresholds (`-E`, `-T`, `--cut_ga|nc|tc`),
//!   sequence/domain E-values, and streaming search of HMM libraries.
//!
//! Bit scores are on hmmsearch's scale (validated against HMMER 3.4; see
//! `VALIDATION.md`), so GA/NC/TC cutoffs from Pfam, PANTHER/EukCC, CheckV and
//! BUSCO can be applied directly.

// Engine API: some accessors and record fields exist for completeness (domtblout fields, alternative thresholds) and are used by tests only.
#![allow(dead_code)]

pub mod alphabet;
pub mod dp;
pub mod model;
pub mod parse;
pub mod pipeline;
pub mod search;

pub use search::{search_hmm_file, search_hmm_file_models, Cutoff, SearchConfig, Targets};
