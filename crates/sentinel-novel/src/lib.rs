//! # sentinel-novel
//!
//! Three **reference-free** validators unique to Sentinel. Where CheckV,
//! CheckM2 and BUSCO depend on large reference databases, these modules judge a
//! contig purely from its own sequence and predicted genes, so they always run —
//! even with no database installed:
//!
//! * [`trace`] — **TRACE**: terminal-repeat & architecture completeness for
//!   viruses/phages.
//! * [`chimera`] — **CHIMERA-SCAN**: compositional-coherence contamination /
//!   chimera detection.
//! * [`sage`] — **SAGE**: single-copy anchor-gene copy-number and codon
//!   coherence for self-consistency-based contamination.

pub mod chimera;
pub mod sage;
pub mod trace;

pub use chimera::{scan as chimera_scan, ChimeraResult};
pub use sage::{sage as run_sage, SageResult};
pub use trace::{trace as run_trace, TraceConfidence, TraceResult};
