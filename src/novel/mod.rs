//! # sentinel-novel
//!
//! Three **reference-free** validators unique to Sentinel. Where CheckV,
//! CheckM2 and BUSCO depend on large reference databases, these modules judge a
//! contig purely from its own sequence and predicted genes, so they always run —
//! even with no database installed:
//!
//! * [`trace`] — **TRACE**: terminal-repeat & architecture completeness for
//!   viruses/phages.
//! * [`chimera`] — **CHIMERA-SCAN**: compositional junctions inside a contig
//!   (inter-species chimeras).
//! * [`sage`] — **SAGE**: contigs of a genome bin that depart from its core
//!   tetranucleotide composition (reference-free contamination signal).

pub mod chimera;
pub mod sage;
pub mod trace;

pub use chimera::scan as chimera_scan;
pub use sage::sage as run_sage;
pub use trace::trace as run_trace;
