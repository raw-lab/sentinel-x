# Calibration scripts (CHIMERA-SCAN, SAGE)

Python/NumPy re-implementations used to calibrate the reference-free modules
(see AUDIT.md §4.6). The Rust code reproduces their numbers exactly.

* `calib_lib.py` — FASTA reading, canonical tetranucleotide profiles.
* `calib3.py` — CHIMERA-SCAN: controls (window-aligned pieces of complete
  chromosomes) vs inter-genome junctions; rates for candidate rules.
* `viral.py` — the junction statistic on CheckV's test contigs.
* `sage_eval.py` — shows the old protein-composition index falling as
  contamination rises.
* `sage2.py`, `sage5.py` — contig-level discordance on simulated MAGs
  (single-linkage and the core-centroid estimator Sentinel-X uses).
* `sage_xcheck.py` — writes three simulated bins for a Rust/Python check.

Inputs: `TEST_BINS` = folder with `TEST1.fna`, `TEST2.fna`, `TEST3.fna`
(CheckM2's test genomes); `TEST_PROTEINS` = folder with `TEST1.ref.faa` …;
`CHECKV_TEST_FNA`, `CHECKV_QUALITY` = CheckV's test set and its
`quality_summary.tsv`. Needs NumPy.
