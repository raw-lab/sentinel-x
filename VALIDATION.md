# Sentinel — validation & benchmarks against CheckV and CheckM2

All comparisons run the **original tools' own code** on the **same inputs**.
Harnesses live in the crate `examples/` directories and were re-run on this
build. Compilation: `cargo build --release --workspace` is clean (zero
warnings); `cargo test --workspace` passes **137/137**.

## 1. CheckV — terminal-repeat detection (57 cases, all identical)

### 1a. Real test set (46 contigs)
CheckV ships `test/test_sequences.fna` (46 contigs) with ground-truth calls in
`test/ground_truth/complete_genomes.tsv`. Three detectors on the identical 46
contigs:

* **Sentinel** — `sentinel_checkv::detect_terminal_repeat`
* **CheckV (live)** — CheckV's own `fetch_dtr` / `fetch_itr`, run from source
* **CheckV (ground truth)** — the shipped TSV

| contig | Sentinel | CheckV (live) | Ground truth |
|---|---|---|---|
| `…MeugFWYF_C6321` | DTR:82  | DTR:82  | DTR:82  |
| `…DCKB1_C1454`    | DTR:39  | DTR:39  | DTR:39  |
| `…DCKB1_C2382`    | DTR:518 | DTR:518 | DTR:518 |
| `…DCKB1_C3067`    | DTR:46  | DTR:46  | DTR:46  |
| `…DCKB1_C717`     | DTR:46  | DTR:46  | DTR:46  |

**Identical on all 46 contigs** (5 DTR + 41 no-repeat), and an exact match to
CheckV's ground-truth type and repeat length.

### 1b. Synthetic stress set (11 cases)
Constructed DTR/ITR/edge cases with known answers:

| case | Sentinel | CheckV |
|---|---|---|
| DTR 21/30/50/100/250 | DTR:{21,30,50,100,250} | identical |
| ITR 21/40/75 | ITR:{22,41,76} | identical |
| random (no repeat) | None | None |
| 15 bp / 20 bp (below 21 min) | None | None |

**Identical on all 11** — including CheckV's boundary-extension behaviour on
ITRs (the repeat extends one base past the seed) and correct rejection of
sub-threshold repeats.

## 2. CheckM2 — full feature vector (21,241 features, element-wise identical)

Same 126 predicted proteins + the same explicit KO set fed to CheckM2's own
`MetadataCalculator` + `KeggCalculator` (assembled in `feature_ordering.json`
order) and to Sentinel's `FeatureSpace::build_vector`:

| block | features | max \|diff\| |
|---|---|---|
| Metadata (AA counts, AALength, CDS) | 22 | 0.00e+00 |
| KO_Genes | 19,999 | 0.00e+00 |
| KO_Pathways | 416 | 0.00e+00 |
| KO_Modules | 757 | 0.00e+00 |
| KO_Categories | 47 | 0.00e+00 |
| **total** | **21,241** | **0.00e+00** |

Every feature matches; 0 of 21,241 differ.

## 3. CheckM2 — LightGBM prediction (6 diverse vectors, bit-identical)

Six feature vectors spanning 0 → 5000 KOs present (completeness 5 → 91%,
contamination 0 → 33%) scored by CheckM2's real models
(`general_model_COMP.gbm`, `model_CONT.gbm`) via Python `lightgbm` and by
Sentinel's `sentinel-gbm`:

| KOs present | CheckM2 comp | Sentinel comp | CheckM2 cont | Sentinel cont |
|---|---|---|---|---|
| 0    | 10.8193228677 | 10.8193228677 | 0.0031008532 | 0.0031008532 |
| 50   | 6.7079036131  | 6.7079036131  | 0.0 | 0.0 |
| 150  | 5.0507329672  | 5.0507329672  | 0.0 | 0.0 |
| 500  | 11.3158223976 | 11.3158223976 | 0.0 | 0.0 |
| 1500 | 18.9511046133 | 18.9511046133 | 0.0 | 0.0 |
| 5000 | 90.7587077065 | 90.7587077065 | 33.1833440251 | 33.1833440251 |

Worst deviation across all completeness + contamination: **0.00e+00**.

---

### Summary

| Path | Comparison | Outcome |
|---|---|---|
| CheckV terminal repeats | vs CheckV code + ground truth, 46 real + 11 synthetic | identical (57/57) |
| CheckM2 feature vector | vs `MetadataCalculator` + `KeggCalculator`, all blocks | identical (21,241/21,241) |
| CheckM2 LightGBM | vs real models, 6 vectors across the full range | bit-identical (0.00e+00) |

The reference databases (CheckV genome/HMM DB, CheckM2 DIAMOND DB, Pfam, BUSCO,
KOfam) are fetched on demand by `sentinel download <db>` and verified by
`sentinel doctor`. The comparisons above isolate the **algorithmic** equivalence
and therefore need no multi-gigabyte database; the database-dependent steps
(CheckV AAI/HMM search, CheckM2 DIAMOND KO annotation) are deterministic lookups
layered on top of the validated scoring code.
