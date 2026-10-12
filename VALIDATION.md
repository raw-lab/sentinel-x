# Sentinel 0.2.0 — validation against the original tools

Each section runs the original tool and Sentinel on the same inputs and
compares their native output files. The results are summarised in
[AUDIT.md §4](AUDIT.md#4-head-to-head-comparisons); this file records the
setup, the commands and the full tables.

## Environment

| tool | version | how it was run |
|---|---|---|
| CheckV | 1.1.1 (pip), with pyrodigal-gv 0.3.2 | `checkv end_to_end` |
| CheckM2 | 1.1.0 (source), LightGBM 4.7.0, Keras 3.8 on the JAX backend | `checkm2 predict` (`--genes` for protein input) |
| BUSCO | 6.1.0 (source) | `busco -m proteins --offline -l <lineage dir>` |
| viralVerify | 1.1 (source) | `viralverify -f contigs.fna --hmm Pfam-subset.hmm` |
| HMMER | 3.4 | `hmmsearch` (the reference for HMM scores) |
| DIAMOND | 2.1.9 (also 2.2.0) | `diamond blastp` with CheckV's and CheckM2's settings |
| Prodigal | 2.6.3 | as called by CheckM2 |

* Sentinel was built with `cargo build --release --locked` on rustc 1.97 and on rustc 1.75 (Ubuntu apt). Both builds produce the same results and pass 126/126 tests.
* Python 3.13 throughout.

### Databases

The sandbox could not reach the upstream servers (§7), so every comparison
uses a **small database built from public data**. Both tools always share
the same one:

* **CheckV:** `checkv-db-v1.5-mini`. It has CheckV's layout (`genome_db/`, `hmm_db/checkv_hmms/`): 3 HMM files, 69 reference genomes from CheckV's own test set, and a 900-cell error table.
* **CheckM2:** a synthetic `uniref100.KO.1.dmnd` (DIAMOND v3, 2.08 M letters). It holds KO-labelled proteins (`<id>~K#####`) from the test genomes plus decoys, so CheckM2's KO path has real hits.
* **Pfam:** a 548-family subset of Pfam-A with NC/GA/TC cutoffs.
* **BUSCO:** three 60-BUSCO lineages in odb10 (`lengths_cutoff` with σ), odb12 (no lengths file) and odb12.2 (median/scaledMAD) formats. Each has `dataset.cfg`, `hmms/` and `scores_cutoff`.

## 1. CheckV

* **Input:** CheckV's `test/test_sequences.fna`, 46 contigs with 36 viruses and 10 proviruses.
* **Comparison:** `scripts/validation/cmp_checkv.py`. Every row and column of the eight CheckV outputs is compared, first as text, then numerically at 1e-9.

### 1a. Strict mode on CheckV's own proteins

```bash
checkv end_to_end test_sequences.fna cv_out -d checkv-db-v1.5-mini
sentinel run -i test_sequences.fna -d virus --db sroot --proteins cv_out/tmp/proteins.faa \
     --checkv-strict --no-novel --outdir s_strict -o s_strict.tsv
python3 scripts/validation/cmp_checkv.py cv_out s_strict/checkv
```

| table | rows | identical | remaining differences |
|---|---|---|---|
| quality_summary.tsv | 46 | 46 | — |
| contamination.tsv | 46 | 46 | — |
| complete_genomes.tsv | 9 | 9 | — |
| tmp/aai.tsv | 150 | 150 | — |
| viruses.fna | 36 | identical | — |
| proviruses.fna | 10 | identical | — |
| completeness.tsv | 46 | 25 | 21 rows differ only in the last digit of `hmm_completeness_lower/upper` (≤ 1e-9) |
| tmp/gene_features.tsv | 1,177 | 1,116 | 61 rows: one printed E-value digit, or a score off by 0.1–0.2 at ≥ 250 bits |

* **The completeness.tsv digits.** CheckV sums HMM completeness over a Python `set`, so the last digit changes with the hash seed. Two CheckV runs with different seeds agree on 26/46 rows and differ in 20, the same kind of difference.
* **The gene_features.tsv scores.** These are float32 summation-order differences in Forward on very high-scoring hits. No cutoff decision changed.

### 1b. Default mode

In default mode Sentinel computes one AAI per contig. CheckV 1.1.1 instead
computes it per fragment, a side effect of the DIAMOND output order (see
AUDIT UP-01).

| run | High | Medium | Low | Not-determined |
|---|---|---|---|---|
| CheckV 1.1.1 | 5 | 28 | 9 | 4 |
| Sentinel `--checkv-strict` | 5 | 28 | 9 | 4 |
| Sentinel default, CheckV's proteins | 19 | 20 | 5 | 2 |
| Sentinel default, own gene calls | 19 | 20 | 5 | 2 |

The mini database contains the test genomes, so using all of a contig's hits
pushes many contigs to ~100 % AAI completeness. Expect a much smaller shift
on the real database.

### 1c. Gene calls

Sentinel's prodigal-gv preset (dynamite-orfs) against pyrodigal-gv:

* 1,328 vs 1,177 genes.
* 949 genes (80.6 %) have identical coordinates; 1,055 (89.6 %) share their stop codon.

## 2. CheckM2

**Inputs** (`bins/`):

* `TEST1`–`TEST3`: CheckM2's bundled test genomes.
* `mixed`: TEST2 + TEST3.
* `partial`: 42 × 50 kb fragments.
* `small`: an 80 kb fragment, which triggers meta mode.

### 2a. Same proteins

```bash
checkm2 predict -i prot_bins/ --genes -x faa -o cm2_genes --database_path db/.../uniref100.KO.1.dmnd
sentinel run -i prot_bins/ -d bacteria --genes --db db --outdir s_genes -o s_genes.tsv
diff cm2_genes/quality_report.tsv s_genes/checkm2/quality_report.tsv     # no output
```

* `quality_report.tsv` is **byte-identical** for 6/6 genomes.
* **Intermediate values** were checked against CheckM2's dumped arrays (`SENTINEL_CHECKM2_DEBUG=<dir>` writes Sentinel's vector and raw predictions):
  * The 21,241-feature vectors are bit-identical: 22 metadata, 19,999 KO counts, 416 pathways, 757 modules, 47 categories.
  * Both LightGBM predictions and the maximum cosine similarity are bit-identical.
  * The specific network's output is within 1–2 float32 ulps, and identical after CheckM2's rounding.

### 2b. Each tool calls its own genes

```bash
checkm2 predict -i bins/ -x fna -o cm2_out --database_path …
sentinel run -i bins/ -d bacteria --db db --outdir s_own -o s_own.tsv
```

| genome | CheckM2 comp / cont | Sentinel comp / cont | model | ttable | CDS (Prodigal / Sentinel) | identical gene coords |
|---|---|---|---|---|---|---|
| TEST1 | 87.76 / 10.12 | 87.75 / 10.12 | specific, both | 11 | 4876 / 4873 | 99.9 % |
| TEST2 | 88.10 / 0.23 | 88.12 / 0.23 | specific, both | 11 | 706 / 711 | 93.8 % |
| TEST3 | 93.26 / 0.77 | 93.27 / 0.76 | specific, both | 11 | 1222 / 1221 | 97.8 % |
| mixed | 100.00 / 23.09 | 100.00 / 23.09 | specific, both | 11 | 1933 / 1933 | 99.0 % |
| partial | 36.62 / 3.87 | 36.69 / 4.31 | specific, both | 11 | 2018 / 2020 | 99.9 % |
| small | 7.18 / 0.00 | 7.33 / 0.00 | specific, both | 11 | 85 / 85 | 87.1 % (meta) |

* Genome size, GC, contig N50, contig count and maximum contig length are identical for every genome.
* Coding density differs by ≤ 0.004.
* Stop codons agree for ≥ 99.9 % of genes in single mode and 97.6 % in meta mode.

## 3. BUSCO

```bash
busco -i genome.faa -m proteins -l lineages/mini_odb12.2 --offline -o busco_mini_odb12.2
sentinel run -i genome.faa -d bacteria --genes --lineage lineages/mini_odb12.2 --outdir s_mini_odb12.2
```

| lineage format | BUSCO summary | Sentinel summary | status + gene identical |
|---|---|---|---|
| odb10 | `C:83.3%[S:73.3%,D:10.0%],F:3.3%,M:13.3%,n:60` | identical | 60/60 |
| odb12 | `C:83.3%[S:73.3%,D:10.0%],F:3.3%,M:13.3%,n:60` | identical | 60/60 |
| odb12.2 | `C:83.3%[S:75.0%,D:8.3%],F:3.3%,M:13.3%,n:60` | identical | 60/60 |

In `full_table.tsv`, 4 scores (all above 1,650 bits) differ by 0.1–0.2 bits;
every length and status matches. The test proteome includes duplicated and
truncated copies, so the Duplicated and Fragmented rules are exercised.

## 4. viralVerify

The input is the same 46 contigs, searched against the 548-family Pfam
subset.

| comparison | identical predictions |
|---|---|
| viralVerify's classifier on Sentinel's proteins and domain table (`scripts/validation/vv_classify_domtbl.py`) vs Sentinel | 46/46, scores identical too |
| end to end, each tool calling its own genes | 45/46 |

The one end-to-end difference is CheckV's provirus `…1106445189670`:
viralVerify says Chromosome (−66.04), Sentinel says Uncertain (−5.12). The
two gene callers find different gene sets on this contig. The embedded
classifier table is byte-identical to viralVerify's
`classifier_table_uniq.txt`.

## 5. DIAMOND-compatible search (`sentinel blastp`)

Run with 2 threads; the outputs are compared with `cmp`.

| database / queries | settings | DIAMOND rows | identical rows | time DIAMOND / Sentinel |
|---|---|---|---|---|
| synthetic CheckM2 DB (2.08 M letters) / 4,876 proteins | CheckM2: `--max-target-seqs 1 --query-cover 80 --subject-cover 80 --id 30 --evalue 1e-05` | 2,490 | 2,490, file identical | 1.0 s / 1.6 s |
| simulated 227,425-sequence DB (75 MB) / 4,869 queries | CheckM2 | 4,658 | 4,658, file identical | 12.0 s / 15.8 s |
| simulated 227,425-sequence DB / 4,869 queries | CheckV: `--evalue 1e-5 --query-cover 50 --subject-cover 50 -k 10000` | 72,638 | 72,622 (99.98 %) | 12.3 s / 25.2 s |

Under CheckV's settings, 16 rows differ:

* **6 rows** are alignments of proteins longer than 1,000 aa. DIAMOND splits the seed diagonals into chained bands and stops early; Sentinel's band covers all of them and finds the higher-scoring alignment (for example 812 vs 803 bits).
* **10 rows** are borderline hits that Sentinel misses.

Sentinel also reports 28 hits that DIAMOND does not. On CheckV's test data
the search is exact: `tmp/aai.tsv` matches 150/150 rows.

## 6. Reference-free modules

These have no reference implementation, so they were calibrated with
simulations built from the three complete test genomes.

* **CHIMERA-SCAN** (5 kb windows, flag at D ≥ 0.12):
  * Whole chromosomes give D = 0.0008, 0.0354 and 0.0040. None is flagged; the old rules flagged all three.
  * 7,808 window-aligned pieces of 20 kb–1 Mb: 74 flagged (0.95 %).
  * 1,200 inter-species junctions with 10–980 kb on each side: 1,145 detected (95.4 %).
* **SAGE:**
  * Bins simulated by fragmenting the genomes into contigs averaging 8, 15 or 40 kb, then adding 5–45 % of contigs from another genome: the estimate lands within −1.2 … +2.0 points of the true contaminant bp fraction.
  * Fragmented single genomes score 0–2.7 %.
  * The `mixed` bin scores 39.61 %. Its 0.73 Mb TEST2 contig is 39.6 % of the bin.

The Rust code gives the same numbers as the Python calibration scripts on
spot checks:

* CHIMERA-SCAN: identical D values on all three chromosomes.
* SAGE: identical values on three simulated bins.

## 7. Database downloads

* **Local mirror test.** A local HTTP mirror reproduced each upstream's protocol: CheckV's `CURRENT_RELEASE.txt`, a DOI redirect, the Zenodo API with a `CONTENTS.json` backpack, Pfam's `md5_checksums`, and BUSCO's `file_versions.tsv`.
  * Sentinel pointed at it with the `SENTINEL_*_URL` variables installed all four databases and passed `doctor --verify`.
  * A second run skipped all four.
  * A wrong SHA-256, an MD5 mismatch and an unknown lineage each failed cleanly and installed nothing.
* **Real servers.** These were not reachable from the build sandbox: the proxy answered 403 for `portal.nersc.gov`, `doi.org`, `zenodo.org`, `ftp.ebi.ac.uk`, `busco-data.ezlab.org` and `busco-data2.ezlab.org`. Downloading the real databases still needs a networked machine:

  ```bash
  sentinel download all --lineage bacteria_odb12.2
  sentinel doctor --verify
  ```
