# Sentinel

```
            .-.
            |█|            ____  _____ _   _ _____ ___ _   _ _____ _
           _|=|_          / ___|| ____| \ | |_   _|_ _| \ | | ____| |
          |     |         \___ \|  _| |  \| | | |  | ||  \| |  _| | |
          | (O) |          ___) | |___| |\  | | |  | || |\  | |___| |___
          |     |         |____/|_____|_| \_| |_| |___|_| \_|_____|_____|
          | (O) |        ────────────────────────────────────────────────
          |_____|         guarding the genome, base by base
```

**Sentinel** (`sentinel-qc`) is a single pure-Rust binary for genome quality
control across viruses, bacteria, archaea and eukaryotes. It re-implements
the scoring of five standard tools, adds three reference-free checks, and
writes one table, one JSON and one HTML report.

* **Reference tools re-implemented:** CheckV, CheckM2, BUSCO, EukCC and viralVerify.
* **Reference-free checks:** TRACE, CHIMERA-SCAN and SAGE.
* **External programs needed at run time:** none. That covers Python, HMMER, DIAMOND, Prodigal and LightGBM.

| validator | reproduces | needs |
|---|---|---|
| **CheckV** | CheckV 1.1.1: AAI and HMM completeness, host contamination and proviruses, complete genomes (DTR/ITR), quality tiers, all CheckV tables | CheckV database |
| **CheckM2** | CheckM2 1.1.0: 21,241 KEGG features, LightGBM general models, neural-network specific model, model selection, `quality_report.tsv` | CheckM2 DIAMOND database (models are embedded) |
| **BUSCO** | BUSCO 6.1.0 scoring for odb10, odb12 and odb12.2 datasets: Single / Duplicated / Fragmented / Missing, `full_table.tsv` | a BUSCO lineage |
| **EukCC** | EukCC 2 completeness/contamination from a marker set (GA cutoffs) | a selected marker HMM set |
| **viralVerify** | viralVerify 1.1 virus / plasmid / chromosome Naive Bayes on Pfam domains (`--cut_nc`) | Pfam-A |
| TRACE | terminal repeats and coding architecture for phages | — |
| CHIMERA-SCAN | compositional junctions inside a contig (inter-species chimeras) | — |
| SAGE | contigs of a bin that depart from its core composition | — |

Built-in engines:

* **HMMER engine** with hmmsearch 3.4 scoring: length model, null2, domains, GA/NC/TC cutoffs, MSV/Viterbi/Forward filters.
* **DIAMOND-compatible protein search** (`protsearch`). Output is byte-identical to DIAMOND 2.1.9 with CheckM2's settings, and it reads `.dmnd` files in place.
* **LightGBM evaluator** for CheckM2's models.
* **[dynamite-orfs](https://crates.io/crates/dynamite-orfs)** gene calling with Prodigal and prodigal-gv presets.

How closely each validator matches the original tool, with numbers, is in
[AUDIT.md](AUDIT.md) and [VALIDATION.md](VALIDATION.md).

## Install

Needs Rust **1.75 or newer**. The apt `rustc`/`cargo` on Ubuntu 24.04 is enough.

```bash
cargo build --release --locked           # → target/release/sentinel (~21 MB, models embedded)
./install.sh                             # same, then installs to ~/.local/bin
PREFIX=/usr/local ./install.sh           # system-wide
```

Sentinel is one crate with one binary. Its features:

| feature | default | effect |
|---|---|---|
| `embed-models` | on | embed CheckM2's LightGBM and neural-network models (about 9 MB compressed) |
| `reads` | on | `--chunk-size`: stream huge FASTA files in chunks with [Cleaver-X](https://crates.io/crates/Cleaver-X) |
| `scale` | on | `cluster-status` for a [hydra-mpp](https://crates.io/crates/hydra-mpp) cluster |
| `rustyomestats` | off | cross-check assembly statistics with [rustyomestats](https://crates.io/crates/rustyomestats) on every run. Needs rustc ≥ 1.88 and pulls in polars, plotters and rust-bio (~220 crates) |

## Quick start

```bash
sentinel download checkv                       # viral reference database
sentinel download checkm2                      # CheckM2 DIAMOND database
sentinel download busco --lineage bacteria_odb12.2
sentinel download pfam                         # for --viralverify
sentinel doctor --verify                       # status, checksums, readiness

# viruses / phages: one row per contig and validator
sentinel run -i phages.fna -d virus -o viral.tsv --outdir viral_out -t 8

# bacterial or archaeal bins: a directory holds one genome per file
sentinel run -i bins/ -d bacteria -o bins.tsv --outdir bins_out -t 16

# add BUSCO
sentinel run -i genome.fna -d bacteria --lineage bacteria_odb12.2 -o qc.tsv

# eukaryotes
sentinel run -i genome.fna -d eukaryote --markers markers.hmm --lineage eukaryota_odb12 -o euk.tsv
```

### Domains

* `-d virus` and `-d mixed` work **per contig**:
  * CheckV, TRACE and CHIMERA-SCAN always run.
  * viralVerify runs with `--viralverify`. It is opt-in because a Pfam-A search is the heaviest step.
* `-d bacteria` and `-d archaea` work **per genome**:
  * CheckM2 and SAGE run, plus BUSCO with `--lineage`.
  * CHIMERA-SCAN runs per contig.
* `-d eukaryote` works **per genome**:
  * EukCC runs with `--markers`, BUSCO with `--lineage`, and SAGE.
  * CHIMERA-SCAN runs per contig.

`auto` is accepted as an alias of `mixed`; Sentinel does not detect the domain.

### Main options

```
-i, --input <FASTA|DIR>   FASTA (.fa/.fna[.gz], gzip detected by content); a directory = one genome per file
-d, --domain <D>          virus | bacteria | archaea | eukaryote | mixed       [mixed]
    --db <DIR>            database root            [--db > $SENTINEL_DB > ~/.sentinel/db]
-o, --output <TSV>        unified table ("-" = stdout)                          [stdout]
    --json <PATH>         the same table as JSON
    --html <PATH>         HTML report        [default: next to -o; never for stdout or /dev/*]
    --outdir <DIR>        each tool's own tables (CheckV, CheckM2, BUSCO, viralVerify, SAGE)
-t, --threads <N>         threads for every stage (0 = all CPUs)                [1]
    --min-len <N>         drop contigs shorter than N
    --no-novel            skip TRACE / CHIMERA-SCAN / SAGE
virus / mixed:   --proteins <FAA>  --checkv-strict  --viralverify  --chunk-size <BYTES>
prokaryotes:     --ttable <N>  --ko <FILE>  --checkm2-comp <GBM> --checkm2-cont <GBM>
genome-level:    --lineage <NAME|DIR>  --genes (protein FASTA input)
eukaryote:       --markers <HMM>
```

`sentinel help` prints everything; `sentinel info` lists the engines and features.

## Outputs

**Unified table** (`-o`): one row per contig or genome **per validator**.

```
contig_id  length  domain  validator  gene_count  completeness  contamination  quality  method  warnings  extra
```

* `quality` is CheckV's tier for viruses and the MIMAG tier (completeness and contamination) for CheckM2 and EukCC.
* `extra` holds `key=value` pairs separated by `;`. Inside values, `;` and `=` are written as `%3B` and `%3D`, and `%` is written as `%25` only when two hex digits follow it.
* For SAGE, `contamination` is the discordant bp % (see below).

**Each tool's own tables** (`--outdir`):

* `checkv/`: `quality_summary.tsv`, `completeness.tsv`, `contamination.tsv`, `complete_genomes.tsv`, `tmp/gene_features.tsv`, `tmp/aai.tsv`, `viruses.fna`, `proviruses.fna`.
* `checkm2/`: `quality_report.tsv` and `protein_files/`.
* `busco/<genome>/`: `full_table.tsv` and a short summary.
* `viralverify/`: `<input>_result_table.csv`.
* `sage/`: `<genome>.discordant_contigs.tsv`.

**HTML report**: one self-contained file with no external assets.

* **Summary panel:** completeness and contamination gauges for a single genome or contig. For many, the quality-tier distribution and median estimates.
* **Statistics strip:** size, contigs, N50, GC, coding density.
* **Results table:** sortable, one row per genome or contig, up to 5,000 rows.
* **Detail cards and markers:** per-validator cards for small inputs, plus marker bars for BUSCO and EukCC.
* **Warnings:** each one tagged with its validator and genome or contig.

### CheckV: default vs `--checkv-strict`

CheckV 1.1.1 groups DIAMOND's output rows by contig on the assumption that
they arrive in input order. DIAMOND ≥ 2.1 sorts them by query length under
CheckV's settings, so CheckV computes AAI completeness per fragment of each
contig and keeps the last fragment.

* **Default:** Sentinel uses all of a contig's hits.
* **`--checkv-strict`:** reproduces CheckV's tables exactly, bug included.

## Databases

| database | `sentinel download` | source | verified by |
|---|---|---|---|
| CheckV | `checkv` | `portal.nersc.gov/CheckV` (`CURRENT_RELEASE.txt`) | SHA-256 and size recorded (none published) |
| CheckM2 | `checkm2` | DOI `10.5281/zenodo.14897628` (DIAMOND DB v3, CheckM2 ≥ 1.1) | Zenodo MD5, `CONTENTS.json`, SHA-256 of the `.dmnd` |
| Pfam-A | `pfam` | EBI `current_release` | `md5_checksums` |
| BUSCO | `busco --lineage NAME` | `busco-data.ezlab.org/v5/data` (backup `busco-data2`) | MD5 from `file_versions.tsv` |

How installs work:

* Files download into `<root>/.staging/`, are verified and unpacked there, then moved into place atomically.
* A `.sentinel-complete.json` marker records each verified install.
* Failed transfers resume. A checksum failure triggers a fresh download.
* `sentinel download all [--lineage …]` fetches everything and skips what is already installed.
* DIAMOND databases are read in place; no conversion is needed. `sentinel convert` extracts a `.dmnd` to FASTA if you want one.
* For mirrors and offline tests, the URLs can be overridden with `SENTINEL_CHECKV_URL`, `SENTINEL_DOI_RESOLVER`, `SENTINEL_ZENODO_API`, `SENTINEL_PFAM_URL` and `SENTINEL_BUSCO_URL`. `SENTINEL_CHECKM2_SHA256` sets the expected `.dmnd` hash.

## The reference-free modules

These run without any database. They judge a sequence by its own
composition and genes, so they report flags and indices, not
reference-calibrated estimates.

* **TRACE** combines a CheckV-filtered terminal repeat (DTR/ITR), a head-to-tail k-mer overlap and the coding architecture (union coding density, largest gap). It returns a label: `complete`, `highly-complete`, `partial` or `fragment`.
* **CHIMERA-SCAN** tiles a contig into 5 kb windows and finds the best split into two segments.
  * It reports the junction strength D = 1 − cos(left mean profile, right mean profile) and flags D ≥ 0.12.
  * Calibration on complete chromosomes and inter-species junctions: 0.95 % false positives, 95.4 % sensitivity.
  * Phage–host junctions are usually invisible to composition; CheckV finds those.
* **SAGE** estimates the share of a bin's sequence that departs from its core tetranucleotide composition (contigs ≥ 5 kb), reported as `discordant_bp_pct`, and lists the discordant contigs.
  * On simulated MAGs it tracks the true contaminant fraction to within −1.2…+2.0 points.
  * It cannot see close relatives, and atypical plasmids or islands count as discordant.

## Other commands

```
sentinel blastp -q q.faa -d db.dmnd --preset checkm2|checkv [-t N]   # DIAMOND-compatible tabular output
sentinel convert db.dmnd out.faa                                       # .dmnd → FASTA
sentinel doctor [--db DIR] [--build] [--verify]
sentinel cluster-status <addr> [port]                                  # feature `scale`
sentinel info | version | banner | help
```

## Differences

The main ones:

* **Gene calls differ from the reference callers.**
  * dynamite-orfs matches Prodigal coordinates for 94–100 % of genes in single mode, 87 % in meta mode, and 81 % against pyrodigal-gv. It does not mask N runs.
  * These differences cause most residual end-to-end differences. `--proteins` and `--genes` accept external gene calls.
* **EukCC** does no tree placement or marker-set selection.
* **BUSCO genome mode** uses Sentinel's gene calls.
* **The parity tests used miniature databases**, because the real hosts were unreachable from the build environment.

## License

Creative Commons Attribution-NonCommercial (CC BY-NC 4.0)

See the LICENSE file for details.

## Credits

Sentinel re-implements published methods. Please cite the original tools:

* **CheckV:** Nayfach et al. 2021, *Nat. Biotechnol.*
* **CheckM2:** Chklovski et al. 2023, *Nat. Methods*
* **BUSCO:** Manni et al. 2021, *Mol. Biol. Evol.*
* **EukCC:** Saary et al. 2020, *Genome Biol.*
* **viralVerify:** Antipov et al. 2020, *Bioinformatics*
* **HMMER:** Eddy 2011
* **DIAMOND:** Buchfink et al. 2021
* **Prodigal:** Hyatt et al. 2010
* **MIMAG:** Bowers et al. 2017
