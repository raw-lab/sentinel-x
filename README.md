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
         /|=====|\        A=T  G≡C
        /_|=====|_\       T=A  C≡G
          |_____|
```

**Sentinel** is a single, (almost) pure-Rust tool that unifies four pillars of
genome quality control —
[**CheckV**](https://bitbucket.org/berkeleylab/checkv) (viruses),
[**CheckM2**](https://github.com/chklovski/CheckM2) (bacteria/archaea),
[**BUSCO**](https://busco.ezlab.org/) (eukaryotes/single-copy orthologs) and
[**ViralVerify**](https://github.com/ablab/viralVerify) (virus/plasmid/chromosome
classification) — behind one CLI, one gene call, and one report. It adds three
**novel reference-free validators** that work with *no database at all*.

It assesses **completeness**, **contamination**, and **structural validity** for
contigs and genomes from every domain of life.

---

## Why

The standard QC stack is four tools in three languages with four installation
procedures, four sets of dependencies, and four output formats. Sentinel folds
them into one statically-linked binary that:

* calls genes **once** (via the bundled `dynamite` engine) and shares the result
  across every validator;
* embeds a **pure-Rust Plan-7 HMMER** engine (`sentinel-hmmer`) and a
  **LightGBM evaluator** (`sentinel-gbm`) so inference needs no Python, no
  `hmmsearch`, and no `lightgbm` runtime;
* degrades gracefully — the parts that need gigabyte reference databases use them
  *when provided*, while the structural signals and all novel modules **always
  run**;
* parallelises with `rayon` locally and can scale across nodes with the bundled
  `hydra-mpp`, with optional GPU kernels via `wgpu`.

---

## Architecture

```
sentinel (CLI)
│
├── sentinel-core         FASTA I/O, Domain routing, gene calling, composition, report
├── sentinel-hmmer        pure-Rust Plan-7 profile HMM search (Forward/Viterbi, log-odds)
├── sentinel-gbm          LightGBM text-model evaluator (exact CheckM2 reproduction)
│
├── sentinel-checkv       viral: DTR/ITR terminal repeats, provirus boundary, completeness
├── sentinel-checkm2      bacterial/archaeal: 21,241-feature KEGG vector → LightGBM
├── sentinel-busco        single-copy orthologs: Single/Duplicated/Fragmented/Missing
├── sentinel-viralverify  Naive-Bayes virus / plasmid / chromosome from Pfam content
│
├── sentinel-novel        TRACE · CHIMERA-SCAN · SAGE  (reference-free, always on)
├── sentinel-gpu          rayon CPU kernels + optional wgpu compute (feature `gpu`)
└── sentinel-db           database registry + downloader + doctor (fetch/extract/verify)

vendored engines:  dynamite-core (ORF calling) · cleaver-core (FASTA/FASTQ/BAM)
                   hydra-mpp-core (distributed scaling) · rustyomestats
```

Every box above is its own crate in the Cargo workspace.

---

## What is validated

Sentinel's inference paths were cross-checked against the original
implementations:

| Component | Check | Result |
|---|---|---|
| `sentinel-gbm` | Predictions vs Python `lightgbm` on the **real** CheckM2 `general_model_COMP.gbm` and `model_CONT.gbm` | **exact to 8 decimals** |
| `sentinel-checkm2` | KO → pathway/module/category aggregation vs CheckM2's `KeggCalculator` | **exact on all 208 non-zero aggregates** |
| `sentinel-checkm2` | Feature-vector length & ordering | **21,241** features, matching CheckM2 |
| `sentinel-hmmer` | Bit-score convention (no-information sequence ≈ 0 bits) | matches HMMER log-odds model |
| `sentinel-checkv` | DTR/ITR detection logic | mirrors CheckV `fetch_dtr`/`fetch_itr` |

The full workspace test-suite (`cargo test --workspace`) passes with **0
failures**.

---

## Install / build

Requires a Rust toolchain (MSRV **1.75**).

```bash
cd sentinel
cargo build --release            # default features: reads + scale
./target/release/sentinel help   # run by path from the build dir
```

### Putting `sentinel` on your PATH

The quickest path is the bundled script:

```bash
./install.sh                     # builds + installs to ~/.local/bin
PREFIX=/usr/local ./install.sh   # system-wide (may need sudo)
```

Or do it manually. After building, the binary is at `target/release/sentinel`.
If you run `sentinel ...` and get `sentinel: command not found`, it just isn't
on your PATH yet. Pick one:

```bash
# (a) install with cargo  ->  ~/.cargo/bin/sentinel  (usually already on PATH)
cargo install --path crates/sentinel-cli

# (b) copy into a user bin dir on your PATH
install -m 0755 target/release/sentinel ~/.local/bin/sentinel

# (c) system-wide (needs sudo)
sudo install -m 0755 target/release/sentinel /usr/local/bin/sentinel
```

If you downloaded the prebuilt `sentinel` binary instead of building, make it
executable first, then do (b) or (c):

```bash
chmod +x sentinel
./sentinel help                  # or move it onto your PATH as above
```

Then `sentinel -h`, `sentinel doctor`, etc. work from anywhere.

Feature flags:

| Feature | Default | Effect |
|---|---|---|
| `reads` | on | transparent gzip input + large-file chunking (`cleaver`) |
| `scale` | on | rayon parallel contigs + `hydra-mpp` cluster status |
| `gpu`   | off | `wgpu` compute backend for pairwise-distance kernels |
| `polars`| off | (in `sentinel-cli`) polars dataframe export |

GPU and polars are off by default to keep the build light; the same APIs run on
CPU without them.

---

## Usage

```bash
sentinel run -i contigs.fna -d auto -o report.tsv
```

```
sentinel run OPTIONS:
  -i, --input <FASTA>        Input nucleotide FASTA (.fa/.fna[.gz])       [required]
  -d, --domain <D>           virus | bacteria | archaea | eukaryote | auto [default: auto]
  -o, --output <TSV>         Write the QC table here                      [default: stdout]
      --json <PATH>          Also write a JSON report
      --min-len <N>          Ignore contigs shorter than N bases          [default: 0]
  -t, --threads <N>          Worker threads                               [default: 1]
      --checkm2-comp <GBM>   CheckM2 completeness LightGBM model (.gbm)
      --checkm2-cont <GBM>   CheckM2 contamination LightGBM model (.gbm)
      --ko <FILE>            KOs present (one KO id per line) for CheckM2
      --gpu                  Use the GPU backend if available (else CPU)
      --no-novel             Skip the novel TRACE/CHIMERA/SAGE modules
      --chunk-size <BYTES>   Pre-split huge inputs (requires `reads`)
      --scale                Parallelise contigs across threads
```

### Examples

Viral contigs (structural CheckV + novel modules, no database needed):

```bash
sentinel run -i phage_contigs.fna -d virus -o viral_qc.tsv --json viral_qc.json
```

Bacterial genome with the real CheckM2 models:

```bash
sentinel run -i genome.fna -d bacteria \
  --checkm2-comp models/general_model_COMP.gbm \
  --checkm2-cont models/model_CONT.gbm \
  --ko genome_kofam_hits.txt \
  -o checkm2_qc.tsv
```

Each contig/genome yields one row **per validator**, e.g.:

```
contig_id          length  domain  validator      gene_count  completeness  contamination  quality          method                     extra
phage_1            23404   virus   checkv          50          NA            0.00           Not-determined   structural                 viral_genes=0;host_genes=0
phage_1            23404   virus   trace           50          NA            NA             High-quality     architecture               trace_confidence=highly-complete;coding_density=0.905
phage_1            23404   virus   chimera-scan    50          NA            0.00           Not-determined   tetranucleotide-coherence  n_windows=4;max_zscore=1.72;gc=0.500
```

---

## The novel modules

These are unique to Sentinel and run **database-free**, judging a sequence purely
from itself and its predicted genes.

### TRACE — Terminal-Repeat & Architecture Completeness
A completeness estimate for viruses/phages that fuses (1) direct/inverted
terminal repeats, (2) head-to-tail k-mer overlap signalling a circular/headful
molecule, and (3) coding architecture (gene density, largest non-coding gap).
Reports `complete` / `highly-complete` / `partial` / `fragment`.

### CHIMERA-SCAN — Compositional-Coherence Contamination
Tiles a contig into windows, computes each window's tetranucleotide signature,
and finds the sharpest composition seam plus the overall dispersion. Flags
cross-domain contamination and mis-assembled chimeras and reports a putative
breakpoint — including the clean two-organism case that defeats simple
centroid-dispersion methods.

### SAGE — Single-copy Anchor-Gene & Codon-Coherence
Asks whether a genome's proteins look like they came from **one** organism: it
measures amino-acid composition coherence around the proteome centroid, flags
compositional-outlier proteins (a database-free contamination proxy), and detects
near-identical single-copy violations.

---

## Databases (managed by `sentinel doctor` / `sentinel download`)

The reference validators need large databases that Sentinel does **not** bundle.
Sentinel manages them under a single root (`$SENTINEL_DB`, default
`~/.sentinel/db`):

```bash
sentinel doctor --build           # create the folder skeleton + manifest.json
sentinel download all             # fetch every reference database at once
sentinel download checkv          # …or fetch one DB:  CheckV genome/HMM DB
sentinel download busco --lineage bacteria_odb10
sentinel doctor                   # verify what's installed + per-validator readiness
```

`sentinel download all` fetches every database in one pass (skipping any already
installed, and continuing past individual failures); BUSCO is included only when
`--lineage` is supplied, since it cannot be fetched without one.

`sentinel doctor` reports the runtime, the database downloader, optional external
tools (DIAMOND, HMMER), and the status of every database (`OK` / `PARTIAL` /
`EMPTY` / `MISSING`), then summarises which validators are ready. The downloader
fetches from the **canonical upstream sources** — CheckV from
`portal.nersc.gov/CheckV`, CheckM2 from its Zenodo record
(`10.5281/zenodo.14897628`), Pfam-A from EBI, BUSCO from `busco-data.ezlab.org` —
extracts the archives (pure-Rust `tar`+`flate2`) and verifies checksums
(`sha2`). Layout, e.g. for CheckV, mirrors the original tool exactly
(`<root>/checkv/<release>/genome_db/…`, `…/hmm_db/…`) so the files are
drop-in.

A run can point at the root and auto-discover models:

```bash
sentinel run -i genome.fna -d bacteria --db ~/.sentinel/db -o qc.tsv
```

| Validator | Needs | Without it |
|---|---|---|
| CheckV (reference completeness) | CheckV HMM DB + reference-length table | terminal-repeat & structural verdict still produced |
| CheckM2 (prediction) | KOfam/DIAMOND DB (KO annotations) + the `.gbm` models | feature vector built; prediction only with models + `--ko` |
| BUSCO | a lineage `odb` dataset | SAGE + composition reported instead |
| ViralVerify | Pfam-A HMM DB (per-contig hits) | classifier table embedded; needs Pfam hits to classify |

What **is** embedded: the ViralVerify Pfam classification table and CheckM2's KEGG
feature-definition JSONs (feature ordering, module definitions, path/category
mapping).

## Validation against CheckV & CheckM2

Sentinel's inference was checked against the original tools' own code on the same
inputs (full details in `VALIDATION.md`):

* **CheckV terminal repeats** — identical to CheckV (and its ground truth) on all
  46 test contigs (DTR type + length).
* **CheckM2 metadata** — all 22 amino-acid/length/CDS features identical on the
  same 126 proteins.
* **CheckM2 KEGG aggregation** — exact on all 208 pathway/module/category
  aggregates.
* **CheckM2 LightGBM** — completeness and contamination bit-for-bit identical
  (diff `0.00e+00`) on the same 21,241-dimension vector.


---

## Status

* 11 Sentinel crates, ~5,500 lines of new Rust, plus four vendored Rust engines.
* `cargo test --workspace` — **all green**.
* Inference cross-validated against CheckM2's Python/LightGBM stack to numerical
  equality.

Licensed GPL-3.0-or-later (inherits the licenses of the wrapped tools).
```
