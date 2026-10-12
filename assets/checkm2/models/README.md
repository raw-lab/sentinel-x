# Bundled CheckM2 models

Compressed copies of CheckM2 1.1.0's pre-trained models, embedded in the
binary when the `embed-models` feature is on (the default):

| file | CheckM2 original | used for |
|---|---|---|
| `general_model_COMP.gbm.z` | `general_model_COMP.gbm` (LightGBM) | general completeness |
| `model_CONT.gbm.z` | `model_CONT.gbm` (LightGBM) | contamination |
| `../specific_model.bin.z` | `specific_model_COMP.keras` + `scaler.sav` | specific (neural-network) completeness |
| `../ref_rsdata.bin.z` | `min_ref_rsdata_v1.npz` | cosine similarity for model selection |

The LightGBM files are zlib-compressed but otherwise unmodified. The Keras
network, the MinMaxScaler and the reference matrix were converted to a flat
little-endian layout by `scripts/validation/make_checkm2_assets.py`, with no
change to any value.

## Provenance

These models are the work of the CheckM2 authors, not of Sentinel-X. Please
cite CheckM2 when you use Sentinel-X's CheckM2 results:

> Chklovski, A., Parks, D.H., Woodcroft, B.J. & Tyson, G.W. (2023).
> *CheckM2: a rapid, scalable and accurate tool for assessing microbial genome
> quality using machine learning.* Nature Methods 20, 1203–1212.
> Source: https://github.com/chklovski/CheckM2

Sentinel-X is distributed under CC-BY-NC-4.0; see the `LICENSE` file at the
repository root.

Without `embed-models`, Sentinel-X has no specific model. It then reports the
general model only, from LightGBM files passed with `--checkm2-comp` and
`--checkm2-cont`.

## KO annotations

The models predict from CheckM2's 21,241-feature vector, and most of those
features are KO counts. Sentinel-X annotates KOs with its DIAMOND-compatible
search against CheckM2's own database (`sentinel-x download checkm2`), or
reads them from `--ko`. Without KO annotations, completeness and
contamination are reported as NA, and the metadata-only predictions are kept
in the extra column.
