# Validation scripts

Helpers used for the head-to-head comparisons in VALIDATION.md. They were
written for the validation workspace and contain its absolute paths; adjust
them before reuse.

* `cmp_checkv.py A B` — compare two CheckV output folders table by table.
* `make_mini_checkv_db.py` — build a CheckV-layout mini database from the test set.
* `make_checkm2_db_realistic.py`, `run_checkm2_synthdb.py`, `make_checkm2_h2h.py`,
  `cmp_checkm2_vectors.py` — synthetic CheckM2 DIAMOND database, runs, vector comparison.
* `make_checkm2_assets.py` — converts CheckM2's Keras model, scaler and reference
  matrix into Sentinel's embedded `assets/checkm2/*.bin.z`.
* `make_mini_busco.py`, `make_mini_pfam.py`, `vv_classify_domtbl.py` — BUSCO/Pfam
  mini datasets and viralVerify's classifier on a given domain table.
* `build_mirror.py`, `mirror_server.py` — local mirror of the database
  upstreams for downloader tests (`SENTINEL_*_URL` overrides).
