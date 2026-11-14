# CHANGELOG

## 2026-11-13
- `tools/check_window_concordance.R`: write an NA `odds_data.tsv` when a replicate pair has insufficient data
  (<4 concordance cells); it is a declared output, so the rule failed with MissingOutputException
  (DDX1_HepG2_4226, encore_edgecases_HepG2).
- `rules/variants_rbpnet.smk` (`download_vep_cache`): `wget` -> `curl -fLO`; the `ensemblorg/ensembl-vep:release_113.4`
  container has no `wget` (job exited 127 in the 2025 revision run). Runtime 1h -> 4h for the ~20 GB download + untar.
- `HANDOFF.md`: state of the 2025 revision full run (configs in
  `charlene_paper_20250521/permanent_data/skipper2025revision`).
