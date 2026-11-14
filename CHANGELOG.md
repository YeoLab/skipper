# CHANGELOG

## 2026-11-13
- `rules/variants_rbpnet.smk` (`download_vep_cache`): `wget` -> `curl -fLO`; the `ensemblorg/ensembl-vep:release_113.4`
  container has no `wget` (job exited 127 in the 2025 revision run). Runtime 1h -> 4h for the ~20 GB download + untar.
- `HANDOFF.md`: state of the 2025 revision full run (configs in
  `charlene_paper_20250521/permanent_data/skipper2025revision`).
