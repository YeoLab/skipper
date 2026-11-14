# CHANGELOG

## 2026-11-14
- `rules/variants_rbpnet.smk` (`download_vep_cache`): `curl -C -` (resume partial tarball) and runtime 4h -> 12h.
  Ensembl FTP gave ~1.5 MB/s; 23 GB tarball was 17 GB after 3 h in encore_HepG2 (job 14279323), would time out.
- `2025revision_configs/`: frozen `config.yaml`, `manifest.csv`, `run_skipper.sh` for the 8 revision runs + README
  (run order, sequential driver).
- ML rules (train_rbpnet.smk x6, benchmark_ml.smk x5, variants_rbpnet.smk, upload_rbpnet_to_hf.smk):
  `NUMBA_CACHE_DIR=$(mktemp -d)` per job instead of shared `/tscc/lustre/ddn/scratch/$USER`. numba 0.60 has no
  `NUMBA_DISABLE_CACHING`, so 2288aed didn't stop concurrent jobs corrupting the shared seqpro cache
  (`LLVM ERROR: Symbol not found: __gufunc__...seqpro...`; 17 train_model, 64 score_fa, 2 validation failures in encore_HepG2).
- `rules/variants_rbpnet.smk` (`score_fa`): GPU slurm resources moved into the rule; removed the profile's
  `set-resources: rbpnet_variants_score_fa` block. Snakemake stripped its outer quotes, so sbatch got
  `--qos=condo-gpu'` -> "Invalid qos specification" (encore_HepG2). Runtime now the rule's 20 min (profile had 10).
- GPU partition `rtx3090` (rules/train_rbpnet.smk x5, profiles/tscc2_snakemake9 `rbpnet_variants_score_fa`) and
  `rtx6000` (rules/benchmark_ml.smk) -> `hotel-v100`: rtx3090 nodes decommissioned 2026-11-05 (sinfo: drain,
  "PBB decommission"), rtx6000 down; GPU jobs pended forever. hotel-v100 accepts qos `condo-gpu` with account csd792
  (V100 nodes; test job ran).

## 2026-11-13
- `tools/check_window_concordance.R`: write an NA `odds_data.tsv` when a replicate pair has insufficient data
  (<4 concordance cells); it is a declared output, so the rule failed with MissingOutputException
  (DDX1_HepG2_4226, encore_edgecases_HepG2).
- `rules/variants_rbpnet.smk` (`download_vep_cache`): `wget` -> `curl -fLO`; the `ensemblorg/ensembl-vep:release_113.4`
  container has no `wget` (job exited 127 in the 2025 revision run). Runtime 1h -> 4h for the ~20 GB download + untar.
- `HANDOFF.md`: state of the 2025 revision full run (configs in
  `charlene_paper_20250521/permanent_data/skipper2025revision`).
