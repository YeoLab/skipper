# CHANGELOG

## 2026-11-14
- GPU memory requests lowered: `rules/train_rbpnet.smk` (x5) mem_mb 160000 -> 32000 (23 completed train_model jobs
  peaked at 1.9-4.7 GB RSS); `rules/variants_rbpnet.smk` score_fa 80000 -> 16000 (178 clean exits peaked at 1.6 GB;
  inputs <=180 MB FASTA scored in batches of 512). hotel-v100 bills mem at 0.2 SU/GB-h, so train_model 62 -> ~37 SU/h,
  score_fa 47 -> ~34 SU/h; smaller requests also let several 1-GPU jobs share a 4-GPU node under MaxNodePerAccount.
- `tools/consult_encode_reference.R`: empty-input guards. Placeholder TSV/PDF when no experiment has reproducible
  windows or REs (as main's consult_encode_reference_windows.R); empty feature table when only one of windows/REs is
  empty (simplify_*/inner_join failed on 0-column tibbles). Tested real (identical to encore_HepG2 output),
  windows-only, RE-only, none.
- Reviewed all tools/ against origin/main (Skipper 2.0.0); ported the robustness fixes this branch is exposed to:
  `tools/finemap_enriched_windows.R` writes a valid empty gzip (0-byte .gz broke downstream `zcat`, e.g. run_homer);
  `tools/annotate_finemapped_regions.py` catches EmptyDataError on that file and writes an empty output;
  `tools/identify_reproducible_windows.R` no-shared-windows case writes the full 30-column header
  (input_sum..q_min were missing; RPS10/RPS21_K562 in encore_edgecases_K562).
  Not ported (change results or this branch's rule I/O): call_enriched_windows threshold/normalization/blacklist
  changes, identify_reproducible_* reproducibility filter + na.rm, parse_gff base_feature_id column, output path layout.
- Empty-input fallbacks aligned with main (Skipper 2.0.0):
  `tools/consult_term_reference.R` empty case now uses main's 15-column `enriched_terms` schema + "0 or 1 significant terms" PDF;
  `tools/identify_reproducible_re.R` writes an empty `reproducible_enriched_re.tsv.gz` when no replicate has enriched REs
  (group_by on 0 columns errored; main has an equivalent fallback). Main's NA dummy row not used: this branch's
  consult_encode_reference.R has no dummy filter. Already guarded on this branch: call_enriched_windows.R,
  identify_reproducible_windows.R, sample_matched_background_by_region.R, finemap_enriched_windows.R,
  annotate_finemapped_regions (rule-level `-s` check).
- `tools/consult_term_reference.R`: write empty `enriched_terms.tsv.gz` + placeholder PDF when there are 0 reproducible
  enriched windows (summarize(table()) errored; RPS21_K562_9990, RPS10_K562_4104 in encore_edgecases_K562).
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
