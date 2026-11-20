# CHANGELOG

## 2026-11-20
- `rules/bigwig.smk` (make_unscaled_bigwig): back to `bedtools genomecov -5` (read 5' ends = truncation sites) from
  `-split` (coverage, since f8f7c84). RBPNet `prep_data.py` trains on these bigwigs; the Her 2026 models used `-5`.
  Coverage-trained models disagree with the paper models (median variant δ r 0.43 on 28 ENCORE models); `-5` retrains
  match within retrain noise. YeoLab/skipper#69. Scaled bigwigs unchanged (coverage). Rebuild runs need their
  unscaled bigwigs and all RBPNet outputs regenerated.

## 2026-11-18
- `2025revision_configs/driver_20261118.sh`: sequential restart of encore_edgecases_K562 then encore_edgecases_HepG2
  (dry runs 221 and 904 jobs; 11/14 failures already fixed by 1e35f13, 6a156c4, 935c192, 944df60). encode3_H_HepG2 /
  encode3_K_K562 dry runs OK (17,415 / 20,288 jobs), not launched. Run records: igvf-proposal `modules/eclip/runs/`.

## 2026-11-16
- Scratch at 19.5T/24T after pytest_igvf removed; driver relaunched 11:23 (encore_HepG2 -> encore_K562).
- `2025revision_configs/driver_20261116.sh`: log run_skipper.sh's exit code (printed `$(date)`'s, always 0).
- encore_HepG2 (817/819) and encore_K562 stopped: Lustre scratch quota full (24T), snakemake source cache can't be written.
- `2025revision_configs/driver_20261116.sh`: sequential restart of encore2_FINAL_HepG2, encore_HepG2, encore_K562 (user).

## 2026-11-15
- encode3_testPUM_K562 finished (346/346): PUM2 top motif TGTAHATA (3'UTR); PUM1 top motif GAAGAAGA (CDS), PUM motif
  only rank 44. Results in HANDOFF.md.
- `rules/train_rbpnet.smk` (train_model): runtime 3h -> 6h base (x attempt). PUM2_K562 (582k training windows, ~3.4 min/epoch)
  hit the 3h limit at epoch ~44; train.py (max_epochs 100, patience 10) does not resume from checkpoints, so a
  timeout wastes the whole attempt. Earlier datasets trained in 25-66 min.
- `2025revision_configs/encode3_testPUM_K562/`: PUM1/PUM2 K562 test (encode3 K manifest subset; config = encode3_K_K562
  with WORKDIR/MANIFEST changed). Started 12:44.
- encode3_test5_HepG2 finished (853/853, no errors): de novo motifs match the known ones for RBFOX2, PTBP1, QKI and U2AF2;
  TIA1 top motif is G-rich (weak U-rich). Results table in HANDOFF.md.
- Stopped the 8-run revision driver and all runs (user); finished outputs kept.
- `2025revision_configs/encode3_test5_HepG2/`: 5-RBP HepG2 test (RBFOX2 ENCSR987FTF, PTBP1, QKI, TIA1, U2AF2)
  subset of `encode3_manifest_H_fastq.csv`; config = encode3_H_HepG2 with WORKDIR/MANIFEST changed. Started 08:42.
- Runs 1 (encore2_FINAL_HepG2) and 3 (encore_HepG2) restarted at 05:37, at the same time as run 5, at the user's request (dry runs: 444 and 1073 jobs).
- `2025revision_configs/*/run_skipper.sh` (and live `$B/<run>/`): also `module unload python3essential`. Its python
  (loaded in the driver's env) stayed ahead of the snakemake conda env, so `rbpnet_variants_variant_analysis` ran
  python3essential's scipy -> `GLIBCXX_3.4.30 not found` (30 jobs, encore_edgecases_K562). Other conda python rules
  in runs 1-5 ran under the same python (succeeded; not rerun).
- Same scripts: `--rerun-triggers mtime`. Guard-only rule edits otherwise re-trigger every finished job via code
  provenance (encore_edgecases_K562 retry: 3260 jobs incl. 1540 score_fa -> 221 = failed jobs + downstream).
- `rules/variants_rbpnet.smk` (fetch_SNP_from_gnomAD_and_roulette, fetch_Clinvar_SNP): empty-finemap guard
  `[ -s ]` -> `zcat | wc -l`; since 5d97e05 the empty finemapped bed.gz is a valid 20-byte gzip, so `-s` passed
  and `bcftools query -R` failed "Failed to read the regions" (RPS10/RPS21_K562).
- `rules/analysis.smk` (run_homer): skip with placeholder html/motif outputs when there are no finemapped windows
  (as on main / Skipper 2.0.0; zcat of the 0-byte sampled background failed).

## 2026-11-14
- All rules: `mem_mb` and `runtime` are `lambda wildcards, attempt: N * attempt`; profile adds `retries: 2` and
  attempt-scaled default-resources (`20000*attempt` MB, `30*attempt` min). A job that runs out of memory/time is
  resubmitted with 2x, then 3x.
- Base mem_mb lowered to ~2x measured peak RSS (sacct, all encore runs since 2026-11-13; rule from SubmitLine) where that
  is <75% of the old request; short jobs with unsampled RSS capped at 16 GB. Largest cuts: se_align_reads 160->67 GB,
  rbpnet_variants_fetch_variant_sequence 80->19, rbpnet_prepare_data 80->34, annotate_finemap 45->19, finemap_windows
  45->25, train_model 32->10, seqlet 32->9, score_fa 16->4, fastqc/trim 16->2. rbpnet_validation stays 32 GB (peak 18).
  PE-only and unrun rules (no measurements) keep their values.
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
