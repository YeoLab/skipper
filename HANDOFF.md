# HANDOFF

## Current state (2026-11-14 05:55)
Full Skipper run (target `all`, branch `dev/brian/2025revision`, created 2026-11-14 from `dev/charlene/2025revision` @ 362dc27) on the 6 rbp-portal manifests
(`/tscc/projects/ps-yeolab5/rbp-portal/enco*/eclip/Skipper/2.0.0/v2/output/*.csv`), split into 8 runs
(manifest x cell line), running ONE AT A TIME in
`B=/tscc/nfs/home/bay001/projects/charlene_paper_20250521/permanent_data/skipper2025revision` (on ps-yeolab3).

| Order | Run | Protocol | Overdispersion | Dry-run jobs |
|---|---|---|---|---|
| 1 | encore2_FINAL_HepG2 | ENCODE4 (SE) | input | 552 |
| 2 | encore_edgecases_HepG2 | SE | input | 3267 |
| 3 | encore_HepG2 | SE | input | 3991 |
| 4 | encore_edgecases_K562 | SE | input | 6344 |
| 5 | encore_K562 | SE | input | 8154 |
| 6 | encode3_H_HepG2 | ENCODE3 (PE) | clip | 17416 |
| 7 | encore2_FINAL_K562 | SE | input | 18109 |
| 8 | encode3_K_K562 | PE | clip | 20289 |

- Driver: inline `setsid nohup bash -c` loop (log `$B/sequential.log`, command documented in `$B/README.md`):
  runs `$B/<run>/run_skipper.sh` in order, then `df -h /tscc/projects/ps-yeolab3`; STOPS if <5T available.
  After the pass, retries once every run whose `supervise.log` doesn't end in `done`.
- Status 00:55: run 1 encore2_FINAL_HepG2 failed (old Lustre quota + VEP bug) -> retry pass.
  Run 2 encore_edgecases_HepG2 stopped at 72% (2362/3267): GPU jobs stuck on decommissioned rtx3090 partition;
  fixed 0f8608c (-> hotel-v100, qos condo-gpu), 18 stuck jobs scancelled, snakemake SIGTERMed -> retry pass.
  Also had check_window_concordance (fixed 1e35f13) + VEP (fixed a078230) errors.
  Run 3 encore_HepG2: score_fa invalid qos fixed 0d4a4b5; then GPU rules failed 'LLVM ERROR ... seqpro' from shared
  numba cache (2288aed's NUMBA_DISABLE_CACHING isn't a numba 0.60 setting) -> per-job NUMBA_CACHE_DIR abcf926.
  Snakemake SIGTERMed 05:45, draining 3 running train_model jobs; then driver starts run 4; run 3 -> retry pass.
  Stale shared cache /tscc/lustre/ddn/scratch/bay001/seqpro_* now unused (delete after run 3 drains). 17T free.
- Each run: `$B/<run>/run_skipper.sh` (retries on "sbatch failed"; `-n` = dry run), `config.yaml`,
  `manifest.csv` (only for mixed-cell manifests, split by `Cells`), logs `supervise.log`, `snakemake.log`,
  `.snakemake/slurm_logs/`.
- Configs from `encode_configs/Skipper_pe_small_test_bay001.yaml` / `encore_configs/Skipper_se_small_test_bay001.yaml`;
  changed only WORKDIR, MANIFEST, cell-matched `gencode.v40.annotation.<cell>_totalrna.gt1.*`, TMPDIR -> `$B/tmp/`,
  VEP_CACHEDIR -> `$B/vep_cache/` (Her et al. 2026: cell-specific >1 TPM transcriptome).
- Monitoring: session cron every 4 h (expires 2026-11-20 or when session ends).

## History
- 15:47 first attempt: 8 concurrent runs in `/tscc/lustre/ddn/scratch/bay001/charlene2025revision`
  (repo `scratch/` symlink). All died ~17:10: Lustre user quota full (24T/24T). Partial outputs (613G), that folder,
  and the repo `scratch/` symlink deleted 2026-11-13.
- 18:05 first driver logged rc=0 for a failed run (`$(date)` reset `$?`); replaced 19:0x.

## Gotchas
- Skipper fixes go on `dev/brian/2025revision` (user, 2026-11-14). Fixes a078230..362dc27 were first pushed to
  `dev/charlene/2025revision` (left as-is); both branches are identical at 362dc27.
- Old shared numba caches in Lustre scratch deleted 2026-11-14.
- charlene_paper_20250521 GitHub repo (byee4) doesn't exist; commit 919dd41 is local-only.
- Lustre quota nearly full (~510G free after cleanup): rules hardcode `NUMBA_CACHE_DIR`/`MPLCONFIGDIR=/tscc/lustre/ddn/scratch/$USER`
  (variants_rbpnet.smk, train_rbpnet.smk, benchmark_ml.smk); conda/singularity prefixes on Lustre (read-only use OK).
  If ML rules fail on writes, free Lustre space or patch rules.
- ps-yeolab3 is shared (1.6P, 18T free at launch); other users' writes count toward the 5T stop.
- "ILP solver ... PULP_CBC_CMD" warning is harmless.
- Run configs/scripts committed in byee4/charlene_paper_20250521 (branch `skipper2025revision`); pipeline fixes here on `dev/brian/2025revision`.
