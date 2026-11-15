# HANDOFF

## Current state (2026-11-15 04:31)
Full Skipper run (target `all`, branch `dev/brian/2025revision`) on the 6 rbp-portal manifests
(`/tscc/projects/ps-yeolab5/rbp-portal/enco*/eclip/Skipper/2.0.0/v2/output/*.csv`), split into 8 runs
(manifest x cell line), ONE AT A TIME in
`B=/tscc/nfs/home/bay001/projects/charlene_paper_20250521/permanent_data/skipper2025revision`.

| Order | Run | Status |
|---|---|---|
| 1 | encore2_FINAL_HepG2 | failed -> retry pass |
| 2 | encore_edgecases_HepG2 | stopped 72% -> retry pass |
| 3 | encore_HepG2 | failed (numba cache, VEP timeout) -> retry pass |
| 4 | encore_edgecases_K562 | ended 22:06 rc=1 at 6123/6344 -> retry pass (dry run: 221 jobs) |
| 5 | encore_K562 | running since 22:06; 04:31: 3731/8153 (46%), no errors, no OOM/timeout retries |
| 6-8 | encode3_H_HepG2, encore2_FINAL_K562, encode3_K_K562 | queued |

- Driver: setsid bash loop (PID 929206, log `$B/sequential.log`): runs `$B/<run>/run_skipper.sh` in order,
  `df -h /tscc/projects/ps-yeolab3` after each, STOPS if <5T. Then retries once every run whose `supervise.log`
  doesn't end in `done`. 15T free at 04:31.
- Frozen run params: `2025revision_configs/<run>/` (keep in sync with live `$B/<run>/`).

## This session (00:31)
- Run 4 failures: consult_term_reference, run_homer, fetch_Clinvar/gnomAD (RPS10/RPS21 empty windows) -> fixed
  (944df60: bcftools guards use `zcat | wc -l`, run_homer skips on empty). 30 variant_analysis failures: python3essential
  module python shadowed the conda env (GLIBCXX error) -> run_skipper.sh unloads python3essential.
- run_skipper.sh adds `--rerun-triggers mtime` so guard-only rule edits don't rerun finished jobs.

## Open issues / next steps
- Run 5 started under the OLD env and rules: its variant_analysis (and any empty-window bcftools/homer jobs) will
  fail -> redone in the retry pass.
- From run 5: base mem ~2x measured peak, mem/runtime scale with attempt (`retries: 2`). Watch for OOM retries.
- GPU rules on hotel-v100 (a40 dropped, user 2026-11-14); jobs pend on account-wide MaxNodePerAccount.
- download_vep_cache resumable (curl -C -), 12h; partial tarball in `$B/vep_cache` kept.
- Delete stale `/tscc/lustre/ddn/scratch/bay001/seqpro_*` (unused).
- Monitoring: session cron every 4 h (expires 2026-11-20 or when session ends).

## Gotchas
- charlene_paper_20250521 repo not used (user). All fixes go here on `dev/brian/2025revision`.
- Run scripts must unload nodejs/galyleo/python3essential before `conda activate snakemake9`; jobs inherit the env.
- Conda python rules in runs 1-5 ran under python3essential's python (succeeded; not rerun).
- Lustre quota nearly full: ML rules write caches under `/tscc/lustre/ddn/scratch/$USER`.
- ps-yeolab3 is shared; other users' writes count toward the 5T stop.
- "ILP solver ... PULP_CBC_CMD" warning is harmless.
