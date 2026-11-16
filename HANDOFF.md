# HANDOFF

## Current state (2026-11-16 09:55) - BLOCKED
`B=/tscc/nfs/home/bay001/projects/charlene_paper_20250521/permanent_data/skipper2025revision`
Driver `$B/driver_20261116.sh` (copy in `2025revision_configs/`): encore2_FINAL_HepG2 -> encore_HepG2 -> encore_K562.
- encore2_FINAL_HepG2 done 07:36 (148/148). CDK9 972 windows, RACK1 346, PABPN1 0; enrichment r 0.20/0.07/0.12
  (all fail >0.5); top motifs UGAUGA-core (C/D-box-like), weak.
- encore_HepG2 stopped 09:48 at 817/819: only `all_ml_variants_output` (`touch`) left; it and encore_K562 fail with
  `OSError: [Errno 122] Disk quota exceeded`. Cause: **Lustre scratch user quota full** (`lfs quota -h -u bay001
  /tscc/lustre/ddn`: 24T/24T). Snakemake source cache (`~/.cache` -> `/tscc/lustre/ddn/scratch/bay001/.cache`),
  `--conda-prefix` and `--singularity-prefix` live there. Not cleaned (needs user OK).
- Next: user frees scratch; then `setsid nohup bash -c "$B/driver_20261116.sh >> $B/sequential.log 2>&1" &`
  (encore2_FINAL_HepG2 is a no-op). No monitoring cron (session ended).
- Tests done: `$B/encode3_test5_HepG2/` (853/853), `$B/encode3_testPUM_K562/` (346/346); configs in `2025revision_configs/`.
- Stopped runs (finished outputs preserved; killed jobs' outputs marked incomplete, redone on restart):

| Run | Last progress |
|---|---|
| 2 encore_edgecases_HepG2 | 2362/3267 |
| 4 encore_edgecases_K562 | 6123/6344 (221 left) |
| 5 encore_K562 | 2203 jobs left |
| 6-8 encode3_H_HepG2, encore2_FINAL_K562, encode3_K_K562 | not started |

  Restart any with `$B/<run>/run_skipper.sh > $B/<run>/supervise.log 2>&1` (one at a time unless the user says otherwise).

## Recent fixes
See CHANGELOG.md (python3essential unload, mtime rerun trigger, empty-window guards, train_model 6h).

## Test results (2026-11-15)
| RBP | repro windows | top HOMER de novo | RBPNet profile r (total>50) | sig GO terms |
|---|---|---|---|---|
| RBFOX2 | 27467 | TGCATG (matches) | 0.69 | 646 |
| PTBP1 | 14178 | TYYCTTTCY (CU-rich, matches) | 0.63 | 368 |
| QKI | 31380 | TACTAAC (matches) | 0.68 | 421 |
| TIA1 | 2681 | G-rich; TTTKGGGG #5 (weak U-rich) | 0.70 | 137 |
| U2AF2 | 58251 | TTYYTYCC (polypyrimidine, matches) | 0.76 | 253 |
| PUM1 (K562) | 24231 | GAAGAAGA, CDS-dominated; UGUAHAUA only #44 (1.7%) | 0.66 | 623 |
| PUM2 (K562) | 6739 | TGTAHATA (matches; 62% of sites, 3'UTR) | 0.55 | 70 |

- vs Her et al. 2026: RBPNet enrichment corr (dlogodds pearson, total>50; paper pass >0.5): RBFOX2 0.56, PTBP1 0.77,
  QKI 0.77, U2AF2 0.68, PUM2 0.69 pass; TIA1 0.50, PUM1 0.46 fail (paper: PUM1 ~0.4, 75% rRNA). Seqlets (97.5th pct):
  RBFOX2 UGCAUG, QKI UACUAAC, PTBP1 CU, U2AF2 polyU+AG|GU, PUM2 UGUAUAUA, TIA1 U-rich (+5'SS), PUM1 GA-rich.

## Next steps
- After the driver: ask the user about encore_edgecases_K562 (221 jobs), encore_edgecases_HepG2 (904) and the 3 unstarted runs.

## Gotchas
- charlene_paper_20250521 repo not used (user). All fixes go here on `dev/brian/2025revision`.
- `pkill -f run_skipper.sh` also matches the calling shell; kill by PID.
- Snakemake 9 SIGTERM drains running jobs; to stop fast, scancel jobs by run UUID (squeue job name) then SIGINT.
- GPU rules on hotel-v100; jobs pend on account-wide MaxNodePerAccount.
- ps-yeolab3 shared; stop if <5T free.
