# HANDOFF

## Current state (2026-11-16 05:45)
8-run revision pass stopped 2026-11-15. **Restarted 2026-11-16 05:43** (user): `$B/driver_20261116.sh`
(copy in `2025revision_configs/`) runs encore2_FINAL_HepG2 (148 jobs) -> encore_HepG2 (819) -> encore_K562 (2203),
one at a time, 5T check before each. Others stay stopped. Monitoring: session cron every 4 h.
`B=/tscc/nfs/home/bay001/projects/charlene_paper_20250521/permanent_data/skipper2025revision`

- **Done: `$B/encode3_test5_HepG2/`** 08:42-11:41, 853/853, no errors, no OOM/timeout retries.
  5 HepG2 RBPs from `encode3_manifest_H_fastq.csv` (ENCODE3 PE, OVERDISPERSION_MODE clip, hepg2 >1 TPM transcriptome):
  RBFOX2_HepG2_ENCSR987FTF ("204"), PTBP1_HepG2_ENCSR384KAN, QKI_HepG2_ENCSR570WLM, TIA1_HepG2_ENCSR623VEQ,
  U2AF2_HepG2_ENCSR202BFN. Config = encode3_H_HepG2 with WORKDIR/MANIFEST changed. Dry run: 853 jobs.
  Frozen params: `2025revision_configs/encode3_test5_HepG2/`.
- **Done: `$B/encode3_testPUM_K562/`** (12:44-22:28, 346/346; PUM2 train_model timed out once at 3h, retry OK): PUM1_K562_ENCSR308YNT, PUM2_K562_ENCSR661ICQ
  (only in encode3 K manifest); config = encode3_K_K562 with WORKDIR/MANIFEST changed; dry run 346 jobs.
  Frozen params: `2025revision_configs/encode3_testPUM_K562/`. Expected motif UGUAHAUA (PUM).
- Stopped runs (finished outputs preserved; killed jobs' outputs marked incomplete, redone on restart):

| Run | Last progress |
|---|---|
| 1 encore2_FINAL_HepG2 | rerun started 05:37 (444 jobs), stopped |
| 2 encore_edgecases_HepG2 | 2362/3267 |
| 3 encore_HepG2 | rerun started 05:37 (1073 jobs), stopped |
| 4 encore_edgecases_K562 | 6123/6344 (221 left) |
| 5 encore_K562 | ~60% of 8153, stopped |
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
