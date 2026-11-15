# HANDOFF

## Current state (2026-11-15 11:45)
8-run revision pass STOPPED by the user (2026-11-15 ~08:20). Replaced by a 5-RBP test run.
`B=/tscc/nfs/home/bay001/projects/charlene_paper_20250521/permanent_data/skipper2025revision`

- **Done: `$B/encode3_test5_HepG2/`** 08:42-11:41, 853/853, no errors, no OOM/timeout retries.
  5 HepG2 RBPs from `encode3_manifest_H_fastq.csv` (ENCODE3 PE, OVERDISPERSION_MODE clip, hepg2 >1 TPM transcriptome):
  RBFOX2_HepG2_ENCSR987FTF ("204"), PTBP1_HepG2_ENCSR384KAN, QKI_HepG2_ENCSR570WLM, TIA1_HepG2_ENCSR623VEQ,
  U2AF2_HepG2_ENCSR202BFN. Config = encode3_H_HepG2 with WORKDIR/MANIFEST changed. Dry run: 853 jobs.
  Frozen params: `2025revision_configs/encode3_test5_HepG2/`.
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
  The old driver (PID 929206) is killed; don't restart it without asking.

## Recent fixes
- run_skipper.sh: unloads nodejs/galyleo/python3essential (module python shadowed conda envs), `--rerun-triggers mtime`.
- Empty-window guards: bcftools rules `zcat | wc -l`, run_homer placeholder (944df60).
- Attempt-scaled mem/runtime in all rules (`retries: 2`); no OOM/timeout retries seen in run 5.

## Test results (2026-11-15)
| RBP | repro windows | top HOMER de novo | RBPNet profile r (total>50) | sig GO terms |
|---|---|---|---|---|
| RBFOX2 | 27467 | TGCATG (matches) | 0.69 | 646 |
| PTBP1 | 14178 | TYYCTTTCY (CU-rich, matches) | 0.63 | 368 |
| QKI | 31380 | TACTAAC (matches) | 0.68 | 421 |
| TIA1 | 2681 | G-rich; TTTKGGGG #5 (weak U-rich) | 0.70 | 137 |
| U2AF2 | 58251 | TTYYTYCC (polypyrimidine, matches) | 0.76 | 253 |

## Next steps
- Ask the user whether to resume the stopped revision runs (and in what order/concurrency).

## Gotchas
- charlene_paper_20250521 repo not used (user). All fixes go here on `dev/brian/2025revision`.
- `pkill -f run_skipper.sh` also matches the calling shell; kill by PID.
- Snakemake 9 SIGTERM drains running jobs; to stop fast, scancel jobs by run UUID (squeue job name) then SIGINT.
- GPU rules on hotel-v100; jobs pend on account-wide MaxNodePerAccount.
- ps-yeolab3 shared; stop if <5T free.
