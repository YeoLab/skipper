# 2025revision_configs

Frozen run parameters (`config.yaml`, `manifest.csv`, `run_skipper.sh` per run) for the full Skipper re-run
(Cell Systems 2026 revision). Live run folders (WORKDIR, outputs, logs):
`/tscc/nfs/home/bay001/projects/charlene_paper_20250521/permanent_data/skipper2025revision/<run>/`.
`manifest.csv` here is a copy of the manifest each `config.yaml` points to. Source manifests:
`/tscc/projects/ps-yeolab5/rbp-portal/enco*/eclip/Skipper/2.0.0/v2/output/*.csv`.

- Pipeline: this repo, branch `dev/brian/2025revision`, target `all` (basic output + RBPNet/variant scoring).
- Environment: `module load singularitypro; module unload nodejs galyleo; conda activate snakemake9`, profile
  `profiles/tscc2_snakemake9`, conda/singularity prefixes `/tscc/lustre/ddn/scratch/bay001/{snakeconda,singularity}`.

## Runs (one folder per manifest x cell line)
Mixed-cell manifests are split by `Cells` (`<run>/manifest.csv`) so each run uses the cell-matched
`gencode.v40.annotation.<cell>_totalrna.gt1` partition (Her et al. 2026: >1 TPM expressed transcriptome).
Overdispersion: `clip` for ENCODE3 (PE), `input` for ENCORE (SE).

| Run | Source manifest | Protocol |
|---|---|---|
| encore2_FINAL_HepG2 / _K562 | encore2/.../Skipper_manifest_FINAL.csv | ENCODE4 (SE) |
| encore_edgecases_HepG2 / _K562 | encore/.../Skipper_manifest_all_20230727_edgecases.csv.wexpt.csv | SE |
| encore_HepG2 | encore/.../Skipper_manifest_HepG2_20230727.csv.wexpt.csv | SE |
| encore_K562 | encore/.../Skipper_manifest_K562_20230727.csv.wexpt.csv | SE |
| encode3_H_HepG2 | encode3/.../encode3_manifest_H_fastq.csv | ENCODE3 (PE) |
| encode3_K_K562 | encode3/.../encode3_manifest_K_fastq.csv | PE |

Each live `<run>/` has `config.yaml`, `run_skipper.sh` (supervised snakemake; `-n` for dry run), outputs in `output/`,
logs `supervise.log`, `snakemake.log`, `.snakemake/slurm_logs/`. Keep the copies here in sync if a run's params change.

## How it was run
Sequentially (one run at a time), stopping if `/tscc/projects/ps-yeolab3` has <5T free (log `sequential.log`):
```
ok() { avail=$(df -B1 --output=avail /tscc/projects/ps-yeolab3 | tail -1); [ $avail -ge 5000000000000 ] || exit 1; }
for r in encore2_FINAL_HepG2 encore_edgecases_HepG2 encore_HepG2 encore_edgecases_K562 encore_K562 \
         encode3_H_HepG2 encore2_FINAL_K562 encode3_K_K562; do
    $r/run_skipper.sh > $r/supervise.log 2>&1; ok
done
# then one retry pass over runs whose supervise.log doesn't end in "done"
```
