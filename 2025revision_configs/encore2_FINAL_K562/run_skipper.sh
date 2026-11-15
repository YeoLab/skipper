#!/bin/bash
# Supervised Skipper run for the config.yaml in this folder: retries when Slurm refuses
# submissions ("sbatch failed"), stops on success or any other failure.
# usage: ./run_skipper.sh [extra snakemake args, e.g. -n for a dry run]
wd=$(dirname $(readlink -f $0))
R=/tscc/projects/ps-yeolab3/bay001/codebase/skipper
TARGETS="all"
module load singularitypro
# nodejs/galyleo/python3essential modules put their own python ahead of snakemake's conda envs on PATH; jobs inherit this env
module unload nodejs galyleo python3essential 2>/dev/null
source ~/.bashrc >/dev/null 2>&1
conda activate snakemake9
unset SLURM_JOB_ID
cd $wd
for attempt in $(seq 1 50); do
    [ -f snakemake.log ] && mv snakemake.log snakemake.$(date +%Y%m%d-%H%M%S).log
    echo "$(date) attempt $attempt"
    snakemake $TARGETS -s $R/Skipper.py \
        --configfile $wd/config.yaml \
        --profile $R/profiles/tscc2_snakemake9 \
        --rerun-triggers mtime \
        --conda-prefix /tscc/lustre/ddn/scratch/bay001/snakeconda \
        --singularity-prefix /tscc/lustre/ddn/scratch/bay001/singularity \
        "$@" > snakemake.log 2>&1 < /dev/null
    rc=$?
    if [ $rc -eq 0 ]; then echo "$(date) done"; exit 0; fi
    if grep -q "sbatch failed" snakemake.log; then
        echo "$(date) sbatch refused (rc=$rc); retrying in 5 min"; sleep 300
    else
        echo "$(date) snakemake failed (rc=$rc); see snakemake.log"; exit $rc
    fi
done
echo "$(date) giving up after 50 attempts"; exit 1
