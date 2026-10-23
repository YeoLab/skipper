from pathlib import Path
import os
import sys
import pandas as pd
locals().update(config)
UTILS_WORKDIR = config.get('UTILS_WORKDIR', '/tscc/lustre/ddn/scratch/bay001')
WD = Path(UTILS_WORKDIR) / 'gnomAD_per_experiment'
workdir: str(WD)

"""
Regenerate output/variants/gnomAD_roulette/{experiment}.total.csv for the selected
RBPNet models, using the main pipeline's variant rules (rules/variants_rbpnet.smk).
Trained models are linked in from model_dir instead of retrained; finemapped windows
are taken from the first of FINEMAP_DIRS (Skipper WORKDIRs) that has them.

snakemake -s utils/score_gnomAD_per_experiment.smk \
    --configfile encode_configs/encode_pe_rules_config_K562_20230929.yaml \
    --profile profiles/tscc2_snakemake9 -n
"""

tabledir = Path('/tscc/nfs/home/hsher/projects/ENCODE/tables/')
metrics = pd.read_csv(tabledir/'model_performance.csv')
selected_models = metrics.loc[metrics['selected'], 'Experiment'].sort_values().tolist()
model_dir = Path('/tscc/nfs/home/hsher/ps-yeolab5/rbpnet_models')
FINEMAP_DIRS = config.get('FINEMAP_DIRS', [
    f'{UTILS_WORKDIR}/finemap_regen/{cohort}' for cohort in ['ENCODE3_K562', 'ENCODE3_HepG2', 'ENCORE_K562', 'ENCORE_HepG2']
] + [
    f'/tscc/nfs/home/hsher/scratch/{cohort}' for cohort in ['ENCODE3_K562', 'ENCODE3_HepG2', 'ENCORE_K562', 'ENCORE_HepG2']
])
EXPERIMENTS = config.get('GNOMAD_EXPERIMENTS', selected_models)

# the imported rules read models and windows from fixed paths under output/, so link them in
os.makedirs(WD / 'output/ml/rbpnet_model', exist_ok=True)
os.makedirs(WD / 'output/finemapping/mapped_sites', exist_ok=True)
missing_finemap = []
for exp in EXPERIMENTS:
    model_link = WD / 'output/ml/rbpnet_model' / exp
    if not model_link.exists():
        model_link.symlink_to(model_dir / exp)
    finemap_link = WD / f'output/finemapping/mapped_sites/{exp}.finemapped_windows.bed.gz'
    if not finemap_link.exists():
        found = [Path(d) / f'output/finemapping/mapped_sites/{exp}.finemapped_windows.bed.gz' for d in FINEMAP_DIRS]
        found = [f for f in found if f.exists()]
        if found:
            finemap_link.symlink_to(found[0])
        else:
            missing_finemap.append(exp)
if missing_finemap:
    sys.stderr.write(f"WARNING: no finemapped windows yet for {len(missing_finemap)} experiments; regenerate them first (see encode_configs/finemap_regen_bay001)\n")

rule all:
    input:
        expand("output/variants/gnomAD_roulette/{experiment_label}.total.csv",
               experiment_label = [e for e in EXPERIMENTS if e not in missing_finemap])

module variants:
    snakefile:
        "../rules/variants_rbpnet.smk"
    config:
        config

use rule fetch_SNP_from_gnomAD_and_roulette, slop_finemap, fetch_peak_sequence, fetch_variant_sequence, join_gnomAD_info from variants

# score every chromosome's ref/alt fasta of an experiment in one GPU job (model loaded once)
# instead of the module's per-chromosome score_fa; job startup dominates otherwise (~1 min each)
REPO_TOOLS = Path(workflow.basedir).parent / 'tools'
CHROMS = list(range(1, 23))
TYPES = ['ref', 'alt']

rule score_experiment:
    input:
        model = "output/ml/rbpnet_model/{experiment_label}/training_done",
        fa = expand("output/variants/gnomAD_roulette/{{experiment_label}}.chr{chr_number}.{type}.fa",
                    chr_number = CHROMS, type = TYPES),
    output:
        score = temp(expand("output/variants/gnomAD_roulette/{{experiment_label}}.chr{chr_number}.{type}.score.csv",
                            chr_number = CHROMS, type = TYPES)),
    wildcard_constraints:
        experiment_label = "[^./]+"
    params:
        pairs = lambda wildcards, input, output: ' '.join(f'{fa} {out}' for fa, out in zip(input.fa, output.score)),
    threads: 1
    resources:
        mem_mb=80000,
        runtime=30,
        slurm_partition="rtx3090",
        slurm_account="csd792",
        slurm_extra="'--qos=condo-gpu' '--gpus=1'",
    container:
        "/tscc/nfs/home/bay001/bin/eugene-tools_0.1.2.sif"
    shell:
        """
        export NUMBA_CACHE_DIR=/tscc/lustre/ddn/scratch/${{USER}}
        export MPLCONFIGDIR=/tscc/lustre/ddn/scratch/${{USER}}
        export PYTHONPATH={RBPNET_PATH}
        python {REPO_TOOLS}/score_fa_batch.py \
            output/ml/rbpnet_model/{wildcards.experiment_label}/ \
            output/variants/gnomAD_roulette/zarr/{wildcards.experiment_label} \
            {params.pairs}
        """
