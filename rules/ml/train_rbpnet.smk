import pandas as pd
locals().update(config)

rule prepare_data:
    input:
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz",
        ip_bigwigs = lambda wildcards: expand("output/bigwigs/unscaled/plus/{replicate_label}.unscaled.plus.bw",
        replicate_label = experiment_to_replicate_labels[wildcards.experiment_label]),
    output:
        "output/ml/rbpnet_data/{experiment_label}/prep_done"
    resources:
        mem_mb=80000,
        runtime=40
    container: 
        "/tscc/nfs/home/bay001/eugene-tools_0.1.2.sif"
    shell:
        """
        export NUMBA_CACHE_DIR=/tscc/lustre/ddn/scratch/${{USER}} # TODO: HARCODED IS BAD
        export MPLCONFIGDIR=/tscc/lustre/ddn/scratch/${{USER}}
        python {RBPNET_PATH}/prep_data.py {CONFIG_PATH} {wildcards.experiment_label} \
            output/ml/rbpnet_data/{wildcards.experiment_label}
        """

rule train_model:
    input:
        "output/ml/rbpnet_data/{experiment_label}/prep_done"
    output:
        model = "output/ml/rbpnet_model/{experiment_label}/training_done",
    resources:
        mem_mb=160000,
        runtime="3h",
        slurm_partition="rtx3090",
        slurm_account="csd792",
        slurm_extra="'--qos=condo-gpu' '--gpus=1'",
    container:
        "/tscc/nfs/home/bay001/eugene-tools_0.1.2.sif"
    shell:
        """
        module load gpu
        export NUMBA_CACHE_DIR=/tscc/lustre/ddn/scratch/${{USER}} # TODO: HARCODED IS BAD
        export MPLCONFIGDIR=/tscc/lustre/ddn/scratch/${{USER}}
        python {RBPNET_PATH}/train.py output/ml/rbpnet_data/{wildcards.experiment_label} \
        output/ml/rbpnet_model/{wildcards.experiment_label}
        """

rule train_original_model:
    input:
        "output/ml/rbpnet_data/{experiment_label}/prep_done"
    output:
        model = "output/ml/rbpnet_model_original/{experiment_label}/training_done",
    resources:
        mem_mb=160000,
        runtime="3h",
        slurm_partition="rtx3090",
        slurm_account="csd792",
        slurm_extra="'--qos=condo-gpu' '--gpus=1'",
    container:
        "/tscc/nfs/home/bay001/eugene-tools_0.1.2.sif"
    shell:
        """
        module load gpu
        export NUMBA_CACHE_DIR=/tscc/lustre/ddn/scratch/${{USER}} # TODO: HARCODED IS BAD
        export MPLCONFIGDIR=/tscc/lustre/ddn/scratch/${{USER}}
        python {RBPNET_PATH}/train_original_rbpnet.py output/ml/rbpnet_data/{wildcards.experiment_label} \
        output/ml/rbpnet_model_original/{wildcards.experiment_label}
        """

rule validation:
    input:
        model = "output/ml/rbpnet_model/{experiment_label}/training_done",
        zarr = "output/ml/rbpnet_data/{experiment_label}/prep_done"
    output:
        validation = "output/ml/rbpnet_model/{experiment_label}/valid/test_data_metric.csv",
    resources:
        mem_mb=160000,
        runtime="1h",
        slurm_partition="rtx3090",
        slurm_account="csd792",
        slurm_extra="'--qos=condo-gpu' '--gpus=1'",
    container:
        "/tscc/nfs/home/bay001/eugene-tools_0.1.2.sif"
    shell:
        """
        module load gpu
        export NUMBA_CACHE_DIR=/tscc/lustre/ddn/scratch/${{USER}}
        export MPLCONFIGDIR=/tscc/lustre/ddn/scratch/${{USER}}
        python {RBPNET_PATH}/validate.py {CONFIG_PATH} \
        {wildcards.experiment_label} \
        output/ml/rbpnet_data/{wildcards.experiment_label} \
        output/ml/rbpnet_model/{wildcards.experiment_label}
        """

rule validate_original_model:
    input:
        model = "output/ml/rbpnet_model_original/{experiment_label}/training_done",
        zarr = "output/ml/rbpnet_data/{experiment_label}/prep_done"
    output:
        validation = "output/ml/rbpnet_model_original/{experiment_label}/valid/test_data_metric.csv",
    resources:
        mem_mb=160000,
        runtime="1h",
        slurm_partition="rtx3090",
        slurm_account="csd792",
        slurm_extra="'--qos=condo-gpu' '--gpus=1'",
    container:
        "/tscc/nfs/home/bay001/eugene-tools_0.1.2.sif"
    shell:
        """
        module load gpu
        export NUMBA_CACHE_DIR=/tscc/lustre/ddn/scratch/${{USER}}
        export MPLCONFIGDIR=/tscc/lustre/ddn/scratch/${{USER}}
        python {RBPNET_PATH}/validate.py {CONFIG_PATH} \
        {wildcards.experiment_label} \
        output/ml/rbpnet_data/{wildcards.experiment_label} \
        output/ml/rbpnet_model_original/{wildcards.experiment_label}
        """

rule seqlet:
    input:
        validation = "output/ml/rbpnet_model/{experiment_label}/valid/test_data_metric.csv",
        model = "output/ml/rbpnet_model/{experiment_label}/training_done",
        zarr = "output/ml/rbpnet_data/{experiment_label}/prep_done"
    output:
        validation = "output/ml/rbpnet_model/{experiment_label}/motif_done",
    resources:
        mem_mb=160000,
        runtime="1h",
        slurm_partition="rtx3090",
        slurm_account="csd792",
        slurm_extra="'--qos=condo-gpu' '--gpus=1'",
    container:
        "/tscc/nfs/home/bay001/eugene-tools_0.1.2.sif"
    shell:
        """
        module load gpu
        export NUMBA_CACHE_DIR=/tscc/lustre/ddn/scratch/${{USER}}
        export MPLCONFIGDIR=/tscc/lustre/ddn/scratch/${{USER}}
        python {RBPNET_PATH}/seqlet.py \
        output/ml/rbpnet_data/{wildcards.experiment_label} \
        output/ml/rbpnet_model/{wildcards.experiment_label}
        """