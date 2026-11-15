import pandas as pd
locals().update(config)

rule prepare_data:
    input:
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz",
        ip_bigwigs = lambda wildcards: expand("output/bigwigs/unscaled/plus/{replicate_label}.unscaled.plus.bw",
        replicate_label = experiment_to_replicate_labels[wildcards.experiment_label]),
    output:
        "output/ml/rbpnet_data/{experiment_label}/prep_done"
    # container:
    #     "docker://brianyee/eugene-tools:0.1.2" # THIS DOCKER IS NOT UPDATED WITH PYAROOW YET? NO SPACE LEFT ON DEVICE PROBLEM. PLUS CHARLENE CAN NEVER PULL CORRECTLY
    resources:
        mem_mb=lambda wildcards, attempt: 34000 * attempt,
        runtime=lambda wildcards, attempt: 40 * attempt
    container: 
        "/tscc/nfs/home/bay001/bin/eugene-tools_0.1.2.sif"
    shell:
        """
        export NUMBA_DISABLE_CACHING=1
        export NUMBA_CACHE_DIR=$(mktemp -d) # per-job: a shared numba cache gets corrupted by concurrent jobs
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
        mem_mb=lambda wildcards, attempt: 10000 * attempt,
        runtime=lambda wildcards, attempt: 180 * attempt,
        slurm_partition="hotel-v100",
        slurm_account="csd792",
        slurm_extra="'--qos=condo-gpu' '--gpus=1'",
    container:
        "/tscc/nfs/home/bay001/bin/eugene-tools_0.1.2.sif"
    # container:
    #     "docker://brianyee/eugene-tools:0.1.2" #NO SPACE LEFT ON DEVICE PROBLEM
    shell:
        """
        module load gpu
        export NUMBA_DISABLE_CACHING=1
        export NUMBA_CACHE_DIR=$(mktemp -d) # per-job: a shared numba cache gets corrupted by concurrent jobs
        export MPLCONFIGDIR=/tscc/lustre/ddn/scratch/${{USER}}
        python {RBPNET_PATH}/train.py output/ml/rbpnet_data/{wildcards.experiment_label} \
        output/ml/rbpnet_model/{wildcards.experiment_label}
        """

# reviewer asked us to benchmark model without additional layer
rule train_original_model:
    input:
        "output/ml/rbpnet_data/{experiment_label}/prep_done"
    output:
        model = "output/ml/rbpnet_model_original/{experiment_label}/training_done",
    resources:
        mem_mb=lambda wildcards, attempt: 32000 * attempt,
        runtime=lambda wildcards, attempt: 360 * attempt,
        slurm_partition="hotel-v100",
        slurm_account="csd792",
        slurm_extra="'--qos=condo-gpu' '--gpus=1'",
    container:
        "/tscc/nfs/home/bay001/bin/eugene-tools_0.1.2.sif"
    # container:
    #     "docker://brianyee/eugene-tools:0.1.2" #NO SPACE LEFT ON DEVICE PROBLEM
    shell:
        """
        module load gpu
        export NUMBA_DISABLE_CACHING=1
        export NUMBA_CACHE_DIR=$(mktemp -d) # per-job: a shared numba cache gets corrupted by concurrent jobs
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
        mem_mb=lambda wildcards, attempt: 32000 * attempt,
        runtime=lambda wildcards, attempt: 60 * attempt,
        slurm_partition="hotel-v100",
        slurm_account="csd792",
        slurm_extra="'--qos=condo-gpu' '--gpus=1'",
    container:
        "/tscc/nfs/home/bay001/bin/eugene-tools_0.1.2.sif"
    # container:
    #     "docker://brianyee/eugene-tools:0.1.2"
    shell:
        """
        module load gpu
        export NUMBA_DISABLE_CACHING=1
        export NUMBA_CACHE_DIR=$(mktemp -d) # per-job: a shared numba cache gets corrupted by concurrent jobs
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
        mem_mb=lambda wildcards, attempt: 32000 * attempt,
        runtime=lambda wildcards, attempt: 60 * attempt,
        slurm_partition="hotel-v100",
        slurm_account="csd792",
        slurm_extra="'--qos=condo-gpu' '--gpus=1'",
    container:
        "/tscc/nfs/home/bay001/bin/eugene-tools_0.1.2.sif"
    # container:
    #     "docker://brianyee/eugene-tools:0.1.2"
    shell:
        """
        module load gpu
        export NUMBA_DISABLE_CACHING=1
        export NUMBA_CACHE_DIR=$(mktemp -d) # per-job: a shared numba cache gets corrupted by concurrent jobs
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
        mem_mb=lambda wildcards, attempt: 9000 * attempt,
        runtime=lambda wildcards, attempt: 60 * attempt,
        slurm_partition="hotel-v100",
        slurm_account="csd792",
        slurm_extra="'--qos=condo-gpu' '--gpus=1'",
    container:
        "/tscc/nfs/home/bay001/bin/eugene-tools_0.1.2.sif"
    # container:
    #     "docker://brianyee/eugene-tools:0.1.2"
    shell:
        """
        module load gpu
        export NUMBA_DISABLE_CACHING=1
        export NUMBA_CACHE_DIR=$(mktemp -d) # per-job: a shared numba cache gets corrupted by concurrent jobs
        export MPLCONFIGDIR=/tscc/lustre/ddn/scratch/${{USER}}
        python {RBPNET_PATH}/seqlet.py \
        output/ml/rbpnet_data/{wildcards.experiment_label} \
        output/ml/rbpnet_model/{wildcards.experiment_label}
        """