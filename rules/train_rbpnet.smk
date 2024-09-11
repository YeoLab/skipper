import pandas as pd
locals().update(config)

rule prepare_data:
    input:
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz",
        ip_bigwigs = lambda wildcards: expand("output/bigwigs/unscaled/plus/{replicate_label}.unscaled.plus.bw",
        replicate_label = experiment_to_replicate_labels[wildcards.experiment_label]),
    output:
        "output/ml/rbpnet_data/{experiment_label}/prep_done"
    params:
        error_file = "stderr/prep_rbpnet.{experiment_label}.err",
        out_file = "stdout/prep_rbpnet.{experiment_label}.out",
        run_time = "40:00",
        memory = "80000",
    # container:
    #     "docker://brianyee/eugene-tools:0.1.2" # THIS DOCKER IS NOT UPDATED WITH PYAROOW YET? NO SPACE LEFT ON DEVICE PROBLEM. PLUS CHARLENE CAN NEVER PULL CORRECTLY
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
    params:
        error_file = "stderr/train_rbpnet.{experiment_label}.err",
        out_file = "stdout/train_rbpnet.{experiment_label}.out",
        run_time = "3:40:00",
        memory = "160000",
        gpu = True
    container:
        "/tscc/nfs/home/bay001/eugene-tools_0.1.2.sif"
    # container:
    #     "docker://brianyee/eugene-tools:0.1.2" #NO SPACE LEFT ON DEVICE PROBLEM
    shell:
        """
        module load gpu
        export NUMBA_CACHE_DIR=/tscc/lustre/ddn/scratch/${{USER}} # TODO: HARCODED IS BAD
        export MPLCONFIGDIR=/tscc/lustre/ddn/scratch/${{USER}}
        python {RBPNET_PATH}/train.py output/ml/rbpnet_data/{wildcards.experiment_label} \
        output/ml/rbpnet_model/{wildcards.experiment_label}
        """

rule validation:
    input:
        model = "output/ml/rbpnet_model/{experiment_label}/training_done",
        zarr = "output/ml/rbpnet_data/{experiment_label}/prep_done"
    output:
        validation = "output/ml/rbpnet_model/{experiment_label}/valid/test_data_metric.csv",
    params:
        error_file = "stderr/validate_rbpnet.{experiment_label}.err",
        out_file = "stdout/validate_rbpnet.{experiment_label}.out",
        run_time = "40:00",
        memory = "160000",
        gpu = True
    container:
        "/tscc/nfs/home/bay001/eugene-tools_0.1.2.sif"
    # container:
    #     "docker://brianyee/eugene-tools:0.1.2"
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
        
rule seqlet:
    input:
        validation = "output/ml/rbpnet_model/{experiment_label}/valid/test_data_metric.csv",
        model = "output/ml/rbpnet_model/{experiment_label}/training_done",
        zarr = "output/ml/rbpnet_data/{experiment_label}/prep_done"
    output:
        validation = "output/ml/rbpnet_model/{experiment_label}/motif_done",
    params:
        error_file = "stderr/seqlet_rbpnet.{experiment_label}.err",
        out_file = "stdout/seqlet_rbpnet.{experiment_label}.out",
        run_time = "40:00",
        memory = "160000",
        gpu = True
    container:
        "/tscc/nfs/home/bay001/eugene-tools_0.1.2.sif"
    # container:
    #     "docker://brianyee/eugene-tools:0.1.2"
    shell:
        """
        module load gpu
        export NUMBA_CACHE_DIR=/tscc/lustre/ddn/scratch/${{USER}}
        export MPLCONFIGDIR=/tscc/lustre/ddn/scratch/${{USER}}
        python {RBPNET_PATH}/seqlet.py \
        output/ml/rbpnet_data/{wildcards.experiment_label} \
        output/ml/rbpnet_model/{wildcards.experiment_label}
        """