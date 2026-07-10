
import pandas as pd
from pathlib import Path
locals().update(config)

rule upload_hf_dataset:
    input:
        parsed_data = "output/ml/rbpnet_data/{experiment_label}/prep_done",
    output:
        "output/ml/rbpnet_data/parquet/{experiment_label}.train.parquet",
        "output/ml/rbpnet_data/parquet/{experiment_label}.test.parquet",
        "output/ml/rbpnet_data/parquet/{experiment_label}.val.parquet",
    # container:
    #     "docker://brianyee/eugene-tools:0.1.2" # THIS DOCKER IS NOT UPDATED WITH PYAROOW YET? NO SPACE LEFT ON DEVICE PROBLEM. PLUS CHARLENE CAN NEVER PULL CORRECTLY
    resources:
        mem_mb=80000,
        runtime=40
    singularity:
        "/tscc/nfs/home/hsher/scratch/singularity/eugene_nt_lora_latest.sif"
    shell:
        """
        export NUMBA_DISABLE_CACHING=1
        export NUMBA_CACHE_DIR=/tscc/lustre/ddn/scratch/${{USER}} # TODO: HARCODED IS BAD
        export MPLCONFIGDIR=/tscc/lustre/ddn/scratch/${{USER}}
        export HF_HOME=/tscc/nfs/home/${{USER}}/scratch/.cache/huggingface
        python {TOOL_DIR}/upload_to_hf_dataset.py {wildcards.experiment_label} output/ml/rbpnet_data/parquet
        """