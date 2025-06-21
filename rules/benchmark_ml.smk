locals().update(config)

rule test_sdata_to_fasta:
    input:
        zarr = "output/ml/rbpnet_data/{experiment_label}/prep_done"
    output:
        fasta = temp("output/ml/rbpnet_data/{experiment_label}/test_data.fasta"),
    resources:
        mem_mb=80000,
        runtime="10m",
    container:
        "/tscc/nfs/home/bay001/eugene-tools_0.1.2.sif"
    shell:
        """
        export NUMBA_CACHE_DIR=/tscc/lustre/ddn/scratch/${{USER}}
        export MPLCONFIGDIR=/tscc/lustre/ddn/scratch/${{USER}}
        python {TOOL_DIR}/sdata_to_fasta.py \
            output/ml/rbpnet_data/{wildcards.experiment_label}/test.zarr \
            {output}
        """

rule score_with_homer:
    input:
        fasta = "output/ml/rbpnet_data/{experiment_label}/test_data.fasta",
        homer_pwm = lambda wildcards: config[f'{wildcards.motif_source}_mapping_df'].query('Experiment == @wildcards.experiment_label')[wildcards.motif_source].values[0],
    output:
        score = "output/ml/benchmark/homer/{motif_source}/{experiment_label}.csv",
        outdir = temp(directory("output/ml/benchmark/homer/{motif_source}/{experiment_label}.temp_output_dir"))
    wildcard_constraints:
        motif_source="(RBNS|SELEX)",
    resources:
        mem_mb=80000,
        runtime="10m",
    container:
        "docker://howardxu520/skipper:Homer_4.11"
    shell:
        """
        findMotifs.pl {input.fasta} fasta \
            {output.outdir} \
            -find {input.homer_pwm} | grep '+' > {output.score} # only keep the sense strand to fasta
        """
rule score_with_internal_homer:
    input:
        fasta = "output/ml/rbpnet_data/{experiment_label}/test_data.fasta",
        homer_pwm = "output/homer/finemapped_results/{experiment_label}/homerMotifs.all.motifs",
    output:
        score = "output/ml/benchmark/homer/internal/{experiment_label}.csv",
        outdir = temp(directory("output/ml/benchmark/homer/internal/{experiment_label}.temp_output_dir"))
    resources:
        mem_mb=80000,
        runtime="40m",
    container:  
        "docker://howardxu520/skipper:Homer_4.11"
    shell:
        """
        # if either sequence or homer PWM is empty
        if [ ! -s {input.fasta} ] || [ ! -s {input.homer_pwm} ]; then
            touch {output.score};
            mkdir {output.outdir};
        else
            findMotifs.pl {input.fasta} fasta \
                {output.outdir} \
                -find {input.homer_pwm} | grep '+' > {output.score};
        fi
        """

use rule score_with_internal_homer as score_with_mcross_homer with:
    input:
        fasta = "output/ml/rbpnet_data/{experiment_label}/test_data.fasta",
        homer_pwm = "output/ctk/ctk_mcross/mcross/{data_types}.{experiment_label}/{data_types}.{experiment_label}.homer"
    output:
        score = "output/ml/benchmark/homer/{data_types}_mcross/{experiment_label}.csv",
        outdir = temp(directory("output/ml/benchmark/homer/{data_types}_mcross/{experiment_label}.temp_output_dir"))

rule calculate_pearson_auprc_for_homer:
    input:
        score = "output/ml/benchmark/homer/{motif_source}/{experiment_label}.csv",
        zarr = "output/ml/rbpnet_data/{experiment_label}/prep_done",
        test_data = "output/ml/rbpnet_data/{experiment_label}/test_data.fasta",
        reproducible_windows = "output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz"
    output:
        pearson_auprc = "output/ml/benchmark/homer/{motif_source}/{experiment_label}.pearson_auprc.csv",
    resources:
        mem_mb=80000,
        runtime="10m",
    wildcard_constraints:
        motif_source="(RBNS|SELEX)",
    container:
        "/tscc/nfs/home/bay001/eugene-tools_0.1.2.sif"
    shell:
        """
        export NUMBA_CACHE_DIR=/tscc/lustre/ddn/scratch/${{USER}}
        export MPLCONFIGDIR=/tscc/lustre/ddn/scratch/${{USER}}
        python {TOOL_DIR}/benchmark_homer.py \
            . \
            {wildcards.experiment_label} \
            {wildcards.motif_source} \
            {RBPNET_PATH} \
            {output.pearson_auprc} 
        """
rule calculate_pearson_auprc_for_internal_homer:
    input:
        score = "output/ml/benchmark/homer/internal/{experiment_label}.csv",
        zarr = "output/ml/rbpnet_data/{experiment_label}/prep_done",
        test_data = "output/ml/rbpnet_data/{experiment_label}/test_data.fasta",
        reproducible_windows = "output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz"
    output:
        pearson_auprc = "output/ml/benchmark/homer/internal/{experiment_label}.pearson_auprc.csv",
    resources:
        mem_mb=80000,
        runtime="10m",
    container:
        "/tscc/nfs/home/bay001/eugene-tools_0.1.2.sif"
    shell:
        """
        export NUMBA_CACHE_DIR=/tscc/lustre/ddn/scratch/${{USER}}
        export MPLCONFIGDIR=/tscc/lustre/ddn/scratch/${{USER}}
        python {TOOL_DIR}/benchmark_homer.py \
            . \
            {wildcards.experiment_label} \
            internal \
            {RBPNET_PATH} \
            {output.pearson_auprc} 
        """

rule calculate_pearson_auprc_for_mcross_homer:
    input:
        score = "output/ml/benchmark/homer/{data_types}_mcross/{experiment_label}.csv",
        zarr = "output/ml/rbpnet_data/{experiment_label}/prep_done",
        test_data = "output/ml/rbpnet_data/{experiment_label}/test_data.fasta",
        reproducible_windows = "output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz"
    output:
        pearson_auprc = "output/ml/benchmark/homer/{data_types}_mcross/{experiment_label}.pearson_auprc.csv",
    resources:
        mem_mb=80000,
        runtime="10m",
    wildcard_constraints:
        data_types="(CIMS|CITS)",
    container:
        "/tscc/nfs/home/bay001/eugene-tools_0.1.2.sif"
    shell:
        """
        export NUMBA_CACHE_DIR=/tscc/lustre/ddn/scratch/${{USER}}
        export MPLCONFIGDIR=/tscc/lustre/ddn/scratch/${{USER}}
        python {TOOL_DIR}/benchmark_homer.py \
            . \
            {wildcards.experiment_label} \
            {wildcards.data_types}_mcross \
            {RBPNET_PATH} \
            {output.pearson_auprc} 
        """