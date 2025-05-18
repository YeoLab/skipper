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