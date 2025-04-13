locals().update(config)
rule multiqc:
    input: 
        trimmed_fastqc = lambda wildcards: 
            [f"output/fastqc/processed/{replicate_label}.trimmed.umi_fastqc.zip" 
            for replicate_label in experiment_to_replicate_labels[wildcards.experiment_label]
            ] if config['protocol']=='ENCODE4' 
            else [f"output/fastqc/processed/{replicate_label}-trimmed-pair1_fastqc.zip" 
            for replicate_label in experiment_to_replicate_labels[wildcards.experiment_label]
            ]+[f"output/fastqc/processed/{replicate_label}-trimmed-pair2_fastqc.zip" 
            for replicate_label in experiment_to_replicate_labels[wildcards.experiment_label]
            ],
        initial_fastqc = lambda wildcards: [f"output/fastqc/initial/{replicate_label}_fastqc.zip"
         for replicate_label in experiment_to_replicate_labels[wildcards.experiment_label]
         ] if config['protocol'] == 'ENCODE4'
         else [f"output/fastqc/initial/{replicate_label}-1_fastqc.zip"
         for replicate_label in experiment_to_replicate_labels[wildcards.experiment_label]
         ]+[f"output/fastqc/initial/{replicate_label}-2_fastqc.zip"
         for replicate_label in experiment_to_replicate_labels[wildcards.experiment_label]
         ],
        star_log = lambda wildcards: [f"output/bams/raw/genome/{replicate_label}.genome.Log.final.out" for replicate_label in experiment_to_replicate_labels[wildcards.experiment_label]],
        fastp = lambda wildcards: [f"output/fastp/{replicate_label}.fastp.json" for replicate_label in experiment_to_replicate_labels[wildcards.experiment_label]] if config['protocol']=='ENCODE4' else [],
        trimmed = lambda wildcards: [f"output/fastqs/trimmed/{replicate_label}-trimmed.log" for replicate_label in experiment_to_replicate_labels[wildcards.experiment_label]] if config['protocol']=='ENCODE4' else []
    output:
        multiqc_results = directory("output/multiqc/{experiment_label}/multiqc_data/"),
        multiqc_plots = directory("output/multiqc/{experiment_label}/multiqc_plots/"),
        multiqc_report = "output/multiqc/{experiment_label}/multiqc_report.html"
    params:
        error_file = "stderr/{experiment_label}.multiqc.err",
        out_file = "stdout/{experiment_label}.multiqc.out",
        run_time = "15:00",
        memory = "4000",
        job_name = "multiqc"
    benchmark: "benchmarks/multiqc/{experiment_label}.multiqc.txt"
    container:
        "docker://jeltje/multiqc:1.6"
    resources:
        mem_mb=4000
    shell:
        """
        ls {input.trimmed_fastqc} {input.initial_fastqc} {input.star_log} {input.fastp} {input.trimmed} > output/multiqc/{wildcards.experiment_label}/files.txt
        multiqc --outdir output/multiqc/{wildcards.experiment_label} -f --export --data-format json --file-list output/multiqc/{wildcards.experiment_label}/files.txt
        """

rule quantify_gc_bias:
    input:
        "output/counts/genome/tables/{experiment_label}.tsv.gz"
    output:
        gc_bias = "output/qc/{experiment_label}.gc_bias.txt"
    params:
        error_file = "stderr/{experiment_label}.gc_bias.err",
        out_file = "stdout/{experiment_label}.gc_bias.out",
        run_time = "5:00",
        memory = "40000",
        job_name = "gc_bias"
    conda:
        "envs/metadensity.yaml"
    resources:
        mem_mb=40000
    shell:
        """
        python {TOOL_DIR}/quantify_gcbias.py {input} {output}
        """
