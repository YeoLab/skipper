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
        mem_mb=8000,
        runtime=60
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
    conda:
        "envs/metadensity.yaml"
    resources:
        mem_mb=40000,
        runtime=5
    shell:
        """
        module purge;
        python {TOOL_DIR}/quantify_gcbias.py {input} {output}
        """

def get_bams(wildcards):
    ''' return a list of final bam for all IP and Input replicate given experiment label '''
    if config['protocol']=='ENCODE':
        return [f"output/bams/dedup/genome_R{INFORMATIVE_READ}/{replicate_label}.genome.Aligned.sort.dedup.R{INFORMATIVE_READ}.bam"
        for replicate_label in experiment_to_replicate_labels[wildcards.experiment_label]
        ]
    else:
        return [f"output/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam"
        for replicate_label in experiment_to_replicate_labels[wildcards.experiment_label]
        ]

rule nread_in_finemapped_regions:
    input:
        bam=lambda wildcards: get_bams(wildcards),
        bed="output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz"
    output:
        sorted_bed="output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.sorted.bed.gz",
        nread_in_finemapped_regions = "output/qc/{experiment_label}.nread_in_finemapped_regions.txt"
    container:
        "docker://howardxu520/skipper:bigwig_1.0"
    params:
        genome=os.path.join(config['STAR_DIR'], 'chrNameLength.txt')
    resources:
        mem_mb=lambda wildcards, attempt: 32000 * (2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 5 * (2 ** (attempt - 1)),
    shell:
        """
        zcat {input.bed} | sort -k1,1 -k2,2 -V | gzip -c > {output.sorted_bed};
        for bam in {input.bam}
        do
            nread_in_peak=$(bedtools coverage \
            -a {output.sorted_bed} \
            -b $bam \
            -sorted \
            -counts \
            -g {params.genome} \
            -s | cut -f 10 | awk '{{sum += $NF}} END {{print sum}}')
            replicate_label=$(basename $bam | cut -d'.' -f1)
            echo -e "$replicate_label\t$nread_in_peak" >> {output.nread_in_finemapped_regions}
        done
        """