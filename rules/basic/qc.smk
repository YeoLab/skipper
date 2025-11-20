locals().update(config)
import pandas as pd

manifest = pd.read_csv(MANIFEST, comment = "#", index_col = False).dropna(subset=['Experiment','Sample'])
experiment_to_sample = dict(zip(manifest["Experiment"], manifest["Sample"]))

rule multiqc:
    input: 
        trimmed_fastqc = lambda wildcards: (
            [f"output/QC/fastqc/processed/{replicate_label}.trimmed.umi_fastqc.zip" 
             for replicate_label in experiment_to_replicate_labels[wildcards.experiment_label]]
            if config['protocol'] == 'ENCODE4' else
            [f"output/QC/fastqc/processed/{replicate_label}-trimmed-pair1_fastqc.zip" 
             for replicate_label in experiment_to_replicate_labels[wildcards.experiment_label]]
            + [f"output/QC/fastqc/processed/{replicate_label}-trimmed-pair2_fastqc.zip" 
               for replicate_label in experiment_to_replicate_labels[wildcards.experiment_label]]
        ),
        initial_fastqc = lambda wildcards: (
            [f"output/QC/fastqc/initial/{replicate_label}_fastqc.zip"
             for replicate_label in experiment_to_replicate_labels[wildcards.experiment_label]]
            if config['protocol'] == 'ENCODE4' else
            [f"output/QC/fastqc/initial/{replicate_label}-1_fastqc.zip"
             for replicate_label in experiment_to_replicate_labels[wildcards.experiment_label]]
            + [f"output/QC/fastqc/initial/{replicate_label}-2_fastqc.zip"
               for replicate_label in experiment_to_replicate_labels[wildcards.experiment_label]]
        ),
        star_log = lambda wildcards: [
            f"output/secondary_results/bams/raw/genome/{replicate_label}.genome.Log.final.out"
            for replicate_label in experiment_to_replicate_labels[wildcards.experiment_label]
        ],
        fastp = lambda wildcards: (
            [f"output/QC/fastp/{replicate_label}.fastp.json"
             for replicate_label in experiment_to_replicate_labels[wildcards.experiment_label]]
            if config['protocol'] == 'ENCODE4' else []
        ),
        trimmed = lambda wildcards: (
            [f"output/secondary_results/fastqs/trimmed/{replicate_label}-trimmed.log"
             for replicate_label in experiment_to_replicate_labels[wildcards.experiment_label]]
            if config['protocol'] == 'ENCODE4' else []
        )
    output:
        multiqc_results = directory("output/multiqc/{experiment_label}/multiqc_data/"),
        multiqc_plots = directory("output/multiqc/{experiment_label}/multiqc_plots/"),
        multiqc_report = "output/multiqc/{experiment_label}/multiqc_report.html"
    benchmark: "benchmarks/multiqc/{experiment_label}.multiqc.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.multiqc.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.multiqc.err",
    conda:
        "envs/multiqc2.yaml"
    resources:
        mem_mb=4000,
        runtime="30m"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting multiqc" | tee -a {log.stdout}

        (ls {input.trimmed_fastqc} {input.initial_fastqc} {input.star_log} {input.fastp} {input.trimmed} \
            > output/multiqc/{wildcards.experiment_label}/files.txt) >> {log.stdout} 2> {log.stderr}

        multiqc \
            --outdir output/multiqc/{wildcards.experiment_label} \
            -f \
            --export \
            --data-format json \
            --file-list output/multiqc/{wildcards.experiment_label}/files.txt \
            >> {log.stdout} 2>> {log.stderr}

        echo "[`date`] Finished multiqc" | tee -a {log.stdout}
        """

