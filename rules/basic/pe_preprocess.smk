import pandas as pd
from functools import reduce
import re
import os
import sys
import glob
from time import sleep
locals().update(config)

rule run_star_genome_generate:
    input:
        gff = "output/gff/filtered.gff3",
        fasta_file = ancient(GENOME),
    output:
        chrom_sizes = CHROM_SIZES,
    params:
        star_dir = STAR_DIR
    threads: 8
    resources:
        mem_mb = 48000,
        runtime = "2h",
    benchmark: "benchmarks/run_star_genome_generate.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/run_star_genome_generate.out",
        stderr = config["WORKDIR"] + "/stderr/run_star_genome_generate.err",
    conda:
        "envs/star.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting star_genome_generate." | tee -a {log.stdout}

        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {params.star_dir} \
            --genomeFastaFiles {input.fasta_file} \
            --sjdbGTFfile {input.gff} \
            --sjdbOverhang 99 \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished star_genome_generate." | tee -a {log.stdout}
        """

rule copy_with_umi:
    input:
        fq_1 = lambda wildcards: config['replicate_label_to_fastq_1'][wildcards.replicate_label],
        fq_2 = lambda wildcards: config['replicate_label_to_fastq_2'][wildcards.replicate_label],
    output:
        fq_1 = temp("output/secondary_results/fastqs/copy/{replicate_label}-1.fastq.gz"), 
        fq_2 = temp("output/secondary_results/fastqs/copy/{replicate_label}-2.fastq.gz"),      
    threads: 2
    resources:
        runtime="2h",
        mem_mb=2000
    benchmark: "benchmarks/umi/unassigned_experiment.{replicate_label}.copy_with_umi.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{replicate_label}.copy_with_umi.out",
        stderr = config["WORKDIR"] + "/stderr/{replicate_label}.copy_with_umi.err",
    shell:
        r"""
        set -euo pipefail
        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting copy_with_umi" | tee -a {log.stdout}

        zcat {input.fq_1} \
            | awk 'NR % 4 != 1 {{print}} NR % 4 == 1 {{split($1,header,":"); print $1 ":" substr(header[1],2,length(header[1]) - 1) }}' \
            | gzip \
            > {output.fq_1} \
            >> {log.stdout} 2> {log.stderr}

        zcat {input.fq_2} \
            | awk 'NR % 4 != 1 {{print}} NR % 4 == 1 {{split($1,header,":"); print $1 ":" substr(header[1],2,length(header[1]) - 1) }}' \
            | gzip \
            > {output.fq_2} \
            >> {log.stdout} 2>> {log.stderr}

        echo "[`date`] Finished copy_with_umi" | tee -a {log.stdout}
        """

rule run_initial_fastqc:
    input:
        r1 = rules.copy_with_umi.output.fq_1,
        r2 = rules.copy_with_umi.output.fq_2
    output:
        report_r1 = "output/QC/fastqc/initial/{replicate_label}-1_fastqc.html",
        zip_file_r1 = "output/QC/fastqc/initial/{replicate_label}-1_fastqc.zip",
        report_r2 = "output/QC/fastqc/initial/{replicate_label}-2_fastqc.html",
        zip_file_r2 = "output/QC/fastqc/initial/{replicate_label}-2_fastqc.zip",
    threads: 2
    benchmark: "benchmarks/fastqc/unassigned_experiment.{replicate_label}.initial_fastqc.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{replicate_label}.run_initial_fastqc.out",
        stderr = config["WORKDIR"] + "/stderr/{replicate_label}.run_initial_fastqc.err",
    conda:
        "envs/fastqc.yaml"
    resources:
        mem_mb=16000,
        runtime="3h"
    params:
        outdir="output/QC/fastqc/initial/"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting run_initial_fastqc" | tee -a {log.stdout}

        fastqc {input.r1} \
            --extract \
            --outdir {params.outdir} \
            -t {threads} \
            >> {log.stdout} 2> {log.stderr}

        fastqc {input.r2} \
            --extract \
            --outdir {params.outdir} \
            -t {threads} \
            >> {log.stdout} 2>> {log.stderr}

        echo "[`date`] Finished run_initial_fastqc" | tee -a {log.stdout}
        """
        
rule trim_fastq_encode:
    input:
        fq_1 = rules.copy_with_umi.output.fq_1,
        fq_2 = rules.copy_with_umi.output.fq_2,
        adapter_1 = lambda wildcards: config['replicate_label_to_adapter_1'][wildcards.replicate_label],
        adapter_2 = lambda wildcards: config['replicate_label_to_adapter_2'][wildcards.replicate_label],
    output:
        fq_1_trimmed = temp("output/secondary_results/fastqs/trimmed/{replicate_label}-trimmed-pair1.fastq.gz"), 
        fq_2_trimmed = temp("output/secondary_results/fastqs/trimmed/{replicate_label}-trimmed-pair2.fastq.gz"), 
        metrics = "output/secondary_results/fastqs/trimmed/{replicate_label}-trimmed.log"
    threads: 8
    resources:
        tmpdir = "/tscc/nfs/home/hsher/scratch/singularity_tmp",
        mem_mb=16000,
        runtime="3h"
    benchmark: "benchmarks/trim/unassigned_experiment.{replicate_label}.trim.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{replicate_label}.trim_fastq_encode.out",
        stderr = config["WORKDIR"] + "/stderr/{replicate_label}.trim_fastq_encode.err",
    conda:
        "envs/skewer.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting trim_fastq_encode" | tee -a {log.stdout}

        skewer \
            -t {threads} \
            -x {input.adapter_1} \
            -y {input.adapter_2} \
            -o output/secondary_results/fastqs/trimmed/{wildcards.replicate_label} \
            -z -r 0.2 -d 0.2 -q 13 -l 20 \
            {input.fq_1} {input.fq_2} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished trim_fastq_encode" | tee -a {log.stdout}
        """

rule run_trimmed_fastqc:
    input:
        r1 = rules.trim_fastq_encode.output.fq_1_trimmed,
        r2 = rules.trim_fastq_encode.output.fq_2_trimmed,
    output:
        report_r1 = "output/QC/fastqc/processed/{replicate_label}-trimmed-pair1_fastqc.html",
        zip_file_r1 = "output/QC/fastqc/processed/{replicate_label}-trimmed-pair1_fastqc.zip",
        report_r2 = "output/QC/fastqc/processed/{replicate_label}-trimmed-pair2_fastqc.html",
        zip_file_r2 = "output/QC/fastqc/processed/{replicate_label}-trimmed-pair2_fastqc.zip",
    threads: 2
    benchmark: "benchmarks/fastqc/unassigned_experiment.{replicate_label}.trimmed_fastqc.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{replicate_label}.run_trimmed_fastqc.out",
        stderr = config["WORKDIR"] + "/stderr/{replicate_label}.run_trimmed_fastqc.err",
    conda:
        "envs/fastqc.yaml"
    resources:
        mem_mb=16000,
        runtime="2h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting run_trimmed_fastqc" | tee -a {log.stdout}

        fastqc {input.r1} \
            --extract \
            --outdir output/QC/fastqc/processed \
            -t {threads} \
        >> {log.stdout} 2> {log.stderr}

        fastqc {input.r2} \
            --extract \
            --outdir output/QC/fastqc/processed \
            -t {threads} \
        >> {log.stdout} 2>> {log.stderr}

        echo "[`date`] Finished run_trimmed_fastqc" | tee -a {log.stdout}
        """
        
rule align_reads_encode:
    input:
        fq_1 = rules.trim_fastq_encode.output.fq_1_trimmed,
        fq_2 = rules.trim_fastq_encode.output.fq_2_trimmed,
        chrom_sizes = CHROM_SIZES,
    output:
        ubam = temp("output/secondary_results/bams/raw/genome/{replicate_label}.genome.Aligned.out.bam"),
        log= "output/secondary_results/bams/raw/genome/{replicate_label}.genome.Log.final.out",
    threads: 8
    params:
        star_sjdb = STAR_DIR,
        outprefix = "output/secondary_results/bams/raw/genome/{replicate_label}.genome.",
        rg = "{replicate_label}"
    benchmark: "benchmarks/align/unassigned_experiment.{replicate_label}.align_reads_genome.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{replicate_label}.align_reads_encode.out",
        stderr = config["WORKDIR"] + "/stderr/{replicate_label}.align_reads_encode.err",
    conda:
        "envs/star.yaml"
    resources:
        mem_mb=160000,
        runtime="2h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting align_reads_encode" | tee -a {log.stdout}

        STAR \
            --alignEndsType EndToEnd \
            --genomeDir {params.star_sjdb} \
            --genomeLoad NoSharedMemory \
            --outBAMcompression 10 \
            --outFileNamePrefix {params.outprefix} \
            --winAnchorMultimapNmax 100 \
            --outFilterMultimapNmax 100 \
            --outFilterMultimapScoreRange 1 \
            --outSAMmultNmax 1 \
            --outMultimapperOrder Random \
            --outFilterScoreMin 10 \
            --outFilterType BySJout \
            --limitOutSJcollapsed 5000000 \
            --outReadsUnmapped None \
            --outSAMattrRGline ID:{wildcards.replicate_label} \
            --outSAMattributes All \
            --outSAMmode Full \
            --outSAMtype BAM Unsorted \
            --outSAMunmapped Within \
            --readFilesCommand zcat \
            --outStd Log \
            --readFilesIn {input.fq_1} {input.fq_2} \
            --runMode alignReads \
            --runThreadN {threads} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished align_reads_encode" | tee -a {log.stdout}
        """

rule sort_bam:
    input:
        bam = "output/secondary_results/bams/raw/{ref}/{replicate_label}.{ref}.Aligned.out.bam",
    output:
        sort = "output/secondary_results/bams/raw/{ref}/{replicate_label}.{ref}.Aligned.sort.bam",
    threads: 2
    benchmark: "benchmarks/sort/{ref}/unassigned_experiment.{replicate_label}.sort_bam.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{replicate_label}.{ref}.sort_bam.out",
        stderr = config["WORKDIR"] + "/stderr/{replicate_label}.{ref}.sort_bam.err",
    conda:
        "envs/bedbam_tools.yaml"
    resources:
        mem_mb = 16000,
        runtime = "1h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting sort_bam" | tee -a {log.stdout}

        samtools sort \
            -T {wildcards.replicate_label} \
            -@ {threads} \
            -o {output.sort} \
            {input.bam} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished sort_bam" | tee -a {log.stdout}
        """

rule index_bams:
    input:
        bam = "output/secondary_results/bams/{round}/{ref}/{replicate_label}.Aligned.{mid}.bam"
    output:
        ibam = "output/secondary_results/bams/{round}/{ref}/{replicate_label}.Aligned.{mid}.bam.bai"
    threads: 2
    benchmark: "benchmarks/index_bam/{round}/{ref}/{mid}/unassigned_experiment.{replicate_label}.index_bam.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{round}.{ref}.{mid}.{replicate_label}.index_bams.out",
        stderr = config["WORKDIR"] + "/stderr/{round}.{ref}.{mid}.{replicate_label}.index_bams.err",
    conda:
        "envs/bedbam_tools.yaml"
    resources:
        mem_mb = 1000,
        runtime = "1h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting index_bams" | tee -a {log.stdout}

        samtools index \
            -@ {threads} \
            {input.bam} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished index_bams" | tee -a {log.stdout}
        """

rule dedup_umi:
    input:
        bam = "output/secondary_results/bams/raw/genome/{replicate_label}.genome.Aligned.sort.bam",
        ibam = "output/secondary_results/bams/raw/genome/{replicate_label}.genome.Aligned.sort.bam.bai"
    output:
        bam_dedup = "output/secondary_results/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam"
    benchmark: "benchmarks/dedup/genome/unassigned_experiment.{replicate_label}.dedup_umi.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{replicate_label}.dedup_umi.out",
        stderr = config["WORKDIR"] + "/stderr/{replicate_label}.dedup_umi.err",
    conda:
        "envs/umicollapse.yaml"
    resources:
        mem_mb = 48000,
        runtime = "2h",
        tmpdir = TMPDIR
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting dedup_umi for {wildcards.replicate_label}" | tee -a{log.stdout}

        JAR=$(dirname $(which umicollapse))/../share/umicollapse*/umicollapse.jar

        java -server -Xms32G -Xmx32G -Xss40M \
            -jar $JAR bam \
            -i {input.bam} \
            -o {output.bam_dedup} \
            --umi-sep : \
            --two-pass \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished dedup_umi for {wildcards.replicate_label}" | tee -a {log.stdout}
        """

rule select_informative_read:
    input:
        bam_combined = "output/secondary_results/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam"
    output:
        bam_informative = "output/secondary_results/bams/dedup/genome_R" + str(INFORMATIVE_READ) + "/{replicate_label}.genome.Aligned.sort.dedup.R" + str(INFORMATIVE_READ) + ".bam"
    benchmark: "benchmarks/select/unassigned_experiment.{replicate_label}.select_informative_read.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{replicate_label}.select_informative_read.out",
        stderr = config["WORKDIR"] + "/stderr/{replicate_label}.select_informative_read.err",
    conda:
        "envs/bedbam_tools.yaml"
    resources:
        mem_mb = 10000,
        runtime = "1h"
    params:
        flag = 64 if UNINFORMATIVE_READ == "1" else 128
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting select_informative_read" | tee -a {log.stdout}

        samtools view \
            -bF {params.flag} \
            {input.bam_combined} \
            > {output.bam_informative} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished select_informative_read" | tee -a {log.stdout}
        """

