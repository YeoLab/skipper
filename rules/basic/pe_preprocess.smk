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
        star_dir = directory(STAR_DIR),
        chrom_sizes = CHROM_SIZES,
    threads: 8
    resources:
        mem_mb = 48000,
        runtime = "2h",
    benchmark: "benchmarks/run_star_genome_generate.txt"
    log: "logs/run_star_genome_generate.log"
    conda:
        "envs/star.yaml"
    shell:        
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting star_genome_generate." | tee {log}

        tmp_gff=tmp/tmp.gff
        zcat {input.gff} > $tmp_gff

        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output.star_dir} \
            --genomeFastaFiles {input.fasta_file} \
            --sjdbGTFfile $tmp_gff \
            --sjdbOverhang 99 \
            2>&1 | tee -a {log}

        rm -f $tmp_gff

        echo "[`date`] Finished star_genome_generate." | tee -a {log}
        """

rule copy_with_umi:
    input:
        fq_1 = lambda wildcards: config['replicate_label_to_fastq_1'][wildcards.replicate_label],
        fq_2 = lambda wildcards: config['replicate_label_to_fastq_2'][wildcards.replicate_label],
    output:
        fq_1 = temp("output/fastqs/copy/{replicate_label}-1.fastq.gz"), 
        fq_2 = temp("output/fastqs/copy/{replicate_label}-2.fastq.gz"),      
    threads: 2
    resources:
        runtime="2h",
        mem_mb=2000
    benchmark: "benchmarks/umi/unassigned_experiment.{replicate_label}.copy_with_umi.txt"
    log: "logs/{replicate_label}.copy_with_umi.log"
    shell:
        r"""
        set -euo pipefail
        echo "[`date`] Starting copy_with_umi" | tee {log}

        zcat {input.fq_1} \
            | awk 'NR % 4 != 1 {print} NR % 4 == 1 {split($1,header,":"); print $1 ":" substr(header[1],2,length(header[1]) - 1) }' \
            | gzip \
            > {output.fq_1} \
            2>&1 | tee -a {log}

        zcat {input.fq_2} \
            | awk 'NR % 4 != 1 {print} NR % 4 == 1 {split($1,header,":"); print $1 ":" substr(header[1],2,length(header[1]) - 1) }' \
            | gzip \
            > {output.fq_2} \
            2>&1 | tee -a {log}

        echo "[`date`] Finished copy_with_umi" | tee -a {log}
        """

rule run_initial_fastqc:
    input:
        r1 = rules.copy_with_umi.output.fq_1,
        r2 = rules.copy_with_umi.output.fq_2
    output:
        report_r1 = "output/fastqc/initial/{replicate_label}-1_fastqc.html",
        zip_file_r1 = "output/fastqc/initial/{replicate_label}-1_fastqc.zip",
        report_r2 = "output/fastqc/initial/{replicate_label}-2_fastqc.html",
        zip_file_r2 = "output/fastqc/initial/{replicate_label}-2_fastqc.zip",
    threads: 2
    benchmark: "benchmarks/fastqc/unassigned_experiment.{replicate_label}.initial_fastqc.txt"
    log: "logs/{replicate_label}.run_initial_fastqc.log"
    conda:
        "envs/fastqc.yaml"
    resources:
        mem_mb=16000,
        runtime="3h"
    params:
        outdir="output/fastqc/initial/"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting run_initial_fastqc" | tee {log}

        fastqc {input.r1} \
            --extract \
            --outdir {params.outdir} \
            -t {threads} \
            2>&1 | tee -a {log}

        fastqc {input.r2} \
            --extract \
            --outdir {params.outdir} \
            -t {threads} \
            2>&1 | tee -a {log}

        echo "[`date`] Finished run_initial_fastqc" | tee -a {log}
        """
        
rule trim_fastq_encode:
    input:
        fq_1 = rules.copy_with_umi.output.fq_1,
        fq_2 = rules.copy_with_umi.output.fq_2,
        adapter_1 = lambda wildcards: config['replicate_label_to_adapter_1'][wildcards.replicate_label],
        adapter_2 = lambda wildcards: config['replicate_label_to_adapter_2'][wildcards.replicate_label],
    output:
        fq_1_trimmed = temp("output/fastqs/trimmed/{replicate_label}-trimmed-pair1.fastq.gz"), 
        fq_2_trimmed = temp("output/fastqs/trimmed/{replicate_label}-trimmed-pair2.fastq.gz"), 
        metrics = "output/fastqs/trimmed/{replicate_label}-trimmed.log"
    threads: 8
    resources:
        tmpdir = "/tscc/nfs/home/hsher/scratch/singularity_tmp",
        mem_mb=16000,
        runtime="3h"
    benchmark: "benchmarks/trim/unassigned_experiment.{replicate_label}.trim.txt"
    log: "logs/{replicate_label}.trim_fastq_encode.log"
    conda:
        "envs/skewer.yaml"
    shell:
        r"""
        set -euo pipefail
        echo "[`date`] Starting trim_fastq_encode" | tee {log}

        skewer \
            -t {threads} \
            -x {input.adapter_1} \
            -y {input.adapter_2} \
            -o output/fastqs/trimmed/{wildcards.replicate_label} \
            -z -r 0.2 -d 0.2 -q 13 -l 20 \
            {input.fq_1} {input.fq_2} \
            2>&1 | tee -a {log}

        echo "[`date`] Finished trim_fastq_encode" | tee -a {log}
        """

rule run_trimmed_fastqc:
    input:
        r1 = rules.trim_fastq_encode.output.fq_1_trimmed,
        r2 = rules.trim_fastq_encode.output.fq_2_trimmed,
    output:
        report_r1 = "output/fastqc/processed/{replicate_label}-trimmed-pair1_fastqc.html",
        zip_file_r1 = "output/fastqc/processed/{replicate_label}-trimmed-pair1_fastqc.zip",
        report_r2 = "output/fastqc/processed/{replicate_label}-trimmed-pair2_fastqc.html",
        zip_file_r2 = "output/fastqc/processed/{replicate_label}-trimmed-pair2_fastqc.zip",
    threads: 2
    benchmark: "benchmarks/fastqc/unassigned_experiment.{replicate_label}.trimmed_fastqc.txt"
    log: "logs/{replicate_label}.run_trimmed_fastqc.log"
    conda:
        "envs/fastqc.yaml"
    resources:
        mem_mb=16000,
        runtime="2h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting run_trimmed_fastqc" | tee {log}

        fastqc {input.r1} \
            --extract \
            --outdir output/fastqc/processed \
            -t {threads} \
            2>&1 | tee -a {log}

        fastqc {input.r2} \
            --extract \
            --outdir output/fastqc/processed \
            -t {threads} \
            2>&1 | tee -a {log}

        echo "[`date`] Finished run_trimmed_fastqc" | tee -a {log}
        """
        
rule align_reads_encode:
    input:
        fq_1 = rules.trim_fastq_encode.output.fq_1_trimmed,
        fq_2 = rules.trim_fastq_encode.output.fq_2_trimmed
        star_sjdb = STAR_DIR
    output:
        ubam = temp("output/bams/raw/genome/{replicate_label}.genome.Aligned.out.bam"),
        log= "output/bams/raw/genome/{replicate_label}.genome.Log.final.out",
    threads: 8
    params:
        outprefix = "output/bams/raw/genome/{replicate_label}.genome.",
        rg = "{replicate_label}"
    benchmark: "benchmarks/align/unassigned_experiment.{replicate_label}.align_reads_genome.txt"
    log: "logs/{replicate_label}.align_reads_encode.log"
    conda:
        "envs/star.yaml"
    resources:
        mem_mb=160000,
        runtime="2h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting align_reads_encode" | tee {log}

        STAR \
            --alignEndsType EndToEnd \
            --genomeDir {input.star_sjdb} \
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
            2>&1 | tee -a {log}

        echo "[`date`] Finished align_reads_encode" | tee -a {log}
        """

rule sort_bam:
    input:
        bam = "output/bams/raw/{ref}/{replicate_label}.{ref}.Aligned.out.bam",
    output:
        sort = "output/bams/raw/{ref}/{replicate_label}.{ref}.Aligned.sort.bam",
    threads: 2
    benchmark: "benchmarks/sort/{ref}/unassigned_experiment.{replicate_label}.sort_bam.txt"
    log: "logs/{replicate_label}.{ref}.sort_bam.log"
    conda:
        "envs/bedbam_tools.yaml"
    resources:
        mem_mb = 16000,
        runtime = "1h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting sort_bam" | tee {log}

        samtools sort \
            -T {wildcards.replicate_label} \
            -@ {threads} \
            -o {output.sort} \
            {input.bam} \
            2>&1 | tee -a {log}

        echo "[`date`] Finished sort_bam" | tee -a {log}
        """

rule index_bams:
    input:
        bam = "output/bams/{round}/{ref}/{replicate_label}.Aligned.{mid}.bam"
    output:
        ibam = "output/bams/{round}/{ref}/{replicate_label}.Aligned.{mid}.bam.bai"
    threads: 2
    benchmark: "benchmarks/index_bam/{round}/{ref}/{mid}/unassigned_experiment.{replicate_label}.index_bam.txt"
    log: "logs/{replicate_label}.{ref}.{mid}.{round}.index_bams.log"
    conda:
        "envs/bedbam_tools.yaml"
    resources:
        mem_mb = 1000,
        runtime = "1h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting index_bams" | tee {log}

        samtools index \
            -@ {threads} \
            {input.bam} \
            2>&1 | tee -a {log}

        echo "[`date`] Finished index_bams" | tee -a {log}
        """

rule dedup_umi_encode:
    input:
        bam = "output/bams/raw/genome/{replicate_label}.genome.Aligned.sort.bam",
        ibam = "output/bams/raw/genome/{replicate_label}.genome.Aligned.sort.bam.bai"
    output:
        bam_dedup = "output/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam"
    params:
        prefix = "output/bams/dedup/genome/{replicate_label}.genome.sort"
    resources:
        mem_mb = 34000,
        runtime = "2h",
        tmpdir = TMPDIR
    benchmark: "benchmarks/dedup/genome/unassigned_experiment.{replicate_label}.dedup_umi.txt"
    log: "logs/{replicate_label}.dedup_umi_encode.log"
    conda:
        "envs/umicollapse.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting dedup_umi_encode" | tee {log}

        java -server -Xms32G -Xmx32G -Xss40M \
            -jar /UMICollapse/umicollapse.jar bam \
            -i {input.bam} \
            -o {output.bam_dedup} \
            --umi-sep : \
            --two-pass \
            2>&1 | tee -a {log}

        echo "[`date`] Finished dedup_umi_encode" | tee -a {log}
        """

rule select_informative_read:
    input:
        bam_combined = "output/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam"
    output:
        bam_informative = "output/bams/dedup/genome_R" + str(INFORMATIVE_READ) + "/{replicate_label}.genome.Aligned.sort.dedup.R" + str(INFORMATIVE_READ) + ".bam"
    benchmark: "benchmarks/select/unassigned_experiment.{replicate_label}.select_informative_read.txt"
    log: "logs/{replicate_label}.select_informative_read.log"
    conda:
        "envs/bedbam_tools.yaml"
    resources:
        mem_mb = 10000,
        runtime = "1h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting select_informative_read" | tee {log}

        samtools view \
            -bF {64 if UNINFORMATIVE_READ == 1 else 128} \
            {input.bam_combined} \
            > {output.bam_informative} \
            2>&1 | tee -a {log}

        echo "[`date`] Finished select_informative_read" | tee -a {log}
        """

rule obtain_unique_reads:
    input:
        rules.select_informative_read.output.bam_informative
    output:
        temp("output/QC/{replicate_label}.uniq_fragments")
    benchmark:
        "benchmarks/{replicate_label}.count_uniq_fragments.txt"
    log: 
        "logs/{replicate_label}.obtain_unique_reads.log"
    conda:
        "envs/bedbam_tools.yaml"
    resources:
        mem_mb=10000,
        runtime="1h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting obtain_unique_reads" | tee {log}

        samtools idxstats {input} \
            | awk -F '\t' '{{s+=$3+$4}} END {{print s}}' \
            > {output} \
            2>&1 | tee -a {log}

        echo "[`date`] Finished obtain_unique_reads" | tee -a {log}
        """

rule obtain_aligned_reads:
    input:
        rules.align_reads_encode.output.ubam
    output:
        "output/QC/{replicate_label}.aligned_reads"
    benchmark:
        "benchmarks/{replicate_label}.count_aligned_reads.txt"
    log: 
        "logs/{replicate_label}.obtain_aligned_reads.log"
    conda:
        "envs/bedbam_tools.yaml"
    resources:
        mem_mb=8000,
        runtime="1h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting obtain_aligned_reads" | tee {log}

        samtools idxstats {input} \
            | awk -F '\t' '{{s+=$3+$4}} END {{print s}}' \
            > {output} \
            2>&1 | tee -a {log}

        echo "[`date`] Finished obtain_aligned_reads" | tee -a {log}
        """