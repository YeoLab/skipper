import pandas as pd
from functools import reduce
import re
import os
import sys
import glob
from time import sleep
locals().update(config)


rule copy_with_umi:
    input:
        fq_1 = lambda wildcards: config['replicate_label_to_fastq_1'][wildcards.replicate_label],
        fq_2 = lambda wildcards: config['replicate_label_to_fastq_2'][wildcards.replicate_label],
    output:
        fq_1 = temp("output/fastqs/copy/{replicate_label}-1.fastq.gz"), #SORT OUT!!
        fq_2 = temp("output/fastqs/copy/{replicate_label}-2.fastq.gz"), #SORT OUT!!        
    threads: 2
    resources:
        runtime="8h",
        mem_mb=8000
    benchmark: "benchmarks/umi/unassigned_experiment.{replicate_label}.copy_with_umi.txt"
    shell:
        "zcat {input.fq_1} | awk 'NR % 4 != 1 {{print}} NR % 4 == 1 {{split($1,header,\":\"); print $1 \":\" substr(header[1],2,length(header[1]) - 1) }}' | gzip > {output.fq_1};"
        "zcat {input.fq_2} | awk 'NR % 4 != 1 {{print}} NR % 4 == 1 {{split($1,header,\":\"); print $1 \":\" substr(header[1],2,length(header[1]) - 1) }}' | gzip > {output.fq_2};"

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
    container:
        "docker://howardxu520/skipper:fastqc_0.12.1"
    resources:
        mem_mb=16000,
        runtime="3h"
    params:
        outdir="output/fastqc/initial/"
    shell:
        "fastqc {input.r1} --extract --outdir {params.outdir} -t {threads};"
        "fastqc {input.r2} --extract --outdir {params.outdir} -t {threads};"
        
rule trim_fastq_encode:
    input:
        # fq_1 = lambda wildcards: replicate_label_to_fastq_1[wildcards.replicate_label],
        # fq_2 = lambda wildcards: replicate_label_to_fastq_2[wildcards.replicate_label],
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
        tmpdir=TMPDIR,
        mem_mb=16000,
        runtime="3h"
    benchmark: "benchmarks/trim/unassigned_experiment.{replicate_label}.trim.txt"
    container:
        "docker://howardxu520/skipper:skewer_0.2.2"
    shell:
        "skewer "
          "-t {threads} "
          "-x {input.adapter_1} "
          "-y {input.adapter_2} "
          "-o output/fastqs/trimmed/{wildcards.replicate_label} " 
          "-z -r 0.2 -d 0.2 -q 13 -l 20 "
          "{input.fq_1} {input.fq_2}"

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
    container:
        "docker://howardxu520/skipper:fastqc_0.12.1"
    resources:
        mem_mb=16000,
        runtime="4h"
    shell:
        "fastqc {input.r1} --extract --outdir output/fastqc/processed -t {threads};"
        "fastqc {input.r2} --extract --outdir output/fastqc/processed -t {threads};"
        
rule align_reads_encode:
    input:
        fq_1 = rules.trim_fastq_encode.output.fq_1_trimmed,
        fq_2 = rules.trim_fastq_encode.output.fq_2_trimmed
    output:
        ubam = temp("output/bams/raw/genome/{replicate_label}.genome.Aligned.out.bam"),
        # unmapped= "output/bams/raw/genome/{replicate_label}.genome.Unmapped.out.mate1",
        log= "output/bams/raw/genome/{replicate_label}.genome.Log.final.out",
    threads: 8
    params:
        star_sjdb = STAR_DIR,
        outprefix = "output/bams/raw/genome/{replicate_label}.genome.",
        rg = "{replicate_label}"
    benchmark: "benchmarks/align/unassigned_experiment.{replicate_label}.align_reads_genome.txt"
    container:
        "docker://howardxu520/skipper:star_2.7.10b"
    resources:
        mem_mb=40000,
        runtime="2h"
    shell:   
        "STAR "
            "--alignEndsType EndToEnd "
            "--genomeDir {params.star_sjdb} "
            "--genomeLoad NoSharedMemory "
            "--outBAMcompression 10 "
            "--outFileNamePrefix {params.outprefix} "
            "--winAnchorMultimapNmax 100 "
            "--outFilterMultimapNmax 100 "
            "--outFilterMultimapScoreRange 1 "
            "--outSAMmultNmax 1 "
            "--outMultimapperOrder Random "
            "--outFilterScoreMin 10 "
            "--outFilterType BySJout "
            "--limitOutSJcollapsed 5000000 "
            "--outReadsUnmapped None "
            "--outSAMattrRGline ID:{wildcards.replicate_label} "
            "--outSAMattributes All "
            "--outSAMmode Full "
            "--outSAMtype BAM Unsorted "
            "--outSAMunmapped Within "
            "--readFilesCommand zcat "
            "--outStd Log "
            "--readFilesIn {input.fq_1} {input.fq_2} "
            "--runMode alignReads "
            "--runThreadN {threads}"

rule sort_bam:
    input:
        bam="output/bams/raw/{ref}/{replicate_label}.{ref}.Aligned.out.bam",
    output:
        sort = "output/bams/raw/{ref}/{replicate_label}.{ref}.Aligned.sort.bam",
    threads: 2
    benchmark: "benchmarks/sort/{ref}/unassigned_experiment.{replicate_label}.sort_bam.txt"
    container:
        "docker://howardxu520/skipper:samtools_1.17"
    resources:
        mem_mb=16000,
        runteim="1h"
    shell:
        "samtools sort -T {wildcards.replicate_label} -@ {threads} -o {output.sort} {input.bam};"

rule index_bams:
    input:
        bam = "output/bams/{round}/{ref}/{replicate_label}.Aligned.{mid}.bam"
    output:
        ibam = "output/bams/{round}/{ref}/{replicate_label}.Aligned.{mid}.bam.bai"
    threads: 2
    benchmark: "benchmarks/index_bam/{round}/{ref}/{mid}/unassigned_experiment.{replicate_label}.index_bam.txt"
    container:
        "docker://howardxu520/skipper:samtools_1.17"
    resources:
        mem_mb=1000,
        runtime=20
    shell:
        "samtools index -@ {threads} {input.bam};"

rule dedup_umi_encode:
    input:
        bam="output/bams/raw/genome/{replicate_label}.genome.Aligned.sort.bam",
        ibam = "output/bams/raw/genome/{replicate_label}.genome.Aligned.sort.bam.bai"
    output:
        bam_dedup="output/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam"
    params:
        prefix='output/bams/dedup/genome/{replicate_label}.genome.sort'
    resources:
        mem_mb=34000,
        runtime="8h"
    benchmark: "benchmarks/dedup/genome/unassigned_experiment.{replicate_label}.dedup_umi.txt"
    container:
        "docker://howardxu520/skipper:umicollapse_1.0.0"
    resources:
        mem_mb=34000
    shell:
        "java -server -Xms32G -Xmx32G -Xss40M -jar /UMICollapse/umicollapse.jar bam "
            "-i {input.bam} -o {output.bam_dedup} --umi-sep : --two-pass"

rule select_informative_read:
    input:
        bam_combined="output/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam"
    output:
        bam_informative="output/bams/dedup/genome_R"+str(INFORMATIVE_READ)+"/{replicate_label}.genome.Aligned.sort.dedup.R"+str(INFORMATIVE_READ)+".bam"
    benchmark: "benchmarks/select/unassigned_experiment.{replicate_label}.select_informative_read.txt"
    container:
        "docker://howardxu520/skipper:samtools_1.17"
    resources:
        mem_mb=10000,
        runtime=30
    shell:
        "samtools view -bF " + str(64 if UNINFORMATIVE_READ == 1 else 128) + " {input.bam_combined} > {output.bam_informative}"

rule obtain_unique_reads:
    input:
        rules.select_informative_read.output.bam_informative
    output:
        temp("output/QC/{replicate_label}.uniq_fragments")
    benchmark:
        "benchmarks/{replicate_label}.count_uniq_fragments.txt"
    container:
        "docker://howardxu520/skipper:samtools_1.17"
    resources:
        mem_mb=10000,
        runtime=30
    shell:
        """
        samtools idxstats {input} | awk -F '\t' '{{s+=$3+$4}}END{{print s}}' > {output}
        """

rule obtain_aligned_reads:
    input:
        rules.align_reads_encode.output.ubam
    output:
        "output/QC/{replicate_label}.aligned_reads"
    benchmark:
        "benchmarks/{replicate_label}.count_aligned_reads.txt"
    container:
        "docker://howardxu520/skipper:samtools_1.17"
    resources:
        mem_mb=8000,
        runtime=30
    shell:
        """
        samtools idxstats {input} | awk -F '\t' '{{s+=$3+$4}}END{{print s}}' > {output}
        """
