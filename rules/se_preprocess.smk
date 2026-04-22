import os
locals().update(config) 
if config.get('AGGRESSIVE_TRIM', False):
    skewer_k = 1
    skewer_m = 'any'
else:
    skewer_k = 2
    skewer_m = 'tail'

# Check to see if path already exists to ensure snakemake does not rebuild costly annotation files. 
star_exists = os.path.exists(CHROM_SIZES)

if not star_exists:
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

rule run_initial_fastqc:
    input:
        fq = lambda wildcards: config['replicate_label_to_fastqs'][wildcards.replicate_label].split(" "),
    output:
        report = "output/QC/fastqc/initial/{replicate_label}_fastqc.html",
        zip_file = "output/QC/fastqc/initial/{replicate_label}_fastqc.zip",
        directory = directory("output/QC/fastqc/initial/{replicate_label}_fastqc")
    threads: 2
    benchmark: "benchmarks/fastqc/unassigned_experiment.{replicate_label}.initial_fastqc.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{replicate_label}.run_initial_fastqc.out",
        stderr = config["WORKDIR"] + "/stderr/{replicate_label}.run_initial_fastqc.err",
    conda:
        "envs/fastqc.yaml"
    resources:
        mem_mb=16000,
        runtime="2h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting run_initial_fastqc for {wildcards.replicate_label}" | tee -a {log.stdout}

        zcat {input.fq} \
          | fastqc stdin:{wildcards.replicate_label} \
              --extract \
              --outdir output/QC/fastqc/initial \
              -t {threads} \
          >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished run_initial_fastqc for {wildcards.replicate_label}" | tee -a {log.stdout}
        """
        
rule trim_fastq:
    input:
        fq = lambda wildcards: config['replicate_label_to_fastqs'][wildcards.replicate_label].split(" "),
        adapter = lambda wildcards: config['replicate_label_to_adapter'][wildcards.replicate_label],
    output:
        fq_trimmed = temp("output/secondary_results/fastqs/trimmed/{replicate_label}-trimmed.fastq.gz"),
        metrics = "output/secondary_results/fastqs/trimmed/{replicate_label}-trimmed.log"
    threads: 8
    params:
        k = skewer_k,
        m = skewer_m
    benchmark: "benchmarks/trim/unassigned_experiment.{replicate_label}.trim.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{replicate_label}.trim_fastq.out",
        stderr = config["WORKDIR"] + "/stderr/{replicate_label}.trim_fastq.err",
    conda:
        "envs/skewer.yaml"
    resources:
        mem_mb=16000,
        runtime="2h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting trim_fastq for {wildcards.replicate_label}" | tee -a {log.stdout}

        zcat {input.fq} \
          | skewer \
              -t {threads} \
              -x {input.adapter} \
              -1 -r 0.2 -d 0.2 -q 13 -l 20 -k {params.k} -m {params.m} - \
          | skewer \
              -t {threads} \
              -x {input.adapter} \
              -o output/secondary_results/fastqs/trimmed/{wildcards.replicate_label} \
              -z -r 0.2 -d 0.2 -q 13 -l 20 -k {params.k} -m {params.m} - \
          >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished trim_fastq for {wildcards.replicate_label}" | tee -a {log.stdout}
        """

rule extract_umi:
    input:
        fq = rules.trim_fastq.output.fq_trimmed,
    output:
        fq_umi = "output/secondary_results/fastqs/umi/{replicate_label}.trimmed.umi.fq.gz",
        json = "output/QC/fastp/{replicate_label}.fastp.json",
        html = "output/QC/fastp/{replicate_label}.fastp.html",
    threads: 8
    params:
        umi_length = config['UMI_SIZE'],
    benchmark: "benchmarks/umi/unassigned_experiment.{replicate_label}.extract_umi.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{replicate_label}.extract_umi.out",
        stderr = config["WORKDIR"] + "/stderr/{replicate_label}.extract_umi.err",
    conda:
        "envs/fastp.yaml"
    resources:
        mem_mb=8000,
        runtime="1h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting extract_umi for {wildcards.replicate_label}" | tee -a {log.stdout}

        fastp \
            -i {input.fq} \
            -o {output.fq_umi} \
            -A \
            -U \
            --umi_len={params.umi_length} \
            --umi_loc=read1 \
            -j {output.json} \
            -h {output.html} \
            -w {threads} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished extract_umi for {wildcards.replicate_label}" | tee -a {log.stdout}
        """

rule run_trimmed_fastqc:
    input:
        rules.extract_umi.output.fq_umi,
    output:
        report = "output/QC/fastqc/processed/{replicate_label}.trimmed.umi_fastqc.html",
        zip_file = "output/QC/fastqc/processed/{replicate_label}.trimmed.umi_fastqc.zip",
    threads: 2
    benchmark: "benchmarks/fastqc/unassigned_experiment.{replicate_label}.trimmed_fastqc.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{replicate_label}.run_trimmed_fastqc.out",
        stderr = config["WORKDIR"] + "/stderr/{replicate_label}.run_trimmed_fastqc.err",
    conda:
        "envs/fastqc.yaml"
    resources:
        mem_mb=16000,
        runtime="3h"
    shell:
        r"""
        set -euo pipefail
        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting run_trimmed_fastqc for {wildcards.replicate_label}" | tee -a {log.stdout}

        fastqc {input} \
            --extract \
            --outdir output/QC/fastqc/processed \
            -t {threads} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished run_trimmed_fastqc for {wildcards.replicate_label}" | tee -a {log.stdout}
        """

rule align_reads:
    input:
        fq = rules.extract_umi.output.fq_umi,
        chrom_sizes = ancient(CHROM_SIZES),
    output:
        ubam = temp("output/secondary_results/bams/raw/genome/{replicate_label}.genome.Aligned.out.bam"),
        log_file = "output/secondary_results/bams/raw/genome/{replicate_label}.genome.Log.final.out",
    threads: 4
    params:
        star_sjdb = STAR_DIR,
        outprefix = "output/secondary_results/bams/raw/genome/{replicate_label}.genome.",
        rg = "{replicate_label}",
    benchmark: "benchmarks/align/unassigned_experiment.{replicate_label}.align_reads_genome.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{replicate_label}.align_reads.out",
        stderr = config["WORKDIR"] + "/stderr/{replicate_label}.align_reads.err",
    conda:
        "envs/star.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: 64000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 120 * (2 ** (attempt - 1)),
        tmpdir=TMPDIR,
        exclusive = True  # <-- tag for exclusivity
    shell:
        r"""
        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting align_reads for {wildcards.replicate_label}." | tee -a {log.stdout}

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
            --limitOutSJcollapsed 10000000 \
            --outReadsUnmapped None \
            --outSAMattrRGline ID:{wildcards.replicate_label} \
            --outSAMattributes All \
            --outSAMmode Full \
            --outSAMtype BAM Unsorted \
            --outSAMunmapped Within \
            --readFilesCommand zcat \
            --outStd Log \
            --readFilesIn {input.fq} \
            --runMode alignReads \
            --runThreadN {threads} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished align_reads for {wildcards.replicate_label}." | tee -a {log.stdout}
        """

rule sort_bam:
    input:
        bam="output/secondary_results/bams/raw/{ref}/{replicate_label}.{ref}.Aligned.out.bam",
    output:
        sort="output/secondary_results/bams/raw/{ref}/{replicate_label}.{ref}.Aligned.sort.bam",
    threads: 4
    benchmark: "benchmarks/sort/{ref}/unassigned_experiment.{replicate_label}.sort_bam.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{replicate_label}.{ref}.sort_bam.out",
        stderr = config["WORKDIR"] + "/stderr/{replicate_label}.{ref}.sort_bam.err",
    conda:
        "envs/bedbam_tools.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: 16000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 120 * (2 ** (attempt - 1)),
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting sort_bam for {wildcards.replicate_label}" | tee -a {log.stdout}

        samtools sort \
            -T {wildcards.replicate_label} \
            -@ {threads} \
            -o {output.sort} \
            {input.bam} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished sort_bam for {wildcards.replicate_label}" | tee -a {log.stdout}
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
        mem_mb=1000,
        runtime="30m"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting index_bams for {wildcards.replicate_label}" | tee -a {log.stdout}

        samtools index \
            -@ {threads} \
            {input.bam} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished index_bams for {wildcards.replicate_label}" | tee -a {log.stdout}
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
        mem_mb=lambda wildcards, attempt: 48000 * (1.5 ** (attempt - 1)),
        java_mem_mb=lambda wildcards, attempt: 32000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 120 * (2 ** (attempt - 1)),
        tmpdir = TMPDIR
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting dedup_umi for {wildcards.replicate_label}" | tee -a {log.stdout}

        JAR=$(dirname $(which umicollapse))/../share/umicollapse*/umicollapse.jar

        java -server -Xms{resources.java_mem_mb}M -Xmx{resources.java_mem_mb}M -Xss40M \
            -jar $JAR bam \
            -i {input.bam} \
            -o {output.bam_dedup} \
            --umi-sep : \
            --two-pass \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished dedup_umi for {wildcards.replicate_label}" | tee -a {log.stdout}
        """
