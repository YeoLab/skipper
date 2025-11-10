locals().update(config) 
if config.get('AGGRESSIVE_TRIM', False):
    skewer_k = 1
    skewer_m = 'any'
else:
    skewer_k = 2
    skewer_m = 'tail'

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
    log: "logs/run_star_genome_generate.log"
    conda:
        "envs/star.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting star_genome_generate." | tee {log}

        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {params.star_dir} \
            --genomeFastaFiles {input.fasta_file} \
            --sjdbGTFfile {input.gff} \
            --sjdbOverhang 99 \
            2>&1 | tee -a {log}

        echo "[`date`] Finished star_genome_generate." | tee -a {log}
        """

rule run_initial_fastqc:
    input:
        fq = lambda wildcards: config['replicate_label_to_fastqs'][wildcards.replicate_label].split(" "),
    output:
        report = "output/fastqc/initial/{replicate_label}_fastqc.html",
        zip_file = "output/fastqc/initial/{replicate_label}_fastqc.zip",
        directory = directory("output/fastqc/initial/{replicate_label}_fastqc")
    threads: 2
    benchmark: "benchmarks/fastqc/unassigned_experiment.{replicate_label}.initial_fastqc.txt"
    log: "logs/{replicate_label}.run_initial_fastqc.log"
    conda:
        "envs/fastqc.yaml"
    resources:
        mem_mb=16000,
        runtime="2h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting run_initial_fastqc for {wildcards.replicate_label}" | tee {log}

        zcat {input.fq} \
          | fastqc stdin:{wildcards.replicate_label} \
              --extract \
              --outdir output/fastqc/initial \
              -t {threads} \
          2>&1 | tee -a {log}

        echo "[`date`] Finished run_initial_fastqc for {wildcards.replicate_label}" | tee -a {log}
        """
        
rule trim_fastq:
    input:
        fq = lambda wildcards: config['replicate_label_to_fastqs'][wildcards.replicate_label].split(" "),
        adapter = lambda wildcards: config['replicate_label_to_adapter'][wildcards.replicate_label],
    output:
        fq_trimmed = temp("output/fastqs/trimmed/{replicate_label}-trimmed.fastq.gz"),
        metrics = "output/fastqs/trimmed/{replicate_label}-trimmed.log"
    threads: 8
    params:
        k = skewer_k,
        m = skewer_m
    benchmark: "benchmarks/trim/unassigned_experiment.{replicate_label}.trim.txt"
    log: "logs/{replicate_label}.trim_fastq.log"
    conda:
        "envs/skewer.yaml"
    resources:
        mem_mb=16000,
        runtime="2h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting trim_fastq for {wildcards.replicate_label}" | tee {log}

        zcat {input.fq} \
          | skewer \
              -t {threads} \
              -x {input.adapter} \
              -1 -r 0.2 -d 0.2 -q 13 -l 20 -k {params.k} -m {params.m} - \
          | skewer \
              -t {threads} \
              -x {input.adapter} \
              -o output/fastqs/trimmed/{wildcards.replicate_label} \
              -z -r 0.2 -d 0.2 -q 13 -l 20 -k {params.k} -m {params.m} - \
          2>&1 | tee -a {log}

        echo "[`date`] Finished trim_fastq for {wildcards.replicate_label}" | tee -a {log}
        """

rule extract_umi:
    input:
        fq = rules.trim_fastq.output.fq_trimmed,
    output:
        fq_umi = "output/fastqs/umi/{replicate_label}.trimmed.umi.fq.gz",
        json = "output/fastp/{replicate_label}.fastp.json",
        html = "output/fastp/{replicate_label}.fastp.html",
    threads: 8
    params:
        umi_length = config['UMI_SIZE'],
    benchmark: "benchmarks/umi/unassigned_experiment.{replicate_label}.extract_umi.txt"
    log: "logs/{replicate_label}.extract_umi.log"
    conda:
        "envs/fastp.yaml"
    resources:
        mem_mb=8000,
        runtime="1h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting extract_umi for {wildcards.replicate_label}" | tee {log}

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
            2>&1 | tee -a {log}

        echo "[`date`] Finished extract_umi for {wildcards.replicate_label}" | tee -a {log}
        """

rule run_trimmed_fastqc:
    input:
        rules.extract_umi.output.fq_umi,
    output:
        report = "output/fastqc/processed/{replicate_label}.trimmed.umi_fastqc.html",
        zip_file = "output/fastqc/processed/{replicate_label}.trimmed.umi_fastqc.zip",
    threads: 2
    benchmark: "benchmarks/fastqc/unassigned_experiment.{replicate_label}.trimmed_fastqc.txt"
    log: "logs/{replicate_label}.run_trimmed_fastqc.log"
    conda:
        "envs/fastqc.yaml"
    resources:
        mem_mb=16000,
        runtime="3h"
    shell:
        r"""
        set -euo pipefail
        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting run_trimmed_fastqc for {wildcards.replicate_label}" | tee {log}

        fastqc {input} \
            --extract \
            --outdir output/fastqc/processed \
            -t {threads} \
            2>&1 | tee -a {log}

        echo "[`date`] Finished run_trimmed_fastqc for {wildcards.replicate_label}" | tee -a {log}
        """

rule align_reads:
    input:
        fq = rules.extract_umi.output.fq_umi,
        chrom_sizes = CHROM_SIZES,
    output:
        ubam = temp("output/bams/raw/genome/{replicate_label}.genome.Aligned.out.bam"),
        log_file = "output/bams/raw/genome/{replicate_label}.genome.Log.final.out",
    threads: 4
    params:
        star_sjdb = STAR_DIR,
        outprefix = "output/bams/raw/genome/{replicate_label}.genome.",
        rg = "{replicate_label}",
    benchmark: "benchmarks/align/unassigned_experiment.{replicate_label}.align_reads_genome.txt"
    log: "logs/{replicate_label}.align_reads.log"
    conda:
        "envs/star.yaml"
    resources:
        mem_mb=64000,
        runtime="2h",
        tmpdir=TMPDIR,
        exclusive = True  # <-- tag for exclusivity
    shell:
        r"""
        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting align_reads for {wildcards.replicate_label}." | tee -a {log}

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
            --runThreadN {threads}

        echo "[`date`] Finished align_reads for {wildcards.replicate_label}." | tee -a {log}
        """

rule sort_bam:
    input:
        bam="output/bams/raw/{ref}/{replicate_label}.{ref}.Aligned.out.bam",
    output:
        sort="output/bams/raw/{ref}/{replicate_label}.{ref}.Aligned.sort.bam",
    threads: 4
    benchmark: "benchmarks/sort/{ref}/unassigned_experiment.{replicate_label}.sort_bam.txt"
    log: "logs/{replicate_label}.{ref}.sort_bam.log"
    conda:
        "envs/bedbam_tools.yaml"
    resources:
        mem_mb=16000,
        runtime="2h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting sort_bam for {wildcards.replicate_label}" | tee {log}

        samtools sort \
            -T {wildcards.replicate_label} \
            -@ {threads} \
            -o {output.sort} \
            {input.bam} \
            2>&1 | tee -a {log}

        echo "[`date`] Finished sort_bam for {wildcards.replicate_label}" | tee -a {log}
        """

rule index_bams:
    input:
        bam = "output/bams/{round}/{ref}/{replicate_label}.Aligned.{mid}.bam"
    output:
        ibam = "output/bams/{round}/{ref}/{replicate_label}.Aligned.{mid}.bam.bai"
    threads: 2
    benchmark: "benchmarks/index_bam/{round}/{ref}/{mid}/unassigned_experiment.{replicate_label}.index_bam.txt"
    log: "logs/{replicate_label}.{round}.{ref}.{mid}.index_bams.log"
    conda:
        "envs/bedbam_tools.yaml"
    resources:
        mem_mb=1000,
        runtime="30m"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting index_bams for {wildcards.replicate_label}" | tee {log}

        samtools index \
            -@ {threads} \
            {input.bam} \
            2>&1 | tee -a {log}

        echo "[`date`] Finished index_bams for {wildcards.replicate_label}" | tee -a {log}
        """

rule dedup_umi:
    input:
        bam = "output/bams/raw/genome/{replicate_label}.genome.Aligned.sort.bam",
        ibam = "output/bams/raw/genome/{replicate_label}.genome.Aligned.sort.bam.bai"
    output:
        bam_dedup = "output/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam"
    benchmark: "benchmarks/dedup/genome/unassigned_experiment.{replicate_label}.dedup_umi.txt"
    log: "logs/{replicate_label}.dedup_umi.log"
    conda:
        "envs/umicollapse.yaml"
    resources:
        mem_mb = 48000,
        runtime = "2h",
        tmpdir = TMPDIR
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting dedup_umi for {wildcards.replicate_label}" | tee {log}

        JAR=$(dirname $(which umicollapse))/../share/umicollapse*/umicollapse.jar

        java -server -Xms32G -Xmx32G -Xss40M \
            -jar $JAR bam \
            -i {input.bam} \
            -o {output.bam_dedup} \
            --umi-sep : \
            --two-pass \
            2>&1 | tee -a {log}

        echo "[`date`] Finished dedup_umi for {wildcards.replicate_label}" | tee -a {log}
        """

rule obtain_unique_reads:
    input:
        rules.dedup_umi.output.bam_dedup
    output:
        "output/QC/{replicate_label}.uniq_fragments"
    benchmark:
        "benchmarks/{replicate_label}.obtain_unique_reads.txt"
    log: 
        "logs/{replicate_label}.obtain_unique_reads.log"
    conda:
        "envs/bedbam_tools.yaml"
    resources:
        mem_mb=8000,
        runtime="1h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting obtain_unique_reads for {wildcards.replicate_label}" | tee {log}

        samtools idxstats {input} \
          | awk -F '\t' '{{s+=$3+$4}} END {{print s}}' \
          > {output} \
          2>&1 | tee -a {log}

        echo "[`date`] Finished obtain_unique_reads for {wildcards.replicate_label}" | tee -a {log}
        """
