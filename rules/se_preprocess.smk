locals().update(config) 
if config.get('AGGRESSIVE_TRIM', False):
    skewer_k = 1
    skewer_m = 'any'
else:
    skewer_k = 2
    skewer_m = 'tail'

rule run_initial_fastqc:
    input:
        fq = lambda wildcards: config['replicate_label_to_fastqs'][wildcards.replicate_label].split(" "),
    output:
        report = "output/fastqc/initial/{replicate_label}_fastqc.html",
        zip_file = "output/fastqc/initial/{replicate_label}_fastqc.zip",
        directory = directory("output/fastqc/initial/{replicate_label}_fastqc")
    threads: 2
    benchmark: "benchmarks/fastqc/unassigned_experiment.{replicate_label}.initial_fastqc.txt"
    container:
        "docker://howardxu520/skipper:fastqc_0.12.1"
    resources:
        mem_mb=16000,
        runtime="6h"
    shell:
        "zcat {input.fq} | fastqc stdin:{wildcards.replicate_label} --extract --outdir output/fastqc/initial -t {threads}"
        
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
    container:
        "docker://howardxu520/skipper:skewer_0.2.2"
    resources:
        mem_mb=16000,
        runtime="6h"
    shell:
        "zcat {input.fq} | skewer "
          "-t {threads} "
          "-x {input.adapter} "
          "-1 -r 0.2 -d 0.2 -q 13 -l 20 -k {params.k} -m {params.m} - "
          "| skewer -t {threads} "
          "-x {input.adapter} "
          "-o output/fastqs/trimmed/{wildcards.replicate_label} "
          "-z -r 0.2 -d 0.2 -q 13 -l 20 -k {params.k} -m {params.m} -"

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
    container:
        "docker://howardxu520/skipper:fastp_0.23.4"
    resources:
        mem_mb=8000,
        runtime=45
    shell:      
        "fastp "
            "-i {input.fq} "
            "-o {output.fq_umi} "
            "-A "
            "-U "
            "--umi_len={params.umi_length} "
            "--umi_loc=read1 "
            "-j output/fastp/{wildcards.replicate_label}.fastp.json "
            "-h output/fastp/{wildcards.replicate_label}.fastp.html "
            "-w {threads}"


rule run_trimmed_fastqc:
    input:
        rules.extract_umi.output.fq_umi,
    output:
        report = "output/fastqc/processed/{replicate_label}.trimmed.umi_fastqc.html",
        zip_file = "output/fastqc/processed/{replicate_label}.trimmed.umi_fastqc.zip",
    threads: 2
    benchmark: "benchmarks/fastqc/unassigned_experiment.{replicate_label}.trimmed_fastqc.txt"
    container:
        "docker://howardxu520/skipper:fastqc_0.12.1"
    resources:
        mem_mb=16000,
        runtime="3h"
    shell:
        "fastqc {input} --extract --outdir output/fastqc/processed -t {threads}"
        
rule align_reads:
    input:
        fq= rules.extract_umi.output.fq_umi,
    output:
        ubam = temp("output/bams/raw/genome/{replicate_label}.genome.Aligned.out.bam"),
        # unmapped= "output/bams/raw/genome/{replicate_label}.genome.Unmapped.out.mate1",
        log= "output/bams/raw/genome/{replicate_label}.genome.Log.final.out",
    threads: 8
    params:
        star_sjdb = config['STAR_DIR'],
        outprefix = "output/bams/raw/genome/{replicate_label}.genome.",
        rg = "{replicate_label}"
    benchmark: "benchmarks/align/unassigned_experiment.{replicate_label}.align_reads_genome.txt"
    container:
        "docker://howardxu520/skipper:star_2.7.10b"
    resources:
        mem_mb=40000,
        runtime=240
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
            "--limitOutSJcollapsed 10000000 "
            "--outReadsUnmapped None "
            "--outSAMattrRGline ID:{wildcards.replicate_label} "
            "--outSAMattributes All "
            "--outSAMmode Full "
            "--outSAMtype BAM Unsorted "
            "--outSAMunmapped Within "
            "--readFilesCommand zcat "
            "--outStd Log "
            "--readFilesIn {input.fq} "
            "--runMode alignReads "
            "--runThreadN {threads}"
        
rule sort_bam:
    input:
        bam="output/bams/raw/{ref}/{replicate_label}.{ref}.Aligned.out.bam",
    output:
        sort="output/bams/raw/{ref}/{replicate_label}.{ref}.Aligned.sort.bam",
    threads: 4
    benchmark: "benchmarks/sort/{ref}/unassigned_experiment.{replicate_label}.sort_bam.txt"
    container:
        "docker://howardxu520/skipper:samtools_1.17"
    resources:
        mem_mb=16000,
        runtime="2h"
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
        runtime=10
    shell:
        "samtools index -@ {threads} {input.bam};"


rule dedup_umi:
    input:
        bam="output/bams/raw/genome/{replicate_label}.genome.Aligned.sort.bam",
        ibam = "output/bams/raw/genome/{replicate_label}.genome.Aligned.sort.bam.bai"
    output:
        bam_dedup="output/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam",
    benchmark: "benchmarks/dedup/genome/unassigned_experiment.{replicate_label}.dedup_umi.txt"
    container:
        "docker://howardxu520/skipper:umicollapse_1.0.0"
    resources:
        mem_mb=34000,
        runtime="3h",
        tmpdir=TMPDIR
    shell:
        "java -server -Xms32G -Xmx32G -Xss40M -Djava.io.tmpdir={resources.tmpdir} -jar /UMICollapse/umicollapse.jar bam "
            "-i {input.bam} -o {output.bam_dedup} --umi-sep : --two-pass"

rule obtain_unique_reads:
    input:
        rules.dedup_umi.output.bam_dedup
    output:
        "output/QC/{replicate_label}.uniq_fragments"
    benchmark:
        "benchmarks/{replicate_label}.count_uniq_fragments.txt"
    container:
        "docker://howardxu520/skipper:samtools_1.17"
    resources:
        mem_mb=8000,
        runtime=30
    shell:
        """
        samtools idxstats {input} | awk -F '\t' '{{s+=$3+$4}}END{{print s}}' > {output}
        """
