locals().update(config) # JAVA_EXE
rule run_initial_fastqc:
    input:
        fq = lambda wildcards: config['replicate_label_to_fastqs'][wildcards.replicate_label].split(" "),
    output:
        report = "output/fastqc/initial/{replicate_label}_fastqc.html",
        zip_file = "output/fastqc/initial/{replicate_label}_fastqc.zip",
        directory = directory("output/fastqc/initial/{replicate_label}_fastqc")
    threads: 2
    params:
        error_file = "stderr/{replicate_label}.fastqc_initial.err",
        out_file = "stdout/{replicate_label}.fastqc_initial.out",
        run_time = "6:00:00",
        job_name = "run_initial_fastqc",
    benchmark: "benchmarks/fastqc/unassigned_experiment.{replicate_label}.initial_fastqc.txt"
    shell:
        "module load fastqc;"
        "less {input.fq} | fastqc stdin:{wildcards.replicate_label} --extract --outdir output/fastqc/initial -t {threads}"
        
rule trim_fastq:
    input:
        fq = lambda wildcards: config['replicate_label_to_fastqs'][wildcards.replicate_label].split(" "),
        adapter = lambda wildcards: config['replicate_label_to_adapter'][wildcards.replicate_label],
    output:
        fq_trimmed = temp("output/fastqs/trimmed/{replicate_label}-trimmed.fastq.gz"),
        metrics = "output/fastqs/trimmed/{replicate_label}-trimmed.log"
    threads: 8
    params:
        run_time = "5:30:00",
        memory = "30000",
        error_file = "stderr/{replicate_label}.trim.err",
        out_file = "stdout/{replicate_label}.trim.out",
        job_name = "trim_fastq"
    benchmark: "benchmarks/trim/unassigned_experiment.{replicate_label}.trim.txt"
    conda: "envs/skewer.yaml"
    shell:
        "less {input.fq} | skewer "
          "-t {threads} "
          "-x {input.adapter} "
          "-o output/fastqs/trimmed/{wildcards.replicate_label} "
          "-z -r 0.2 -d 0.2 -q 13 -l 20 -"

rule extract_umi:
    input:
        fq = rules.trim_fastq.output.fq_trimmed,
    output:
        fq_umi = "output/fastqs/umi/{replicate_label}.trimmed.umi.fq.gz",
        json = "output/fastp/{replicate_label}.fastp.json",
        html = "output/fastp/{replicate_label}.fastp.html",
    threads: 8
    params:
        error_file = "stderr/{replicate_label}.extract_umi.err",
        out_file = "stdout/{replicate_label}.extract_umi.out",
        run_time = "45:00",
        memory = "10000",
        job_name = "extract_umi",
        umi_length = config['UMI_SIZE'],
    benchmark: "benchmarks/umi/unassigned_experiment.{replicate_label}.extract_umi.txt"
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
    params:
        outdir="output/fastqc/processed/",
        run_time = "03:00:00",
        memory = "15000",
        error_file = "stderr/{replicate_label}.run_trimmed_fastqc.err",
        out_file = "stdout/{replicate_label}.run_trimmed_fastqc.out",
        job_name = "run_trimmed_fastqc"
    benchmark: "benchmarks/fastqc/unassigned_experiment.{replicate_label}.trimmed_fastqc.txt"
    shell:
        "module load fastqc;"
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
        error_file = "stderr/{replicate_label}.align_reads_genome.err",
        out_file = "stdout/{replicate_label}.align_reads_genome.out",
        run_time = "04:00:00",
        memory = "50000",
        job_name = "align_reads",
        star_sjdb = config['STAR_DIR'],
        outprefix = "output/bams/raw/genome/{replicate_label}.genome.",
        rg = "{replicate_label}"
    benchmark: "benchmarks/align/unassigned_experiment.{replicate_label}.align_reads_genome.txt"
    shell:  
        "module load star;"      
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
        sort = "output/bams/raw/{ref}/{replicate_label}.{ref}.Aligned.sort.bam",
    threads: 4
    params:
        error_file = "stderr/{ref}_{replicate_label}.sort_bam.err",
        out_file = "stdout/{ref}_{replicate_label}.sort_bam.out",
        run_time = "02:00:00",
        memory = "20000",
        job_name = "sortbam",
    benchmark: "benchmarks/sort/{ref}/unassigned_experiment.{replicate_label}.sort_bam.txt"
    shell:
        "set +eu;"
        "module load samtools/1.16;"
        "samtools sort -T {wildcards.replicate_label} -@ {threads} -o {output.sort} {input.bam};"
        

rule index_bams:
    input:
        bam = "output/bams/{round}/{ref}/{replicate_label}.Aligned.{mid}.bam"
    output:
        ibam = "output/bams/{round}/{ref}/{replicate_label}.Aligned.{mid}.bam.bai"
    threads: 2
    params:
        error_file = "stderr/{round}_{ref}_{mid}_{replicate_label}.index_bams.err",
        out_file = "stdout/{round}_{ref}_{mid}_{replicate_label}.index_bams.out",
        run_time = "10:00",
        memory = "1000",
        job_name = "index_bam"
    benchmark: "benchmarks/index_bam/{round}/{ref}/{mid}/unassigned_experiment.{replicate_label}.index_bam.txt"
    shell:
        "set +eu;"
        "module load samtools/1.16;"
        "samtools index -@ {threads} {input.bam};"


rule dedup_umi:
    input:
        bam="output/bams/raw/genome/{replicate_label}.genome.Aligned.sort.bam",
        ibam = "output/bams/raw/genome/{replicate_label}.genome.Aligned.sort.bam.bai"
    output:
        bam_dedup="output/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam"
    params:
        error_file = "stderr/{replicate_label}.dedup_umi.err",
        out_file = "stdout/{replicate_label}.dedup_umi.out",
        run_time = "12:00:00",
        memory = "60000",
        job_name = "dedup_bam",
        prefix='output/bams/dedup/genome/{replicate_label}.genome.sort'
    benchmark: "benchmarks/dedup/genome/unassigned_experiment.{replicate_label}.dedup_umi.txt"
    shell:
        "{JAVA_EXE} -server -Xms8G -Xmx8G -Xss20M -jar {UMICOLLAPSE_DIR}/umicollapse.jar bam "
            "-i {input.bam} -o {output.bam_dedup} --umi-sep : --two-pass"

rule obtain_unique_reads:
    input:
        rules.dedup_umi.output.bam_dedup
    output:
        "output/QC/{replicate_label}.uniq_fragments"
    params:
        error_file = "stderr/{replicate_label}.count_uniq_fragments.txt",
        out_file = "stdout/{replicate_label}.count_uniq_fragments.txt",
        run_time = "5:00",
        memory = "10000",
        job_name = "count_uniq_fragments",
    benchmark:
        "benchmarks/{replicate_label}.count_uniq_fragments.txt"
    shell:
        """
        module load samtools
        samtools idxstats {input} | awk -F '\t' '{{s+=$3+$4}}END{{print s}}' > {output}
        """