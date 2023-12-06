import pandas as pd
from functools import reduce
import re
import os
import sys
import glob
from time import sleep
locals().update(config)
# example command
# snakemake -kps encode_pe_rules_HepG2.py -w 25 -j 30 --cluster "qsub -e {params.error_file} -o {params.out_file} -l walltime={params.run_time} -l nodes=1:ppn={threads} -q home-yeo" >> encode_pe_rules_HepG2.log 2>&1

# include: "encode_pe_rules_config_HepG2_20230620.py"

# if not os.path.exists("stderr"): os.makedirs("stderr")
# if not os.path.exists("stdout"): os.makedirs("stdout")

# if EXE_DIR not in sys.path: os.environ["PATH"] = EXE_DIR + os.pathsep + os.environ["PATH"]

# if OVERDISPERSION_MODE not in ["clip","input"]:
#     raise Exception("Overdispersion must be calculated using 'clip' or 'input' samples")

# manifest = pd.read_csv(MANIFEST, index_col = False)
# manifest["Input_replicate_label"] = [str(sample) + "_IN_" + str(replicate) for replicate, sample in zip(manifest.Input_replicate.tolist(),manifest.Sample.tolist())]
# manifest["CLIP_replicate_label"] = [str(sample) + "_IP_" + str(replicate) for replicate, sample in zip(manifest.CLIP_replicate.tolist(),manifest.Sample.tolist())]

# input_replicates = manifest.loc[:,manifest.columns.isin(["Input_replicate_label","Input_fastq_1","Input_fastq_2","Input_bam","Input_adapter_1", "Input_adapter_2"])].drop_duplicates()
# clip_replicates = manifest.loc[:,manifest.columns.isin(["CLIP_replicate_label","CLIP_fastq_1","CLIP_fastq_2","CLIP_bam","CLIP_adapter_1", "CLIP_adapter_2"])].drop_duplicates()


# if len(input_replicates) != len(input_replicates[["Input_replicate_label"]].drop_duplicates()) or \
#     len(clip_replicates) != len(clip_replicates[["CLIP_replicate_label"]].drop_duplicates()):
#     raise Exception("Manifest files are not consistent across replicates")


# input_replicate_labels = input_replicates.Input_replicate_label.tolist()
# clip_replicate_labels = clip_replicates.CLIP_replicate_label.tolist()
# replicate_labels = pd.Series(input_replicate_labels + clip_replicate_labels)

# if all(bam in manifest.columns.tolist() for bam in ["Input_bam", "CLIP_bam"]):
#     replicate_label_to_bams = dict(zip(input_replicate_labels + clip_replicate_labels, input_replicate_labels.Input_bam.tolist() + clip_replicate_labels.CLIP_bam.tolist()))    
# else:
#     replicate_label_to_bams = dict(zip(input_replicate_labels + clip_replicate_labels, ["output/bams/dedup/genome/" + replicate_label + ".genome.Aligned.sort.dedup.bam" for replicate_label in input_replicate_labels + clip_replicate_labels] ))




# rule all:
#     input:
#         expand("output/fastqc/initial/{replicate_label}-1_fastqc.html", replicate_label = replicate_labels), 
#         expand("output/fastqc/initial/{replicate_label}-2_fastqc.html", replicate_label = replicate_labels), 
#         expand("output/fastqc/processed/{replicate_label}-trimmed-pair1_fastqc.html", replicate_label = replicate_labels), 
#         expand("output/fastqc/processed/{replicate_label}-trimmed-pair2_fastqc.html", replicate_label = replicate_labels),
#         expand("output/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam", replicate_label = replicate_labels), 
#         expand("output/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam.bai", replicate_label = replicate_labels), 
#         expand("output/bams/dedup/genome_R"+str(INFORMATIVE_READ)+"/{replicate_label}.genome.Aligned.sort.dedup.R"+str(INFORMATIVE_READ)+".bam", replicate_label = replicate_labels)
#     output:
#         "encode_prep.txt"
#     threads: 1
#     params:
#         error_file = "stderr/all.err",
#         out_file = "stdout/all.out",
#         run_time = "00:04:00",
#         memory = "20",
#         job_name = "all"
#     shell:
#         "echo $(date)  > {output};"
#         "echo created by Evan Boyle and the Yeo lab >> {output}"


rule copy_with_umi:
    input:
        fq_1 = lambda wildcards: config['replicate_label_to_fastq_1'][wildcards.replicate_label],
        fq_2 = lambda wildcards: config['replicate_label_to_fastq_2'][wildcards.replicate_label],
    output:
        fq_1 = temp("output/fastqs/copy/{replicate_label}-1.fastq.gz"), #SORT OUT!!
        fq_2 = temp("output/fastqs/copy/{replicate_label}-2.fastq.gz"), #SORT OUT!!        
    threads: 2
    params:
        run_time = "6:00:00",
        error_file = "stderr/{replicate_label}.copy_with_umi.err",
        out_file = "stdout/{replicate_label}.copy_with_umi.out",
        job_name = "copy_with_umi"
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
    params:
        outdir="output/fastqc/initial/",
        run_time = "03:00:00",
        memory = "15000",
        error_file = "stderr/{replicate_label}.run_initial_fastqc.err",
        out_file = "stdout/{replicate_label}.run_initial_fastqc.out",
        job_name = "run_initial_fastqc"
    benchmark: "benchmarks/fastqc/unassigned_experiment.{replicate_label}.initial_fastqc.txt"
    shell:
        "module load fastqc;"
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
    params:
        run_time = "3:30:00",
        memory = "15000",
        error_file = "stderr/{replicate_label}.trim.err",
        out_file = "stdout/{replicate_label}.trim.out",
        job_name = "trim_fastq"
    benchmark: "benchmarks/trim/unassigned_experiment.{replicate_label}.trim.txt"
    conda: "envs/skewer.yaml"
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
        error_file = "stderr/{replicate_label}.align_reads_genome.err",
        out_file = "stdout/{replicate_label}.align_reads_genome.out",
        run_time = "02:00:00",
        memory = "40000",
        job_name = "align_reads",
        star_sjdb = STAR_DIR,
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
    params:
        error_file = "stderr/{ref}_{replicate_label}.sort_bam.err",
        out_file = "stdout/{ref}_{replicate_label}.sort_bam.out",
        run_time = "00:30:00",
        memory = "10000",
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
        run_time = "20:00",
        memory = "1000",
        job_name = "index_bam"
    benchmark: "benchmarks/index_bam/{round}/{ref}/{mid}/unassigned_experiment.{replicate_label}.index_bam.txt"
    shell:
        "set +eu;"
        "module load samtools/1.16;"
        "samtools index -@ {threads} {input.bam};"

rule dedup_umi_encode:
    input:
        bam="output/bams/raw/genome/{replicate_label}.genome.Aligned.sort.bam",
        ibam = "output/bams/raw/genome/{replicate_label}.genome.Aligned.sort.bam.bai"
    output:
        bam_dedup="output/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam"
    params:
        error_file = "stderr/{replicate_label}.dedup_umi.err",
        out_file = "stdout/{replicate_label}.dedup_umi.out",
        run_time = "8:00:00",
        memory = "32000",
        job_name = "dedup_bam",
        prefix='output/bams/dedup/genome/{replicate_label}.genome.sort'
    benchmark: "benchmarks/dedup/genome/unassigned_experiment.{replicate_label}.dedup_umi.txt"
    shell:
        "{JAVA_EXE} -server -Xms8G -Xmx8G -Xss20M -jar {UMICOLLAPSE_DIR}/umicollapse.jar bam "
            "-i {input.bam} -o {output.bam_dedup} --umi-sep : --two-pass"

rule select_informative_read:
    input:
        bam_combined="output/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam"
    output:
        bam_informative="output/bams/dedup/genome_R"+str(INFORMATIVE_READ)+"/{replicate_label}.genome.Aligned.sort.dedup.R"+str(INFORMATIVE_READ)+".bam"
    params:
        error_file = "stderr/{replicate_label}.select_informative_read.err",
        out_file = "stdout/{replicate_label}.select_informative_read.out",
        run_time = "10:00",
        memory = "10000",
        job_name = "select_informative_read",
        # prefix='output/bams/dedup/genome/{replicate_label}.genome.sort'
    benchmark: "benchmarks/select/unassigned_experiment.{replicate_label}.select_informative_read.txt"
    shell:
        "set +eu;"
        "module load samtools/1.16;"
        "samtools view -bF " + str(64 if UNINFORMATIVE_READ == 1 else 128) + " {input.bam_combined} > {output.bam_informative}"

rule obtain_unique_reads:
    input:
        rules.select_informative_read.output.bam_informative
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