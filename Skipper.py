import pandas as pd
from functools import reduce
import os
import sys
import glob

# example command
# snakemake -kps Skipper.py -w 15 -j 30 --cluster "qsub -e {params.error_file} -o {params.out_file} -l walltime={params.run_time} -l nodes=1:ppn={threads} -q home-yeo" 

include: "Skipper_config.py"

if not os.path.exists("stderr"): os.makedirs("stderr")
if not os.path.exists("stdout"): os.makedirs("stdout")

if EXE_DIR not in sys.path: os.environ["PATH"] = EXE_DIR + os.pathsep + os.environ["PATH"]

if OVERDISPERSION_MODE not in ["clip","input"]:
    raise Exception("Overdispersion must be calculated using 'clip' or 'input' samples")

manifest = pd.read_csv(MANIFEST, comment = "#", index_col = False)
manifest["Input_replicate_label"] = [str(sample) + "_IN_" + str(replicate) for replicate, sample in zip(manifest.Input_replicate.tolist(),manifest.Sample.tolist())]
manifest["CLIP_replicate_label"] = [str(sample) + "_IP_" + str(replicate) for replicate, sample in zip(manifest.CLIP_replicate.tolist(),manifest.Sample.tolist())]

input_replicates = manifest.loc[:,manifest.columns.isin(["Input_replicate_label","Input_fastq","Input_bam","Input_adapter"])].drop_duplicates()
clip_replicates = manifest.loc[:,manifest.columns.isin(["CLIP_replicate_label","CLIP_fastq","CLIP_bam","CLIP_adapter"])].drop_duplicates()

if len(input_replicates) != len(input_replicates[["Input_replicate_label"]].drop_duplicates()) or \
    len(clip_replicates) != len(clip_replicates[["CLIP_replicate_label"]].drop_duplicates()):
    raise Exception("Manifest files are not consistent across replicates")

input_replicate_labels = input_replicates.Input_replicate_label.tolist()
clip_replicate_labels = clip_replicates.CLIP_replicate_label.tolist()
replicate_labels = pd.Series(input_replicate_labels + clip_replicate_labels)

if all(bam in manifest.columns.tolist() for bam in ["Input_bam", "CLIP_bam"]):
    replicate_label_to_bams = dict(zip(input_replicate_labels + clip_replicate_labels, input_replicates.Input_bam.tolist() + clip_replicates.CLIP_bam.tolist()))    
else:
    replicate_label_to_bams = dict(zip(input_replicate_labels + clip_replicate_labels, ["output/bams/dedup/genome/" + replicate_label + ".genome.Aligned.sort.dedup.bam" for replicate_label in input_replicate_labels + clip_replicate_labels] ))

experiment_labels = pd.Series(manifest.Experiment.drop_duplicates().tolist())
experiment_data = manifest.groupby("Experiment").agg({"CLIP_replicate_label": list, "Input_replicate_label" : list})

if "Input_fastq" in manifest.columns:
    replicate_label_to_fastqs = dict(zip(input_replicate_labels + clip_replicate_labels, input_replicates.Input_fastq.tolist() + clip_replicates.CLIP_fastq.tolist()))
    replicate_label_to_adapter = dict(zip(input_replicate_labels + clip_replicate_labels, input_replicates.Input_adapter.tolist() + clip_replicates.CLIP_adapter.tolist()))

overdispersion_replicate_lookup = dict(zip(manifest.CLIP_replicate_label.tolist(), manifest.Input_replicate_label.tolist() if OVERDISPERSION_MODE == "input" else manifest.CLIP_replicate_label.tolist()))
clip_to_input_replicate_label = dict(zip(manifest.CLIP_replicate_label.tolist(), manifest.Input_replicate_label.tolist()))
experiment_to_replicate_labels = dict(zip(experiment_data.index.tolist(), [reduce(lambda agg, x: agg if x in agg else agg + [x], inputs, []) + clips for inputs, clips in zip(experiment_data.Input_replicate_label, experiment_data.CLIP_replicate_label)]))
experiment_to_clip_replicate_labels = dict(zip(experiment_data.index.tolist(), experiment_data.CLIP_replicate_label))

experiment_to_input_replicate_labels = {}
for experiment_label, label_list in zip(experiment_data.index, experiment_data.Input_replicate_label):
    experiment_to_input_replicate_labels[experiment_label] = {}
    for entry in label_list:
        replicates = set()
        for other_entry in label_list:
            if other_entry != entry:
                replicates.add(other_entry)
        experiment_to_input_replicate_labels[experiment_label].update({entry : list(replicates)})

rule all:
    input:
        expand("output/fastqc/initial/{replicate_label}_fastqc.html", replicate_label = replicate_labels), 
        expand("output/fastqc/processed/{replicate_label}.trimmed.umi_fastqc.html", replicate_label = replicate_labels), 
        expand("output/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam", replicate_label = replicate_labels), 
        expand("output/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam.bai", replicate_label = replicate_labels), 
        # expand("output/bigwigs/plus/{replicate_label}.plus.bw", replicate_label = replicate_labels),
        expand("output/counts/repeats/vectors/{replicate_label}.counts", replicate_label = replicate_labels),
        expand("output/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_windows.tsv.gz", zip, experiment_label = manifest.Experiment, clip_replicate_label = manifest.CLIP_replicate_label),
        expand("output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz", experiment_label = manifest.Experiment),
        expand("output/figures/enrichment_reproducibility/{experiment_label}.enrichment_reproducibility.pdf", experiment_label = manifest.Experiment),
        expand("output/counts/repeats/tables/family/{experiment_label}.tsv.gz", experiment_label = manifest.Experiment),
        expand("output/reproducible_enriched_re/{experiment_label}.reproducible_enriched_re.tsv.gz", experiment_label = manifest.Experiment),
        expand("output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz", experiment_label = manifest.Experiment),
        expand("output/homer/finemapped_results/{experiment_label}/homerResults.html", experiment_label = manifest.Experiment),
        expand("output/gene_sets/{experiment_label}.enriched_terms.tsv.gz", experiment_label = manifest.Experiment),
        "output/figures/tsne/skipper.tsne_query.pdf",
    output:
        "land_ho.txt"
    threads: 1
    params:
        error_file = "stderr/all.err",
        out_file = "stdout/all.out",
        run_time = "00:04:00",
        memory = "20",
        job_name = "all"
    shell:
        "echo $(date)  > {output};"
        "echo created by Evan Boyle and the Yeo lab >> {output}"

rule parse_gff:
    input:
        gff = GFF,
        rankings = ACCESSION_RANKINGS,
    output:
        partition = PARTITION,
        feature_annotations = FEATURE_ANNOTATIONS,
    threads: 4
    params:
        error_file = "stderr/parse_gff.err",
        out_file = "stdout/parse_gff.out",
        run_time = "3:00:00",
        job_name = "parse_gff"
    benchmark: "benchmarks/parse_gff.txt"
    shell:        
        "{R_EXE} --vanilla {TOOL_DIR}/parse_gff.R {input.gff} {input.rankings} {output.partition} {output.feature_annotations}"

rule run_initial_fastqc:
    input:
        fq = lambda wildcards: replicate_label_to_fastqs[wildcards.replicate_label].split(" "),
    output:
        report = "output/fastqc/initial/{replicate_label}_fastqc.html",
        zip_file = "output/fastqc/initial/{replicate_label}_fastqc.zip",
        directory = directory("output/fastqc/initial/{replicate_label}_fastqc")
    threads: 1
    params:
        error_file = "stderr/{replicate_label}.fastqc_initial.err",
        out_file = "stdout/{replicate_label}.fastqc_initial.out",
        run_time = "3:00:00",
        job_name = "run_initial_fastqc"
    benchmark: "benchmarks/fastqc/unassigned_experiment.{replicate_label}.initial_fastqc.txt"
    shell:
        "module load fastqc;"
        "less {input.fq} | fastqc stdin:{wildcards.replicate_label} --extract --outdir output/fastqc/initial -t {threads}"
        
rule trim_fastq:
    input:
        fq = lambda wildcards: replicate_label_to_fastqs[wildcards.replicate_label].split(" "),
        adapter = lambda wildcards: replicate_label_to_adapter[wildcards.replicate_label],
    output:
        fq_trimmed = temp("output/fastqs/trimmed/{replicate_label}-trimmed.fastq.gz"),
        metrics = "output/fastqs/trimmed/{replicate_label}-trimmed.log"
    threads: 8
    params:
        run_time = "3:30:00",
        memory = "15000",
        error_file = "stderr/{replicate_label}.trim.err",
        out_file = "stdout/{replicate_label}.trim.out",
        job_name = "trim_fastq"
    benchmark: "benchmarks/trim/unassigned_experiment.{replicate_label}.trim.txt"
    shell:
        "less {input.fq} | skewer "
          "-t {threads} "
          "-x {input.adapter} "
          "-o output/fastqs/trimmed/{wildcards.replicate_label} "
          "-z -r 0.2 -d 0.2 -q 13 -l 20 -"

rule extract_umi:
    input:
        fq = "output/fastqs/trimmed/{replicate_label}-trimmed.fastq.gz",
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
        umi_length = UMI_SIZE,
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
        "output/fastqs/umi/{replicate_label}.trimmed.umi.fq.gz",
    output:
        report = "output/fastqc/processed/{replicate_label}.trimmed.umi_fastqc.html",
        zip_file = "output/fastqc/processed/{replicate_label}.trimmed.umi_fastqc.zip",
    threads: 1
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
        fq= "output/fastqs/umi/{replicate_label}.trimmed.umi.fq.gz",
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
            "--readFilesIn {input.fq} "
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
        run_time = "1:00:00",
        memory = "10000",
        job_name = "dedup_bam",
        prefix='output/bams/dedup/genome/{replicate_label}.genome.sort'
    benchmark: "benchmarks/dedup/genome/unassigned_experiment.{replicate_label}.dedup_umi.txt"
    shell:
        "{JAVA_EXE} -server -Xms8G -Xmx8G -Xss20M -jar {UMICOLLAPSE_DIR}/umicollapse.jar bam "
            "-i {input.bam} -o {output.bam_dedup} --umi-sep : --two-pass"

rule make_bigwig:
    input:
        CHROM_SIZES,
        bam = lambda wildcards: replicate_label_to_bams[wildcards.replicate_label],
    output:
        bg_plus = temp("output/bedgraphs/plus/{replicate_label}.plus.bg"),
        bg_minus = temp("output/bedgraphs/minus/{replicate_label}.plus.bg"),
        bw_plus = "output/bigwigs/plus/{replicate_label}.plus.bw",
        bw_minus = "output/bigwigs/minus/{replicate_label}.minus.bw",
    params:
        error_file = "stderr/{replicate_label}.make_bigwig.err",
        out_file = "stdout/{replicate_label}.make_bigwig.out",
        run_time = "40:00",
        memory = "1000",
        job_name = "make_bigwig"
    benchmark: "benchmarks/bigwigs/unassigned_experiment.{replicate_label}.make_bigwig.txt"
    shell:
        "bedtools genomecov -5 -strand + -bg -ibam {input.bam} | awk '($1 != \"chrEBV\")' | sort -k1,1 -k2,2n > {output.bg_plus};"
        "bedtools genomecov -5 -strand - -bg -ibam {input.bam} | awk '($1 != \"chrEBV\")' | sort -k1,1 -k2,2n > {output.bg_minus};"
        "{TOOL_DIR}/bedGraphToBigWig {output.bg_plus} {CHROM_SIZES} {output.bw_plus};" 
        "{TOOL_DIR}/bedGraphToBigWig {output.bg_minus} {CHROM_SIZES} {output.bw_minus};" 

rule uniq_repeats:
    input:
        repeatmasker = REPEAT_TABLE,
        genome = GENOME
    output:
        sorted_bed = temp("repeats.sort.temp.bed.gz"),
        unique_repeats = REPEAT_TABLE.replace(".tsv", ".sort.unique.bed")
    params:
        error_file = "stderr/calc_partition_nuc.err",
        out_file = "stdout/calc_partition_nuc.out",
        run_time = "40:00",
        memory = "1000",
        job_name = "uniq_repeats_nuc"
    benchmark: "benchmarks/uniq_repeats.txt"
    shell:
        "zcat {REPEAT_TABLE} | awk -v OFS=\"\\t\" '{{print $6,$7,$8,$11 \":\" name_count[$11]++, $2, $10,$11,$12,$13}} "
            "$13 == \"L1\" || $13 == \"Alu\" {{$11 = $11 \"_AS\"; $12 = $12 \"_AS\"; $13 = $13 \"_AS\"; "
            "if($10 == \"+\") {{$10 = \"-\"}} else {{$10 = \"+\"}}; print $6,$7,$8,$11 \":\" name_count[$11]++, $2, $10,$11,$12,$13}}' | "
            "tail -n +2 | bedtools sort -i - | gzip > {output.sorted_bed}; "
        "bedtools coverage -s -d -a {output.sorted_bed} -b {output.sorted_bed}  | awk -v OFS=\"\\t\" "
            "'$NF >1 {{print $1,$2+$(NF-1)-1,$2+$(NF-1),$4,$5,$6}}' | "
            "bedtools sort -i - | "
            "bedtools merge -c 4,5,6 -o distinct -s -i - | "
            "bedtools subtract -s -a {output.sorted_bed} -b - | "
            "bedtools nuc -s -fi {input.genome} -bed -  | awk -v OFS=\"\\t\" 'NR > 1 {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11}}' | "
            "gzip -c > {output.unique_repeats}"


rule quantify_repeats:
    input:
        CHROM_SIZES,
        bam = lambda wildcards: replicate_label_to_bams[wildcards.replicate_label],
        repeats = rules.uniq_repeats.output.unique_repeats
    output:
        counts = "output/counts/repeats/vectors/{replicate_label}.counts"
    params:
        error_file = "stderr/{replicate_label}.quantify_repeats.err",
        out_file = "stdout/{replicate_label}.quantify_repeats.out",
        run_time = "15:00",
        memory = "20000",
        job_name = "dedup_bam",
        prefix='output/bams/dedup/genome/{replicate_label}.genome.sort'
    benchmark: "benchmarks/repeats/unassigned_experiment.{replicate_label}.quantify_repeats.txt"
    shell:
        "bedtools bamtobed -i {input.bam} | awk '($1 != \"chrEBV\") && ($4 !~ \"/{UNINFORMATIVE_READ}$\")' | "
            "bedtools flank -s -l 1 -r 0 -g {CHROM_SIZES} -i - | "
            "bedtools shift -p 1 -m -1 -g {CHROM_SIZES} -i - | "
            "bedtools sort -i - | "
            "bedtools coverage -s -counts -a {input.repeats} -b - | "
            "awk 'BEGIN {{print \"{wildcards.replicate_label}\"}} {{print $NF}}' > {output.counts}"

rule make_repeat_count_tables:
    input:
        unique_repeats = rules.uniq_repeats.output.unique_repeats,
        replicate_counts = lambda wildcards: expand("output/counts/repeats/vectors/{replicate_label}.counts", replicate_label = experiment_to_replicate_labels[wildcards.experiment_label]),
    output:
        name_table = "output/counts/repeats/tables/name/{experiment_label}.tsv.gz",
        class_table = "output/counts/repeats/tables/class/{experiment_label}.tsv.gz",
        family_table = "output/counts/repeats/tables/family/{experiment_label}.tsv.gz",
    params:
        error_file = "stderr/{experiment_label}.make_repeat_count_tables.err",
        out_file = "stdout/{experiment_label}.make_repeat_count_tables.out",
        run_time = "00:15:00",
        cores = "1",
        memory = "200",
        job_name = "make_repeat_count_tables"
    benchmark: "benchmarks/counts/{experiment_label}.all_replicates.make_repeat_count_table.txt"
    shell:
        "echo \"repeat_name\" | paste - {input.replicate_counts} | sed -n '1p' | gzip > {output.name_table};"
        "echo \"repeat_class\" | paste - {input.replicate_counts} | sed -n '1p' | gzip > {output.class_table};"
        "echo \"repeat_family\" | paste - {input.replicate_counts} | sed -n '1p' | gzip > {output.family_table};"
        "paste <(zcat {input.unique_repeats} | awk -v OFS=\"\\t\" 'BEGIN {{print \"repeat_name\";}} {{print $7}}') {input.replicate_counts} | "
            "awk -v OFS=\"\\t\" 'NR > 1 {{for(i = 2; i <= NF; i++) {{tabulation[$1][i] += $i}} }} END {{for(name in tabulation) {{ printf name; for(i = 2; i <= NF; i++) {{printf \"\\t\" tabulation[name][i]}} print \"\";}} }}' | sort -k 1,1 | gzip >> {output.name_table};"
        "paste <(zcat {input.unique_repeats} | awk -v OFS=\"\\t\" 'BEGIN {{print \"repeat_class\";}} {{print $8}}') {input.replicate_counts} | "
            "awk -v OFS=\"\\t\" 'NR > 1 {{for(i = 2; i <= NF; i++) {{tabulation[$1][i] += $i}} }} END {{for(name in tabulation) {{ printf name; for(i = 2; i <= NF; i++) {{printf \"\\t\" tabulation[name][i]}} print \"\";}} }}' | sort -k 1,1 | gzip >> {output.class_table};"
        "paste <(zcat {input.unique_repeats} | awk -v OFS=\"\\t\" 'BEGIN {{print \"repeat_family\";}} {{print $9}}') {input.replicate_counts} | "
            "awk -v OFS=\"\\t\" 'NR > 1 {{for(i = 2; i <= NF; i++) {{tabulation[$1][i] += $i}} }} END {{for(name in tabulation) {{ printf name; for(i = 2; i <= NF; i++) {{printf \"\\t\" tabulation[name][i]}} print \"\";}} }}' | sort -k 1,1 | gzip >> {output.family_table};"

rule fit_clip_betabinomial_re_model:
    input:
        table = "output/counts/repeats/tables/name/{experiment_label}.tsv.gz",
    output:
        coef = "output/clip_model_coef_re/{experiment_label}.{clip_replicate_label}.tsv",
        # plot = lambda wildcards: expand("output/figures/clip_distributions/{{experiment_label}}.{{clip_replicate_label}}.{other_label}.clip_distribution.pdf", other_label = experiment_to_input_replicate_labels[wildcards.experiment_label][wildcards.Input_replicate_label])
    params:
        error_file = "stderr/{experiment_label}.{clip_replicate_label}.fit_clip_betabinomial_re_model.err",
        out_file = "stdout/{experiment_label}.{clip_replicate_label}.fit_clip_betabinomial_re_model.out",
        run_time = "10:00",
        memory = "2000",
        job_name = "fit_clip_betabinomial_re_model"
    benchmark: "benchmarks/fit_clip_betabinomial_re_model/{experiment_label}.{clip_replicate_label}.fit_clip.txt"
    shell:
        "{R_EXE} --vanilla {TOOL_DIR}/fit_clip_betabinom_re.R {input.table} {wildcards.experiment_label} {wildcards.clip_replicate_label}"

rule fit_input_betabinomial_re_model:
    input:
        table = "output/counts/repeats/tables/name/{experiment_label}.tsv.gz",
    output:
        coef = "output/input_model_coef_re/{experiment_label}.{input_replicate_label}.tsv",
        # plot = lambda wildcards: expand("output/figures/input_distributions/{{experiment_label}}.{{input_replicate_label}}.{other_label}.input_distribution.pdf", other_label = experiment_to_input_replicate_labels[wildcards.experiment_label][wildcards.Input_replicate_label])
    params:
        error_file = "stderr/{experiment_label}.{input_replicate_label}.fit_input_betabinomial_re_model.err",
        out_file = "stdout/{experiment_label}.{input_replicate_label}.fit_input_betabinomial_re_model.out",
        run_time = "10:00",
        memory = "2000",
        job_name = "fit_input_betabinomial_re_model"
    benchmark: "benchmarks/fit_input_betabinomial_re_model/{experiment_label}.{input_replicate_label}.fit_input.txt"
    shell:
        "{R_EXE} --vanilla {TOOL_DIR}/fit_input_betabinom_re.R {input.table} {wildcards.experiment_label} {wildcards.input_replicate_label}"

rule call_enriched_re:
    input:
        table = "output/counts/repeats/tables/name/{experiment_label}.tsv.gz",
        repeats = rules.uniq_repeats.output.unique_repeats,
        parameters = lambda wildcards: "output/" + OVERDISPERSION_MODE + "_model_coef_re/{experiment_label}." + overdispersion_replicate_lookup[wildcards.clip_replicate_label] + ".tsv",
    output:
        "output/figures/clip_scatter_re/{experiment_label}.{clip_replicate_label}.clip_test_distribution.pdf",
        "output/enriched_re/{experiment_label}.{clip_replicate_label}.enriched_re.tsv.gz"
    params:
        input_replicate_label = lambda wildcards: clip_to_input_replicate_label[wildcards.clip_replicate_label],
        error_file = "stderr/{experiment_label}.{clip_replicate_label}.call_enriched_re.err",
        out_file = "stdout/{experiment_label}.{clip_replicate_label}.call_enriched_re.out",
        run_time = "00:25:00",
        memory = "3000",
        job_name = "call_enriched_re"
    benchmark: "benchmarks/call_enriched_re/{experiment_label}.{clip_replicate_label}.call_enriched_re.txt"
    shell:
        "{R_EXE} --vanilla {TOOL_DIR}/call_enriched_re.R {input.table} {input.repeats} {input.parameters} {params.input_replicate_label} {wildcards.clip_replicate_label} {wildcards.experiment_label}.{wildcards.clip_replicate_label}"

rule find_reproducible_enriched_re:
    input:
        windows = lambda wildcards: expand("output/enriched_re/{{experiment_label}}.{clip_replicate_label}.enriched_re.tsv.gz", clip_replicate_label = experiment_to_clip_replicate_labels[wildcards.experiment_label])
    output:
        reproducible_windows = "output/reproducible_enriched_re/{experiment_label}.reproducible_enriched_re.tsv.gz",
    params:
        error_file = "stderr/{experiment_label}.find_reproducible_enriched_re.err",
        out_file = "stdout/{experiment_label}.find_reproducible_enriched_re.out",
        run_time = "5:00",
        memory = "1000",
        job_name = "find_reproducible_enriched_re"
    benchmark: "benchmarks/find_reproducible_enriched_re/{experiment_label}.all_replicates.reproducible.txt"
    shell:
        "{R_EXE} --vanilla {TOOL_DIR}/identify_reproducible_re.R output/enriched_re/ {wildcards.experiment_label}"
        
rule partition_bam_reads:
    input:
        CHROM_SIZES,
        bam = lambda wildcards: replicate_label_to_bams[wildcards.replicate_label],
        region_partition = PARTITION,
    output:
        counts = "output/counts/genome/vectors/{replicate_label}.counts",
    params:
        error_file = "stderr/{replicate_label}.partition_bam_reads.err",
        out_file = "stdout/{replicate_label}.partition_bam_reads.out",
        run_time = "20:00",
        cores = "1",
        memory = "10000",
        job_name = "partition_bam_reads"
    benchmark: "benchmarks/counts/unassigned_experiment.{replicate_label}.partition_bam_reads.txt"
    shell:
        "bedtools bamtobed -i {input.bam} | awk '($1 != \"chrEBV\") && ($4 !~ \"/{UNINFORMATIVE_READ}$\")' | "
        "bedtools flank -s -l 1 -r 0 -g {CHROM_SIZES} -i - | "
        "bedtools shift -p 1 -m -1 -g {CHROM_SIZES} -i - | "
        "bedtools sort -i - | "
        "bedtools coverage -counts -s -a {input.region_partition} -b - | cut -f 7 | "
        "awk 'BEGIN {{print \"{wildcards.replicate_label}\"}} {{print}}' > {output.counts};"
        
rule calc_partition_nuc:
    input:
        partition = PARTITION,
        genome = GENOME
    output:
        nuc = PARTITION.replace(".bed", ".nuc")
    params:
        error_file = "stderr/calc_partition_nuc.err",
        out_file = "stdout/calc_partition_nuc.out",
        run_time = "00:10:00",
        memory = "1000",
        job_name = "calc_partition_nuc"
    benchmark: "benchmarks/partition_nuc.txt"
    shell:
        "bedtools nuc -s -fi {input.genome} -bed {input.partition} | gzip -c > {output.nuc}"

rule make_genome_count_table:
    input:
        partition = rules.calc_partition_nuc.output.nuc,
        replicate_counts = lambda wildcards: expand("output/counts/genome/vectors/{replicate_label}.counts", replicate_label = experiment_to_replicate_labels[wildcards.experiment_label]),
    output:
        count_table = "output/counts/genome/tables/{experiment_label}.tsv.gz",
    params:
        error_file = "stderr/{experiment_label}.make_count_table.err",
        out_file = "stdout/{experiment_label}.make_count_table.out",
        run_time = "00:05:00",
        cores = "1",
        memory = "200",
        job_name = "make_genome_count_table"
    benchmark: "benchmarks/counts/{experiment_label}.all_replicates.make_genome_count_table.txt"
    shell:
        "paste <(zcat {input.partition} | awk -v OFS=\"\\t\" 'BEGIN {{print \"chr\\tstart\\tend\\tname\\tscore\\tstrand\\tgc\"}} NR > 1 {{print $1,$2,$3,$4,$5,$6,$8}}' ) {input.replicate_counts} | gzip -c > {output.count_table};"

rule fit_input_betabinomial_model:
    input:
        table = "output/counts/genome/tables/{experiment_label}.tsv.gz"
    output:
        coef = "output/input_model_coef/{experiment_label}.{input_replicate_label}.tsv",
        # plot = lambda wildcards: expand("output/figures/input_distributions/{{experiment_label}}.{{input_replicate_label}}.{other_label}.input_distribution.pdf", other_label = experiment_to_input_replicate_labels[wildcards.experiment_label][wildcards.Input_replicate_label])
    params:
        error_file = "stderr/{experiment_label}.{input_replicate_label}.fit_input_betabinom.err",
        out_file = "stdout/{experiment_label}.{input_replicate_label}.fit_input_betabinom.out",
        run_time = "1:00:00",
        memory = "10000",
        job_name = "fit_input_betabinomial_model"
    benchmark: "benchmarks/betabinomial/{experiment_label}.{input_replicate_label}.fit_input.txt"
    shell:
        "{R_EXE} --vanilla {TOOL_DIR}/fit_input_betabinom.R {input.table} {wildcards.experiment_label} {wildcards.input_replicate_label}"

rule fit_clip_betabinomial_model:
    input:
        table = "output/counts/genome/tables/{experiment_label}.tsv.gz"
    output:
        coef = "output/clip_model_coef/{experiment_label}.{clip_replicate_label}.tsv",
        # plot = lambda wildcards: expand("output/figures/clip_distributions/{{experiment_label}}.{{clip_replicate_label}}.{other_label}.clip_distribution.pdf", other_label = experiment_to_input_replicate_labels[wildcards.experiment_label][wildcards.Input_replicate_label])
    params:
        error_file = "stderr/{experiment_label}.{clip_replicate_label}.fit_clip_betabinomial_model.err",
        out_file = "stdout/{experiment_label}.{clip_replicate_label}.fit_clip_betabinomial_model.out",
        run_time = "1:00:00",
        memory = "1000",
        job_name = "fit_clip_betabinomial_model"
    benchmark: "benchmarks/fit_clip_betabinomial_model/{experiment_label}.{clip_replicate_label}.fit_clip.txt"
    shell:
        "{R_EXE} --vanilla {TOOL_DIR}/fit_clip_betabinom.R {input.table} {wildcards.experiment_label} {wildcards.clip_replicate_label}"

rule call_enriched_windows:
    input:
        feature_annotations = FEATURE_ANNOTATIONS,
        accession_rankings = ACCESSION_RANKINGS,
        table = "output/counts/genome/tables/{experiment_label}.tsv.gz",
        parameters = lambda wildcards: "output/" + OVERDISPERSION_MODE + "_model_coef/{experiment_label}." + overdispersion_replicate_lookup[wildcards.clip_replicate_label] + ".tsv",
        # parameters = lambda wildcards: "output/clip_model_coef/{experiment_label}.{wildcards.clip_replicate_label}.tsv",
    output:
        "output/threshold_scan/{experiment_label}.{clip_replicate_label}.threshold_data.tsv",
        "output/tested_windows/{experiment_label}.{clip_replicate_label}.tested_windows.tsv.gz",
        "output/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_windows.tsv.gz",
        "output/enrichment_summaries/{experiment_label}.{clip_replicate_label}.enriched_window_feature_summary.tsv",
        "output/enrichment_summaries/{experiment_label}.{clip_replicate_label}.enriched_window_transcript_summary.tsv",
        "output/enrichment_summaries/{experiment_label}.{clip_replicate_label}.enriched_window_gene_summary.tsv",
        "output/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_fractions_feature_data.tsv",
        "output/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds_feature_data.tsv",
        "output/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds_transcript_data.tsv",
        "output/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds_feature_gc_data.tsv",
        "output/figures/threshold_scan/{experiment_label}.{clip_replicate_label}.threshold_scan.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_coverage.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_rates.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_counts.linear.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_counts.log10.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_odds.feature.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_odds.all_transcript_types.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_odds.select_transcript_types.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_counts.per_gene_feature.pdf",
        "output/figures/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_fractions.feature.pdf",
        "output/figures/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds.feature.pdf",
        "output/figures/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds.all_transcript_types.pdf",
        "output/figures/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds.feature_gc.pdf"
    params:
        input_replicate_label = lambda wildcards: clip_to_input_replicate_label[wildcards.clip_replicate_label],
        error_file = "stderr/{experiment_label}.{clip_replicate_label}.call_enriched_windows.err",
        out_file = "stdout/{experiment_label}.{clip_replicate_label}.call_enriched_windows.out",
        run_time = "00:45:00",
        memory = "6000",
        job_name = "call_enriched_windows"
    benchmark: "benchmarks/call_enriched_windows/{experiment_label}.{clip_replicate_label}.call_enriched_windows.txt"
    shell:
        "{R_EXE} --vanilla {TOOL_DIR}/call_enriched_windows.R {input.table} {input.accession_rankings} {input.feature_annotations} {input.parameters} {params.input_replicate_label} {wildcards.clip_replicate_label} {wildcards.experiment_label}.{wildcards.clip_replicate_label}"

rule check_window_concordance:
    input:
        windows = lambda wildcards: expand("output/tested_windows/{{experiment_label}}.{clip_replicate_label}.tested_windows.tsv.gz", clip_replicate_label = experiment_to_clip_replicate_labels[wildcards.experiment_label])
    output:
        "output/figures/enrichment_reproducibility/{experiment_label}.enrichment_reproducibility.pdf",
        "output/enrichment_reproducibility/{experiment_label}.enrichment_reproducibility.tsv"
    params:
        error_file = "stderr/{experiment_label}.check_window_concordance.err",
        out_file = "stdout/{experiment_label}.check_window_concordance.out",
        run_time = "0:15:00",
        memory = "1000",
        job_name = "check_window_concordance"
    benchmark: "benchmarks/check_window_concordance/{experiment_label}.all_replicates.concordance.txt"
    shell:
        "{R_EXE} --vanilla {TOOL_DIR}/check_window_concordance.R output/tested_windows {wildcards.experiment_label} " + (BLACKLIST if BLACKLIST is not None else "") 

rule find_reproducible_enriched_windows:
    input:
        windows = lambda wildcards: expand("output/enriched_windows/{{experiment_label}}.{clip_replicate_label}.enriched_windows.tsv.gz", clip_replicate_label = experiment_to_clip_replicate_labels[wildcards.experiment_label])
    output:
        reproducible_windows = "output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz",
        linear_bar = "output/figures/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_window_counts.linear.pdf",
        log_bar = "output/figures/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_window_counts.log10.pdf"
    params:
        error_file = "stderr/{experiment_label}.find_reproducible_enriched_windows.err",
        out_file = "stdout/{experiment_label}.find_reproducible_enriched_windows.out",
        run_time = "5:00",
        memory = "2000",
        job_name = "find_reproducible_enriched_windows"
    benchmark: "benchmarks/find_reproducible_enriched_windows/{experiment_label}.all_replicates.reproducible.txt"
    shell:
        "{R_EXE} --vanilla {TOOL_DIR}/identify_reproducible_windows.R output/enriched_windows/ {wildcards.experiment_label} " + (BLACKLIST if BLACKLIST is not None else "") 

rule sample_background_windows_by_region:
    input:
        enriched_windows = "output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz",
        all_windows = FEATURE_ANNOTATIONS,
    output:
        variable_windows = "output/homer/region_matched_background/variable/{experiment_label}.sampled_variable_windows.bed.gz",
        fixed_windows = "output/homer/region_matched_background/fixed/{experiment_label}.sampled_fixed_windows.bed.gz"
    params:
        error_file = "stderr/{experiment_label}.sample_background_windows_by_region.err",
        out_file = "stdout/{experiment_label}.sample_background_windows_by_region.out",
        run_time = "10:00",
        memory = "3000",
        job_name = "sample_background_windows"
    benchmark: "benchmarks/sample_background_windows_by_region/{experiment_label}.sample_background_windows_by_region.txt"
    shell:
        "{R_EXE} --vanilla {TOOL_DIR}/sample_matched_background_by_region.R {input.enriched_windows} {input.all_windows} 75 output/homer/region_matched_background {wildcards.experiment_label};"

rule get_nt_coverage:
    input:
        windows = "output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz",
        clip_bams = lambda wildcards: [replicate_label_to_bams[clip_replicate_label] for clip_replicate_label in experiment_to_clip_replicate_labels[wildcards.experiment_label]],
        input_bams = lambda wildcards: [replicate_label_to_bams[input_replicate_label] for input_replicate_label in experiment_to_input_replicate_labels[wildcards.experiment_label]],         
    output:
        nt_census = temp("output/finemapping/nt_coverage/{experiment_label}.nt_census.bed"),
        nt_input_counts = temp("output/finemapping/nt_coverage/{experiment_label}.nt_coverage.input.counts"),
        nt_clip_counts = temp("output/finemapping/nt_coverage/{experiment_label}.nt_coverage.clip.counts"),
        nt_coverage = "output/finemapping/nt_coverage/{experiment_label}.nt_coverage.bed"
    params:
        error_file = "stderr/{experiment_label}.get_nt_coverage.err",
        out_file = "stdout/{experiment_label}.get_nt_coverage.out",
        run_time = "1:00:00",
        memory = "15000",
        job_name = "get_nt_coverage"
    benchmark: "benchmarks/get_nt_coverage/{experiment_label}.all_replicates.reproducible.txt"
    shell:
        "set +eu;"
        "module load samtools/1.16;"
        "zcat {input.windows} | tail -n +2 | sort -k1,1 -k2,2n | awk -v OFS=\"\t\" '{{print $1, $2 -37, $3+37,$4,$5,$6}}' | "
            "bedtools merge -i - -s -c 6 -o distinct | awk -v OFS=\"\t\" '{{for(i=$2;i< $3;i++) {{print $1,i,i+1,\"MW:\" NR \":\" i - $2,0,$4, NR}} }}' > {output.nt_census}; "
        "samtools cat {input.input_bams} | bedtools intersect -s -wa -a - -b {output.nt_census} | "
            "bedtools bamtobed -i - | awk '($1 != \"chrEBV\") && ($4 !~ \"/{UNINFORMATIVE_READ}$\")' | "
            "bedtools flank -s -l 1 -r 0 -g {CHROM_SIZES} -i - | "
            "bedtools shift -p 1 -m -1 -g {CHROM_SIZES} -i - | "
            "bedtools sort -i - | "
            "bedtools coverage -counts -s -a {output.nt_census} -b - | awk '{{print $NF}}' > {output.nt_input_counts};"
        "samtools cat {input.clip_bams} | bedtools intersect -s -wa -a - -b {output.nt_census} | "
            "bedtools bamtobed -i - | awk '($1 != \"chrEBV\") && ($4 !~ \"/{UNINFORMATIVE_READ}$\")' | "
            "bedtools flank -s -l 1 -r 0 -g {CHROM_SIZES} -i - | "
            "bedtools shift -p 1 -m -1 -g {CHROM_SIZES} -i - | "
            "bedtools sort -i - | "
            "bedtools coverage -counts -s -a {output.nt_census} -b - | awk '{{print $NF}}' > {output.nt_clip_counts};"
        "paste {output.nt_census} {output.nt_input_counts} {output.nt_clip_counts} > {output.nt_coverage}"

rule finemap_windows:
    input:
        nt_coverage = "output/finemapping/nt_coverage/{experiment_label}.nt_coverage.bed",        
    output:
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz"
    params:
        error_file = "stderr/{experiment_label}.finemap_windows.err",
        out_file = "stdout/{experiment_label}.finemap_windows.out",
        run_time = "1:00:00",
        memory = "10000",
        job_name = "finemap_windows"
    benchmark: "benchmarks/finemap_windows/{experiment_label}.all_replicates.reproducible.txt"
    shell:
        "{R_EXE} --vanilla {TOOL_DIR}/finemap_enriched_windows.R {input.nt_coverage} output/finemapping/mapped_sites/ {wildcards.experiment_label}"

rule run_homer:
    input:
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz",
        background = "output/homer/region_matched_background/fixed/{experiment_label}.sampled_fixed_windows.bed.gz",
        genome = GENOME
    output:
        report = "output/homer/finemapped_results/{experiment_label}/homerResults.html"
    params:
        error_file = "stderr/{experiment_label}.run_homer.err",
        out_file = "stdout/{experiment_label}.run_homer.out",
        run_time = "40:00",
        memory = "2000",
        job_name = "run_homer"
    benchmark: "benchmarks/run_homer/{experiment_label}.all_replicates.reproducible.txt"
    shell:
        "findMotifsGenome.pl <(less {input.finemapped_windows} | awk -v OFS=\"\t\" '{{print $4 \":\"$9,$1,$2+1,$3,$6}}') "
            "{input.genome} output/homer/finemapped_results/{wildcards.experiment_label} -preparsedDir output/homer/preparsed -size given -rna -nofacts -S 20 -len 5,6,7,8,9 -nlen 1 "
            "-bg <(zcat {input.background} | awk -v OFS=\"\t\" '{{print $4,$1,$2+1,$3,$6}}') "

rule consult_encode_reference:
    input:
        enriched_windows = lambda wildcards: expand("output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz", experiment_label = experiment_labels),
        enriched_re = lambda wildcards: expand("output/reproducible_enriched_re/{experiment_label}.reproducible_enriched_re.tsv.gz", experiment_label = experiment_labels),
        encode_references = lambda wildcards: expand(TOOL_DIR + "/{reference}.reference.tsv", reference = ["encode3_feature_summary", "encode3_eclip_enrichment", "encode3_class_assignment"])
    output:
        tsne_coordinates = "output/tsne/skipper.tsne_query.tsv",
        tsne_plot = "output/figures/tsne/skipper.tsne_query.pdf"
    params:
        error_file = "stderr/skipper.consult_encode_reference.err",
        out_file = "stdout/skipper.consult_encode_reference.out",
        run_time = "10:00",
        memory = "1000",
        job_name = "consult_encode_reference"
    benchmark: "benchmarks/consult_encode_reference/skipper.txt"
    shell:
        "{R_EXE} --vanilla {TOOL_DIR}/consult_encode_reference.R output/reproducible_enriched_windows output/reproducible_enriched_re {TOOL_DIR} skipper "

rule consult_term_reference:
    input:
        enriched_windows = "output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz",
        gene_sets = GENE_SETS,
        gene_set_reference = GENE_SET_REFERENCE,
        gene_set_distance = GENE_SET_DISTANCE
    output:
        enrichment_results = "output/gene_sets/{experiment_label}.enriched_terms.tsv.gz",
        enrichment_plot = "output/figures/gene_sets/{experiment_label}.clustered_top_terms.pdf"
    params:
        error_file = "stderr/{experiment_label}.consult_term_reference.err",
        out_file = "stdout/{experiment_label}.consult_term_reference.out",
        run_time = "15:00",
        memory = "1000",
        job_name = "consult_term_reference"
    benchmark: "benchmarks/consult_term_reference/{experiment_label}.all_replicates.reproducible.txt"
    shell:
        "mkdir -p output/gene_sets/;{R_EXE} --vanilla {TOOL_DIR}/consult_term_reference.R {input.enriched_windows} {input.gene_sets} {input.gene_set_reference} {input.gene_set_distance} {wildcards.experiment_label} "
